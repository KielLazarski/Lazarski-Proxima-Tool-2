"""Shared HTTP layer.

Two things the first real run on a corporate Windows machine exposed:

1. POSTing JSON without a Content-Type header. Open Targets answered 415 and
   gnomAD answered 400. Both are GraphQL endpoints that require the header;
   the sandbox never reached them, so the bug survived until now.

2. Corporate TLS interception. A company proxy re-signs HTTPS with its own
   root certificate, and Python on Windows does not reliably pick that up,
   so verification fails with "self-signed certificate in certificate chain".
   The certificate is legitimately installed in the Windows store -- Python
   just is not looking there.

   The fix loads the Windows certificate store explicitly. It never disables
   verification: turning verification off to make an error go away is how you
   end up trusting anything on the network.
"""

from __future__ import annotations

import json
import os
import pathlib
import platform
import ssl
import urllib.error
import urllib.request
from typing import Any

TIMEOUT = 30
_CTX: ssl.SSLContext | None = None


def ssl_context() -> ssl.SSLContext:
    """Default trust, plus the Windows store when running on Windows."""
    global _CTX
    if _CTX is not None:
        return _CTX

    ctx = ssl.create_default_context()
    if platform.system() == "Windows":
        loaded = 0
        for store in ("ROOT", "CA"):
            try:
                certs = ssl.enum_certificates(store)
            except Exception:
                continue
            for cert, encoding, _trust in certs:
                if encoding != "x509_asn":
                    continue
                try:
                    ctx.load_verify_locations(cadata=ssl.DER_cert_to_PEM_cert(cert))
                    loaded += 1
                except ssl.SSLError:
                    pass  # one unparseable cert should not sink the whole store
        ctx._tpd_windows_certs = loaded  # type: ignore[attr-defined]

    # Escape hatch for a corporate CA that is not in the Windows store: ask IT
    # for the proxy's root certificate as a .pem and point TPD_EXTRA_CA at it.
    #
    #     set TPD_EXTRA_CA=C:\path\to\company-root.pem
    #
    # This ADDS a trusted root. It never turns verification off -- a bench that
    # trusts anything on the network is a bench you cannot trust the output of.
    extra = os.environ.get("TPD_EXTRA_CA")
    if extra:
        p = pathlib.Path(extra)
        if p.is_file():
            try:
                ctx.load_verify_locations(cafile=str(p))
                ctx._tpd_extra_ca = str(p)  # type: ignore[attr-defined]
            except ssl.SSLError as exc:
                ctx._tpd_extra_ca = f"FAILED to load {p}: {exc}"  # type: ignore[attr-defined]
        else:
            ctx._tpd_extra_ca = f"TPD_EXTRA_CA set but not a file: {extra}"  # type: ignore[attr-defined]

    _CTX = ctx
    return ctx


def request(url: str, payload: dict | None = None) -> Any:
    """GET, or POST JSON when payload is given. Returns parsed JSON."""
    data = None
    headers = {"User-Agent": "tpd-bench/1.0", "Accept": "application/json"}
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"

    req = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(req, timeout=TIMEOUT, context=ssl_context()) as r:
        return json.load(r)
