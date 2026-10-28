"""Check every data source and validate the assumptions made about its schema.

Reachability is the easy half. The half that actually breaks a sweep is a
source that responds fine while the field this code reads has been renamed,
nested differently, or quietly emptied -- that produces a clean-looking run
built on nothing.

So each check fetches one real record and asserts the specific fields the
adapters depend on. Run it, paste the report back, and the failures say exactly
what to fix.

    python3 -m tpd_bench.features.doctor

Python 3.9+. No dependencies, no installs -- standard library only.
"""

from __future__ import annotations

import json
import platform
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field

from . import http

TIMEOUT = 30


@dataclass
class Check:
    name: str
    url: str
    ok: bool = False
    reachable: bool = False
    detail: str = ""
    missing: list[str] = field(default_factory=list)


def _fetch(url: str, payload: dict | None = None) -> tuple[bool, object, str]:
    try:
        return True, http.request(url, payload), "HTTP 200"
    except urllib.error.HTTPError as e:
        body = ""
        try:
            body = " | " + e.read().decode()[:160].replace("\n", " ")
        except Exception:
            pass
        return False, None, f"HTTP {e.code}{body}"
    except Exception as e:  # noqa: BLE001 - report anything, never crash the doctor
        return False, None, f"{type(e).__name__}: {e}"


def _dig(obj, *path):
    for p in path:
        if obj is None:
            return None
        obj = obj[p] if isinstance(p, int) else (obj or {}).get(p)
    return obj


def _identify_interceptor(host: str) -> str:
    """Name whoever is re-signing this connection.

    Deliberately reads the certificate WITHOUT trusting it, and transfers no
    data over the connection. The point is to tell you which certificate to
    ask IT for -- not to work around the failure.
    """
    import socket

    try:
        ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
        with socket.create_connection((host, 443), timeout=10) as sock:
            with ctx.wrap_socket(sock, server_hostname=host) as ss:
                der = ss.getpeercert(binary_form=True)
        if not der:
            return "could not read the presented certificate"
        txt = ssl.DER_cert_to_PEM_cert(der)
        # Parse the issuer out of the PEM without a third-party x509 library
        import base64
        raw = base64.b64decode("".join(
            l for l in txt.splitlines() if "CERTIFICATE" not in l
        ))
        printable = "".join(chr(b) if 32 <= b < 127 else "\x00" for b in raw)
        words = [w for w in printable.split("\x00") if len(w) >= 4]
        head = words[:14]
        return "certificate presented by: " + " / ".join(head[:8])
    except Exception as exc:  # noqa: BLE001
        return f"could not inspect the certificate ({type(exc).__name__})"


# ------------------------------------------------------------------ checks

def check_uniprot() -> Check:
    url = (
        "https://rest.uniprot.org/uniprotkb/search?query=accession:P10275"
        "&fields=accession,gene_primary,cc_subcellular_location,sequence&format=json&size=1"
    )
    c = Check("UniProt (universe, compartment)", url)
    ok, data, detail = _fetch(url)
    c.reachable, c.detail = ok, detail
    if not ok:
        return c
    res = _dig(data, "results", 0)
    if not res:
        c.missing.append("results[0]")
        return c
    has_loc = any(
        cm.get("commentType") == "SUBCELLULAR LOCATION" for cm in res.get("comments", [])
    )
    if not has_loc:
        c.missing.append("comments[].commentType == 'SUBCELLULAR LOCATION'")
    if not _dig(res, "sequence", "value"):
        c.missing.append("sequence.value")
    c.ok = not c.missing
    return c


def check_opentargets() -> Check:
    url = "https://api.platform.opentargets.org/api/v4/graphql"
    c = Check("Open Targets (genetics)", url)
    q = {
        "query": "query($g:String!){target(ensemblId:$g){approvedSymbol "
                 "associatedDiseases{rows{score datatypeScores{id score}}}}}",
        "variables": {"g": "ENSG00000169047"},
    }
    ok, data, detail = _fetch(url, payload=q)
    c.reachable, c.detail = ok, detail
    if not ok:
        return c
    rows = _dig(data, "data", "target", "associatedDiseases", "rows")
    if rows is None:
        c.missing.append("data.target.associatedDiseases.rows")
    elif rows and "datatypeScores" not in rows[0]:
        c.missing.append("rows[].datatypeScores")
    c.ok = not c.missing
    return c


def check_chembl() -> Check:
    url = (
        "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
        "?target_chembl_id=CHEMBL5255&limit=1"
    )
    c = Check("ChEMBL (ligandability, inhibitor status)", url)
    ok, data, detail = _fetch(url)
    c.reachable, c.detail = ok, detail
    if not ok:
        return c
    mechs = _dig(data, "mechanisms")
    if mechs is None:
        c.missing.append("mechanisms")
    elif mechs and "max_phase" not in mechs[0]:
        c.missing.append("mechanisms[].max_phase")
    c.ok = not c.missing
    return c


def check_clinicaltrials() -> Check:
    url = (
        "https://clinicaltrials.gov/api/v2/studies?"
        + urllib.parse.urlencode({"query.cond": "atopic dermatitis", "pageSize": "1"})
    )
    c = Check("ClinicalTrials.gov v2 (competition)", url)
    ok, data, detail = _fetch(url)
    c.reachable, c.detail = ok, detail
    if not ok:
        if "CERTIFICATE_VERIFY" in detail:
            c.detail += " | " + _identify_interceptor("clinicaltrials.gov")
        return c
    s = _dig(data, "studies", 0)
    if s is None:
        c.missing.append("studies[0]")
    elif _dig(s, "protocolSection", "designModule") is None:
        c.missing.append("protocolSection.designModule (phases live here)")
    c.ok = not c.missing
    return c


def check_alphafold() -> Check:
    url = "https://alphafold.ebi.ac.uk/api/prediction/P10275"
    c = Check("AlphaFold DB (degron scoring)", url)
    ok, data, detail = _fetch(url)
    c.reachable, c.detail = ok, detail
    if not ok:
        return c
    rec = data[0] if isinstance(data, list) and data else None
    if not rec:
        c.missing.append("prediction[0]")
    else:
        # the degron scorer needs coordinates and per-residue confidence
        if not (rec.get("cifUrl") or rec.get("pdbUrl")):
            c.missing.append("cifUrl or pdbUrl")
        c.detail += " | model available"
    c.ok = not c.missing
    return c


def check_gnomad() -> Check:
    """Needed for direction of effect, the gap that currently discounts all genetics."""
    url = "https://gnomad.broadinstitute.org/api"
    c = Check("gnomAD (direction of effect -- currently missing)", url)
    # The schema exposes gnomad_constraint (and legacy exac_constraint); plain
    # "constraint" does not exist. Try the modern field, fall back to legacy.
    for field in ("gnomad_constraint", "exac_constraint"):
        q = {
            "query": "query($g:String!){gene(gene_symbol:$g,reference_genome:GRCh38)"
                     "{gene_id %s{pLI oe_lof}}}" % field,
            "variables": {"g": "PCSK9"},
        }
        ok, data, detail = _fetch(url, payload=q)
        c.reachable, c.detail = ok, f"{detail} (field: {field})"
        if ok and _dig(data, "data", "gene", field) is not None:
            c.ok = True
            return c
    if ok:
        c.missing.append("gnomad_constraint and exac_constraint both absent/empty")
    return c


CHECKS = [
    check_uniprot,
    check_opentargets,
    check_chembl,
    check_clinicaltrials,
    check_alphafold,
    check_gnomad,
]


def run() -> int:
    print("=" * 74)
    print("TPD BENCH -- SOURCE DOCTOR")
    print("=" * 74)
    print(f"python {platform.python_version()} on {platform.system()}")
    ctx = http.ssl_context()
    n_certs = getattr(ctx, "_tpd_windows_certs", None)
    if n_certs is not None:
        print(f"loaded {n_certs} certificates from the Windows store")
    extra = getattr(ctx, "_tpd_extra_ca", None)
    print(f"TPD_EXTRA_CA: {extra}" if extra else "TPD_EXTRA_CA: not set")
    if sys.version_info < (3, 9):
        print("!! needs Python 3.9 or newer")
        return 2
    print()

    results = [fn() for fn in CHECKS]
    for c in results:
        if c.ok:
            mark, note = "OK  ", c.detail
        elif c.reachable:
            mark, note = "SCHEMA", f"{c.detail} | missing: {', '.join(c.missing)}"
        else:
            mark, note = "FAIL", c.detail
        print(f"  [{mark:6s}] {c.name}")
        print(f"            {note}")

    good = sum(c.ok for c in results)
    print()
    print(f"{good}/{len(results)} sources usable")
    if good == len(results):
        print("\nAll clear. Next: python3 -m tpd_bench.features.pipeline")
    else:
        print(
            "\nPaste this whole report back. SCHEMA failures mean the source answered\n"
            "but a field moved -- those are quick fixes. FAIL means blocked or down."
        )
    return 0 if good == len(results) else 1


if __name__ == "__main__":
    sys.exit(run())
