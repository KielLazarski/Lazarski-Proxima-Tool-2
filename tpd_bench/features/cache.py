"""On-disk response cache, so a sweep survives interruption.

A 12k-protein sweep is tens of thousands of requests against rate-limited
public APIs. It will be interrupted. Without a cache every restart begins from
zero, which in practice means it never finishes.

Keyed on the URL plus the request body. Deliberately never expires: a sweep
should be reproducible, and a cache that silently refreshes mid-run produces a
dataset assembled from two different days. Delete the directory to refresh.
"""

from __future__ import annotations

import hashlib
import json
import pathlib

DIR = pathlib.Path(__file__).parent / ".cache"


def _key(url: str, body: bytes | None) -> pathlib.Path:
    h = hashlib.sha256(url.encode() + (body or b"")).hexdigest()[:32]
    return DIR / h[:2] / f"{h}.json"


def get(url: str, body: bytes | None = None):
    p = _key(url, body)
    if p.exists():
        try:
            return json.loads(p.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            p.unlink(missing_ok=True)  # corrupt entry, refetch
    return None


def put(url: str, data, body: bytes | None = None) -> None:
    p = _key(url, body)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(data), encoding="utf-8")
    tmp.replace(p)  # atomic, so a kill mid-write cannot leave a partial entry


def stats() -> tuple[int, float]:
    if not DIR.exists():
        return 0, 0.0
    files = list(DIR.rglob("*.json"))
    return len(files), sum(f.stat().st_size for f in files) / 1e6
