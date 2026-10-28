"""Assemble Target records for the whole universe, then score them.

Offline this reports what it would do and stops at the first blocked source.
It never emits a partial universe that looks complete -- a sweep missing a
feature layer produces confident nominations from absent data, which is worse
than no sweep.
"""

from __future__ import annotations

import sys

from .universe import NetworkUnavailable, classify_location, human_proteome

LAYERS = [
    ("universe", "UniProt", "subcellular location, sequence, active sites"),
    ("genetics", "Open Targets", "association class; DIRECTION still unresolved"),
    ("expression", "HPA / GTEx", "disease vs normal fold change"),
    ("pharmacology", "ChEMBL", "ligandability and inhibitor status (weakest layer)"),
    ("competition", "ClinicalTrials.gov", "most advanced competing modality"),
    ("structure", "AlphaFold DB", "models for the degron scorer"),
    ("turnover", "proteomics atlases", "half-life; patchy, needs inference"),
]


def preflight() -> int:
    print("=" * 74)
    print("PROTEOME SWEEP -- preflight")
    print("=" * 74)
    for name, source, what in LAYERS:
        print(f"  {name:14s} {source:22s} {what}")
    print()
    try:
        n = 0
        placed = 0
        for entry in human_proteome():
            n += 1
            if classify_location(entry):
                placed += 1
            if n >= 500:
                break
        print(f"reachable: sampled {n} entries, {placed} with a usable location")
        return 0
    except NetworkUnavailable as exc:
        print(f"BLOCKED: {exc}")
        print()
        print("Nothing was written. Run this where the open-source APIs are")
        print("reachable; the bench continues to use the 36 curated records")
        print("until a sweep completes.")
        return 1


if __name__ == "__main__":
    sys.exit(preflight())
