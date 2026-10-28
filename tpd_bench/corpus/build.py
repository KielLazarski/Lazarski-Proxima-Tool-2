"""Build the corpus and report what it found."""

from __future__ import annotations

import pathlib

from .discordance import build, write_csv
from .seed import COMPARISONS, RESISTANCE

OUT = pathlib.Path(__file__).with_name("scaffold_evidence.csv")


def run() -> int:
    entries = build(RESISTANCE, COMPARISONS)
    write_csv(entries, OUT)

    print("=" * 74)
    print("SCAFFOLDING DISCORDANCE CORPUS")
    print("=" * 74)
    for e in entries:
        mark = "REJECT" if e.sign == "protective" else ""
        print(f"  tier {e.tier}  {e.gene:8s} sign={e.sign:11s} {mark}")
        print(f"          {e.basis}")
        if e.confounders:
            print(f"          demoted by: {', '.join(e.confounders)}")
    print()
    n1 = sum(1 for e in entries if e.tier == 1)
    prot = [e.gene for e in entries if e.sign == "protective"]
    print(f"{len(entries)} genes | {n1} at tier 1 | {len(prot)} rejected on sign: {prot}")
    print(f"written to {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
