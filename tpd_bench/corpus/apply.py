"""Replace hand-entered scaffold features with corpus-derived ones.

This closes the loop: corpus -> features -> scoring. It also exposes drift.
Where the corpus disagrees with what was typed in by hand, the corpus wins and
the difference is printed, because a hand-entered value the corpus cannot
justify was doing work it had not earned.
"""

from __future__ import annotations

import csv
import dataclasses
import pathlib

from ..schema import ScaffoldSign, Target

CSV = pathlib.Path(__file__).with_name("scaffold_evidence.csv")

SIGN = {
    "pathogenic": ScaffoldSign.PATHOGENIC,
    "protective": ScaffoldSign.PROTECTIVE,
    "unknown": ScaffoldSign.UNKNOWN,
}


def load() -> dict[str, tuple[int, ScaffoldSign]]:
    if not CSV.exists():
        # Generated, so a fresh clone does not have it. Build it rather than
        # making every entry point fail on a missing file.
        import io
        from contextlib import redirect_stdout

        from .build import run as build_corpus

        with redirect_stdout(io.StringIO()):
            build_corpus()
    with CSV.open(encoding="utf-8") as fh:
        return {
            r["gene"]: (int(r["tier"]), SIGN[r["sign"]])
            for r in csv.DictReader(fh)
        }


def apply(targets: list[Target], verbose: bool = True) -> list[Target]:
    corpus = load()
    out, drift = [], []
    for t in targets:
        hit = corpus.get(t.gene)
        if not hit:
            # No corpus entry means no evidence, which is not the same as the
            # tier someone typed in. Clear it.
            if t.scaffold_evidence_tier:
                drift.append((t.gene, f"tier {t.scaffold_evidence_tier} -> 0 (no corpus entry)"))
                t = dataclasses.replace(
                    t, scaffold_evidence_tier=0, scaffold_sign=ScaffoldSign.NONE
                )
            out.append(t)
            continue
        tier, sign = hit
        if (tier, sign) != (t.scaffold_evidence_tier, t.scaffold_sign):
            drift.append((
                t.gene,
                f"tier {t.scaffold_evidence_tier}->{tier}, "
                f"sign {t.scaffold_sign.value}->{sign.value}",
            ))
        out.append(dataclasses.replace(t, scaffold_evidence_tier=tier, scaffold_sign=sign))

    if verbose and drift:
        print("corpus drift (corpus overrides hand-entry):")
        for gene, msg in drift:
            print(f"  {gene:10s} {msg}")
        print()
    return out
