"""Calibration harness.

Scores every target from features alone, compares against the held labels, and
prints every disagreement. The disagreements are the product -- a report that
says 100% is a report that was fitted, not calibrated.
"""

from __future__ import annotations

import csv
import pathlib
import sys
from collections import Counter

from .scoring import Verdict, score
from .corpus.apply import apply as apply_corpus
from .targets import BY_ID, TARGETS as RAW_TARGETS

# Scaffold features come from the corpus, never from hand-entry.
TARGETS = apply_corpus(RAW_TARGETS, verbose=False)

LABELS = pathlib.Path(__file__).with_name("calibration_labels.csv")

# Verdicts that all mean "do not take this forward". Distinguishing between
# them matters for diagnosis, but a reject-for-the-right-family still counts
# as directionally correct, so both are reported.
REJECT_FAMILY = {
    Verdict.REJECT_TARGET_INVALID,
    Verdict.REJECT_WRONG_MODALITY,
    Verdict.REJECT_INHIBITION_SUFFICIENT,
    Verdict.REJECT_SCAFFOLD_PROTECTIVE,
    Verdict.REJECT_MODALITY_TOO_WEAK,
}


def load_labels() -> dict[str, tuple[str, str]]:
    with LABELS.open(encoding="utf-8") as fh:
        return {r["target_id"]: (r["expected_verdict"], r["rationale"]) for r in csv.DictReader(fh)}


def run() -> int:
    labels = load_labels()
    missing = set(labels) ^ {t.target_id for t in TARGETS}
    if missing:
        print(f"!! label/target mismatch: {sorted(missing)}")
        return 2

    exact = directional = 0
    disagreements = []
    rows = []

    for t in TARGETS:
        r = score(t)
        expected, rationale = labels[t.target_id]
        got = r.verdict.value

        is_exact = got == expected
        both_reject = r.verdict in REJECT_FAMILY and expected in {v.value for v in REJECT_FAMILY}
        is_directional = is_exact or both_reject

        exact += is_exact
        directional += is_directional
        rows.append((t, r, expected, is_exact, is_directional))
        if not is_exact:
            disagreements.append((t, r, expected, rationale, both_reject))

    n = len(TARGETS)
    print("=" * 78)
    print(f"TPD GLUE BENCH -- CALIBRATION  ({n} cases)")
    print("=" * 78)
    print(f"exact verdict match     {exact}/{n}  ({exact/n:.0%})")
    print(f"directionally correct   {directional}/{n}  ({directional/n:.0%})")
    print()

    # per-area
    print("by disease area")
    areas = sorted({t.disease_area for t in TARGETS})
    for a in areas:
        sub = [x for x in rows if x[0].disease_area == a]
        e = sum(1 for x in sub if x[3])
        print(f"  {a:26s} {e}/{len(sub)} exact")
    print()

    # Two ranked lists, never one. The tracks are judged against different
    # delivery bars, so ranking them against each other would be meaningless.
    from .schema import Track
    keep = (Verdict.TOP, Verdict.MID, Verdict.NOMINATE_FOR_VALIDATION)
    for track in (Track.INTRACELLULAR, Track.EXTRACELLULAR):
        print(f"{track.value} track -- nominations")
        noms = sorted(
            [x for x in rows if x[1].track is track and x[1].verdict in keep],
            key=lambda x: -x[1].composite,
        )
        if not noms:
            print("  none")
        for t, r, *_ in noms:
            print(
                f"  {r.composite:5.1f}  {r.gene:12s} {t.indication:28s} "
                f"V={r.validity.score:5.1f} M={r.modality.score:5.1f}  {r.verdict.value}"
            )
            for f in r.flags:
                print(f"         ! {f}")
        print()

    # E3 opportunities are reported, never scored.
    opps = [t for t in TARGETS if t.e3_restricted_candidate]
    if opps:
        print("restricted-E3 ligand-discovery opportunities (not scored)")
        for t in opps:
            print(f"  {t.gene:12s} {t.e3_restricted_candidate}")
        print()

    print(f"disagreements ({len(disagreements)})")
    if not disagreements:
        print("  none -- which is itself suspicious; check that labels are not leaking into features")
    for t, r, expected, rationale, both_reject in disagreements:
        kind = "same family" if both_reject else "DIRECTIONAL MISS"
        print(f"\n  [{kind}] {t.target_id} {r.gene} / {t.indication}")
        print(f"    expected : {expected}")
        print(f"    got      : {r.verdict.value}   (V={r.validity.score} M={r.modality.score})")
        print(f"    label says: {rationale}")
        for f in r.validity.flags:
            print(f"      V: {f}")
        for f in r.modality.flags:
            print(f"      M: {f}")

    print()
    print("verdict distribution:", dict(Counter(r.verdict.value for _, r, *_ in rows)))
    return 0 if directional == n else 1


if __name__ == "__main__":
    sys.exit(run())
