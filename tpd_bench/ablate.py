"""Ablation diagnostic.

The calibration score is partly fitted -- labels were written while watching the
engine. This is not. It asks a question labels cannot answer: which terms are
actually doing work?

A term whose removal changes no verdict is dead weight dressed up as rigour.
A term whose removal changes most verdicts is load-bearing, and the whole bench
rests on whether that one judgement is right.
"""

from __future__ import annotations

import copy
import dataclasses

from .schema import (
    ClinicalOutcome,
    CompetitorStage,
    GeneticsClass,
    InhibitorStatus,
    MRStatus,
    ScaffoldSign,
)
from .scoring import score
from .targets import TARGETS

# Each ablation neutralises one feature across every target.
ABLATIONS: dict[str, dict] = {
    "human genetics": {"genetics": GeneticsClass.NONE},
    "MR status": {"mr_status": MRStatus.NONE},
    "clinical outcome": {"clinical_outcome": ClinicalOutcome.UNKNOWN},
    "approved drug on pathway": {"approved_drug_on_pathway": False},
    "disease expression": {"disease_expression_fold": 1.0},
    "network redundancy": {"network_redundancy": 0.0},
    "scaffold evidence": {"scaffold_evidence_tier": 0, "scaffold_sign": ScaffoldSign.NONE},
    "kinase-dead phenocopy": {"kinase_dead_phenocopies_ko": False},
    "inhibitor status": {"inhibitor_status": InhibitorStatus.NONE_TRIED},
    "compensatory upregulation": {"upregulated_under_inhibition": False},
    "paralog selectivity": {"paralog_selectivity_blocked": False},
    "degron features": {"degron_feature_score": 0.0},
    "protein turnover": {"protein_half_life_h": None},
    "catalytic advantage": {"catalytic_advantage": False},
    "competitor landscape": {
        "competitor_stage": CompetitorStage.NONE,
        "competitor_dosing_interval_days": None,
    },
}


def run() -> int:
    baseline = {t.target_id: score(t).verdict for t in TARGETS}
    n = len(TARGETS)

    print("=" * 74)
    print("ABLATION -- how many verdicts move when each term is removed")
    print("=" * 74)

    results = []
    for name, patch in ABLATIONS.items():
        moved = []
        for t in TARGETS:
            t2 = dataclasses.replace(copy.deepcopy(t), **patch)
            if score(t2).verdict is not baseline[t.target_id]:
                moved.append(t.gene)
        results.append((len(moved), name, moved))

    results.sort(reverse=True)
    for count, name, moved in results:
        bar = "#" * count
        note = ""
        if count == 0:
            note = "   <- changes nothing; either dead weight or untested by this set"
        elif count >= n * 0.3:
            note = "   <- load-bearing; the bench rests on this judgement"
        print(f"  {count:3d}/{n}  {name:28s} {bar}{note}")
        if 0 < count <= 6:
            print(f"          {', '.join(moved)}")

    dead = [name for c, name, _ in results if c == 0]
    print()
    if dead:
        print("Terms with no effect on any verdict:")
        for d in dead:
            print(f"  - {d}")
        print(
            "\nThese are not necessarily wrong. They may matter on targets this set\n"
            "does not contain -- which is itself a gap in the set, not a defence of\n"
            "the term. Either find cases that exercise them, or cut them."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
