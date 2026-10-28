"""Structure-based degron feature scoring.

The glue-first path needs a glue-amenability number for every protein in the
universe, and no public resource provides one. This computes it from structure.

WHAT IS BEING SCORED
A CRBN molecular glue works by inducing a protein-protein interface, so the
question is whether a surface-exposed, structurally ordered loop can sit in the
CRBN tri-tryptophan pocket. The classical answer is a beta-hairpin degron with
a key glycine at the turn -- the "G-loop" rule, from the C2H2 zinc fingers.

WHY THIS IS NOT JUST THE G-LOOP RULE
The G-loop rule is actively being broken. CRBN glues reach VAV1 through an SH3
RT-loop degron with no G-loop, and work mining the CRBN target space has
revised neosubstrate recognition more broadly. A binary G-loop gate would
systematically miss the targets the field is currently discovering.

So this scores a graded set of features, of which the G-loop is the
highest-weighted but not a requirement. `RULES_VERSION` exists because this
will need revisiting -- treat a stale version as a reason to recompute, not as
settled fact.

INPUT
Per-residue data from an AlphaFold model: sequence, secondary structure,
relative solvent accessibility, and pLDDT. Low pLDDT means the model does not
know where the loop is, so confidence propagates into the score rather than
being silently ignored.
"""

from __future__ import annotations

from dataclasses import dataclass, field

RULES_VERSION = "2026-10-28"

# Minimum model confidence for a residue to count toward a degron claim.
# Below this, AlphaFold is not asserting a conformation.
PLDDT_FLOOR = 70.0

# Relative solvent accessibility above which a residue counts as exposed.
RSA_EXPOSED = 0.25


@dataclass
class Residue:
    index: int
    aa: str
    ss: str          # H helix | E strand | C coil/turn
    rsa: float       # 0-1 relative solvent accessibility
    plddt: float


@dataclass
class DegronHit:
    start: int
    end: int
    motif: str
    score: float
    features: list[str] = field(default_factory=list)


def _turn_regions(res: list[Residue]) -> list[tuple[int, int]]:
    """Coil stretches flanked by strands: candidate beta-hairpin turns."""
    out = []
    i = 0
    while i < len(res):
        if res[i].ss == "C":
            j = i
            while j < len(res) and res[j].ss == "C":
                j += 1
            # a hairpin turn is short, and bracketed by strand on both sides
            if 2 <= (j - i) <= 6 and i > 0 and j < len(res):
                if res[i - 1].ss == "E" and res[j].ss == "E":
                    out.append((i, j - 1))
            i = j
        else:
            i += 1
    return out


def _exposed(r: Residue) -> bool:
    return r.rsa >= RSA_EXPOSED and r.plddt >= PLDDT_FLOOR


def score_protein(res: list[Residue]) -> tuple[float, list[DegronHit]]:
    """Return a 0-1 glue-amenability score and the hits behind it.

    The score is the best single hit, not a sum: a protein with one good degron
    is glue-amenable. Ten mediocre loops are not a substitute for one real one,
    and summing would reward large proteins for being large.
    """
    hits: list[DegronHit] = []

    for start, end in _turn_regions(res):
        turn = res[start : end + 1]
        if not any(_exposed(r) for r in turn):
            continue

        feats: list[str] = []
        s = 0.35  # base: an exposed, ordered beta-hairpin turn

        seq = "".join(r.aa for r in turn)

        # Classical G-loop: glycine at the turn. Highest-weighted single
        # feature, because it is the best-evidenced one.
        if "G" in seq:
            s += 0.30
            feats.append("G-loop glycine at turn")

        # A second glycine gives the backbone the flexibility the pocket wants.
        if seq.count("G") >= 2:
            s += 0.05
            feats.append("second glycine")

        # Aromatic adjacent to the turn stacks against the tri-Trp pocket.
        flank = res[max(0, start - 3) : start] + res[end + 1 : end + 4]
        if any(r.aa in "FWY" for r in flank if _exposed(r)):
            s += 0.12
            feats.append("aromatic flanking the turn")

        # A surface lysine within reach is needed for ubiquitin transfer.
        near = res[max(0, start - 12) : min(len(res), end + 13)]
        if any(r.aa == "K" and _exposed(r) for r in near):
            s += 0.10
            feats.append("accessible lysine within transfer range")
        else:
            s -= 0.15
            feats.append("no accessible lysine nearby -- transfer may fail")

        # Model confidence scales the whole claim rather than gating it.
        conf = min(1.0, sum(r.plddt for r in turn) / (len(turn) * 90.0))
        s *= conf
        if conf < 0.85:
            feats.append(f"low model confidence (x{conf:.2f})")

        hits.append(
            DegronHit(
                start=res[start].index,
                end=res[end].index,
                motif=seq,
                score=round(max(0.0, min(1.0, s)), 3),
                features=feats,
            )
        )

    # Non-canonical route: an exposed, ordered loop in an SH3-like context,
    # which is how CRBN glues reach VAV1 without a G-loop. Scored lower than a
    # true hairpin because the evidence base is thinner, but not excluded.
    for i in range(len(res) - 8):
        w = res[i : i + 9]
        if all(r.ss == "C" for r in w[2:7]) and sum(_exposed(r) for r in w) >= 5:
            conf = min(1.0, sum(r.plddt for r in w) / (len(w) * 90.0))
            s = 0.40 * conf
            if any(r.aa == "K" and _exposed(r) for r in w):
                s += 0.08
            hits.append(
                DegronHit(
                    start=w[0].index,
                    end=w[-1].index,
                    motif="".join(r.aa for r in w),
                    score=round(max(0.0, min(1.0, s)), 3),
                    features=["non-canonical exposed loop (RT-loop-like)"],
                )
            )

    if not hits:
        return 0.0, []
    hits.sort(key=lambda h: -h.score)
    return hits[0].score, hits[:5]
