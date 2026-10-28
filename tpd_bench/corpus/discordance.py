"""Scaffolding discordance corpus.

Two pipelines feeding one field, `scaffold_evidence_tier`, plus the thing that
actually matters: `scaffold_sign`.

  Tier 1  human resistance mutations -- somatic variants that ablate catalytic
          function while signalling persists. Small, decisive, and the reason
          BTK is the strongest case in the set.
  Tier 2/3 knockout-vs-inhibitor phenotype join -- broad, systematic, and noisy
          enough that most raw hits are artefacts rather than scaffolds.

THE FALSE-DISCORDANCE PROBLEM
A naive join asks "does the knockout phenotype differ from the inhibitor
phenotype?" and finds a great many hits, nearly all of them wrong. Four
confounders produce a difference with no scaffold involved:

  1. The inhibitor never engaged the target in vivo. No phenotype because no
     pharmacology, not because of a scaffold.
  2. The inhibitor hit something else as well. A phenotype the knockout lacks
     is an off-target effect, which is discordance pointing the other way.
  3. Constitutive knockout removes the protein through development. The adult
     phenotype may be developmental compensation, not adult target biology.
  4. Timescale. Chronic genetic absence versus acute pharmacological block are
     different experiments even when both are clean.

Every candidate is checked against all four, and anything failing a check is
demoted rather than dropped -- the demotion is informative.

AND THE SIGN
Discordance alone is not a degradation argument. RIPK1 has textbook tier-2
discordance and the scaffold is protective: removing the protein is worse than
inhibiting it. The pipeline therefore classifies direction, and a PROTECTIVE
sign is a rejection, not a nomination.
"""

from __future__ import annotations

import csv
import pathlib
from dataclasses import asdict, dataclass, field
from enum import Enum


class Confounder(str, Enum):
    NO_TARGET_ENGAGEMENT = "no_target_engagement"
    INHIBITOR_OFF_TARGET = "inhibitor_off_target"
    DEVELOPMENTAL_COMPENSATION = "developmental_compensation"
    TIMESCALE_MISMATCH = "timescale_mismatch"
    UNMATCHED_MODEL = "unmatched_model"


class Sign(str, Enum):
    PATHOGENIC = "pathogenic"      # removing the protein helps beyond inhibition
    PROTECTIVE = "protective"      # removing the protein is worse -- reject
    UNKNOWN = "unknown"


@dataclass
class ResistanceEvidence:
    """Tier 1. A human somatic variant that keeps the protein, loses the enzyme."""

    gene: str
    variant: str
    catalytic_impact: str          # ablates | impairs | none
    signalling_persists: bool
    source: str

    def qualifies(self) -> bool:
        return self.catalytic_impact in ("ablates", "impairs") and self.signalling_persists


@dataclass
class PhenotypeComparison:
    """Tier 2/3. One knockout-vs-inhibitor comparison, with its controls."""

    gene: str
    model: str
    ko_phenotype: str
    inhibitor_phenotype: str
    ko_worse_than_inhibitor: bool
    # controls
    target_engagement_shown: bool
    inhibitor_selectivity_known: bool
    conditional_ko: bool
    matched_model: bool
    source: str
    # Whether the two arms reached the SAME outcome. Phenotype strings differ
    # in wording all the time; only an annotator can say if they differ in
    # substance. Defaulting this to "they differ" invents scaffolds.
    outcomes_equivalent: bool = False

    def confounders(self) -> list[Confounder]:
        out = []
        if not self.target_engagement_shown:
            out.append(Confounder.NO_TARGET_ENGAGEMENT)
        if not self.inhibitor_selectivity_known:
            out.append(Confounder.INHIBITOR_OFF_TARGET)
        if not self.conditional_ko:
            out.append(Confounder.DEVELOPMENTAL_COMPENSATION)
        if not self.matched_model:
            out.append(Confounder.UNMATCHED_MODEL)
        return out

    def sign(self) -> Sign:
        if self.outcomes_equivalent or self.ko_phenotype == self.inhibitor_phenotype:
            return Sign.UNKNOWN
        return Sign.PROTECTIVE if self.ko_worse_than_inhibitor else Sign.PATHOGENIC

    def tier(self) -> int:
        """Tier 2 only when every control passes. Each failed control demotes."""
        c = len(self.confounders())
        if c == 0:
            return 2
        if c <= 2:
            return 3
        return 4


@dataclass
class CorpusEntry:
    gene: str
    tier: int
    sign: str
    basis: str
    confounders: list[str] = field(default_factory=list)
    source: str = ""


def build(
    resistance: list[ResistanceEvidence],
    comparisons: list[PhenotypeComparison],
) -> list[CorpusEntry]:
    """Tier 1 wins where it exists; a comparison can only confirm or downgrade."""
    out: dict[str, CorpusEntry] = {}

    for c in comparisons:
        e = CorpusEntry(
            gene=c.gene,
            tier=c.tier(),
            sign=c.sign().value,
            basis=f"{c.ko_phenotype} (KO) vs {c.inhibitor_phenotype} (inhibitor), {c.model}",
            confounders=[x.value for x in c.confounders()],
            source=c.source,
        )
        prev = out.get(c.gene)
        if prev is None or e.tier < prev.tier:
            out[c.gene] = e

    by_gene: dict[str, list[ResistanceEvidence]] = {}
    for r in resistance:
        if r.qualifies():
            by_gene.setdefault(r.gene, []).append(r)

    for gene, variants in by_gene.items():
        r = variants[0]
        listed = ", ".join(f"{v.variant} ({v.catalytic_impact})" for v in variants)
        prev = out.get(gene)
        # A protective sign from a clean comparison is not overridden by a
        # resistance variant: the direction matters more than the tier.
        if prev and prev.sign == Sign.PROTECTIVE.value and prev.tier <= 2:
            prev.basis += f"; tier-1 variants {listed} noted but sign stands"
            continue
        out[gene] = CorpusEntry(
            gene=gene,
            tier=1,
            sign=Sign.PATHOGENIC.value,
            basis=f"{len(variants)} variant(s): {listed}; signalling persists",
            source=r.source,
        )

    return sorted(out.values(), key=lambda e: (e.tier, e.gene))


def write_csv(entries: list[CorpusEntry], path: pathlib.Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(entries[0])))
        w.writeheader()
        for e in entries:
            row = asdict(e)
            row["confounders"] = "|".join(row["confounders"])
            w.writerow(row)
