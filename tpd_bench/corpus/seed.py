"""Verified seed data for the discordance corpus.

Every entry here was checked against a primary or near-primary source during
the 2026-10-28 research pass. This is the offline cache: the sandbox cannot
reach cBioPortal, IMPC or UniProt, so live mining runs in your environment via
`sources.py`. Nothing in this file is inferred or filled in.

Where a control could not be confirmed from the source, it is recorded as
False. That demotes the tier, which is the correct outcome -- an unverified
control is not a passed control.
"""

from __future__ import annotations

from .discordance import PhenotypeComparison, ResistanceEvidence

RESISTANCE: list[ResistanceEvidence] = [
    ResistanceEvidence(
        gene="BTK", variant="L528W", catalytic_impact="ablates",
        signalling_persists=True,
        source="Montoya et al. Science 2024; Nurix ICML 2023 poster",
    ),
    ResistanceEvidence(
        gene="BTK", variant="T474I", catalytic_impact="impairs",
        signalling_persists=True,
        source="Montoya et al. Science 2024",
    ),
    ResistanceEvidence(
        gene="BTK", variant="V416L", catalytic_impact="impairs",
        signalling_persists=True,
        source="Montoya et al. Science 2024",
    ),
]

COMPARISONS: list[PhenotypeComparison] = [
    # The case that defines the sign problem.
    PhenotypeComparison(
        gene="RIPK1", model="mouse, kinase-dead knock-in vs full knockout",
        ko_phenotype="postnatal lethality, dysregulated cell death, systemic inflammation",
        inhibitor_phenotype="viable, phenotypically normal",
        ko_worse_than_inhibitor=True,
        target_engagement_shown=True, inhibitor_selectivity_known=True,
        conditional_ko=True, matched_model=True,
        source="Front Immunol 2021;12:676946; Mannion et al. Immunity 2024",
    ),
    # The negative control: no discordance at all.
    PhenotypeComparison(
        gene="HPK1", model="Jurkat and primary human T cells",
        ko_phenotype="loss of TCR suppression",
        inhibitor_phenotype="loss of TCR suppression",
        ko_worse_than_inhibitor=False, outcomes_equivalent=True,
        target_engagement_shown=True, inhibitor_selectivity_known=True,
        conditional_ko=False, matched_model=True,
        source="Wang et al. Front Immunol 2025;16:1449106",
    ),
    PhenotypeComparison(
        gene="IRAK4", model="mouse lupus background, kinase-dead knock-in",
        ko_phenotype="reduced splenomegaly, proteinuria, ANA, glomerulonephritis",
        inhibitor_phenotype="reduced disease (kinase inhibition alone sufficient in this model)",
        ko_worse_than_inhibitor=False, outcomes_equivalent=True,
        target_engagement_shown=True, inhibitor_selectivity_known=False,
        conditional_ko=False, matched_model=True,
        source="Murphy et al. Eur J Immunol 2017;47(5):880-891",
    ),
    # Degrader-vs-inhibitor only: no genetic arm, so controls fail by design.
    PhenotypeComparison(
        gene="PTK2", model="A549 and related lines, degrader vs inhibitor",
        ko_phenotype="wound healing -53pct, invasion -65pct at 50 nM (degrader)",
        inhibitor_phenotype="no effect at matched concentration (defactinib)",
        ko_worse_than_inhibitor=False,
        target_engagement_shown=True, inhibitor_selectivity_known=False,
        conditional_ko=False, matched_model=False,
        source="Pogasha & Fletcher RSC Med Chem 2025;16:2958-2973; Cromm JACS 2018",
    ),
    PhenotypeComparison(
        gene="AURKA", model="neuroblastoma lines, degrader vs inhibitor",
        ko_phenotype="S-phase arrest; TPX2-associated spindle pool removed (JB170)",
        inhibitor_phenotype="G2/M arrest; both pools inhibited (alisertib)",
        ko_worse_than_inhibitor=False,
        target_engagement_shown=True, inhibitor_selectivity_known=True,
        conditional_ko=False, matched_model=True,
        source="Pogasha & Fletcher RSC Med Chem 2025",
    ),
]
