"""Feature schema for the TPD glue bench.

Design rule: everything in Target is an OBJECTIVE, MINEABLE feature — something
you could in principle populate from a database, a paper, or a structure, without
knowing whether the target is a good idea. Ground-truth verdicts live in
calibration_labels.csv and are never read by the scorer.

If you find yourself wanting to add a field like `is_good_target`, stop: that is
a label, not a feature.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


# ---------------------------------------------------------------- enumerations


class Subcellular(str, Enum):
    """Where the protein sits. Gates intracellular degradation absolutely."""

    NUCLEAR = "nuclear"
    CYTOSOLIC = "cytosolic"
    MEMBRANE = "membrane"
    SECRETED = "secreted"


class GeneticsClass(str, Enum):
    """Strength of human genetic evidence, ordered."""

    NONE = "none"
    GWAS_ASSOCIATION = "gwas_association"
    RARE_VARIANT_ASSOCIATION = "rare_variant_association"
    LOF_PROTECTIVE = "lof_protective"
    LOF_PROTECTIVE_REPLICATED = "lof_protective_replicated"


class MRStatus(str, Enum):
    """Mendelian randomisation status.

    MATCHED_ENGAGEMENT is the distinction that separates Lp-PLA2 (where the
    genetic estimate was standardised to the reduction the drug actually
    achieves, and predicted the trial) from IL6R and LPA (where it was not).
    """

    NONE = "none"
    SUPPORTED = "supported"
    SUPPORTED_MATCHED_ENGAGEMENT = "supported_matched_engagement"
    REFUTED = "refuted"


class ClinicalOutcome(str, Enum):
    """Hardest clinical evidence available for this target-indication pair."""

    UNKNOWN = "unknown"
    PHASE1_BIOMARKER = "phase1_biomarker"
    SURROGATE_APPROVAL = "surrogate_approval"
    HARD_ENDPOINT_POSITIVE = "hard_endpoint_positive"
    FAILED_FUTILITY = "failed_futility"
    FAILED_HARM = "failed_harm"
    DISCONTINUED_STRATEGIC = "discontinued_strategic"


class InhibitorStatus(str, Enum):
    """What existing small-molecule pharmacology does to this target."""

    NONE_TRIED = "none_tried"
    NO_TRACTABLE_INHIBITOR = "no_tractable_inhibitor"
    FAILED_EFFICACY = "failed_efficacy"
    FAILED_TOXICITY = "failed_toxicity"
    WORKS_TOLERATED = "works_tolerated"
    WORKS_NARROW_WINDOW = "works_narrow_window"


class ScaffoldSign(str, Enum):
    """Direction of the non-catalytic contribution.

    The RIPK1 lesson: evidence that knockout differs from inhibition is not
    enough. PROTECTIVE means removing the protein is worse than inhibiting it.
    """

    NONE = "none"
    UNKNOWN = "unknown"
    PATHOGENIC = "pathogenic"
    PROTECTIVE = "protective"


class Ligandability(str, Enum):
    """Bifunctional tractability cascade, in the order specified."""

    KNOWN_LIGAND = "known_ligand"
    REPORTED_POCKET = "reported_pocket"
    PREDICTED_POCKET = "predicted_pocket"
    NONE = "none"


class Track(str, Enum):
    """Which degradation chemistry a target routes to.

    Set by compartment, not chosen. Nuclear and cytosolic proteins go to the
    ubiquitin-proteasome track (glues, PROTACs); membrane and secreted proteins
    go to the lysosomal track (LYTAC-style conjugates). The two are judged
    against different delivery bars, so they are ranked as separate lists.
    """

    INTRACELLULAR = "intracellular"
    EXTRACELLULAR = "extracellular"


class CompetitorStage(str, Enum):
    """Most advanced competing modality against this target-indication pair.

    Most validated secreted targets are already owned by siRNA or antibodies.
    An extracellular nomination has to beat what exists, not just work.
    """

    NONE = "none"
    PRECLINICAL = "preclinical"
    CLINICAL = "clinical"
    PHASE3 = "phase3"
    APPROVED = "approved"


class TissueOfAction(str, Enum):
    """Determines which delivery gate the candidate is judged against."""

    SYSTEMIC = "systemic"
    HEPATIC = "hepatic"
    LOCAL = "local"


# scaffold evidence tiers -> confidence multiplier
SCAFFOLD_TIER_WEIGHT = {
    1: 1.00,  # human resistance mutations, or clinical failure + genetic validation
    2: 0.75,  # inhibitor vs knockout, in vivo
    3: 0.50,  # inhibitor vs knockdown, or degrader vs inhibitor only
    4: 0.30,  # structural PPI evidence, predictive
    0: 0.00,  # none
}


@dataclass
class Target:
    """One target-indication pair. Validity is indication-specific by design."""

    target_id: str
    gene: str
    indication: str
    disease_area: str  # immunology_inflammation | cardiometabolic

    # --- where it is -------------------------------------------------------
    subcellular: Subcellular
    tissue_of_action: TissueOfAction = TissueOfAction.SYSTEMIC

    # --- axis 1 inputs: is the target real? --------------------------------
    genetics: GeneticsClass = GeneticsClass.NONE
    genetics_direction_supports_removal: bool = False
    mr_status: MRStatus = MRStatus.NONE
    approved_drug_on_pathway: bool = False
    clinical_outcome: ClinicalOutcome = ClinicalOutcome.UNKNOWN
    disease_expression_fold: float = 1.0  # diseased vs matched normal
    pathway_node_distance: int = 0  # 0 = the validated node itself
    network_redundancy: float = 0.0  # 0 = unique, 1 = highly degenerate class

    # --- axis 2 inputs: is degradation the right tool? ---------------------
    scaffold_evidence_tier: int = 0
    scaffold_sign: ScaffoldSign = ScaffoldSign.NONE
    kinase_dead_phenocopies_ko: bool = False
    inhibitor_status: InhibitorStatus = InhibitorStatus.NONE_TRIED
    upregulated_under_inhibition: bool = False
    paralog_selectivity_blocked: bool = False
    ligandability: Ligandability = Ligandability.NONE
    degron_feature_score: float = 0.0  # 0-1, glue amenability
    protein_half_life_h: Optional[float] = None
    # gnomAD LOEUF: upper bound of the observed/expected LoF ratio. Below ~0.35
    # the gene is strongly LoF-intolerant, meaning humans carrying one broken
    # copy are already affected -- a direct argument against removing the
    # protein pharmacologically. This is NOT direction of effect; it is a
    # safety gate that the bench previously had nothing to fill.
    lof_constraint_loeuf: Optional[float] = None

    # --- extracellular track inputs ----------------------------------------
    # Sub-stoichiometric clearance vs occupancy-driven binding. The core
    # mechanistic argument for a degrader over an antibody or siRNA.
    catalytic_advantage: bool = False
    # Binding fails but removal works: aggregates, lipoprotein particles,
    # complexes that neutralisation cannot address.
    clears_non_neutralisable: bool = False
    competitor_stage: CompetitorStage = CompetitorStage.NONE
    competitor_dosing_interval_days: Optional[int] = None

    # --- E3 opportunity (reported, never scored) ---------------------------
    # Restricted-expression ligase that would suit this target's cell type.
    # Emitted as a ligand-discovery opportunity; carries no scoring weight,
    # because oral delivery pins the intracellular track to CRBN today.
    e3_restricted_candidate: Optional[str] = None

    # Fields this record never looked up. Scoring renormalises these out
    # rather than treating them as zero. Curated records leave it empty; the
    # sweep fills it with whatever layer it could not populate. Absent
    # evidence and negative evidence are different things, and conflating
    # them is what crushed the first three sweeps.
    unknown_fields: frozenset = frozenset()

    # --- provenance --------------------------------------------------------
    notes: str = ""
    sources: list[str] = field(default_factory=list)

    @property
    def is_intracellular(self) -> bool:
        return self.subcellular in (Subcellular.NUCLEAR, Subcellular.CYTOSOLIC)

    @property
    def track(self) -> Track:
        return Track.INTRACELLULAR if self.is_intracellular else Track.EXTRACELLULAR
