"""Two-axis scoring for the TPD glue bench.

The two axes are never collapsed into a single number without both being shown.
A target can be real and undegradable, or degradable and not real, and those are
different failure modes that need different responses.

Hard gates produce VERDICTS, not just low scores. "This is not a degradation
problem" has to be a first-class output or the bench will recommend a degrader
for every target it likes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from .schema import (
    SCAFFOLD_TIER_WEIGHT,
    CompetitorStage,
    Track,
    ClinicalOutcome,
    GeneticsClass,
    InhibitorStatus,
    Ligandability,
    MRStatus,
    ScaffoldSign,
    Target,
    TissueOfAction,
)


I_WORKS_TOLERATED = InhibitorStatus.WORKS_TOLERATED


class Verdict(str, Enum):
    TOP = "top"
    MID = "mid"
    REJECT_TARGET_INVALID = "reject_target_invalid"
    REJECT_WRONG_MODALITY = "reject_wrong_modality"
    # Right compartment, degradable in principle, but the modality case is
    # simply too thin. Distinct from wrong_modality, which means unreachable.
    REJECT_MODALITY_TOO_WEAK = "reject_modality_too_weak"
    REJECT_INHIBITION_SUFFICIENT = "reject_inhibition_sufficient"
    REJECT_SCAFFOLD_PROTECTIVE = "reject_scaffold_protective"
    # Strong modality case, validity below the commit bar but not refuted.
    # A hypothesis engine that cannot emit this can only re-nominate targets
    # someone has already validated -- which is not a hypothesis engine.
    NOMINATE_FOR_VALIDATION = "nominate_for_validation"
    FLAG_CONTESTED = "flag_contested"
    FLAG_STRATEGIC_DISCONTINUATION = "flag_strategic_discontinuation"


@dataclass
class Profile:
    """Disease-area weight profile. Selected by the toggle before a run."""

    name: str
    w_genetics: float
    w_pathway_drug: float
    w_clinical: float
    w_expression: float
    # MR supported but NOT matched to achievable engagement is discounted.
    # ZEUS is why this is not 1.0.
    mr_unmatched_discount: float = 0.60
    top_threshold: float = 55.0
    mid_threshold: float = 35.0
    strong_modality: float = 55.0
    # A weak modality axis cannot be rescued by a strong validity axis, so it
    # must not reach MID through the composite either.
    modality_floor: float = 20.0
    validation_floor: float = 12.0


PROFILES: dict[str, Profile] = {
    # Human genetics is a near-causal oracle here, so it carries more weight --
    # but only when the instrument is matched to achievable target engagement.
    "cardiometabolic": Profile(
        name="cardiometabolic",
        w_genetics=0.45,
        w_pathway_drug=0.10,
        w_clinical=0.35,
        w_expression=0.10,
    ),
    # Genetics gives loci, not mechanisms. An approved biologic on the axis is
    # usually the stronger evidence, so it is weighted up and genetics down.
    "immunology_inflammation": Profile(
        name="immunology_inflammation",
        w_genetics=0.25,
        w_pathway_drug=0.30,
        w_clinical=0.35,
        w_expression=0.10,
    ),
    # Mechanism controls that exist to grade the modality axis. Balanced
    # weights so the validity axis neither rescues nor sinks them.
    "cross_cutting": Profile(
        name="cross_cutting",
        w_genetics=0.30,
        w_pathway_drug=0.25,
        w_clinical=0.35,
        w_expression=0.10,
    ),
}


GENETICS_POINTS = {
    GeneticsClass.NONE: 0.0,
    GeneticsClass.GWAS_ASSOCIATION: 25.0,
    GeneticsClass.RARE_VARIANT_ASSOCIATION: 50.0,
    GeneticsClass.LOF_PROTECTIVE: 80.0,
    GeneticsClass.LOF_PROTECTIVE_REPLICATED: 100.0,
}

CLINICAL_POINTS = {
    # UNKNOWN is handled by renormalisation, not by a score. Scoring "nobody
    # has tested this" as weak evidence makes the bench structurally unable to
    # nominate a novel target, which is the one thing it exists to do.
    ClinicalOutcome.UNKNOWN: None,
    ClinicalOutcome.PHASE1_BIOMARKER: 50.0,
    ClinicalOutcome.SURROGATE_APPROVAL: 70.0,
    ClinicalOutcome.HARD_ENDPOINT_POSITIVE: 100.0,
    ClinicalOutcome.FAILED_FUTILITY: 10.0,
    ClinicalOutcome.FAILED_HARM: 0.0,
    # A program dropped for competitive or portfolio reasons says nothing about
    # the target. Scoring it as a failure teaches that good targets are bad.
    ClinicalOutcome.DISCONTINUED_STRATEGIC: 55.0,
}


@dataclass
class AxisResult:
    score: float
    flags: list[str] = field(default_factory=list)
    # Stable machine-readable reasons. Verdict logic reads these, never the
    # human-readable flag strings -- matching on prose is how you get a bench
    # that breaks when someone fixes a typo.
    reasons: set[str] = field(default_factory=set)
    # Per-term contributions, so a slider can reweight the score without
    # re-running the bench. None means the term was never looked up and is
    # renormalised out rather than counted as zero.
    terms: dict = field(default_factory=dict)


@dataclass
class Result:
    target_id: str
    gene: str
    indication: str
    validity: AxisResult
    modality: AxisResult
    verdict: Verdict
    track: Track
    composite: float
    flags: list[str] = field(default_factory=list)


# --------------------------------------------------------------- axis 1


def score_validity(t: Target, p: Profile) -> AxisResult:
    flags: list[str] = []
    reasons: set[str] = set()

    # MR refutation is the most reliable rejection signal in the whole set.
    # HDL-C, Lp-PLA2, sPLA2 and CRP all have strong observational signal and
    # null causal signal, and all four produced null or harmful trials.
    if t.mr_status is MRStatus.REFUTED:
        return AxisResult(
            0.0,
            ["mr_refuted: causal signal null despite observational association"],
            {"validity_refuted"},
        )

    genetics = GENETICS_POINTS[t.genetics]
    if genetics > 0 and not t.genetics_direction_supports_removal:
        if "direction" in t.unknown_fields:
            # Unknown is not contrary. Flag it; do not discount for it.
            flags.append("direction of effect unknown -- not established either way")
        else:
            genetics *= 0.4
            flags.append("genetic direction does not support target removal")

    if t.mr_status is MRStatus.SUPPORTED:
        genetics *= p.mr_unmatched_discount
        flags.append(
            "MR supported but not standardised to achievable engagement -- discounted (see ZEUS)"
        )
    elif t.mr_status is MRStatus.SUPPORTED_MATCHED_ENGAGEMENT:
        flags.append("MR standardised to achievable engagement -- full weight")

    if t.genetics is GeneticsClass.NONE:
        flags.append("no human genetic support -- flag, not veto")

    pathway = 100.0 if t.approved_drug_on_pathway else 0.0
    if pathway and t.pathway_node_distance > 0:
        # IL-17 sits one node from IL-23 on the same axis and carries the
        # opposite sign. Distance from the validated node is a real discount.
        pathway *= max(0.0, 1.0 - 0.45 * t.pathway_node_distance)
        flags.append(
            f"target is {t.pathway_node_distance} node(s) from the validated one -- discounted"
        )

    expression = min(100.0, max(0.0, (t.disease_expression_fold - 1.0) * 40.0))

    # Absent evidence is renormalised out; it is not evidence against.
    terms = []
    if "genetics" not in t.unknown_fields:
        terms.append((p.w_genetics, genetics))
    if "pathway" not in t.unknown_fields:
        terms.append((p.w_pathway_drug, pathway))
    else:
        flags.append("pathway drug status not looked up -- renormalised out")
    if "expression" not in t.unknown_fields:
        terms.append((p.w_expression, expression))
    else:
        flags.append("disease expression not looked up -- renormalised out")
    clinical = CLINICAL_POINTS[t.clinical_outcome]
    if clinical is None:
        flags.append("no clinical data yet -- renormalised out, not penalised")
    else:
        terms.append((p.w_clinical, clinical))

    breakdown = {
        "genetics": None if "genetics" in t.unknown_fields else genetics,
        "pathway": None if "pathway" in t.unknown_fields else pathway,
        "expression": None if "expression" in t.unknown_fields else expression,
        "clinical": clinical,
    }
    total_w = sum(w for w, _ in terms)
    score = sum(w * v for w, v in terms) / total_w if total_w else 0.0
    # Renormalising is right, but a score resting on one layer is thin even
    # when it is high. Say so rather than letting 100.0 imply four agreeing
    # sources when it means one.
    if len(terms) <= 1:
        flags.append(
            f"score rests on {len(terms)} evidence layer -- treat as provisional"
        )

    if t.network_redundancy > 0:
        score *= 1.0 - 0.5 * t.network_redundancy
        if t.network_redundancy >= 0.5:
            flags.append("highly redundant target class -- structural penalty")

    if t.clinical_outcome is ClinicalOutcome.FAILED_HARM:
        flags.append("trial showed harm -- validity refuted regardless of other evidence")
        score = 0.0
        reasons.add("validity_refuted")
    if t.clinical_outcome is ClinicalOutcome.SURROGATE_APPROVAL:
        flags.append("surrogate approval only -- outcome unconfirmed")
    if t.clinical_outcome is ClinicalOutcome.DISCONTINUED_STRATEGIC:
        flags.append("strategic discontinuation -- not scored as scientific failure")

    return AxisResult(round(min(100.0, score), 1), flags, reasons, breakdown)


# --------------------------------------------------------------- axis 2


def score_modality(t: Target) -> AxisResult:
    flags: list[str] = []

    # Compartment routes rather than rejects. Membrane and secreted proteins
    # go to the lysosomal track, judged against its own bar.
    if not t.is_intracellular:
        return _score_extracellular(t)

    # Hard gate: RIPK1. Scaffold evidence exists and points the wrong way.
    if t.scaffold_sign is ScaffoldSign.PROTECTIVE:
        return AxisResult(
            0.0,
            ["scaffold is protective -- removing the protein is harmful (see RIPK1)"],
            {"scaffold_protective"},
        )

    # Hard gate: HPK1. Kinase-dead phenocopies knockout, so there is nothing
    # for degradation to add.
    if t.kinase_dead_phenocopies_ko:
        return AxisResult(
            10.0,
            ["kinase-dead phenocopies knockout -- inhibition is sufficient"],
            {"inhibition_sufficient"},
        )

    # Hard gate: a tolerated inhibitor that already works, with no scaffold
    # evidence to argue otherwise. TYK2 is the case.
    if t.inhibitor_status is I_WORKS_TOLERATED and t.scaffold_evidence_tier == 0:
        return AxisResult(
            5.0,
            ["a tolerated inhibitor already works and no scaffold evidence exists"],
            {"inhibition_sufficient"},
        )

    score = 0.0
    reasons: set[str] = set()

    # --- entry criterion 1: scaffolding ---------------------------------
    tier_w = SCAFFOLD_TIER_WEIGHT.get(t.scaffold_evidence_tier, 0.0)
    if tier_w and t.scaffold_sign is ScaffoldSign.PATHOGENIC:
        score += 40.0 * tier_w
        flags.append(f"scaffold evidence tier {t.scaffold_evidence_tier} (x{tier_w:.2f})")
    elif tier_w and t.scaffold_sign is ScaffoldSign.UNKNOWN:
        score += 20.0 * tier_w
        flags.append("scaffold sign unknown -- half credit pending direction")

    # --- entry criterion 2: compensatory upregulation --------------------
    if t.upregulated_under_inhibition:
        score += 18.0
        flags.append("target accumulates under inhibition -- degradation addresses this")

    # --- entry criterion 3: ligand-silent --------------------------------
    if t.inhibitor_status is InhibitorStatus.NO_TRACTABLE_INHIBITOR:
        score += 24.0
        flags.append("no tractable inhibitor -- degradation is the only route")
        if t.ligandability is Ligandability.NONE:
            # Nothing to inhibit and no pocket to find: the purest glue case,
            # and exactly the TF-shaped target this bench exists to surface.
            score += 10.0
            flags.append("no pocket and no inhibitor -- pure glue play, no competing modality")
    elif t.inhibitor_status is InhibitorStatus.FAILED_EFFICACY:
        score += 8.0
        flags.append("inhibitor engaged target without efficacy -- weak degradation argument")
    elif t.inhibitor_status is InhibitorStatus.WORKS_TOLERATED:
        score -= 25.0
        flags.append("a tolerated inhibitor already works -- degradation adds little")
    elif t.inhibitor_status is InhibitorStatus.WORKS_NARROW_WINDOW:
        score += 10.0
        flags.append("inhibitor works but window is narrow -- selectivity opportunity")

    # --- entry criterion 4: selectivity-limited --------------------------
    if t.paralog_selectivity_blocked:
        score += 15.0
        flags.append("paralog selectivity unreachable by inhibitors")

    # --- glue-first tractability ----------------------------------------
    glue = t.degron_feature_score * 32.0
    if t.degron_feature_score >= 0.5:
        flags.append(f"glue-amenable surface (degron score {t.degron_feature_score:.2f})")
    else:
        # bifunctional fallback: the ligandability cascade
        ladder = {
            Ligandability.KNOWN_LIGAND: 18.0,
            Ligandability.REPORTED_POCKET: 12.0,
            Ligandability.PREDICTED_POCKET: 6.0,
            Ligandability.NONE: 0.0,
        }
        glue = max(glue, ladder[t.ligandability])
        flags.append(f"bifunctional fallback: {t.ligandability.value}")
        if t.ligandability is Ligandability.NONE:
            flags.append("no ligand route -- route to ligand-discovery queue")
    score += glue

    # --- turnover: slow resynthesis favours degradation ------------------
    if t.protein_half_life_h is not None and t.protein_half_life_h >= 24:
        score += 8.0
        flags.append(f"slow resynthesis (t-half {t.protein_half_life_h:.0f} h)")

    # --- safety: human LoF intolerance ------------------------------------
    if t.lof_constraint_loeuf is not None:
        if t.lof_constraint_loeuf < 0.35:
            score -= 25.0
            flags.append(
                f"strongly LoF-intolerant in humans (LOEUF {t.lof_constraint_loeuf:.2f}) "
                "-- heterozygous loss is already disease-causing"
            )
            reasons.add("lof_intolerant")
        elif t.lof_constraint_loeuf > 1.0:
            score += 5.0
            flags.append(f"LoF-tolerant in humans (LOEUF {t.lof_constraint_loeuf:.2f})")

    # --- delivery gate ----------------------------------------------------
    if t.tissue_of_action is TissueOfAction.HEPATIC:
        score += 6.0
        flags.append("hepatic restriction available -- delivery gate relaxed")
    elif t.tissue_of_action is TissueOfAction.LOCAL:
        score += 4.0
        flags.append("local delivery available -- delivery gate relaxed")

    return AxisResult(round(max(0.0, min(100.0, score)), 1), flags, reasons)


def _score_extracellular(t: Target) -> AxisResult:
    """Lysosomal-targeting track. Different chemistry, different bar.

    Two rulings shape this: a nomination needs a catalytic-turnover argument,
    and crowded competition is penalised hard. Without both, the track would
    hand back PCSK9 -- a target with two approved antibodies and an approved
    siRNA -- which is not a nomination, it is a reading of the news.
    """
    flags: list[str] = [f"{t.subcellular.value} protein -- lysosomal degradation track"]
    reasons: set[str] = {"extracellular_track"}
    score = 30.0  # base: the compartment is reachable by this chemistry

    if t.catalytic_advantage:
        score += 28.0
        flags.append("catalytic turnover vs occupancy-driven binding")
    if t.clears_non_neutralisable:
        score += 22.0
        flags.append("clears species that binding cannot neutralise")

    # Gate, not a penalty. The ruling was that an extracellular nomination
    # must beat what exists; a target with no capability argument does not get
    # to limp through on a high validity score.
    if not (t.catalytic_advantage or t.clears_non_neutralisable):
        return AxisResult(
            0.0,
            flags + ["no capability argument over an antibody or siRNA -- gated out"],
            reasons | {"no_displacement_argument"},
        )

    penalty = {
        CompetitorStage.NONE: 0.0,
        CompetitorStage.PRECLINICAL: 5.0,
        CompetitorStage.CLINICAL: 15.0,
        CompetitorStage.PHASE3: 30.0,
        CompetitorStage.APPROVED: 40.0,
    }[t.competitor_stage]
    if penalty:
        score -= penalty
        flags.append(f"competing modality at {t.competitor_stage.value} (-{penalty:.0f})")
        reasons.add("crowded")
        # A long-interval competitor also guts the dosing-frequency argument.
        if (
            t.catalytic_advantage
            and t.competitor_dosing_interval_days
            and t.competitor_dosing_interval_days >= 90
        ):
            score -= 15.0
            flags.append(
                f"competitor already dosed every {t.competitor_dosing_interval_days} days "
                "-- turnover argument does not buy an interval"
            )

    return AxisResult(round(max(0.0, min(100.0, score)), 1), flags, reasons)


# --------------------------------------------------------------- combine


def score(t: Target, profile_name: str | None = None) -> Result:
    p = PROFILES[profile_name or t.disease_area]
    v = score_validity(t, p)
    m = score_modality(t)

    flags: list[str] = []
    verdict: Verdict

    # Verdict precedence. A refuted target is reported as refuted even when the
    # compartment also rules it out: "this target is harmful" transfers to
    # other modalities, "we can't reach it" does not.
    if "validity_refuted" in v.reasons:
        verdict = Verdict.REJECT_TARGET_INVALID
    elif "scaffold_protective" in m.reasons:
        verdict = Verdict.REJECT_SCAFFOLD_PROTECTIVE
    elif "inhibition_sufficient" in m.reasons:
        verdict = Verdict.REJECT_INHIBITION_SUFFICIENT
    elif m.score == 0.0:
        verdict = Verdict.REJECT_WRONG_MODALITY
    elif m.score < p.modality_floor:
        verdict = (
            Verdict.REJECT_WRONG_MODALITY
            if "wrong_compartment" in m.reasons
            else Verdict.REJECT_MODALITY_TOO_WEAK
        )
        if "no_displacement_argument" in m.reasons:
            flags.append("extracellular target with no advantage over existing modalities")
        elif "crowded" in m.reasons:
            flags.append("displacement argument does not survive the competitive landscape")
    elif t.clinical_outcome is ClinicalOutcome.DISCONTINUED_STRATEGIC:
        verdict = Verdict.FLAG_STRATEGIC_DISCONTINUATION
    elif (
        v.score < p.mid_threshold
        and t.mr_status is MRStatus.SUPPORTED
        and t.clinical_outcome is ClinicalOutcome.FAILED_FUTILITY
    ):
        verdict = Verdict.FLAG_CONTESTED
        flags.append(
            "MR support stands and the trial failure is explicable by threshold and power "
            "-- resolve before scoring either way"
        )
    elif v.score < p.mid_threshold:
        if m.score >= p.strong_modality and v.score >= p.validation_floor:
            verdict = Verdict.NOMINATE_FOR_VALIDATION
            flags.append(
                "degradation is the right tool; validity needs a perturbation experiment, "
                "not more mining"
            )
        else:
            verdict = Verdict.REJECT_TARGET_INVALID
    else:
        # Geometric mean: a target must clear BOTH axes. A high score on one
        # cannot rescue a low score on the other, which is the whole point.
        combined = (v.score * m.score) ** 0.5
        verdict = Verdict.TOP if combined >= p.top_threshold else Verdict.MID

    # ZEUS ruling: genetics alone can reach TOP, but the exposure is named.
    if verdict is Verdict.TOP and t.clinical_outcome in (
        ClinicalOutcome.UNKNOWN,
        ClinicalOutcome.PHASE1_BIOMARKER,
    ):
        flags.append(
            "unanchored by perturbation data -- genetics-driven nomination, track how it ages"
        )

    composite = round((v.score * m.score) ** 0.5, 1)
    return Result(t.target_id, t.gene, t.indication, v, m, verdict, t.track, composite, flags)
