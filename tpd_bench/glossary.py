"""Plain-language and detailed readings of every flag the bench emits.

Each entry has three parts: what it means in one sentence, what is actually
being measured, and where the underlying data or claim comes from. A score
nobody can interrogate is a score nobody should act on.

Matching is by substring against the flag text, longest pattern first, so a
specific flag wins over a general one.
"""

from __future__ import annotations

# Sources referenced below, defined once.
SRC = {
    "gnomad": ("gnomAD v4 constraint metrics", "https://gnomad.broadinstitute.org"),
    "ot": ("Open Targets Platform association scores", "https://platform.opentargets.org"),
    "af": ("AlphaFold Protein Structure Database", "https://alphafold.ebi.ac.uk"),
    "uniprot": ("UniProt subcellular location annotations", "https://www.uniprot.org"),
    "chembl": ("ChEMBL mechanism and activity data", "https://www.ebi.ac.uk/chembl/"),
    "ctgov": ("ClinicalTrials.gov registry", "https://clinicaltrials.gov"),
    "crbn": ("Mining the CRBN target space, Science 2025",
             "https://www.science.org/doi/10.1126/science.adt6736"),
    "vav1": ("CRBN glues reach VAV1 via an SH3 RT-loop degron, not a G-loop",
             "https://www.researchgate.net/publication/392621760"),
    "ripk1": ("RIPK1 kinase-dead knock-in versus knockout, Front Immunol 2021;12:676946",
              "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.676946/full"),
    "btk": ("BTK resistance mutants retain signalling without kinase activity, Science 2024",
            "https://www.nurixtx.com/wp-content/uploads/2023/07/NX-2127-Scaffholding-Function-ICML-2023-poster-1.pdf"),
    "hpk1": ("Kinase-dead HPK1 phenocopies knockout, Front Immunol 2025;16:1449106",
             "https://pubmed.ncbi.nlm.nih.gov/39981246/"),
    "zeus": ("ZEUS: IL-6 inhibition returned HR 0.99 despite strong Mendelian randomisation",
             "https://ecgwaves.com/articles/zeus-and-lpa-horizon-trials-why-ziltivekimab-and-pelacarsen-failed"),
    "lppla2": ("Lp-PLA2 genetic estimate matched the trial at equal target engagement, "
               "Eur J Prev Cardiol 2017", "https://pubmed.ncbi.nlm.nih.gov/27940953"),
    "oral": ("Physicochemical determinants of oral absorption for degraders, J Med Chem 2023",
             "https://www.arvinas.com/wp-content/uploads/2024/07/Physicochemical-Property-Determinants-of-Oral-Absorption-for-PROTAC-Protein-Degraders.pdf"),
    "design": ("Bench design decision, not an external finding", ""),
}


def _e(simple: str, detail: str, src: str) -> dict:
    name, url = SRC[src]
    return {"simple": simple, "detail": detail, "source": name, "url": url}


# Longest pattern wins, so order here does not matter.
GLOSSARY: dict[str, dict] = {
    "direction of effect unknown": _e(
        "We know this gene is linked to the disease, but not whether removing the protein would help.",
        "Association scores say a gene matters; they do not say whether the associated variants "
        "reduce or increase protein function. Degradation only follows from a protective "
        "loss-of-function. The score is left intact rather than discounted, because unknown is "
        "not the same as contrary.",
        "ot"),
    "genetic direction does not support target removal": _e(
        "The genetics point the other way: losing this protein looks unhelpful.",
        "Where the protective variants increase rather than reduce function, removing the protein "
        "argues against the human evidence. Genetics is discounted to 40% of its weight.",
        "ot"),
    "mr_refuted": _e(
        "The statistical genetics say this marker tracks the disease without causing it.",
        "Mendelian randomisation finds no causal effect despite a strong observational "
        "association. HDL-C, Lp-PLA2, sPLA2 and CRP all had this shape, and all produced null "
        "or harmful trials.",
        "lppla2"),
    "MR supported but not standardised": _e(
        "The genetic case is real but has not been matched to the dose a drug can actually reach.",
        "A causal estimate only predicts a trial when it is standardised to the target engagement "
        "the drug achieves. ZEUS had multi-instrument support and still returned a null, so "
        "unmatched estimates are discounted.",
        "zeus"),
    "MR standardised to achievable engagement": _e(
        "The genetic estimate was matched to the effect a drug can actually deliver.",
        "Standardising the causal estimate to achievable target engagement is what separates a "
        "prediction from an association. Done properly it predicted the Lp-PLA2 trials almost exactly.",
        "lppla2"),
    "no human genetic support": _e(
        "No human genetic evidence links this gene to this disease.",
        "Recorded as a flag rather than a veto. Several approved drugs had no genetic nomination, "
        "so absence of genetics lowers confidence without disqualifying a target.",
        "ot"),
    "no clinical data yet": _e(
        "Nobody has run a trial on this target yet, so there is nothing to score.",
        "The term is removed from the weighted average and the remaining weights are "
        "renormalised. Scoring 'untested' as weak evidence would make the tool unable to "
        "nominate anything new.",
        "design"),
    "pathway drug status not looked up": _e(
        "We have not checked whether an approved drug already hits this pathway.",
        "Not looked up in this run. Renormalised out rather than scored as absent.",
        "design"),
    "disease expression not looked up": _e(
        "We have not checked whether this protein is enriched in diseased tissue.",
        "Not looked up in this run. Renormalised out rather than scored as absent.",
        "design"),
    "score rests on": _e(
        "This score comes from a single kind of evidence, so treat it as provisional.",
        "Renormalising missing terms is correct, but a high score built on one layer is thin. "
        "The number is arithmetically right and epistemically weak.",
        "design"),
    "trial showed harm": _e(
        "A trial of this target made patients worse.",
        "Validity is set to zero regardless of other evidence. Harm transfers across modalities, "
        "so it is not something a different drug design fixes.",
        "ctgov"),
    "surrogate approval only": _e(
        "Approved on a biomarker, not on whether patients actually do better.",
        "Obeticholic acid hit its histologic endpoint and was still rejected on benefit-risk. "
        "Surrogate approval is kept as a distinct, weaker state.",
        "ctgov"),
    "strategic discontinuation": _e(
        "The programme stopped for business reasons, not because the science failed.",
        "Scored near neutral. Treating a portfolio decision as a scientific failure would teach "
        "the tool that good targets are bad.",
        "ctgov"),
    "highly redundant target class": _e(
        "Other proteins can do this one's job, so removing it may change nothing.",
        "Chemokine receptors are the worst case: roughly 50 ligands across 20 promiscuous "
        "receptors. Redundancy is a property of the class, not of the experiment.",
        "design"),
    "node(s) from the validated one": _e(
        "The proven target is nearby in the pathway, but not this protein.",
        "Adjacent nodes can carry opposite signs. IL-23 blockade helps in Crohn's disease and "
        "IL-17 blockade harms, one step apart on the same axis.",
        "ot"),
    "scaffold is protective": _e(
        "Removing this protein is worse than blocking it. Do not degrade it.",
        "RIPK1 kinase-dead mice are normal while knockouts die postnatally: the non-catalytic "
        "role holds a checkpoint against cell death. Evidence of a scaffold is not an argument "
        "for degradation until its direction is known.",
        "ripk1"),
    "kinase-dead phenocopies knockout": _e(
        "Switching the enzyme off does everything removing the protein would. An inhibitor is enough.",
        "In HPK1, kinase-dead knock-ins match knockouts and a degrader matches an inhibitor "
        "across TCR activation readouts. There is no scaffold left for degradation to address.",
        "hpk1"),
    "a tolerated inhibitor already works": _e(
        "A safe, effective inhibitor exists. Degradation would have to beat it.",
        "Deucravacitinib solved TYK2 selectivity allosterically with no class boxed warning. "
        "Where inhibition is already clean, degradation adds cost without a clear gain.",
        "chembl"),
    "scaffold evidence tier": _e(
        "This protein has a job beyond its enzyme activity, which an inhibitor would not touch.",
        "Tiered by how much the comparison can be trusted. Human resistance mutations that "
        "ablate catalysis while signalling continues are the strongest tier; structural "
        "interaction evidence the weakest.",
        "btk"),
    "scaffold sign unknown": _e(
        "The protein has a non-enzyme role, but we do not know if removing it helps or hurts.",
        "Half credit pending direction. The sign matters more than the existence of the scaffold.",
        "design"),
    "target accumulates under inhibition": _e(
        "Blocking this protein makes the cell produce more of it. Removal avoids that.",
        "Compensatory upregulation blunts an inhibitor over time. Degradation removes the "
        "protein rather than occupying it, so the feedback has nothing to push against.",
        "chembl"),
    "no tractable inhibitor": _e(
        "Nothing can be made to block this protein, so degradation is the only route.",
        "Transcription factors and adaptors mostly lack the deep pocket an inhibitor needs. "
        "An induced interface does not require one.",
        "chembl"),
    "no pocket and no inhibitor": _e(
        "No pocket to inhibit and no competing approach. The purest case for a molecular glue.",
        "Nothing to inhibit and no bifunctional route either, which leaves an induced-proximity "
        "glue as the only option worth pursuing.",
        "crbn"),
    "inhibitor engaged target without efficacy": _e(
        "A drug hit this protein and patients did not improve.",
        "Weak support for degradation. Target engagement without benefit usually means the "
        "target is wrong, not that the modality was.",
        "ctgov"),
    "inhibitor works but window is narrow": _e(
        "An inhibitor works but the safe dose range is tight. Selectivity could widen it.",
        "A narrow therapeutic index is a selectivity problem, which induced proximity can "
        "sometimes address where the binding site cannot.",
        "chembl"),
    "paralog selectivity unreachable": _e(
        "Near-identical sister proteins share the same pocket, so an inhibitor cannot tell them apart.",
        "Degraders can discriminate on surface complementarity rather than pocket shape. "
        "SMARCA2 versus SMARCA4 is the worked example.",
        "crbn"),
    "glue-amenable surface": _e(
        "The protein has an exposed loop a molecular glue could grip.",
        "Scored from the AlphaFold model: an exposed, ordered beta-hairpin turn, a glycine at "
        "the turn, an aromatic nearby, and a reachable lysine for ubiquitin transfer. Graded "
        "rather than gated, because CRBN glues also reach proteins through non-canonical loops.",
        "af"),
    "no ligand route": _e(
        "No pocket and no glue surface found. Chemistry would have to be invented first.",
        "Routed to a ligand-discovery queue rather than rejected. Absence of chemistry is not "
        "absence of biology.",
        "af"),
    "bifunctional fallback": _e(
        "No glue surface, so this would need a two-headed molecule built on a known pocket.",
        "Falls back through the ligandability cascade: a known ligand, then a reported pocket, "
        "then a pocket predicted from structure.",
        "chembl"),
    "strongly LoF-intolerant": _e(
        "People born with one broken copy of this gene are already ill. Removing it is risky.",
        "LOEUF below 0.35 means strong selection against loss of function in healthy humans. "
        "This is a safety signal, separate from whether the target drives disease.",
        "gnomad"),
    "LoF-tolerant in humans": _e(
        "People missing one copy of this gene are healthy, so removing the protein looks tolerable.",
        "A high LOEUF means loss of function is not selected against. Reassuring about safety; "
        "it says nothing about efficacy.",
        "gnomad"),
    "slow resynthesis": _e(
        "The cell replaces this protein slowly, so one dose keeps working longer.",
        "Slow turnover favours degradation and disfavours occupancy-driven inhibition, since "
        "the protein does not come straight back.",
        "design"),
    "hepatic restriction available": _e(
        "The drug can be aimed at the liver, which limits exposure elsewhere.",
        "First-pass metabolism and ASGPR-mediated uptake make hepatic restriction close to free, "
        "so the delivery bar is relaxed for liver targets.",
        "design"),
    "local delivery available": _e(
        "The drug can be delivered where it is needed rather than throughout the body.",
        "Inhaled, topical, gut-restricted and intra-articular routes each buy selectivity "
        "without requiring it from the molecule.",
        "design"),
    "outside intracellular degradation": _e(
        "This protein sits outside the cell, where the proteasome cannot reach it.",
        "Routed to the lysosomal track instead, which uses a different chemistry and is judged "
        "against a different delivery bar.",
        "uniprot"),
    "lysosomal degradation track": _e(
        "Outside the cell, so this would need a lysosome-targeting conjugate rather than a glue.",
        "Membrane and secreted proteins are reachable, but by antibody or glycopeptide "
        "conjugates that are not orally dosed.",
        "uniprot"),
    "no capability argument": _e(
        "An antibody or siRNA already does this job as well. Degradation adds nothing.",
        "An extracellular nomination has to offer something binding cannot: catalytic turnover, "
        "or clearance of a species that neutralisation leaves behind.",
        "ctgov"),
    "catalytic turnover": _e(
        "One drug molecule can clear many copies, unlike an antibody that must occupy each one.",
        "Sub-stoichiometric clearance is the core mechanistic argument for a degrader over an "
        "occupancy-driven binder.",
        "design"),
    "competing modality at": _e(
        "Someone already has a drug against this target at this stage.",
        "Penalised by how advanced the competitor is. An approved twice-yearly siRNA is a much "
        "harder thing to displace than a preclinical programme.",
        "ctgov"),
    "does not buy an interval": _e(
        "The existing drug is already dosed rarely, so turnover gains nothing.",
        "PCSK9 is the example: the catalytic argument is real, but inclisiran is already "
        "twice-yearly, so there is no dosing interval left to win.",
        "ctgov"),
    "degradation is the right tool": _e(
        "Degradation clearly fits, but whether the target drives disease is still unproven.",
        "A distinct verdict rather than a rejection. What is needed next is a perturbation "
        "experiment in disease-relevant cells, not more database mining.",
        "design"),
    "unanchored by perturbation data": _e(
        "This nomination rests on genetics with no experiment behind it. Watch how it ages.",
        "Genetics-driven nominations are allowed through, but the exposure is named so it can "
        "be tracked rather than forgotten.",
        "zeus"),
    "MR support stands and the trial failure": _e(
        "The trial failed, but it may have been too small or dosed at the wrong threshold to tell.",
        "Lp(a)HORIZON missed its endpoint while the genetic case stood. The causal estimate "
        "predicts benefit per absolute unit lowered, and at the entry threshold used the "
        "expected effect was likely below what the trial could detect. Held as contested "
        "rather than scored either way.",
        "zeus"),
    "clears species that binding cannot neutralise": _e(
        "Some targets cannot be neutralised by binding them \u2014 only by removing them.",
        "Aggregates, lipoprotein particles and complexes can keep functioning with an antibody "
        "attached. Removal addresses what occupancy does not.",
        "design"),
    "displacement argument does not survive": _e(
        "The advantage over existing drugs disappears once the competition is accounted for.",
        "A capability argument has to beat what is already approved or in Phase 3, not merely "
        "exist in principle.",
        "ctgov"),
    "oral absorption": _e(
        "Degraders are large molecules, and oral dosing constrains how they can be built.",
        "Unsatisfied hydrogen-bond donors are the binding constraint; molecular weight and "
        "polar surface area relax further than the Rule of 5 would suggest.",
        "oral"),
}


def explain(flag: str) -> dict | None:
    """Longest matching pattern wins, so a specific flag beats a general one."""
    hits = [k for k in GLOSSARY if k.lower() in flag.lower()]
    if not hits:
        return None
    return GLOSSARY[max(hits, key=len)]
