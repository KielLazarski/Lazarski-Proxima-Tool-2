"""Proteome universe and open-source feature adapters.

Universe: every human protein UniProt places in the nucleus or cytosol, plus
membrane and secreted proteins routed to the extracellular track. Roughly 12k
either way, and every feature layer has to cover it before a nomination from
this sweep means anything.

Open sources only, per the data-stack ruling:
  UniProt            subcellular location, sequence, active sites
  Open Targets       genetic association, direction of effect
  ChEMBL             ligandability, clinical phase of known binders
  ClinicalTrials.gov competitor landscape
  AlphaFold DB       structures for the degron scorer
  Human Protein Atlas / GTEx   expression

THE WEAK LINK
`inhibitor_status` moves more verdicts than any other term in the bench, and it
is the one that degrades worst under automation. The distinctions that matter
-- an inhibitor that engaged its target and failed, versus one that never
engaged, versus one that works and is tolerated -- are clinical facts written
in prose. ChEMBL knows that binders exist and at what phase; it does not know
why a programme stopped.

`derive_inhibitor_status` therefore returns a status AND a confidence, and
anything below `CONFIDENCE_FLOOR` is forced to NONE_TRIED rather than guessed.
A wrong confident answer on this field moves more verdicts than a wrong answer
anywhere else, so the bias is deliberately toward admitting ignorance.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from typing import Any, Iterator

from . import cache, http
from ..schema import (
    CompetitorStage,
    GeneticsClass,
    InhibitorStatus,
    Ligandability,
    Subcellular,
)

TIMEOUT = 30
CONFIDENCE_FLOOR = 0.5
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"


class NetworkUnavailable(RuntimeError):
    """Raised rather than returning empty, so a blocked run cannot look clean."""


def _get(url: str) -> Any:
    hit = cache.get(url)
    if hit is not None:
        return hit
    try:
        data = http.request(url)
        cache.put(url, data)
        return data
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise NetworkUnavailable(
            f"could not reach {urllib.parse.urlparse(url).netloc} "
            "-- run outside the sandbox"
        ) from exc


# --------------------------------------------------------------- universe

# UniProt subcellular vocabulary -> our four compartments. Order matters:
# a protein annotated both nuclear and secreted is treated as secreted,
# because the extracellular pool is what a degrader would have to reach.
# Topology phrases that mean the protein genuinely spans the plasma membrane
# with an extracellular face. Anything else that merely mentions "membrane" is
# usually a cytosolic protein that associates with one.
TRANSMEMBRANE = (
    "Single-pass", "Multi-pass", "Type I membrane", "Type II membrane",
    "Type III membrane", "Type IV membrane",
)
# Membrane association from the cytoplasmic side. Still intracellular, still
# reachable by the proteasome.
PERIPHERAL = ("Peripheral membrane protein", "Lipid-anchor", "Anchor")


def human_proteome() -> Iterator[dict]:
    """Stream reviewed human entries with location and sequence."""
    base = (
        "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:9606"
        "&fields=accession,gene_primary,cc_subcellular_location,sequence,ft_act_site,xref_ensembl"
        "&format=json&size=500"
    )
    url = base
    while url:
        data = _get(url)
        for e in data.get("results", []):
            yield e
        url = data.get("nextUrl")  # cursor pagination


def proteome_by_gene(genes: list[str]) -> Iterator[dict]:
    """Fetch named genes only. For checking whether a layer works at all:
    50 arbitrary accessions tell you nothing about the genetics layer, because
    obscure proteins genuinely have no disease associations."""
    for g in genes:
        q = urllib.parse.quote(f'gene:{g} AND reviewed:true AND organism_id:9606')
        url = (f"https://rest.uniprot.org/uniprotkb/search?query={q}"
               "&fields=accession,gene_primary,cc_subcellular_location,sequence,"
               "ft_act_site,xref_ensembl&format=json&size=1")
        try:
            data = _get(url)
        except NetworkUnavailable:
            raise
        for e in data.get("results", []):
            yield e


def classify_location(entry: dict) -> Subcellular | None:
    """Route a protein to a compartment, or None when UniProt does not say.

    The first version of this looked for the word "Membrane" anywhere and sent
    the protein to the extracellular track. That misrouted 46% of a test sweep:
    ARHGAP10, SNAP23 and HAX1 are cytosolic proteins that associate with a
    membrane, and the proteasome reaches all of them. Only a genuine
    transmembrane topology with an extracellular face belongs on the lysosomal
    track.

    Precedence: secreted, then true transmembrane, then nucleus, then cytosol,
    then any remaining membrane. An unplaced protein returns None and is never
    assumed cytosolic.
    """
    locations: list[str] = []
    topologies: list[str] = []
    for c in entry.get("comments", []):
        if c.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for loc in c.get("subcellularLocations", []):
            v = (loc.get("location") or {}).get("value")
            if v:
                locations.append(v)
            t = (loc.get("topology") or {}).get("value")
            if t:
                topologies.append(t)

    if not locations:
        return None

    blob = "; ".join(locations)
    topo = "; ".join(topologies)

    if "Secreted" in blob:
        return Subcellular.SECRETED

    # Genuine transmembrane protein with an extracellular face
    if any(k in topo for k in TRANSMEMBRANE) and (
        "Cell membrane" in blob or "Plasma membrane" in blob
    ):
        return Subcellular.MEMBRANE

    if "Nucleus" in blob:
        return Subcellular.NUCLEAR
    if "Cytoplasm" in blob or "Cytosol" in blob:
        return Subcellular.CYTOSOLIC

    # Peripheral or lipid-anchored association is cytoplasmic-side
    if any(k in topo for k in PERIPHERAL):
        return Subcellular.CYTOSOLIC

    # An internal membrane (mitochondrion, ER, Golgi, endosome) is intracellular.
    # Case-insensitive: UniProt writes "Endoplasmic reticulum membrane".
    low = blob.lower()
    if "membrane" in low and "cell membrane" not in low:
        return Subcellular.CYTOSOLIC
    if "cell membrane" in low:
        return Subcellular.MEMBRANE
    return None


def ensembl_gene_id(entry: dict) -> str | None:
    """ENSG id from the UniProt cross-references.

    Open Targets is keyed on Ensembl gene ids. The first sweep passed it
    UniProt accessions, which match nothing, so every genetics lookup came
    back empty and every validity score was zero.
    """
    for x in entry.get("uniProtKBCrossReferences", []):
        if x.get("database") != "Ensembl":
            continue
        for p in x.get("properties", []):
            if p.get("key") == "GeneId":
                return (p.get("value") or "").split(".")[0] or None
    return None


# --------------------------------------------------------------- genetics

def target_associations(ensembl_id: str) -> dict[str, dict]:
    """All disease associations for a gene, keyed by disease id.

    Deliberately unfiltered. The server-side filter argument is called `Bs`,
    which is undocumented-looking enough that relying on it is a liability --
    and fetching once per gene costs fewer requests than once per
    gene-and-indication anyway. Disease ids are MONDO, not EFO.
    """
    q = """query($g:String!,$i:Int!,$s:Int!){ target(ensemblId:$g){
             associatedDiseases(page:{index:$i,size:$s}){
               count rows{ score datatypeScores{ id score } disease{ id } } } } }"""
    out: dict[str, dict] = {}
    size, index = 500, 0
    while index < 6:  # hard ceiling; no gene needs 3000 associations
        payload = {"query": q, "variables": {"g": ensembl_id, "i": index, "s": size}}
        key = f"{OT_URL}#assoc#{ensembl_id}#{index}"
        hit = cache.get(key)
        if hit is None:
            hit = http.request(OT_URL, payload)
            cache.put(key, hit)
        ad = (((hit.get("data") or {}).get("target") or {}).get("associatedDiseases") or {})
        rows = ad.get("rows") or []
        for r in rows:
            did = (r.get("disease") or {}).get("id")
            if did:
                out[did] = {
                    "score": r.get("score", 0.0),
                    "datatypes": {d["id"]: d["score"] for d in r.get("datatypeScores", [])},
                }
        if len(rows) < size or len(out) >= (ad.get("count") or 0):
            break
        index += 1
    return out


def open_targets_genetics(ensembl_id: str, disease_id: str) -> tuple[GeneticsClass, bool]:
    """Association class for one gene-disease pair, plus direction.

    Direction of effect is still unresolved: the association score says a gene
    matters, not whether protective variants are loss-of-function. Degradation
    only follows from that, so this returns False and genetics stays
    discounted. gnomAD LOEUF covers a different question (safety), not this one.
    """
    assoc = target_associations(ensembl_id).get(disease_id)
    if not assoc:
        return GeneticsClass.NONE, False

    rare = assoc["datatypes"].get("genetic_association", 0.0)
    overall = assoc["score"]

    if rare >= 0.7:
        cls = GeneticsClass.LOF_PROTECTIVE_REPLICATED
    elif rare >= 0.5:
        cls = GeneticsClass.LOF_PROTECTIVE
    elif rare >= 0.3:
        cls = GeneticsClass.RARE_VARIANT_ASSOCIATION
    elif overall > 0:
        cls = GeneticsClass.GWAS_ASSOCIATION
    else:
        cls = GeneticsClass.NONE
    return cls, False


# --------------------------------------------------------------- pharmacology

@dataclass
class InhibitorCall:
    status: InhibitorStatus
    ligandability: Ligandability
    confidence: float
    basis: str


def derive_inhibitor_status(chembl_target_id: str) -> InhibitorCall:
    """Best-effort inhibitor status from ChEMBL, with explicit confidence.

    What ChEMBL supports well: whether potent binders exist, and the highest
    clinical phase reached. What it cannot tell you: whether a programme
    stopped for efficacy, toxicity or portfolio reasons -- which is exactly the
    distinction the bench leans on hardest.

    So a drug that reached phase 4 is read as WORKS_TOLERATED with reasonable
    confidence, and anything that stalled earlier returns low confidence and is
    downgraded to NONE_TRIED by the caller.
    """
    data = _get(
        "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
        f"?target_chembl_id={chembl_target_id}&limit=200"
    )
    mechs = data.get("mechanisms", [])
    acts = _get(
        "https://www.ebi.ac.uk/chembl/api/data/activity.json"
        f"?target_chembl_id={chembl_target_id}&standard_type=IC50&limit=1"
    )
    has_potent = acts.get("page_meta", {}).get("total_count", 0) > 0

    ligand = Ligandability.KNOWN_LIGAND if has_potent else Ligandability.NONE

    if not mechs:
        return InhibitorCall(
            InhibitorStatus.NO_TRACTABLE_INHIBITOR if not has_potent
            else InhibitorStatus.NONE_TRIED,
            ligand,
            0.6 if not has_potent else 0.4,
            "no mechanism-of-action entries in ChEMBL",
        )

    phase = max((m.get("max_phase") or 0) for m in mechs)
    if phase >= 4:
        return InhibitorCall(
            InhibitorStatus.WORKS_TOLERATED, ligand, 0.75,
            f"approved drug against this target (max_phase {phase})",
        )
    # Everything between is genuinely ambiguous: a phase-2 compound may have
    # failed, stalled, or still be running. ChEMBL does not distinguish them.
    return InhibitorCall(
        InhibitorStatus.NONE_TRIED, ligand, 0.25,
        f"clinical compounds exist (max_phase {phase}) but outcome is not in ChEMBL",
    )


def apply_confidence(call: InhibitorCall) -> InhibitorStatus:
    """Admit ignorance below the floor rather than guess the heaviest term."""
    return call.status if call.confidence >= CONFIDENCE_FLOOR else InhibitorStatus.NONE_TRIED


# --------------------------------------------------------------- competition

PHASE_TO_STAGE = {
    0: CompetitorStage.NONE,
    1: CompetitorStage.CLINICAL,
    2: CompetitorStage.CLINICAL,
    3: CompetitorStage.PHASE3,
    4: CompetitorStage.APPROVED,
}


class CompetitorDataUnavailable(RuntimeError):
    """Raised when the competitive landscape cannot be established.

    Returning NONE here would be a silent, directional error: no competition
    means no penalty, which makes a crowded target look attractive. The
    extracellular track rejects almost everything on this term, so a quiet
    default would invert its output. Better to stop the sweep.
    """


def competitor_stage(gene: str, condition: str) -> CompetitorStage:
    """Most advanced interventional trial against this target-indication pair.

    Counts siRNA, ASO and antibody interventions, since those are what an
    extracellular degrader would have to displace.
    """
    url = (
        "https://clinicaltrials.gov/api/v2/studies?"
        + urllib.parse.urlencode(
            {
                "query.intr": gene,
                "query.cond": condition,
                "filter.overallStatus": "RECRUITING|ACTIVE_NOT_RECRUITING|COMPLETED",
                "pageSize": "100",
            }
        )
    )
    try:
        data = _get(url)
    except NetworkUnavailable as exc:
        raise CompetitorDataUnavailable(
            f"cannot establish competition for {gene}/{condition}; refusing to "
            "assume the field is clear"
        ) from exc
    best = 0
    for s in data.get("studies", []):
        phases = s.get("protocolSection", {}).get("designModule", {}).get("phases", [])
        for p in phases:
            digits = "".join(c for c in p if c.isdigit())
            if digits:
                best = max(best, int(digits))
    return PHASE_TO_STAGE.get(best, CompetitorStage.NONE)
