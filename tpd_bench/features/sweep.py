"""The sweep. Builds Target records from live sources and scores them.

    py -m tpd_bench.features.sweep --limit 50     # smoke test, minutes
    py -m tpd_bench.features.sweep                # full proteome, hours

Start with --limit. A full sweep is tens of thousands of requests against
rate-limited public APIs and there is no point discovering a bad field mapping
on protein 9,000.

Everything caches, so a second run resumes rather than restarting, and the
smoke test's work is not thrown away when you go full scale.

WHAT THIS DOES NOT DO
Direction of effect is still unresolved. Open Targets gives an association
score, not whether the associated variants are protective loss-of-function or
protective gain-of-function, and degradation only follows from the first. Every
record therefore carries direction_supports_removal=False, which discounts
genetics to 40%. gnomAD fills a different hole -- a LoF-intolerance safety gate
-- but it does not fill this one.

Until that is closed, read genetics-driven ranks as systematically conservative.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from dataclasses import replace

from ..scoring import Verdict, score
from ..schema import (
    GeneticsClass,
    InhibitorStatus,
    Ligandability,
    Subcellular,
    Target,
    TissueOfAction,
    Track,
)
from . import cache, http
from .degron import score_protein
from .structure import to_residues
from .universe import (
    NetworkUnavailable,
    apply_confidence,
    classify_location,
    derive_inhibitor_status,
    ensembl_gene_id,
    human_proteome,
    proteome_by_gene,
    open_targets_genetics,
)

STATS: dict = {"genetics_errors": 0, "no_efo": 0, "no_ensembl": 0}

OT = "https://api.platform.opentargets.org/api/v4/graphql"
GNOMAD = "https://gnomad.broadinstitute.org/api"
AF = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# Indication per disease area. Resolved to EFO ids at runtime rather than
# hard-coded, because a wrong id fails silently as "no association".
INDICATIONS = {
    "immunology_inflammation": ["atopic dermatitis", "systemic lupus erythematosus"],
    "cardiometabolic": ["non-alcoholic steatohepatitis", "coronary artery disease"],
}


def resolve_disease(name: str) -> str | None:
    """Disease name to id. Open Targets returns MONDO ids, not EFO."""
    q = {
        "query": 'query($s:String!){search(queryString:$s,entityNames:["disease"])'
                 "{hits{id name entity}}}",
        "variables": {"s": name},
    }
    key = f"{OT}#efo#{name}"
    hit = cache.get(key)
    if hit is None:
        try:
            hit = http.request(OT, q)
        except Exception:
            return None
        cache.put(key, hit)
    hits = (((hit or {}).get("data") or {}).get("search") or {}).get("hits") or []
    for h in hits:
        if h.get("entity") == "disease":
            return h.get("id")
    return None


def loeuf(gene: str) -> float | None:
    """gnomAD LoF-intolerance. None when the gene is not covered."""
    q = {
        "query": "query($g:String!){gene(gene_symbol:$g,reference_genome:GRCh38)"
                 "{gnomad_constraint{oe_lof_upper}}}",
        "variables": {"g": gene},
    }
    key = f"{GNOMAD}#{gene}"
    hit = cache.get(key)
    if hit is None:
        try:
            hit = http.request(GNOMAD, q)
        except Exception:
            return None
        cache.put(key, hit)
    c = ((((hit or {}).get("data") or {}).get("gene") or {}).get("gnomad_constraint") or {})
    return c.get("oe_lof_upper")


def degron_for(acc: str) -> float:
    """Fetch the AlphaFold model and score its best degron. 0.0 if unavailable."""
    key = f"{AF.format(acc=acc)}#degron"
    hit = cache.get(key)
    if hit is not None:
        return hit.get("score", 0.0)
    try:
        meta = http.request(AF.format(acc=acc))
        rec = meta[0] if isinstance(meta, list) and meta else None
        url = rec and (rec.get("pdbUrl") or rec.get("cifUrl"))
        if not url or not url.endswith(".pdb"):
            cache.put(key, {"score": 0.0, "why": "no pdb model"})
            return 0.0
        import urllib.request as ur
        req = ur.Request(url, headers={"User-Agent": "tpd-bench/1.0"})
        with ur.urlopen(req, timeout=60, context=http.ssl_context()) as r:
            text = r.read().decode()
        res = to_residues(text)
        s, _hits = score_protein(res) if res else (0.0, [])
        cache.put(key, {"score": s})
        return s
    except Exception:
        cache.put(key, {"score": 0.0, "why": "fetch or parse failed"})
        return 0.0


def build_target(entry: dict, area: str, indication: str, efo: str | None) -> Target | None:
    acc = entry.get("primaryAccession")
    genes = entry.get("genes") or []
    gene = (genes[0].get("geneName", {}) or {}).get("value") if genes else None
    if not acc or not gene:
        return None
    loc = classify_location(entry)
    if loc is None:
        return None  # unplaced proteins are not scored, never assumed cytosolic

    gclass, direction = (GeneticsClass.NONE, False)
    ensg = ensembl_gene_id(entry)
    if efo and not ensg:
        STATS["no_ensembl"] += 1
    if efo and ensg:
        try:
            gclass, direction = open_targets_genetics(ensg, efo)
        except Exception as exc:
            # The first sweep swallowed this and produced 140 rows of
            # validity 0.0 with no indication anything had gone wrong.
            STATS["genetics_errors"] += 1
            STATS.setdefault("genetics_error_example", f"{type(exc).__name__}: {exc}")
    else:
        STATS["no_efo"] += 1

    t = Target(
        target_id=acc,
        gene=gene,
        indication=indication,
        disease_area=area,
        subcellular=loc,
        tissue_of_action=TissueOfAction.HEPATIC if area == "cardiometabolic" else TissueOfAction.SYSTEMIC,
        genetics=gclass,
        genetics_direction_supports_removal=direction,
        # The sweep populates genetics, structure and constraint. It does not
        # look up pathway drugs, disease expression, or direction of effect.
        unknown_fields=frozenset({"pathway", "expression", "direction"}),
        lof_constraint_loeuf=loeuf(gene),
    )
    if t.track is Track.INTRACELLULAR:
        t = replace(t, degron_feature_score=degron_for(acc))
    return t


def run(limit: int | None, area: str | None, genes: list[str] | None = None) -> int:
    areas = [area] if area else list(INDICATIONS)
    print(f"sweeping: {', '.join(areas)}")
    if limit:
        print(f"--limit {limit}: smoke test, not a real sweep\n")

    plan = []
    for a in areas:
        for name in INDICATIONS[a]:
            efo = resolve_disease(name)
            print(f"  {a:26s} {name:34s} id: {efo or 'UNRESOLVED'}")
            plan.append((a, name, efo))
    print()

    rows, seen, t0 = [], 0, time.time()
    try:
        source = proteome_by_gene(genes) if genes else human_proteome()
        for entry in source:
            seen += 1
            for a, name, efo in plan:
                t = build_target(entry, a, name, efo)
                if t is None:
                    continue
                r = score(t)
                rows.append((t, r))
            if seen % 25 == 0:
                print(f"  {seen} proteins | {len(rows)} scored | {time.time()-t0:.0f}s")
            if limit and seen >= limit:
                break
    except NetworkUnavailable as exc:
        print(f"\nBLOCKED: {exc}")
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted -- cached work is kept, rerun to resume")

    out = "sweep_results.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([
            "accession", "gene", "indication", "area", "track", "compartment",
            "validity", "modality", "composite", "verdict", "loeuf",
            "degron", "genetics", "flags",
        ])
        for t, r in rows:
            w.writerow([
                t.target_id, t.gene, t.indication, t.disease_area, r.track.value,
                t.subcellular.value, r.validity.score, r.modality.score,
                r.composite, r.verdict.value,
                "" if t.lof_constraint_loeuf is None else round(t.lof_constraint_loeuf, 3),
                round(t.degron_feature_score, 3), t.genetics.value,
                " | ".join(r.validity.flags + r.modality.flags),
            ])

    # A column that is constant across every row carries no information. Say so
    # loudly: the first sweep wrote 140 rows of validity 0.0 and looked fine.
    problems = []
    if rows:
        if all(r.validity.score == 0.0 for _t, r in rows):
            problems.append(
                "VALIDITY IS 0.0 FOR EVERY ROW -- nothing can be nominated. "
                "The validity layer is supplying no data."
            )
        if all(t.genetics is GeneticsClass.NONE for t, _r in rows):
            problems.append("GENETICS IS 'none' FOR EVERY ROW -- Open Targets returned nothing.")
        intra = [t for t, _r in rows if t.track is Track.INTRACELLULAR]
        if intra and all(t.degron_feature_score == 0.0 for t in intra):
            problems.append("ALL DEGRON SCORES ARE 0.0 -- structures are not parsing.")
        if all(t.lof_constraint_loeuf is None for t, _r in rows):
            problems.append("LOEUF IS EMPTY FOR EVERY ROW -- gnomAD returned nothing.")

    if STATS["genetics_errors"]:
        problems.append(
            f"{STATS['genetics_errors']} genetics lookups FAILED "
            f"(example: {STATS.get('genetics_error_example')})"
        )
    if STATS["no_ensembl"]:
        problems.append(
            f"{STATS['no_ensembl']} proteins had no Ensembl cross-reference, so "
            "genetics could not be looked up for them"
        )
    if STATS["no_efo"]:
        problems.append(f"{STATS['no_efo']} lookups skipped: disease could not be resolved to an EFO id")

    if problems:
        print("\n" + "!" * 68)
        print("THE SWEEP RAN BUT THE OUTPUT IS NOT USABLE")
        print("!" * 68)
        for p in problems:
            print(f"  - {p}")
        print("!" * 68)

    keep = (Verdict.TOP, Verdict.MID, Verdict.NOMINATE_FOR_VALIDATION)
    noms = sorted([x for x in rows if x[1].verdict in keep], key=lambda x: -x[1].composite)
    print(f"\n{seen} proteins seen, {len(rows)} scored, {len(noms)} nominated")
    print(f"written to {out}\n")
    for t, r in noms[:25]:
        print(f"  {r.composite:5.1f}  {t.gene:12s} {t.indication:30s} "
              f"V={r.validity.score:5.1f} M={r.modality.score:5.1f}  {r.verdict.value}")
    if not noms:
        print("  no nominations -- expected on a small --limit sample")
    n_cached, mb = cache.stats()
    print(f"\ncache: {n_cached} responses, {mb:.1f} MB")
    print("\nREMINDER: direction of effect is unresolved, so every genetics score")
    print("here is discounted to 40%. Ranks are conservative, not wrong.")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=None, help="stop after N proteins")
    p.add_argument("--area", choices=list(INDICATIONS), default=None)
    p.add_argument("--genes", default=None,
                   help="comma-separated gene symbols instead of the whole proteome")
    a = p.parse_args()
    sys.exit(run(a.limit, a.area, [g.strip() for g in a.genes.split(',')] if a.genes else None))
