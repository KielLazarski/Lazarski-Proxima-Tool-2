"""Live data sources for the discordance corpus.

None of these reach anything from the Claude sandbox -- its network allowlist
covers package registries and GitHub only, so every call here returns an
explicit NetworkUnavailable rather than a silent empty list. An empty corpus
that looks like a successful run is worse than an error.

Run these in your own environment. The seed cache in `seed.py` is what the
bench uses offline.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from typing import Any

from .discordance import PhenotypeComparison, ResistanceEvidence

TIMEOUT = 30


class NetworkUnavailable(RuntimeError):
    """Raised instead of returning empty, so a blocked run cannot look clean."""


def _get(url: str) -> Any:
    try:
        with urllib.request.urlopen(url, timeout=TIMEOUT) as r:
            return json.load(r)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        raise NetworkUnavailable(
            f"could not reach {url.split('/')[2]} -- run this outside the sandbox"
        ) from exc


# --------------------------------------------------------------- tier 1

def catalytic_residues(uniprot_acc: str) -> set[int]:
    """Active-site and ATP-binding residue positions for a protein.

    A resistance mutation only counts as tier-1 evidence if it sits somewhere
    that plausibly ablates catalysis. Position alone is weak -- it narrows the
    candidate list for review, it does not settle the question.
    """
    data = _get(
        f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.json?fields=ft_act_site,ft_binding"
    )
    out: set[int] = set()
    for f in data.get("features", []):
        if f.get("type") in ("Active site", "Binding site"):
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value", start)
            if start:
                out.update(range(int(start), int(end) + 1))
    return out


def recurrent_mutations(gene: str, min_count: int = 3) -> list[dict]:
    """Recurrent somatic variants in a gene, from cBioPortal's public API.

    Recurrence under inhibitor pressure is the signal. A variant seen once is
    noise; a variant that keeps reappearing in treated patients is selection.
    """
    studies = _get("https://www.cbioportal.org/api/studies")
    ids = [s["studyId"] for s in studies][:50]
    hits: dict[str, int] = {}
    for sid in ids:
        try:
            muts = _get(
                f"https://www.cbioportal.org/api/molecular-profiles/{sid}_mutations"
                f"/mutations?sampleListId={sid}_all&projection=DETAILED"
            )
        except NetworkUnavailable:
            raise
        except Exception:
            continue
        for m in muts:
            if m.get("gene", {}).get("hugoGeneSymbol") == gene:
                p = m.get("proteinChange")
                if p:
                    hits[p] = hits.get(p, 0) + 1
    return [{"variant": v, "count": c} for v, c in hits.items() if c >= min_count]


def mine_resistance(gene: str, uniprot_acc: str) -> list[ResistanceEvidence]:
    """Candidates only. `signalling_persists` cannot be mined and is set False.

    That field is the whole tier-1 claim and it comes from functional work, not
    a database. Everything this returns needs a human to confirm it before the
    bench treats it as tier 1 -- which is why nothing here auto-promotes.
    """
    sites = catalytic_residues(uniprot_acc)
    out = []
    for m in recurrent_mutations(gene):
        pos = "".join(c for c in m["variant"] if c.isdigit())
        if pos and int(pos) in sites:
            out.append(
                ResistanceEvidence(
                    gene=gene,
                    variant=m["variant"],
                    catalytic_impact="impairs",
                    signalling_persists=False,  # requires manual confirmation
                    source=f"cBioPortal recurrent (n={m['count']}); UniProt {uniprot_acc}",
                )
            )
    return out


# --------------------------------------------------------------- tier 2/3

def impc_phenotypes(gene: str) -> list[dict]:
    """Mouse knockout phenotypes from IMPC.

    Note that IMPC knockouts are overwhelmingly constitutive, so anything built
    from them alone fails the developmental-compensation control and is capped
    at tier 3. Conditional data has to come from MGI or the literature.
    """
    data = _get(
        "https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select"
        f"?q=marker_symbol:{gene}&rows=200&wt=json"
    )
    return data.get("response", {}).get("docs", [])


def build_comparisons(gene: str) -> list[PhenotypeComparison]:
    """Skeleton. The inhibitor arm has no public structured source.

    This is the real obstacle to automating tier 2: knockout phenotypes are in
    a database, inhibitor phenotypes are in prose. Until the inhibitor side is
    curated, the join is half a join.
    """
    raise NotImplementedError(
        "inhibitor-arm phenotypes are not available from any structured source; "
        "curate them or extract from literature before enabling this path"
    )
