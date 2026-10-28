"""Diagnose the genetics layer by asking Open Targets what it actually accepts.

STAT6, IRF5, PNPLA3 and PCSK9 all returned "none". Those have well-known
associations, so the query is wrong rather than the data being absent. Rather
than guess at argument names, this introspects the schema and prints the real
ones, then runs each step separately so the failing one is obvious.

    py -m tpd_bench.features.test_genetics
"""

from __future__ import annotations

import json

from . import http

OT = "https://api.platform.opentargets.org/api/v4/graphql"


def show(label: str, obj, limit: int = 600) -> None:
    txt = json.dumps(obj, indent=2) if not isinstance(obj, str) else obj
    print(f"  {label}:")
    for line in txt[:limit].splitlines():
        print(f"    {line}")
    if len(txt) > limit:
        print("    ...")


def step(name: str, fn):
    print(f"\n--- {name} ---")
    try:
        return fn()
    except Exception as exc:  # noqa: BLE001
        print(f"  FAILED: {type(exc).__name__}: {exc}")
        body = getattr(exc, "read", None)
        if body:
            try:
                print(f"  body: {body().decode()[:400]}")
            except Exception:
                pass
        return None


def run() -> int:
    print("=" * 70)
    print("OPEN TARGETS GENETICS DIAGNOSTIC")
    print("=" * 70)

    # 1. What arguments does Target.associatedDiseases actually take?
    def introspect():
        q = {"query": """
          query { __type(name:"Target"){ fields{
            name args{ name type{ name kind ofType{ name kind } } }
          }}}"""}
        d = http.request(OT, q)
        fields = (((d.get("data") or {}).get("__type") or {}).get("fields") or [])
        for f in fields:
            if f["name"] in ("associatedDiseases", "approvedSymbol", "id"):
                args = ", ".join(
                    f"{a['name']}: {a['type'].get('name') or a['type'].get('kind')}"
                    for a in f.get("args", [])
                ) or "(no args)"
                print(f"  Target.{f['name']}({args})")
        return fields
    step("1. Target field signatures", introspect)

    # 2. What does an Association row expose?
    def assoc_fields():
        q = {"query": '''query { __type(name:"AssociatedDisease"){ fields{ name } } }'''}
        d = http.request(OT, q)
        t = (d.get("data") or {}).get("__type")
        if not t:
            print("  no type named AssociatedDisease; trying DiseaseAssociation")
            q = {"query": '''query { __type(name:"DiseaseAssociation"){ fields{ name } } }'''}
            d = http.request(OT, q)
            t = (d.get("data") or {}).get("__type")
        print("  fields:", [f["name"] for f in (t or {}).get("fields", [])])
        return t
    step("2. Association row fields", assoc_fields)

    # 3. Does the disease search work, and what id comes back?
    def search():
        q = {
            "query": 'query($s:String!){search(queryString:$s,entityNames:["disease"])'
                     "{hits{id name entity}}}",
            "variables": {"s": "atopic dermatitis"},
        }
        d = http.request(OT, q)
        if d.get("errors"):
            show("errors", d["errors"])
        hits = (((d.get("data") or {}).get("search") or {}).get("hits") or [])[:3]
        for h in hits:
            print(f"  {h.get('id')}  {h.get('entity')}  {h.get('name')}")
        return hits[0]["id"] if hits else None
    efo = step("3. Disease search -> EFO id", search)

    # 4. The actual association query, unfiltered first
    def assoc_unfiltered():
        q = {
            "query": """query($g:String!){ target(ensemblId:$g){
                 approvedSymbol
                 associatedDiseases{ count rows{ score disease{ id name } } } } }""",
            "variables": {"g": "ENSG00000166888"},  # STAT6
        }
        d = http.request(OT, q)
        if d.get("errors"):
            show("errors", d["errors"])
            return None
        t = (d.get("data") or {}).get("target") or {}
        ad = t.get("associatedDiseases") or {}
        print(f"  symbol: {t.get('approvedSymbol')}  total associations: {ad.get('count')}")
        for r in (ad.get("rows") or [])[:5]:
            print(f"    {r['disease']['id']:16s} {r['score']:.3f}  {r['disease']['name']}")
        return ad
    step("4. STAT6 associations (unfiltered)", assoc_unfiltered)

    # 5. The filtered form the sweep uses -- this is the suspect
    def assoc_filtered():
        if not efo:
            print("  skipped: no EFO id from step 3")
            return None
        q = {
            "query": """query($g:String!,$d:String!){ target(ensemblId:$g){
                 associatedDiseases(efoIds:[$d]){ count rows{ score } } } }""",
            "variables": {"g": "ENSG00000166888", "d": efo},
        }
        d = http.request(OT, q)
        if d.get("errors"):
            show("errors", d["errors"])
            print("\n  ^ If this names a bad argument, that is the bug.")
            return None
        print("  filtered query succeeded:", d.get("data"))
        return d
    step("5. Filtered by EFO id (what the sweep does)", assoc_filtered)

    print("\n" + "=" * 70)
    print("Send me everything above.")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
