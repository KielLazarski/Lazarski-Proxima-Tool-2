"""Build the viewer: one self-contained HTML file, no server, no network.

    python3 -m tpd_bench.viewer

Data is generated from the live bench at build time rather than hand-copied,
so the page cannot drift out of step with the scoring code. The bundled demo
run in demo/ is embedded as raw CSV so the page opens on real output; dropping
a sweep CSV on it replaces that through the same parser.
"""

from __future__ import annotations

import csv
import io
import json
import pathlib
from contextlib import redirect_stdout

from .ablate import ABLATIONS
from .glossary import GLOSSARY
from .corpus.apply import apply as apply_corpus
from .scoring import PROFILES, Verdict, score
from .schema import Track
from .targets import TARGETS as RAW

# index.html at the repo root, not a path relative to the current directory:
# the build must land in the same place whatever folder it is run from.
#
# It is called index.html, and it is committed rather than ignored, so that
# GitHub Pages serves the tool itself. Pages looks for index.html at the site
# root and falls back to rendering README.md when there isn't one -- which is
# what a visitor got while this was a build artifact.
#
# The cost is a generated file under version control, which can go stale. Run
# OPEN_VIEWER.py (or python3 -m tpd_bench.viewer) after changing anything the
# page displays, and commit the result; git status will show it.
OUT = pathlib.Path(__file__).resolve().parent.parent / "index.html"


def scored_rows() -> list[dict]:
    rows = []
    labels = {}
    lp = pathlib.Path(__file__).with_name("calibration_labels.csv")
    with lp.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            labels[r["target_id"]] = (r["expected_verdict"], r["rationale"])
    for t in apply_corpus(RAW, verbose=False):
        r = score(t)
        exp, why = labels.get(t.target_id, ("", ""))
        rows.append({
            "id": t.target_id, "gene": t.gene, "indication": t.indication,
            "area": t.disease_area, "track": r.track.value,
            "compartment": t.subcellular.value,
            "v": r.validity.score, "m": r.modality.score, "c": r.composite,
            "verdict": r.verdict.value, "expected": exp, "why": why,
            "agree": r.verdict.value == exp,
            "vflags": r.validity.flags, "mflags": r.modality.flags,
            "flags": r.flags, "notes": t.notes,
            "genetics": t.genetics.value,
            "degron": t.degron_feature_score,
            "scaffold": t.scaffold_evidence_tier,
            "sign": t.scaffold_sign.value,
            "terms": r.validity.terms,
            # Hard gates do not depend on the validity weights, so a slider
            # must not be able to overturn them.
            "gated": r.verdict.value not in (
                "top", "mid", "nominate_for_validation", "reject_target_invalid"),
            "refuted": "validity_refuted" in r.validity.reasons,
            # Applied after the weighted mean: a degenerate target class is
            # penalised whatever the weights say.
            "vmul": 1.0 - 0.5 * t.network_redundancy,
        })
    return rows


def profiles() -> dict:
    return {k: {"w_genetics": p.w_genetics, "w_pathway": p.w_pathway_drug,
                "w_clinical": p.w_clinical, "w_expression": p.w_expression,
                "top": p.top_threshold, "mid": p.mid_threshold,
                "strong_modality": p.strong_modality,
                "validation_floor": p.validation_floor,
                "modality_floor": p.modality_floor}
            for k, p in PROFILES.items()}


def glossary_patterns() -> list:
    """The glossary's own pattern table, longest pattern first.

    Embedding patterns rather than a dict keyed on the exact flag strings seen
    at build time means the page can explain a flag it has never seen -- which
    matters, because flags carry numbers. 'LoF-tolerant in humans (LOEUF 1.26)'
    is one key per gene, and a user's own sweep will produce values no build
    could have enumerated. The page mirrors explain(): longest match wins.
    """
    return [[k, v] for k, v in sorted(GLOSSARY.items(), key=lambda kv: -len(kv[0]))]


def ablation_rows() -> list[dict]:
    import copy
    import dataclasses
    tg = apply_corpus(RAW, verbose=False)
    base = {t.target_id: score(t).verdict for t in tg}
    out = []
    for name, patch in ABLATIONS.items():
        moved = []
        for t in tg:
            t2 = dataclasses.replace(copy.deepcopy(t), **patch)
            if score(t2).verdict is not base[t.target_id]:
                moved.append(t.gene)
        out.append({"term": name, "n": len(moved), "genes": sorted(set(moved))})
    return sorted(out, key=lambda r: -r["n"])


def demo_sweep() -> dict:
    """The bundled sweep run, embedded so the page opens with real output.

    The raw CSV text is embedded rather than parsed rows, so the demo goes
    through exactly the same parser as a file the user drops in. One code
    path, so the demo cannot render differently from a real run.
    """
    p = pathlib.Path(__file__).resolve().parent.parent / "demo" / "sweep_results_demo.csv"
    if not p.exists():
        return {}
    return {
        "name": p.name,
        "text": p.read_text(encoding="utf-8"),
        "when": "28 October 2026",
        "scope": "6 genes against 4 indications",
    }


def corpus_rows() -> list[dict]:
    p = pathlib.Path(__file__).parent / "corpus" / "scaffold_evidence.csv"
    if not p.exists():
        from .corpus.build import run as build
        with redirect_stdout(io.StringIO()):
            build()
    with p.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Degradation target bench</title>
<style>
:root{
  --ink:#16202B; --paper:#F7F8F6; --panel:#FFFFFF; --rule:#D3D8D4;
  --dim:#5E6B74; --faint:#8A94A0;
  --take:#1F6F5C; --validate:#B07A12; --reject:#8A94A0; --flagc:#6B4FA8;
  --link:#0E6E8C; --grid:#E4E8E4;
  --ui:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  --num:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.5 var(--ui);
  -webkit-font-smoothing:antialiased}
header{padding:22px 28px 0;max-width:1480px;margin:0 auto}
h1{font-size:23px;font-weight:650;letter-spacing:-.012em;margin:0 0 3px}
.sub{color:var(--dim);font-size:14px;margin:0 0 18px;max-width:72ch}
main{max-width:1480px;margin:0 auto;padding:0 28px 60px}

nav{display:flex;gap:2px;border-bottom:1px solid var(--rule);margin-bottom:20px}
nav button{appearance:none;border:0;background:none;font:inherit;color:var(--dim);
  padding:9px 15px;cursor:pointer;border-bottom:2px solid transparent;margin-bottom:-1px}
nav button:hover{color:var(--ink)}
nav button[aria-selected=true]{color:var(--ink);border-bottom-color:var(--ink);font-weight:600}
nav button:focus-visible{outline:2px solid var(--link);outline-offset:-2px}

.split{display:grid;grid-template-columns:minmax(380px,460px) 1fr;gap:26px;align-items:start}
@media(max-width:1020px){.split{grid-template-columns:1fr}}

.plotwrap{background:var(--panel);border:1px solid var(--rule);padding:16px 16px 10px}
.plotcap{font-size:13px;color:var(--dim);margin:10px 2px 0;line-height:1.45}
svg{display:block;width:100%;height:auto}
.pt{cursor:pointer}
.pt circle{transition:r .12s}
.pt:hover circle{r:7}
.pt.on circle{stroke:var(--ink);stroke-width:2.5}

.controls{display:flex;flex-wrap:wrap;gap:8px;margin:0 0 12px}
select,input[type=search]{font:inherit;padding:6px 9px;border:1px solid var(--rule);
  background:var(--panel);color:var(--ink);border-radius:0}
input[type=search]{min-width:190px}
select:focus-visible,input:focus-visible{outline:2px solid var(--link);outline-offset:1px}

table{width:100%;border-collapse:collapse;font-size:14px;background:var(--panel)}
th{text-align:left;font-weight:600;font-size:12.5px;color:var(--dim);
  border-bottom:1px solid var(--rule);padding:8px 10px;position:sticky;top:0;background:var(--panel);
  cursor:pointer;white-space:nowrap}
th:hover{color:var(--ink)}
td{padding:7px 10px;border-bottom:1px solid var(--grid);vertical-align:top}
tbody tr{cursor:pointer}
tbody tr:hover{background:#F0F3F0}
tbody tr.on{background:#E8EFEC;box-shadow:inset 3px 0 0 var(--ink)}
.n{font-family:var(--num);font-size:13px;text-align:right;font-variant-numeric:tabular-nums}
.gene{font-weight:600}
.scroll{max-height:60vh;overflow:auto;border:1px solid var(--rule)}

.v{display:inline-block;font-size:12px;padding:1px 7px;border:1px solid;border-radius:999px;white-space:nowrap}
.v.take{color:var(--take);border-color:var(--take)}
.v.validate{color:var(--validate);border-color:var(--validate)}
.v.reject{color:var(--reject);border-color:var(--reject)}
.v.flag{color:var(--flagc);border-color:var(--flagc)}

.detail{background:var(--panel);border:1px solid var(--rule);padding:18px 20px;margin-top:18px}
.detail h3{margin:0 0 2px;font-size:19px}
.detail .where{color:var(--dim);font-size:13.5px;margin:0 0 16px}
.axes{display:grid;grid-template-columns:1fr 1fr;gap:20px}
@media(max-width:700px){.axes{grid-template-columns:1fr}}
.axis h4{margin:0 0 8px;font-size:13px;color:var(--dim);font-weight:600}
.bar{height:7px;background:var(--grid);margin:0 0 10px;position:relative}
.bar i{display:block;height:100%;background:var(--ink)}
.score{font-family:var(--num);font-size:27px;line-height:1;font-weight:500}
.score small{font-size:13px;color:var(--faint);font-family:var(--ui);margin-left:3px}
.why{list-style:none;padding:0;margin:10px 0 0;font-size:13.5px;color:var(--dim)}
.why li{padding:3px 0 3px 15px;position:relative;line-height:1.45}
.why li:before{content:"";position:absolute;left:0;top:11px;width:7px;height:1px;background:var(--faint)}
.why li.warn{color:#8A5A00}
.disagree{margin-top:16px;padding:11px 13px;background:#FBF6E9;border-left:3px solid var(--validate);font-size:13.5px}

.drop{border:1.5px dashed var(--rule);padding:26px;text-align:center;color:var(--dim);
  background:var(--panel);margin-bottom:18px}
.drop.over{border-color:var(--link);color:var(--link);background:#F2F9FB}
.drop b{color:var(--ink)}
.empty{padding:40px 10px;text-align:center;color:var(--faint)}
.ab{display:flex;align-items:center;gap:11px;padding:5px 0;font-size:14px}
.ab .nm{flex:0 0 220px}
.ab .track{flex:1;height:11px;background:var(--grid)}
.ab .track i{display:block;height:100%;background:var(--ink)}
.ab.zero .nm,.ab.zero .cnt{color:var(--faint)}
.ab .cnt{font-family:var(--num);font-size:13px;width:46px;text-align:right}
.ab .who{flex:0 0 100%;font-size:12.5px;color:var(--faint);padding-left:231px;margin-top:-3px}
.note{font-size:13.5px;color:var(--dim);max-width:76ch;line-height:1.55}
.note strong{color:var(--ink);font-weight:600}

#tip{position:fixed;pointer-events:none;z-index:50;background:var(--ink);color:#fff;
  padding:7px 10px;font-size:12.5px;line-height:1.4;max-width:240px;opacity:0;transition:opacity .1s}
#tip.on{opacity:1}
#tip b{font-weight:650}
#tip .co{font-family:var(--num);color:#B8C4C0;font-size:11.5px}
.ptools{display:flex;gap:6px;align-items:center;margin:8px 2px 0;font-size:12.5px;color:var(--dim)}
.ptools button{font:inherit;border:1px solid var(--rule);background:var(--panel);
  padding:3px 9px;cursor:pointer;color:var(--ink)}
.ptools button:hover{border-color:var(--ink)}
.plot{cursor:grab;touch-action:none}
.plot.drag{cursor:grabbing}
.wts{background:var(--panel);border:1px solid var(--rule);padding:14px 16px;margin-bottom:14px}
.wts summary{cursor:pointer;font-size:13.5px;font-weight:600;list-style:none}
.wts summary::-webkit-details-marker{display:none}
.wts summary:before{content:"\25b8 ";color:var(--faint)}
.wts[open] summary:before{content:"\25be "}
.wrow{display:grid;grid-template-columns:150px 1fr 48px;gap:10px;align-items:center;margin-top:10px;font-size:13.5px}
.wrow input{width:100%}
.wrow .wv{font-family:var(--num);font-size:12.5px;text-align:right;color:var(--dim)}
.wts .foot{margin-top:13px;font-size:12.5px;color:var(--faint);line-height:1.5}
.wts .foot button{font:inherit;border:1px solid var(--rule);background:var(--panel);
  padding:2px 9px;cursor:pointer;color:var(--ink);margin-right:7px}
.expl{display:flex;gap:0;border:1px solid var(--rule);margin-left:auto}
.expl button{font:inherit;font-size:12.5px;border:0;background:var(--panel);color:var(--dim);
  padding:4px 11px;cursor:pointer}
.expl button[aria-pressed=true]{background:var(--ink);color:#fff}
.dhead{display:flex;align-items:baseline;gap:12px;flex-wrap:wrap}
.why li .src{display:block;font-size:12px;color:var(--faint);margin-top:2px}
.why li .src a{color:var(--link)}
.why li.det{padding-bottom:9px}
@media(prefers-reduced-motion:reduce){*{transition:none!important}}
</style></head><body>
<div id="tip" role="status" aria-live="polite"></div>
<header>
  <h1>Degradation target bench</h1>
  <p class="sub">Every candidate is judged twice &mdash; whether the target is real, and whether
  degradation is the right tool. The two scores are kept apart, because a real target you
  cannot reach and a reachable target that does not drive disease fail for different reasons.</p>
</header>
<main>
  <nav role="tablist">
    <button role="tab" data-tab="noms" aria-selected="true">Nominations</button>
    <button role="tab" data-tab="cal" aria-selected="false">Calibration set</button>
    <button role="tab" data-tab="corp" aria-selected="false">Scaffolding corpus</button>
    <button role="tab" data-tab="diag" aria-selected="false">What moves verdicts</button>
  </nav>
  <section id="noms"></section>
  <section id="cal" hidden></section>
  <section id="corp" hidden></section>
  <section id="diag" hidden></section>
</main>
<script>
const CAL = __CAL__, ABL = __ABL__, CORP = __CORP__, PROF = __PROF__, GPAT = __GLOSS__;
const DEMO = __DEMO__;

/* Mirrors explain() in glossary.py: longest matching pattern wins, so a
   specific flag beats a general one. Patterns, not exact flags, because flags
   carry numbers -- a sweep nobody has run yet still gets a cited explanation. */
const _gcache = new Map();
function gloss(f){
  if(_gcache.has(f)) return _gcache.get(f);
  const lf = f.toLowerCase();
  let hit = null;
  for(const [k,v] of GPAT){ if(lf.includes(k.toLowerCase())){ hit = v; break; } }
  _gcache.set(f, hit);
  return hit;
}
const $ = (s,r)=> (r||document).querySelector(s);
const el = (t,a,kids)=>{const n=document.createElement(t);
  for(const k in (a||{})) k==='html'?n.innerHTML=a[k]:(k==='text'?n.textContent=a[k]:n.setAttribute(k,a[k]));
  (kids||[]).forEach(c=>n.appendChild(c)); return n;};

const FAMILY = v => v==='top'||v==='mid' ? 'take'
  : v==='nominate_for_validation' ? 'validate'
  : v.startsWith('flag') ? 'flag' : 'reject';
const LABEL = {top:'Take forward', mid:'Worth a look',
  nominate_for_validation:'Needs an experiment',
  reject_wrong_modality:'Cannot reach it', reject_modality_too_weak:'Weak modality case',
  reject_inhibition_sufficient:'Inhibition suffices', reject_scaffold_protective:'Scaffold is protective',
  reject_target_invalid:'Target not supported', flag_contested:'Contested',
  flag_strategic_discontinuation:'Dropped for portfolio reasons'};
const lab = v => LABEL[v] || v.replace(/_/g,' ');

/* ---------- the scatter: the whole design in one picture ---------- */
const TIP = () => document.getElementById('tip');
function showTip(e, r){
  const t = TIP();
  t.innerHTML = `<b>${r.gene}</b><br>${r.indication}<br>`
    + `<span class="co">relevance ${r.v.toFixed(1)} · modality ${r.m.toFixed(1)}`
    + `${isFinite(r.c)?' · combined '+r.c.toFixed(1):''}</span><br>${lab(r.verdict)}`;
  t.classList.add('on');
  const pad = 14, w = t.offsetWidth, h = t.offsetHeight;
  let x = e.clientX + pad, y = e.clientY + pad;
  if (x + w > innerWidth - 8) x = e.clientX - w - pad;
  if (y + h > innerHeight - 8) y = e.clientY - h - pad;
  t.style.left = x + 'px'; t.style.top = y + 'px';
}
const hideTip = () => TIP().classList.remove('on');

function plot(rows, sel, onPick, zoom){
  const W=440,H=400,P={l:48,r:16,t:26,b:42}, iw=W-P.l-P.r, ih=H-P.t-P.b;
  const X=v=>P.l+(v/100)*iw, Y=v=>P.t+ih-(v/100)*ih;
  const ns='http://www.w3.org/2000/svg';
  const mk=(t,a,kids)=>{const n=document.createElementNS(ns,t);
    for(const k in a) n.setAttribute(k,a[k]); (kids||[]).forEach(c=>n.appendChild(c)); return n;};
  const svg=mk('svg',{viewBox:`0 0 ${W} ${H}`,class:'plot',role:'img',
    'aria-label':'Therapeutic relevance against modality fit'});

  /* axis furniture stays outside the zoom group so labels keep their size */
  const t=(x,y,str,o)=>{const n=mk('text',Object.assign({x,y,fill:'#8A94A0','font-size':10.5},o||{}));
    n.textContent=str; return n;};

  const clip=mk('clipPath',{id:'pclip'},[mk('rect',{x:P.l,y:P.t,width:iw,height:ih})]);
  svg.appendChild(mk('defs',{},[clip]));
  const zg=mk('g',{'clip-path':'url(#pclip)'});
  const g=mk('g',{transform:`translate(${zoom.tx} ${zoom.ty}) scale(${zoom.k})`});
  zg.appendChild(g); svg.appendChild(zg);

  g.appendChild(mk('rect',{x:X(0),y:Y(100),width:X(55)-X(0),height:Y(35)-Y(100),fill:'#FBF6E9'}));
  g.appendChild(mk('rect',{x:X(55),y:Y(35),width:X(100)-X(55),height:Y(0)-Y(35),fill:'#FBF6E9'}));
  g.appendChild(mk('rect',{x:X(55),y:Y(100),width:X(100)-X(55),height:Y(35)-Y(100),fill:'#EBF3EF'}));
  for(const v of [0,25,50,75,100]){
    g.appendChild(mk('line',{x1:X(v),y1:P.t,x2:X(v),y2:P.t+ih,stroke:'#E4E8E4','vector-effect':'non-scaling-stroke'}));
    g.appendChild(mk('line',{x1:P.l,y1:Y(v),x2:P.l+iw,y2:Y(v),stroke:'#E4E8E4','vector-effect':'non-scaling-stroke'}));
  }
  g.appendChild(mk('line',{x1:X(55),y1:P.t,x2:X(55),y2:P.t+ih,stroke:'#B6BEB8',
    'stroke-dasharray':'3 3','vector-effect':'non-scaling-stroke'}));
  g.appendChild(mk('line',{x1:P.l,y1:Y(35),x2:P.l+iw,y2:Y(35),stroke:'#B6BEB8',
    'stroke-dasharray':'3 3','vector-effect':'non-scaling-stroke'}));
  if(zoom.k < 1.6){
    g.appendChild(t(X(27),Y(97)+4,'relevant, but out of reach',{'text-anchor':'middle',fill:'#9C7B2E'}));
    g.appendChild(t(X(78),Y(5),'degradable, not driving disease',{'text-anchor':'middle',fill:'#9C7B2E'}));
    g.appendChild(t(X(78),Y(97)+4,'take forward',{'text-anchor':'middle',fill:'#2F7A65'}));
  }

  const C={take:'#1F6F5C',validate:'#B07A12',reject:'#8A94A0',flag:'#6B4FA8'};
  const rr = 5/Math.sqrt(zoom.k);
  rows.forEach(r=>{
    const f=FAMILY(r.verdict);
    const on = sel&&sel.id===r.id&&sel.indication===r.indication;
    const gp=mk('g',{class:'pt'+(on?' on':''),tabindex:'0',role:'button',
      'aria-label':`${r.gene}, ${r.indication}, relevance ${r.v.toFixed(0)}, modality ${r.m.toFixed(0)}, ${lab(r.verdict)}`});
    gp.appendChild(mk('circle',{cx:X(r.m),cy:Y(r.v),r:rr,fill:C[f],
      'fill-opacity':f==='reject'?.42:.85,stroke:'#fff',
      'stroke-width':1/zoom.k,'vector-effect':'non-scaling-stroke'}));
    gp.addEventListener('pointerenter',e=>showTip(e,r));
    gp.addEventListener('pointermove',e=>showTip(e,r));
    gp.addEventListener('pointerleave',hideTip);
    gp.addEventListener('click',e=>{e.stopPropagation();showTip(e,r);onPick(r);});
    gp.addEventListener('focus',()=>{const b=gp.getBoundingClientRect();
      showTip({clientX:b.left+b.width/2,clientY:b.top},r);});
    gp.addEventListener('blur',hideTip);
    gp.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();onPick(r);}});
    g.appendChild(gp);
  });

  svg.appendChild(t(P.l-8,Y(0)+4,'0',{'text-anchor':'end'}));
  svg.appendChild(t(P.l-8,Y(100)+4,'100',{'text-anchor':'end'}));
  svg.appendChild(t(P.l,H-10,'0'));
  svg.appendChild(t(P.l+iw,H-10,'100',{'text-anchor':'end'}));
  svg.appendChild(t(P.l+iw/2,H-10,'Is degradation the right tool?',{'text-anchor':'middle',fill:'#5E6B74','font-size':11.5}));
  const yl=mk('text',{x:14,y:P.t+ih/2,fill:'#5E6B74','font-size':11.5,'text-anchor':'middle',
    transform:`rotate(-90 14 ${P.t+ih/2})`}); yl.textContent='Is the target therapeutically relevant?';
  svg.appendChild(yl);

  /* wheel to zoom about the cursor, drag to pan */
  svg.addEventListener('wheel',e=>{
    e.preventDefault();
    const pt=svg.createSVGPoint(); pt.x=e.clientX; pt.y=e.clientY;
    const loc=pt.matrixTransform(svg.getScreenCTM().inverse());
    const f=e.deltaY<0?1.18:1/1.18;
    const k=Math.max(1,Math.min(14,zoom.k*f));
    const real=k/zoom.k;
    zoom.tx = loc.x-(loc.x-zoom.tx)*real;
    zoom.ty = loc.y-(loc.y-zoom.ty)*real;
    zoom.k = k;
    if(k===1){zoom.tx=0;zoom.ty=0;}
    zoom.redraw();
  },{passive:false});
  let drag=null;
  svg.addEventListener('pointerdown',e=>{
    if(e.target.closest('.pt')) return;
    drag={x:e.clientX,y:e.clientY,tx:zoom.tx,ty:zoom.ty};
    svg.classList.add('drag'); svg.setPointerCapture(e.pointerId);
  });
  svg.addEventListener('pointermove',e=>{
    if(!drag) return;
    const s=svg.getBoundingClientRect().width/W;
    zoom.tx=drag.tx+(e.clientX-drag.x)/s; zoom.ty=drag.ty+(e.clientY-drag.y)/s;
    zoom.redraw();
  });
  const stop=()=>{drag=null;svg.classList.remove('drag');};
  svg.addEventListener('pointerup',stop); svg.addEventListener('pointercancel',stop);
  svg.addEventListener('pointerleave',()=>{stop();hideTip();});
  return svg;
}

function detail(r, mode, onMode){
  const bar=(v)=>el('div',{class:'bar'},[el('i',{style:`width:${Math.max(0,Math.min(100,v))}%`})]);
  const line=(f)=>{
    const g=gloss(f);
    const warn=/unknown|provisional|renormalis|not looked up|rests on/.test(f);
    if(mode==='simple'&&g) return el('li',{text:g.simple,class:warn?'warn':''});
    const li=el('li',{class:(warn?'warn ':'')+'det'});
    li.appendChild(document.createTextNode(g?g.detail:f));
    if(g&&g.source){
      const sp=el('span',{class:'src'});
      sp.appendChild(document.createTextNode('Source: '));
      if(g.url){const a=el('a',{href:g.url,target:'_blank',rel:'noopener',text:g.source});sp.appendChild(a);}
      else sp.appendChild(document.createTextNode(g.source));
      li.appendChild(sp);
    }
    return li;
  };
  const list=(fs)=>el('ul',{class:'why'},fs.map(line));

  const toggle=el('div',{class:'expl',role:'group','aria-label':'Explanation detail'});
  [['simple','Plain'],['detail','Detailed']].forEach(([k,t])=>{
    const b=el('button',{text:t,'aria-pressed':String(mode===k)});
    b.addEventListener('click',()=>onMode(k)); toggle.appendChild(b);
  });

  const box=el('div',{class:'detail'},[
    el('div',{class:'dhead'},[el('h3',{text:r.gene}),toggle]),
    el('p',{class:'where',text:`${r.indication} — ${r.compartment}, ${r.track} track`}),
    el('div',{class:'axes'},[
      el('div',{class:'axis'},[el('h4',{text:'Is the target therapeutically relevant?'}),
        el('div',{class:'score',html:r.v.toFixed(1)+'<small>/100</small>'}),bar(r.v),list(r.vflags)]),
      el('div',{class:'axis'},[el('h4',{text:'Is degradation the right tool?'}),
        el('div',{class:'score',html:r.m.toFixed(1)+'<small>/100</small>'}),bar(r.m),list(r.mflags)]),
    ])
  ]);
  if(r.flags&&r.flags.length) box.appendChild(list(r.flags));
  if(r.notes) box.appendChild(el('p',{class:'note',style:'margin-top:16px',text:r.notes}));
  if(r.expected && !r.agree) box.appendChild(el('div',{class:'disagree',
    html:`The bench said <b>${lab(r.verdict)}</b>; the calibration set expects
      <b>${lab(r.expected)}</b>. ${r.why}`}));
  return box;
}

/* ---------- reweighting ---------- */
const WDEF = {w_genetics:'Human genetics', w_pathway:'Approved drug on the pathway',
  w_clinical:'Clinical outcome', w_expression:'Disease expression'};
const TDEF = {top:'Take forward above', mid:'Worth a look above',
  strong_modality:'Strong modality above', modality_floor:'Reject modality below'};

function rescore(r, w, th, touched){
  if(!r.terms) return r;                       /* a run without term columns */
  /* Until a slider is moved, each row uses the weight profile for its own
     disease area -- cardiometabolic leans on genetics, immunology does not.
     Moving any slider overrides every area at once, which is the point of
     asking "what if we weighted this differently". */
  if(!touched){
    const pr = PROF[r.area];
    if(pr) w = {w_genetics:pr.w_genetics,w_pathway:pr.w_pathway,
                w_clinical:pr.w_clinical,w_expression:pr.w_expression};
  }
  const parts=[['genetics',w.w_genetics],['pathway',w.w_pathway],
    ['expression',w.w_expression],['clinical',w.w_clinical]]
    .filter(([k])=>r.terms[k]!==null&&r.terms[k]!==undefined);
  const tot=parts.reduce((a,[k,wt])=>a+wt,0);
  let v=tot?parts.reduce((a,[k,wt])=>a+wt*r.terms[k],0)/tot:0;
  v *= (r.vmul===undefined?1:r.vmul);
  /* A trial that showed harm zeroes relevance regardless of weighting, and
     a refuted causal case does the same. No slider should soften that. */
  if(r.refuted) v=0;
  const c=Math.sqrt(Math.max(0,v)*Math.max(0,r.m));
  let verdict=r.verdict;
  /* Hard gates and refutations do not depend on these weights, so a slider
     must not be able to overturn them. */
  if(!r.gated && !r.refuted){
    verdict = v<th.mid
      ? (r.m>=th.strong_modality&&v>=th.validation_floor?'nominate_for_validation':'reject_target_invalid')
      : (c>=th.top?'top':'mid');
  }
  return Object.assign({},r,{v:Math.round(v*10)/10,c:Math.round(c*10)/10,verdict});
}

function weightPanel(w, th, base, onChange, canReweight, touched){
  const d=el('details',{class:'wts'});
  d.appendChild(el('summary',{text:'Weights and thresholds'}));
  if(!canReweight){
    d.appendChild(el('p',{class:'foot',text:
      'This run has no per-term columns, so its scores cannot be reweighted here. '
      +'Re-run the sweep with the current version to enable the sliders.'}));
    return d;
  }
  if(!touched) d.appendChild(el('p',{class:'foot',style:'margin:8px 0 0',text:
    'Each disease area currently uses its own weights. Moving any slider applies '
    + 'one set of weights to everything.'}));
  const row=(label,val,min,max,step,set)=>{
    const out=el('span',{class:'wv',text:(+val).toFixed(2)});
    const i=el('input',{type:'range',min,max,step,value:val,'aria-label':label});
    i.addEventListener('input',e=>{out.textContent=(+e.target.value).toFixed(2);set(+e.target.value);});
    return el('div',{class:'wrow'},[el('label',{text:label}),i,out]);
  };
  Object.entries(WDEF).forEach(([k,label])=>
    d.appendChild(row(label,w[k],0,1,0.05,v=>{w[k]=v;onChange();})));
  Object.entries(TDEF).forEach(([k,label])=>
    d.appendChild(row(label,th[k],0,100,1,v=>{th[k]=v;onChange();})));
  const foot=el('div',{class:'foot'});
  const reset=el('button',{text:'Reset to defaults'});
  reset.addEventListener('click',()=>{Object.assign(w,base.w);Object.assign(th,base.th);onChange(true);});
  d.setAttribute('open','');
  foot.appendChild(reset);
  foot.appendChild(document.createTextNode(
    'Weights change how the relevance score is combined, not what the evidence says. '
    + 'Hard rejections — wrong compartment, protective scaffold, an inhibitor that already '
    + 'works — are not affected, because no weighting makes those reasons go away.'));
  d.appendChild(foot);
  return d;
}

function tableOf(rows, sort, onSort, sel, onPick){
  const cols=[['gene','Gene'],['indication','Indication'],['c','Score'],
    ['v','Target'],['m','Modality'],['verdict','Verdict']];
  const thead=el('tr',{},cols.map(([k,t])=>el('th',{text:t+(sort.k===k?(sort.d>0?' ▲':' ▼'):''),
    class:['c','v','m'].includes(k)?'n':''})));
  thead.querySelectorAll('th').forEach((th,i)=>th.addEventListener('click',()=>onSort(cols[i][0])));
  const body=el('tbody',{},rows.map(r=>{
    const tr=el('tr',{class:(sel&&sel.id===r.id&&sel.indication===r.indication)?'on':''},[
      el('td',{text:r.gene,class:'gene'}),
      el('td',{text:r.indication}),
      el('td',{text:r.c.toFixed(1),class:'n'}),
      el('td',{text:r.v.toFixed(0),class:'n'}),
      el('td',{text:r.m.toFixed(0),class:'n'}),
      el('td',{},[el('span',{class:'v '+FAMILY(r.verdict),text:lab(r.verdict)})]),
    ]);
    tr.addEventListener('click',()=>onPick(r)); return tr;}));
  return el('div',{class:'scroll'},[el('table',{},[el('thead',{},[thead]),body])]);
}

/* ---------- views ---------- */
function browser(host, rows, caption, profileKey){
  const base0 = PROF[profileKey] || PROF.cardiometabolic;
  const base = {w:{w_genetics:base0.w_genetics,w_pathway:base0.w_pathway,
                   w_clinical:base0.w_clinical,w_expression:base0.w_expression},
                th:{top:base0.top,mid:base0.mid,strong_modality:base0.strong_modality,
                    modality_floor:base0.modality_floor,validation_floor:base0.validation_floor}};
  const w=Object.assign({},base.w), th=Object.assign({},base.th);
  const canReweight = rows.some(r=>r.terms);
  let sel=null, sort={k:'c',d:-1}, area='', verdict='', q='', mode='simple', touched=false;
  const zoom={k:1,tx:0,ty:0,redraw:()=>view(true)};

  function view(keepScroll){
    const scored = canReweight ? rows.map(r=>rescore(r,w,th,touched)) : rows;
    const f=scored.filter(r=>(!area||r.area===area)&&(!verdict||FAMILY(r.verdict)===verdict)
      &&(!q||(r.gene+' '+r.indication).toLowerCase().includes(q)));
    f.sort((a,b)=>{const x=a[sort.k],y=b[sort.k];
      return (typeof x==='number'?x-y:String(x).localeCompare(String(y)))*sort.d;});
    if(sel){const m=f.find(r=>r.id===sel.id&&r.indication===sel.indication); if(m) sel=m;}
    const y=window.scrollY;
    host.textContent='';

    host.appendChild(weightPanel(w,th,base,(reset)=>{touched=!reset;view(true);},canReweight,touched));

    const areas=[...new Set(rows.map(r=>r.area))].sort();
    const ctl=el('div',{class:'controls'});
    const mkSel=(opts,cur,on,none)=>{const sl=el('select');
      sl.appendChild(el('option',{value:'',text:none}));
      opts.forEach(o=>sl.appendChild(el('option',{value:o[0],text:o[1],...(cur===o[0]?{selected:''}:{})})));
      sl.addEventListener('change',e=>on(e.target.value)); return sl;};
    ctl.appendChild(mkSel(areas.map(a=>[a,a.replace(/_/g,' ')]),area,v=>{area=v;view();},'All disease areas'));
    ctl.appendChild(mkSel([['take','Take forward'],['validate','Needs an experiment'],
      ['flag','Flagged'],['reject','Rejected']],verdict,v=>{verdict=v;view();},'All verdicts'));
    const sq=el('input',{type:'search',placeholder:'Find a gene',value:q});
    sq.addEventListener('input',e=>{q=e.target.value.toLowerCase();view(true);});
    ctl.appendChild(sq);
    ctl.appendChild(el('span',{style:'align-self:center;color:var(--dim);font-size:13.5px',
      text:`${f.length} of ${rows.length}`}));

    const pick=r=>{sel=r;view(true);
      const d=$('#detailhost'); if(d) d.scrollIntoView({block:'nearest',behavior:'smooth'});};
    const left=el('div',{class:'plotwrap'},[plot(f,sel,pick,zoom)]);
    const tools=el('div',{class:'ptools'});
    const zb=(t,fn)=>{const b=el('button',{text:t});b.addEventListener('click',fn);return b;};
    tools.appendChild(zb('−',()=>{zoom.k=Math.max(1,zoom.k/1.4);if(zoom.k===1){zoom.tx=0;zoom.ty=0;}view(true);}));
    tools.appendChild(zb('+',()=>{zoom.k=Math.min(14,zoom.k*1.4);view(true);}));
    tools.appendChild(zb('Reset view',()=>{zoom.k=1;zoom.tx=0;zoom.ty=0;view(true);}));
    tools.appendChild(el('span',{text:zoom.k>1.02?`${zoom.k.toFixed(1)}× — drag to move`
      :'Scroll to zoom, drag to move'}));
    left.appendChild(tools);
    left.appendChild(el('p',{class:'plotcap',text:caption}));

    host.appendChild(el('div',{class:'split'},[left,el('div',{},[ctl,tableOf(f,sort,k=>{
      sort = sort.k===k?{k,d:-sort.d}:{k,d:typeof f[0]?.[k]==='number'?-1:1}; view(true);},sel,pick)])]));
    const dh=el('div',{id:'detailhost'}); host.appendChild(dh);
    if(sel) dh.appendChild(detail(sel,mode,m=>{mode=m;view(true);}));
    if(keepScroll) window.scrollTo(0,y);
  }
  view();
}

/* Nominations: opens on the bundled demo run, replaceable by any sweep CSV */
function nomsView(){
  const host=$('#noms'); host.textContent='';
  if(DEMO && DEMO.text) host.appendChild(el('p',{class:'note',style:'margin:0 0 14px',html:
    `Opening on a <b>bundled demo run</b> &mdash; ${DEMO.scope}, swept against the live databases on
     ${DEMO.when}. Real output, not a mock-up, but a hand-picked slice rather than a proteome sweep.
     Worth a look: <b>PCSK9 / coronary artery disease</b> scores validity 100 and composite 0. The
     axes multiply, so the best-validated target here is still the wrong thing to degrade.`}));
  const drop=el('div',{class:'drop',html:
    'Drop <b>sweep_results.csv</b> here to browse a run, or <b>click to choose a file</b>.<br>'
    +'<span style="font-size:13px">Nothing is uploaded &mdash; the file is read in this page.</span>'});
  const inp=el('input',{type:'file',accept:'.csv',style:'display:none'});
  drop.addEventListener('click',()=>inp.click());
  ['dragenter','dragover'].forEach(e=>drop.addEventListener(e,ev=>{ev.preventDefault();drop.classList.add('over');}));
  ['dragleave','drop'].forEach(e=>drop.addEventListener(e,ev=>{ev.preventDefault();drop.classList.remove('over');}));
  drop.addEventListener('drop',ev=>{if(ev.dataTransfer.files[0]) read(ev.dataTransfer.files[0]);});
  inp.addEventListener('change',e=>{if(e.target.files[0]) read(e.target.files[0]);});
  host.appendChild(drop); host.appendChild(inp);
  const body=el('div'); host.appendChild(body);
  if(DEMO && DEMO.text) load(DEMO.name, DEMO.text, true);
  else body.appendChild(el('p',{class:'empty',text:'No run loaded yet.'}));

  function read(file){
    const fr=new FileReader();
    fr.onload=()=>load(file.name, fr.result, false);
    fr.readAsText(file);
  }
  /* One path for the bundled run and for a dropped file, so the demo cannot
     render differently from the real thing. */
  function load(name, text, isDemo){
    let rows;
    try{ rows=parse(text); }
    catch(err){ body.textContent=''; body.appendChild(el('p',{class:'empty',
      text:'That file could not be read as a sweep result. Expected the columns written by RUN_SWEEP.py.'})); return; }
    if(!rows.length){ body.textContent=''; body.appendChild(el('p',{class:'empty',text:'That file has no rows.'})); return; }
    const dead=['v','m','degron'].filter(k=>rows.every(r=>r[k]===rows[0][k]));
    body.textContent='';
    if(dead.length) body.appendChild(el('div',{class:'disagree',
      html:`Every row has the same value for <b>${dead.join(', ')}</b>. A column that never
        varies carries no information &mdash; that layer probably failed rather than agreeing.`}));
    drop.innerHTML = isDemo
      ? '<b>'+name+'</b> — bundled demo run, '+rows.length+' rows. Drop your own '
        +'<b>sweep_results.csv</b> here to replace it.'
      : '<b>'+name+'</b> — '+rows.length+' rows. Drop another file to replace it.';
    browser(body, rows, 'Each point is one target in one indication. Points in the shaded corners '
      +'are the traps: strong on one axis, weak on the other.', 'cardiometabolic');
  }
  function parse(txt){
    const lines=txt.replace(/\r/g,'').split('\n').filter(Boolean);
    const head=split(lines[0]);
    const need=['gene','validity','modality','verdict'];
    if(!need.every(n=>head.includes(n))) throw new Error('missing columns');
    return lines.slice(1).map((l,i)=>{
      const c=split(l), o={}; head.forEach((h,j)=>o[h]=c[j]);
      const fl=(o.flags||'').split('|').map(s=>s.trim()).filter(Boolean);
      return {id:o.accession||String(i), gene:o.gene, indication:o.indication,
        area:o.area, track:o.track, compartment:o.compartment,
        v:+o.validity||0, m:+o.modality||0, c:+o.composite||0,
        verdict:o.verdict, vflags:fl, mflags:[], flags:[], notes:'',
        degron:+o.degron||0, expected:'', why:'', agree:true};
    });
  }
  function split(l){ /* handles the quoted flags column */
    const out=[]; let cur='', q=false;
    for(const ch of l){ if(ch==='"'){q=!q;} else if(ch===','&&!q){out.push(cur);cur='';} else cur+=ch; }
    out.push(cur); return out;
  }
}

function calView(){
  const host=$('#cal'); host.textContent='';
  const bad=CAL.filter(r=>!r.agree).length;
  host.appendChild(el('p',{class:'note',style:'margin:0 0 16px',html:
    `37 cases with known answers, labelled on both axes. The bench agrees with
     <b>${CAL.length-bad} of ${CAL.length}</b>. Treat that with suspicion: labels were revised
     while the scoring was being tuned, so the figure is partly fitted. The tab on the right
     is the test labels cannot fit.`}));
  const body=el('div'); host.appendChild(body);
  browser(body, CAL, 'The calibration set as the bench scores it. Shaded corners are the '
    +'traps — cases designed to catch a scorer that trusts one axis too much.');
}

function corpView(){
  const host=$('#corp'); host.textContent='';
  host.appendChild(el('p',{class:'note',style:'margin:0 0 16px',html:
    `Evidence that a protein's non-catalytic role matters, tiered by how much the comparison
     can be trusted. A protein can show a clean difference between knockout and inhibitor and
     still be the wrong thing to degrade: <b>RIPK1</b> does, and its scaffold is protective, so
     the corpus rejects it.`}));
  const t=el('table',{},[el('thead',{},[el('tr',{},
    ['Tier','Gene','Direction','Basis','Demoted by'].map(h=>el('th',{text:h})))]),
    el('tbody',{},CORP.map(r=>el('tr',{},[
      el('td',{text:r.tier,class:'n'}),
      el('td',{text:r.gene,class:'gene'}),
      el('td',{},[el('span',{class:'v '+(r.sign==='protective'?'reject':r.sign==='unknown'?'flag':'take'),
        text:r.sign==='protective'?'Protective — do not degrade'
          :r.sign==='unknown'?'Direction unclear':'Pathogenic'})]),
      el('td',{text:r.basis,style:'color:var(--dim);font-size:13px;max-width:46ch'}),
      el('td',{text:(r.confounders||'').split('|').filter(Boolean).map(s=>s.replace(/_/g,' ')).join(', ')||'—',
        style:'color:var(--faint);font-size:12.5px'}),
    ])))]);
  host.appendChild(el('div',{class:'scroll'},[t]));
}

function diagView(){
  const host=$('#diag'); host.textContent='';
  host.appendChild(el('p',{class:'note',style:'margin:0 0 18px',html:
    `Each scoring term is switched off in turn to see how many verdicts move. Labels cannot
     influence this, which is why it is the number to trust. A term that moves nothing is
     either doing no work or has no case in the set that exercises it.`}));
  const max=Math.max(...ABL.map(a=>a.n))||1;
  ABL.forEach(a=>{
    const row=el('div',{class:'ab'+(a.n?'':' zero')},[
      el('div',{class:'nm',text:a.term}),
      el('div',{class:'track'},[el('i',{style:`width:${(a.n/max)*100}%`})]),
      el('div',{class:'cnt',text:a.n?`${a.n}/${CAL.length}`:'none'}),
    ]);
    host.appendChild(row);
    if(a.n&&a.n<=6) host.appendChild(el('div',{class:'ab'},[el('div',{class:'who',text:a.genes.join(', ')})]));
  });
  const top=ABL[0];
  host.appendChild(el('p',{class:'note',style:'margin-top:22px',html:
    `<b>${top.term}</b> moves more verdicts than anything else — more than genetics, clinical
     outcome or scaffold evidence. It is also the hardest field to populate automatically,
     which is worth knowing before trusting a ranked list built on it.`}));
}

/* ---------- tabs ---------- */
const views={noms:nomsView,cal:calView,corp:corpView,diag:diagView};
let built={};
document.querySelectorAll('nav button').forEach(b=>b.addEventListener('click',()=>{
  document.querySelectorAll('nav button').forEach(x=>x.setAttribute('aria-selected',String(x===b)));
  document.querySelectorAll('main section').forEach(s=>s.hidden = s.id!==b.dataset.tab);
  if(!built[b.dataset.tab]){ views[b.dataset.tab](); built[b.dataset.tab]=1; }
}));
nomsView(); built.noms=1;
</script></body></html>
"""


def _js(obj) -> str:
    """JSON for embedding inside a <script> block. '</' is escaped so a string
    containing '</script>' cannot close the tag early."""
    return json.dumps(obj).replace("</", "<\\/")


def build() -> pathlib.Path:
    html = (HTML
            .replace("__CAL__", _js(scored_rows()))
            .replace("__ABL__", _js(ablation_rows()))
            .replace("__CORP__", _js(corpus_rows()))
            .replace("__PROF__", _js(profiles()))
            .replace("__GLOSS__", _js(glossary_patterns()))
            .replace("__DEMO__", _js(demo_sweep())))
    OUT.write_text(html, encoding="utf-8")
    return OUT


if __name__ == "__main__":
    p = build()
    print(f"wrote {p.resolve()}  ({p.stat().st_size/1024:.0f} KB)")
    print("Double-click it. No server, no network, nothing to install.")
