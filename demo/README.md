# Demo data

Real sweep output, bundled so the viewer has something to show without anyone
waiting on a sweep. Open `tpd_viewer.html` and the Nominations tab opens on
`sweep_results_demo.csv` already loaded.

Nothing here is synthetic. Both files are verbatim output from
`RUN_SWEEP.py`, written on 28 October 2026 against the live public APIs
(UniProt, Open Targets, gnomAD, AlphaFold DB). If you re-run the sweep today
the numbers may differ: the sources are living databases.

---

## `sweep_results_demo.csv` — the run the viewer opens with

24 rows: 6 genes scored against all 4 indications.

| | |
|---|---|
| Genes | STAT6, IRF5, IL23R, PNPLA3, HSD17B13, PCSK9 |
| Indications | atopic dermatitis, systemic lupus erythematosus, non-alcoholic steatohepatitis, coronary artery disease |
| Nominated | 3 (`mid`) |
| Rejected | 21 |

**This is a 6-gene slice, not a proteome sweep.** It was run to confirm the
feature pipeline after the compartment and genetics fixes, so the gene list was
chosen by hand. It demonstrates the scoring; it is not a result.

Four things in it are worth pointing at during a demo:

- **STAT6 / atopic dermatitis → `mid`.** Intracellular, replicated protective
  loss-of-function genetics, degron score 0.87. The axis dupilumab validated,
  reached at a node no antibody can.
- **PNPLA3 / NASH → `mid`.** Same genetics class, and the hepatic restriction
  flag relaxes the delivery gate — the one place where the oral requirement
  softens because the tissue of action is the liver.
- **PCSK9 / coronary artery disease → `reject_wrong_modality`, despite validity
  100.** The best-validated target in the file scores composite 0.0. Validity
  and modality are combined by geometric mean, so a perfect score on one axis
  cannot rescue a zero on the other. PCSK9 is secreted, inclisiran already
  doses twice a year, and there is no displacement argument. This is the single
  most useful cell in the file: it is the bench refusing to reward a famous
  target.
- **Cross-area pairings.** PNPLA3 against lupus, PCSK9 against atopic
  dermatitis. The sweep scores every protein against every indication and most
  pairings come back invalid, which is what it should do.

Every row carries `direction of effect unknown` and `score rests on 1 evidence
layer`. That is accurate, not a gap in the demo file: pathway, expression and
direction of effect are not yet wired, and terms with no data are renormalised
out rather than scored as zero. Read the genetics-driven ranks as
systematically conservative.

---

## `sweep_results_failed_run.csv` — a broken run, kept on purpose

160 rows, validity 0.0 in every one of them.

This is the sweep from earlier the same day, before two bugs were fixed: a
swallowed exception in the genetics lookup, and a compartment classifier that
matched the string "Membrane" anywhere in a UniProt annotation and so routed
cytosolic proteins to the extracellular track.

It is here because it looks completely normal. 160 rows, every column
populated, no error anywhere. Drag it onto the Nominations tab and the viewer
says so: *"Every row has the same value for v. A column that never varies
carries no information — that layer probably failed rather than agreeing."*

Worth 30 seconds in a demo. A pipeline that cannot tell you it failed is worse
than one that crashes.

---

## Using your own run

`RUN_SWEEP.py` writes `sweep_results.csv` to the repo root. Drag it onto the
Nominations tab — it replaces the demo and nothing leaves the page. The viewer
reads the file locally; there is no server and no network call in it.
