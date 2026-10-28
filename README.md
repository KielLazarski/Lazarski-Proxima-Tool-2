# TPD Glue Bench

A database-mining tool that nominates degradation targets in immunology,
inflammation and cardiometabolic disease.

<!-- Replace GITHUB-USERNAME below with your handle, then delete this comment. -->
**[Open the tool](https://GITHUB-USERNAME.github.io/Lazarski-Proxima-Interview-Tool-2/)**
— runs in the browser, nothing to install, opens on a real sweep.

Outside oncology there is no essentiality screen to borrow, and the
disease-biology target lists are picked over. This scores candidates on two
independent axes and returns a reasoned verdict rather than a single number.

## Quick start

Python 3.9+. No dependencies — standard library only.

```bash
python3 OPEN_VIEWER.py       # browse the analysis in a browser
python3 START_HERE.py        # runs every check, writes report.txt
python3 RUN_SWEEP.py         # scores targets against live databases
```

`OPEN_VIEWER.py` builds a single self-contained HTML page and opens it. No
server, no network, nothing to install. It opens on a **bundled demo sweep**
(`demo/`), so there is something real to look at without waiting on a run, and
also shows the calibration set as the bench scores it, the scaffolding corpus,
and which scoring terms actually move verdicts. Drop your own
`sweep_results.csv` on it to replace the demo.

The scatter is the point: target validity on one axis, modality fit on the
other, with the two trap corners shaded. Click a point to see every reason
behind its score, each one cited, with a plain/detailed toggle. The weight
sliders re-score live; hard rejections are deliberately immune to them.

## Demoing it

Open the viewer and go to **PCSK9 / coronary artery disease** in the bundled
run. Validity 100 — the best-validated target in the file — and composite 0.
The axes are combined by geometric mean, so a perfect score on one cannot
rescue a zero on the other, and PCSK9 is secreted with inclisiran already
dosing twice a year. That one cell is the whole argument for scoring two axes
separately.

`demo/sweep_results_failed_run.csv` is a second, deliberately broken run. Drop
it on the page and the viewer says the validity column never varies, so that
layer failed rather than agreed. See `demo/README.md`.

`START_HERE.py` needs no network. `RUN_SWEEP.py` needs access to UniProt,
Open Targets, gnomAD, AlphaFold and ChEMBL.

`index.html` is the built page, committed rather than ignored because it is
what GitHub Pages serves. Rebuild and commit it after changing anything the
page displays — `git status` will show it. `.nojekyll` next to it stops Pages
putting the file through Jekyll.

## The two axes

| Axis | Question |
| --- | --- |
| **Target validity** | Does this protein causally drive the disease? |
| **Modality fit** | Is degradation the right tool, versus inhibition or an antibody? |

They are never collapsed into one number. A target can be real and
undegradable (PCSK9) or degradable and not real (PPARα in MASH). Those are
different failure modes needing different responses. The composite is a
geometric mean, so a strong score on one axis cannot rescue a weak score on
the other.

## Verdicts

| Verdict | Meaning |
| --- | --- |
| `top` / `mid` | Take it forward |
| `nominate_for_validation` | Modality case is strong, validity unproven. Run a perturbation experiment. |
| `reject_wrong_modality` | Real target, but degradation cannot reach it |
| `reject_modality_too_weak` | Reachable, but the modality case does not survive the competitive landscape |
| `reject_inhibition_sufficient` | A tolerated inhibitor already works, or kinase-dead phenocopies knockout |
| `reject_scaffold_protective` | Scaffold evidence exists and points the wrong way |
| `reject_target_invalid` | Target refuted or validity below bar |
| `flag_contested` | Trial failed but the genetic case stands; resolve before scoring |
| `flag_strategic_discontinuation` | Programme dropped for portfolio reasons, not science |

## Layout

```
tpd_bench/
  schema.py              feature definitions — objective, mineable facts only
  scoring.py             both axes, disease-area weight profiles, hard gates
  targets.py             curated records for the calibration cases
  calibration_labels.csv ground truth — never read by the scorer
  calibrate.py           harness; prints every disagreement
  ablate.py              which scoring terms actually move verdicts
  corpus/                scaffolding discordance corpus
  features/              live data adapters, degron scorer, proteome sweep
  glossary.py            cited explanation for every flag the page can show
  viewer.py              builds the self-contained HTML page
demo/
  sweep_results_demo.csv       real sweep output, bundled so the viewer opens full
  sweep_results_failed_run.csv a broken run, kept to show the dead-column check
  README.md                    provenance and what to point at
docs/
  TPD_Glue_Bench_Status.docx   two-page status summary
  calibration_set_annotated.csv
  frozen_predictions.md        predictions recorded before their readouts
  design_notes.md
```

## The one architectural rule

Features and labels live in separate files, and `scoring.py` never imports the
labels. If a feature was chosen to make the answer come out right, the
calibration report is theatre.

## How to judge it

The harness reports 36/36 on the calibration set. **Discount that.** Labels
were revised while tuning, so it is partly fitted.

Two things are not:

- `python3 -m tpd_bench.ablate` measures which scoring terms actually change
  verdicts. Labels cannot fit it. It currently shows `inhibitor status` as the
  most load-bearing term in the whole bench — worth deciding you believe before
  building further on it.
- `docs/frozen_predictions.md` holds predictions recorded before their
  readouts, with the result that would falsify the platform rather than just
  one molecule.

## Known limits

- **Direction of effect is unresolved.** Open Targets reports that a gene is
  associated with a disease, not whether the associated variants are protective
  loss-of-function — and degradation only follows from that. Genetics-driven
  ranks are provisional and flagged in the output.
- Only one validity layer is populated in the sweep, so scores rest on genetics
  alone.
- All three evidence layers are observational; none distinguishes *associated*
  from *required*.
- Tier-2 scaffolding evidence cannot be mined at scale: knockout phenotypes are
  in databases, inhibitor phenotypes are in prose.
- Expression and network data (GTEx, STRING) are specified, not wired.
- Scored on pre-July-2026 evidence, IL6R would have passed. It had replicated
  protective loss-of-function genetics and MR support, and ZEUS still returned
  HR 0.99. No configuration of the current validity layer rejects it in advance.
