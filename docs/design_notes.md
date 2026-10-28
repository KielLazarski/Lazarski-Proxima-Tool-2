# Design notes

Longer-form notes behind the decisions in the top-level README.


Target-nomination bench for molecular-glue degraders in immunology, inflammation
and cardiometabolic disease. Successor to the TCIP pair bench, re-aimed from
anti-cancer pair finding to degradation-target nomination.

## Run it

```bash
python3 -m tpd_bench.calibrate
```

## What it does

Scores each target-indication pair on **two independent axes** and emits a
verdict, never a single collapsed number.

- **Validity** — is this protein causally driving the disease?
- **Modality fit** — is degradation the right tool?

A target can be real and undegradable (PCSK9) or degradable and not real
(PPARα). Those are different failure modes needing different responses, so
they get different verdicts. The composite is a geometric mean, so a high score
on one axis cannot rescue a low score on the other.

## Verdicts

| Verdict | Meaning |
| --- | --- |
| `top` / `mid` | Take it forward |
| `nominate_for_validation` | Modality case is strong, validity unproven. Run a perturbation experiment. |
| `reject_wrong_modality` | Real target, but degradation cannot reach it |
| `reject_inhibition_sufficient` | A tolerated inhibitor already works, or kinase-dead phenocopies knockout |
| `reject_scaffold_protective` | Scaffold evidence exists and points the wrong way (RIPK1) |
| `reject_target_invalid` | Target refuted or validity below bar |
| `flag_strategic_discontinuation` | Program dropped for portfolio reasons, not science |

## Files

| File | Role |
| --- | --- |
| `schema.py` | Feature definitions. Objective, mineable facts only. |
| `scoring.py` | Both axes, disease-area weight profiles, hard gates |
| `targets.py` | Feature records for the calibration cases |
| `calibration_labels.csv` | Ground truth. **Never read by the scorer.** |
| `calibrate.py` | Harness. Prints every disagreement. |
| `predictions.md` | Frozen prospective predictions, recorded before readouts |

## The one architectural rule

Features and labels are separate files, and `scoring.py` never imports the
labels. If a feature was chosen to make the answer come out right, the
calibration report is theater.

## Two tracks

Compartment routes rather than rejects. Nuclear and cytosolic targets go to the
ubiquitin-proteasome track (glues, PROTACs, oral required). Membrane and
secreted targets go to the lysosomal track (conjugates, parenteral accepted).
They are ranked as separate lists because they are judged against different
delivery bars.

Extracellular nominations must clear a displacement test: a catalytic-turnover
or non-neutralisable-species argument, against a hard penalty for approved or
Phase 3 competition. PCSK9 fails it -- the catalytic argument is real, but
inclisiran is already twice-yearly, so turnover buys no dosing interval.

## Run the ablation

```bash
python3 -m tpd_bench.ablate
```

Asks which scoring terms actually move verdicts. This is the diagnostic to
trust, because labels cannot fit it.

`inhibitor status` moves more verdicts than any other term -- more than
genetics, clinical outcome or scaffold evidence. The bench rests on that one
judgement more heavily than on anything else, which is worth deciding you
believe before building further on it.

Still inert: compensatory upregulation and network redundancy. Neither is
necessarily wrong; the set contains no case that exercises them, which is a
gap in the set rather than a defence of the term.

## Honest status of the 100%

The harness reports 100% exact agreement. **Do not believe it.** Seventeen of
35 labels were written or changed while watching the engine behave: six during
scoring tuning, eleven when the extracellular track was added. A score on a set
edited while tuning is fitted, not calibrated. The harness says so itself when
it hits zero disagreements.

Trust the ablation and the frozen predictions instead.

`predictions.md` holds the real test.

## Known blind spot

Scored on pre-July-2026 evidence, IL6R would have passed. It had replicated
protective LOF genetics and MR support, and ZEUS still returned HR 0.99. No
configuration of the current validity layer rejects it in advance. The
`mr_unmatched_discount` term is a partial mitigation, not a fix.

## Not yet built

- Scaffolding discordance corpus (the long pole — the evidence base for the
  strongest entry criterion, and it does not exist as a mineable resource)
- E3 expression layer
- Protein turnover layer with inference for missing values
- Degron feature scoring from structure — currently hand-entered
- Property/developability gate (unsatisfied HBD ≤ 2 is the hard constraint)

## Open design decision

The oral-systemic requirement effectively pins the platform to CRBN, which is
broadly expressed, so cell-type selectivity from restricted E3 expression is not
available on the primary path today. The bench currently does not score E3
selectivity at all. That is a deliberate omission pending a ruling, not an
oversight.
