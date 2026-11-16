# kfre-transport

Open, tested building blocks for one real-world-evidence question: **does the 4-variable Kidney
Failure Risk Equation (KFRE) stay calibrated when it moves between health systems that differ in
who gets a urine albumin (uACR) test?**

> Research code only. **Information only, not medical advice, not a medical device.** It is not a
> clinical calculator and must not be used for decisions about any person.
> Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. (mv2@munimtech.com).
> No human expert has reviewed it yet. Code Apache-2.0; text CC BY 4.0.

## Why this question

- In Stockholm's SCREAM data only 26,229 of 116,158 adults with CKD stages 3-4 had a uACR result,
  and among them the 4-variable KFRE under-predicted (observed/expected 2.11)
  (Cleary et al., BMJ Med 2026;5:e001950, PMID 41658049).
- The accompanying editorial (Walker & Lees, BMJ Med 2026;5:e002559, PMID 41658048) describes
  under-estimation in patients with uACR data and over-estimation in those without, and suggests
  unmeasured factors drive both the decision to test and the risk.
- If so, a model that is calibrated at one health system can fail at another purely because the
  two test different people. Complete-case analysis hides this. NIH NCATS's Real-World Evidence
  Innovation Challenge names kidney disease as a priority area and asks for methods whose
  reproducibility, transportability and stability can be measured across health systems.

## What is here (version 0.1.0)

| Module | What it does | Checked by |
|---|---|---|
| `kfretransport/equations.py` | CKD-EPI 2021 (race-free) and 2009 creatinine eGFR; 4-variable KFRE (2- and 5-year; North American and non-North American baselines); KDIGO G and A categories; unit conversions | hand-worked unit tests (`test_equations.py`) |
| `kfretransport/metrics.py` | O/E, Brier, C-statistic, logistic calibration slope/intercept and calibration-in-the-large, DerSimonian-Laird pooling with I^2 and a 95% prediction interval for a new site, random-effects meta-regression with Knapp-Hartung inference | hand-worked unit tests (`test_metrics.py`) |
| `kfretransport/sim.py` | synthetic multi-site study of informative uACR testing (4 scenarios) | property tests (`test_sim.py`); the full study is pre-registered and queued, **not run yet** |

Python 3.9+, standard library only, no network, no data files.

```
python -m unittest discover -v
python -m kfretransport.sim --k 200 --n 10000 --seed 20261116 --out-dir out
```

## Status

- Unit tests: hand-traced; they must pass before any commit to this repository.
- Simulation study v1: pre-registered in
  [findings/sim-prereg-v1-2026-11-16.md](findings/sim-prereg-v1-2026-11-16.md) with seven
  predictions (S1-S7) written before any run. Results will be published here whatever they show.
- Not done: any real data, any OMOP/N3C code, competing risk of death, censoring, expert review.

## Formulas and their sources

KFRE 4-variable: risk = 1 - S0^exp(L), where

    L = -0.2201 (age/10 - 7.036) + 0.2467 (male - 0.5642) - 0.5567 (eGFR/5 - 7.222) + 0.4510 (ln uACR[mg/g] - 5.137)

and S0 = 0.9750 (2 y) / 0.9240 (5 y) North
American, 0.9832 / 0.9365 elsewhere (Tangri et al., JAMA 2011 and JAMA 2016; as reproduced in the
CDC supplement "Item S1. Calculation of 5-Year Kidney Failure Risk"). Some secondary sources give
0.9751 for the North American 2-year baseline or drop the sex centering; this code keeps 0.9750 and
the centered form. eGFR: Inker et al., NEJM 2021 (CKD-EPI 2021); Levey et al., Ann Intern Med 2009.
Please check every formula against the primary papers before relying on it, and open an issue if
anything is wrong.

## Licence

Apache-2.0 (code). Text and findings CC BY 4.0.
