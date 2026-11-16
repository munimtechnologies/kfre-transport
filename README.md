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

## What is here (version 0.2.0)

| Module | What it does | Checked by |
|---|---|---|
| `kfretransport/equations.py` | CKD-EPI 2021 (race-free) and 2009 creatinine eGFR; 4-variable KFRE (2- and 5-year; North American and non-North American baselines); KDIGO G and A categories; unit conversions | hand-worked unit tests (`test_equations.py`) |
| `kfretransport/metrics.py` | O/E, Brier, C-statistic, logistic calibration slope/intercept and calibration-in-the-large, DerSimonian-Laird pooling with I^2 and a 95% prediction interval for a new site, random-effects meta-regression with Knapp-Hartung inference | hand-worked unit tests (`test_metrics.py`) |
| `kfretransport/sim.py` | synthetic multi-site study of informative uACR testing (4 scenarios) | property tests (`test_sim.py`); the full study is pre-registered and queued, **not run yet** |
| `kfretransport/conceptsets.py` (new) | concept sets written as LOINC / SNOMED / ICD10CM / UCUM **codes**, not OMOP ids; unit table; plausible ranges | each code checked on 2026-11-16 against tx.fhir.org (LOINC 2.83, SNOMED CT International 2026-07-01, ICD-10-CM 2026-04-01) |
| `kfretransport/omop.py` (new) | reads an OMOP CDM v5.4 extract (CSV or Athena-style TSV), resolves the concept sets through the site's own vocabulary, builds the KFRE cohort and outcome, writes a site summary with no person-level rows, and reports vocabulary drift | hand-traced tests on an invented 7-person extract (`test_omop.py`) |

Python 3.9+, standard library only, no network.

```
python -m unittest discover -v
python -m kfretransport.sim --k 200 --n 10000 --seed 20261116 --out-dir out
python -m kfretransport.omop --cdm-dir path/to/omop_csv --out site-summary.json
```

### OMOP cohort definitions (pre-declared in `omop.py`)

- eGFR is always recomputed from serum creatinine (LOINC 2160-0, 38483-4, 14682-9, 59826-8;
  mg/dL, umol/L or mmol/L) with CKD-EPI 2021. The site's stored eGFR is never used.
- Index date: the first eGFR in [10, 60) that is at least 90 days into an unbroken run of eGFR
  below 60. Adults only; sex must be recorded as male or female, because the KFRE needs it.
- Excluded: any kidney replacement therapy (KRT) record on or before the index date.
- uACR (LOINC 9318-7, 14959-1, 32294-1): the latest value in the 365 days up to the index
  date, converted to mg/g. People without one stay in the cohort as "untested".
- Outcome: the first KRT record after the index date (dialysis, transplant, dialysis
  dependence, transplant status), matched on standard **or** source concepts (ICD10CM Z99.2,
  Z94.0). Death first is a competing event. Follow-up ends at the end of the observation period.
- The output has counts and rates only, never person-level rows.

## Vocabulary drift: SNOMED retired "End-stage renal disease" in June 2026

SNOMED CT International inactivated 46177005 "End-stage renal disease" and its subconcepts on
2026-06-01, folding them into "Chronic kidney disease stage 5". The OHDSI vocabulary release of
2026-08-29 still uses SNOMED 2026-02-01, so OMOP sites will meet this change at their next
vocabulary refresh, each at a different time. That is why the KFRE outcome here is built from
KRT records and never from the ESRD diagnosis. `omop.vocabulary_drift()` also reports where each
site stands. Details and sources:
[findings/snomed-esrd-retirement-2026-11-16.md](findings/snomed-esrd-retirement-2026-11-16.md).

## Status

- Unit tests: hand-traced. They must pass before any commit to this repository.
- Simulation study v1: pre-registered in
  [findings/sim-prereg-v1-2026-11-16.md](findings/sim-prereg-v1-2026-11-16.md) with seven
  predictions (S1-S7) written before any run. The results will be published here whatever
  they show.
- OMOP code: tested only on invented records. It has **not** yet run on any OMOP database,
  synthetic or real.
- Not done: any real data, censoring-aware estimators (for example, Aalen-Johansen), small-cell
  suppression rules for a given data enclave, a container, and expert review.

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
