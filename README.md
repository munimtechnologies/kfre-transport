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
- In 291,155 insured US adults with CKD (2015-2019), 63% of those with diabetes had albuminuria
  testing, against 25% of those without (Chu et al., BMC Prim Care 2022;23:299, PMID 36434513).
- If testing is this selective, a model that is calibrated at one health system can fail at
  another purely because the two test different people. Complete-case analysis hides this.
  NIH NCATS's Real-World Evidence Innovation Challenge names kidney disease as a priority area
  and asks for methods whose reproducibility, transportability and stability can be measured
  across health systems.
- Death competes with kidney failure. In 30,801 Alberta adults with stage 4 CKD more people died
  than reached kidney failure, and treating death as censoring overstated 5-year kidney-failure
  risk by about 7% (Ravani et al., JASN 2019;30:2219-2227, PMID 31540963). A transportability
  study must therefore say which "observed risk" it compares the KFRE with.

## New in 0.2.3: a closed-form selection-calibration curve

If the part of kidney-failure risk that the KFRE's predictors miss is correlated with who gets
tested, then among the tested

    O/E (tested) = Phi(Phi^-1(p) + theta) / p      relative to the whole cohort's O/E,

where p is the site's testing rate and theta is a single number: the covariance between the
unexplained log-hazard and the testing propensity (Gaussian unexplained risk, probit testing;
exact for the hazard, and for the risk when the outcome is rare). Its consequences:

- Selection on the KFRE's own predictors does not miscalibrate the KFRE.
- Less testing means more apparent under-prediction among the tested, and the relationship is
  curved, not linear.
- The untested get the mirror image.

If all of Stockholm's O/E of 2.11 came from selective testing, theta would be 0.69. Derivation,
tables, limits and the prior-art search are in
[findings/selection-curve-2026-11-17.md](findings/selection-curve-2026-11-17.md). It is a
pre-specified secondary analysis
([findings/protocol-addendum-2-2026-11-17.md](findings/protocol-addendum-2-2026-11-17.md)).
It has **not** been run on any real data.

## What is here (version 0.2.3)

| Module | What it does | Checked by |
|---|---|---|
| `kfretransport/equations.py` | CKD-EPI 2021 (race-free) and 2009 creatinine eGFR; 4-variable KFRE (2- and 5-year; North American and non-North American baselines); KDIGO G and A categories; unit conversions | hand-worked unit tests (`test_equations.py`) |
| `kfretransport/metrics.py` | O/E, Brier, C-statistic, logistic calibration slope/intercept and calibration-in-the-large, DerSimonian-Laird pooling with I^2 and a 95% prediction interval for a new site, random-effects meta-regression with Knapp-Hartung inference, Aalen-Johansen cumulative incidence (death competing), 1 - Kaplan-Meier (death censored), time-to-event O/E; delta-method (Greenwood-type) variance of the Aalen-Johansen estimate and `se_log_oe_dm` | hand-worked unit tests (`test_metrics.py`), including checks that the variance equals the multinomial variance without censoring and Greenwood's variance without competing events |
| `kfretransport/pool.py` | pools site summaries from `omop.py` into the pre-declared primary and secondary analyses (random-effects O/E and meta-regression on the site's testing rate) and masks small cells in published per-site rows | hand-worked unit tests (`test_pool.py`) |
| `kfretransport/selection.py` (new) | the selection-calibration curve for tested and untested people, theta from one site's O/E, its exact risk-scale version by quadrature, and a multi-site fit of theta with a profile-likelihood interval | unit tests (`test_selection.py`): hand values from normal tables, exact identities, agreement of the quadrature with the closed form for rare outcomes, recovery of a known theta |
| `kfretransport/nhanes.py` | pre-registered public-data study: applies US testing patterns to NHANES 2007-2018, where everyone has a uACR, and measures what they do to KFRE scoring and uACR imputation | hand-worked unit tests on invented records (`test_nhanes.py`); **not run on NHANES yet** |
| `kfretransport/albuminuria.py` | estimated uACR from a urine protein-to-creatinine ratio (Sumida et al. 2020 crude equation), for a sensitivity analysis only | tests that it reproduces five figures the paper reports (`test_albuminuria.py`); see provenance note below |
| `kfretransport/sim.py` | synthetic multi-site study of informative uACR testing (4 scenarios) | property tests (`test_sim.py`); the full study is pre-registered and queued, **not run yet** |
| `kfretransport/conceptsets.py` | concept sets written as LOINC / SNOMED / ICD10CM / UCUM **codes**, not OMOP ids; unit tables (also uPCR); plausible ranges | each code checked on 2026-11-16 against tx.fhir.org (LOINC 2.83, SNOMED CT International 2026-07-01, ICD-10-CM 2026-04-01) |
| `kfretransport/omop.py` | reads an OMOP CDM v5.4 extract (CSV or Athena-style TSV), resolves the concept sets through the site's own vocabulary, builds the KFRE cohort and outcome, writes a site summary with no person-level rows (complete-case, Aalen-Johansen and Kaplan-Meier O/E side by side), and reports vocabulary drift | hand-traced tests on an invented 7-person extract (`test_omop.py`) |

Python 3.9+, standard library only, no network (except `nhanes.py`, which downloads public
NHANES files and needs pandas to read them).

```
python -m unittest discover -v
python -m kfretransport.sim --k 200 --n 10000 --seed 20261116 --out-dir out
python -m kfretransport.omop --cdm-dir path/to/omop_csv --out site-a.json
python -m kfretransport.omop --cdm-dir path/to/omop_csv --out site-a-pcr.json --pcr-conversion
python -m kfretransport.pool site-a.json site-b.json site-c.json --out pooled.json --min-cell 20
python -m kfretransport.nhanes --data-dir nhanes-data --out-dir out
python -m kfretransport.selection --theta 0.1 0.3 0.69 --rates 0.1 0.25 0.4 0.63 0.9
```

## Pre-registrations and protocol

- [findings/protocol-addendum-1-2026-11-17.md](findings/protocol-addendum-1-2026-11-17.md):
  the primary multi-site analysis, fixed before any real data. It uses log O/E with
  Aalen-Johansen O and the delta-method variance, random effects, and meta-regression on the
  site's measured-uACR testing rate. The predicted slope is negative.
- [findings/protocol-addendum-2-2026-11-17.md](findings/protocol-addendum-2-2026-11-17.md):
  a secondary analysis that fits the selection-calibration curve across sites; prediction
  theta > 0.
- [findings/sim-prereg-v1-2026-11-16.md](findings/sim-prereg-v1-2026-11-16.md): simulation
  study, seven predictions (S1-S7), queued.
- [findings/nhanes-testing-prereg-v1-2026-11-17.md](findings/nhanes-testing-prereg-v1-2026-11-17.md):
  NHANES study, four predictions (H1-H4), queued. Its primary estimand is the
  `selection_ratio`: how much more imputing uACR from tested people misstates the untested
  group's KFRE under US testing rates than under non-selective testing.

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

### Three ways to count "observed"

`site_summary()` reports, for each horizon (2 and 5 years):

| Key | Who is counted | Observed risk |
|---|---|---|
| `oe` | people with a prediction **and** a known horizon outcome (complete case) | share with KRT by the horizon (death first counts as no KRT); people censored earlier are left out |
| `tte_aj` | **everyone** with a prediction | Aalen-Johansen cumulative incidence of KRT, death competing |
| `tte_km` | **everyone** with a prediction | 1 - Kaplan-Meier, death treated as censoring |

The expected risk is the mean KFRE prediction of the same people. Each time-to-event entry has
two standard errors of log O/E: `se_log_oe` (Poisson, 1/sqrt(events)) and `se_log_oe_dm`
(delta-method variance of the Aalen-Johansen or Kaplan-Meier estimate, E held fixed; the form
printed in Marubini and Valsecchi 1995 and in the R package cifmodeling's formula vignette). The
primary analysis uses `tte_aj` with `se_log_oe_dm`.

### Small cells

`pool.py --min-cell 20` (the default) masks counts from 1 to 19 in published per-site rows. It
also masks anything that would give them back: the testing rate, a count recoverable by
subtraction, and the Poisson variance (which is 1/events). This follows N3C's rule that no
result about fewer than 20 people leaves its enclave (N3C's download guidance at n3c.ncats.nih.gov
and the N3C governance paper, PMC10789985). Each enclave's own rules come first.

### uPCR to uACR conversion (sensitivity analysis only, off by default)

With `--pcr-conversion`, people with no uACR in the window but a urine protein-to-creatinine
ratio (LOINC 2890-2) get an estimated uACR from the crude equation of Sumida et al., Ann Intern
Med 2020;173:426-435 (PMID 32658569; 919,383 adults, CKD Prognosis Consortium):

    ACR = exp(5.3920 + 0.3072 ln(min(PCR/50, 1)) + 1.5793 ln(max(min(PCR/500, 1), 0.1)) + 1.1266 ln(max(PCR/500, 1)))

with PCR and ACR in mg/g. "Tested" still means a **measured** uACR, so the testing rate does not
change. Provenance, stated plainly: the agents could not open the paper's PDF. They took the
coefficients from a search-engine extract of the paper's Table 2. They then checked that the
equation reproduces five figures reported for the crude model (PCR 150 -> ACR 33, PCR 500 -> 220,
ACR 30 <-> PCR 142, ACR 300 <-> PCR 660, and 2.99-fold per doubling of PCR between 50 and 500).
The slope below PCR 50 mg/g is not checked by any reported figure, and the paper says the
association is inconsistent there. The adjusted equation and the dipstick conversion are not
implemented. Please check against the paper before any use, and open an issue if anything
differs.

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
- Simulation study v1 and NHANES study v1: pre-registered, queued, **not run yet**. Results
  will be published here whatever they show.
- Selection-calibration curve (0.2.3): an analytic result plus tested code; not yet fitted to
  any real or simulated multi-site data.
- OMOP code: tested only on invented records. It has **not** yet run on any OMOP database,
  synthetic or real. A smoke test on OHDSI's public synthetic Synthea27Nj set is queued.
- Not done: any real health-system data; a container; expert review.

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
