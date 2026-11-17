# Pre-registration: what US-style uACR testing does to the KFRE (NHANES 2007-2018), v1

Dated 2026-11-17. Written and published **before any NHANES data were downloaded or read**
by this project. Code: `kfretransport/nhanes.py` (version 0.2.2), tests `test_nhanes.py`.
Results will be published in this folder whatever they show. Any change will be a new,
dated addendum; this file will not be edited.

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it yet. Research only: information only, not medical advice. Text CC BY 4.0.

## Why

In routine care, only some people with CKD have a urine albumin-to-creatinine ratio (uACR),
and the 4-variable Kidney Failure Risk Equation (KFRE) needs one. In 291,155 insured US adults
with CKD (2015-2019), 63% of those with diabetes had albuminuria testing, against 25% of
those without (Chu et al., BMC Prim Care 2022;23:299, PMID 36434513, PMC9700954). In
Stockholm, the KFRE under-predicted among the minority with a uACR (Cleary et al., BMJ Med
2026, PMID 41658049), and an editorial described over-estimation among those without one
(Walker & Lees, PMID 41658048). Health-system data cannot show what the untested people's
uACR really was. NHANES can: every examined adult gave a urine sample, so the KFRE inputs
are known for everyone, tested in routine care or not.

Prior work used NHANES to describe the KFRE risk of US adults with CKD (e.g. Chu et al., AJKD
2020, PMID 32305206: awareness by KFRE risk). On 2026-11-17 a Europe PMC search ("kidney
failure risk equation" OR KFRE, AND NHANES: 39 records) found no study that applies
real-world testing patterns to NHANES to measure (a) how many high-KFRE-risk adults would go
unscored or (b) how far imputing uACR from tested people misstates the untested group's risk.

## Data (public)

- NHANES 2007-2008 to 2017-2018 (six cycles): DEMO, BIOPRO (serum creatinine LBXSCR; the
  BIOPRO_E notes say 2007 values were adjusted to the IDMS-standardized 2008 method), ALB_CR
  (URDACT, mg/g), DIQ (DIQ010), BPQ (BPQ020), from wwwn.cdc.gov.
- NCHS public-use linked mortality files, follow-up to 2019-12-31 (ELIGSTAT, MORTSTAT,
  UCOD_LEADING; 9 = nephritis, nephrotic syndrome and nephrosis).

## Cohort

Age 20+, not pregnant, MEC exam weight > 0, serum creatinine and uACR present, CKD-EPI 2021
eGFR in [15, 60) from a single measurement (CKD G3-G4 at one visit). Weights: WTMEC2YR / 6.
Diagnosed diabetes = DIQ010 = 1; diagnosed hypertension = BPQ020 = 1 (self-report, used as
what a clinician would know). 5-year KFRE, 4-variable, North American baseline.

## Testing policies (probability that a person has a uACR in routine care)

| Policy | Diagnosed diabetes | No diagnosed diabetes |
|---|---|---|
| P0 universal | 1 | 1 |
| **P1 US claims (primary)** | 0.63 | 0.25 (Chu 2022) |
| P2 guideline | 1 | 1 if diagnosed hypertension, else 0 |
| P3 diabetes only | 1 | 0 |
| P1q curve (exploratory) | 0.63 | q = 0.05, 0.10, 0.25, 0.50, 0.75, 1.0 |

Estimates use expected weights (survey weight x probability), so nothing is random.

## Estimands (survey-weighted; delete-one-PSU jackknife SE; 95% CI with t on PSUs - strata df)

- **selection_ratio (primary)**: the untested group's mean 5-year KFRE with uACR imputed
  from tested people in the same cell (eGFR 15-29/30-44/45-59 x age <65/65+ x sex; weighted
  mean of log uACR; fallback to the eGFR band), divided by their mean KFRE with their own
  measured uACR (= imputation_ratio), and then divided by the same ratio under
  non-selective testing (imputation_ratio_uniform). The uniform reference is needed
  because cell-mean imputation is biased even with random testing: the KFRE is not linear
  in log uACR. selection_ratio isolates the part caused by who gets tested.
- high_scored_share: share of adults with 5-year KFRE >= 5% (KDIGO 2024 practice point
  2.2.1 uses 3-5% for referral) who would have a uACR, and so a KFRE.
- high_unscored_weighted: the estimated number of US adults behind the unscored share.
- tested_all_ratio: mean KFRE among tested / mean KFRE among everyone.
- kdeath_unscored_share (exploratory): share of kidney-cause deaths through 2019 that
  happened among people who would have been untested.

## Predictions, written before any data were read

- **H1 (primary):** under P1, selection_ratio > 1, with a 95% CI lower bound above 1.
  The prior guess for its size is 1.05-1.35.
- H2: under P1, high_scored_share (5% threshold) < 0.60. Prior guess about 0.45.
- H3: under P1, tested_all_ratio > 1, with a 95% CI excluding 1.
- H4: under P2, high_scored_share (5% threshold) >= 0.85.
- Exploratory, no prediction: P3, the 3% threshold, the q curve, kidney-cause deaths.

## What this can and cannot show

It can show how much of the "over-estimation in the untested" pattern one measurable driver
of testing (diagnosed diabetes) produces by itself, in a nationally representative sample,
and how many high-risk US adults a diabetes-driven testing pattern leaves without a KFRE.
It cannot show KFRE calibration. The public files have no kidney-failure (KRT) outcome, and
kidney-cause death is not kidney failure. Other limits: a single eGFR and uACR, with no
3-month chronicity, so G3a is over-counted, especially in older adults. Self-reported
diagnoses stand in for what clinicians know. The testing rates come from one insured
2015-2019 population. Real testing also depends on things NHANES cannot see, such as
nephrology care and clinician concern. Few G4 participants.

## The "wow" if it works

A first nationally representative measurement showing that US testing patterns alone make
uACR-imputed kidney-failure risk wrong for the people who were never tested, by a stated
amount. It would also show how many high-risk US adults a diabetes-centred testing pattern
leaves unscored. This is the mechanism a multi-site KFRE study (protocol addendum 1) has to
correct for.
