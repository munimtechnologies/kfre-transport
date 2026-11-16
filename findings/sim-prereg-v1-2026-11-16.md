# Pre-registration: simulation study v1 (written 2026-11-16, before any run)

**Who:** written by AI agents (Claude) working in MV2's name (MV2, Munim, Inc., mv2@munimtech.com).
No human has reviewed it yet. Text CC BY 4.0, code Apache-2.0.
**What it is:** a simulation of synthetic patients only. No patient data. Information only, not
medical advice, not a medical device. It does not test the KFRE in any real population.

## Question

When health systems differ in WHO gets a urine albumin (uACR) test, does the 4-variable Kidney
Failure Risk Equation (KFRE) look miscalibrated among the tested patients in a way that tracks each
site's testing rate? And can a pre-declared site-level analysis (random-effects meta-regression of
log observed/expected on the site testing rate) detect that pattern while staying quiet when
testing is uninformative?

Why it matters: in Stockholm's SCREAM data only 26,229 of 116,158 adults with CKD G3-4 had a uACR
(Cleary et al., BMJ Med 2026, PMID 41658049), and the 4-variable KFRE there had O/E 2.11. The
accompanying editorial (Walker & Lees, BMJ Med 2026, PMID 41658048) describes "underestimation in
patients with uACR data and overestimation in those without such data" and suggests that unmeasured
factors drive both testing and risk. If that is the mechanism, a model's calibration at a new
health system depends on that system's testing practice, which is exactly the kind of cross-site
instability the NIH NCATS Real-World Evidence Innovation Challenge asks entrants to measure.

## Generator (fixed in code, kfretransport/sim.py, version 0.1.0)

Per patient: age ~ N(70, 10) clipped to [30, 95]; male ~ Bernoulli(0.5); eGFR ~ Uniform(15, 45);
latent severity U ~ N(0, 1), not recorded; ln uACR (mg/g) = 4.0 + 1.0 U + N(0, 1.2^2).
True 5-year risk = KFRE 4-variable linear predictor with the TRUE uACR, plus g U, North American
5-year baseline 0.9240. Event ~ Bernoulli(true risk). Tested with probability expit(a_s + d U),
a_s = logit of the site's target rate; 8 sites with targets 10%, 20%, ..., 80%; 10,000 patients per
site; 200 replicates per scenario; seed 20261116.

| Scenario | d (testing follows U) | g (risk beyond KFRE inputs) |
|---|---|---|
| A_null | 0 | 0 |
| B_hidden_risk_random_testing | 0 | 0.5 |
| C_informative_testing | 1.0 | 0.5 |
| D_testing_follows_acr_only | 1.0 | 0 |

## Analysis (fixed in code)

Per site: O/E among the tested (complete case); O/E among the untested with uACR imputed at the
site's geometric mean of tested values; oracle O/E (everyone, true uACR). Across sites (sites with
at least one event among the tested): y = log O/E (tested), variance 1/O; DerSimonian-Laird pooling
with I^2 and a 95% prediction interval; random-effects meta-regression of y on the realized site
testing rate (method-of-moments tau^2, Knapp-Hartung t test with k-2 df, two-sided 5%).
Aggregate O/E for a scenario = (sum of O over sites) / (sum of E over sites), median over replicates.

## Predictions (graded exactly as written; sim.json "grades" is a convenience copy)

- **S1** In A, the meta-regression rejects (slope = 0) in at most 8% of analysed replicates.
- **S2** In B, it rejects in at most 8%.
- **S3** In D, it rejects in at most 8%.
- **S4** In C, it rejects in at least 50%, and the median slope is negative (sites that test fewer
  people show higher O/E among the tested).
- **S5** In C, the median aggregate O/E among the tested is above 1 AND the median aggregate O/E
  among the untested (imputed uACR) is below 1 (the pattern the editorial describes).
- **S6** In A and in D, the median aggregate O/E among the tested lies in [0.9, 1.1] (selection on
  U alone does not miscalibrate the tested when U carries no risk beyond the model's inputs).
- **S7** In B, the median aggregate O/E among the tested AND the median oracle O/E are both above 1.1.

## What each outcome would mean

- S1-S3 pass and S4 passes: the site-level moderator analysis is a usable, error-controlled check
  to pre-declare in a multi-site protocol (e.g. an NCATS RWE Challenge entry).
- S4 fails: with 8 sites the test is under-powered at this effect size; a protocol must then use
  more sites or a patient-level model of testing, and must say so.
- S1, S2 or S3 fails: the test raises false alarms; it must not be used as specified.
- S6 fails in D: the generator or the KFRE code is wrong; nothing else is interpreted until fixed.

## Limits stated before the run

One simple generator; one effect size; testing depends only on U; no competing risk of death; no
censoring; uniform eGFR 15-45; O/E is the only calibration measure (no calibration slope, to keep
the run small); variance of log O/E approximated by 1/O. Passing says the method behaves as
intended in this synthetic world only. Real-data behaviour is unknown.
