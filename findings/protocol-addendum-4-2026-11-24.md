# Protocol addendum 4 (2026-11-24): the held-out transport test

Status: **frozen on 2026-11-24, before any real health-system data has been seen.** No site
summary, real or synthetic, has been fitted with this test. Later changes go in new, dated
addenda; this file is not edited.

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it yet. Research only; information only, not medical advice, not a medical device.
Text CC BY 4.0; code (`kfretransport/transport.py`, release 0.2.7) Apache-2.0.

## Why this addendum

Addenda 1 and 2 test whether testing rate and KFRE calibration are *associated* across sites.
The practical claim is stronger: that a health system's uACR testing rate **predicts**, before
any outcome is seen there, how far the KFRE must be recalibrated for its tested patients (the
"transport formula"). A fitted slope cannot show that. Prediction of sites the fit has not
seen can. This addendum fixes that test in advance.

## Data and model

- The same site summaries and inclusion rules as the primary analysis (addendum 1): 5-year
  horizon, observed risk by Aalen-Johansen with death competing (`tte_aj`), y_k = log O/E among
  the tested, v_k = `se_log_oe_dm` squared, p_k = the share of the cohort with a measured uACR in
  the 365 days before index.
- Curve model (addendum 2): y_k = beta0 + log[Phi(Phi^-1(p_k) + theta) / p_k] + e_k,
  e_k ~ N(0, v_k + tau2), fitted by `selection.fit_sites` (theta grid -1 to 3, step 0.001).
- Null model: y_k = beta0 + e_k (theta = 0), the random-effects mean.
- tau2: the DerSimonian-Laird between-site variance of the null model on all included sites
  (`metrics.py`), used for **both** models, so the weights favour neither. If the curve explains
  part of the heterogeneity, this makes the curve's intervals conservative. That is stated in
  advance and accepted.

## The test (`transport.leave_one_out`)

Each site in turn is left out, both models are fitted on the others, and the left-out site's
log O/E is predicted from its testing rate (curve) or without it (null). For the curve, the
prediction variance is the linearised fit variance plus tau2 plus the site's own v_k.

- **T1 (primary for this addendum): skill > 0 and the curve wins at more than half the sites.**
  skill = 1 - SSE_curve / SSE_null, where SSE weights each site's squared prediction error by
  1 / (v_k + tau2). "Wins" means |error_curve| < |error_null|. Both conditions must hold.
- **T2: honest intervals.** At least 85% of left-out sites fall inside the curve's 95%
  prediction interval.
- **T3 (descriptive): stability.** The share of sites whose leave-one-out theta lies inside the
  all-site 95% profile interval of theta, and the largest change in theta from dropping one site.

## Prospective form (Phase 2)

Data partners that join after a fit has been published are predicted **before** their outcomes
are pooled: the prediction (log O/E and 95% interval from `transport.predict`) is committed to the
public repository with a timestamp, and graded the same way once their summaries arrive.

## What would make it fail, said in advance

- If theta is near 0, the curve spends a parameter for nothing and skill will be about 0 or
  negative. T1 then fails, which is the correct answer.
- Sites with the lowest testing rates carry the most information about theta (0.2.4). Leaving one
  out makes the fit extrapolate, so their prediction errors are the largest. This is part of the
  test, not a reason to drop them.
- A Gaussian selection model may misfit even when testing drives calibration. Then T2 fails
  before T1 does.
- Fewer than 6 included sites: the test is reported as not estimable and is not graded.

## Place in the protocol

Secondary to addendum 1's primary meta-regression. A pass here cannot rescue a failed primary,
and a failure here does not overturn a passed primary: it would mean the association is real but
not yet useful for prediction.
