# Protocol addendum 3: within-system nested-window test (2026-11-18)

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it. Text CC BY 4.0. Information only, not medical advice. Written and published
**before any real health-system data have been read by this project**, before any multi-site
result, and before any result from the queued simulation, smoke and NHANES jobs. Addenda 1 and 2
and the earlier pre-registrations are unchanged. Any later change will be a new, dated addendum.

The design, the model and the formulas are in
[nested-window-design-2026-11-18.md](nested-window-design-2026-11-18.md). The code is
`kfretransport.nested` (version 0.2.5; CLI `kfre-nested`).

## Why

The cross-site test (addenda 1 and 2) compares health systems, so anything else that differs
between them (case mix, how kidney failure is recorded, baseline risk) can mimic or hide the
testing-rate effect. The evidence map of 2026-11-17 found that published validations cannot
answer the question (2 qualifying studies of about 45 screened). This addendum adds a test that
runs **inside one system**, where those differences are held fixed: only the uACR lookback
window changes.

## Fixed choices

| Item | Fixed choice |
|---|---|
| Cohort | `omop.build_cohort` as in addendum 1 (unchanged); the cohort does not depend on the window |
| Windows (primary) | W = 90, 180 and 365 days before index. The three bands partition exactly the tested set of the primary analysis (365 days), so no uACR older than one year is used |
| Windows (secondary) | W = 90, 365 and 730 days (adds a band of people whose latest uACR is one to two years old) |
| Bands | band 1: latest uACR 0-90 days before index; band j: in (W_(j-1), W_j]. Disjoint people, so the band estimates are independent; no covariance term |
| Band estimate | log O/E at 5 years (primary; 2 years secondary): O = Aalen-Johansen cumulative incidence of KRT with death competing, E = mean KFRE of the band (North American baseline at US sites); variance `se_log_oe_dm^2` |
| Small bands | a band with fewer than 20 KRT events by the horizon is merged with the next band by dropping its upper window (the last band merges with the one before it); repeated until every band qualifies. A site needs at least 2 bands left (`nested.bands_from_rows`) |
| Model | log O/E_kj = beta0_k + log{[Phi(z_kj + theta) - Phi(z_k(j-1) + theta)] / (P_kj - P_k(j-1))}, z = Phi^-1(cumulative testing rate P), z_0 = -inf. One theta, one intercept per site. Weighted least squares with known variances; theta on the grid [-1, 3], step 0.001; 95% profile interval (chi2(1)) (`nested.fit_within`) |
| Disclosure | published per-band rows follow `nested.public_report`: counts from 1 to 19 and everything that would give them back are blanked |

## Predictions (stated now)

- **N1 (primary for this addendum):** theta_W > 0 with a 95% profile interval that excludes 0,
  fitted over all usable sites together (or at the one site, if only one runs). In words: within a
  system, people whose albumin test was more recent show more under-prediction by the KFRE, in the
  amount the selection curve gives.
- **N2 (descriptive):** the number of usable sites in which band 1 has a higher O/E than the
  last band.
- **N3 (consistency, needs the cross-site fit of addendum 2 with at least 3 sites):** z =
  (theta_W - theta_B) / sqrt(se_W^2 + se_B^2), with each se read off its 95% interval as
  half-width / 1.96 (`nested.consistency`). |z| <= 1.96 counts as consistent.

How the results will be read:

| N1 | N3 | Reading |
|---|---|---|
| holds | consistent | one selection parameter explains both within-system and between-system miscalibration; the transport formula of addendum 2 is supported |
| holds | inconsistent | selection is real inside systems, but between-system differences also have other causes |
| fails | consistent or not estimable | no within-system evidence for threshold-type selection; see the mechanism check below |

N1 is separate from the primary hypothesis of addendum 1. It cannot rescue a failed addendum-1
primary, and a failed N1 does not overturn a passed one. The secondary windows and the 2-year
horizon are reported next to the primary choices and cannot replace them.

## Pre-specified checks (no prediction)

- Goodness of fit: q_min against chi2 with (bands - sites - 1) degrees of freedom.
- **Mechanism.** If people are tested at a rate tied to their risk (the Poisson recency model in
  the design note) rather than once by a threshold, the band gradient is flat for short windows
  and theta_W comes out smaller than theta_B. The fitted band factors will be set next to both
  shapes, and the result will say which one fits better. It will not be called a test.
- **Stale values.** Older bands are scored with older uACR results. If albuminuria tends to rise
  over time in people who progress, stale values under-predict, which raises O/E in older bands
  and works **against** N1. Lowering it is also possible (for example, after treatment starts).
  The direction is not assumed; the secondary windows show how much it matters.

## Status on this date

Nothing has run on real data. The code is tested on invented data only (hand-worked unit tests).
