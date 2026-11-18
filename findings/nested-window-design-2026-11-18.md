# Testing the selection curve inside one health system: the nested-window design (2026-11-18)

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it. Text CC BY 4.0; code Apache-2.0 (`kfretransport/nested.py`, version 0.2.5).
Information only, not medical advice. Every number below was worked out by hand and is marked
that way. The unit tests check the identities and hand values against the code. No real or
synthetic health-system data were used.

## The question this serves

kfre-transport asks whether a health system's urine albumin (uACR) testing practice decides
whether the Kidney Failure Risk Equation (KFRE) is accurate there. The selection curve of
[selection-curve-2026-11-17.md](selection-curve-2026-11-17.md) says that if testing follows risk
the equation's predictors do not capture, the tested people's O/E is Phi(z_p + theta) / p times
the whole cohort's, where p is the testing rate and z_p = Phi^-1(p). Comparing many systems tests
that, but systems differ in many other ways. [power-design-2026-11-17.md](power-design-2026-11-17.md)
shows how much that between-site variance (tau^2) costs: at theta = 0.3 and tau^2 = 0.05, about 58
sites of 50 events each. The 2026-11-17 evidence map found that published validations cannot
supply those sites.

## The design

Keep one system, one cohort, one equation, one baseline and one estimator. Change only the uACR
lookback window W. Let P_W be the share of the cohort with a measured uACR in the W days before
the index date. A person's latest value inside a window is the same in every longer window, so
the tested sets are nested and P_W rises with W. Score **disjoint recency bands**:

- band 1: latest uACR 0 to W_1 days before index;
- band j: latest uACR in (W_(j-1), W_j] days before index.

The bands hold different people, so their O/E estimates are independent. Nested sets would share
people and need a covariance term. Bands do not. The primary windows (90, 180 and 365 days)
split exactly the tested set of the primary multi-site analysis (365 days).

## The band curve

Use the model of the selection-curve note, with one latent testing propensity V for all windows
("threshold model"): a person is tested within W when V > c_W, and P_W = P(V > c_W). The hazard
multiplier is exp(gU - g^2/2), and V = rho U + sqrt(1 - rho^2) e. Exponential tilting shifts V by
theta = g rho, so the expected multiplier in the band of testing ranks (P_(j-1), P_j] is exactly
(on the hazard scale)

    band factor_j = [Phi(z_j + theta) - Phi(z_(j-1) + theta)] / (P_j - P_(j-1)),   z_0 = -inf.

Checks, all in `test_nested.py`:

- theta = 0 gives 1 for every band;
- band 1 is the tested-curve value Phi(z_1 + theta) / P_1;
- the share-weighted band factors add up to the tested-curve value at the widest window
  (the sum telescopes);
- for theta > 0 the factors fall from band 1 to band J; for theta < 0 they rise.

Within a site, each band's log O/E is beta0_k + log band factor_j + error. beta0_k (the O/E the
site would have if everyone were tested) is free for each site. So theta is identified only
by contrasts **inside** sites, and between-site differences in case mix, outcome capture or
baseline risk cancel. `nested.fit_within` fits one theta across any number of sites, with one
intercept per site, by the same profile-likelihood grid as `selection.fit_sites`.

## How much one system can see (hand-computed illustration)

The information about theta from one site is the weighted centred sum of squares of the band
scores d log(band factor)/d theta = [phi(z_j) - phi(z_(j-1))] / (P_j - P_(j-1)) at theta = 0,
with weights 1 / variance (about the number of KRT events in the band).

Illustration. The inputs are assumptions, not data: cumulative testing rates 0.15, 0.30 and
0.40 at the three windows, and 90, 60 and 30 KRT events in the three bands (180 in all).

| Band | Ranks | Score (by hand) | Events (assumed) |
|---|---|---|---|
| 1 | 0-0.15 | 0.233159 / 0.15 = 1.5544 | 90 |
| 2 | 0.15-0.30 | (0.347693 - 0.233159) / 0.15 = 0.7636 | 60 |
| 3 | 0.30-0.40 | (0.386343 - 0.347693) / 0.10 = 0.3865 | 30 |

Weighted mean score 1.0961; information 40.64; standard error of theta 0.157 (checked by
`test_illustration_by_hand`). Information grows in proportion to the number of events. For a
two-sided 5% test at 80% power (linearised at theta = 0), the required information is
(2.802 / theta)^2. That is 16.5 at theta = 0.69, about 73 events with this band pattern, and 87.2
at theta = 0.3, about 386 events. The cross-site design under tau^2 = 0.05 needs about 13 and 58
sites of 50 events (650 and 2,900 events) for the same power. One large system, or a few, could
test the curve. These are linearised hand values. They ignore the curve's bending at larger
theta, any extra within-site variance, and stale uACR values (below).

## What the shape says about how people are tested

The threshold model is not the only way testing can follow risk. In a **recency (Poisson)
model**, each person is tested at the times of a Poisson process whose rate is rate0 * exp(kappa
U). Band j then holds people whose most recent test was (W_(j-1), W_j] days back.
`nested.poisson_band_factors` computes these band factors by quadrature. Two facts follow:

- For a very short window, the tested group's factor tends to E[exp(gU - g^2/2) lambda] /
  E[lambda] = exp(g kappa), the same for every short window. The gradient is flat at small P.
  The unit test checks this limit: exp(0.3) = 1.3499 for g = 0.5, kappa = 0.6.
- The threshold curve is steepest at small P. As P goes to 0, Phi(z_P + theta) / P grows
  without bound for theta > 0.

So within one system, the shape of O/E across recency bands separates "some people are never
tested" (threshold) from "everyone is tested, the sick more often" (recency). The two mechanisms
imply different things. Under the threshold model, the bias in the tested group grows without
bound as fewer people are tested. Under the recency model, the bias among the most recently
tested tends to the finite value exp(g kappa), however short the window. Addendum 3 reports
which shape fits better as a check, not as a test.

## Limits

- **Stale values.** Older bands are scored with older uACR results. The direction of the bias
  is not known in advance. Rising albuminuria in people who progress would raise O/E in older
  bands, against the prediction. Falling albuminuria after treatment would do the opposite. The
  secondary windows (90/365/730) show how much this matters.
- **Mixing predictor patterns.** The curve is exact within one pattern of the equation's
  predictors. Pooled over patterns with different testing rates, it holds only approximately,
  the same caveat as the cross-site curve.
- **One latent propensity.** The threshold model assumes the same V decides testing in every
  window. The consistency check against the cross-site theta (addendum 3, N3) and the
  goodness-of-fit statistic test this assumption.
- **Hazard scale.** The band factors are exact for hazards and close for risks when the outcome
  is rare, as with the selection curve. `selection.risk_scale_factors` quantifies the gap for
  the tested curve.

## Prior art (searched 2026-11-17 in Europe PMC, as recorded in the evidence-map note)

No KFRE validation was found that compares calibration across uACR recency or across testing
windows within one system. The band curve is a direct consequence of the truncated-normal
selection model (Heckman 1979) applied to nested thresholds. It is a design for this study, not
a claimed new piece of statistics.
