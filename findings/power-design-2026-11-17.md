# How many health systems does it take to see the testing effect? (2026-11-17)

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it. Text CC BY 4.0, code Apache-2.0. Information only, not medical advice. This is a
design calculation. It uses no patient data and makes no claim that the effect exists.

## The question

[selection-curve-2026-11-17.md](selection-curve-2026-11-17.md) gives a one-parameter curve.
Under selective uACR testing, a site that tests a share p of its eligible patients sees,
among the tested,

    log O/E = beta0 + log Phi(Phi^-1(p) + theta) - log p.

A multi-site study (protocol addenda 1 and 2) asks whether theta > 0. Before anyone spends
an N3C data-use agreement on it, the question is: **how many sites, of what size, with
what spread of testing rates, are needed to detect a given theta?**
`kfretransport/power.py` (0.2.4) answers it. This note gives the reasoning and
hand-computed values. The exact table comes from the job `kfre-power-v1`: [RESULT SLOT].

## Result 1: the right covariate is the inverse Mills ratio, not the testing rate

Expanding the curve to first order in theta gives

    log O/E_k = beta0 + theta * m(p_k) + ...,   m(p) = phi(Phi^-1(p)) / p,

where m is the inverse Mills ratio at the testing threshold (Heckman 1979, Econometrica
47:153). With known site variances, and with beta0 profiled out, the information about
theta is

    I = sum_k w_k (m(p_k) - m_bar_w)^2,   w_k = 1 / (v_k + tau^2).

Three consequences:

- **Sites that test few people carry the information.** m rises steeply as p falls. By hand
  (normal tables): m(0.7) = 0.497, m(0.5) = 0.798, m(0.3) = 1.159, m(0.1) = 1.755. A
  network of sites that all test 40-60% of patients can barely see theta, however many
  sites it has.
- **The spread of testing rates matters as much as the number of sites.** Simpson's rule
  over 13 points gives the variance of m(p) for p uniform on [0.1, 0.7] as 0.1160. For p
  uniform on [0.2, 0.5] it is 0.0293. The narrower network needs about **4 times as many
  sites** for the same power.
- **The pre-registered linear primary loses little here.** Addendum 1 regresses log O/E on
  p itself. Its efficiency relative to the curve is the squared correlation between m(p)
  and p. For p uniform on [0.1, 0.7] that is 0.98 by hand (Cov = -0.0584, Var p = 0.03,
  Var m = 0.1160). The loss grows when some sites test under 10%. In every design,
  Cauchy-Schwarz guarantees that the linear test's noncentrality never exceeds the curve
  test's. `test_power.py` checks this.

## Result 2: sites needed for 80% power (two-sided 5%)

Linearised, for K equal-sized sites with testing rates spread evenly over [lo, hi]:

    K ~ (1.960 + 0.842)^2 (v + tau^2) / (theta^2 Var m) = 7.849 (v + tau^2) / (theta^2 Var m).

Here v is one site's sampling variance of log O/E. It is roughly 1/(number of kidney-failure
events), so v = 0.02 means about 50 events per site. tau^2 is the between-site variance that
testing does **not** explain. The curve bends, which costs information at larger theta. By
hand, from the curve's exact values at p = 0.1 and p = 0.7, the cost is about 8% at
theta = 0.3, 14% at 0.5 and 19% at 0.69 (for 0.2 and 0.1, 5.5% and 3% are interpolated).
The values below include that cost.

Hand values, rates spread over [0.1, 0.7], v = 0.02 (about 50 events per site), rounded up:

| theta | tau^2 = 0 | tau^2 = 0.02 | tau^2 = 0.05 | tau^2 = 0.10 |
|---|---|---|---|---|
| 0.69 (all of SCREAM's O/E 2.11 is selection) | 4 | 8 | 13 | 22 |
| 0.5 | 7 | 13 | 22 | 38 |
| 0.3 | 17 | 33 | 58 | 99 |
| 0.2 | 36 | 72 | 126 | 215 |
| 0.1 | 140 | 279 | 488 | 836 |

For rates spread over [0.2, 0.5], multiply by about 4. Evenly spaced sites do a few percent
better than these continuous-spread values. Exact values for every cell, and for
v = 0.005-0.04 and the spread [0.05, 0.9], come from `python -m kfretransport.power`:
[RESULT SLOT, kfre-power-v1].

What this means for a study in N3C, which had 84 contributing sites in August 2024
("As of August 2024, N3C contains data on over 22.9 million patients from 84 contributing
sites": Walters et al., JAMIA 2025;32(2):391-397, PMID 39579089; the number may have changed
since):

- If testing explains a large share of SCREAM-sized miscalibration (theta 0.5 or more),
  20-40 sites are enough even with substantial residual heterogeneity.
- theta = 0.3 is detectable only if residual heterogeneity is modest (tau^2 up to about
  0.05) **and** the network includes sites that test 10-20% of eligible patients.
- theta = 0.1 is out of reach for any realistic network. A null result from 84 sites
  would therefore bound theta, not rule it out. The protocol should report the interval,
  never "no effect".
- Once sites have about 50 events, residual heterogeneity dominates. Adding events per site
  helps less than reducing tau^2. Candidate reductions: adjusting for case mix, stratifying
  by region, and the NHANES-based prior on theta.

## Checks

- `test_power.py` uses only exact identities and hand values: power equals alpha at
  theta = 0; m is the derivative of the log curve (finite differences); a two-point design
  gives information 39.58 by hand; linear noncentrality <= curve noncentrality; the small-theta
  noncentrality equals theta^2 times the information; the search for K returns the first K
  that reaches the target power. It also checks that a seeded Monte Carlo run of the actual
  profile-likelihood test (the `fit_sites` method) gives about 5% rejections at theta = 0 and
  agrees with the analytic power within 0.12 for one design.
- Larger Monte Carlo checks (2,000 replicates at six design points, including 84 sites) are
  in the job `kfre-power-v1`: [RESULT SLOT].

## Limits

- Known site variances and normal errors on log O/E. Small sites with few events break
  this. The real analysis uses the delta-method variance from `metrics.py`.
- The Gaussian-probit selection model of the curve. Any other selection mechanism changes
  the shape. The linear primary does not depend on it.
- The linear test here is a z-test with known variances. Addendum 1's Knapp-Hartung
  inference is a little less powerful with few sites.
- Equal-sized sites with evenly spread testing rates. Real networks are lumpier. Pass the
  real rates and variances to `curve_power()`.
- The tau^2 grid is a set of scenarios, not an estimate. No published source gives the
  residual between-site variance of KFRE log O/E after accounting for testing.
