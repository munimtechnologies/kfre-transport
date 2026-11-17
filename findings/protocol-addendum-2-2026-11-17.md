# Protocol addendum 2: selection-curve secondary analysis (2026-11-17)

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it. Text CC BY 4.0. Information only, not medical advice. Written before any real
health-system data, any multi-site result and any result from the simulation or NHANES studies
existed. Addendum 1 (the primary analysis) is unchanged. The earlier pre-registrations are
unchanged.

## What is added

A pre-specified **secondary** analysis for multi-site runs of `omop.py` and `pool.py`, using
`kfretransport.selection.fit_sites` (v0.2.3), described in
[selection-curve-2026-11-17.md](selection-curve-2026-11-17.md).

- Input: the primary site estimates from addendum 1. Each site gives its 5-year log O/E (O from
  Aalen-Johansen with death competing, E the mean KFRE), the variance `se_log_oe_dm^2`, and the
  site's measured-uACR testing rate p_k. Sites must have 0 < p_k < 1 and at least one KRT
  event.
- Model: log O/E_k = beta0 + log Phi(Phi^-1(p_k) + theta) - log p_k + error. The weights are
  1 / (variance_k + tau^2), where tau^2 is the DerSimonian-Laird estimate from the primary pooled
  fit. theta is fitted on the grid [-1, 3] with step 0.001.
- Reported: theta with its profile-likelihood 95% interval; beta0 and exp(beta0) (the O/E the
  sites would share if everyone were tested); and whether theta or the interval reaches the
  edge of the grid.
- Minimum: 3 sites with at least 2 distinct testing rates. Otherwise the analysis is reported as
  not estimable.

## Prediction

- **C1:** theta > 0, with the 95% interval excluding 0. This is the same direction as the
  negative slope predicted in addendum 1.
- C1 is secondary. It cannot rescue a failed primary, and a failed C1 does not overturn a
  passed primary.

## Exploratory (no prediction)

- Fit the same curve to the site-level output of simulation v1 (scenarios A-D). Then report how
  close the fitted theta comes to zero in scenarios A, B and D, and whether scenario C gives a
  positive theta. The v1 generator is not the Gaussian probit model, so this tests robustness,
  not the formula.
- Compare the weighted residual sum of squares of the curve with that of the linear
  meta-regression from addendum 1, fitted to the same sites.
