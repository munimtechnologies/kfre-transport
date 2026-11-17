# A closed-form selection-calibration curve for risk equations that need a lab test (2026-11-17)

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. (mv2@munimtech.com). No
human expert has reviewed this yet. Text CC BY 4.0; code (`kfretransport/selection.py`) Apache-2.0.
Research only. Information only, not medical advice, not a medical device.

## Question

The 4-variable Kidney Failure Risk Equation (KFRE) needs a urine albumin-to-creatinine ratio
(uACR). In routine care many people with CKD never get one: 63% of insured US adults with CKD
and diabetes had albuminuria testing, against 25% of those without (Chu et al., BMC Prim Care
2022, PMID 36434513). In Stockholm only 26,229 of 116,158 adults with CKD stages 3-4 had a uACR,
and among them the KFRE under-predicted kidney failure, observed/expected (O/E) 2.11 (Cleary et
al., BMJ Med 2026, PMID 41658049). **How much miscalibration can selective testing alone
produce, and how does it depend on a health system's testing rate?**

## Result

Let U ~ N(0, 1) be the part of a person's log-hazard of kidney failure that the equation's
predictors (age, sex, eGFR, uACR) do not explain, scaled so the hazard is multiplied by
exp(gU - g^2/2) (mean 1 over everyone). Let a person be tested when a standard-normal propensity
V = rho U + sqrt(1 - rho^2) e exceeds a threshold, with e independent noise, so that a share p of
the cohort is tested. Write z_p = Phi^-1(p) and theta = g * rho.

Under the measure tilted by exp(gU - g^2/2), U ~ N(g, 1) and V ~ N(theta, 1), so

    O/E among the tested   = Phi(z_p + theta) / p          (relative to the whole cohort's O/E)
    O/E among the untested = (1 - Phi(z_p + theta)) / (1 - p)

These are exact for the expected hazard and hold for the risk when the outcome is rare.
`risk_scale_factors()` computes the exact risk-scale version by quadrature. For common
outcomes and theta > 0, the tested factor always lies between 1 and the closed form. The reason: risk
f(x) = 1 - S^x is concave in the hazard multiplier x, so f(x)/x falls as x rises, and testing
shifts weight toward larger x (monotone likelihood ratio). A unit test checks this property
at one set of values.

What follows from the formula:

1. **One number decides it.** Only theta matters: the covariance between the unexplained part of
   the log-hazard and the standardized testing propensity. Testing that depends only on the
   equation's own predictors (for example on diabetes, eGFR or the true uACR itself) does not
   miscalibrate the equation in the tested. Within each predictor pattern, the curve applies
   with that pattern's own testing rate.
2. **Less testing, more apparent under-prediction, and the effect is not linear.** For theta > 0,
   tested O/E falls monotonically toward 1 as p rises. log O/E is steepest at low testing rates.
   The linear meta-regression in protocol addendum 1 is therefore a first-order approximation;
   this curve is a one-parameter alternative. It also gives the mirror image in the untested
   and a whole-cohort O/E, which is what transport between health systems needs.
3. **The two groups always average out:** p x tested + (1 - p) x untested = 1.

## Table (hand-computed from printed normal tables, about +/-0.002; the code computes them exactly)

Tested O/E relative to the whole cohort:

| theta | p = 0.10 | 0.25 | 0.40 | 0.63 | 0.90 |
|---|---|---|---|---|---|
| 0.10 | 1.187 | 1.131 | 1.098 | 1.059 | 1.018 |
| 0.30 | 1.632 | 1.416 | 1.297 | 1.169 | 1.048 |
| 0.69 | 2.771 | 2.025 | 1.672 | 1.344 | 1.084 |

Untested (true risk relative to what an equation calibrated in the whole cohort predicts):

| theta | p = 0.10 | 0.25 | 0.40 | 0.63 | 0.90 |
|---|---|---|---|---|---|
| 0.10 | 0.979 | 0.956 | 0.935 | 0.900 | 0.836 |
| 0.30 | 0.930 | 0.861 | 0.802 | 0.713 | 0.569 |
| 0.69 | 0.803 | 0.658 | 0.552 | 0.415 | 0.243 |

Between testing rates of 25% and 63% (Chu et al.'s no-diabetes and diabetes rates), tested log
O/E falls by about 0.050 per 10 percentage points at theta = 0.3. That is an O/E ratio of about
0.95 per 10 points. At theta = 0.69 the fall is about 0.108 per 10 points (ratio about 0.90).

## What it says about Stockholm (an illustration, not a finding about SCREAM)

With p = 26,229 / 116,158 = 0.2258 and tested O/E 2.11, `theta_from_oe` gives
**theta = 0.694**. That is the selection strength needed if the *whole* 2.11 came from testing
and the KFRE were calibrated in the whole Stockholm cohort. It would imply that untested people
carry about 0.68 times the risk the equation predicts for them, so at the same predictor values
the tested would carry about 3.1 times the risk of the untested. If the whole cohort's O/E were
itself 1.3 (for other reasons, such as a different baseline hazard), the selection needed falls
to theta = 0.41. Cleary et al. do not report O/E in the untested, so this cannot be checked
from their abstract. Walker and Lees's editorial (PMID 41658048) describes over-estimation in
patients without uACR data, which is the direction the curve predicts.

## Why it may be new, and its limits

The mathematics is standard: exponential tilting of a bivariate normal, the same algebra as
Heckman-type selection models (Heckman, Econometrica 1979). We found no paper that uses it to
link a health system's **testing rate** to the **calibration of a risk equation** among the
tested. Our Europe PMC searches on 2026-11-17 combined informative or selective testing,
Heckman and inverse Mills with calibration, external validation and KFRE. Nor did we find one
that turns the link into a one-parameter transport model across sites. Please tell us if one
exists.

Limits:
- The model assumes Gaussian unexplained log-hazard and probit testing. Real selection may have
  heavier tails.
- It ignores competing death; the multi-site analysis handles that separately with
  Aalen-Johansen.
- theta is not observed directly: it has to be estimated from several sites (`fit_sites`), and
  sites also differ for other reasons. That is what the random-effects term and beta0 absorb.
- Nothing here has been run on real data. The curve is a model of a mechanism. It does not show
  the mechanism is at work.

## How to use

```
python -m kfretransport.selection --theta 0.1 0.3 0.69 --rates 0.1 0.25 0.4 0.63 0.9
```

`fit_sites(rates, log_oe, variances, tau2)` fits beta0 and theta across sites by weighted least
squares, with a profile-likelihood 95% interval for theta. Protocol addendum 2 makes it a
pre-specified secondary analysis.
