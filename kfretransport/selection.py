# SPDX-License-Identifier: Apache-2.0
"""Selection-calibration curve: how far selective uACR testing alone can make a
risk equation look miscalibrated in the people who were tested.

Model (stated in full in findings/selection-curve-2026-11-17.md):
- U ~ N(0, 1) is the part of a person's log-hazard of kidney failure that the
  equation's own predictors do not explain; the hazard is multiplied by
  exp(g U - g^2 / 2), so the multiplier averages 1 over everyone.
- A person is tested when V = rho U + sqrt(1 - rho^2) e > c, e ~ N(0, 1)
  independent of U, and the share tested is p = P(V > c), so c = -z_p with
  z_p = Phi^-1(p).
- Exponential tilting gives, exactly on the hazard scale (and for the risk
  when the outcome is rare):

      O/E among the tested   = Phi(z_p + theta) / p
      O/E among the untested = (1 - Phi(z_p + theta)) / (1 - p)

  with theta = g * rho, relative to the O/E of the whole cohort. Only the
  product theta matters. Selection on the equation's own predictors does not
  enter: within each predictor pattern the curve applies with that pattern's p.

Research code only. Information only, not medical advice, not a medical device.
Made with AI agents (Claude) working in MV2's name (MV2, Munim, Inc.).
"""

import argparse
import json
import math
from statistics import NormalDist

_N = NormalDist()
CHI2_1_95 = _N.inv_cdf(0.975) ** 2


def _check_rate(p):
    if not (0.0 < p < 1.0):
        raise ValueError("testing rate p must be strictly between 0 and 1")


def tested_factor(p, theta):
    """O/E among the tested divided by O/E in the whole cohort: Phi(z_p + theta) / p."""
    _check_rate(p)
    return _N.cdf(_N.inv_cdf(p) + theta) / p


def untested_factor(p, theta):
    """O/E among the untested divided by O/E in the whole cohort."""
    _check_rate(p)
    return (1.0 - _N.cdf(_N.inv_cdf(p) + theta)) / (1.0 - p)


def theta_from_oe(p, oe_ratio):
    """The theta at which the curve gives tested O/E = oe_ratio at testing rate p
    (whole-cohort O/E taken as 1): Phi^-1(oe_ratio * p) - Phi^-1(p)."""
    _check_rate(p)
    q = oe_ratio * p
    if not (0.0 < q < 1.0):
        raise ValueError("oe_ratio * p must be strictly between 0 and 1")
    return _N.inv_cdf(q) - _N.inv_cdf(p)


def log_slope(p, theta):
    """d/dp of log tested_factor(p, theta): phi(z+theta) / (Phi(z+theta) phi(z)) - 1/p."""
    _check_rate(p)
    z = _N.inv_cdf(p)
    a = z + theta
    return _N.pdf(a) / (_N.cdf(a) * _N.pdf(z)) - 1.0 / p


def risk_scale_factors(p, g, rho, base_risk, half_width=8.0, steps=4000):
    """Exact tested and untested factors on the RISK scale, by quadrature.

    Each person's risk is 1 - (1 - base_risk) ** exp(g u - g^2 / 2). The
    reference E is the mean risk over everyone, i.e. what an equation
    calibrated in the whole cohort predicts for this predictor pattern.
    Returns tested, untested, tested_share (the quadrature's own p) and
    closed_form (the rare-outcome curve tested_factor(p, g * rho)).
    """
    _check_rate(p)
    if not (0.0 < base_risk < 1.0):
        raise ValueError("base_risk must be strictly between 0 and 1")
    if not (-1.0 < rho < 1.0):
        raise ValueError("rho must be strictly between -1 and 1")
    if steps < 2 or half_width <= 0:
        raise ValueError("need steps >= 2 and half_width > 0")
    z = _N.inv_cdf(p)
    s = math.sqrt(1.0 - rho * rho)
    log_surv = math.log1p(-base_risk)
    h = 2.0 * half_width / steps
    mass = 0.0
    mass_t = 0.0
    risk_all = 0.0
    risk_t = 0.0
    for i in range(steps + 1):
        u = -half_width + i * h
        w = 0.5 if i == 0 or i == steps else 1.0
        dens = w * math.exp(-0.5 * u * u)
        mult = math.exp(g * u - 0.5 * g * g)
        risk = -math.expm1(mult * log_surv)
        pt = _N.cdf((z + rho * u) / s)
        mass += dens
        mass_t += dens * pt
        risk_all += dens * risk
        risk_t += dens * risk * pt
    mean_all = risk_all / mass
    tested = (risk_t / mass_t) / mean_all
    untested = ((risk_all - risk_t) / (mass - mass_t)) / mean_all
    return {
        "tested": tested,
        "untested": untested,
        "tested_share": mass_t / mass,
        "mean_risk": mean_all,
        "closed_form": tested_factor(p, g * rho),
    }


def fit_sites(rates, log_oe, variances, tau2=0.0, lo=-1.0, hi=3.0, step=0.001):
    """Fit log O/E_k = beta0 + log tested_factor(p_k, theta) across sites.

    Weighted least squares with weights 1 / (variance_k + tau2), beta0 profiled
    out in closed form, theta on a grid. With known variances this is the
    normal likelihood, so the 95% interval is the set of grid thetas whose
    weighted residual sum of squares is within chi2(1, 0.95) of the minimum.
    """
    k = len(rates)
    if k != len(log_oe) or k != len(variances):
        raise ValueError("rates, log_oe and variances must have the same length")
    if k < 3:
        raise ValueError("need at least 3 sites")
    if len(set(rates)) < 2:
        raise ValueError("need at least 2 distinct testing rates")
    if step <= 0 or hi <= lo:
        raise ValueError("need step > 0 and hi > lo")
    for p in rates:
        _check_rate(p)
    weights = []
    for v in variances:
        if v + tau2 <= 0:
            raise ValueError("variance + tau2 must be positive")
        weights.append(1.0 / (v + tau2))
    sw = sum(weights)
    n = int(round((hi - lo) / step))
    profile = []
    for i in range(n + 1):
        th = lo + i * step
        resid = [y - math.log(tested_factor(p, th)) for p, y in zip(rates, log_oe)]
        b0 = sum(w * r for w, r in zip(weights, resid)) / sw
        q = sum(w * (r - b0) ** 2 for w, r in zip(weights, resid))
        profile.append((th, b0, q))
    best = min(profile, key=lambda t: t[2])
    inside = [t[0] for t in profile if t[2] <= best[2] + CHI2_1_95]
    ci = [min(inside), max(inside)]
    edge = step / 2.0
    return {
        "theta": best[0],
        "beta0": best[1],
        "oe_all": math.exp(best[1]),
        "q_min": best[2],
        "ci95_profile": ci,
        "theta_at_grid_edge": best[0] <= lo + edge or best[0] >= hi - edge,
        "ci_at_grid_edge": ci[0] <= lo + edge or ci[1] >= hi - edge,
        "k": k,
        "tau2": tau2,
    }


def curve_table(thetas, rates):
    rows = []
    for th in thetas:
        for p in rates:
            rows.append({
                "theta": th,
                "p": p,
                "tested": round(tested_factor(p, th), 4),
                "untested": round(untested_factor(p, th), 4),
            })
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description="Selection-calibration curve (research only).")
    ap.add_argument("--theta", type=float, nargs="+", default=[0.1, 0.3, 0.69])
    ap.add_argument("--rates", type=float, nargs="+", default=[0.1, 0.25, 0.4, 0.63, 0.9])
    args = ap.parse_args(argv)
    print(json.dumps(curve_table(args.theta, args.rates), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
