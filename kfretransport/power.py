# SPDX-License-Identifier: Apache-2.0
"""How many health systems a multi-site study needs to see the selection effect.

Design calculations for the selection-calibration curve of selection.py:

    log O/E_k = beta0 + log Phi(Phi^-1(p_k) + theta) - log p_k + error_k,
    error_k ~ N(0, v_k + tau2),

where p_k is site k's uACR testing rate, v_k the sampling variance of its log
O/E and tau2 the between-site variance the curve does not explain.

- To first order in theta, log O/E_k = beta0 + theta * m(p_k), with
  m(p) = phi(Phi^-1(p)) / p, the inverse Mills ratio of the testing threshold
  (Heckman 1979). So the information about theta is the weighted spread of
  m(p_k) across sites, and sites that test few people carry the most.
- Curve test (likelihood ratio, theta = 0 against the curve, beta0 free):
  noncentrality = weighted centred sum of squares of log tested_factor(p_k, theta).
- Linear test (slope of log O/E on p, the primary analysis of protocol addendum 1,
  here with known variances): noncentrality = (S_pm)^2 / S_pp. By Cauchy-Schwarz it
  never exceeds the curve test's noncentrality.
- Power of a two-sided level-alpha test with noncentrality L:
  Phi(sqrt(L) - z) + Phi(-sqrt(L) - z), z = Phi^-1(1 - alpha / 2).

monte_carlo_power() checks these approximations by simulation, using the same
profile-likelihood grid as selection.fit_sites().

Research code only. Information only, not medical advice, not a medical device.
Made with AI agents (Claude) working in MV2's name (MV2, Munim, Inc.).
"""

import argparse
import json
import math
import random
from statistics import NormalDist

from kfretransport.selection import tested_factor

_N = NormalDist()


def _weights(rates, variances, tau2):
    if len(rates) != len(variances):
        raise ValueError("rates and variances must have the same length")
    if len(rates) < 3:
        raise ValueError("need at least 3 sites")
    if len(set(rates)) < 2:
        raise ValueError("need at least 2 distinct testing rates")
    for p in rates:
        if not (0.0 < p < 1.0):
            raise ValueError("testing rate p must be strictly between 0 and 1")
    weights = []
    for v in variances:
        if v + tau2 <= 0:
            raise ValueError("variance + tau2 must be positive")
        weights.append(1.0 / (v + tau2))
    return weights


def _centred_ss(weights, values):
    sw = sum(weights)
    mean = sum(w * x for w, x in zip(weights, values)) / sw
    return sum(w * (x - mean) ** 2 for w, x in zip(weights, values))


def mills(p, theta=0.0):
    """d/dtheta of log tested_factor(p, theta) = phi(z_p + theta) / Phi(z_p + theta).
    At theta = 0 this is phi(z_p) / p."""
    if not (0.0 < p < 1.0):
        raise ValueError("testing rate p must be strictly between 0 and 1")
    a = _N.inv_cdf(p) + theta
    return _N.pdf(a) / _N.cdf(a)


def theta_information(rates, variances, tau2=0.0, theta=0.0):
    """Fisher information for theta with beta0 profiled out (known variances)."""
    w = _weights(rates, variances, tau2)
    return _centred_ss(w, [mills(p, theta) for p in rates])


def theta_se(rates, variances, tau2=0.0, theta=0.0):
    """Large-sample standard error of the fitted theta: 1 / sqrt(information)."""
    return 1.0 / math.sqrt(theta_information(rates, variances, tau2, theta))


def curve_noncentrality(rates, variances, tau2, theta):
    w = _weights(rates, variances, tau2)
    return _centred_ss(w, [math.log(tested_factor(p, theta)) for p in rates])


def linear_noncentrality(rates, variances, tau2, theta):
    w = _weights(rates, variances, tau2)
    mu = [math.log(tested_factor(p, theta)) for p in rates]
    sw = sum(w)
    pbar = sum(wi * p for wi, p in zip(w, rates)) / sw
    sxx = sum(wi * (p - pbar) ** 2 for wi, p in zip(w, rates))
    sxy = sum(wi * (p - pbar) * m for wi, p, m in zip(w, rates, mu))
    return sxy * sxy / sxx


def power_from_noncentrality(nc, alpha=0.05):
    if not (0.0 < alpha < 1.0):
        raise ValueError("alpha must be strictly between 0 and 1")
    z = _N.inv_cdf(1.0 - alpha / 2.0)
    s = math.sqrt(max(nc, 0.0))
    return _N.cdf(s - z) + _N.cdf(-s - z)


def curve_power(rates, variances, tau2, theta, alpha=0.05):
    return power_from_noncentrality(curve_noncentrality(rates, variances, tau2, theta), alpha)


def linear_power(rates, variances, tau2, theta, alpha=0.05):
    return power_from_noncentrality(linear_noncentrality(rates, variances, tau2, theta), alpha)


def even_rates(k, lo, hi):
    """k testing rates spread evenly from lo to hi (both included)."""
    if k < 2:
        raise ValueError("need k >= 2")
    if not (0.0 < lo < hi < 1.0):
        raise ValueError("need 0 < lo < hi < 1")
    return [lo + (hi - lo) * i / (k - 1) for i in range(k)]


def sites_needed(theta, site_variance, tau2=0.0, lo=0.1, hi=0.7, target=0.8,
                 alpha=0.05, test="curve", kmax=2000):
    """Smallest number of equal-sized sites, testing rates evenly spread over
    [lo, hi], whose power reaches target. None if kmax sites are not enough."""
    if test == "curve":
        fn = curve_power
    elif test == "linear":
        fn = linear_power
    else:
        raise ValueError("test must be 'curve' or 'linear'")
    if not (0.0 < target < 1.0):
        raise ValueError("target must be strictly between 0 and 1")
    for k in range(3, kmax + 1):
        if fn(even_rates(k, lo, hi), [site_variance] * k, tau2, theta, alpha) >= target:
            return k
    return None


def monte_carlo_power(rates, variances, tau2, theta, beta0=0.0, reps=400, seed=1,
                      alpha=0.05, lo=-1.0, hi=3.0, step=0.01):
    """Share of simulated studies whose profile-likelihood 95% interval for
    theta (the fit_sites method, on the grid [lo, hi]) excludes 0."""
    w = _weights(rates, variances, tau2)
    if not (lo < 0.0 < hi) or step <= 0:
        raise ValueError("need lo < 0 < hi and step > 0")
    if reps < 1:
        raise ValueError("need reps >= 1")
    crit = _N.inv_cdf(1.0 - alpha / 2.0) ** 2
    n = int(round((hi - lo) / step))
    grid = [lo + i * step for i in range(n + 1)]
    i0 = min(range(len(grid)), key=lambda i: abs(grid[i]))
    sw = sum(w)
    rows = []
    for th in grid:
        row = [math.log(tested_factor(p, th)) for p in rates]
        rw = sum(wi * r for wi, r in zip(w, row))
        rww = sum(wi * r * r for wi, r in zip(w, row))
        rows.append((row, rw, rww))
    sd = [math.sqrt(v + tau2) for v in variances]
    mu = [beta0 + math.log(tested_factor(p, theta)) for p in rates]
    rng = random.Random(seed)
    hits = 0
    for _ in range(reps):
        y = [m + rng.gauss(0.0, s) for m, s in zip(mu, sd)]
        sy = sum(wi * yi for wi, yi in zip(w, y))
        syy = sum(wi * yi * yi for wi, yi in zip(w, y))
        qs = []
        for row, rw, rww in rows:
            cross = sum(wi * yi * r for wi, yi, r in zip(w, y, row))
            sr = sy - rw
            qs.append(syy - 2.0 * cross + rww - sr * sr / sw)
        if qs[i0] - min(qs) > crit:
            hits += 1
    return hits / reps


def design_table(thetas, site_variances, tau2s, spreads, target=0.8, alpha=0.05, kmax=2000):
    rows = []
    for lo, hi in spreads:
        for th in thetas:
            for v in site_variances:
                for t2 in tau2s:
                    rows.append({
                        "rates": [lo, hi],
                        "theta": th,
                        "site_variance": v,
                        "tau2": t2,
                        "sites_curve": sites_needed(th, v, t2, lo, hi, target, alpha, "curve", kmax),
                        "sites_linear": sites_needed(th, v, t2, lo, hi, target, alpha, "linear", kmax),
                    })
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description="Sites needed to detect the selection curve (research only).")
    ap.add_argument("--theta", type=float, nargs="+", default=[0.1, 0.2, 0.3, 0.5, 0.69])
    ap.add_argument("--site-variance", type=float, nargs="+", default=[0.005, 0.01, 0.02, 0.04])
    ap.add_argument("--tau2", type=float, nargs="+", default=[0.0, 0.02, 0.05, 0.1])
    ap.add_argument("--spread", type=float, nargs=2, action="append",
                    help="lowest and highest site testing rate (repeatable)")
    ap.add_argument("--target", type=float, default=0.8)
    ap.add_argument("--kmax", type=int, default=2000)
    args = ap.parse_args(argv)
    spreads = args.spread or [[0.1, 0.7], [0.2, 0.5]]
    print(json.dumps(design_table(args.theta, args.site_variance, args.tau2, spreads,
                                  args.target, 0.05, args.kmax), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
