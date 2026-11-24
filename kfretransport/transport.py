# SPDX-License-Identifier: Apache-2.0
"""Transport formula: predict a NEW health system's KFRE calibration from its
uACR testing rate alone, and test that prediction on held-out systems (new in
0.2.7).

The selection-calibration curve (selection.py) says that among the tested

    log O/E_k = beta0 + log tested_factor(p_k, theta) + error_k,
    error_k ~ N(0, v_k + tau2),

where p_k is site k's testing rate, v_k the sampling variance of its log O/E
and tau2 the between-site variance the curve does not explain. Fitted on some
sites (selection.fit_sites), it predicts a new site's log O/E from p_new before
any outcome has been seen there.

Prediction variance = var_fit + tau2 + v_new. var_fit comes from the
Gauss-Newton (linearised) covariance of (beta0, theta) at the fitted theta, with
Jacobian rows (1, s_k), s_k = d/dtheta log tested_factor(p_k, theta)
= phi(z_k + theta) / Phi(z_k + theta), z_k = Phi^-1(p_k).

The comparator ("null") ignores testing: a random-effects mean of the other
sites (theta fixed at 0), prediction variance 1/sum(w) + tau2 + v_new.

leave_one_out() is the held-out test of addendum 4
(findings/protocol-addendum-4-2026-11-24.md): each site is predicted from the
others, and the curve's weighted squared prediction error is compared with the
null's (skill = 1 - SSE_curve / SSE_null), together with the coverage of the
curve's 95% prediction intervals.

tau2 is an input here (pre-specified or estimated elsewhere), not estimated.

Research code only. Information only, not medical advice, not a medical device.
Made with AI agents (Claude) working in MV2's name (MV2, Munim, Inc.).
"""

import argparse
import json
import math
from statistics import NormalDist

from kfretransport.selection import fit_sites, tested_factor

_N = NormalDist()
Z975 = _N.inv_cdf(0.975)


def _check_rate(p):
    if not (0.0 < p < 1.0):
        raise ValueError("testing rate p must be strictly between 0 and 1")


def theta_score(p, theta):
    """d/dtheta of log tested_factor(p, theta) = phi(z_p + theta) / Phi(z_p + theta)."""
    _check_rate(p)
    a = _N.inv_cdf(p) + theta
    return _N.pdf(a) / _N.cdf(a)


def _weights(variances, tau2):
    w = []
    for v in variances:
        if not v + tau2 > 0.0:
            raise ValueError("variance + tau2 must be positive")
        w.append(1.0 / (v + tau2))
    return w


def fit_covariance(rates, variances, theta, tau2=0.0):
    """Linearised covariance [[var b0, cov], [cov, var theta]] of (beta0, theta)
    at theta: the inverse of J^T W J, J_k = (1, theta_score(p_k, theta)),
    W = diag(1 / (v_k + tau2))."""
    if len(rates) != len(variances):
        raise ValueError("rates and variances must have the same length")
    w = _weights(variances, tau2)
    a = b = c = 0.0
    for p, wk in zip(rates, w):
        s = theta_score(p, theta)
        a += wk
        b += wk * s
        c += wk * s * s
    det = a * c - b * b
    if not det > 1e-12 * a * c:
        raise ValueError("no information about theta: the testing rates must differ")
    return [[c / det, -b / det], [-b / det, a / det]]


def predict(rates, log_oe, variances, p_new, v_new=0.0, tau2=0.0, step=0.001):
    """Fit the curve on the given sites and predict a new site's log O/E among
    its tested from its testing rate p_new. v_new is the new site's own
    sampling variance (0 for the true value)."""
    _check_rate(p_new)
    if v_new < 0.0 or tau2 < 0.0:
        raise ValueError("v_new and tau2 must be >= 0")
    fit = fit_sites(rates, log_oe, variances, tau2=tau2, step=step)
    theta = fit["theta"]
    cov = fit_covariance(rates, variances, theta, tau2)
    s_new = theta_score(p_new, theta)
    mean = fit["beta0"] + math.log(tested_factor(p_new, theta))
    var_fit = cov[0][0] + 2.0 * s_new * cov[0][1] + s_new * s_new * cov[1][1]
    var_pred = var_fit + tau2 + v_new
    sd = math.sqrt(var_pred)
    return {
        "p_new": p_new,
        "log_oe": mean,
        "oe": math.exp(mean),
        "var_fit": var_fit,
        "var_pred": var_pred,
        "pi95": [math.exp(mean - Z975 * sd), math.exp(mean + Z975 * sd)],
        "theta": theta,
        "beta0": fit["beta0"],
        "theta_at_grid_edge": fit["theta_at_grid_edge"],
    }


def predict_null(log_oe, variances, v_new=0.0, tau2=0.0):
    """Prediction that ignores testing: weighted mean of the sites' log O/E."""
    if len(log_oe) != len(variances) or not log_oe:
        raise ValueError("need log_oe and variances of the same, non-zero length")
    if v_new < 0.0 or tau2 < 0.0:
        raise ValueError("v_new and tau2 must be >= 0")
    w = _weights(variances, tau2)
    sw = sum(w)
    mean = sum(wk * y for wk, y in zip(w, log_oe)) / sw
    var_fit = 1.0 / sw
    var_pred = var_fit + tau2 + v_new
    sd = math.sqrt(var_pred)
    return {
        "log_oe": mean,
        "oe": math.exp(mean),
        "var_fit": var_fit,
        "var_pred": var_pred,
        "pi95": [math.exp(mean - Z975 * sd), math.exp(mean + Z975 * sd)],
    }


def leave_one_out(rates, log_oe, variances, tau2=0.0, step=0.001):
    """Held-out test: predict each site from all the others, by the curve and
    by the null. Errors are observed minus predicted log O/E; each site's
    squared error is weighted by 1 / (v_k + tau2)."""
    k = len(rates)
    if k != len(log_oe) or k != len(variances):
        raise ValueError("rates, log_oe and variances must have the same length")
    if k < 4:
        raise ValueError("need at least 4 sites (3 left to fit the curve)")
    w = _weights(variances, tau2)
    rows = []
    sse_c = sse_n = 0.0
    covered = 0
    for i in range(k):
        r = rates[:i] + rates[i + 1:]
        y = log_oe[:i] + log_oe[i + 1:]
        v = variances[:i] + variances[i + 1:]
        cur = predict(r, y, v, rates[i], v_new=variances[i], tau2=tau2, step=step)
        nul = predict_null(y, v, v_new=variances[i], tau2=tau2)
        e_c = log_oe[i] - cur["log_oe"]
        e_n = log_oe[i] - nul["log_oe"]
        z_c = e_c / math.sqrt(cur["var_pred"])
        z_n = e_n / math.sqrt(nul["var_pred"])
        sse_c += w[i] * e_c * e_c
        sse_n += w[i] * e_n * e_n
        if abs(z_c) <= Z975:
            covered += 1
        rows.append({
            "site": i,
            "rate": rates[i],
            "observed": log_oe[i],
            "pred_curve": cur["log_oe"],
            "pred_null": nul["log_oe"],
            "err_curve": e_c,
            "err_null": e_n,
            "z_curve": z_c,
            "z_null": z_n,
            "theta_without_site": cur["theta"],
            "theta_at_grid_edge": cur["theta_at_grid_edge"],
        })
    return {
        "n_sites": k,
        "tau2": tau2,
        "sse_curve": sse_c,
        "sse_null": sse_n,
        "skill": (1.0 - sse_c / sse_n) if sse_n > 0.0 else None,
        "coverage95_curve": covered / k,
        "mean_z_curve": sum(row["z_curve"] for row in rows) / k,
        "rows": rows,
    }


def _read_sites(path):
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    sites = data["sites"] if isinstance(data, dict) else data
    rates = [float(s["rate"]) for s in sites]
    log_oe = [float(s["log_oe"]) for s in sites]
    variances = [float(s["variance"]) for s in sites]
    return rates, log_oe, variances


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Transport formula and held-out test (research only, not medical advice).")
    ap.add_argument("--sites", required=True,
                    help='JSON: [{"rate": p, "log_oe": y, "variance": v}, ...] (site summaries)')
    ap.add_argument("--tau2", type=float, default=0.0)
    ap.add_argument("--predict", type=float, nargs="*", default=[],
                    help="testing rates of new sites to predict")
    ap.add_argument("--loo", action="store_true", help="run the leave-one-site-out test")
    args = ap.parse_args(argv)
    rates, log_oe, variances = _read_sites(args.sites)
    out = {"predictions": [predict(rates, log_oe, variances, p, tau2=args.tau2)
                           for p in args.predict]}
    if args.loo:
        out["leave_one_out"] = leave_one_out(rates, log_oe, variances, tau2=args.tau2)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
