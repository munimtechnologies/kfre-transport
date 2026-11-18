# SPDX-License-Identifier: Apache-2.0
"""Within-system nested-window design (new in 0.2.5).

One health system, one cohort, one baseline and one method. Only the uACR
lookback window W changes. Let P_W be the share of the cohort with a measured
uACR in the W days up to the index date. P_W rises with W. The person's latest
value inside the window is the same in every window that contains it, so the
tested sets are nested. The design scores DISJOINT recency bands:

    band 1: latest uACR 0..W_1 days before index
    band j: latest uACR in (W_(j-1), W_j] days before index

Disjoint bands hold different people, so their log O/E estimates are
independent. No covariance term is needed (nested sets would share people).

Threshold model (selection.py, one latent testing propensity V for all
windows): tested within W  <=>  V > c_W, with P_W = P(V > c_W). Exponential
tilting then gives, exactly on the hazard scale and relative to the whole
cohort's O/E,

    band factor_j = [Phi(z_j + theta) - Phi(z_(j-1) + theta)] / (P_j - P_(j-1))

with z_j = Phi^-1(P_j), z_0 = -inf, and the same theta = g * rho as the
cross-site curve. Band 1 is tested_factor(P_1, theta). For theta > 0 the band
factors fall from band 1 to band J.

Recency (Poisson) model, for contrast: each person is tested at times of a
Poisson process with rate r0 * exp(kappa U). The band of a person is the one
holding their most recent test. Then the band factors follow no closed form
(poisson_band_factors(), by quadrature). For very short windows the tested
group's factor tends to exp(g * kappa) whatever the window, so the gradient is
flat at small P. The threshold model's gradient is steepest at small P. The
shape of the within-site gradient therefore says which kind of selection a
system has.

Pre-registered use: findings/protocol-addendum-3-2026-11-18.md.

Research code only. Information only, not medical advice, not a medical device.
Made with AI agents (Claude) working in MV2's name (MV2, Munim, Inc.).
"""

import argparse
import json
import math
from statistics import NormalDist

from kfretransport.metrics import oe_time_to_event
from kfretransport.omop import HORIZON_DAYS, build_cohort, read_cdm

_N = NormalDist()
_SQRT2 = math.sqrt(2.0)
_SQRT2PI = math.sqrt(2.0 * math.pi)
Z975 = _N.inv_cdf(0.975)
CHI2_1_95 = Z975 ** 2
PRIMARY_WINDOWS = (90, 180, 365)
SECONDARY_WINDOWS = (90, 365, 730)
MIN_BAND_EVENTS = 20
MIN_CELL = 20


# ---------------------------------------------------------------- the curve
def _phi_cdf(x):
    """Standard normal CDF through erfc (accurate in the lower tail; +-inf ok)."""
    return 0.5 * math.erfc(-x / _SQRT2)


def _phi_pdf(x):
    if math.isinf(x):
        return 0.0
    return math.exp(-0.5 * x * x) / _SQRT2PI


def _edge_z(p):
    if p <= 0.0:
        return -math.inf
    if p >= 1.0:
        return math.inf
    return _N.inv_cdf(p)


def _mass(za, zb, theta):
    """Phi(zb + theta) - Phi(za + theta), using upper tails when both are > 0."""
    a = za + theta
    b = zb + theta
    if a > 0.0:
        return _phi_cdf(-a) - _phi_cdf(-b)
    return _phi_cdf(b) - _phi_cdf(a)


def _check_band(lo, hi):
    if not (0.0 <= lo < hi <= 1.0):
        raise ValueError("a band needs 0 <= lower rate < upper rate <= 1")


def bands_of(cum_rates):
    """[(P_(j-1), P_j)] from cumulative testing rates P_1 < ... < P_J, P_0 = 0."""
    if not cum_rates:
        raise ValueError("need at least one cumulative testing rate")
    out = []
    lo = 0.0
    for hi in cum_rates:
        _check_band(lo, hi)
        out.append((lo, hi))
        lo = hi
    return out


def band_factor(lo, hi, theta):
    """O/E in the band of testing ranks (lo, hi], relative to the whole cohort."""
    _check_band(lo, hi)
    return _mass(_edge_z(lo), _edge_z(hi), theta) / (hi - lo)


def band_factors(cum_rates, theta):
    return [band_factor(lo, hi, theta) for lo, hi in bands_of(cum_rates)]


def band_score(lo, hi, theta=0.0):
    """d/dtheta of log band_factor(lo, hi, theta):
    [phi(z_hi + theta) - phi(z_lo + theta)] / [Phi(z_hi + theta) - Phi(z_lo + theta)].
    With lo = 0 this is the inverse Mills ratio of power.mills()."""
    _check_band(lo, hi)
    za, zb = _edge_z(lo), _edge_z(hi)
    den = _mass(za, zb, theta)
    if den <= 0.0:
        raise ValueError("band has no probability mass at this theta")
    return (_phi_pdf(zb + theta) - _phi_pdf(za + theta)) / den


# ---------------------------------------------------------------- fitting
def _prepare(site, need_y=True):
    rates = list(site["cum_rates"])
    v = list(site["variances"])
    y = list(site["log_oe"]) if need_y else [0.0] * len(rates)
    if not (len(rates) == len(v) == len(y)):
        raise ValueError("cum_rates, log_oe and variances must have the same length")
    if len(rates) < 2:
        raise ValueError("each site needs at least 2 bands")
    bands = bands_of(rates)
    w = []
    for vi in v:
        if not vi > 0.0:
            raise ValueError("every band variance must be > 0")
        w.append(1.0 / vi)
    edges = [_edge_z(0.0)] + [_edge_z(p) for p in rates]
    widths = [hi - lo for lo, hi in bands]
    return edges, widths, y, w, bands


def fit_within(sites, lo=-1.0, hi=3.0, step=0.001):
    """Common theta, one intercept beta0_k per site (profiled out in closed
    form), weights 1 / variance:

        log O/E_kj = beta0_k + log band_factor(P_k(j-1), P_kj, theta) + error

    theta on a grid; the 95% interval is the set of grid values whose weighted
    residual sum of squares is within chi2(1, 0.95) of the minimum (the
    fit_sites method). Only within-site contrasts inform theta."""
    if not sites:
        raise ValueError("need at least one site")
    if step <= 0 or hi <= lo:
        raise ValueError("need step > 0 and hi > lo")
    prepared = [_prepare(s) for s in sites]
    n = int(round((hi - lo) / step))
    profile = []
    for i in range(n + 1):
        th = lo + i * step
        q = 0.0
        b0s = []
        for edges, widths, y, w, _ in prepared:
            resid = []
            for j in range(len(y)):
                m = _mass(edges[j], edges[j + 1], th)
                if m <= 0.0:
                    resid = None
                    break
                resid.append(y[j] - math.log(m / widths[j]))
            if resid is None:
                q = math.inf
                b0s.append(None)
                continue
            sw = sum(w)
            b0 = sum(wi * r for wi, r in zip(w, resid)) / sw
            q += sum(wi * (r - b0) ** 2 for wi, r in zip(w, resid))
            b0s.append(b0)
        profile.append((th, b0s, q))
    best = min(profile, key=lambda t: t[2])
    inside = [t[0] for t in profile if t[2] <= best[2] + CHI2_1_95]
    ci = [min(inside), max(inside)]
    edge = step / 2.0
    n_bands = sum(len(p[2]) for p in prepared)
    return {
        "theta": best[0],
        "beta0": best[1],
        "oe_all": [math.exp(b) for b in best[1]],
        "q_min": best[2],
        "df": n_bands - len(prepared) - 1,
        "ci95_profile": ci,
        "theta_at_grid_edge": best[0] <= lo + edge or best[0] >= hi - edge,
        "ci_at_grid_edge": ci[0] <= lo + edge or ci[1] >= hi - edge,
        "n_sites": len(prepared),
        "n_bands": n_bands,
    }


def within_information(sites, theta=0.0):
    """Fisher information for theta from within-site contrasts (known
    variances, beta0_k profiled out): sum over sites of the weighted centred
    sum of squares of band_score."""
    info = 0.0
    for s in sites:
        _, _, _, w, bands = _prepare(s, need_y=False)
        x = [band_score(a, b, theta) for a, b in bands]
        sw = sum(w)
        mean = sum(wi * xi for wi, xi in zip(w, x)) / sw
        info += sum(wi * (xi - mean) ** 2 for wi, xi in zip(w, x))
    return info


def within_se(sites, theta=0.0):
    info = within_information(sites, theta)
    if info <= 0.0:
        raise ValueError("no within-site information about theta")
    return 1.0 / math.sqrt(info)


def consistency(theta_w, ci_w, theta_b, ci_b):
    """z for the difference between the within-site theta and the cross-site
    theta (selection.fit_sites). Standard errors are read off the profile 95%
    intervals as half-width / 1.96; the two estimates are treated as
    independent (within-site contrasts vs site-level totals)."""
    se_w = (ci_w[1] - ci_w[0]) / (2.0 * Z975)
    se_b = (ci_b[1] - ci_b[0]) / (2.0 * Z975)
    se = math.sqrt(se_w * se_w + se_b * se_b)
    if se <= 0.0:
        raise ValueError("need intervals with positive width")
    z = (theta_w - theta_b) / se
    return {"diff": theta_w - theta_b, "se": se, "z": z, "consistent": abs(z) <= Z975}


# ---------------------------------------------------------------- recency model
def poisson_band_factors(windows, rate0, kappa, g, half_width=8.0, steps=4000):
    """Band factors (hazard scale, relative to the whole cohort) when each
    person is tested at the times of a Poisson process with rate
    rate0 * exp(kappa * U) per day, U ~ N(0, 1), and the hazard multiplier is
    exp(g U - g^2 / 2). Band j holds people whose most recent test before
    index is in (W_(j-1), W_j] days. Trapezoid quadrature over U.
    Returns {"cum_rates": [...], "factors": [...]}."""
    ws = list(windows)
    if not ws or any(w <= 0 for w in ws) or any(b <= a for a, b in zip(ws, ws[1:])):
        raise ValueError("windows must be positive and increasing")
    if rate0 <= 0:
        raise ValueError("rate0 must be > 0")
    if steps < 2 or half_width <= 0:
        raise ValueError("need steps >= 2 and half_width > 0")
    k = len(ws)
    h = 2.0 * half_width / steps
    mass = 0.0
    tilt = 0.0
    band_mass = [0.0] * k
    band_tilt = [0.0] * k
    for i in range(steps + 1):
        u = -half_width + i * h
        wt = 0.5 if i == 0 or i == steps else 1.0
        dens = wt * math.exp(-0.5 * u * u)
        mult = math.exp(g * u - 0.5 * g * g)
        lam = rate0 * math.exp(kappa * u)
        prev = 1.0
        for j, w in enumerate(ws):
            cur = math.exp(-lam * w)
            pj = prev - cur
            band_mass[j] += dens * pj
            band_tilt[j] += dens * pj * mult
            prev = cur
        mass += dens
        tilt += dens * mult
    ref = tilt / mass
    cum = []
    c = 0.0
    for j in range(k):
        c += band_mass[j]
        cum.append(c / mass)
    factors = [(band_tilt[j] / band_mass[j]) / ref if band_mass[j] > 0 else None
               for j in range(k)]
    return {"cum_rates": cum, "factors": factors}


# ---------------------------------------------------------------- cohorts
def split_bands(rows_by_window):
    """rows_by_window[j]: omop.build_cohort rows built with uACR lookback
    window W_j (windows increasing). Band j holds the people first tested at
    window j, with their rows from that window's run."""
    if not rows_by_window:
        raise ValueError("need at least one window")
    ids = [r["person_id"] for r in rows_by_window[0]]
    for rows in rows_by_window[1:]:
        if [r["person_id"] for r in rows] != ids:
            raise ValueError("every window must give the same cohort in the same order")
    first = {}
    for j, rows in enumerate(rows_by_window):
        for r in rows:
            pid = r["person_id"]
            if r["tested"]:
                if pid not in first:
                    first[pid] = j
            elif pid in first:
                raise ValueError("windows are not nested: tested in a shorter window only")
    bands = [[r for r in rows if first.get(r["person_id"]) == j]
             for j, rows in enumerate(rows_by_window)]
    cum = []
    c = 0
    for b in bands:
        c += len(b)
        cum.append(c)
    return {"n_cohort": len(ids), "bands": bands, "cum_tested": cum}


def band_summary(rows, years):
    """Aalen-Johansen O/E (death competing) for one band at a horizon."""
    p_key = "kfre%d" % years
    rows = [r for r in rows if r[p_key] is not None]
    if not rows:
        return {"n": 0, "events": 0, "oe": None, "se_log_oe": None, "se_log_oe_dm": None,
                "observed_risk": None, "expected_risk": None}
    s = oe_time_to_event([r["time_days"] for r in rows], [r["event"] for r in rows],
                         [r[p_key] for r in rows], HORIZON_DAYS[years], cause=1, competing=True)
    return {"n": s["n"], "events": s["events"], "oe": s["oe"], "se_log_oe": s["se_log_oe"],
            "se_log_oe_dm": s["se_log_oe_dm"], "observed_risk": s["observed_risk"],
            "expected_risk": s["expected_risk"]}


def bands_from_rows(rows_by_window, windows, years, min_events=MIN_BAND_EVENTS):
    """Band summaries for one horizon. A band with fewer than min_events KRT
    events by the horizon is merged with the next band by dropping its upper
    window (the last band merges with the one before it, so the widest window,
    and the tested set, stays the same). Repeats until every band qualifies or
    one window is left. The site is usable when >= 2 bands remain, each with
    >= max(min_events, 1) events and a positive variance; then theta is fitted
    within this site alone."""
    windows = list(windows)
    if len(windows) != len(rows_by_window):
        raise ValueError("one row list per window")
    by_w = dict(zip(windows, rows_by_window))
    kept = list(windows)
    dropped = []
    while True:
        split = split_bands([by_w[w] for w in kept])
        summaries = [band_summary(b, years) for b in split["bands"]]
        small = [j for j, s in enumerate(summaries) if s["events"] < min_events]
        if not small or len(kept) == 1:
            break
        j = small[0]
        w = kept[j] if j < len(kept) - 1 else kept[j - 1]
        kept.remove(w)
        dropped.append(w)
    n = split["n_cohort"]
    cum_rates = [c / n for c in split["cum_tested"]] if n else []
    bands = []
    prev = 0
    for w, s in zip(kept, summaries):
        row = {"window_from_days": prev, "window_to_days": w}
        row.update(s)
        bands.append(row)
        prev = w
    usable = (len(kept) >= 2 and all(
        b["events"] >= max(min_events, 1) and b["se_log_oe_dm"] is not None
        and b["se_log_oe_dm"] > 0 for b in bands))
    fit = None
    if usable:
        fit = fit_within([{"cum_rates": cum_rates,
                           "log_oe": [math.log(b["oe"]) for b in bands],
                           "variances": [b["se_log_oe_dm"] ** 2 for b in bands]}])
    return {"years": years, "windows_requested": list(windows), "windows": kept,
            "windows_dropped": dropped, "n_cohort": n, "cum_tested": split["cum_tested"],
            "cum_rates": cum_rates, "bands": bands, "usable": usable, "fit": fit}


def site_bands(cdm, windows=PRIMARY_WINDOWS, years=(2, 5), region="north_america",
               min_events=MIN_BAND_EVENTS):
    windows = sorted(set(int(w) for w in windows))
    if len(windows) < 2 or windows[0] <= 0:
        raise ValueError("need at least 2 distinct positive windows")
    rows_by_window = [build_cohort(cdm, region=region, acr_lookback_days=w)["rows"]
                      for w in windows]
    return [bands_from_rows(rows_by_window, windows, y, min_events) for y in years]


def public_report(result, min_cell=MIN_CELL):
    """Copy of a bands_from_rows result safe to share: any band with a count
    from 1 to min_cell - 1 is blanked, and then every cumulative count and rate
    too (they would give band counts back by subtraction)."""
    def small(x):
        return x is not None and 0 < x < min_cell
    out = dict(result)
    bands = []
    any_masked = False
    for b in result["bands"]:
        b = dict(b)
        if small(b["n"]) or small(b["events"]):
            any_masked = True
            for key in ("n", "events", "oe", "se_log_oe", "se_log_oe_dm", "observed_risk",
                        "expected_risk"):
                b[key] = None
            b["masked"] = True
        else:
            b["masked"] = False
        bands.append(b)
    out["bands"] = bands
    if small(result["n_cohort"]):
        out["n_cohort"] = None
        any_masked = True
    if any_masked or any(small(c) for c in result["cum_tested"]):
        out["cum_tested"] = None
        out["cum_rates"] = None
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Within-site nested-window KFRE analysis from an "
                                 "OMOP CDM CSV extract (research only; not medical advice)")
    ap.add_argument("--cdm-dir", required=True)
    ap.add_argument("--out", required=True, help="JSON file to write (summary only)")
    ap.add_argument("--windows", type=int, nargs="+", default=list(PRIMARY_WINDOWS))
    ap.add_argument("--region", default="north_america",
                    choices=["north_america", "non_north_america"])
    ap.add_argument("--min-events", type=int, default=MIN_BAND_EVENTS)
    ap.add_argument("--min-cell", type=int, default=MIN_CELL)
    a = ap.parse_args(argv)
    cdm = read_cdm(a.cdm_dir)
    results = site_bands(cdm, a.windows, (2, 5), a.region, a.min_events)
    report = {"design": "nested-window v1 (protocol addendum 3)",
              "horizons": [public_report(r, a.min_cell) for r in results]}
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1, sort_keys=True)
    return report


if __name__ == "__main__":
    main()
