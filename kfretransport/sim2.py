# SPDX-License-Identifier: Apache-2.0
"""Simulation v2 (new in 0.2.6): person-level uACR testing histories.

Question: when the analysis code of this package is run end to end on
person-level time-to-event data (censoring, death competing), does it recover
the selection parameter theta of the selection curve? How does it behave when
testing follows a different mechanism, or none? Pre-registration:
findings/sim2-prereg-2026-11-19.md. Everything is synthetic; no patient data.

Per person (all independent draws):
- U ~ N(0, 1): the part of the log-hazard of kidney failure (KRT) that the
  equation's predictors miss. x ~ N(0, 1): the equation's own linear predictor.
- Predicted 5-year risk p = 1 - exp(-H), H = H0 * exp(SX * x).
- True KRT hazard per day: (H / 1826) * exp(g U - g^2 / 2), constant in time,
  so the prediction is right on average over U (hazard scale).
- Death hazard per day: -ln(1 - DEATH_5Y) / 1826 * exp(g_death U - g_death^2 / 2).
- Follow-up ends at C ~ Uniform(365, 3650) days. Observed: the first of KRT
  (event 1), death (event 2) or end of follow-up (event 0).
- Days since the latest uACR before index:
  threshold: V = rho U + beta_x x + sqrt(1 - rho^2 - beta_x^2) e, e ~ N(0, 1).
    Tested within window W_j  <=>  V > -Phi^-1(share * FRACTION_j). The person
    gets the midpoint of their band (45, 135 or 272.5 days); only the band
    matters. theta = g * rho.
  recency: the time since the latest test is Exponential with rate
    r0 * exp(kappa U - kappa^2 / 2), r0 = -ln(1 - share) / 365 (the Poisson
    recency model of nested.py). kappa = 0 is no selection at all.
- The prediction does not change with the window (no stale uACR values here).

Per site: nested.bands_from_rows (the code a health system would run: windows
90, 180, 365; Aalen-Johansen O/E per disjoint band; bands with < 20 events
merged), the Aalen-Johansen O/E of everyone tested within 365 days, and the
same bands with 1 - Kaplan-Meier (death censored) for comparison.
Per replicate (8 sites with different testing intensity): theta_W from
nested.fit_within over the usable sites, theta_B from selection.fit_sites on
the sites' 365-day tested groups, and nested.consistency between them.

Research code only. Information only, not medical advice, not a medical device.
Made with AI agents (Claude) working in MV2's name (MV2, Munim, Inc.).
"""

import argparse
import json
import math
import os
import random
from statistics import NormalDist

from kfretransport import CODE_VERSION
from kfretransport.metrics import oe_time_to_event
from kfretransport.nested import MIN_BAND_EVENTS, bands_from_rows, consistency, fit_within, split_bands
from kfretransport.selection import fit_sites

_N = NormalDist()
YEARS = 5
DAYS = 1826                      # omop.HORIZON_DAYS[5]
WINDOWS = (90, 180, 365)
THRESHOLD_FRACTIONS = (0.45, 0.70, 1.0)
SITE_SHARES = (0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70)
N_PER_SITE = 20000
H0 = 0.04
SX = 0.8
DEATH_5Y = 0.20
FOLLOW_UP = (365.0, 3650.0)

SCENARIOS = {
    "T69": {"mechanism": "threshold", "g": 0.8, "rho": 0.8625, "beta_x": 0.0, "g_death": 0.0,
            "theta": 0.69},
    "T30": {"mechanism": "threshold", "g": 0.8, "rho": 0.375, "beta_x": 0.0, "g_death": 0.0,
            "theta": 0.30},
    "N0": {"mechanism": "recency", "g": 0.8, "kappa": 0.0, "g_death": 0.0, "theta": 0.0},
    "R69": {"mechanism": "recency", "g": 0.8, "kappa": 0.8625, "g_death": 0.0, "theta": None},
    "TX": {"mechanism": "threshold", "g": 0.8, "rho": 0.8625, "beta_x": 0.3, "g_death": 0.0,
           "theta": 0.69},
    "TD": {"mechanism": "threshold", "g": 0.8, "rho": 0.8625, "beta_x": 0.0, "g_death": 0.5,
           "theta": 0.69},
}


def _midpoints(windows):
    out = []
    lo = 0.0
    for w in windows:
        out.append((lo + w) / 2.0)
        lo = float(w)
    return out


def simulate_site(rng, n, share, sc, windows=WINDOWS, fractions=THRESHOLD_FRACTIONS):
    """Returns a list of (days_since_latest_uacr, time_days, event, predicted_5y_risk).
    days_since_latest_uacr is math.inf for nobody-tested-within-the-last-window
    under the threshold mechanism."""
    if not (0.0 < share < 1.0):
        raise ValueError("share must be strictly between 0 and 1")
    g = sc["g"]
    gd = sc.get("g_death", 0.0)
    lam_death = -math.log(1.0 - DEATH_5Y) / DAYS
    mech = sc["mechanism"]
    if mech == "threshold":
        if len(fractions) != len(windows):
            raise ValueError("one fraction per window")
        rho = sc["rho"]
        bx = sc.get("beta_x", 0.0)
        resid = 1.0 - rho * rho - bx * bx
        if resid <= 0.0:
            raise ValueError("need rho^2 + beta_x^2 < 1")
        s = math.sqrt(resid)
        cuts = [-_N.inv_cdf(share * f) for f in fractions]
        mids = _midpoints(windows)
    elif mech == "recency":
        kappa = sc["kappa"]
        rate0 = -math.log(1.0 - share) / windows[-1]
    else:
        raise ValueError("unknown mechanism: %r" % (mech,))
    people = []
    for _ in range(n):
        u = rng.gauss(0.0, 1.0)
        x = rng.gauss(0.0, 1.0)
        e = rng.gauss(0.0, 1.0)
        h5 = H0 * math.exp(SX * x)
        p = -math.expm1(-h5)
        lam1 = h5 / DAYS * math.exp(g * u - 0.5 * g * g)
        lam2 = lam_death * math.exp(gd * u - 0.5 * gd * gd)
        t1 = rng.expovariate(lam1)
        t2 = rng.expovariate(lam2)
        c = rng.uniform(FOLLOW_UP[0], FOLLOW_UP[1])
        if mech == "threshold":
            v = rho * u + bx * x + s * e
            days = math.inf
            for cut, mid in zip(cuts, mids):
                if v > cut:
                    days = mid
                    break
        else:
            days = rng.expovariate(rate0 * math.exp(kappa * u - 0.5 * kappa * kappa))
        if t1 <= t2 and t1 <= c:
            t, ev = t1, 1
        elif t2 <= c:
            t, ev = t2, 2
        else:
            t, ev = c, 0
        people.append((days, t, ev, p))
    return people


def rows_by_window(people, windows=WINDOWS):
    """The rows omop.build_cohort would give with each uACR lookback window."""
    out = []
    for w in windows:
        rows = []
        for i, (days, t, ev, p) in enumerate(people):
            tested = days <= w
            rows.append({"person_id": i, "tested": tested, "time_days": t, "event": ev,
                         "kfre5": p if tested else None})
        out.append(rows)
    return out


def _tte(rows, competing):
    return oe_time_to_event([r["time_days"] for r in rows], [r["event"] for r in rows],
                            [r["kfre5"] for r in rows], DAYS, cause=1, competing=competing)


def analyse_site(people, windows=WINDOWS, min_events=MIN_BAND_EVENTS):
    rbw = rows_by_window(people, windows)
    res = bands_from_rows(rbw, windows, YEARS, min_events)
    n = len(people)
    out = {"n": n, "usable": res["usable"], "windows": res["windows"],
           "full": res["windows"] == list(windows), "cum_rates": res["cum_rates"],
           "band_oe": [b["oe"] for b in res["bands"]],
           "band_events": [b["events"] for b in res["bands"]],
           "band_var": [(b["se_log_oe_dm"] ** 2) if b["se_log_oe_dm"] else None
                        for b in res["bands"]],
           "site_theta": res["fit"]["theta"] if res["fit"] else None}
    widest = [r for r in rbw[-1] if r["tested"]]
    out["share_tested"] = len(widest) / n if n else 0.0
    out["oe_tested"] = None
    out["var_tested"] = None
    if widest:
        s = _tte(widest, True)
        if s["se_log_oe_dm"]:
            out["oe_tested"] = s["oe"]
            out["var_tested"] = s["se_log_oe_dm"] ** 2
    split = split_bands(rbw)
    km_oe, km_var, ok = [], [], True
    for b in split["bands"]:
        if not b:
            ok = False
            break
        s = _tte(b, False)
        if s["events"] < max(min_events, 1) or not s["se_log_oe_dm"]:
            ok = False
            break
        km_oe.append(s["oe"])
        km_var.append(s["se_log_oe_dm"] ** 2)
    out["km"] = None
    if ok and len(km_oe) >= 2:
        out["km"] = {"cum_rates": [c / n for c in split["cum_tested"]],
                     "log_oe": [math.log(o) for o in km_oe], "variances": km_var}
    return out


def run_replicate(rng, sc, n=N_PER_SITE, shares=SITE_SHARES, min_events=MIN_BAND_EVENTS):
    sites = [analyse_site(simulate_site(rng, n, s, sc), WINDOWS, min_events) for s in shares]
    rep = {"usable_sites": sum(1 for s in sites if s["usable"]),
           "share_tested": [s["share_tested"] for s in sites],
           "site_cum_rates": [s["cum_rates"] if s["full"] else None for s in sites],
           "site_band_oe": [s["band_oe"] if s["full"] else None for s in sites],
           "theta_w": None, "ci_w": None, "q_min": None, "df": None, "edge_w": None,
           "theta_b": None, "ci_b": None, "z": None, "consistent": None, "theta_w_km": None}
    within = [{"cum_rates": s["cum_rates"], "log_oe": [math.log(o) for o in s["band_oe"]],
               "variances": s["band_var"]} for s in sites if s["usable"]]
    if within:
        fw = fit_within(within)
        rep.update({"theta_w": fw["theta"], "ci_w": fw["ci95_profile"], "q_min": fw["q_min"],
                    "df": fw["df"], "edge_w": fw["theta_at_grid_edge"] or fw["ci_at_grid_edge"]})
    cross = [s for s in sites if s["oe_tested"] and s["var_tested"]
             and 0.0 < s["share_tested"] < 1.0]
    if len(cross) >= 3 and len(set(s["share_tested"] for s in cross)) >= 2:
        fb = fit_sites([s["share_tested"] for s in cross],
                       [math.log(s["oe_tested"]) for s in cross],
                       [s["var_tested"] for s in cross])
        rep.update({"theta_b": fb["theta"], "ci_b": fb["ci95_profile"]})
    km = [s["km"] for s in sites if s["km"]]
    if km:
        rep["theta_w_km"] = fit_within(km)["theta"]
    if rep["theta_w"] is not None and rep["theta_b"] is not None:
        try:
            c = consistency(rep["theta_w"], rep["ci_w"], rep["theta_b"], rep["ci_b"])
            rep["z"] = c["z"]
            rep["consistent"] = c["consistent"]
        except ValueError:
            pass
    return rep


# ---------------------------------------------------------------- summaries
def chi2_95(df):
    """Upper 5% point of chi-squared, Wilson-Hilferty approximation (df >= 1)."""
    if df < 1:
        raise ValueError("df must be >= 1")
    a = 2.0 / (9.0 * df)
    return df * (1.0 - a + _N.inv_cdf(0.95) * math.sqrt(a)) ** 3


def _vals(xs):
    return [x for x in xs if x is not None]


def _mean(xs):
    v = _vals(xs)
    return sum(v) / len(v) if v else None


def _sd(xs):
    v = _vals(xs)
    if len(v) < 2:
        return None
    m = sum(v) / len(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def _share(flags):
    v = _vals(flags)
    return sum(1 for f in v if f) / len(v) if v else None


def _covers(ci, t):
    if ci is None or t is None:
        return None
    return ci[0] <= t <= ci[1]


def _excludes_zero(ci):
    if ci is None:
        return None
    return ci[0] > 0.0 or ci[1] < 0.0


def _misfit(rep):
    if rep.get("q_min") is None or rep.get("df") is None or rep["df"] < 1:
        return None
    return rep["q_min"] > chi2_95(rep["df"])


def _w_below_b(rep):
    if rep.get("theta_w") is None or rep.get("theta_b") is None:
        return None
    return rep["theta_w"] < rep["theta_b"]


def _site_means(reps, key):
    n_sites = len(reps[0][key]) if reps else 0
    out = []
    for i in range(n_sites):
        rows = [r[key][i] for r in reps if r[key][i] is not None]
        if not rows:
            out.append(None)
            continue
        out.append([sum(row[j] for row in rows) / len(rows) for j in range(len(rows[0]))])
    return out


def summarize(reps, theta):
    return {
        "replicates": len(reps),
        "theta_true": theta,
        "analysed_within": len(_vals([r["theta_w"] for r in reps])),
        "analysed_cross": len(_vals([r["theta_b"] for r in reps])),
        "mean_theta_w": _mean([r["theta_w"] for r in reps]),
        "sd_theta_w": _sd([r["theta_w"] for r in reps]),
        "coverage_w": _share([_covers(r["ci_w"], theta) for r in reps]),
        "reject0_w": _share([_excludes_zero(r["ci_w"]) for r in reps]),
        "mean_theta_b": _mean([r["theta_b"] for r in reps]),
        "sd_theta_b": _sd([r["theta_b"] for r in reps]),
        "coverage_b": _share([_covers(r["ci_b"], theta) for r in reps]),
        "reject0_b": _share([_excludes_zero(r["ci_b"]) for r in reps]),
        "share_w_below_b": _share([_w_below_b(r) for r in reps]),
        "consistency_reject": _share([(not r["consistent"]) if r["consistent"] is not None
                                      else None for r in reps]),
        "misfit_share": _share([_misfit(r) for r in reps]),
        "mean_theta_w_km": _mean([r["theta_w_km"] for r in reps]),
        "edge_share": _share([r["edge_w"] for r in reps]),
        "mean_usable_sites": _mean([r["usable_sites"] for r in reps]),
        "min_usable_sites": min(r["usable_sites"] for r in reps) if reps else None,
        "mean_share_tested": [_mean([r["share_tested"][i] for r in reps])
                              for i in range(len(reps[0]["share_tested"]))] if reps else [],
        "mean_site_cum_rates": _site_means(reps, "site_cum_rates"),
        "mean_site_band_oe": _site_means(reps, "site_band_oe"),
    }


def grade(s):
    """Grades for P1-P7 of findings/sim2-prereg-2026-11-19.md. The prose there
    is the authority; this is a convenience copy."""
    def ok(c):
        return "PASS" if c else "FAIL"

    def near(x, t, tol):
        return x is not None and abs(x - t) <= tol

    def within(x, lo, hi):
        return x is not None and lo <= x <= hi

    t69, t30, n0, r69, td = s["T69"], s["T30"], s["N0"], s["R69"], s["TD"]
    return {
        "P1": ok(near(t69["mean_theta_w"], 0.69, 0.10) and within(t69["coverage_w"], 0.85, 0.99)),
        "P2": ok(near(t30["mean_theta_w"], 0.30, 0.10) and within(t30["coverage_w"], 0.85, 0.99)),
        "P3": ok(within(n0["reject0_w"], 0.0, 0.08) and near(n0["mean_theta_w"], 0.0, 0.05)),
        "P4": ok(r69["mean_theta_w"] is not None and r69["mean_theta_b"] is not None
                 and r69["mean_theta_w"] < r69["mean_theta_b"]
                 and within(r69["share_w_below_b"], 0.80, 1.0)),
        "P5": ok(within(t69["consistency_reject"], 0.0, 0.10)
                 and within(t30["consistency_reject"], 0.0, 0.10)),
        "P6": ok(near(t69["mean_theta_b"], 0.69, 0.10) and within(t69["coverage_b"], 0.85, 0.99)
                 and near(t30["mean_theta_b"], 0.30, 0.10)
                 and within(t30["coverage_b"], 0.85, 0.99)),
        "P7": ok(t69["mean_theta_w"] is not None and td["mean_theta_w"] is not None
                 and t69["mean_theta_w"] < 0.69 and td["mean_theta_w"] < t69["mean_theta_w"]),
    }


def _compact(rep):
    keys = ("usable_sites", "theta_w", "ci_w", "q_min", "df", "theta_b", "ci_b", "z",
            "theta_w_km")
    return {k: rep[k] for k in keys}


def run_study(reps, seed, n=N_PER_SITE, names=None, shares=SITE_SHARES):
    names = list(names) if names else list(SCENARIOS)
    out = {}
    for name in names:
        sc = SCENARIOS[name]
        rng = random.Random("%d:%s" % (seed, name))
        rs = [run_replicate(rng, sc, n, shares) for _ in range(reps)]
        summary = summarize(rs, sc["theta"])
        summary["params"] = sc
        summary["replicate_rows"] = [_compact(r) for r in rs]
        out[name] = summary
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Simulation v2: person-level testing histories "
                                 "(research only; synthetic data)")
    ap.add_argument("--reps", type=int, default=200, help="replicates per scenario")
    ap.add_argument("--n", type=int, default=N_PER_SITE, help="people per site")
    ap.add_argument("--seed", type=int, default=20261119)
    ap.add_argument("--scenarios", nargs="+", default=None, choices=sorted(SCENARIOS))
    ap.add_argument("--out-dir", default="out")
    a = ap.parse_args(argv)
    res = run_study(a.reps, a.seed, a.n, a.scenarios)
    out = {"code_version": CODE_VERSION, "design": "sim2 (findings/sim2-prereg-2026-11-19.md)",
           "reps": a.reps, "n_per_site": a.n, "seed": a.seed, "site_shares": list(SITE_SHARES),
           "windows": list(WINDOWS), "threshold_fractions": list(THRESHOLD_FRACTIONS),
           "h0": H0, "sx": SX, "death_5y": DEATH_5Y, "follow_up_days": list(FOLLOW_UP),
           "scenarios": res}
    if set(SCENARIOS) <= set(res):
        out["grades"] = grade(res)
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, "sim2.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(json.dumps(out.get("grades", {}), sort_keys=True))
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
