# SPDX-License-Identifier: Apache-2.0
"""Pool per-site summaries written by kfretransport.omop into the multi-site
analysis pre-declared in findings/protocol-addendum-1-2026-11-17.md.

Primary analysis (per horizon, 2 and 5 years):
- per site: log O/E, where O is the Aalen-Johansen cumulative incidence of
  kidney replacement therapy by the horizon (death competing) among everyone
  with a KFRE prediction and E is their mean predicted risk ("tte_aj"), with
  the delta-method variance (metrics.cumulative_incidence_variance);
- DerSimonian-Laird random-effects pooling of log O/E, with a 95% prediction
  interval for a new site;
- random-effects meta-regression of log O/E on the site's measured-uACR
  testing rate (Knapp-Hartung, t with k-2 df). The question: does a site's
  testing practice predict how far the KFRE is off there?

Secondary analyses: the same with the Poisson variance, with 1 - Kaplan-Meier
(death censored), and complete-case O/E (Poisson variance only).

Small cells: with min_cell (N3C does not allow results about fewer than 20
people to leave its enclave), counts from 1 to min_cell - 1 are masked in the
published per-site rows, and so is anything that would give them back (the
testing rate, a count recoverable by subtraction, the Poisson variance 1/events).
Pooled results are computed from the unmasked values inside the enclave.

Research code only. Information only, not medical advice, not a medical device.
"""
import argparse
import json
import math
import os

from kfretransport.metrics import dersimonian_laird, meta_regression

MEASURES = ("tte_aj", "tte_km", "oe")
VARIANCES = ("delta", "poisson")
PRIMARY = ("tte_aj", "delta")
SECONDARY = (("tte_aj", "poisson"), ("tte_km", "delta"), ("oe", "poisson"))


def _summary_for(result, years):
    for s in result.get("summary", []):
        if s.get("years") == years:
            return s
    raise ValueError("no %d-year summary in this site result" % years)


def site_estimate(result, years, measure="tte_aj", variance="delta", site=None):
    """One site's log O/E and its variance from an omop.py result dict.
    status: ok, no_prediction, no_events, no_variance."""
    if measure not in MEASURES or variance not in VARIANCES:
        raise ValueError("unknown measure or variance")
    if measure == "oe" and variance == "delta":
        raise ValueError("complete-case O/E has only the Poisson variance")
    s = _summary_for(result, years)
    row = {"site": site, "years": years, "measure": measure, "variance": variance,
           "n_cohort": s.get("n_cohort"), "n_tested": s.get("n_tested"),
           "testing_rate": s.get("testing_rate"), "events": None, "oe": None,
           "log_oe": None, "var_log_oe": None, "status": "ok"}
    if measure == "oe":
        oe, events = s.get("oe"), s.get("events_scored")
        se = (1.0 / math.sqrt(events)) if events else None
        if oe is None:
            row["status"] = "no_prediction"
            return row
    else:
        m = s.get(measure)
        if m is None:
            row["status"] = "no_prediction"
            return row
        oe, events = m.get("oe"), m.get("events")
        se = m.get("se_log_oe_dm") if variance == "delta" else m.get("se_log_oe")
    row["events"], row["oe"] = events, oe
    if not events or not oe or oe <= 0:
        row["status"] = "no_events"
        return row
    if se is None or se <= 0:
        row["status"] = "no_variance"
        return row
    row["log_oe"] = math.log(oe)
    row["var_log_oe"] = se * se
    return row


def pool(rows):
    """Random-effects pooling and meta-regression on the testing rate over the
    rows whose status is ok."""
    ok = [r for r in rows if r["status"] == "ok"]
    excluded = {}
    for r in rows:
        if r["status"] != "ok":
            excluded[r["status"]] = excluded.get(r["status"], 0) + 1
    out = {"k_sites": len(rows), "k_used": len(ok), "excluded": excluded,
           "random_effects": None, "meta_regression": None, "notes": []}
    if len(ok) < 2:
        out["notes"].append("fewer than 2 usable sites: no pooling")
        return out
    y = [r["log_oe"] for r in ok]
    v = [r["var_log_oe"] for r in ok]
    dl = dersimonian_laird(y, v)
    out["random_effects"] = {
        "pooled_log_oe": dl["pooled"], "se": dl["se"], "tau2": dl["tau2"],
        "q": dl["q"], "i2": dl["i2"], "pooled_oe": math.exp(dl["pooled"]),
        "ci_oe": [math.exp(dl["ci"][0]), math.exp(dl["ci"][1])],
        "pi_oe": [math.exp(dl["pi"][0]), math.exp(dl["pi"][1])] if dl["pi"] else None}
    x = [r["testing_rate"] for r in ok]
    if len(ok) < 3:
        out["notes"].append("fewer than 3 usable sites: no meta-regression")
    elif any(xi is None for xi in x):
        out["notes"].append("a usable site has no testing rate: no meta-regression")
    else:
        try:
            mr = meta_regression(y, v, x)
        except ValueError as exc:
            out["notes"].append("no meta-regression: %s" % exc)
        else:
            out["meta_regression"] = {
                "moderator": "testing_rate (share of the cohort with a measured uACR, 0-1)",
                "b0": mr["b0"], "b1": mr["b1"], "se_b1": mr["se_b1"], "t": mr["t"],
                "df": mr["df"], "ci_b1": list(mr["ci_b1"]), "reject_5pct": mr["reject"],
                "tau2": mr["tau2"], "q_res": mr["q_res"],
                "oe_ratio_per_10_points": math.exp(0.1 * mr["b1"])}
    return out


def mask_row(row, min_cell):
    """Published copy of a site row with small cells masked (see module doc).
    Zero counts are not masked."""
    out = dict(row)
    if not min_cell:
        return out
    tag = "<%d" % min_cell

    def small(v):
        return isinstance(v, int) and 0 < v < min_cell

    n, t = row.get("n_cohort"), row.get("n_tested")
    untested = (n - t) if isinstance(n, int) and isinstance(t, int) else None
    for key in ("n_cohort", "n_tested", "events"):
        if small(row.get(key)):
            out[key] = tag
    if small(untested) and out.get("n_tested") != tag:
        out["n_tested"] = tag
    if out.get("n_cohort") == tag or out.get("n_tested") == tag:
        out["testing_rate"] = tag
    if out.get("events") == tag and row.get("variance") == "poisson":
        out["var_log_oe"] = tag
    return out


def run_plan(named_results, years, min_cell=None):
    """Primary and secondary analyses for one horizon. named_results: list of
    (site name, omop.py result dict)."""
    plan = {"years": years, "min_cell": min_cell, "analyses": []}
    for role, combos in (("primary", (PRIMARY,)), ("secondary", SECONDARY)):
        for measure, variance in combos:
            rows = [site_estimate(res, years, measure, variance, site=name)
                    for name, res in named_results]
            plan["analyses"].append({"role": role, "measure": measure, "variance": variance,
                                     "sites": [mask_row(r, min_cell) for r in rows],
                                     "result": pool(rows)})
    return plan


def main(argv=None):
    ap = argparse.ArgumentParser(description="Pool kfre-omop site summaries "
                                 "(research only; not medical advice)")
    ap.add_argument("summaries", nargs="+", help="site summary JSON files from kfre-omop")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-cell", type=int, default=20,
                    help="mask published counts from 1 to this minus 1 (0 = off)")
    a = ap.parse_args(argv)
    named = []
    for path in a.summaries:
        with open(path, encoding="utf-8") as f:
            named.append((os.path.splitext(os.path.basename(path))[0], json.load(f)))
    result = {"plan": [run_plan(named, y, a.min_cell or None) for y in (2, 5)],
              "pcr_conversion": sorted({bool(r.get("pcr_conversion")) for _, r in named})}
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1, sort_keys=True)
    return result


if __name__ == "__main__":
    main()
