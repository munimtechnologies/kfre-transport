# SPDX-License-Identifier: Apache-2.0
"""Simulation study: does WHO GETS A URINE ALBUMIN TEST make the 4-variable
KFRE look miscalibrated at some health systems, and can a pre-declared
site-level meta-regression on the testing rate tell this apart from noise?

Everything here is synthetic. No patient data. The generator is deliberately
simple and is described in findings/sim-prereg-v1-2026-11-16.md.

Per patient: age ~ N(70, 10) clipped to [30, 95]; male ~ Bernoulli(0.5);
eGFR ~ Uniform(15, 45); a latent severity U ~ N(0, 1) that the record does not
hold; ln uACR = 4.0 + 1.0*U + N(0, 1.2^2).
True 5-year risk: KFRE linear predictor (with the TRUE uACR) + g*U, North
American 5-year baseline. Event ~ Bernoulli(true risk).
Tested for uACR with probability expit(a_site + d*U).

Per site: KFRE O/E among the tested (complete case), O/E among the untested
with uACR imputed at the site's geometric mean of tested values, and the
oracle O/E (everyone, true uACR).
Across sites: log O/E (tested) pooled by DerSimonian-Laird, and a random-
effects meta-regression of log O/E (tested) on the site testing rate
(Knapp-Hartung t test, k-2 df). Variance of log O/E taken as 1/O.

Research code only. Information only, not medical advice, not a medical device.
"""
import argparse
import json
import math
import os
import random
import statistics

from kfretransport import CODE_VERSION
from kfretransport.equations import kfre4_linear_predictor, risk_from_lp
from kfretransport.metrics import dersimonian_laird, logit, meta_regression

YEARS = 5
REGION = "north_america"
SITE_RATES = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80]

SCENARIOS = {
    "A_null": {"d": 0.0, "g": 0.0},
    "B_hidden_risk_random_testing": {"d": 0.0, "g": 0.5},
    "C_informative_testing": {"d": 1.0, "g": 0.5},
    "D_testing_follows_acr_only": {"d": 1.0, "g": 0.0},
}


def simulate_site(rng, n, test_logit, d, g):
    """Returns a list of (age, male, egfr, ln_acr, event, tested)."""
    rows = []
    for _ in range(n):
        age = min(95.0, max(30.0, rng.gauss(70.0, 10.0)))
        male = rng.random() < 0.5
        egfr = rng.uniform(15.0, 45.0)
        u = rng.gauss(0.0, 1.0)
        ln_acr = 4.0 + 1.0 * u + rng.gauss(0.0, 1.2)
        lp = kfre4_linear_predictor(age, male, egfr, math.exp(ln_acr)) + g * u
        event = rng.random() < risk_from_lp(lp, YEARS, REGION)
        z = test_logit + d * u
        p_test = 1.0 / (1.0 + math.exp(-z)) if z >= 0 else math.exp(z) / (1.0 + math.exp(z))
        tested = rng.random() < p_test
        rows.append((age, male, egfr, ln_acr, event, tested))
    return rows


def _pred(age, male, egfr, ln_acr):
    return risk_from_lp(kfre4_linear_predictor(age, male, egfr, math.exp(ln_acr)), YEARS, REGION)


def analyze_site(rows):
    tested = [r for r in rows if r[5]]
    untested = [r for r in rows if not r[5]]
    out = {"n": len(rows), "n_tested": len(tested), "rate": len(tested) / len(rows)}
    out["o_tested"] = sum(1 for r in tested if r[4])
    out["e_tested"] = sum(_pred(r[0], r[1], r[2], r[3]) for r in tested)
    out["o_all"] = sum(1 for r in rows if r[4])
    out["e_oracle"] = sum(_pred(r[0], r[1], r[2], r[3]) for r in rows)
    if tested and untested:
        imp = sum(r[3] for r in tested) / len(tested)
        out["o_untested"] = sum(1 for r in untested if r[4])
        out["e_untested"] = sum(_pred(r[0], r[1], r[2], imp) for r in untested)
    else:
        out["o_untested"] = out["e_untested"] = None
    return out


def run_replicate(rng, n_per_site, d, g, rates=SITE_RATES):
    sites = [analyze_site(simulate_site(rng, n_per_site, logit(rate), d, g)) for rate in rates]
    usable = [s for s in sites if s["o_tested"] > 0 and s["e_tested"] > 0]
    rep = {"sites": len(sites), "usable": len(usable)}
    o_t = sum(s["o_tested"] for s in sites)
    e_t = sum(s["e_tested"] for s in sites)
    o_u = sum(s["o_untested"] for s in sites if s["o_untested"] is not None)
    e_u = sum(s["e_untested"] for s in sites if s["e_untested"] is not None)
    rep["oe_tested"] = o_t / e_t if e_t > 0 else None
    rep["oe_untested"] = o_u / e_u if e_u > 0 else None
    rep["oe_oracle"] = sum(s["o_all"] for s in sites) / sum(s["e_oracle"] for s in sites)
    rep["rates"] = [s["rate"] for s in sites]
    rep["site_oe_tested"] = [s["o_tested"] / s["e_tested"] if s["e_tested"] > 0 else None
                             for s in sites]
    if len(usable) >= 3:
        y = [math.log(s["o_tested"] / s["e_tested"]) for s in usable]
        v = [1.0 / s["o_tested"] for s in usable]
        x = [s["rate"] for s in usable]
        mr = meta_regression(y, v, x)
        dl = dersimonian_laird(y, v)
        rep.update({"b1": mr["b1"], "reject": mr["reject"], "i2": dl["i2"],
                    "pooled_oe_dl": math.exp(dl["pooled"]),
                    "pi_oe": [math.exp(dl["pi"][0]), math.exp(dl["pi"][1])]})
    else:
        rep.update({"b1": None, "reject": None, "i2": None, "pooled_oe_dl": None, "pi_oe": None})
    return rep


def _median(values):
    vals = [x for x in values if x is not None]
    return statistics.median(vals) if vals else None


def summarize(reps):
    analysed = [r for r in reps if r["reject"] is not None]
    k = len(analysed)
    return {
        "replicates": len(reps),
        "analysed": k,
        "reject_rate": (sum(1 for r in analysed if r["reject"]) / k) if k else None,
        "median_b1": _median([r["b1"] for r in analysed]),
        "median_oe_tested": _median([r["oe_tested"] for r in reps]),
        "median_oe_untested": _median([r["oe_untested"] for r in reps]),
        "median_oe_oracle": _median([r["oe_oracle"] for r in reps]),
        "median_i2": _median([r["i2"] for r in analysed]),
        "median_pooled_oe_dl": _median([r["pooled_oe_dl"] for r in analysed]),
        "median_pi_low": _median([r["pi_oe"][0] for r in analysed]),
        "median_pi_high": _median([r["pi_oe"][1] for r in analysed]),
        "median_site_rates": [_median([r["rates"][i] for r in reps])
                              for i in range(len(reps[0]["rates"]))] if reps else [],
        "median_site_oe_tested": [_median([r["site_oe_tested"][i] for r in reps])
                                  for i in range(len(reps[0]["site_oe_tested"]))] if reps else [],
        "min_usable_sites": min(r["usable"] for r in reps) if reps else None,
    }


def grade(summary):
    """Grades for the seven predictions in findings/sim-prereg-v1-2026-11-16.md.
    The prose there is the authority; this is a convenience copy."""
    a = summary["A_null"]
    b = summary["B_hidden_risk_random_testing"]
    c = summary["C_informative_testing"]
    dd = summary["D_testing_follows_acr_only"]

    def ok(cond):
        return "PASS" if cond else "FAIL"

    def safe(x, f):
        return x is not None and f(x)

    return {
        "S1": ok(safe(a["reject_rate"], lambda r: r <= 0.08)),
        "S2": ok(safe(b["reject_rate"], lambda r: r <= 0.08)),
        "S3": ok(safe(dd["reject_rate"], lambda r: r <= 0.08)),
        "S4": ok(safe(c["reject_rate"], lambda r: r >= 0.50) and safe(c["median_b1"], lambda x: x < 0)),
        "S5": ok(safe(c["median_oe_tested"], lambda x: x > 1.0)
                 and safe(c["median_oe_untested"], lambda x: x < 1.0)),
        "S6": ok(safe(a["median_oe_tested"], lambda x: 0.9 <= x <= 1.1)
                 and safe(dd["median_oe_tested"], lambda x: 0.9 <= x <= 1.1)),
        "S7": ok(safe(b["median_oe_tested"], lambda x: x > 1.1)
                 and safe(b["median_oe_oracle"], lambda x: x > 1.1)),
    }


def run_study(k, seed, n_per_site, scenarios=None):
    scenarios = scenarios or SCENARIOS
    rng = random.Random(seed)
    results = {}
    for name, par in scenarios.items():
        reps = [run_replicate(rng, n_per_site, par["d"], par["g"]) for _ in range(k)]
        results[name] = summarize(reps)
        results[name]["params"] = par
    return results


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--k", type=int, default=200, help="replicates per scenario")
    ap.add_argument("--n", type=int, default=10000, help="patients per site")
    ap.add_argument("--seed", type=int, default=20261116)
    ap.add_argument("--out-dir", default="out")
    args = ap.parse_args(argv)
    summary = run_study(args.k, args.seed, args.n)
    out = {"code_version": CODE_VERSION, "k": args.k, "n_per_site": args.n,
           "seed": args.seed, "site_rates_target": SITE_RATES, "years": YEARS,
           "region": REGION, "scenarios": summary, "grades": grade(summary)}
    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, "sim.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(json.dumps(out["grades"], sort_keys=True))
    print("wrote", path)


if __name__ == "__main__":
    main()
