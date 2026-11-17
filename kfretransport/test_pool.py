# SPDX-License-Identifier: Apache-2.0
# Every expected value is worked out by hand in the comments.
import json
import math
import os
import tempfile
import unittest

from kfretransport.pool import main, mask_row, pool, run_plan, site_estimate


def result(oe, se_dm, se_p, events, n_cohort, n_tested, years=2, cc_oe=None, cc_events=0):
    """A minimal omop.py-style site result with one horizon."""
    tte = {"n": n_tested, "events": events, "oe": oe, "observed_risk": oe * 0.1,
           "expected_risk": 0.1, "se_log_oe": se_p}
    if se_dm != "absent":
        tte["se_log_oe_dm"] = se_dm
    return {"pcr_conversion": False,
            "summary": [{"years": years, "n_cohort": n_cohort, "n_tested": n_tested,
                         "testing_rate": n_tested / n_cohort, "oe": cc_oe,
                         "events_scored": cc_events, "tte_aj": tte, "tte_km": tte}]}


class TestSiteEstimate(unittest.TestCase):
    def test_delta_and_poisson(self):
        r = result(2.0, 0.5, 0.25, 16, 1000, 400)
        row = site_estimate(r, 2, site="A")
        self.assertEqual(row["status"], "ok")
        self.assertAlmostEqual(row["log_oe"], math.log(2.0), places=12)
        self.assertAlmostEqual(row["var_log_oe"], 0.25, places=12)   # 0.5^2
        self.assertAlmostEqual(row["testing_rate"], 0.4, places=12)
        row = site_estimate(r, 2, variance="poisson")
        self.assertAlmostEqual(row["var_log_oe"], 0.0625, places=12)  # 0.25^2

    def test_complete_case_poisson_only(self):
        r = result(2.0, 0.5, 0.25, 16, 1000, 400, cc_oe=1.5, cc_events=4)
        row = site_estimate(r, 2, measure="oe", variance="poisson")
        self.assertAlmostEqual(row["var_log_oe"], 0.25, places=12)    # 1/4
        self.assertAlmostEqual(row["log_oe"], math.log(1.5), places=12)
        with self.assertRaises(ValueError):
            site_estimate(r, 2, measure="oe", variance="delta")

    def test_statuses(self):
        self.assertEqual(site_estimate(result(0.0, None, None, 0, 100, 50), 2)["status"],
                         "no_events")
        # a 0.2.1 summary has no se_log_oe_dm: usable only with the Poisson variance
        old = result(1.2, "absent", 0.5, 4, 100, 50)
        self.assertEqual(site_estimate(old, 2)["status"], "no_variance")
        self.assertEqual(site_estimate(old, 2, variance="poisson")["status"], "ok")
        none = result(1.0, 0.5, 0.5, 4, 100, 50)
        none["summary"][0]["tte_aj"] = None
        self.assertEqual(site_estimate(none, 2)["status"], "no_prediction")
        with self.assertRaises(ValueError):
            site_estimate(none, 5)


class TestPool(unittest.TestCase):
    def test_exact_line_in_testing_rate(self):
        # log O/E 0, 1, 2 with variance 1 at testing rates 0.2, 0.5, 0.8.
        # DL: fixed mean 1, Q = 2 = k-1 -> tau2 0; pooled O/E = e; se = sqrt(1/3)
        # meta-regression: exact line, b1 = 1/0.3 = 3.3333333, b0 = -2/3;
        # O/E ratio per 10 points of testing = exp(0.33333) = 1.3956124
        rows = [site_estimate(result(math.exp(y), 1.0, 1.0, 50, 1000, t), 2)
                for y, t in ((0.0, 200), (1.0, 500), (2.0, 800))]
        out = pool(rows)
        self.assertEqual(out["k_used"], 3)
        re = out["random_effects"]
        self.assertAlmostEqual(re["pooled_oe"], math.e, places=9)
        self.assertAlmostEqual(re["tau2"], 0.0, places=12)
        self.assertAlmostEqual(re["se"], 0.5773503, places=6)
        self.assertAlmostEqual(re["ci_oe"][1], math.exp(1.0 + 1.959964 * 0.5773503), places=4)
        mr = out["meta_regression"]
        self.assertAlmostEqual(mr["b1"], 10.0 / 3.0, places=9)
        self.assertAlmostEqual(mr["b0"], -2.0 / 3.0, places=9)
        self.assertAlmostEqual(mr["oe_ratio_per_10_points"], 1.3956124, places=6)
        self.assertTrue(mr["reject_5pct"])

    def test_meta_regression_not_significant(self):
        # rates 0.25, 0.25, 0.75, 0.75 and log O/E 0, 1, 1, 2 (variance 1):
        # in units of (rate - 0.25)/0.5 this is the 0,0,1,1 / 0,1,1,2 case of
        # test_metrics (b1 1, se 0.7071068, t 1.4142136), so here b1 = 2,
        # se_b1 = 1.4142136, t = 1.4142136 < t(2) = 4.302653
        rows = [site_estimate(result(math.exp(y), 1.0, 1.0, 50, 100, t), 2)
                for y, t in ((0.0, 25), (1.0, 25), (1.0, 75), (2.0, 75))]
        mr = pool(rows)["meta_regression"]
        self.assertAlmostEqual(mr["b1"], 2.0, places=9)
        self.assertAlmostEqual(mr["se_b1"], 1.4142136, places=6)
        self.assertAlmostEqual(mr["t"], 1.4142136, places=6)
        self.assertFalse(mr["reject_5pct"])

    def test_too_few_sites_and_exclusions(self):
        rows = [site_estimate(result(1.0, 0.5, 0.5, 4, 100, 50), 2),
                site_estimate(result(0.0, None, None, 0, 100, 50), 2)]
        out = pool(rows)
        self.assertEqual((out["k_sites"], out["k_used"]), (2, 1))
        self.assertEqual(out["excluded"], {"no_events": 1})
        self.assertIsNone(out["random_effects"])
        # two usable sites: pooling but no meta-regression
        rows[1] = site_estimate(result(2.0, 0.5, 0.5, 4, 100, 60), 2)
        out = pool(rows)
        self.assertIsNotNone(out["random_effects"])
        self.assertIsNone(out["meta_regression"])
        # same testing rate everywhere: no spread, no meta-regression
        rows = [site_estimate(result(math.exp(y), 1.0, 1.0, 50, 100, 50), 2)
                for y in (0.0, 1.0, 2.0)]
        out = pool(rows)
        self.assertIsNone(out["meta_regression"])
        self.assertIn("no meta-regression: moderator has no spread", out["notes"])


class TestMasking(unittest.TestCase):
    def test_small_counts_and_what_gives_them_back(self):
        row = site_estimate(result(1.5, 0.4, 0.5, 3, 500, 12), 2, variance="poisson")
        m = mask_row(row, 20)
        self.assertEqual((m["n_cohort"], m["n_tested"], m["events"]), (500, "<20", "<20"))
        self.assertEqual(m["testing_rate"], "<20")
        self.assertEqual(m["var_log_oe"], "<20")      # Poisson variance is 1/events
        self.assertEqual(row["n_tested"], 12)          # the input is not changed
        row = site_estimate(result(1.5, 0.4, 0.5, 3, 500, 12), 2, variance="delta")
        self.assertAlmostEqual(mask_row(row, 20)["var_log_oe"], 0.16, places=12)

    def test_untested_complement(self):
        # 100 in the cohort, 95 tested: 5 untested would follow by subtraction
        row = site_estimate(result(1.0, 0.3, 0.3, 40, 100, 95), 2)
        m = mask_row(row, 20)
        self.assertEqual((m["n_cohort"], m["n_tested"], m["testing_rate"]), (100, "<20", "<20"))
        self.assertEqual(m["events"], 40)

    def test_off(self):
        row = site_estimate(result(1.0, 0.3, 0.3, 4, 10, 5), 2)
        self.assertEqual(mask_row(row, None), row)


class TestPlanAndCli(unittest.TestCase):
    def test_plan_and_cli(self):
        sites = [("s%d" % i, result(math.exp(y), 1.0, 1.0, 50, 1000, t))
                 for i, (y, t) in enumerate(((0.0, 200), (1.0, 500), (2.0, 800)))]
        plan = run_plan(sites, 2, min_cell=20)
        roles = [(a["role"], a["measure"], a["variance"]) for a in plan["analyses"]]
        self.assertEqual(roles[0], ("primary", "tte_aj", "delta"))
        self.assertEqual(len(roles), 4)
        self.assertAlmostEqual(plan["analyses"][0]["result"]["random_effects"]["pooled_oe"],
                               math.e, places=9)
        # complete case: no cc O/E in these results -> nobody usable
        self.assertEqual(plan["analyses"][3]["result"]["k_used"], 0)
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for name, res in sites:
                # the CLI needs both horizons; reuse the 2-year block as 5-year
                five = dict(res["summary"][0], years=5)
                res = dict(res, summary=[res["summary"][0], five])
                p = os.path.join(d, name + ".json")
                with open(p, "w", encoding="utf-8") as f:
                    json.dump(res, f)
                paths.append(p)
            out = os.path.join(d, "pooled.json")
            main(paths + ["--out", out])
            with open(out, encoding="utf-8") as f:
                data = json.load(f)
        self.assertEqual([p["years"] for p in data["plan"]], [2, 5])
        self.assertEqual(data["plan"][0]["min_cell"], 20)
        self.assertEqual(data["plan"][0]["analyses"][0]["sites"][0]["site"], "s0")
        self.assertEqual(data["pcr_conversion"], [False])


if __name__ == "__main__":
    unittest.main()
