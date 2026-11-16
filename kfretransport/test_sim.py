# SPDX-License-Identifier: Apache-2.0
import math
import random
import unittest

from kfretransport.equations import kfre4
from kfretransport.sim import (SCENARIOS, analyze_site, grade, run_replicate,
                               run_study, simulate_site, summarize)


class TestSimulateSite(unittest.TestCase):
    def test_shape_and_ranges(self):
        rows = simulate_site(random.Random(1), 200, 0.0, 0.0, 0.0)
        self.assertEqual(len(rows), 200)
        for age, male, egfr, ln_acr, event, tested in rows:
            self.assertTrue(30.0 <= age <= 95.0)
            self.assertTrue(15.0 <= egfr <= 45.0)
            self.assertIn(male, (True, False))
            self.assertIn(event, (True, False))
            self.assertIn(tested, (True, False))
            self.assertTrue(math.isfinite(ln_acr))

    def test_deterministic(self):
        a = simulate_site(random.Random(7), 50, -1.0, 1.0, 0.5)
        b = simulate_site(random.Random(7), 50, -1.0, 1.0, 0.5)
        self.assertEqual(a, b)

    def test_testing_logit_extremes(self):
        hi = simulate_site(random.Random(3), 300, 20.0, 0.0, 0.0)
        lo = simulate_site(random.Random(3), 300, -20.0, 0.0, 0.0)
        self.assertEqual(sum(1 for r in hi if r[5]), 300)
        self.assertEqual(sum(1 for r in lo if r[5]), 0)


class TestAnalyzeSite(unittest.TestCase):
    def test_counts_and_imputation(self):
        rows = [
            (60.0, True, 30.0, 2.0, True, True),
            (70.0, False, 20.0, 4.0, False, True),
            (80.0, True, 40.0, 5.0, True, False),
            (65.0, False, 25.0, 1.0, False, False),
        ]
        s = analyze_site(rows)
        self.assertEqual((s["n"], s["n_tested"], s["rate"]), (4, 2, 0.5))
        self.assertEqual((s["o_tested"], s["o_untested"], s["o_all"]), (1, 1, 2))
        e_t = kfre4(60.0, True, 30.0, math.exp(2.0)) + kfre4(70.0, False, 20.0, math.exp(4.0))
        self.assertAlmostEqual(s["e_tested"], e_t, places=12)
        # untested are imputed at the mean tested ln uACR = (2 + 4) / 2 = 3
        e_u = kfre4(80.0, True, 40.0, math.exp(3.0)) + kfre4(65.0, False, 25.0, math.exp(3.0))
        self.assertAlmostEqual(s["e_untested"], e_u, places=12)
        e_o = e_t + kfre4(80.0, True, 40.0, math.exp(5.0)) + kfre4(65.0, False, 25.0, math.exp(1.0))
        self.assertAlmostEqual(s["e_oracle"], e_o, places=12)

    def test_nobody_untested(self):
        s = analyze_site([(60.0, True, 30.0, 2.0, True, True)])
        self.assertIsNone(s["o_untested"])
        self.assertIsNone(s["e_untested"])


class TestStudy(unittest.TestCase):
    def test_replicate_keys(self):
        rep = run_replicate(random.Random(5), 400, 1.0, 0.5)
        self.assertEqual(rep["sites"], 8)
        self.assertEqual(len(rep["rates"]), 8)
        for key in ("oe_tested", "oe_untested", "oe_oracle", "b1", "reject", "i2"):
            self.assertIn(key, rep)

    def test_small_study_runs_and_grades(self):
        res = run_study(2, 11, 300)
        self.assertEqual(set(res), set(SCENARIOS))
        for name, s in res.items():
            self.assertEqual(s["replicates"], 2)
            self.assertEqual(len(s["median_site_rates"]), 8)
        g = grade(res)
        self.assertEqual(sorted(g), ["S1", "S2", "S3", "S4", "S5", "S6", "S7"])
        for v in g.values():
            self.assertIn(v, ("PASS", "FAIL"))

    def test_summarize_skips_unanalysed(self):
        reps = [
            {"reject": True, "b1": -1.0, "oe_tested": 1.2, "oe_untested": 0.9, "oe_oracle": 1.0,
             "i2": 0.5, "pooled_oe_dl": 1.1, "pi_oe": [0.8, 1.5], "rates": [0.1, 0.2],
             "site_oe_tested": [1.3, 1.1], "usable": 2},
            {"reject": None, "b1": None, "oe_tested": 1.0, "oe_untested": None, "oe_oracle": 1.0,
             "i2": None, "pooled_oe_dl": None, "pi_oe": None, "rates": [0.3, 0.4],
             "site_oe_tested": [None, 1.0], "usable": 1},
        ]
        s = summarize(reps)
        self.assertEqual((s["replicates"], s["analysed"]), (2, 1))
        self.assertEqual(s["reject_rate"], 1.0)
        self.assertAlmostEqual(s["median_oe_tested"], 1.1)
        self.assertAlmostEqual(s["median_oe_untested"], 0.9)
        self.assertAlmostEqual(s["median_site_oe_tested"][0], 1.3)
        self.assertAlmostEqual(s["median_site_oe_tested"][1], 1.05)
        self.assertEqual(s["min_usable_sites"], 1)


if __name__ == "__main__":
    unittest.main()
