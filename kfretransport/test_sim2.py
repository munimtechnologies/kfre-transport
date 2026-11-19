# SPDX-License-Identifier: Apache-2.0
# Tests for sim2.py (0.2.6): generator properties, a hand-traced site, determinism,
# and the summary / grading arithmetic on invented replicate rows (working in comments).
import math
import random
import unittest

from kfretransport.sim2 import (SCENARIOS, analyse_site, chi2_95, grade, rows_by_window,
                                run_replicate, simulate_site, summarize)


def _share_within(people, w):
    return sum(1 for d, _, _, _ in people if d <= w) / len(people)


class TestGenerator(unittest.TestCase):
    def test_threshold_shares(self):
        # share 0.4 and fractions 0.45 / 0.70 / 1.0 -> 0.18, 0.28, 0.40 tested within
        # 90, 180, 365 days (binomial sd at n = 4000 is under 0.008)
        people = simulate_site(random.Random(1), 4000, 0.4, SCENARIOS["T69"])
        for w, target in ((90, 0.18), (180, 0.28), (365, 0.40)):
            self.assertAlmostEqual(_share_within(people, w), target, delta=0.03)
        # everyone tested sits at a band midpoint; the rest are never tested within 365 days
        self.assertEqual(set(d for d, _, _, _ in people), {45.0, 135.0, 272.5, math.inf})

    def test_recency_without_selection(self):
        # kappa = 0: P(tested within W) = 1 - 0.7 ** (W / 365); W = 90 gives
        # 1 - exp(-0.356675 * 0.246575) = 1 - exp(-0.087947) = 0.084191 (by hand)
        people = simulate_site(random.Random(2), 4000, 0.3, SCENARIOS["N0"])
        self.assertAlmostEqual(_share_within(people, 90), 0.084191, delta=0.03)
        self.assertAlmostEqual(_share_within(people, 365), 0.30, delta=0.03)
        self.assertTrue(all(0.0 < d < math.inf for d, _, _, _ in people))

    def test_outcomes_are_well_formed(self):
        people = simulate_site(random.Random(3), 2000, 0.25, SCENARIOS["TD"])
        self.assertTrue(all(0.0 < t <= 3650.0 for _, t, _, _ in people))
        self.assertEqual(set(e for _, _, e, _ in people) - {0, 1, 2}, set())
        self.assertTrue(all(0.0 < p < 1.0 for _, _, _, p in people))
        # KRT, death and censoring all occur
        self.assertEqual(set(e for _, _, e, _ in people), {0, 1, 2})

    def test_bad_parameters(self):
        with self.assertRaises(ValueError):
            simulate_site(random.Random(0), 10, 0.0, SCENARIOS["T69"])
        with self.assertRaises(ValueError):
            simulate_site(random.Random(0), 10, 0.3, dict(SCENARIOS["T69"], beta_x=0.6))
        with self.assertRaises(ValueError):
            simulate_site(random.Random(0), 10, 0.3, dict(SCENARIOS["T69"], mechanism="other"))

    def test_scenarios_theta(self):
        self.assertAlmostEqual(SCENARIOS["T69"]["g"] * SCENARIOS["T69"]["rho"], 0.69, places=12)
        self.assertAlmostEqual(SCENARIOS["T30"]["g"] * SCENARIOS["T30"]["rho"], 0.30, places=12)


# (days since latest uACR, time, event, predicted risk)
HAND = ([(45.0, 100.0, 1, 0.25), (45.0, 200.0, 0, 0.25), (45.0, 300.0, 1, 0.25),
         (45.0, 2000.0, 0, 0.25)]
        + [(135.0, 150.0, 1, 0.5), (135.0, 400.0, 0, 0.5), (135.0, 500.0, 1, 0.5),
           (135.0, 1900.0, 0, 0.5)]
        + [(math.inf, 500.0, 0, 0.3), (math.inf, 500.0, 0, 0.3)])


class TestHandSite(unittest.TestCase):
    def test_rows(self):
        rbw = rows_by_window(HAND, (90, 180))
        self.assertEqual([r["tested"] for r in rbw[0]], [True] * 4 + [False] * 6)
        self.assertEqual([r["tested"] for r in rbw[1]], [True] * 8 + [False] * 2)
        self.assertIsNone(rbw[0][4]["kfre5"])
        self.assertEqual(rbw[1][4]["kfre5"], 0.5)

    def test_site(self):
        s = analyse_site(HAND, (90, 180), min_events=2)
        # Band 1 (Aalen-Johansen, 1826 days): t=100 1 of 4 -> .25, S .75; t=200 censored;
        # t=300 1 of 2 -> .25 + .75/2 = .625; O/E = .625 / .25 = 2.5. Band 2: 1.25.
        self.assertTrue(s["usable"])
        self.assertTrue(s["full"])
        self.assertEqual(s["cum_rates"], [0.4, 0.8])
        self.assertAlmostEqual(s["band_oe"][0], 2.5, places=12)
        self.assertAlmostEqual(s["band_oe"][1], 1.25, places=12)
        self.assertEqual(s["band_events"], [2, 2])
        # All 8 tested: t100 1/8 -> .125 (S .875); t150 .875/7 -> .25 (S .75); 200 censored;
        # t300 .75/5 -> .40 (S .6); 400 censored; t500 .6/3 -> .60. E = 3/8. O/E = 1.6.
        self.assertAlmostEqual(s["share_tested"], 0.8, places=12)
        self.assertAlmostEqual(s["oe_tested"], 1.6, places=9)
        # no deaths, so 1 - Kaplan-Meier equals Aalen-Johansen
        self.assertEqual(s["km"]["cum_rates"], [0.4, 0.8])
        self.assertAlmostEqual(s["km"]["log_oe"][0], math.log(2.5), places=12)
        self.assertAlmostEqual(s["km"]["log_oe"][1], math.log(1.25), places=12)
        # same bands as test_nested: theta near 0.53
        self.assertGreater(s["site_theta"], 0.4)
        self.assertLess(s["site_theta"], 0.7)

    def test_death_is_competing_not_censoring(self):
        # turn the censoring at t=200 in band 1 into a death: Aalen-Johansen drops to .5
        # (S after 200 is .75 * 2/3 = .5, then t300 1 of 2 -> .25 + .5/2 = .5); 1 - KM
        # treats the death as censoring and keeps .625. So AJ O/E 2.0, KM O/E 2.5.
        people = list(HAND)
        people[1] = (45.0, 200.0, 2, 0.25)
        s = analyse_site(people, (90, 180), min_events=2)
        self.assertAlmostEqual(s["band_oe"][0], 2.0, places=12)
        self.assertAlmostEqual(s["km"]["log_oe"][0], math.log(2.5), places=12)


class TestReplicate(unittest.TestCase):
    def test_small_replicate_is_deterministic(self):
        a = run_replicate(random.Random(7), SCENARIOS["T69"], 3000, (0.2, 0.4, 0.6))
        b = run_replicate(random.Random(7), SCENARIOS["T69"], 3000, (0.2, 0.4, 0.6))
        self.assertEqual(a, b)
        self.assertEqual(len(a["share_tested"]), 3)
        self.assertLessEqual(a["usable_sites"], 3)
        for key in ("theta_w", "ci_w", "theta_b", "ci_b", "z", "theta_w_km", "q_min"):
            self.assertIn(key, a)


REP_A = {"theta_w": 0.6, "ci_w": [0.5, 0.8], "q_min": 30.0, "df": 15, "edge_w": False,
         "theta_b": 0.7, "ci_b": [0.6, 0.8], "z": -1.0, "consistent": True, "theta_w_km": 0.62,
         "usable_sites": 8, "share_tested": [0.15, 0.30],
         "site_cum_rates": [[0.10, 0.20], None], "site_band_oe": [[2.0, 1.5], None]}
REP_B = {"theta_w": 0.8, "ci_w": [0.75, 0.9], "q_min": 10.0, "df": 15, "edge_w": False,
         "theta_b": 0.6, "ci_b": [0.5, 0.7], "z": 2.5, "consistent": False, "theta_w_km": None,
         "usable_sites": 7, "share_tested": [0.17, 0.32],
         "site_cum_rates": [[0.12, 0.22], [0.3, 0.5]], "site_band_oe": [[2.2, 1.3], [1.4, 1.1]]}


class TestSummary(unittest.TestCase):
    def test_chi2(self):
        # Wilson-Hilferty by hand: df 15 -> 24.985 (exact 24.996); df 30 -> 43.767 (43.773)
        self.assertAlmostEqual(chi2_95(15), 24.996, delta=0.05)
        self.assertAlmostEqual(chi2_95(30), 43.773, delta=0.05)
        with self.assertRaises(ValueError):
            chi2_95(0)

    def test_summarize(self):
        s = summarize([REP_A, REP_B], 0.69)
        self.assertAlmostEqual(s["mean_theta_w"], 0.7, places=12)
        self.assertAlmostEqual(s["sd_theta_w"], math.sqrt(0.02), places=12)
        self.assertEqual(s["coverage_w"], 0.5)          # [.5,.8] covers .69, [.75,.9] does not
        self.assertEqual(s["reject0_w"], 1.0)
        self.assertAlmostEqual(s["mean_theta_b"], 0.65, places=12)
        self.assertEqual(s["coverage_b"], 1.0)
        self.assertEqual(s["share_w_below_b"], 0.5)     # .6 < .7 yes, .8 < .6 no
        self.assertEqual(s["consistency_reject"], 0.5)
        self.assertEqual(s["misfit_share"], 0.5)        # 30 > 24.98, 10 not
        self.assertAlmostEqual(s["mean_theta_w_km"], 0.62, places=12)
        self.assertEqual((s["mean_usable_sites"], s["min_usable_sites"]), (7.5, 7))
        for got, want in zip(s["mean_share_tested"], [0.16, 0.31]):
            self.assertAlmostEqual(got, want, places=12)
        for got, want in zip(s["mean_site_cum_rates"][0], [0.11, 0.21]):
            self.assertAlmostEqual(got, want, places=12)
        self.assertEqual(s["mean_site_cum_rates"][1], [0.3, 0.5])
        for got, want in zip(s["mean_site_band_oe"][0], [2.1, 1.4]):
            self.assertAlmostEqual(got, want, places=12)
        self.assertEqual((s["analysed_within"], s["analysed_cross"]), (2, 2))

    def test_grade(self):
        good = {
            "T69": {"mean_theta_w": 0.65, "coverage_w": 0.93, "consistency_reject": 0.05,
                    "mean_theta_b": 0.68, "coverage_b": 0.94},
            "T30": {"mean_theta_w": 0.28, "coverage_w": 0.95, "consistency_reject": 0.04,
                    "mean_theta_b": 0.30, "coverage_b": 0.95},
            "N0": {"reject0_w": 0.05, "mean_theta_w": 0.01},
            "R69": {"mean_theta_w": 0.2, "mean_theta_b": 0.4, "share_w_below_b": 0.9},
            "TD": {"mean_theta_w": 0.6},
        }
        self.assertEqual(set(grade(good).values()), {"PASS"})
        bad = {k: dict(v) for k, v in good.items()}
        bad["T69"]["coverage_w"] = 0.80
        bad["TD"]["mean_theta_w"] = 0.66
        bad["R69"]["share_w_below_b"] = 0.7
        g = grade(bad)
        self.assertEqual((g["P1"], g["P4"], g["P7"]), ("FAIL", "FAIL", "FAIL"))
        self.assertEqual((g["P2"], g["P3"], g["P5"], g["P6"]), ("PASS",) * 4)
        bad["N0"]["mean_theta_w"] = None
        self.assertEqual(grade(bad)["P3"], "FAIL")


if __name__ == "__main__":
    unittest.main()
