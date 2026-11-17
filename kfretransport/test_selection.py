# SPDX-License-Identifier: Apache-2.0
import math
import unittest
from statistics import NormalDist

from kfretransport.selection import (fit_sites, log_slope, main, risk_scale_factors,
                                     tested_factor, theta_from_oe, untested_factor)

N = NormalDist()


class TestCurve(unittest.TestCase):
    def test_no_selection_is_one(self):
        # z_0.5 = 0, Phi(0) = 0.5
        self.assertAlmostEqual(tested_factor(0.5, 0.0), 1.0, places=12)
        for p in (0.1, 0.25, 0.63, 0.9):
            self.assertAlmostEqual(tested_factor(p, 0.0), 1.0, places=9)
            self.assertAlmostEqual(untested_factor(p, 0.0), 1.0, places=9)

    def test_hand_value_half(self):
        # p = 0.5, theta = Phi^-1(0.75): tested 0.75 / 0.5, untested 0.25 / 0.5
        th = N.inv_cdf(0.75)
        self.assertAlmostEqual(tested_factor(0.5, th), 1.5, places=9)
        self.assertAlmostEqual(untested_factor(0.5, th), 0.5, places=9)

    def test_hand_values_table(self):
        # hand-computed from printed normal tables (findings note)
        self.assertAlmostEqual(tested_factor(0.25, 0.3), 1.416, places=2)
        self.assertAlmostEqual(tested_factor(0.10, 0.69), 2.771, places=2)
        self.assertAlmostEqual(tested_factor(0.63, 0.69), 1.344, places=2)

    def test_total_identity(self):
        for p in (0.05, 0.25, 0.5, 0.63, 0.95):
            for th in (-0.5, 0.2, 0.69, 1.5):
                t = tested_factor(p, th)
                u = untested_factor(p, th)
                self.assertAlmostEqual(p * t + (1 - p) * u, 1.0, places=12)

    def test_monotone_in_rate(self):
        rates = [0.05 + 0.05 * i for i in range(19)]
        vals = [tested_factor(p, 0.4) for p in rates]
        for a, b in zip(vals, vals[1:]):
            self.assertGreater(a, b)
        for v in vals:
            self.assertGreater(v, 1.0)
        for p in rates:
            self.assertLess(untested_factor(p, 0.4), 1.0)
            self.assertLess(tested_factor(p, -0.4), 1.0)

    def test_theta_roundtrip_and_scream(self):
        for p in (0.1, 0.4, 0.8):
            for th in (-0.3, 0.25, 0.9):
                self.assertAlmostEqual(theta_from_oe(p, tested_factor(p, th)), th, places=9)
        # SCREAM: 26,229 of 116,158 tested, O/E 2.11 (Cleary et al. 2026); hand value 0.6937
        self.assertAlmostEqual(theta_from_oe(26229 / 116158, 2.11), 0.694, places=2)

    def test_log_slope(self):
        self.assertAlmostEqual(log_slope(0.3, 0.0), 0.0, places=12)
        h = 1e-5
        for p, th in ((0.25, 0.3), (0.6, 0.69), (0.1, -0.2)):
            fd = (math.log(tested_factor(p + h, th)) - math.log(tested_factor(p - h, th))) / (2 * h)
            self.assertAlmostEqual(log_slope(p, th), fd, places=6)
        self.assertLess(log_slope(0.4, 0.5), 0.0)

    def test_errors(self):
        with self.assertRaises(ValueError):
            tested_factor(0.0, 0.1)
        with self.assertRaises(ValueError):
            tested_factor(1.0, 0.1)
        with self.assertRaises(ValueError):
            theta_from_oe(0.5, 2.0)
        with self.assertRaises(ValueError):
            fit_sites([0.2, 0.4], [0.1, 0.0], [0.01, 0.01])
        with self.assertRaises(ValueError):
            fit_sites([0.3, 0.3, 0.3], [0.1, 0.0, 0.2], [0.01, 0.01, 0.01])
        with self.assertRaises(ValueError):
            risk_scale_factors(0.3, 0.5, 1.0, 0.1)


class TestRiskScale(unittest.TestCase):
    def test_quadrature_share_and_identity(self):
        r = risk_scale_factors(0.3, 0.8, 0.5, 0.2)
        self.assertAlmostEqual(r["tested_share"], 0.3, places=9)
        s = r["tested_share"]
        self.assertAlmostEqual(s * r["tested"] + (1 - s) * r["untested"], 1.0, places=12)

    def test_no_correlation_is_one(self):
        r = risk_scale_factors(0.4, 1.0, 0.0, 0.3)
        self.assertAlmostEqual(r["tested"], 1.0, places=9)
        self.assertAlmostEqual(r["untested"], 1.0, places=9)

    def test_rare_outcome_matches_closed_form(self):
        for p, g, rho in ((0.25, 0.8, 0.5), (0.6, 1.2, 0.7), (0.1, 0.5, -0.4)):
            r = risk_scale_factors(p, g, rho, 1e-7)
            self.assertAlmostEqual(r["tested"], tested_factor(p, g * rho), places=4)
            self.assertAlmostEqual(r["untested"], untested_factor(p, g * rho), places=4)

    def test_common_outcome_is_attenuated(self):
        # risk saturates, so the risk-scale factor lies between 1 and the hazard-scale curve
        r = risk_scale_factors(0.25, 1.0, 0.6, 0.3)
        self.assertGreater(r["tested"], 1.0)
        self.assertLess(r["tested"], r["closed_form"])


class TestFit(unittest.TestCase):
    def test_recovers_exact_curve(self):
        rates = [0.1, 0.25, 0.4, 0.63, 0.9]
        y = [0.1 + math.log(tested_factor(p, 0.5)) for p in rates]
        f = fit_sites(rates, y, [0.01] * 5)
        self.assertAlmostEqual(f["theta"], 0.5, delta=1e-6)
        self.assertAlmostEqual(f["beta0"], 0.1, delta=1e-6)
        self.assertAlmostEqual(f["oe_all"], math.exp(0.1), delta=1e-6)
        self.assertLess(f["ci95_profile"][0], 0.5)
        self.assertGreater(f["ci95_profile"][1], 0.5)
        self.assertFalse(f["theta_at_grid_edge"])
        self.assertFalse(f["ci_at_grid_edge"])
        self.assertEqual(f["k"], 5)

    def test_flat_sites_give_theta_near_zero(self):
        rates = [0.2, 0.4, 0.6, 0.8]
        f = fit_sites(rates, [0.05] * 4, [0.02] * 4)
        self.assertAlmostEqual(f["theta"], 0.0, delta=1e-6)
        self.assertAlmostEqual(f["beta0"], 0.05, delta=1e-6)
        self.assertLess(f["ci95_profile"][0], 0.0)
        self.assertGreater(f["ci95_profile"][1], 0.0)

    def test_main_runs(self):
        self.assertEqual(main(["--theta", "0.3", "--rates", "0.25", "0.63"]), 0)


if __name__ == "__main__":
    unittest.main()
