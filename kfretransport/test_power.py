# SPDX-License-Identifier: Apache-2.0
import math
import unittest
from statistics import NormalDist

from kfretransport.power import (curve_noncentrality, curve_power, design_table, even_rates,
                                 linear_noncentrality, linear_power, main, mills,
                                 monte_carlo_power, power_from_noncentrality, sites_needed,
                                 theta_information, theta_se)
from kfretransport.selection import tested_factor

N = NormalDist()


class TestMills(unittest.TestCase):
    def test_hand_values(self):
        # phi(0) / 0.5 = 0.3989423 / 0.5
        self.assertAlmostEqual(mills(0.5), 0.7978846, places=6)
        # phi(1.2815516) = 0.175498 (normal tables); / 0.1
        self.assertAlmostEqual(mills(0.1), 1.75498, places=3)
        # phi(0.5244005) = 0.347693; / 0.7
        self.assertAlmostEqual(mills(0.7), 0.496704, places=4)

    def test_is_derivative_of_log_curve(self):
        h = 1e-5
        for p, th in ((0.1, 0.0), (0.3, 0.4), (0.8, -0.3)):
            fd = (math.log(tested_factor(p, th + h)) - math.log(tested_factor(p, th - h))) / (2 * h)
            self.assertAlmostEqual(mills(p, th), fd, places=6)

    def test_decreasing_in_rate(self):
        vals = [mills(0.05 * i) for i in range(1, 20)]
        for a, b in zip(vals, vals[1:]):
            self.assertGreater(a, b)


class TestAnalyticPower(unittest.TestCase):
    def test_information_two_point_design(self):
        # two sites at p = 0.1 and two at p = 0.7, variance 0.04 each (weight 25):
        # info = 4 * 25 * ((m(0.1) - m(0.7)) / 2)^2 = 100 * 0.629138^2 = 39.58 by hand
        rates = [0.1, 0.1, 0.7, 0.7]
        info = theta_information(rates, [0.04] * 4)
        exact = 100.0 * ((mills(0.1) - mills(0.7)) / 2.0) ** 2
        self.assertAlmostEqual(info, exact, places=9)
        self.assertAlmostEqual(info, 39.58, delta=0.05)
        self.assertAlmostEqual(theta_se(rates, [0.04] * 4), 1.0 / math.sqrt(exact), places=9)
        # tau2 adds to every variance: weight 1 / 0.05 = 20
        self.assertAlmostEqual(theta_information(rates, [0.04] * 4, tau2=0.01), exact * 0.8, places=9)

    def test_null_power_is_alpha(self):
        rates = even_rates(8, 0.1, 0.7)
        v = [0.02] * 8
        self.assertAlmostEqual(curve_power(rates, v, 0.05, 0.0), 0.05, places=9)
        self.assertAlmostEqual(linear_power(rates, v, 0.05, 0.0), 0.05, places=9)
        self.assertAlmostEqual(power_from_noncentrality(0.0, 0.1), 0.1, places=12)

    def test_power_hand_value(self):
        # noncentrality (1.96 + 0.8416)^2 gives power 0.80 (plus a negligible lower tail)
        nc = (N.inv_cdf(0.975) + N.inv_cdf(0.8)) ** 2
        self.assertAlmostEqual(power_from_noncentrality(nc), 0.8, places=4)

    def test_linear_never_beats_curve(self):
        for lo, hi in ((0.1, 0.7), (0.2, 0.5), (0.05, 0.95)):
            rates = even_rates(9, lo, hi)
            for th in (-0.4, 0.1, 0.3, 0.69, 1.2):
                lin = linear_noncentrality(rates, [0.02] * 9, 0.03, th)
                cur = curve_noncentrality(rates, [0.02] * 9, 0.03, th)
                self.assertLessEqual(lin, cur * (1 + 1e-12) + 1e-15)

    def test_small_theta_matches_information(self):
        rates = even_rates(7, 0.1, 0.7)
        v = [0.02] * 7
        th = 0.001
        nc = curve_noncentrality(rates, v, 0.0, th)
        approx = th * th * theta_information(rates, v, 0.0, 0.0)
        self.assertAlmostEqual(nc / approx, 1.0, delta=0.01)

    def test_monotone(self):
        rates = even_rates(10, 0.1, 0.7)
        v = [0.02] * 10
        p1 = curve_power(rates, v, 0.05, 0.1)
        p2 = curve_power(rates, v, 0.05, 0.2)
        p3 = curve_power(rates, v, 0.05, 0.4)
        self.assertLess(p1, p2)
        self.assertLess(p2, p3)
        self.assertGreater(curve_power(rates, v, 0.0, 0.2), p2)
        # a narrower spread of testing rates gives less power
        self.assertLess(curve_power(even_rates(10, 0.2, 0.5), v, 0.05, 0.2), p2)


class TestSitesNeeded(unittest.TestCase):
    def test_bracket(self):
        k = sites_needed(0.3, 0.02, 0.05)
        self.assertIsNotNone(k)
        self.assertGreater(k, 3)
        good = curve_power(even_rates(k, 0.1, 0.7), [0.02] * k, 0.05, 0.3)
        bad = curve_power(even_rates(k - 1, 0.1, 0.7), [0.02] * (k - 1), 0.05, 0.3)
        self.assertGreaterEqual(good, 0.8)
        self.assertLess(bad, 0.8)
        # linearised hand value: 7.849 * 0.07 / (0.09 * 0.1160) = 52.6 sites
        # (variance of m(p) over p uniform on [0.1, 0.7] by Simpson's rule);
        # evenly spaced sites do a little better than a continuous spread, and the
        # curve's bend at theta = 0.3 costs about 8% of the information (hand: ~55)
        self.assertGreater(k, 40)
        self.assertLess(k, 75)

    def test_linear_needs_at_least_as_many(self):
        for th in (0.3, 0.69):
            kc = sites_needed(th, 0.02, 0.05)
            kl = sites_needed(th, 0.02, 0.05, test="linear")
            self.assertGreaterEqual(kl, kc)

    def test_none_when_too_few(self):
        self.assertIsNone(sites_needed(0.05, 0.04, 0.1, kmax=10))

    def test_errors(self):
        with self.assertRaises(ValueError):
            sites_needed(0.3, 0.02, test="cubic")
        with self.assertRaises(ValueError):
            even_rates(5, 0.5, 0.2)
        with self.assertRaises(ValueError):
            theta_information([0.2, 0.4], [0.01, 0.01])
        with self.assertRaises(ValueError):
            theta_information([0.3, 0.3, 0.3], [0.01] * 3)
        with self.assertRaises(ValueError):
            mills(1.0)

    def test_table_and_main(self):
        rows = design_table([0.69], [0.02], [0.05], [[0.1, 0.7]])
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["sites_curve"], sites_needed(0.69, 0.02, 0.05))
        self.assertEqual(main(["--theta", "0.69", "--site-variance", "0.02", "--tau2", "0.05",
                               "--spread", "0.1", "0.7"]), 0)


class TestMonteCarlo(unittest.TestCase):
    def test_null_rejection_near_alpha(self):
        rates = even_rates(12, 0.1, 0.7)
        r = monte_carlo_power(rates, [0.02] * 12, 0.0, 0.0, reps=300, seed=7, lo=-1.0, hi=1.5)
        self.assertGreater(r, 0.01)
        self.assertLess(r, 0.11)

    def test_matches_analytic(self):
        rates = even_rates(20, 0.1, 0.7)
        v = [0.02] * 20
        analytic = curve_power(rates, v, 0.0, 0.15)
        self.assertGreater(analytic, 0.25)
        self.assertLess(analytic, 0.7)
        mc = monte_carlo_power(rates, v, 0.0, 0.15, beta0=0.2, reps=300, seed=11, lo=-1.0, hi=1.5)
        self.assertAlmostEqual(mc, analytic, delta=0.12)

    def test_deterministic(self):
        rates = even_rates(6, 0.1, 0.7)
        a = monte_carlo_power(rates, [0.02] * 6, 0.0, 0.3, reps=20, seed=3, lo=-1.0, hi=1.5)
        b = monte_carlo_power(rates, [0.02] * 6, 0.0, 0.3, reps=20, seed=3, lo=-1.0, hi=1.5)
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
