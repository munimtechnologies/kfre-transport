# SPDX-License-Identifier: Apache-2.0
"""Tests for transport.py (hand values and exact identities).
Made with AI agents (Claude) working in MV2's name (MV2, Munim, Inc.)."""

import contextlib
import io
import json
import math
import os
import tempfile
import unittest
from statistics import NormalDist

from kfretransport.selection import tested_factor
from kfretransport.transport import (
    Z975,
    fit_covariance,
    leave_one_out,
    main,
    predict,
    predict_null,
    theta_score,
)

RATES = [0.1, 0.25, 0.4, 0.55, 0.7]
V = [0.01, 0.01, 0.01, 0.01, 0.01]


def curve_data(theta, b0=0.1, noise=None):
    noise = noise or [0.0] * len(RATES)
    return [b0 + math.log(tested_factor(p, theta)) + e for p, e in zip(RATES, noise)]


class TestScore(unittest.TestCase):
    def test_hand_value_at_half(self):
        # phi(0) / Phi(0) = 0.39894228 / 0.5
        self.assertAlmostEqual(theta_score(0.5, 0.0), 0.7978845608, places=9)

    def test_hand_value_at_z_one(self):
        # phi(1) / Phi(1) = 0.24197072 / 0.84134475 = 0.28760
        p = NormalDist().cdf(1.0)
        self.assertAlmostEqual(theta_score(p, 0.0), 0.28760, places=4)

    def test_is_derivative_of_log_curve(self):
        p, th, h = 0.2, 0.4, 1e-5
        fd = (math.log(tested_factor(p, th + h)) - math.log(tested_factor(p, th - h))) / (2 * h)
        self.assertAlmostEqual(theta_score(p, th), fd, places=7)

    def test_rate_checked(self):
        with self.assertRaises(ValueError):
            theta_score(0.0, 0.1)


class TestCovariance(unittest.TestCase):
    def test_matches_weighted_regression_identities(self):
        rates, var, th, tau2 = [0.1, 0.4, 0.7], [0.01, 0.02, 0.04], 0.3, 0.01
        w = [1.0 / (v + tau2) for v in var]
        s = [theta_score(p, th) for p in rates]
        sw = sum(w)
        sbar = sum(wi * si for wi, si in zip(w, s)) / sw
        sxx = sum(wi * (si - sbar) ** 2 for wi, si in zip(w, s))
        cov = fit_covariance(rates, var, th, tau2)
        self.assertAlmostEqual(cov[1][1], 1.0 / sxx, places=10)
        self.assertAlmostEqual(cov[0][1], -sbar / sxx, places=10)
        self.assertAlmostEqual(cov[1][0], cov[0][1], places=12)
        self.assertAlmostEqual(cov[0][0], 1.0 / sw + sbar * sbar / sxx, places=10)

    def test_two_point_hand_value(self):
        # s1 = 0.79788, s2 = 0.28760, w = 100 each: var theta = 1 / (50 (s1 - s2)^2) = 0.0768
        cov = fit_covariance([0.5, NormalDist().cdf(1.0)], [0.01, 0.01], 0.0)
        self.assertAlmostEqual(cov[1][1], 0.0768, places=4)

    def test_equal_rates_have_no_information(self):
        with self.assertRaises(ValueError):
            fit_covariance([0.3, 0.3, 0.3], [0.01, 0.02, 0.03], 0.5)


class TestPredict(unittest.TestCase):
    def test_recovers_exact_curve(self):
        y = curve_data(0.5)
        res = predict(RATES, y, V, 0.3)
        self.assertAlmostEqual(res["theta"], 0.5, places=9)
        self.assertAlmostEqual(res["beta0"], 0.1, places=9)
        self.assertAlmostEqual(res["log_oe"], 0.1 + math.log(tested_factor(0.3, 0.5)), places=9)
        self.assertFalse(res["theta_at_grid_edge"])

    def test_variance_parts_and_interval(self):
        y = curve_data(0.5)
        res = predict(RATES, y, V, 0.3, v_new=0.02, tau2=0.05)
        self.assertGreater(res["var_fit"], 0.0)
        self.assertAlmostEqual(res["var_pred"] - res["var_fit"], 0.07, places=12)
        sd = math.sqrt(res["var_pred"])
        self.assertAlmostEqual(math.log(res["pi95"][1]) - res["log_oe"], Z975 * sd, places=9)
        self.assertAlmostEqual(res["log_oe"] - math.log(res["pi95"][0]), Z975 * sd, places=9)

    def test_shift_invariance(self):
        y = curve_data(0.5, noise=[0.02, -0.02, 0.02, -0.02, 0.02])
        a = predict(RATES, y, V, 0.3)
        b = predict(RATES, [yi + 0.7 for yi in y], V, 0.3)
        self.assertAlmostEqual(a["theta"], b["theta"], places=9)
        self.assertAlmostEqual(b["log_oe"] - a["log_oe"], 0.7, places=9)
        self.assertAlmostEqual(b["var_pred"], a["var_pred"], places=12)

    def test_new_rate_checked(self):
        with self.assertRaises(ValueError):
            predict(RATES, curve_data(0.5), V, 1.0)


class TestNull(unittest.TestCase):
    def test_hand_value(self):
        # weights 100, 100, 50; mean = (0 + 20 + 20) / 250 = 0.16; var_fit = 1/250
        res = predict_null([0.0, 0.2, 0.4], [0.01, 0.01, 0.02], v_new=0.03)
        self.assertAlmostEqual(res["log_oe"], 0.16, places=12)
        self.assertAlmostEqual(res["var_fit"], 0.004, places=12)
        self.assertAlmostEqual(res["var_pred"], 0.034, places=12)


class TestLeaveOneOut(unittest.TestCase):
    def test_exact_curve_is_predicted_perfectly(self):
        res = leave_one_out(RATES, curve_data(0.5), V)
        self.assertEqual(res["n_sites"], 5)
        for row in res["rows"]:
            self.assertAlmostEqual(row["err_curve"], 0.0, places=7)
            self.assertAlmostEqual(row["theta_without_site"], 0.5, places=9)
        self.assertAlmostEqual(res["skill"], 1.0, places=6)
        self.assertEqual(res["coverage95_curve"], 1.0)
        self.assertGreater(res["sse_null"], 0.0)

    def test_noisy_curve_beats_null(self):
        y = curve_data(0.69, noise=[0.02, -0.02, 0.02, -0.02, 0.02])
        res = leave_one_out(RATES, y, V)
        self.assertGreater(res["skill"], 0.5)
        self.assertLess(res["sse_curve"], res["sse_null"])

    def test_null_errors_by_hand_definition(self):
        y = curve_data(0.5)
        res = leave_one_out(RATES, y, V)
        row = res["rows"][0]
        others = y[1:]
        self.assertAlmostEqual(row["pred_null"], sum(others) / 4.0, places=12)
        self.assertAlmostEqual(row["err_null"], y[0] - sum(others) / 4.0, places=12)

    def test_needs_four_sites(self):
        with self.assertRaises(ValueError):
            leave_one_out(RATES[:3], curve_data(0.5)[:3], V[:3])


class TestMain(unittest.TestCase):
    def test_cli_runs(self):
        y = curve_data(0.5)
        sites = [{"rate": p, "log_oe": yi, "variance": v} for p, yi, v in zip(RATES, y, V)]
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "sites.json")
            with open(path, "w", encoding="utf-8") as fh:
                json.dump({"sites": sites}, fh)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = main(["--sites", path, "--loo", "--predict", "0.3"])
        self.assertEqual(rc, 0)
        out = json.loads(buf.getvalue())
        self.assertEqual(out["leave_one_out"]["n_sites"], 5)
        self.assertAlmostEqual(out["predictions"][0]["theta"], 0.5, places=9)


if __name__ == "__main__":
    unittest.main()
