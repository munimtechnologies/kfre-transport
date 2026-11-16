# SPDX-License-Identifier: Apache-2.0
# Every expected value is worked out by hand in the comments.
import math
import unittest

from kfretransport.metrics import (brier, c_statistic, calibration_in_the_large,
                                   calibration_slope, dersimonian_laird,
                                   meta_regression, observed_expected, t975)

L3 = math.log(3.0)  # logit(0.75) = -logit(0.25) = 1.0986123
P8 = [0.25] * 4 + [0.75] * 4


class TestSimpleMetrics(unittest.TestCase):
    def test_oe(self):
        self.assertAlmostEqual(observed_expected([1, 0, 0, 1], [0.5, 0.25, 0.25, 0.5]), 2 / 1.5)

    def test_brier(self):
        self.assertAlmostEqual(brier([1, 0], [0.8, 0.4]), 0.1)

    def test_c_statistic(self):
        # events 0.35, 0.8 vs non-events 0.1, 0.4: 3 of 4 pairs ordered
        self.assertAlmostEqual(c_statistic([0, 0, 1, 1], [0.1, 0.4, 0.35, 0.8]), 0.75)
        self.assertAlmostEqual(c_statistic([0, 1], [0.5, 0.5]), 0.5)
        self.assertIsNone(c_statistic([1, 1], [0.2, 0.3]))

    def test_t975(self):
        self.assertEqual(t975(1), 12.706205)
        self.assertEqual(t975(30), 2.042272)
        self.assertAlmostEqual(t975(60), 2.000298, delta=0.001)


class TestCalibration(unittest.TestCase):
    def test_perfect_two_groups(self):
        # observed rates 1/4 and 3/4 equal the predictions: a = 0, b = 1.
        # weights 0.1875 x 8: h00 = 1.5, h01 = 0, h11 = 1.5 L^2, det = 2.25 L^2
        # se_a = sqrt(1.5/2.25) = 0.8164966; se_b = 0.8164966 / L = 0.7432072
        a, b, se_a, se_b = calibration_slope([1, 0, 0, 0, 1, 1, 1, 0], P8)
        self.assertAlmostEqual(a, 0.0, places=9)
        self.assertAlmostEqual(b, 1.0, places=9)
        self.assertAlmostEqual(se_a, 0.8164966, places=5)
        self.assertAlmostEqual(se_b, 0.7432072, places=5)

    def test_flattened_two_groups(self):
        # observed 2/4 and 3/4 at logits -L and +L:
        # a - bL = 0 and a + bL = L  ->  a = L/2 = 0.5493061, b = 0.5
        # h00 = 1.75, h01 = -0.25 L, h11 = 1.75 L^2, det = 3 L^2
        # se_a = sqrt(1.75/3) = 0.7637626; se_b = sqrt(1.75/(3 L^2)) = 0.6952067
        a, b, se_a, se_b = calibration_slope([1, 1, 0, 0, 1, 1, 1, 0], P8)
        self.assertAlmostEqual(a, L3 / 2, places=7)
        self.assertAlmostEqual(b, 0.5, places=7)
        self.assertAlmostEqual(se_a, 0.7637626, places=5)
        self.assertAlmostEqual(se_b, 0.6952067, places=4)

    def test_calibration_in_the_large(self):
        # predictions 0.5, observed 3/4: a = logit(0.75) = ln 3; se = sqrt(1/(4*0.1875))
        a, se = calibration_in_the_large([1, 1, 1, 0], [0.5] * 4)
        self.assertAlmostEqual(a, L3, places=9)
        self.assertAlmostEqual(se, 1.1547005, places=6)


class TestPooling(unittest.TestCase):
    def test_dl_homogeneous(self):
        # fixed mean 1, Q = 2 = k-1 -> tau2 0, I2 0; se = sqrt(1/3) = 0.5773503
        # PI: 1 -/+ 12.706205 * 0.5773503 = 1 -/+ 7.3359309
        r = dersimonian_laird([0.0, 1.0, 2.0], [1.0, 1.0, 1.0])
        self.assertAlmostEqual(r["pooled"], 1.0)
        self.assertAlmostEqual(r["tau2"], 0.0)
        self.assertAlmostEqual(r["i2"], 0.0)
        self.assertAlmostEqual(r["se"], 0.5773503, places=6)
        self.assertAlmostEqual(r["pi"][0], -6.3359309, places=4)
        self.assertAlmostEqual(r["pi"][1], 8.3359309, places=4)

    def test_dl_heterogeneous(self):
        # Q = 4 + 0 + 4 = 8, C = 3 - 3/3 = 2, tau2 = (8-2)/2 = 3, I2 = 6/8
        # RE weights 1/4: pooled 2, se = sqrt(4/3) = 1.1547005
        # PI half-width 12.706205 * sqrt(3 + 4/3) = 26.4500749
        r = dersimonian_laird([0.0, 2.0, 4.0], [1.0, 1.0, 1.0])
        self.assertAlmostEqual(r["q"], 8.0)
        self.assertAlmostEqual(r["tau2"], 3.0)
        self.assertAlmostEqual(r["i2"], 0.75)
        self.assertAlmostEqual(r["pooled"], 2.0)
        self.assertAlmostEqual(r["se"], 1.1547005, places=6)
        self.assertAlmostEqual(r["pi"][1] - 2.0, 26.4500749, places=3)

    def test_dl_two_sites_has_no_pi(self):
        self.assertIsNone(dersimonian_laird([0.0, 1.0], [1.0, 1.0])["pi"])


class TestMetaRegression(unittest.TestCase):
    def test_exact_line(self):
        # y = 1 + 2x exactly: b0 = 1, b1 = 2, residuals 0, tau2 0
        r = meta_regression([1.0, 3.0, 5.0], [1.0, 1.0, 1.0], [0.0, 1.0, 2.0])
        self.assertAlmostEqual(r["b0"], 1.0)
        self.assertAlmostEqual(r["b1"], 2.0)
        self.assertAlmostEqual(r["tau2"], 0.0)
        self.assertTrue(r["reject"])

    def test_no_excess_heterogeneity(self):
        # x = 0,0,1,1; y = 0,1,1,2: b0 0.5, b1 1, Q_res 1 < df 2 -> tau2 0
        # (X'WX)^-1[1,1] = 1; KH s2 = 1/2; se_b1 = 0.7071068; t = 1.4142136 < 4.302653
        r = meta_regression([0.0, 1.0, 1.0, 2.0], [1.0] * 4, [0.0, 0.0, 1.0, 1.0])
        self.assertAlmostEqual(r["b0"], 0.5)
        self.assertAlmostEqual(r["b1"], 1.0)
        self.assertAlmostEqual(r["q_res"], 1.0)
        self.assertAlmostEqual(r["tau2"], 0.0)
        self.assertAlmostEqual(r["se_b1"], 0.7071068, places=6)
        self.assertAlmostEqual(r["t"], 1.4142136, places=6)
        self.assertEqual(r["df"], 2)
        self.assertFalse(r["reject"])
        self.assertAlmostEqual(r["ci_b1"][1], 1.0 + 4.302653 * 0.7071068, places=5)

    def test_excess_heterogeneity(self):
        # y = 0,2,1,3: b1 1, Q_res 4, tr = 2, denom 2, tau2 = (4-2)/2 = 1
        # RE weights 1/2: (X'W*X)^-1[1,1] = 2; KH s2 = 1; se_b1 = sqrt 2
        r = meta_regression([0.0, 2.0, 1.0, 3.0], [1.0] * 4, [0.0, 0.0, 1.0, 1.0])
        self.assertAlmostEqual(r["tau2"], 1.0)
        self.assertAlmostEqual(r["b1"], 1.0)
        self.assertAlmostEqual(r["se_b1"], math.sqrt(2.0), places=9)

    def test_needs_spread(self):
        with self.assertRaises(ValueError):
            meta_regression([0.0, 1.0, 2.0], [1.0] * 3, [0.5, 0.5, 0.5])


if __name__ == "__main__":
    unittest.main()
