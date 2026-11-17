# SPDX-License-Identifier: Apache-2.0
# Every expected value is worked out by hand in the comments.
import math
import unittest

from kfretransport.metrics import (brier, c_statistic, calibration_in_the_large,
                                   calibration_slope, cumulative_incidence,
                                   cumulative_incidence_variance,
                                   dersimonian_laird, km_risk, meta_regression,
                                   observed_expected, oe_time_to_event, t975)

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


# Five people: KRT at 2, death at 3, censored 4, KRT at 5, censored 6.
T5 = [2, 3, 4, 5, 6]
E5 = [1, 2, 0, 1, 0]


class TestTimeToEvent(unittest.TestCase):
    def test_aalen_johansen(self):
        # t=2: 5 at risk, 1 KRT -> CIF 1/5 = 0.2; S = 4/5
        # t=3: 4 at risk, 1 death -> CIF unchanged; S = 0.8 * 3/4 = 0.6
        # t=4: censored, 3 -> 2 at risk
        # t=5: 2 at risk, 1 KRT -> CIF += 0.6 * 1/2 = 0.3 -> 0.5
        cif, n = cumulative_incidence(T5, E5, 10)
        self.assertAlmostEqual(cif, 0.5, places=12)
        self.assertEqual(n, 2)
        cif, n = cumulative_incidence(T5, E5, 4.5)
        self.assertAlmostEqual(cif, 0.2, places=12)
        self.assertEqual(n, 1)
        # an event exactly at the horizon counts
        self.assertAlmostEqual(cumulative_incidence(T5, E5, 5)[0], 0.5, places=12)
        # cause 2 (death): only t=3, S(3-) = 0.8, 4 at risk -> 0.8 * 1/4 = 0.2
        self.assertAlmostEqual(cumulative_incidence(T5, E5, 10, cause=2)[0], 0.2, places=12)

    def test_kaplan_meier_death_censored(self):
        # death at 3 treated as censoring: t=2 S = 4/5; t=5 2 at risk, S = 0.8 * 1/2
        # risk = 1 - 0.4 = 0.6 (> 0.5 from Aalen-Johansen)
        risk, n = km_risk(T5, E5, 10)
        self.assertAlmostEqual(risk, 0.6, places=12)
        self.assertEqual(n, 2)

    def test_ties(self):
        # KRT and death both at t=1 with 4 at risk: CIF 1/4, S = 1 - 2/4 = 0.5;
        # t=3: 2 at risk (times 3 and 4), CIF += 0.5 * 1/2 -> 0.5
        cif, _ = cumulative_incidence([1, 1, 3, 4], [1, 2, 1, 0], 5)
        self.assertAlmostEqual(cif, 0.5, places=12)
        # KM: t=1 S = 3/4; the death leaves; t=3 2 at risk, S = 0.375 -> risk 0.625
        self.assertAlmostEqual(km_risk([1, 1, 3, 4], [1, 2, 1, 0], 5)[0], 0.625, places=12)
        # someone censored at the same time as an event is still at risk then
        self.assertAlmostEqual(cumulative_incidence([2, 2], [1, 0], 3)[0], 0.5, places=12)

    def test_oe_time_to_event(self):
        # predictions all 0.25: E = 0.25. AJ: O = 0.5 -> O/E 2.0; KM: 0.6 -> 2.4
        # 2 KRT events: se(ln O/E) = 1/sqrt(2) = 0.7071068
        r = oe_time_to_event(T5, E5, [0.25] * 5, 10)
        self.assertAlmostEqual(r["oe"], 2.0, places=12)
        self.assertAlmostEqual(r["observed_risk"], 0.5, places=12)
        self.assertAlmostEqual(r["expected_risk"], 0.25, places=12)
        self.assertEqual((r["n"], r["events"]), (5, 2))
        self.assertAlmostEqual(r["se_log_oe"], 0.7071068, places=6)
        r = oe_time_to_event(T5, E5, [0.25] * 5, 10, competing=False)
        self.assertAlmostEqual(r["oe"], 2.4, places=12)
        self.assertEqual(r["events"], 2)
        r = oe_time_to_event(T5, E5, [0.25] * 5, 1)
        self.assertEqual((r["oe"], r["events"], r["se_log_oe"]), (0.0, 0, None))
        self.assertIsNone(r["se_log_oe_dm"])
        with self.assertRaises(ValueError):
            oe_time_to_event(T5, E5, [0.25] * 4, 10)


class TestAalenJohansenVariance(unittest.TestCase):
    def test_worked_example(self):
        # event times (Y, d, d1, S-, F after): t=2 (5,1,1,1.0,0.2), t=3 (4,1,0,0.8,0.2),
        # t=5 (2,1,1,0.6,0.5); F(10) = 0.5
        # term 1: 0.3^2*1/(5*4) + 0.3^2*1/(4*3) + 0 = 0.0045 + 0.0075 = 0.012
        # term 2: 1*1*4/125 + 0 + 0.36*1*1/8 = 0.032 + 0.045 = 0.077
        # term 3: -2 * (0.3*1*1/25 + 0 + 0) = -0.024
        # var = 0.065; se = 0.2549510; se(ln O) = 0.2549510 / 0.5 = 0.5099020
        cif, var, n = cumulative_incidence_variance(T5, E5, 10)
        self.assertAlmostEqual(cif, 0.5, places=12)
        self.assertAlmostEqual(var, 0.065, places=12)
        self.assertEqual(n, 2)
        r = oe_time_to_event(T5, E5, [0.25] * 5, 10)
        self.assertAlmostEqual(r["var_observed_risk"], 0.065, places=12)
        self.assertAlmostEqual(r["se_log_oe_dm"], 0.5099020, places=6)

    def test_no_censoring_is_multinomial(self):
        # 10 people: t=1 2 KRT + 3 deaths; t=2 1 KRT; 4 followed past the horizon.
        # F = 0.2 + 0.5*1/5 = 0.3; multinomial variance 0.3*0.7/10 = 0.021
        # term 1: 0.1^2 * 5/(10*5) = 0.001; term 2: 2*8/1000 + 0.25*1*4/125 = 0.024;
        # term 3: -2 * 0.1 * 1 * 2/100 = -0.004; total 0.021
        times = [1] * 5 + [2] + [9] * 4
        events = [1, 1, 2, 2, 2, 1, 0, 0, 0, 0]
        cif, var, n = cumulative_incidence_variance(times, events, 5)
        self.assertAlmostEqual(cif, 0.3, places=12)
        self.assertAlmostEqual(var, 0.021, places=12)
        self.assertEqual(n, 3)

    def test_no_competing_is_greenwood(self):
        # death censored (1 - KM): S = 0.4, Greenwood 0.16 * (1/(5*4) + 1/(2*1)) = 0.088
        # delta form: 0.4^2/20 + (4/125 + 0.64/8) - 2*0.4/25 = 0.008 + 0.112 - 0.032
        r = oe_time_to_event(T5, E5, [0.25] * 5, 10, competing=False)
        self.assertAlmostEqual(r["var_observed_risk"], 0.088, places=12)
        # sqrt(0.088) = 0.2966479; / 0.6 = 0.4944132
        self.assertAlmostEqual(r["se_log_oe_dm"], 0.4944132, places=6)

    def test_edges(self):
        # one of two at risk has KRT, the other censored then: 1*1*1/8 = 0.125
        self.assertAlmostEqual(cumulative_incidence_variance([2, 2], [1, 0], 3)[1], 0.125,
                               places=12)
        # everyone at risk has the event: Y = d, F = 1, variance 0
        cif, var, _ = cumulative_incidence_variance([1, 1], [1, 1], 2)
        self.assertEqual((cif, var), (1.0, 0.0))
        with self.assertRaises(ValueError):
            cumulative_incidence_variance([], [], 1)


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
