# SPDX-License-Identifier: Apache-2.0
# Expected values below were worked out by hand (series for ln, products of
# exp terms); the working is in the comments so anyone can re-check it.
import unittest

from kfretransport.equations import (acr_mg_mmol_to_mg_g, albuminuria_category,
                                     creatinine_umol_to_mg_dl, egfr_ckdepi_2009,
                                     egfr_ckdepi_2021, gfr_category, kfre4,
                                     kfre4_linear_predictor, risk_from_lp)


class TestEgfr(unittest.TestCase):
    def test_2021_male_at_kappa(self):
        # Scr/kappa = 1, so only 142 * 0.9938**40 remains.
        # ln 0.9938 = -0.0062193; x40 = -0.2487720; exp = 0.7797577; x142 = 110.7256
        self.assertAlmostEqual(egfr_ckdepi_2021(0.9, 40, female=False), 110.7256, delta=0.01)

    def test_2021_female_at_kappa_is_male_times_1012(self):
        f = egfr_ckdepi_2021(0.7, 40, female=True)
        m = egfr_ckdepi_2021(0.9, 40, female=False)
        self.assertAlmostEqual(f, 112.0543, delta=0.01)
        self.assertAlmostEqual(f / m, 1.012, places=9)

    def test_2021_female_60_scr1(self):
        # 142 * 0.7**1.2 (0.6518049) * 0.9938**60 (0.6885560) * 1.012 = 64.4950
        self.assertAlmostEqual(egfr_ckdepi_2021(1.0, 60, female=True), 64.495, delta=0.02)

    def test_2009_male_at_kappa(self):
        # 141 * 0.993**40: ln 0.993 = -0.00702461; x40 = -0.2809844; exp = 0.7550401
        self.assertAlmostEqual(egfr_ckdepi_2009(0.9, 40, female=False), 106.4607, delta=0.02)

    def test_2009_female_black_factors(self):
        base = egfr_ckdepi_2009(0.9, 40, female=False)
        fb = egfr_ckdepi_2009(0.7, 40, female=True, black=True)
        self.assertAlmostEqual(fb / base, 1.018 * 1.159 * 141.0 / 141.0, places=9)

    def test_higher_creatinine_lowers_egfr(self):
        a = egfr_ckdepi_2021(1.0, 60, female=False)
        b = egfr_ckdepi_2021(2.0, 60, female=False)
        self.assertLess(b, a)

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            egfr_ckdepi_2021(0.0, 50, female=True)
        with self.assertRaises(ValueError):
            egfr_ckdepi_2021(1.0, 17, female=True)

    def test_unit_conversion(self):
        self.assertAlmostEqual(creatinine_umol_to_mg_dl(88.4), 1.0, places=12)
        self.assertAlmostEqual(acr_mg_mmol_to_mg_g(3.0), 26.52, places=9)


class TestCategories(unittest.TestCase):
    def test_gfr(self):
        got = [gfr_category(x) for x in (95, 90, 89.9, 60, 59.9, 45, 44.9, 30, 29.9, 15, 14.9)]
        self.assertEqual(got, ["G1", "G1", "G2", "G2", "G3a", "G3a", "G3b", "G3b", "G4", "G4", "G5"])

    def test_albuminuria(self):
        got = [albuminuria_category(x) for x in (5, 29.9, 30, 300, 300.1)]
        self.assertEqual(got, ["A1", "A1", "A2", "A2", "A3"])


class TestKfre(unittest.TestCase):
    def test_linear_predictor_hand_case(self):
        # 80-year-old male, eGFR 20, uACR 900 mg/g:
        # -0.2201*(8-7.036)    = -0.2121764
        # +0.2467*(1-0.5642)   = +0.1075119
        # -0.5567*(4-7.222)    = +1.7936874
        # +0.4510*(ln900-5.137) = +0.4510*1.6653948 = +0.7510930
        # sum = 2.4401159
        self.assertAlmostEqual(kfre4_linear_predictor(80, True, 20, 900), 2.4401159, places=6)

    def test_risk_hand_case_north_america(self):
        # exp(2.4401159) = 11.474371
        # 2 y: 1 - exp(11.474371 * ln 0.975) = 1 - exp(-0.2905059) = 0.2521149
        # 5 y: 1 - exp(11.474371 * ln 0.924) = 1 - exp(-0.9069710) = 0.5962547
        self.assertAlmostEqual(kfre4(80, True, 20, 900, years=2), 0.25211, delta=0.0002)
        self.assertAlmostEqual(kfre4(80, True, 20, 900, years=5), 0.59625, delta=0.0002)

    def test_lp_zero_gives_one_minus_baseline(self):
        self.assertAlmostEqual(risk_from_lp(0.0, 5, "north_america"), 0.076, places=12)
        self.assertAlmostEqual(risk_from_lp(0.0, 2, "non_north_america"), 0.0168, places=12)

    def test_monotone(self):
        self.assertGreater(kfre4(60, False, 30, 300), kfre4(60, False, 30, 30))
        self.assertLess(kfre4(60, False, 45, 300), kfre4(60, False, 30, 300))
        self.assertGreater(kfre4(60, True, 30, 300), kfre4(60, False, 30, 300))
        self.assertGreater(kfre4(60, False, 30, 300, years=5), kfre4(60, False, 30, 300, years=2))

    def test_acr_must_be_positive(self):
        with self.assertRaises(ValueError):
            kfre4(60, False, 30, 0.0)


if __name__ == "__main__":
    unittest.main()
