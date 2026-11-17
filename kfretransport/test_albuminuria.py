# SPDX-License-Identifier: Apache-2.0
# The crude Sumida 2020 equation must reproduce the figures the paper reports
# for it. Hand working in the comments.
import math
import unittest

from kfretransport.albuminuria import acr_from_pcr


class TestSumidaCrude(unittest.TestCase):
    def test_knot_500(self):
        # every log term is 0 at PCR 500: ACR = exp(5.392) = 219.64, reported as 220
        self.assertAlmostEqual(acr_from_pcr(500.0), math.exp(5.392), places=9)
        self.assertEqual(round(acr_from_pcr(500.0)), 220)

    def test_pcr_150(self):
        # 5.392 + 1.5793 ln 0.3 = 5.392 - 1.9014342 = 3.4905658 -> 32.80, reported 33
        self.assertAlmostEqual(acr_from_pcr(150.0), 32.8045, delta=0.01)
        self.assertEqual(round(acr_from_pcr(150.0)), 33)

    def test_reported_thresholds(self):
        # ACR 30 <-> PCR 142: 5.392 + 1.5793 ln 0.284 = 3.4040065 -> 30.08
        self.assertAlmostEqual(acr_from_pcr(142.0), 30.0, delta=0.2)
        # ACR 300 <-> PCR 660: 5.392 + 1.1266 ln 1.32 = 5.7047799 -> 300.30
        self.assertAlmostEqual(acr_from_pcr(660.0), 300.0, delta=1.0)

    def test_fold_change_per_doubling(self):
        # 2^1.5793 = exp(1.0946874) = 2.98823, reported as 2.99 between 50 and 500
        self.assertAlmostEqual(acr_from_pcr(200.0) / acr_from_pcr(100.0), 2.99, delta=0.005)
        # above 500: 2^1.1266
        self.assertAlmostEqual(acr_from_pcr(2000.0) / acr_from_pcr(1000.0),
                               2.0 ** 1.1266, places=9)

    def test_low_segment_and_continuity(self):
        # at PCR 50 the middle term floors at ln 0.1: exp(5.392 + 1.5793 ln 0.1)
        self.assertAlmostEqual(acr_from_pcr(50.0), math.exp(5.392 + 1.5793 * math.log(0.1)),
                               places=9)
        # continuous at both knots
        self.assertAlmostEqual(acr_from_pcr(50.0 * (1 - 1e-9)), acr_from_pcr(50.0), places=6)
        self.assertAlmostEqual(acr_from_pcr(500.0 * (1 + 1e-9)), acr_from_pcr(500.0), places=6)
        # below 50 only the 0.3072 slope moves: ratio for halving = 2^-0.3072
        self.assertAlmostEqual(acr_from_pcr(20.0) / acr_from_pcr(40.0), 2.0 ** -0.3072, places=9)

    def test_rejects_non_positive(self):
        for bad in (0.0, -1.0, None):
            with self.assertRaises(ValueError):
                acr_from_pcr(bad)


if __name__ == "__main__":
    unittest.main()
