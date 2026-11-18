# SPDX-License-Identifier: Apache-2.0
# Tests for nested.py (0.2.5). Exact identities first, then hand values (working in the
# comments), then the cohort path on the invented 7-person OMOP extract of test_omop.py.
import json
import math
import os
import tempfile
import unittest
from statistics import NormalDist

from kfretransport.equations import egfr_ckdepi_2021, kfre4
from kfretransport.nested import (band_factor, band_factors, band_score, bands_from_rows,
                                  bands_of, consistency, fit_within, main,
                                  poisson_band_factors, public_report, site_bands,
                                  split_bands, within_information, within_se)
from kfretransport.omop import read_cdm
from kfretransport.selection import tested_factor
from kfretransport.test_omop import FILES, MEASUREMENT, write_fixture

N = NormalDist()

# Person 7's only uACR moved from 2023-04-01 (30 days before index 2023-05-01) to
# 2022-03-01: 365 days to 2023-03-01, + 31 + 30 = 426 days before index.
MOVED = dict(FILES)
MOVED["measurement.csv"] = MEASUREMENT.replace("24,7,1003,2023-04-01,10,2003,mg/g",
                                               "24,7,1003,2022-03-01,10,2003,mg/g")


class TestCurve(unittest.TestCase):
    def test_no_selection_gives_one(self):
        for lo, hi in ((0.0, 0.2), (0.1, 0.35), (0.6, 0.8), (0.9, 1.0)):
            self.assertAlmostEqual(band_factor(lo, hi, 0.0), 1.0, places=12)

    def test_bands_add_up_to_the_tested_curve(self):
        # sum_j (P_j - P_(j-1)) * factor_j telescopes to Phi(z_PJ + theta)
        cum, th = [0.12, 0.2, 0.31], 0.7
        f = band_factors(cum, th)
        widths = [hi - lo for lo, hi in bands_of(cum)]
        total = sum(w * x for w, x in zip(widths, f))
        self.assertAlmostEqual(total, N.cdf(N.inv_cdf(0.31) + 0.7), places=12)
        self.assertAlmostEqual(total, 0.31 * tested_factor(0.31, 0.7), places=12)

    def test_first_band_is_tested_factor(self):
        for p, th in ((0.05, 0.69), (0.3, 0.3), (0.7, -0.4)):
            self.assertAlmostEqual(band_factor(0.0, p, th), tested_factor(p, th), places=12)

    def test_whole_range_is_one(self):
        self.assertAlmostEqual(band_factor(0.0, 1.0, 0.9), 1.0, places=12)

    def test_order_follows_sign_of_theta(self):
        cum = [0.1, 0.3, 0.6, 0.9]
        down = band_factors(cum, 0.5)
        up = band_factors(cum, -0.5)
        for a, b in zip(down, down[1:]):
            self.assertGreater(a, b)
        for a, b in zip(up, up[1:]):
            self.assertLess(a, b)

    def test_upper_tail_band_is_finite(self):
        # both edges far in the upper tail after the shift; the erfc form keeps it positive
        self.assertGreater(band_factor(0.97, 0.99, 3.0), 0.0)

    def test_bad_bands(self):
        with self.assertRaises(ValueError):
            band_factor(0.3, 0.3, 0.1)
        with self.assertRaises(ValueError):
            bands_of([0.2, 0.1])
        with self.assertRaises(ValueError):
            bands_of([0.2, 1.2])


class TestScore(unittest.TestCase):
    def test_is_derivative_of_log_band(self):
        h = 1e-5
        for lo, hi, th in ((0.1, 0.3, 0.4), (0.0, 0.2, 0.0), (0.6, 0.95, -0.3)):
            fd = (math.log(band_factor(lo, hi, th + h)) - math.log(band_factor(lo, hi, th - h))) / (2 * h)
            self.assertAlmostEqual(band_score(lo, hi, th), fd, places=6)

    def test_first_band_is_mills(self):
        for p, th in ((0.1, 0.0), (0.4, 0.5)):
            a = N.inv_cdf(p) + th
            self.assertAlmostEqual(band_score(0.0, p, th), N.pdf(a) / N.cdf(a), places=10)

    def test_hand_values(self):
        # phi(z_0.15) = phi(1.0364334) = 0.233160; phi(z_0.30) = phi(0.5244005) = 0.347693;
        # phi(z_0.40) = phi(0.2533471) = 0.386342 (all by hand)
        self.assertAlmostEqual(band_score(0.0, 0.15), 0.233160 / 0.15, places=3)          # 1.55440
        self.assertAlmostEqual(band_score(0.15, 0.30), (0.347693 - 0.233160) / 0.15, places=3)  # 0.76355
        self.assertAlmostEqual(band_score(0.30, 0.40), (0.386342 - 0.347693) / 0.10, places=3)  # 0.38649


class TestInformation(unittest.TestCase):
    def test_two_equal_bands(self):
        # weights 25 and 25: centred SS = 25 * 2 * ((m1 - m2) / 2)^2 = 12.5 (m1 - m2)^2
        site = {"cum_rates": [0.2, 0.5], "variances": [0.04, 0.04]}
        m1, m2 = band_score(0.0, 0.2), band_score(0.2, 0.5)
        self.assertAlmostEqual(within_information([site]), 12.5 * (m1 - m2) ** 2, places=10)

    def test_illustration_by_hand(self):
        # P = 0.15, 0.30, 0.40 with 90, 60 and 30 events (variance 1/events):
        # scores 1.55440, 0.76355, 0.38649; weighted mean 197.3039 / 180 = 1.096133;
        # info = 90 * 0.458267^2 + 60 * 0.332580^2 + 30 * 0.709643^2 = 40.645; se = 0.15686
        site = {"cum_rates": [0.15, 0.30, 0.40], "variances": [1 / 90.0, 1 / 60.0, 1 / 30.0]}
        self.assertAlmostEqual(within_information([site]), 40.645, delta=0.05)
        self.assertAlmostEqual(within_se([site]), 0.15686, delta=0.001)

    def test_sites_add(self):
        a = {"cum_rates": [0.2, 0.5], "variances": [0.04, 0.04]}
        b = {"cum_rates": [0.1, 0.3, 0.6], "variances": [0.02, 0.03, 0.05]}
        self.assertAlmostEqual(within_information([a, b]),
                               within_information([a]) + within_information([b]), places=10)


class TestFit(unittest.TestCase):
    def test_recovers_theta_one_site(self):
        cum = [0.1, 0.25, 0.4]
        y = [0.2 + math.log(f) for f in band_factors(cum, 0.5)]
        fit = fit_within([{"cum_rates": cum, "log_oe": y, "variances": [0.01, 0.02, 0.03]}])
        self.assertAlmostEqual(fit["theta"], 0.5, places=9)
        self.assertAlmostEqual(fit["beta0"][0], 0.2, places=9)
        self.assertLess(fit["q_min"], 1e-12)
        self.assertLess(fit["ci95_profile"][0], 0.5)
        self.assertGreater(fit["ci95_profile"][1], 0.5)
        self.assertEqual((fit["df"], fit["n_sites"], fit["n_bands"]), (1, 1, 3))
        self.assertFalse(fit["theta_at_grid_edge"])

    def test_recovers_theta_two_sites_own_intercepts(self):
        sites = []
        for b0, cum in ((0.3, [0.05, 0.15, 0.3]), (-0.2, [0.3, 0.5, 0.7])):
            y = [b0 + math.log(f) for f in band_factors(cum, 0.8)]
            sites.append({"cum_rates": cum, "log_oe": y, "variances": [0.02, 0.02, 0.04]})
        fit = fit_within(sites)
        self.assertAlmostEqual(fit["theta"], 0.8, places=9)
        self.assertAlmostEqual(fit["beta0"][0], 0.3, places=9)
        self.assertAlmostEqual(fit["beta0"][1], -0.2, places=9)
        self.assertAlmostEqual(fit["oe_all"][1], math.exp(-0.2), places=9)
        self.assertEqual(fit["df"], 3)

    def test_bad_input(self):
        with self.assertRaises(ValueError):
            fit_within([{"cum_rates": [0.3], "log_oe": [0.1], "variances": [0.01]}])
        with self.assertRaises(ValueError):
            fit_within([{"cum_rates": [0.3, 0.2], "log_oe": [0.1, 0.0], "variances": [0.01, 0.01]}])
        with self.assertRaises(ValueError):
            fit_within([{"cum_rates": [0.2, 0.3], "log_oe": [0.1], "variances": [0.01, 0.01]}])
        with self.assertRaises(ValueError):
            fit_within([{"cum_rates": [0.2, 0.3], "log_oe": [0.1, 0.0], "variances": [0.01, 0.0]}])
        with self.assertRaises(ValueError):
            fit_within([])

    def test_consistency(self):
        # se_w = 0.392 / 3.919928 = 0.1, se_b = 0.784 / 3.919928 = 0.2 (to 6 places)
        c = consistency(0.5, [0.304, 0.696], 0.2, [-0.192, 0.592])
        self.assertAlmostEqual(c["se"], math.sqrt(0.1 ** 2 + 0.2 ** 2), places=4)
        self.assertAlmostEqual(c["z"], 0.3 / math.sqrt(0.05), places=4)    # 1.3416
        self.assertTrue(c["consistent"])
        self.assertFalse(consistency(1.0, [0.9, 1.1], 0.0, [-0.1, 0.1])["consistent"])


class TestPoisson(unittest.TestCase):
    def test_no_link_to_risk_gives_one(self):
        r = poisson_band_factors([30, 365, 730], 0.003, 0.0, 0.7)
        for f in r["factors"]:
            self.assertAlmostEqual(f, 1.0, places=9)
        for p, w in zip(r["cum_rates"], [30, 365, 730]):
            self.assertAlmostEqual(p, 1.0 - math.exp(-0.003 * w), places=9)
        for f in poisson_band_factors([30, 365], 0.003, 1.0, 0.0)["factors"]:
            self.assertAlmostEqual(f, 1.0, places=9)

    def test_short_window_limit(self):
        # E[exp(gU - g^2/2) * lambda] / E[lambda] = exp(g kappa) = exp(0.3) = 1.349859
        r = poisson_band_factors([0.0001], 1.0, 0.6, 0.5)
        self.assertAlmostEqual(r["factors"][0], math.exp(0.3), delta=1e-3)

    def test_recent_band_is_riskier(self):
        r = poisson_band_factors([30, 365, 3650], 0.01, 1.0, 0.5)
        f = r["factors"]
        self.assertGreater(f[0], f[1])
        self.assertGreater(f[1], f[2])
        self.assertLess(r["cum_rates"][-1], 1.0)

    def test_bad_windows(self):
        with self.assertRaises(ValueError):
            poisson_band_factors([365, 90], 0.01, 1.0, 0.5)


def _row(pid, tested, time=0, event=0, p=None):
    return {"person_id": pid, "tested": tested, "time_days": time, "event": event,
            "kfre5": p, "kfre2": p}


class TestSplitAndFit(unittest.TestCase):
    def rows(self):
        # windows 30 and 365. Band 1 = ids 1-4 (tested in both), band 2 = ids 5-8
        # (tested only at 365), ids 9-10 never tested.
        b1 = [(100, 1), (200, 0), (300, 1), (2000, 0)]
        b2 = [(150, 1), (400, 0), (500, 1), (1900, 0)]
        w30, w365 = [], []
        for i, (t, e) in enumerate(b1):
            w30.append(_row(i + 1, True, t, e, 0.25))
            w365.append(_row(i + 1, True, t, e, 0.25))
        for i, (t, e) in enumerate(b2):
            w30.append(_row(i + 5, False, t, e))
            w365.append(_row(i + 5, True, t, e, 0.5))
        for pid in (9, 10):
            w30.append(_row(pid, False, 500, 0))
            w365.append(_row(pid, False, 500, 0))
        return [w30, w365]

    def test_split(self):
        s = split_bands(self.rows())
        self.assertEqual(s["n_cohort"], 10)
        self.assertEqual(s["cum_tested"], [4, 8])
        self.assertEqual([r["person_id"] for r in s["bands"][1]], [5, 6, 7, 8])

    def test_bands_and_fit(self):
        # Band 1, 5 y (1826 d): t=100 1 of 4 -> CIF .25, S .75; 200 censored; t=300 1 of 2
        # -> CIF .25 + .75/2 = .625; 2000 is past the horizon. O/E = .625 / .25 = 2.5.
        # Band 2 the same shape with mean prediction .5: O/E = 1.25. Rates .4 and .8.
        r = bands_from_rows(self.rows(), [30, 365], 5, min_events=2)
        self.assertEqual(r["cum_rates"], [0.4, 0.8])
        self.assertAlmostEqual(r["bands"][0]["oe"], 2.5, places=12)
        self.assertAlmostEqual(r["bands"][1]["oe"], 1.25, places=12)
        self.assertEqual((r["bands"][1]["window_from_days"], r["bands"][1]["window_to_days"]),
                         (30, 365))
        self.assertTrue(r["usable"])
        # band-1 / band-2 factor ratio 2 needs theta near 0.53 (by hand: 1.91 at 0.5, 2.19 at 0.6)
        self.assertGreater(r["fit"]["theta"], 0.4)
        self.assertLess(r["fit"]["theta"], 0.7)
        self.assertEqual(r["fit"]["df"], 0)

    def test_not_nested(self):
        rows = self.rows()
        rows[1][0] = _row(1, False, 100, 1)
        with self.assertRaises(ValueError):
            split_bands(rows)

    def test_different_cohorts(self):
        rows = self.rows()
        rows[1] = rows[1][:-1]
        with self.assertRaises(ValueError):
            split_bands(rows)


class TestCohortPath(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        write_fixture(self.tmp.name, MOVED)
        self.cdm = read_cdm(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_three_bands_no_merging(self):
        # Cohort 1, 2, 6, 7. Latest uACR before index: person 6 17 days (05-15 -> 06-01),
        # person 1 50 days (03-01 -> 04-20), person 7 426 days; person 2 none.
        # Windows 30, 60, 730 -> bands [6], [1], [7]; cumulative 1, 2, 3 of 4.
        y2, y5 = site_bands(self.cdm, [730, 30, 60], min_events=0)
        self.assertEqual(y5["windows"], [30, 60, 730])
        self.assertEqual(y5["cum_tested"], [1, 2, 3])
        self.assertEqual(y5["cum_rates"], [0.25, 0.5, 0.75])
        b = y5["bands"]
        # band 1: person 6 dies on day 365 (competing): no KRT, O/E 0
        self.assertEqual((b[0]["n"], b[0]["events"], b[0]["oe"]), (1, 0, 0.0))
        # band 2: person 1, dialysis day 407: CIF 1, O/E = 1 / KFRE5 (uACR 265.2 mg/g)
        k1 = kfre4(70, True, egfr_ckdepi_2021(2.0, 70, female=False), 265.2, years=5)
        self.assertAlmostEqual(b[1]["oe"], 1.0 / k1, places=9)
        # band 3: person 7, scored with the 426-day-old 10 mg/g; censored day 245
        k7 = kfre4(73, True, egfr_ckdepi_2021(2.0, 73, female=False), 10.0, years=5)
        self.assertAlmostEqual(b[2]["expected_risk"], k7, places=12)
        self.assertEqual(b[2]["events"], 0)
        self.assertFalse(y5["usable"])
        self.assertIsNone(y5["fit"])
        self.assertEqual(y2["years"], 2)

    def test_merging(self):
        # min 1 event: band [6] has none -> drop 30; bands [1, 6] (1 event) and [7] (none)
        # -> last band merges back: drop 60; one window left, not usable.
        r = site_bands(self.cdm, [30, 60, 730], years=(5,), min_events=1)[0]
        self.assertEqual(r["windows"], [730])
        self.assertEqual(r["windows_dropped"], [30, 60])
        self.assertEqual((r["bands"][0]["n"], r["bands"][0]["events"]), (3, 1))
        self.assertFalse(r["usable"])

    def test_public_report_masks(self):
        r = site_bands(self.cdm, [30, 60, 730], years=(5,), min_events=0)[0]
        p = public_report(r)
        self.assertTrue(all(b["masked"] for b in p["bands"]))
        self.assertIsNone(p["bands"][1]["oe"])
        self.assertIsNone(p["n_cohort"])
        self.assertIsNone(p["cum_rates"])
        # the original is untouched
        self.assertEqual(r["bands"][1]["n"], 1)

    def test_cli(self):
        out = os.path.join(self.tmp.name, "nested.json")
        main(["--cdm-dir", self.tmp.name, "--out", out, "--windows", "30", "60", "730",
              "--min-events", "0"])
        with open(out, encoding="utf-8") as f:
            text = f.read()
        data = json.loads(text)
        self.assertNotIn("person_id", text)
        self.assertEqual([h["years"] for h in data["horizons"]], [2, 5])
        self.assertEqual(data["horizons"][1]["windows"], [30, 60, 730])
        self.assertTrue(all(b["masked"] for b in data["horizons"][1]["bands"]))

    def test_needs_two_windows(self):
        with self.assertRaises(ValueError):
            site_bands(self.cdm, [365])


if __name__ == "__main__":
    unittest.main()
