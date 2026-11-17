# SPDX-License-Identifier: Apache-2.0
# Every expected value is worked out by hand in the comments. Invented records
# only; no NHANES data are read by these tests. Stdlib only (no pandas).
import math
import unittest

from kfretransport.equations import egfr_ckdepi_2021, kfre4
from kfretransport.nhanes import (_num, jackknife, make_row, parse_mortality_line,
                                  policy_estimates, policy_prob)

NAN = float("nan")
BASE = {"SEQN": 41475.0, "RIDAGEYR": 70.0, "RIAGENDR": 1.0, "RIDEXPRG": NAN,
        "WTMEC2YR": 6000.0, "SDMVSTRA": 80.0, "SDMVPSU": 2.0, "LBXSCR": 2.0,
        "URDACT": 300.0, "DIQ010": 1.0, "BPQ020": 2.0}


def rec(**kw):
    r = dict(BASE)
    r.update(kw)
    return r


class TestMakeRow(unittest.TestCase):
    def test_valid(self):
        # man, 70, creatinine 2.0: CKD-EPI 2021 = 142 * (2/0.9)^-1.2 * 0.9938^70
        # = 142 * 0.38358 * 0.64706 = 35.2 -> G3b, kept
        row, why = make_row(rec(), 6, "2007-2008")
        self.assertIsNone(why)
        egfr = egfr_ckdepi_2021(2.0, 70.0, False)
        self.assertAlmostEqual(row["egfr"], egfr, places=12)
        self.assertAlmostEqual(row["egfr"], 35.2, delta=0.1)
        self.assertAlmostEqual(row["w"], 1000.0, places=9)        # 6000 / 6 cycles
        self.assertEqual((row["seqn"], row["stratum"], row["psu"]), (41475, 80, 2))
        self.assertTrue(row["dm"])
        self.assertFalse(row["htn"])
        self.assertTrue(row["male"])
        self.assertAlmostEqual(row["kfre5"], kfre4(70.0, True, egfr, 300.0, 5), places=12)
        self.assertAlmostEqual(row["kfre2"], kfre4(70.0, True, egfr, 300.0, 2), places=12)

    def test_exclusions(self):
        cases = [(rec(RIDAGEYR=19.0), "age_under_20"),
                 (rec(RIAGENDR=NAN), "sex_missing"),
                 (rec(RIAGENDR=2.0, RIDEXPRG=1.0), "pregnant"),
                 (rec(WTMEC2YR=5.397605e-79), "no_exam_weight"),   # SAS XPORT zero
                 (rec(LBXSCR=NAN), "no_creatinine"),
                 (rec(URDACT=None), "no_uacr"),
                 (rec(LBXSCR=0.8, RIDAGEYR=40.0), "egfr_not_15_to_59"),
                 (rec(LBXSCR=8.0), "egfr_not_15_to_59")]
        for r, expected in cases:
            self.assertEqual(make_row(r, 6), (None, expected))

    def test_num(self):
        self.assertIsNone(_num(NAN))
        self.assertIsNone(_num("x"))
        self.assertEqual(_num(5.397605e-79), 0.0)
        self.assertEqual(_num(2), 2.0)


def person(w, dm, htn, kfre5, age=70.0, male=True, egfr=40.0, uacr=100.0, kdeath=None):
    return {"w": w, "dm": dm, "htn": htn, "kfre5": kfre5, "age": age, "male": male,
            "egfr": egfr, "uacr": uacr, "kdeath": kdeath, "stratum": 1, "psu": 1}


class TestPolicyEstimates(unittest.TestCase):
    def setUp(self):
        self.rows = [person(1.0, True, True, 0.10, kdeath=False),
                     person(1.0, False, True, 0.10, kdeath=True),
                     person(2.0, False, False, 0.01)]

    def test_us_claims(self):
        # P1: p = 0.63, 0.25, 0.25; W = 4; tested weights 0.63, 0.25, 0.50 (T = 1.38)
        r = policy_estimates(self.rows, policy_prob("P1_us_claims"), 0.05)
        self.assertAlmostEqual(r["tested_share"], 0.345, places=12)       # 1.38 / 4
        self.assertAlmostEqual(r["high_share"], 0.5, places=12)           # 2 / 4
        self.assertAlmostEqual(r["high_scored_share"], 0.44, places=12)   # 0.88 / 2
        self.assertAlmostEqual(r["high_unscored_weighted"], 1.12, places=12)  # 0.37 + 0.75
        self.assertAlmostEqual(r["mean_kfre_all"], 0.055, places=12)      # 0.22 / 4
        # tested mean (0.063 + 0.025 + 0.005) / 1.38 = 0.0673913; / 0.055 = 1.2252964
        self.assertAlmostEqual(r["mean_kfre_tested"], 0.0673913, places=6)
        self.assertAlmostEqual(r["tested_all_ratio"], 1.2252964, places=6)
        # kidney deaths: only person 2 (weight 1), untested weight 0.75
        self.assertAlmostEqual(r["kdeath_unscored_share"], 0.75, places=12)

    def test_guideline_and_diabetes_only(self):
        r = policy_estimates(self.rows, policy_prob("P2_guideline"), 0.05)
        self.assertAlmostEqual(r["tested_share"], 0.5, places=12)
        self.assertAlmostEqual(r["high_scored_share"], 1.0, places=12)
        self.assertAlmostEqual(r["high_unscored_weighted"], 0.0, places=12)
        self.assertAlmostEqual(r["tested_all_ratio"], 0.1 / 0.055, places=12)  # 1.8181818
        r = policy_estimates(self.rows, policy_prob("P3_diabetes_only"), 0.05)
        self.assertAlmostEqual(r["high_scored_share"], 0.5, places=12)
        self.assertAlmostEqual(r["high_unscored_weighted"], 1.0, places=12)
        self.assertAlmostEqual(r["kdeath_unscored_share"], 1.0, places=12)

    def test_threshold_and_universal(self):
        r = policy_estimates(self.rows, policy_prob("P1_us_claims"), 0.005)
        self.assertAlmostEqual(r["high_scored_share"], 0.345, places=12)  # everyone high
        r = policy_estimates(self.rows, policy_prob("P0_universal"), 0.05)
        self.assertAlmostEqual(r["tested_all_ratio"], 1.0, places=12)
        self.assertIsNone(r["imputation_ratio"])       # nobody untested
        self.assertIsNone(r["selection_ratio"])
        q = policy_estimates(self.rows, policy_prob("P1q", 0.25), 0.05)
        self.assertAlmostEqual(q["tested_share"], 0.345, places=12)       # same as P1
        with self.assertRaises(ValueError):
            policy_prob("P9")

    def test_imputation_geometric_mean(self):
        # one cell (70, man, eGFR 40): tested donors (diabetes) uACR 1000 and 10,
        # equal weights -> imputed log uACR = (ln 1000 + ln 10)/2 = ln 100.
        # The untested person (no diabetes) really has uACR 10.
        k = lambda a: kfre4(70.0, True, 40.0, a, 5)  # noqa: E731
        rows = [person(1.0, True, False, k(1000.0), uacr=1000.0),
                person(1.0, True, False, k(10.0), uacr=10.0),
                person(1.0, False, False, k(10.0), uacr=10.0)]
        r = policy_estimates(rows, policy_prob("P3_diabetes_only"), 0.05)
        self.assertAlmostEqual(r["imputation_ratio"], k(100.0) / k(10.0), places=12)
        self.assertGreater(r["imputation_ratio"], 1.0)
        self.assertEqual((r["imputation_fallbacks"], r["imputation_no_donor"]), (0, 0))
        # uniform reference: all three are donors and all three are imputed with
        # mean log uACR (ln 1000 + 2 ln 10)/3 = ln 46.4159;
        # ratio = 3 k(46.4159) / (k(1000) + 2 k(10)), below 1 because the KFRE is
        # not linear in log uACR (Jensen). selection_ratio = P3 ratio / uniform ratio.
        m = math.exp((math.log(1000.0) + 2.0 * math.log(10.0)) / 3.0)
        self.assertAlmostEqual(m, 46.4158883, places=6)
        ref = 3.0 * k(m) / (k(1000.0) + 2.0 * k(10.0))
        self.assertAlmostEqual(r["imputation_ratio_uniform"], ref, places=12)
        self.assertLess(ref, 1.0)
        self.assertAlmostEqual(r["selection_ratio"], (k(100.0) / k(10.0)) / ref, places=12)

    def test_imputation_fallback_and_no_donor(self):
        # untested man aged 50 has no donor in his cell (donors are 70): falls
        # back to the eGFR band 30-44 (donor uACR 100 -> imputed 100).
        # untested person with eGFR 20 has no donor in band 15-29: skipped.
        k50 = kfre4(50.0, True, 40.0, 10.0, 5)
        rows = [person(1.0, True, False, kfre4(70.0, True, 40.0, 100.0, 5), uacr=100.0),
                person(1.0, False, False, k50, age=50.0, uacr=10.0),
                person(1.0, False, False, 0.5, egfr=20.0, uacr=10.0)]
        r = policy_estimates(rows, policy_prob("P3_diabetes_only"), 0.05)
        self.assertEqual((r["imputation_fallbacks"], r["imputation_no_donor"]), (1, 1))
        self.assertAlmostEqual(r["imputation_ratio"],
                               kfre4(50.0, True, 40.0, 100.0, 5) / k50, places=12)


def mean_stat(rows):
    def stat(w):
        s = sum(w)
        return {"mean": (sum(wi * r["x"] for wi, r in zip(w, rows)) / s) if s > 0 else None,
                "label": "not a number"}
    return stat


class TestJackknife(unittest.TestCase):
    def test_one_stratum_two_psus(self):
        # x = 1 and 3: mean 2; drop PSU 1 -> 3, drop PSU 2 -> 1;
        # var = (1/2) * (1 + 1) = 1 -> se 1 (= sd / sqrt(n) for two values)
        rows = [{"w": 1.0, "stratum": 1, "psu": 1, "x": 1.0},
                {"w": 1.0, "stratum": 1, "psu": 2, "x": 3.0}]
        theta, se, info = jackknife(rows, mean_stat(rows))
        self.assertAlmostEqual(theta["mean"], 2.0, places=12)
        self.assertAlmostEqual(se["mean"], 1.0, places=12)
        self.assertNotIn("label", se)
        self.assertEqual((info["strata"], info["psus"], info["df"]), (1, 2, 1))

    def test_single_psu_stratum_is_skipped(self):
        # A: PSUs x=1, x=3; B: one PSU x=2. mean 2. Drop A1: weights 0, 2, 1 ->
        # 8/3; drop A2: 2, 0, 1 -> 4/3. var = (1/2)(4/9 + 4/9) = 4/9 -> se 2/3
        rows = [{"w": 1.0, "stratum": 1, "psu": 1, "x": 1.0},
                {"w": 1.0, "stratum": 1, "psu": 2, "x": 3.0},
                {"w": 1.0, "stratum": 2, "psu": 1, "x": 2.0}]
        theta, se, info = jackknife(rows, mean_stat(rows))
        self.assertAlmostEqual(se["mean"], 2.0 / 3.0, places=12)
        self.assertEqual(info["strata_with_one_psu"], 1)
        self.assertEqual(info["df"], 1)     # 3 PSUs - 2 strata


class TestMortality(unittest.TestCase):
    def test_fixed_width(self):
        # columns: SEQN 1-6, ELIGSTAT 15, MORTSTAT 16, UCOD_LEADING 17-19
        line = "041475" + " " * 8 + "1" + "1" + "009" + "00" + " " * 21 + "100100"
        self.assertEqual(parse_mortality_line(line),
                         {"seqn": 41475, "eligstat": 1, "mortstat": 1, "ucod_leading": 9})
        line = "041476" + " " * 8 + "2" + " " + "   "
        self.assertEqual(parse_mortality_line(line),
                         {"seqn": 41476, "eligstat": 2, "mortstat": None, "ucod_leading": None})


if __name__ == "__main__":
    unittest.main()
