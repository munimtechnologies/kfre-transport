# SPDX-License-Identifier: Apache-2.0
# Tests for omop.py on a tiny, invented OMOP extract (7 invented people; no real
# data). Concept ids in the fixture are made up on purpose: omop.py must find
# everything through vocabulary codes. Expected values were worked out by hand;
# the working is in the comments.
import json
import os
import tempfile
import unittest

from kfretransport.equations import egfr_ckdepi_2021, kfre4
from kfretransport.omop import (age_on, build_cohort, main, parse_date, read_cdm,
                                read_table_text, site_summary, vocabulary_drift)

CONCEPT = """concept_id,concept_name,domain_id,vocabulary_id,concept_class_id,standard_concept,concept_code,valid_start_date,valid_end_date,invalid_reason
1001,Creatinine SerPl mass,Measurement,LOINC,Lab Test,S,2160-0,1970-01-01,2099-12-31,
1002,Creatinine SerPl moles,Measurement,LOINC,Lab Test,S,14682-9,1970-01-01,2099-12-31,
1003,Albumin/Creat Ur mass ratio,Measurement,LOINC,Lab Test,S,9318-7,1970-01-01,2099-12-31,
1004,Albumin/Creat Ur ratio,Measurement,LOINC,Lab Test,S,32294-1,1970-01-01,2099-12-31,
2001,milligram per deciliter,Unit,UCUM,Unit,S,mg/dL,1970-01-01,2099-12-31,
2002,micromole per liter,Unit,UCUM,Unit,S,umol/L,1970-01-01,2099-12-31,
2003,milligram per gram,Unit,UCUM,Unit,S,mg/g,1970-01-01,2099-12-31,
2004,milligram per millimole,Unit,UCUM,Unit,S,mg/mmol,1970-01-01,2099-12-31,
3001,Hemodialysis,Procedure,SNOMED,Procedure,S,302497006,1970-01-01,2099-12-31,
3002,Transplant of kidney,Procedure,SNOMED,Procedure,S,70536003,1970-01-01,2099-12-31,
3003,End-stage renal disease,Condition,SNOMED,Disorder,S,46177005,1970-01-01,2099-12-31,
3004,Dependence on renal dialysis,Observation,SNOMED,Clinical Finding,S,105502003,1970-01-01,2099-12-31,
3005,Invented home hemodialysis child,Procedure,SNOMED,Procedure,S,X-HOME-HD,1970-01-01,2099-12-31,
4001,End stage renal disease,Condition,ICD10CM,5-char billing code,,N18.6,1970-01-01,2099-12-31,
4002,Dependence on renal dialysis,Observation,ICD10CM,4-char billing code,,Z99.2,1970-01-01,2099-12-31,
"""

ANCESTOR = """ancestor_concept_id,descendant_concept_id,min_levels_of_separation,max_levels_of_separation
3001,3001,0,0
3001,3005,1,1
"""

PERSON = """person_id,gender_concept_id,year_of_birth,month_of_birth,day_of_birth
1,8507,1950,1,1
2,8532,1960,6,15
3,8507,1945,3,1
4,0,1950,1,1
5,8507,1970,1,1
6,8507,1950,1,1
7,8507,1950,1,1
"""

MEASUREMENT = """measurement_id,person_id,measurement_concept_id,measurement_date,value_as_number,unit_concept_id,unit_source_value
1,1,1001,2020-01-10,2.0,2001,mg/dL
2,1,1001,2020-02-10,2.0,2001,mg/dL
3,1,1002,2020-04-20,176.8,2002,umol/L
4,1,1003,2019-03-01,999,2003,mg/g
5,1,1004,2020-03-01,30,2004,mg/mmol
6,1,1003,2020-05-01,5,2003,mg/g
7,2,1001,2019-01-05,1.5,2001,mg/dL
8,2,1001,2019-03-01,1.0,2001,mg/dL
9,2,1001,2019-06-01,1.5,2001,mg/dL
10,2,1001,2019-10-01,1.4,2001,mg/dL
11,2,1001,2019-10-01,1.6,2001,mg/dL
12,3,1001,2020-01-01,2.0,2001,mg/dL
13,3,1001,2020-06-01,2.0,2001,mg/dL
14,5,1001,2020-01-01,2.0,2001,mg/dL
15,5,1001,2020-01-31,2.0,2001,mg/dL
16,5,1001,2020-02-15,1.0,2001,mg/dL
17,6,1001,2018-01-01,2.0,2001,mg/dL
18,6,1001,2018-06-01,2.0,2001,mg/dL
19,6,1003,2018-05-01,50,2003,mg/g
20,6,1003,2018-05-15,60,0,mg/g
21,6,1003,2018-05-20,70,0,xyz
22,7,1001,2023-01-01,2.0,2001,mg/dL
23,7,1001,2023-05-01,2.0,2001,mg/dL
24,7,1003,2023-04-01,10,2003,mg/g
25,4,1001,2020-01-01,2.0,2001,mg/dL
"""

CONDITION = """condition_occurrence_id,person_id,condition_concept_id,condition_start_date,condition_source_concept_id
1,2,0,2022-01-01,4002
2,7,3003,2023-06-01,4001
"""

PROCEDURE = """procedure_occurrence_id,person_id,procedure_concept_id,procedure_date,procedure_source_concept_id
1,1,3005,2021-06-01,0
"""

OBSERVATION = """observation_id,person_id,observation_concept_id,observation_date,observation_source_concept_id
1,3,3004,2019-12-01,0
"""

DEATH = """person_id,death_date
2,2023-01-01
6,2019-06-01
"""

PERIOD = """observation_period_id,person_id,observation_period_start_date,observation_period_end_date
1,1,2015-01-01,2025-12-31
2,2,2018-01-01,2024-12-31
3,3,2015-01-01,2025-12-31
4,5,2015-01-01,2025-12-31
5,6,2010-01-01,2025-12-31
6,7,2020-01-01,2024-01-01
"""

FILES = {"CONCEPT.csv": CONCEPT, "concept_ancestor.csv": ANCESTOR, "person.csv": PERSON,
         "measurement.csv": MEASUREMENT, "condition_occurrence.csv": CONDITION,
         "procedure_occurrence.csv": PROCEDURE, "observation.csv": OBSERVATION,
         "death.csv": DEATH, "observation_period.csv": PERIOD}


def write_fixture(folder, files=FILES):
    for name, text in files.items():
        with open(os.path.join(folder, name), "w", encoding="utf-8", newline="") as f:
            f.write(text)


class FixtureCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        write_fixture(self.tmp.name)
        self.cdm = read_cdm(self.tmp.name)
        self.cohort = build_cohort(self.cdm)
        self.by_id = {r["person_id"]: r for r in self.cohort["rows"]}

    def tearDown(self):
        self.tmp.cleanup()


class TestReading(unittest.TestCase):
    def test_dates(self):
        self.assertEqual(parse_date("2020-04-20").isoformat(), "2020-04-20")
        self.assertEqual(parse_date("2020-04-20 13:45:00").isoformat(), "2020-04-20")
        self.assertEqual(parse_date("20200420").isoformat(), "2020-04-20")
        self.assertIsNone(parse_date(""))
        with self.assertRaises(ValueError):
            parse_date("20/04/2020")

    def test_tsv_with_bom_and_quotes(self):
        # Athena files are tab-separated and may hold bare quote characters.
        rows = read_table_text(chr(0xFEFF) + 'CONCEPT_ID\tconcept_name\n7\tsays "hi"\n')
        self.assertEqual(rows, [{"concept_id": "7", "concept_name": 'says "hi"'}])

    def test_age_on(self):
        self.assertEqual(age_on(parse_date("1960-06-15"), parse_date("2019-06-14")), 58)
        self.assertEqual(age_on(parse_date("1960-06-15"), parse_date("2019-06-15")), 59)


class TestCohort(FixtureCase):
    def test_who_is_in_and_out(self):
        # 1, 2, 6, 7 qualify. 3: dialysis dependence recorded 2019-12-01, before its
        # index 2020-06-01. 4: gender 0. 5: low 2020-01-01 and 01-31 (30 days),
        # then eGFR about 91.7 on 02-15 ends the run, so no confirmation.
        self.assertEqual(sorted(self.by_id), [1, 2, 6, 7])
        self.assertEqual(self.cohort["excluded"], {"sex_not_male_or_female": 1,
                                                   "krt_on_or_before_index": 1,
                                                   "no_confirmed_egfr_10_to_60": 1})

    def test_person1_index_units_and_uacr(self):
        r = self.by_id[1]
        # run starts 2020-01-10; 02-10 is 31 days in; 04-20 (176.8 umol/L = 2.0
        # mg/dL) is 101 days in -> index. Age 70.
        self.assertEqual(r["index_date"], "2020-04-20")
        self.assertEqual(r["age"], 70)
        self.assertTrue(r["male"])
        # 142 * (2/0.9)^-1.2 (0.383579) * 0.9938^70 (0.647037) = 35.2430
        self.assertAlmostEqual(r["egfr"], 35.2430, delta=0.01)
        self.assertAlmostEqual(r["egfr"], egfr_ckdepi_2021(2.0, 70, female=False), places=9)
        # uACR: 2019-03-01 is 416 days before index (out of the 365-day window),
        # 2020-05-01 is after index; 2020-03-01: 30 mg/mmol * 8.84 = 265.2 mg/g.
        self.assertTrue(r["tested"])
        self.assertAlmostEqual(r["uacr_mg_g"], 265.2, places=9)

    def test_person1_outcome_found_through_descendant(self):
        r = self.by_id[1]
        # procedure 3005 is a descendant of Hemodialysis (3001) in CONCEPT_ANCESTOR.
        # 2020-04-20 -> 2021-04-20 = 365 days, -> 05-20 = 30, -> 06-01 = 12: 407.
        self.assertEqual((r["event"], r["time_days"], r["outcome_path"]), (1, 407, "standard"))
        self.assertEqual((r["y2"], r["y5"]), (1, 1))

    def test_person1_kfre(self):
        r = self.by_id[1]
        # L = -0.2201*(7.0-7.036) + 0.2467*(1-0.5642) - 0.5567*(35.2430/5-7.222)
        #     + 0.4510*(ln 265.2 - 5.137)
        #   = 0.0079236 + 0.1075119 + 0.0965320 + 0.4510*0.4434851 (0.2000118) = 0.4119793
        # exp(L) = 1.509803; 2 y: 1 - exp(1.509803 * ln 0.975) = 1 - exp(-0.0382249) = 0.037504
        self.assertAlmostEqual(r["kfre2"], 0.037504, delta=0.0002)
        self.assertAlmostEqual(r["kfre5"], kfre4(70, True, r["egfr"], 265.2, years=5), places=12)

    def test_person2_reset_sameday_mean_and_source_path(self):
        r = self.by_id[2]
        # 01-05 low (eGFR ~40.1), 03-01 eGFR ~65.3 resets, 06-01 low again (age 58),
        # 10-01 is 122 days later -> index; two same-day values 1.4 and 1.6 -> 1.5.
        self.assertEqual(r["index_date"], "2019-10-01")
        self.assertEqual(r["age"], 59)
        self.assertFalse(r["male"])
        self.assertAlmostEqual(r["egfr"], egfr_ckdepi_2021(1.5, 59, female=True), places=9)
        self.assertFalse(r["tested"])
        self.assertIsNone(r["kfre2"])
        # Z99.2 found only through condition_source_concept_id (standard id is 0).
        # 2019-10-01 -> 2022-01-01 = 92 + 366 + 365 = 823 days; death came later.
        self.assertEqual((r["event"], r["time_days"], r["outcome_path"]), (1, 823, "source"))
        self.assertEqual((r["y2"], r["y5"]), (0, 1))

    def test_person6_unit_text_fallback_and_death(self):
        r = self.by_id[6]
        # uACR 05-01 (50, mg/g concept), 05-15 (60, unit from text 'mg/g'),
        # 05-20 (70, unit 'xyz': dropped) -> latest usable is 60.
        self.assertEqual(r["index_date"], "2018-06-01")
        self.assertAlmostEqual(r["uacr_mg_g"], 60.0, places=12)
        # death 2019-06-01, 365 days after index: competing event, counts as y = 0.
        self.assertEqual((r["event"], r["time_days"]), (2, 365))
        self.assertEqual((r["y2"], r["y5"]), (0, 0))
        self.assertEqual(self.cohort["counts"], {"uacr_unit_from_text": 1,
                                                 "uacr_unknown_units": {"xyz": 1}})

    def test_person7_censored_and_esrd_diagnosis_is_not_krt(self):
        r = self.by_id[7]
        # An ESRD diagnosis (N18.6 -> 46177005) is not dialysis or transplant.
        # Follow-up ends 2024-01-01: 245 days after 2023-05-01.
        self.assertEqual((r["event"], r["time_days"]), (0, 245))
        self.assertEqual((r["y2"], r["y5"]), (None, None))

    def test_concept_notes(self):
        notes = self.cohort["notes"]
        for missing in ("missing LOINC 38483-4", "missing SNOMED 71192002",
                        "missing CPT4 1012740", "missing ICD10CM Z94.0"):
            self.assertIn(missing, notes)
        self.assertFalse([n for n in notes if "2160-0" in n or "Z99.2" in n])

    def test_site_summary(self):
        s2 = site_summary(self.cohort["rows"], 2)
        # tested: 1, 6, 7 of 4; scored at 2 y: 1 (y=1) and 6 (y=0); 7 has no outcome.
        self.assertEqual((s2["n_cohort"], s2["n_tested"], s2["n_scored"], s2["events_scored"]),
                         (4, 3, 2, 1))
        self.assertAlmostEqual(s2["testing_rate"], 0.75, places=12)
        exp = self.by_id[1]["kfre2"] + self.by_id[6]["kfre2"]
        self.assertAlmostEqual(s2["expected_scored"], exp, places=12)
        self.assertAlmostEqual(s2["oe"], 1.0 / exp, places=12)
        self.assertEqual((s2["n_untested_with_outcome"], s2["events_untested"]), (1, 0))
        s5 = site_summary(self.cohort["rows"], 5)
        self.assertEqual((s5["n_untested_with_outcome"], s5["events_untested"]), (1, 1))


class TestVocabularyDrift(unittest.TestCase):
    def run_drift(self, files):
        with tempfile.TemporaryDirectory() as d:
            write_fixture(d, files)
            cdm = read_cdm(d)
            return vocabulary_drift(cdm), build_cohort(cdm)

    def test_before_refresh(self):
        drift, _ = self.run_drift(FILES)
        w = {x["code"]: x for x in drift["watch"]}
        self.assertTrue(w["SNOMED 46177005"]["standard"])
        self.assertEqual(w["SNOMED 46177005"]["condition_records"], 1)
        self.assertFalse(w["SNOMED 433146000"]["in_site_vocabulary"])
        self.assertEqual(drift["icd10cm_n18_6_maps_to"], {"SNOMED 46177005": 1})

    def test_after_refresh_outcome_unchanged(self):
        # Simulated refresh: ESRD becomes non-standard ('U' = upgraded/replaced)
        # and N18.6 now lands on CKD stage 5 (433146000).
        files = dict(FILES)
        files["CONCEPT.csv"] = CONCEPT.replace(
            "3003,End-stage renal disease,Condition,SNOMED,Disorder,S,46177005,1970-01-01,2099-12-31,",
            "3003,End-stage renal disease,Condition,SNOMED,Disorder,,46177005,1970-01-01,2026-06-01,U"
        ) + "3006,Chronic kidney disease stage 5,Condition,SNOMED,Disorder,S,433146000,1970-01-01,2099-12-31,\n"
        files["condition_occurrence.csv"] = CONDITION.replace("2,7,3003,", "2,7,3006,")
        drift, cohort = self.run_drift(files)
        w = {x["code"]: x for x in drift["watch"]}
        self.assertFalse(w["SNOMED 46177005"]["standard"])
        self.assertEqual(w["SNOMED 46177005"]["invalid_reason"], "U")
        self.assertEqual(w["SNOMED 46177005"]["condition_records"], 0)
        self.assertEqual(w["SNOMED 433146000"]["condition_records"], 1)
        self.assertEqual(drift["icd10cm_n18_6_maps_to"], {"SNOMED 433146000": 1})
        # The KRT-based outcome does not use ESRD, so the cohort is identical.
        with tempfile.TemporaryDirectory() as d:
            write_fixture(d)
            before = build_cohort(read_cdm(d))
        self.assertEqual(cohort["rows"], before["rows"])


class TestCli(unittest.TestCase):
    def test_writes_summary_without_person_rows(self):
        with tempfile.TemporaryDirectory() as d:
            write_fixture(d)
            out = os.path.join(d, "out.json")
            main(["--cdm-dir", d, "--out", out])
            with open(out, encoding="utf-8") as f:
                text = f.read()
        data = json.loads(text)
        self.assertNotIn("person_id", text)
        self.assertEqual([s["years"] for s in data["summary"]], [2, 5])
        self.assertEqual(data["excluded"]["krt_on_or_before_index"], 1)


if __name__ == "__main__":
    unittest.main()
