# SPDX-License-Identifier: Apache-2.0
"""Concept sets for the KFRE transportability study, written as vocabulary CODES
(LOINC, SNOMED, ICD10CM, UCUM), never as site-specific or release-specific OMOP
concept_ids. omop.py resolves each code through the site's own CONCEPT table,
and reports every code it cannot find, or finds non-standard or invalid.

Each code below was checked on 2026-11-16 against the HL7 FHIR terminology
server tx.fhir.org (LOINC 2.83; SNOMED CT International 2026-07-01;
ICD-10-CM 2026-04-01). The display names are the ones it returned. Check them
again before any real use.

Research code only. Information only, not medical advice, not a medical device.
"""

# (vocabulary_id, concept_code, include_descendants, display)
SERUM_CREATININE = [
    ("LOINC", "2160-0", False, "Creatinine [Mass/volume] in Serum or Plasma"),
    ("LOINC", "38483-4", False, "Creatinine [Mass/volume] in Blood"),
    ("LOINC", "14682-9", False, "Creatinine [Moles/volume] in Serum or Plasma"),
    ("LOINC", "59826-8", False, "Creatinine [Moles/volume] in Blood"),
]

# Urine albumin-to-creatinine ratio. 9318-7 is usually reported in ug/mg,
# which is numerically equal to mg/g; 32294-1 is usually reported in mg/mmol.
UACR = [
    ("LOINC", "9318-7", False, "Albumin/Creatinine [Mass Ratio] in Urine"),
    ("LOINC", "14959-1", False, "Microalbumin/Creatinine [Mass Ratio] in Urine"),
    ("LOINC", "32294-1", False, "Albumin/Creatinine [Ratio] in Urine"),
]

# Urine protein-to-creatinine ratio: used only by the pre-declared sensitivity
# analysis that converts PCR to an estimated ACR (albuminuria.py, 0.2.1).
UPCR = [("LOINC", "2890-2", False, "Protein/Creatinine [Mass Ratio] in Urine")]
# Dipstick protein: kept for reference; no conversion is implemented.
URINE_PROTEIN_DIPSTICK = [("LOINC", "5804-0", False, "Protein [Mass/volume] in Urine by Test strip")]

# Kidney failure, as the KFRE defines it: start of kidney replacement therapy
# (dialysis or transplant). Matched on STANDARD concepts...
KRT_STANDARD = [
    ("SNOMED", "302497006", True, "Hemodialysis"),
    ("SNOMED", "71192002", True, "Peritoneal dialysis"),
    ("SNOMED", "70536003", True, "Transplant of kidney"),
    ("SNOMED", "105502003", True, "Dependence on renal dialysis"),
    ("SNOMED", "737295003", True, "Transplanted kidney present"),
    ("SNOMED", "714152005", True, "Chronic kidney disease stage 5 on dialysis"),
    ("SNOMED", "714153000", True, "Chronic kidney disease stage 5 with transplant"),
    # CPT4 classification header used by OHDSI LEGEND's published dialysis
    # concept set (US sites code dialysis with CPT4, which stays standard).
    ("CPT4", "1012740", True, "Dialysis Services and Procedures"),
]
# ...and, separately, on SOURCE concepts, so that a record whose mapping to a
# standard concept changed or failed is still found.
KRT_SOURCE = [
    ("ICD10CM", "Z99.2", False, "Dependence on renal dialysis"),
    ("ICD10CM", "Z94.0", False, "Kidney transplant status"),
]

# A diagnosis-coded kidney-failure outcome, used ONLY as a sensitivity analysis.
# It is fragile: see VOCABULARY_WATCH.
DIAGNOSIS_KIDNEY_FAILURE_STANDARD = [
    ("SNOMED", "46177005", True, "End-stage renal disease"),
    ("SNOMED", "714152005", True, "Chronic kidney disease stage 5 on dialysis"),
    ("SNOMED", "714153000", True, "Chronic kidney disease stage 5 with transplant"),
]
DIAGNOSIS_KIDNEY_FAILURE_SOURCE = [
    ("ICD10CM", "N18.6", False, "End stage renal disease"),
]

# SNOMED CT International 2026-06-01 inactivated 'End-stage renal disease'
# (46177005) and its subconcepts as "outdated", with 'replaced by' links towards
# 'Chronic kidney disease stage 5' (433146000) (SNOMED International briefing
# note, 2025-11-26). It also added 1392348006. The OHDSI vocabulary release of
# 2026-08-29 still uses SNOMED International 2026-02-01, so OMOP sites will see
# the change only after a later vocabulary refresh. omop.vocabulary_drift reports
# the state of each code below in a site's CONCEPT table.
VOCABULARY_WATCH = [
    ("SNOMED", "46177005", "End-stage renal disease", "inactive from 2026-06-01"),
    ("SNOMED", "236434000", "End stage renal failure without renal replacement therapy", "inactive from 2026-06-01"),
    ("SNOMED", "236435004", "End stage renal failure on dialysis", "inactive from 2026-06-01"),
    ("SNOMED", "236436003", "End stage renal failure with renal transplant", "inactive from 2026-06-01"),
    ("SNOMED", "433146000", "Chronic kidney disease stage 5", "active"),
    ("SNOMED", "714152005", "Chronic kidney disease stage 5 on dialysis", "active"),
    ("SNOMED", "714153000", "Chronic kidney disease stage 5 with transplant", "active"),
    ("SNOMED", "1392348006", "Chronic kidney disease stage 5 without renal replacement therapy", "new 2026-06-01"),
]

# Units, keyed by the UCUM code (OMOP stores UCUM codes as the concept_code of
# the unit concept) or by the raw unit text, compared lower-case without spaces.
# Factor converts the value to mg/dL (creatinine) or mg/g (uACR, uPCR).
CREATININE_UNITS = {
    "mg/dl": 1.0,
    "umol/l": 1.0 / 88.4,
    "mmol/l": 1000.0 / 88.4,
}
UACR_UNITS = {
    "mg/g": 1.0,
    "mg/g{creat}": 1.0,
    "ug/mg": 1.0,
    "ug/mg{creat}": 1.0,
    "mg/mmol": 8.84,
    "mg/mmol{creat}": 8.84,
    "g/mol": 8.84,
    "g/mol{creat}": 8.84,
    "mg/mg": 1000.0,
    "mg/mg{creat}": 1000.0,
}
# uPCR is often reported in mg/mg (for example 0.2) or g/g in the US.
UPCR_UNITS = dict(UACR_UNITS)
UPCR_UNITS.update({"g/g": 1000.0, "g/g{creat}": 1000.0})

# Plausible ranges after conversion; values outside are dropped and counted.
CREATININE_RANGE_MG_DL = (0.2, 20.0)
UACR_RANGE_MG_G = (0.0, 50000.0)   # lower bound exclusive (the KFRE takes ln uACR)
UPCR_RANGE_MG_G = (0.0, 100000.0)  # lower bound exclusive (the conversion takes logs)

# OMOP's fixed gender concepts.
GENDER_MALE = 8507
GENDER_FEMALE = 8532
