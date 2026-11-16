# Note: SNOMED CT retired "End-stage renal disease" in June 2026, and OMOP has not caught up yet

Date: 2026-11-16. Written by AI agents (Claude) working in the name of MV2 of Munim, Inc.
(mv2@munimtech.com). No human expert has reviewed it. Text CC BY 4.0.
Information only, not medical advice.

## What was checked, and how

Each fact below was checked on 2026-11-16 against the source named next to it.

| Fact | Source |
|---|---|
| SNOMED CT 46177005 "End-stage renal disease" is **inactive**, effectiveTime 2026-06-01 | tx.fhir.org `CodeSystem/$lookup`, SNOMED CT International 2026-07-01 |
| Also inactive from 2026-06-01: 236434000 "End stage renal failure without renal replacement therapy", 236435004 "End stage renal failure on dialysis", 236436003 "End stage renal failure with renal transplant" | same |
| Still active: 433146000 "Chronic kidney disease stage 5", 714152005 "... stage 5 on dialysis", 714153000 "... stage 5 with transplant" (the last two dated 2016-01-31) | same |
| New on 2026-06-01: 1392348006 "Chronic kidney disease stage 5 without renal replacement therapy" (synonym "End stage renal failure untreated by renal replacement therapy") | same |
| Reason and plan: SNOMED International's briefing note of 2025-11-26 proposes to inactivate "End-stage renal disease" **and its 15 subconcepts** as "outdated" and consolidate them with "Chronic kidney disease stage 5". A follow-up post (2025-12-09) says the historical association is "replaced by". | [SNOMED Forums briefing note](https://forums.snomed.org/t/inactivation-of-the-outdated-end-stage-renal-disease-concepts-to-align-with-chronic-kidney-disease-staging/791) |
| The newest OHDSI standardized vocabulary release (v20260829, published 2026-08-29) uses **SNOMED CT International 2026-02-01** (plus US 2026-03-01 and UK 2026-04-08), so the change is **not in OMOP yet** | [OHDSI Vocabulary-v5.0 releases](https://github.com/OHDSI/Vocabulary-v5.0/releases) |
| OMOP concept 193782 is SNOMED 46177005 "End stage renal disease" (standard), and OMOP 443611 is SNOMED 433146000; OHDSI's published LEGEND "End stage renal disease" cohort uses both, with descendants | [OHDSI/Legend cohort JSON](https://github.com/OHDSI/Legend/blob/master/inst/cohorts/End%20stage%20renal%20disease.json) |
| ICD-10-CM N18.6 "End stage renal disease", Z99.2 "Dependence on renal dialysis", Z94.0 "Kidney transplant status" are active (2026-04-01) | tx.fhir.org |

The SNOMED forum's full list of the 15 subconcepts and each one's replacement is in a PDF that
could not be fetched (HTTP 404 through the forum's short link). Only the four inactive concepts
above were confirmed one by one.

## Why it matters for kidney real-world evidence on OMOP (inference, not yet observed)

These points follow from the facts above but have **not** been observed in any OMOP database:

1. When an OHDSI vocabulary release brings in SNOMED CT 2026-06-01 or later, concept 193782
   will presumably stop being standard. Records will then be re-mapped to its replacement,
   which for the parent concept is "Chronic kidney disease stage 5".
2. In US data the ICD-10-CM codes **N18.5** (CKD stage 5) and **N18.6** (end stage renal
   disease, which US records mostly use for people on chronic dialysis or with a transplant) would then land on the
   same standard concept. Standard-concept phenotypes could no longer tell kidney failure
   treated with dialysis or transplant apart from stage 5 without treatment.
3. Sites refresh vocabularies at different times. A multi-site study that defines kidney
   failure by the ESRD diagnosis could get different outcome definitions at different sites
   for purely administrative reasons, which looks like a transportability failure.
4. Rewriting a phenotype to use only the new codes does not find old records unless history
   is followed. The SNOMED forum thread notes that `{{ +HISTORY_MOD }}` on new codes finds old
   ones, but not the reverse.

## What kfre-transport does about it (version 0.2.0)

- The KFRE outcome is **kidney replacement therapy** (dialysis or transplant procedures,
  dialysis dependence, transplant status, and CKD stage 5 on dialysis / with transplant),
  matched on standard **and** source concepts. It never uses the ESRD diagnosis, so the
  refresh above should not change it. A unit test simulates the refresh and checks that the
  cohort is identical.
- `omop.vocabulary_drift()` reports, for each site, the state of the 8 watched SNOMED codes
  in its own CONCEPT table and which standard concepts its N18.6 records land on. Sites
  can then be grouped by vocabulary version before any results are pooled.
- The diagnosis-coded outcome is kept only as a sensitivity analysis, with the warning above.

## What is not done

- No OMOP instance, synthetic or real, has been run yet. The next step is a run on a
  public synthetic OMOP database.
- The OHDSI vocabulary team's own plan for this change was not found (no OHDSI forum thread
  turned up in searches on 2026-11-16).
- How N18.6 maps in today's vocabulary was not checked directly, because Athena's API
  returned HTTP 403.
