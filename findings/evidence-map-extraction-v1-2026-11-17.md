# Evidence map v1, extraction and result: do published real-world KFRE validations already show the testing effect? (2026-11-17)

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it. Text CC BY 4.0. Information only, not medical advice, not a medical device.

This note carries out the pre-registration
[evidence-map-prereg-v1-2026-11-17.md](evidence-map-prereg-v1-2026-11-17.md), which was
committed to this repository on 2026-11-17 at 16:50:47 UTC (commit 5a98703), before any
abstract or full text was opened. The searches and extraction below were run on 2026-11-17
between about 20:40 and 21:30 UTC. The machine-readable table is
[evidence-map-extraction-v1-2026-11-17.csv](evidence-map-extraction-v1-2026-11-17.csv).

## Result in one paragraph

**Not estimable, as the pre-registration requires when fewer than 5 studies qualify.** Of
about 45 published KFRE validations in routine data that the three fixed searches surfaced,
only **two** report both things the question needs: the share of eligible patients who had a
urine albumin test (p) and a numeric observed-versus-predicted comparison among the people
scored. Neither is fully clean (one lacks any uncertainty; the other's p has the wrong
denominator). Seven more report p but no numeric calibration, and several report numeric
calibration but no p. So the primary fit (theta) and the meta-regression were not run, E1 and
E2 are **not tested**, and no power calculation was made (it needs at least 5 studies). The
published literature, as it stands, cannot answer whether testing practice predicts KFRE
miscalibration; a multi-site study that measures p by design (protocol addenda 1-2) is needed.

## Searches (run exactly as fixed, Europe PMC REST, 2026-11-17)

| # | Query | Hits |
|---|---|---|
| 1 | `("kidney failure risk equation") AND (validation OR calibration) AND ("electronic health record" OR "primary care" OR "real-world" OR "health system")` | 228 |
| 2 | `("kidney failure risk equation" OR KFRE) AND (albuminuria OR ACR OR "albumin-to-creatinine") AND (missing OR available OR measured OR testing)` | 476 |
| 3 | records citing PMID 26757465 with "validation" in title or abstract (`CITES:26757465_MED AND (TITLE:validation OR ABSTRACT:validation)`) | 120 (of 584 citing records) |

How screening was done, stated plainly: all pages were requested in one call per query
(pageSize 1000, JSON). An automated fetch-and-summarise step read the titles. For query 1 it
listed every title; for queries 2 and 3 it listed the titles that looked like KFRE
validations, KFRE use in routine data, or albuminuria-testing studies, with instructions to be
inclusive. A handful of records were cut at page-chunk boundaries (about 3). Title screening by
an automated reader can miss studies; this is a limitation of v1.

## Screening and eligibility

Titles kept for abstract review: about 45 KFRE validations or applications in a defined
population. Excluded at title or abstract stage, by the fixed rules:

- Nephrology-clinic or referred cohorts that report no pre-uACR denominator: Kwek 2022
  (36043308), Kang 2020 (32737361), da Silva 2023 (36529658), Larrarte 2024 (38405300),
  Escamilla-Cabrera 2025 (40451701), Temprado Collado 2026 (42453745), Prouvot 2025 (40927378),
  Ali 2021 (34030639), Hundemer 2020/2021 (32928746, 34901569), Thanabalasingam 2022/2024
  (35445190, 38562968), Ramspek 2021 (33685974), Lennartz 2016 (26787778), Creon 2026
  (42422567), Lanot 2026 (42724472), Atiquzzaman 2024 (38327282, advanced-CKD program).
- Transplant, paediatric, disease-specific or protocol-only papers: e.g. Masoud 2025
  (41218935), Ali 2021 (33563222), Tangri 2020 (32549052), Nogueira Machado 2024 (39624526),
  Engole 2026 (42554043), Phannajit 2026 (42152438, protocol).
- Not the KFRE, or a new model only: Black 2025 (41214777, incident-CKD models), Al-Wahsh 2021
  (33944922, competing-risk refit; abstract gives no pre-uACR denominator).
- Duplicate population: Maher 2023 (37798742) uses the same 35,539-person Leicester cohort as
  Major 2019 (31693662).
- Research protocol measuring uACR in everyone: UK Biobank arm of Walker 2026 (41288315).

## Extraction (every number quoted from the paper; full text read through Europe PMC)

### A. Studies with both p and a numeric observed-versus-predicted comparison

| Study | Setting | p (quoted source) | KFRE, horizon, baseline | Observed vs predicted (quoted) | Variance of ln O/E | Issue |
|---|---|---|---|---|---|---|
| Cleary 2026, BMJ Med, PMID 41658049 (SCREAM, Stockholm) | routine primary and outpatient care, CKD 3-4, 2010-2018 | "26 658 (23.0%) patients with a urine albumin to creatinine ratio test" of "116 158 patients [who] met the cohort entry criteria": **p = 0.2295** | 4-variable, non-North American, 5-year implied | "an observed-expected event probability ratio of 2.11"; observed risk "accounting for censoring", deaths censored (Kaplan-Meier) | **not accessible**: no CI, no event count among the 26,229 scored, no observed proportion in the text read | known in advance (the reason for the question) |
| Bravo-Zuniga 2024, BMJ Open, PMID 38184316 (EsSalud primary care, Lima) | 17 primary-care centres, CKD 3a-4, 2013-2017 | "Among the 22 744 patients with CKD screened ... only 13 890 had complete ACR data, of which 7519 patients were eligible due to a diagnosis of CKD 3a-4": 13,890 / 22,744 = 0.611 is across **all CKD stages**, not the eGFR-eligible 3a-4 group the rule requires | 4-variable, non-North American, 2- and 5-year | Table 2, 3a-4: 2-year predicted 0.96%, observed 1.52% (1.24-1.79), **O/E 1.57 (1.39-1.76)**; 5-year predicted 3.18%, observed 3.37% (2.95-3.8), **O/E 1.06 (0.93-1.19)**; Aalen-Johansen | 5-year 0.003955; 2-year 0.003625 (from the CIs) | p denominator does not match the rule (stage-specific denominator only in a figure image) |

### B. Studies that report p but no numeric calibration of the original KFRE

| Study | p (quoted) | What calibration they report (quoted) |
|---|---|---|
| Whitlock 2017, Can J Kidney Health Dis, PMID 28491341 (Manitoba lab database) | "Only 11.7% of stage 3-5 CKD patients had simultaneous ACR measurement" | none; C statistic 0.900 only; 151 of 1512 reached kidney failure in 5 years |
| Naranjo 2021, PMID 35369014 (OptumLabs EHR, 39 organisations) | "The proportion with ACR testing was 19% within the previous 3 years" (abstract) | "Calibration of the four-variable KFRE was acceptable in each group" (abstract); full text not accessible (HTTP 500) |
| Stewart 2026, BJGP, PMID 41344814 (Greater Manchester Care Record) | "71 383 patients with CKD stages 3-5 ... 30 197 (42.3%) had uACR measured in the index year" | no O/E or mean predicted risk; RRT by KFRE band 0.1% / 2.0% / 20.1% (crude), 234 RRT events among the scored by our addition |
| Walker 2025, BMC Nephrol, PMID 40038631 (Glasgow) | "uACR was available in 44.5% of people" (10,874 adults with CKD) | "calibration was also good, but with over-prediction of risk amongst higher risk groups"; outcome includes sustained eGFR < 15 |
| Bravo-Zuniga 2025, BMC Nephrol, PMID 41350644 (EsSalud national, 30,031 CKD G3-4) | "Only 38.4% (11,540) had complete data for the four variables"; "most missing data attributed to ... ACR measurements (60.7%)" | O/E 1.84 (1.7-1.99) at 2 years and 1.06 (1-1.13) at 5 years, but computed on **all 30,031 after multiple imputation of ACR**, so it is not the O/E among those tested; outcome includes ICD-10 N18.5; overlaps the 2024 Lima cohort |
| Goubar 2024, BMJ Open Diab Res Care, PMID 39277182 (two London hospitals, T2D) | "ACR at baseline was available for 81% of individuals" | ICI / E50 / E90 only, UK-KFRE, imputed cohort; "KFRE 'under-predicted' ESKD risk overall" |
| Wang 2019, BMC Nephrol, PMID 31801468 (Singapore primary care) | excluded "patients without ACR (5-year, n = 2403" of 19,857 with eGFR < 60: p = 0.879, but ACR came from a dipstick conversion | median bias (observed minus predicted) 4.5% (-1.4 to 5.9%) at 5 years; no O/E |

### C. Numeric calibration but no pre-uACR denominator in the text read

Major 2019 (31693662; calibration slope 1.080 (1.031-1.129) at 5 years, overprediction in
lower-risk groups; denominator is the whole practice population); Bundo-Luque 2024 (39944820;
"overestimation of RRT risk", denominator only in a figure image); Irish 2023 (37076122;
"predicted risk being more than the observed"; denominators inconsistent); Pan 2024
(38689834; only "312 of the 4,612 participants" had a measured uACR, the rest dipstick-
converted; "calibration was acceptable"); Walker 2026 (41288315; SCREAM O/E 1.09 (1.01-1.18) at
2 years and 1.08 (1.02-1.14) at 5 years with the UK-calibrated KFRE; selection counts only in
supplementary figures, not read).

## Pre-registered analysis

- Qualifying studies: 2 (Cleary 2026 with unknown variance; Bravo-Zuniga 2024 with a
  mismatched denominator). Fewer than 5, so **the analysis is not estimable**.
- Power (fixed to be reported first): not computed, because the design needs at least 5 studies.
- Primary (theta), secondary (Knapp-Hartung meta-regression on p and on m(p)) and the
  sensitivity analyses: **not run**.
- E1 (theta > 0) and E2 (slope < 0): **not tested.** This is not a null result; it is an
  absence of usable published evidence.

## Exploratory observation (not pre-registered; suggestive only)

The same health system appears twice. Both 2026 papers validate the 4-variable KFRE in
Stockholm's SCREAM data at 5 years, but on different scored subsets:

- Cleary et al.: everyone with CKD 3-4 by creatinine who had a uACR within 12 months before to
  3 months after the index date (23% of the eligible); non-North American KFRE; deaths censored.
  O/E **2.11**.
- Walker et al.: people with creatinine **and cystatin C** measured the same day, eGFR < 60 by
  any equation, and any quantifiable albuminuria or proteinuria within 12 months (uPCR and
  dipstick converted); UK-calibrated KFRE; observed risk by Aalen-Johansen (the paper is not
  fully explicit for the O/E tables). O/E **1.08** (1.02-1.14).

The two baselines can be put on one scale. Assuming the UK-calibrated KFRE differs from the
non-North American one only in its baseline survival (5-year 0.9570 versus 0.9365; Major et al.
2019 Table 2, and as Goubar et al. 2024 describe it), a person's UK prediction is
1 - (1 - r)^k with k = ln 0.9570 / ln 0.9365 = 0.66994, where r is the non-North American
prediction. At r = 5% the UK prediction is 3.378% (ratio 1.480); at 20%, 13.886% (1.440); in
the low-risk limit the ratio is 1.4927. So Walker's 1.08 corresponds to about **0.73 to 0.80**
on the non-North American baseline Cleary used. Going the other way, death-censored observed
risk overstates Aalen-Johansen risk in older cohorts (in Bravo-Zuniga 2024, 5-year observed
risk was 4.24% censoring deaths against 3.37% with Aalen-Johansen, a factor 1.26), so Cleary's
2.11 might be about 1.7 on Walker's method. The two scored subsets of one health system then
differ in apparent calibration by a factor of roughly **2 to 3**, not explained by the
baselines (which widen the gap) and only partly by the observed-risk method.

What this does and does not show: it is consistent with the hypothesis that who gets scored
changes how calibrated the KFRE looks within one system. But the two subsets differ in more
than testing (cystatin C ordering, eGFR definition, PCR and dipstick conversion, index-date
rules), so it is not a test of the testing-rate mechanism, and the hand rescaling rests on the
stated assumption. Walker et al.'s supplement reports the non-North American KFRE for SCREAM
directly and would replace the rescaling; it was not accessible here.

## What follows for the study

1. Published validations almost never report the pre-uACR denominator next to calibration, so
   the multi-site study must measure p itself. kfre-transport's `omop.py` already keeps the
   untested people in the cohort and reports p per site.
2. The SCREAM contrast suggests a cheap, sharper design inside one system: score nested subsets
   that differ only in testing (for example uACR within 3, 12 and 24 months) and check whether O/E
   tracks the subset's testing rate as Phi(Phi^-1(p) + theta) / p predicts.
3. A proposed reporting item for KFRE validations in routine data (eligible N before requiring
   uACR, tested N, and O/E in both groups) would make a future evidence map possible.

## Limits

Title screening was automated; three search pages were cut at chunk boundaries; only open full
text and abstracts were read; supplements and figure images were not; one full text (Naranjo
2021) was not accessible. Quoted numbers are as printed, including the papers' own
inconsistencies noted in the CSV.
