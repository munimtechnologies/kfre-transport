# Pre-registration: do published real-world KFRE validations already show the testing effect? (v1, 2026-11-17)

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it. Text CC BY 4.0. Information only, not medical advice.

**Status when written:** the agents had run one Europe PMC title search (listed below) and had
seen titles only. They had not opened any abstract or full text to extract numbers. The
SCREAM figures (Cleary et al. 2026: 26,229 of 116,158 tested, O/E 2.11) were already known and
are the reason for the question, so SCREAM enters the analysis but is flagged as known in advance.
Extraction starts only after this file is public in the repository. It is never edited
afterwards; any change is a dated addendum.

## Question

Across published external validations of the 4-variable KFRE in routinely collected health
data, does the share of eligible patients who had a uACR (the testing rate p) predict how far
the KFRE under-predicts among those tested? The selection-calibration curve
([selection-curve-2026-11-17.md](selection-curve-2026-11-17.md)) predicts that log O/E rises as
p falls, through the inverse Mills ratio m(p) = phi(Phi^-1(p)) / p.

This is a study-level evidence map. Studies differ in country, baseline, horizon, outcome
definition and case mix, so between-study heterogeneity will be large and power will be low
([power-design-2026-11-17.md](power-design-2026-11-17.md)). It is a first look with public
evidence, not the multi-site test of protocol addenda 1-2.

## Search (fixed)

1. Europe PMC REST, run on the day of extraction (and its date recorded):
   `("kidney failure risk equation") AND (validation OR calibration) AND ("electronic health
   record" OR "primary care" OR "real-world" OR "health system")`, all pages.
2. Europe PMC: `("kidney failure risk equation" OR KFRE) AND (albuminuria OR ACR OR "albumin-to-creatinine")
   AND (missing OR available OR measured OR testing)`, all pages.
3. Europe PMC records that cite Tangri et al. JAMA 2016 (PMID 26757465) and have
   "validation" in the title or abstract.
Any date of publication up to the extraction date; any language with an English abstract.

## Eligibility (fixed)

A study is included if all of these hold:

- It applies the published 4-variable KFRE (2- or 5-year, either baseline) without refitting
  the coefficients. A recalibrated version counts only if the paper also reports the
  original equation's calibration.
- The population is defined from routine data (EHR, primary-care or laboratory databases,
  registries, insurance claims with labs) by eGFR or CKD stage, not by a research protocol
  that measured uACR in everyone.
- The paper reports, or lets the reader compute, **both**: (a) the testing rate p = the number
  with a uACR (or the albuminuria measure the KFRE was computed from) divided by the number
  eligible before requiring that measure, with 0 < p < 1; and (b) observed and predicted risk
  of kidney failure at the same horizon in the people scored. This may be O/E,
  calibration-in-the-large, or observed and mean predicted proportions.

Excluded: development-only papers; papers where the denominator before requiring uACR is
not reported; nephrology-clinic cohorts that report no such denominator; duplicate
populations (the largest or most recent report is kept).

## Extraction (fixed)

Per study: first author, year, PMID, country, setting, years, number eligible, number
tested (p), horizon (5-year preferred, else 2-year), baseline used (North American or
non-North American), outcome definition, whether death was treated as competing or censored,
observed and predicted risk or O/E, and its uncertainty. Variance of log O/E, in this order
of preference: from a reported 95% interval, ((ln U - ln L) / 3.92)^2; else 1 / (observed
events); else 1 / (observed proportion x number scored). Every number is quoted with the
sentence or table it came from. Only open full text or abstracts readable through Europe
PMC are used; anything unreadable is listed as "not accessible", never guessed.

## Analysis (fixed)

- **Primary:** `kfretransport.selection.fit_sites` on (p_k, ln O/E_k, var_k), with tau^2
  taken as the DerSimonian-Laird estimate of the random-effects pool of ln O/E
  (`metrics.py`). Reported: theta, its profile 95% interval, beta0, and whether either hits
  the grid edge.
- **Secondary:** random-effects meta-regression of ln O/E on p with Knapp-Hartung
  inference (`metrics.py`), and on m(p).
- **Sensitivity:** 5-year horizon only; region-matched baseline only (North American
  baseline in North America, the other baseline elsewhere); SCREAM left out.
- **Power, reported first:** before reading the fitted result, the agents report
  `power.curve_power` at theta = 0.3 and theta = 0.69 for the extracted p_k, var_k and tau^2.

## Prediction (fixed)

- **E1:** the primary theta estimate is > 0.
- **E2:** the meta-regression slope of ln O/E on p is < 0.
- No prediction is made about the interval excluding 0, because power is expected to be low.
- If fewer than 5 studies qualify, the analysis is reported as not estimable. The extraction
  table is published anyway.
- Results are published whatever they show. A null or reversed result is reported as such and
  is not explained away. Study-level confounding (country, case mix, outcome definition) can
  create or hide a relation of this kind, so even a positive result is only suggestive. It
  would not support any claim about individual patients.

## Search already run before this file (titles only)

Europe PMC, 2026-11-17, query 1 above, first page of 100 of 228 hits, titles only.
