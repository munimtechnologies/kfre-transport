# Protocol addendum 1: primary analysis for multi-site (real-data) runs

Dated 2026-11-17, written and published **before any real patient data have been read** by
this project. It fixes the primary measure, the variance, the moderator and the model for
every future multi-site run of `kfretransport.omop` + `kfretransport.pool` (version 0.2.2
or later). It does not change the simulation pre-registration
([sim-prereg-v1-2026-11-16.md](sim-prereg-v1-2026-11-16.md)), which stays frozen. Any later
change will be a new, dated addendum; nothing here will be edited after real data are seen.

Made with AI agents (Claude) working in the name of MV2 of Munim, Inc. No human expert has
reviewed it yet. Research only: information only, not medical advice. Text CC BY 4.0.

## Question

Across health systems, does the share of people with CKD G3-G4 who have a measured urine
albumin-to-creatinine ratio (uACR) predict how far the 4-variable KFRE is from the observed
risk of kidney replacement therapy (KRT) there?

## Primary analysis (per horizon, 2 and 5 years; 5 years is primary)

| Item | Fixed choice |
|---|---|
| Site cohort | `omop.build_cohort` definitions as published in 0.2.0-0.2.2 (eGFR recomputed with CKD-EPI 2021; index = first eGFR in [10, 60) at least 90 days into an unbroken run below 60; adults; KRT before index excluded) |
| KFRE | 4-variable, North American baseline for US sites |
| Site estimate | `tte_aj`: log(O/E), O = Aalen-Johansen cumulative incidence of KRT by the horizon with death as a competing event, among **everyone with a prediction**; E = their mean predicted risk |
| Site variance | delta-method (Greenwood-type) variance of the Aalen-Johansen estimate (`se_log_oe_dm`, new in 0.2.2), E treated as fixed |
| Moderator | the site's testing rate: share of the cohort with a **measured** uACR in the 365 days up to the index date (uPCR conversion never counts as tested) |
| Pooling | DerSimonian-Laird random effects on log O/E; 95% CI and 95% prediction interval for a new site |
| Test | random-effects meta-regression of log O/E on testing rate, method-of-moments tau^2, Knapp-Hartung SE, t with k-2 df, two-sided 5% |
| Minimum | pooling needs >= 2 usable sites, meta-regression >= 3 with different testing rates; sites with 0 KRT events by the horizon are excluded and counted |
| Disclosure | counts from 1 to 19 masked in published per-site rows, plus anything that would give them back (`pool.mask_row`, `--min-cell 20`); this follows N3C's rule that no result about fewer than 20 people leaves its enclave |

**Primary hypothesis (direction stated now):** the meta-regression slope b1 is **negative**:
sites that test fewer people (whose tested people are therefore more selected) show higher
O/E (more under-prediction) among the people they can score. This is the pattern one Swedish
system showed (only 26,229 of 116,158 adults with CKD G3-G4 had a uACR; O/E 2.11 among them;
Cleary et al., BMJ Med 2026, PMID 41658049). Only a multi-site test can tell whether it holds
across systems. A null or positive slope will be reported the same way.

## Secondary analyses (reported next to the primary, never instead of it)

1. Same, with the Poisson variance 1/events.
2. 1 - Kaplan-Meier (death censored) instead of Aalen-Johansen. Expected to give higher O/E
   (Ravani et al., JASN 2019, PMID 31540963).
3. Complete-case O/E (people with a known horizon outcome only), Poisson variance.
4. uPCR-to-uACR conversion on (`--pcr-conversion`; Sumida et al. 2020 crude equation),
   testing rate unchanged.
5. 2-year horizon.

## What would count as support

The primary hypothesis is supported only if b1 < 0 with a 95% CI that excludes 0 in the
primary analysis. Secondary analyses can describe how robust it is but cannot replace a
failed primary.

## Status on this date

No site has run this code on real data. A smoke run on OHDSI's 27-person synthetic Synthea27Nj
set and a pre-registered simulation are queued. A separate, pre-registered public-data study
of the testing mechanism in NHANES is described in
[nhanes-testing-prereg-v1-2026-11-17.md](nhanes-testing-prereg-v1-2026-11-17.md).
