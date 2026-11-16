# SPDX-License-Identifier: Apache-2.0
"""Published kidney equations, written out in plain Python (stdlib only).

Sources (check them yourself before any use):
- CKD-EPI 2021 creatinine equation, race-free (Inker et al., NEJM 2021).
- CKD-EPI 2009 creatinine equation (Levey et al., Ann Intern Med 2009).
- 4-variable Kidney Failure Risk Equation (Tangri et al., JAMA 2011; recalibrated
  baselines in Tangri et al., JAMA 2016, "Multinational assessment of accuracy of
  equations for predicting risk of kidney failure: a meta-analysis").
  Coefficients and centering as reproduced in the CDC supplement "Item S1.
  Calculation of 5-Year Kidney Failure Risk" (stacks.cdc.gov/view/cdc/94475).
  Note: some secondary sources give the North American 2-year baseline as 0.9751
  instead of 0.9750, and some omit the sex centering (0.5642); this code keeps the
  centered form and 0.9750.

Research code only. Information only, not medical advice, not a medical device.
"""
import math

# KFRE 4-variable: linear predictor terms (coefficient, centering)
KFRE4_AGE = (-0.2201, 7.036)    # per 10 years: age/10 - 7.036
KFRE4_MALE = (0.2467, 0.5642)   # male indicator - 0.5642
KFRE4_EGFR = (-0.5567, 7.222)   # per 5 mL/min/1.73m2: eGFR/5 - 7.222
KFRE4_LNACR = (0.4510, 5.137)   # ln(uACR in mg/g) - 5.137

# Baseline survival S0 by (region, years)
KFRE4_BASELINE = {
    ("north_america", 2): 0.9750,
    ("north_america", 5): 0.9240,
    ("non_north_america", 2): 0.9832,
    ("non_north_america", 5): 0.9365,
}

MG_PER_MMOL_TO_MG_PER_G = 8.84      # uACR: 1 mg/mmol = 8.84 mg/g (approximately)
UMOL_PER_L_PER_MG_DL = 88.4         # serum creatinine: 1 mg/dL = 88.4 umol/L


def creatinine_umol_to_mg_dl(scr_umol_l):
    return scr_umol_l / UMOL_PER_L_PER_MG_DL


def acr_mg_mmol_to_mg_g(acr_mg_mmol):
    return acr_mg_mmol * MG_PER_MMOL_TO_MG_PER_G


def egfr_ckdepi_2021(scr_mg_dl, age, female):
    """Race-free CKD-EPI 2021 creatinine eGFR in mL/min/1.73m2."""
    if scr_mg_dl <= 0 or age < 18:
        raise ValueError("needs creatinine > 0 and age >= 18")
    kappa = 0.7 if female else 0.9
    alpha = -0.241 if female else -0.302
    r = scr_mg_dl / kappa
    value = 142.0 * min(r, 1.0) ** alpha * max(r, 1.0) ** -1.200 * 0.9938 ** age
    if female:
        value *= 1.012
    return value


def egfr_ckdepi_2009(scr_mg_dl, age, female, black=False):
    """CKD-EPI 2009 creatinine eGFR (with its race term) in mL/min/1.73m2.
    Kept only as a pre-declared sensitivity analysis: most KFRE validation work
    used it, and the choice of equation changes KFRE calibration."""
    if scr_mg_dl <= 0 or age < 18:
        raise ValueError("needs creatinine > 0 and age >= 18")
    kappa = 0.7 if female else 0.9
    alpha = -0.329 if female else -0.411
    r = scr_mg_dl / kappa
    value = 141.0 * min(r, 1.0) ** alpha * max(r, 1.0) ** -1.209 * 0.993 ** age
    if female:
        value *= 1.018
    if black:
        value *= 1.159
    return value


def gfr_category(egfr):
    """KDIGO GFR category G1-G5."""
    if egfr >= 90:
        return "G1"
    if egfr >= 60:
        return "G2"
    if egfr >= 45:
        return "G3a"
    if egfr >= 30:
        return "G3b"
    if egfr >= 15:
        return "G4"
    return "G5"


def albuminuria_category(acr_mg_g):
    """KDIGO albuminuria category A1-A3 (uACR in mg/g)."""
    if acr_mg_g < 30:
        return "A1"
    if acr_mg_g <= 300:
        return "A2"
    return "A3"


def kfre4_linear_predictor(age, male, egfr, acr_mg_g):
    if acr_mg_g <= 0:
        raise ValueError("uACR must be > 0 (the equation takes its logarithm)")
    return (KFRE4_AGE[0] * (age / 10.0 - KFRE4_AGE[1])
            + KFRE4_MALE[0] * ((1.0 if male else 0.0) - KFRE4_MALE[1])
            + KFRE4_EGFR[0] * (egfr / 5.0 - KFRE4_EGFR[1])
            + KFRE4_LNACR[0] * (math.log(acr_mg_g) - KFRE4_LNACR[1]))


def risk_from_lp(lp, years=5, region="north_america"):
    s0 = KFRE4_BASELINE[(region, years)]
    return 1.0 - s0 ** math.exp(lp)


def kfre4(age, male, egfr, acr_mg_g, years=5, region="north_america"):
    """4-variable KFRE risk of kidney failure within `years` (2 or 5)."""
    return risk_from_lp(kfre4_linear_predictor(age, male, egfr, acr_mg_g), years, region)
