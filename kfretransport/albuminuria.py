# SPDX-License-Identifier: Apache-2.0
"""Estimated urine albumin-to-creatinine ratio (ACR) from a urine
protein-to-creatinine ratio (PCR), for a pre-declared SENSITIVITY analysis only.

Source: Sumida K, Nadkarni GN, Grams ME, et al. Conversion of urine
protein-creatinine ratio or urine dipstick protein to urine albumin-creatinine
ratio for use in chronic kidney disease screening and prognosis: an individual
participant-based meta-analysis. Ann Intern Med 2020;173(6):426-435.
doi:10.7326/M20-0529, PMID 32658569, PMC7780415 (CKD Prognosis Consortium,
919,383 adults, 33 cohorts).

Crude equation (PCR and ACR in mg/g, natural logs, linear spline in log PCR with
knots at 50 and 500 mg/g):

    ACR = exp(5.3920 + 0.3072 ln(min(PCR/50, 1))
                     + 1.5793 ln(max(min(PCR/500, 1), 0.1))
                     + 1.1266 ln(max(PCR/500, 1)))

Provenance, honestly stated: the primary PDF could not be opened by the agents
that wrote this code. The coefficients were taken from a search-engine extract
of Table 2 and then checked against five figures that the paper's text and
Table 3 report for the crude model (as quoted by NephJC): PCR 150 -> ACR 33,
PCR 500 -> ACR 220, ACR 30 <-> PCR 142, ACR 300 <-> PCR 660, and a 2.99-fold
ACR increase per doubling of PCR between 50 and 500 mg/g. All five are
reproduced (see test_albuminuria.py). The lowest segment (PCR < 50 mg/g, slope
0.3072) is not checked by any reported figure, and the paper itself says the
PCR-ACR association is inconsistent there. The adjusted model (sex, diabetes,
hypertension) and the dipstick conversion are NOT implemented. Check against
the paper before any real use.

Research code only. Information only, not medical advice, not a medical device.
"""
import math

SUMIDA_2020_CRUDE = (5.3920, 0.3072, 1.5793, 1.1266)
LOW_PCR_MG_G = 50.0     # below this the paper reports an inconsistent PCR-ACR link


def acr_from_pcr(pcr_mg_g, coefficients=SUMIDA_2020_CRUDE):
    """Estimated ACR (mg/g) from PCR (mg/g) with the Sumida 2020 crude equation."""
    if pcr_mg_g is None or pcr_mg_g <= 0:
        raise ValueError("PCR must be > 0 mg/g")
    b0, b1, b2, b3 = coefficients
    x = pcr_mg_g
    return math.exp(b0 + b1 * math.log(min(x / 50.0, 1.0))
                    + b2 * math.log(max(min(x / 500.0, 1.0), 0.1))
                    + b3 * math.log(max(x / 500.0, 1.0)))
