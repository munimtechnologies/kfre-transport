# SPDX-License-Identifier: Apache-2.0
"""What US-style uACR testing does to the KFRE, measured in NHANES 2007-2018,
where every examined adult had a urine albumin-to-creatinine ratio (uACR).

Pre-registered in findings/nhanes-testing-prereg-v1-2026-11-17.md before any
NHANES data were read. Public data only (NCHS NHANES and the NCHS public-use
linked mortality files through 2019-12-31).

Cohort: adults 20+, not pregnant, with an exam weight, serum creatinine and a
uACR, whose CKD-EPI 2021 eGFR is in [15, 60) (CKD G3-G4 by a single eGFR).
Testing policies give each person a probability of having had a uACR test in
routine care:
  P0 universal: 1
  P1 us_claims: 0.63 with diagnosed diabetes, 0.25 without (Chu et al., BMC
     Prim Care 2022;23:299, PMID 36434513: insured US adults with CKD, 2015-2019)
  P2 guideline: 1 with diagnosed diabetes or diagnosed hypertension, else 0
  P3 diabetes_only: 1 with diagnosed diabetes, else 0
  P1q (exploratory curve): 0.63 with diagnosed diabetes, q without
Estimates use expected weights (survey weight x probability), i.e. the
large-population limit of random testing within each stratum, so nothing is
random and every number is reproducible.

Standard errors: delete-one-PSU jackknife over NHANES masked variance units
(SDMVSTRA, SDMVPSU); 95% CI with t on (PSUs - strata) degrees of freedom.

Pandas is needed only to read the SAS XPORT files; everything else is stdlib.
Research only. Information only, not medical advice, not a medical device.
"""
import argparse
import json
import math
import os
import time
import urllib.request

from kfretransport import CODE_VERSION
from kfretransport.equations import egfr_ckdepi_2021, kfre4
from kfretransport.metrics import t975

CYCLES = (("2007-2008", "2007", "E"), ("2009-2010", "2009", "F"), ("2011-2012", "2011", "G"),
          ("2013-2014", "2013", "H"), ("2015-2016", "2015", "I"), ("2017-2018", "2017", "J"))
TABLES = ("DEMO", "BIOPRO", "ALB_CR", "DIQ", "BPQ")
NHANES_URL = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s_%s.xpt"
MORT_URL = ("https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/"
            "NHANES_%s_MORT_2019_PUBLIC.dat")
KEEP = ("SEQN", "RIDAGEYR", "RIAGENDR", "RIDEXPRG", "WTMEC2YR", "SDMVSTRA", "SDMVPSU",
        "LBXSCR", "URDACT", "DIQ010", "BPQ020")
US_CLAIMS_DM, US_CLAIMS_NO_DM = 0.63, 0.25
THRESHOLDS = (0.05, 0.03)          # KDIGO 2024 practice point 2.2.1: 3-5% at 5 years
Q_CURVE = (0.05, 0.10, 0.25, 0.50, 0.75, 1.0)


def policy_prob(name, q=None):
    """Testing probability function for a named policy."""
    if name == "P0_universal":
        return lambda r: 1.0
    if name == "P1_us_claims":
        return lambda r: US_CLAIMS_DM if r["dm"] else US_CLAIMS_NO_DM
    if name == "P2_guideline":
        return lambda r: 1.0 if (r["dm"] or r["htn"]) else 0.0
    if name == "P3_diabetes_only":
        return lambda r: 1.0 if r["dm"] else 0.0
    if name == "P1q":
        return lambda r: US_CLAIMS_DM if r["dm"] else q
    raise ValueError("unknown policy " + name)


POLICIES = ("P0_universal", "P1_us_claims", "P2_guideline", "P3_diabetes_only")


def _num(v):
    """None for missing / NaN; SAS XPORT zeros can arrive as 5.4e-79."""
    if v is None:
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if v != v:
        return None
    return 0.0 if abs(v) < 1e-30 else v


def make_row(r, n_cycles, cycle=""):
    """One merged NHANES record -> (cohort row or None, exclusion reason or None)."""
    age = _num(r.get("RIDAGEYR"))
    if age is None or age < 20:
        return None, "age_under_20"
    sex = _num(r.get("RIAGENDR"))
    if sex not in (1.0, 2.0):
        return None, "sex_missing"
    if _num(r.get("RIDEXPRG")) == 1.0:
        return None, "pregnant"
    w = _num(r.get("WTMEC2YR"))
    if not w or w <= 0:
        return None, "no_exam_weight"
    scr = _num(r.get("LBXSCR"))
    if not scr or scr <= 0:
        return None, "no_creatinine"
    acr = _num(r.get("URDACT"))
    if not acr or acr <= 0:
        return None, "no_uacr"
    male = sex == 1.0
    egfr = egfr_ckdepi_2021(scr, age, not male)
    if not 15.0 <= egfr < 60.0:
        return None, "egfr_not_15_to_59"
    row = {"seqn": int(_num(r.get("SEQN")) or 0), "cycle": cycle, "w": w / n_cycles,
           "stratum": int(_num(r.get("SDMVSTRA")) or 0), "psu": int(_num(r.get("SDMVPSU")) or 0),
           "age": age, "male": male, "egfr": egfr, "uacr": acr,
           "dm": _num(r.get("DIQ010")) == 1.0, "htn": _num(r.get("BPQ020")) == 1.0,
           "kfre5": kfre4(age, male, egfr, acr, 5, "north_america"),
           "kfre2": kfre4(age, male, egfr, acr, 2, "north_america"),
           "died": None, "kdeath": None}
    return row, None


def cell_of(row):
    """Imputation cell: eGFR band x age band x sex."""
    e = row["egfr"]
    band = "15-29" if e < 30 else ("30-44" if e < 45 else "45-59")
    return band, ("65+" if row["age"] >= 65 else "<65"), row["male"]


def imputation_ratio(rows, wt, wu):
    """Impute log uACR of the untested (weights wu) as the wt-weighted mean log
    uACR of tested donors in the same cell (fallback: same eGFR band). Returns
    (sum wu * KFRE5 with imputed uACR / sum wu * KFRE5 with measured uACR,
    number of fallbacks, number with no donor)."""
    s_cell, n_cell, s_band, n_band = {}, {}, {}, {}
    for x, r in zip(wt, rows):
        if x > 0:
            c = cell_of(r)
            lg = math.log(r["uacr"])
            s_cell[c] = s_cell.get(c, 0.0) + x * lg
            n_cell[c] = n_cell.get(c, 0.0) + x
            s_band[c[0]] = s_band.get(c[0], 0.0) + x * lg
            n_band[c[0]] = n_band.get(c[0], 0.0) + x
    num = den = 0.0
    fallbacks = no_donor = 0
    for x, r in zip(wu, rows):
        if x <= 0:
            continue
        c = cell_of(r)
        if n_cell.get(c, 0.0) > 0:
            mu = s_cell[c] / n_cell[c]
        elif n_band.get(c[0], 0.0) > 0:
            mu = s_band[c[0]] / n_band[c[0]]
            fallbacks += 1
        else:
            no_donor += 1
            continue
        num += x * kfre4(r["age"], r["male"], r["egfr"], math.exp(mu), 5, "north_america")
        den += x * r["kfre5"]
    return (num / den if den > 0 else None), fallbacks, no_donor


def policy_estimates(rows, prob, threshold=0.05, weights=None):
    """Expected-weight estimates for one testing policy. Returns a dict:
    tested_share, high_scored_share (share of KFRE >= threshold people who have
    a test), high_unscored_weighted (weighted count of them without a test),
    mean_kfre_tested / mean_kfre_all and their ratio, imputation_ratio (mean
    5-y KFRE of the untested with uACR imputed from tested people in the same
    cell, over their mean KFRE with their own measured uACR), and
    kdeath_unscored_share (share of kidney-cause deaths among the untested)."""
    w = [r["w"] for r in rows] if weights is None else weights
    p = [prob(r) for r in rows]
    W = sum(w)
    out = {"tested_share": None, "high_share": None, "high_scored_share": None,
           "high_unscored_weighted": None, "mean_kfre_all": None, "mean_kfre_tested": None,
           "tested_all_ratio": None, "imputation_ratio": None, "imputation_fallbacks": 0,
           "imputation_no_donor": 0, "imputation_ratio_uniform": None,
           "selection_ratio": None, "kdeath_unscored_share": None}
    if W <= 0:
        return out
    wt = [wi * pi for wi, pi in zip(w, p)]
    wu = [wi * (1.0 - pi) for wi, pi in zip(w, p)]
    hi = [r["kfre5"] >= threshold for r in rows]
    W_hi = sum(wi for wi, h in zip(w, hi) if h)
    T = sum(wt)
    out["tested_share"] = T / W
    out["high_share"] = W_hi / W
    if W_hi > 0:
        out["high_scored_share"] = sum(x for x, h in zip(wt, hi) if h) / W_hi
        out["high_unscored_weighted"] = sum(x for x, h in zip(wu, hi) if h)
    out["mean_kfre_all"] = sum(wi * r["kfre5"] for wi, r in zip(w, rows)) / W
    if T > 0:
        out["mean_kfre_tested"] = sum(x * r["kfre5"] for x, r in zip(wt, rows)) / T
        out["tested_all_ratio"] = out["mean_kfre_tested"] / out["mean_kfre_all"]
    ratio, fb, nd = imputation_ratio(rows, wt, wu)
    out["imputation_ratio"], out["imputation_fallbacks"], out["imputation_no_donor"] = ratio, fb, nd
    # reference: random (non-selective) testing, where donors and untested people
    # both carry the survey weight. Cell-mean imputation of log uACR is biased even
    # then (the KFRE is not linear in log uACR), so selection is judged against it.
    ref, _, _ = imputation_ratio(rows, w, w)
    out["imputation_ratio_uniform"] = ref
    if ratio is not None and ref:
        out["selection_ratio"] = ratio / ref
    kd = [r.get("kdeath") for r in rows]
    W_kd = sum(wi for wi, k in zip(w, kd) if k)
    if W_kd > 0:
        out["kdeath_unscored_share"] = sum(x for x, k in zip(wu, kd) if k) / W_kd
    return out


def jackknife(rows, stat):
    """Delete-one-PSU jackknife (JKn). stat(weights) -> dict of numbers or None.
    Returns (estimates, standard errors, info)."""
    base = [r["w"] for r in rows]
    theta = stat(base)
    strata = {}
    for r in rows:
        strata.setdefault(r["stratum"], set()).add(r["psu"])
    acc = {k: 0.0 for k, v in theta.items() if isinstance(v, (int, float)) and v is not None
           and not isinstance(v, bool)}
    missing = {k: 0 for k in acc}
    single = 0
    for h, psus in sorted(strata.items()):
        n_h = len(psus)
        if n_h < 2:
            single += 1
            continue
        for j in sorted(psus):
            wts = [0.0 if (r["stratum"] == h and r["psu"] == j)
                   else (wb * n_h / (n_h - 1.0) if r["stratum"] == h else wb)
                   for r, wb in zip(rows, base)]
            rep = stat(wts)
            for k in acc:
                if rep.get(k) is None:
                    missing[k] += 1
                else:
                    acc[k] += (n_h - 1.0) / n_h * (rep[k] - theta[k]) ** 2
    se = {k: math.sqrt(v) for k, v in acc.items()}
    n_psu = sum(len(s) for s in strata.values())
    info = {"strata": len(strata), "psus": n_psu, "df": n_psu - len(strata),
            "strata_with_one_psu": single, "replicates_missing": missing}
    return theta, se, info


def with_ci(theta, se, df):
    out = {}
    tq = t975(df) if df >= 1 else None
    for k, v in theta.items():
        if k in se and tq is not None:
            out[k] = {"est": v, "se": se[k], "ci95": [v - tq * se[k], v + tq * se[k]]}
        else:
            out[k] = {"est": v}
    return out


def parse_mortality_line(line):
    """NCHS public-use linked mortality file (2019), fixed width:
    SEQN 1-6, ELIGSTAT 15, MORTSTAT 16, UCOD_LEADING 17-19."""
    def f(a, b):
        s = line[a - 1:b].strip()
        return int(s) if s.isdigit() else None
    return {"seqn": f(1, 6), "eligstat": f(15, 15), "mortstat": f(16, 16),
            "ucod_leading": f(17, 19)}


# ------------------------------------------------------------- I/O (job only)
def download(data_dir, tries=4):
    os.makedirs(data_dir, exist_ok=True)
    urls = []
    for cyc, year, suf in CYCLES:
        for t in TABLES:
            urls.append(NHANES_URL % (year, t, suf))
        urls.append(MORT_URL % cyc.replace("-", "_"))
    got = []
    for u in urls:
        dest = os.path.join(data_dir, u.rsplit("/", 1)[1])
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            got.append(dest)
            continue
        for k in range(tries):
            try:
                req = urllib.request.Request(u, headers={
                    "User-Agent": "kfre-transport research (MV2, Munim, Inc.; mv2@munimtech.com)"})
                with urllib.request.urlopen(req, timeout=120) as resp, open(dest, "wb") as fh:
                    fh.write(resp.read())
                got.append(dest)
                break
            except Exception as exc:  # noqa: BLE001 - report and retry
                print("download failed (%d/%d) %s: %s" % (k + 1, tries, u, exc))
                time.sleep(30 * (k + 1))
    return urls, got


def load(data_dir):
    """Merged records per cycle (dicts keyed by NHANES variable names) and the
    mortality status by SEQN."""
    import pandas as pd  # only the job needs pandas
    records = []
    for cyc, _year, suf in CYCLES:
        merged = {}
        for t in TABLES:
            path = os.path.join(data_dir, "%s_%s.xpt" % (t, suf))
            df = pd.read_sas(path, format="xport")
            cols = [c for c in df.columns if c in KEEP]
            for rec in df[cols].to_dict("records"):
                s = _num(rec.get("SEQN"))
                if s is None:
                    continue
                if t == "DEMO":
                    merged[int(s)] = dict(rec)
                elif int(s) in merged:
                    merged[int(s)].update({k: v for k, v in rec.items() if k != "SEQN"})
        for rec in merged.values():
            rec["_cycle"] = cyc
            records.append(rec)
    mort = {}
    for cyc, _y, _s in CYCLES:
        path = os.path.join(data_dir, "NHANES_%s_MORT_2019_PUBLIC.dat" % cyc.replace("-", "_"))
        if os.path.exists(path):
            with open(path, encoding="ascii", errors="replace") as fh:
                for line in fh:
                    m = parse_mortality_line(line)
                    if m["seqn"] is not None:
                        mort[m["seqn"]] = m
    return records, mort


def analyze(records, mort):
    rows, excluded = [], {}
    for rec in records:
        row, why = make_row(rec, len(CYCLES), rec.get("_cycle", ""))
        if row is None:
            excluded[why] = excluded.get(why, 0) + 1
            continue
        m = mort.get(row["seqn"])
        if m and m["eligstat"] == 1 and m["mortstat"] in (0, 1):
            row["died"] = m["mortstat"] == 1
            row["kdeath"] = m["mortstat"] == 1 and m["ucod_leading"] == 9
        rows.append(row)
    W = sum(r["w"] for r in rows)
    desc = {"n_unweighted": len(rows), "weighted_adults": W,
            "n_by_cycle": {}, "dm_share": sum(r["w"] for r in rows if r["dm"]) / W,
            "htn_share": sum(r["w"] for r in rows if r["htn"]) / W,
            "dm_or_htn_share": sum(r["w"] for r in rows if r["dm"] or r["htn"]) / W,
            "mortality_linked": sum(1 for r in rows if r["died"] is not None),
            "deaths": sum(1 for r in rows if r["died"]),
            "kidney_deaths": sum(1 for r in rows if r["kdeath"])}
    for r in rows:
        desc["n_by_cycle"][r["cycle"]] = desc["n_by_cycle"].get(r["cycle"], 0) + 1
    for label, sel in (("dm", lambda r: r["dm"]), ("no_dm", lambda r: not r["dm"])):
        sub = [r for r in rows if sel(r)]
        ws = sum(r["w"] for r in sub)
        desc["mean_log_uacr_" + label] = (sum(r["w"] * math.log(r["uacr"]) for r in sub) / ws
                                          if ws > 0 else None)
    results = {"code_version": CODE_VERSION, "excluded": excluded, "descriptive": desc,
               "policies": {}, "q_curve": {}}
    info = None
    for thr in THRESHOLDS:
        for name in POLICIES:
            prob = policy_prob(name)
            theta, se, info = jackknife(
                rows, lambda wts, prob=prob, thr=thr: policy_estimates(rows, prob, thr, wts))
            results["policies"]["%s@%g" % (name, thr)] = with_ci(theta, se, info["df"])
    for q in Q_CURVE:
        prob = policy_prob("P1q", q)
        theta, se, info = jackknife(
            rows, lambda wts, prob=prob: policy_estimates(rows, prob, 0.05, wts))
        results["q_curve"]["%g" % q] = with_ci(theta, se, info["df"])
    results["design"] = info
    return results


def main(argv=None):
    ap = argparse.ArgumentParser(description="NHANES uACR testing-policy study "
                                 "(research only; not medical advice)")
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--no-download", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(a.out_dir, exist_ok=True)
    report = {"code_version": CODE_VERSION}
    if not a.no_download:
        urls, got = download(a.data_dir)
        report["files_expected"], report["files_downloaded"] = len(urls), len(got)
    records, mort = load(a.data_dir)
    report["records"], report["mortality_records"] = len(records), len(mort)
    report.update(analyze(records, mort))
    with open(os.path.join(a.out_dir, "nhanes-results.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=1, sort_keys=True)
    return report


if __name__ == "__main__":
    main()
