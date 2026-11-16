# SPDX-License-Identifier: Apache-2.0
"""Build the KFRE validation cohort from an OMOP CDM v5.4 extract (CSV or TSV
files, one per table), standard library only, no network.

Pre-declared definitions (version 0.2.0; changes need a dated note):
- Concept sets are vocabulary codes (conceptsets.py), resolved through the
  site's own CONCEPT (and CONCEPT_ANCESTOR) table. Unknown, non-standard or
  invalid codes are reported, never silently dropped.
- eGFR is always recomputed from serum creatinine with CKD-EPI 2021 (never the
  site's stored eGFR). Same-day creatinine values are averaged.
- Index date: the first date whose eGFR is in [10, 60) and that comes at least
  90 days after the start of an unbroken run of eGFR < 60 results (a result
  >= 60 ends the run). Age >= 18 at every result used; sex must be male or
  female (OMOP 8507/8532) because the KFRE needs it.
- Exclusion: any kidney-replacement-therapy (KRT) record on or before index.
- uACR: the latest value in [index - 365 days, index]; same-day values averaged.
  People without one are kept as 'untested' (the study's point).
- Outcome: first KRT record after index (dialysis or transplant, matched on
  standard OR source concepts). Death before KRT is a competing event; follow-up
  ends at the end of the observation period that contains the index date.
- Horizon outcome (KFRE convention): y = 1 if KRT within the horizon; y = 0 if
  followed through the horizon without KRT or if death came first; None (left
  out) if follow-up ended earlier without either.

Research code only. Information only, not medical advice, not a medical device.
"""
import argparse
import csv
import datetime
import io
import json
import os

from kfretransport import conceptsets as cs
from kfretransport.equations import egfr_ckdepi_2021, kfre4
from kfretransport.metrics import observed_expected

HORIZON_DAYS = {2: 730, 5: 1826}
TABLES = ("person", "measurement", "condition_occurrence", "procedure_occurrence",
          "observation", "death", "observation_period", "concept", "concept_ancestor")


# ---------------------------------------------------------------- reading
def parse_date(text):
    s = (text or "").strip()
    if not s:
        return None
    if len(s) >= 10 and s[4] == "-" and s[7] == "-":
        return datetime.date(int(s[0:4]), int(s[5:7]), int(s[8:10]))
    if len(s) == 8 and s.isdigit():
        return datetime.date(int(s[0:4]), int(s[4:6]), int(s[6:8]))
    raise ValueError("unrecognised date: %r" % s)


def parse_int(text):
    s = (text or "").strip()
    if not s:
        return None
    return int(float(s))


def parse_float(text):
    s = (text or "").strip()
    if not s:
        return None
    return float(s)


def read_table_text(text):
    """Rows of a CSV or TSV table as dicts with lower-case column names.
    Tab-separated if the header has more tabs than commas (Athena files)."""
    if text.startswith(chr(0xFEFF)):    # byte-order mark
        text = text[1:]
    header = text.split("\n", 1)[0]
    if header.count("\t") > header.count(","):
        reader = csv.DictReader(io.StringIO(text), delimiter="\t", quoting=csv.QUOTE_NONE)
    else:
        reader = csv.DictReader(io.StringIO(text))
    rows = []
    for row in reader:
        rows.append({(k or "").strip().lower(): (v or "") for k, v in row.items()})
    return rows


def read_cdm(folder):
    """{table: rows} for every known table found in `folder` (file name
    <table>.csv / .tsv / .txt, any case). Missing tables are empty lists."""
    found = {}
    for name in os.listdir(folder):
        stem, ext = os.path.splitext(name)
        if ext.lower() in (".csv", ".tsv", ".txt") and stem.lower() in TABLES:
            with open(os.path.join(folder, name), encoding="utf-8", newline="") as f:
                found[stem.lower()] = read_table_text(f.read())
    return {t: found.get(t, []) for t in TABLES}


# ---------------------------------------------------------------- vocabulary
class Vocabulary:
    def __init__(self, concept_rows, ancestor_rows=()):
        self.by_id = {}
        self.by_code = {}
        for r in concept_rows:
            cid = parse_int(r.get("concept_id"))
            if cid is None:
                continue
            self.by_id[cid] = r
            key = (r.get("vocabulary_id", "").strip(), r.get("concept_code", "").strip())
            self.by_code[key] = cid
        self.descendants = {}
        for r in ancestor_rows:
            a = parse_int(r.get("ancestor_concept_id"))
            d = parse_int(r.get("descendant_concept_id"))
            if a is not None and d is not None:
                self.descendants.setdefault(a, set()).add(d)

    def code(self, concept_id):
        r = self.by_id.get(concept_id)
        if r is None:
            return None
        return r.get("vocabulary_id", "").strip(), r.get("concept_code", "").strip()

    def is_standard(self, concept_id):
        r = self.by_id.get(concept_id)
        return r is not None and r.get("standard_concept", "").strip() == "S"

    def invalid_reason(self, concept_id):
        r = self.by_id.get(concept_id)
        return "" if r is None else r.get("invalid_reason", "").strip()


def resolve(vocab, entries, expect_standard=True):
    """Concept ids for a list of (vocabulary_id, code, descendants, display).
    Returns (set of ids, list of notes)."""
    ids = set()
    notes = []
    for vocabulary_id, code, descendants, _display in entries:
        cid = vocab.by_code.get((vocabulary_id, code))
        if cid is None:
            notes.append("missing %s %s" % (vocabulary_id, code))
            continue
        reason = vocab.invalid_reason(cid)
        if reason:
            notes.append("invalid(%s) %s %s" % (reason, vocabulary_id, code))
        if expect_standard and not vocab.is_standard(cid):
            notes.append("non-standard %s %s" % (vocabulary_id, code))
        ids.add(cid)
        if descendants:
            ids |= vocab.descendants.get(cid, set())
    return ids, notes


def unit_key(vocab, unit_concept_id, unit_source_value):
    """The unit as a lower-case string without spaces: the UCUM code of the unit
    concept if there is one, else the raw unit text. Returns (key, used_text)."""
    if unit_concept_id:
        c = vocab.code(unit_concept_id)
        if c is not None and c[1]:
            return c[1].replace(" ", "").lower(), False
    return (unit_source_value or "").replace(" ", "").lower(), True


# ---------------------------------------------------------------- helpers
def age_on(birth, day):
    return day.year - birth.year - ((day.month, day.day) < (birth.month, birth.day))


def _birth(row):
    y = parse_int(row.get("year_of_birth"))
    if y is None:
        return None
    m = parse_int(row.get("month_of_birth")) or 7    # mid-year if unknown
    d = parse_int(row.get("day_of_birth")) or 1
    return datetime.date(y, m, d)


def _lab_series(rows, ids, units, value_range, vocab, counts, label):
    """{person_id: {date: [converted values]}} for measurement rows in `ids`."""
    lo, hi = value_range
    out = {}
    for r in rows:
        if parse_int(r.get("measurement_concept_id")) not in ids:
            continue
        value = parse_float(r.get("value_as_number"))
        day = parse_date(r.get("measurement_date"))
        pid = parse_int(r.get("person_id"))
        if value is None or day is None or pid is None:
            counts[label + "_no_value"] = counts.get(label + "_no_value", 0) + 1
            continue
        key, used_text = unit_key(vocab, parse_int(r.get("unit_concept_id")),
                                  r.get("unit_source_value"))
        factor = units.get(key)
        if factor is None:
            unknown = counts.setdefault(label + "_unknown_units", {})
            unknown[key or "(blank)"] = unknown.get(key or "(blank)", 0) + 1
            continue
        if used_text:
            counts[label + "_unit_from_text"] = counts.get(label + "_unit_from_text", 0) + 1
        v = value * factor
        if label == "uacr":
            plausible = lo < v <= hi     # uACR must be > 0 (the KFRE takes its log)
        else:
            plausible = lo <= v <= hi
        if not plausible:
            counts[label + "_implausible"] = counts.get(label + "_implausible", 0) + 1
            continue
        out.setdefault(pid, {}).setdefault(day, []).append(v)
    return out


def _events(cdm, std_ids, src_ids):
    """{person_id: sorted list of (date, path)} from condition, procedure and
    observation rows; path is 'standard', 'source' or 'both'."""
    specs = (("condition_occurrence", "condition_concept_id", "condition_source_concept_id",
              "condition_start_date"),
             ("procedure_occurrence", "procedure_concept_id", "procedure_source_concept_id",
              "procedure_date"),
             ("observation", "observation_concept_id", "observation_source_concept_id",
              "observation_date"))
    out = {}
    for table, c_std, c_src, c_date in specs:
        for r in cdm.get(table, []):
            s = parse_int(r.get(c_std)) in std_ids
            o = parse_int(r.get(c_src)) in src_ids
            if not (s or o):
                continue
            day = parse_date(r.get(c_date))
            pid = parse_int(r.get("person_id"))
            if day is None or pid is None:
                continue
            path = "both" if (s and o) else ("standard" if s else "source")
            out.setdefault(pid, []).append((day, path))
    for v in out.values():
        v.sort()
    return out


# ---------------------------------------------------------------- cohort
def build_cohort(cdm, region="north_america", egfr_range=(10.0, 60.0), confirm_days=90,
                 acr_lookback_days=365):
    """Returns {"rows": [...], "excluded": {...}, "counts": {...}, "notes": [...]}."""
    vocab = Vocabulary(cdm.get("concept", []), cdm.get("concept_ancestor", []))
    notes = []
    creat_ids, n = resolve(vocab, cs.SERUM_CREATININE)
    notes += n
    acr_ids, n = resolve(vocab, cs.UACR)
    notes += n
    krt_std, n = resolve(vocab, cs.KRT_STANDARD)
    notes += n
    krt_src, n = resolve(vocab, cs.KRT_SOURCE, expect_standard=False)
    notes += n

    counts = {}
    creat = _lab_series(cdm.get("measurement", []), creat_ids, cs.CREATININE_UNITS,
                        cs.CREATININE_RANGE_MG_DL, vocab, counts, "creatinine")
    acr = _lab_series(cdm.get("measurement", []), acr_ids, cs.UACR_UNITS,
                      cs.UACR_RANGE_MG_G, vocab, counts, "uacr")
    krt = _events(cdm, krt_std, krt_src)
    deaths = {}
    for r in cdm.get("death", []):
        pid, day = parse_int(r.get("person_id")), parse_date(r.get("death_date"))
        if pid is not None and day is not None:
            deaths[pid] = min(day, deaths.get(pid, day))
    periods = {}
    for r in cdm.get("observation_period", []):
        pid = parse_int(r.get("person_id"))
        a = parse_date(r.get("observation_period_start_date"))
        b = parse_date(r.get("observation_period_end_date"))
        if pid is not None and a is not None and b is not None:
            periods.setdefault(pid, []).append((a, b))

    lo, hi = egfr_range
    rows = []
    excluded = {}

    def drop(reason):
        excluded[reason] = excluded.get(reason, 0) + 1

    for p in sorted(cdm.get("person", []), key=lambda r: parse_int(r.get("person_id")) or 0):
        pid = parse_int(p.get("person_id"))
        gender = parse_int(p.get("gender_concept_id"))
        if gender not in (cs.GENDER_MALE, cs.GENDER_FEMALE):
            drop("sex_not_male_or_female")
            continue
        male = gender == cs.GENDER_MALE
        birth = _birth(p)
        if birth is None:
            drop("no_birth_year")
            continue
        index = index_egfr = None
        run_start = None
        for day in sorted(creat.get(pid, {})):
            age = age_on(birth, day)
            if age < 18:
                continue
            values = creat[pid][day]
            egfr = egfr_ckdepi_2021(sum(values) / len(values), age, female=not male)
            if egfr >= hi:
                run_start = None
                continue
            if run_start is None:
                run_start = day
            if (day - run_start).days >= confirm_days and egfr >= lo:
                index, index_egfr = day, egfr
                break
        if index is None:
            drop("no_confirmed_egfr_10_to_60")
            continue
        events = krt.get(pid, [])
        if any(d <= index for d, _ in events):
            drop("krt_on_or_before_index")
            continue
        period = [ab for ab in periods.get(pid, []) if ab[0] <= index <= ab[1]]
        if not period:
            drop("index_outside_observation_period")
            continue
        end = max(b for _, b in period)
        acr_value = None
        window = [d for d in acr.get(pid, {}) if 0 <= (index - d).days <= acr_lookback_days]
        if window:
            last = max(window)
            acr_value = sum(acr[pid][last]) / len(acr[pid][last])
        after = [(d, path) for d, path in events if d > index]
        kf = after[0] if after else None
        death = deaths.get(pid)
        if kf is not None and kf[0] <= end and (death is None or kf[0] <= death):
            event, stop, path = 1, kf[0], kf[1]
        elif death is not None and death <= end:
            event, stop, path = 2, death, ""
        else:
            event, stop, path = 0, end, ""
        age = age_on(birth, index)
        row = {"person_id": pid, "index_date": index.isoformat(), "age": age, "male": male,
               "egfr": index_egfr, "tested": acr_value is not None, "uacr_mg_g": acr_value,
               "event": event, "time_days": (stop - index).days, "outcome_path": path}
        for years, h in HORIZON_DAYS.items():
            if event == 1 and row["time_days"] <= h:
                y = 1
            elif row["time_days"] >= h or event == 2:
                y = 0
            else:
                y = None
            row["y%d" % years] = y
            row["kfre%d" % years] = (kfre4(age, male, index_egfr, acr_value, years, region)
                                     if acr_value is not None else None)
        rows.append(row)
    return {"rows": rows, "excluded": excluded, "counts": counts, "notes": notes}


def site_summary(rows, years):
    """Testing rate, and O/E of the KFRE among tested people with a horizon outcome."""
    y_key, p_key = "y%d" % years, "kfre%d" % years
    tested = [r for r in rows if r["tested"]]
    scored = [r for r in tested if r[y_key] is not None]
    untested = [r for r in rows if not r["tested"] and r[y_key] is not None]
    out = {"years": years, "n_cohort": len(rows), "n_tested": len(tested),
           "testing_rate": (len(tested) / len(rows)) if rows else None,
           "n_scored": len(scored), "events_scored": sum(r[y_key] for r in scored),
           "expected_scored": sum(r[p_key] for r in scored),
           "n_untested_with_outcome": len(untested),
           "events_untested": sum(r[y_key] for r in untested)}
    out["oe"] = (observed_expected([r[y_key] for r in scored], [r[p_key] for r in scored])
                 if scored else None)
    return out


def vocabulary_drift(cdm):
    """State of every VOCABULARY_WATCH code in this site's CONCEPT table, and
    which standard concepts this site's ICD10CM N18.6 records landed on."""
    vocab = Vocabulary(cdm.get("concept", []), ())
    watch = []
    for vocabulary_id, code, display, snomed_status in cs.VOCABULARY_WATCH:
        cid = vocab.by_code.get((vocabulary_id, code))
        n = 0
        if cid is not None:
            n = sum(1 for r in cdm.get("condition_occurrence", [])
                    if parse_int(r.get("condition_concept_id")) == cid)
        watch.append({"code": "%s %s" % (vocabulary_id, code), "display": display,
                      "snomed_2026_07": snomed_status, "in_site_vocabulary": cid is not None,
                      "standard": cid is not None and vocab.is_standard(cid),
                      "invalid_reason": vocab.invalid_reason(cid) if cid is not None else "",
                      "condition_records": n})
    n186 = vocab.by_code.get(("ICD10CM", "N18.6"))
    landed = {}
    if n186 is not None:
        for r in cdm.get("condition_occurrence", []):
            if parse_int(r.get("condition_source_concept_id")) == n186:
                c = vocab.code(parse_int(r.get("condition_concept_id")))
                key = "%s %s" % c if c else "unmapped"
                landed[key] = landed.get(key, 0) + 1
    return {"watch": watch, "icd10cm_n18_6_maps_to": landed}


def main(argv=None):
    ap = argparse.ArgumentParser(description="KFRE cohort from an OMOP CDM CSV extract "
                                 "(research only; not medical advice)")
    ap.add_argument("--cdm-dir", required=True)
    ap.add_argument("--out", required=True, help="JSON file to write (summary only)")
    ap.add_argument("--region", default="north_america",
                    choices=["north_america", "non_north_america"])
    a = ap.parse_args(argv)
    cdm = read_cdm(a.cdm_dir)
    cohort = build_cohort(cdm, region=a.region)
    result = {"excluded": cohort["excluded"], "counts": cohort["counts"],
              "concept_notes": cohort["notes"],
              "summary": [site_summary(cohort["rows"], y) for y in (2, 5)],
              "vocabulary_drift": vocabulary_drift(cdm)}
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=1, sort_keys=True)
    return result


if __name__ == "__main__":
    main()
