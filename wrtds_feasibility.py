"""
WRTDS Feasibility Diagnostic
=============================
Checks every site in your NEAA site list against NWIS discharge records
and reports which sites qualify for WRTDS flow-normalized trend analysis.

WRTDS requirements (following Hirsch & De Cicco 2015):
  - >= 100 paired water quality + discharge observations
  - >= 10 years of overlapping WQ and discharge record
  - Daily discharge record with < 20% gaps over the WQ period
  - Discharge gauge is the same site (or a nearby mapped gauge)

Outputs:
  - wrtds_feasibility.csv   : full results table, one row per site
  - wrtds_feasibility.txt   : human-readable summary report
  - Console summary         : tier breakdown and recommended approach

Usage:
    conda activate neaa
    python wrtds_feasibility.py

Edit SITE_LIST_PATH and PROJECT_DIR below before running.
"""

import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from tqdm import tqdm

try:
    import dataretrieval.nwis as nwis
    import dataretrieval.wqp as wqp
except ImportError:
    sys.exit("Install dataretrieval:  pip install dataretrieval")

# =============================================================================
# CONFIGURATION — edit these to match your setup
# =============================================================================
PROJECT_DIR    = Path.home() / "Desktop" / "NEAA"
SITE_LIST_PATH = PROJECT_DIR / "WQP_sites_manual_srp_update6.xlsx"
OUTPUT_DIR     = PROJECT_DIR / "neaa_data_exports"
EXCLUDE_SITES  = ["USGS-02049500"]

# WRTDS thresholds
MIN_WQ_OBS         = 100    # minimum paired WQ observations
MIN_RECORD_YRS     = 10     # minimum overlapping record length (years)
MAX_GAP_PCT        = 20.0   # maximum % missing daily Q over WQ period
PREFERRED_WQ_OBS   = 200    # preferred threshold for confident WRTDS
PREFERRED_REC_YRS  = 20     # preferred record length

# WQ characteristics to check
WQ_CHARACTERISTICS = ["pH", "Alkalinity", "Temperature, water", "Specific conductance"]

# =============================================================================
# HELPERS
# =============================================================================

def load_sites(path: Path) -> list[str]:
    df = pd.read_excel(path, header=None)
    sites = df.iloc[:, 0].dropna().astype(str).tolist()
    return [s for s in sites if s not in EXCLUDE_SITES]


def get_wq_summary(site: str) -> dict:
    """
    Download WQ data for one site and return summary statistics.
    Uses characteristicName (not pCode) — matches pipeline download.
    """
    try:
        df, _ = wqp.get_results(
            siteid=site,
            characteristicName=WQ_CHARACTERISTICS,
        )
        if df is None or df.empty:
            return {"wq_n": 0, "wq_start": None, "wq_end": None,
                    "wq_span_yr": 0, "wq_chars": ""}

        df = df[df["ActivityMediaName"] == "Water"].copy()
        df["date"] = pd.to_datetime(df["ActivityStartDate"], errors="coerce")
        df = df.dropna(subset=["date"])

        # Count only rows that have a numeric result value
        df["val"] = pd.to_numeric(df["ResultMeasureValue"], errors="coerce")
        df = df.dropna(subset=["val"])

        chars_found = df["CharacteristicName"].unique().tolist()
        n           = len(df)
        start       = df["date"].min()
        end         = df["date"].max()
        span_yr     = (end - start).days / 365.25 if n > 0 else 0

        return {
            "wq_n":       n,
            "wq_start":   start.date() if pd.notna(start) else None,
            "wq_end":     end.date()   if pd.notna(end)   else None,
            "wq_span_yr": round(span_yr, 1),
            "wq_chars":   ", ".join(chars_found),
        }
    except Exception as e:
        return {"wq_n": 0, "wq_start": None, "wq_end": None,
                "wq_span_yr": 0, "wq_chars": f"ERROR: {e}"}


def get_discharge_summary(site: str, wq_start, wq_end) -> dict:
    """
    Download daily discharge from NWIS and return summary statistics
    over the period that overlaps with the water quality record.
    """
    site_no = site.replace("USGS-", "")
    empty = {
        "q_record_start": None, "q_record_end": None,
        "q_record_yrs": 0,
        "q_overlap_start": None, "q_overlap_end": None,
        "q_overlap_yrs": 0,
        "q_gap_pct": 100.0,
        "q_n_days": 0,
        "q_mean_cfs": np.nan,
        "q_median_cfs": np.nan,
        "q_cv": np.nan,
        "q_source": "none",
    }

    try:
        df, _ = nwis.get_dv(
            sites=site_no,
            parameterCd="00060",       # daily mean discharge, ft3/s
            start="1900-01-01",
            end=datetime.today().strftime("%Y-%m-%d"),
        )

        if df is None or df.empty:
            return {**empty, "q_source": "no_nwis_gauge"}

        # Flatten MultiIndex if present
        if isinstance(df.columns, pd.MultiIndex):
            df.columns = ["_".join(c).strip("_") for c in df.columns]

        # Find the discharge column
        q_col = next(
            (c for c in df.columns if "00060" in c and "cd" not in c.lower()),
            None
        )
        if q_col is None:
            return {**empty, "q_source": "no_q_column"}

        df.index = pd.to_datetime(df.index)
        # NWIS returns tz-aware timestamps; WQ dates are tz-naive.
        # Strip timezone so date comparisons work — values are unchanged.
        if hasattr(df.index, 'tz') and df.index.tz is not None:
            df.index = df.index.tz_localize(None)
        df = df[[q_col]].rename(columns={q_col: "Q_cfs"})
        df["Q_cfs"] = pd.to_numeric(df["Q_cfs"], errors="coerce")

        q_start = df.index.min().date()
        q_end   = df.index.max().date()
        q_yrs   = (df.index.max() - df.index.min()).days / 365.25

        # Compute overlap with WQ period
        if wq_start is None or wq_end is None:
            return {**empty, "q_source": "no_wq_dates"}

        overlap_start = max(pd.Timestamp(wq_start), df.index.min())
        overlap_end   = min(pd.Timestamp(wq_end),   df.index.max())

        if overlap_start >= overlap_end:
            return {**empty,
                    "q_record_start": q_start, "q_record_end": q_end,
                    "q_record_yrs": round(q_yrs, 1),
                    "q_source": "no_overlap"}

        # Slice to overlap period
        q_overlap = df.loc[overlap_start:overlap_end, "Q_cfs"]
        overlap_yrs = (overlap_end - overlap_start).days / 365.25

        # Expected days vs actual non-null days
        expected_days  = (overlap_end - overlap_start).days + 1
        actual_days    = q_overlap.notna().sum()
        gap_pct        = (1 - actual_days / expected_days) * 100

        q_valid = q_overlap.dropna()

        return {
            "q_record_start":  q_start,
            "q_record_end":    q_end,
            "q_record_yrs":    round(q_yrs, 1),
            "q_overlap_start": overlap_start.date(),
            "q_overlap_end":   overlap_end.date(),
            "q_overlap_yrs":   round(overlap_yrs, 1),
            "q_gap_pct":       round(gap_pct, 1),
            "q_n_days":        int(actual_days),
            "q_mean_cfs":      round(q_valid.mean(), 2) if len(q_valid) else np.nan,
            "q_median_cfs":    round(q_valid.median(), 2) if len(q_valid) else np.nan,
            "q_cv":            round(q_valid.std() / q_valid.mean(), 3)
                               if len(q_valid) > 1 and q_valid.mean() > 0 else np.nan,
            "q_source":        "nwis_dv",
        }

    except Exception as e:
        return {**empty, "q_source": f"ERROR: {e}"}


def classify_wrtds(row: pd.Series) -> tuple[str, str]:
    """
    Assign a WRTDS tier and plain-English recommendation to each site.

    Tier 1 — Full WRTDS:       meets all preferred thresholds
    Tier 2 — WRTDS possible:   meets minimum thresholds
    Tier 3 — SK with Q covar:  not enough for WRTDS; seasonal Kendall + Q covariate
    Tier 4 — TS only:          minimal data; Theil-Sen only, no flow normalization
    Tier X — Insufficient:     too few observations for any trend analysis
    """
    n        = row["wq_n"]
    span     = row["wq_span_yr"]
    overlap  = row["q_overlap_yrs"]
    gap      = row["q_gap_pct"]
    q_src    = row["q_source"]
    has_q    = q_src == "nwis_dv" and overlap > 0

    # Tier X — can't do any trend analysis
    if n < 25 or span < 5:
        return "X — insufficient data", (
            f"Only {n} WQ obs over {span:.1f} yr. "
            "Need >= 25 obs and >= 5 yr for any trend analysis."
        )

    # Tier 1 — full WRTDS, preferred thresholds
    if (has_q and n >= PREFERRED_WQ_OBS and overlap >= PREFERRED_REC_YRS
            and gap <= MAX_GAP_PCT):
        return "1 — Full WRTDS", (
            f"{n} WQ obs, {overlap:.1f} yr discharge overlap, "
            f"{gap:.1f}% gaps. Ideal for WRTDS."
        )

    # Tier 2 — WRTDS possible but marginal
    if (has_q and n >= MIN_WQ_OBS and overlap >= MIN_RECORD_YRS
            and gap <= MAX_GAP_PCT):
        notes = []
        if n < PREFERRED_WQ_OBS:
            notes.append(f"only {n} WQ obs (preferred >= {PREFERRED_WQ_OBS})")
        if overlap < PREFERRED_REC_YRS:
            notes.append(f"only {overlap:.1f} yr overlap (preferred >= {PREFERRED_REC_YRS})")
        return "2 — WRTDS marginal", (
            f"Meets minimums but {'; '.join(notes)}. "
            "WRTDS feasible, interpret cautiously."
        )

    # Tier 3 — seasonal Kendall with discharge covariate
    if n >= 25 and span >= MIN_RECORD_YRS:
        if not has_q:
            reason = f"no NWIS gauge ({q_src})"
        elif overlap < MIN_RECORD_YRS:
            reason = f"discharge overlap only {overlap:.1f} yr"
        elif gap > MAX_GAP_PCT:
            reason = f"discharge record has {gap:.1f}% gaps"
        else:
            reason = f"only {n} WQ obs"
        return "3 — Seasonal Kendall + Q", (
            f"Can't meet WRTDS minimums ({reason}). "
            "Use seasonal Kendall with discharge covariate instead."
        )

    # Tier 4 — Theil-Sen only
    if n >= 25:
        return "4 — Theil-Sen only", (
            f"{n} WQ obs but only {span:.1f} yr record. "
            "Too short for flow normalization; use Theil-Sen trend only."
        )

    return "X — insufficient data", (
        f"{n} WQ obs. Below minimum threshold for trend analysis."
    )


# =============================================================================
# MAIN
# =============================================================================
def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*65}")
    print("  WRTDS Feasibility Diagnostic")
    print(f"  {datetime.now():%Y-%m-%d %H:%M}")
    print(f"{'='*65}\n")

    # Load site list
    if not SITE_LIST_PATH.exists():
        sys.exit(f"Site list not found: {SITE_LIST_PATH}")
    sites = load_sites(SITE_LIST_PATH)
    print(f"  Checking {len(sites)} sites...\n")

    records = []
    for site in tqdm(sites, desc="Checking sites"):
        wq   = get_wq_summary(site)
        q    = get_discharge_summary(site, wq["wq_start"], wq["wq_end"])
        row  = {"site": site, **wq, **q}
        tier, note = classify_wrtds(pd.Series(row))
        row["wrtds_tier"] = tier
        row["note"]       = note
        records.append(row)

    df = pd.DataFrame(records)

    # ── Sort by tier then WQ obs descending ──────────────────────────────────
    df["tier_sort"] = df["wrtds_tier"].str[0]
    df = df.sort_values(["tier_sort", "wq_n"], ascending=[True, False])
    df = df.drop(columns="tier_sort")

    # ── Save CSV ──────────────────────────────────────────────────────────────
    csv_path = OUTPUT_DIR / "wrtds_feasibility.csv"
    df.to_csv(csv_path, index=False)

    # ── Tier counts ───────────────────────────────────────────────────────────
    tier_counts = df["wrtds_tier"].value_counts().sort_index()

    t1 = df[df["wrtds_tier"].str.startswith("1")]
    t2 = df[df["wrtds_tier"].str.startswith("2")]
    t3 = df[df["wrtds_tier"].str.startswith("3")]
    t4 = df[df["wrtds_tier"].str.startswith("4")]
    tx = df[df["wrtds_tier"].str.startswith("X")]

    # ── Console summary ───────────────────────────────────────────────────────
    print(f"\n{'='*65}")
    print("  RESULTS SUMMARY")
    print(f"{'='*65}")
    print(f"  Total sites checked:          {len(df)}")
    print(f"  Tier 1 — Full WRTDS:          {len(t1)}")
    print(f"  Tier 2 — WRTDS marginal:      {len(t2)}")
    print(f"  Tier 3 — Seasonal Kendall+Q:  {len(t3)}")
    print(f"  Tier 4 — Theil-Sen only:      {len(t4)}")
    print(f"  Tier X — Insufficient:        {len(tx)}")
    print(f"\n  WRTDS feasible (Tier 1+2):    {len(t1)+len(t2)} sites")
    print(f"  Flow-normalized any method:   {len(t1)+len(t2)+len(t3)} sites")

    # ── Print Tier 1 sites ────────────────────────────────────────────────────
    if not t1.empty:
        print(f"\n{'─'*65}")
        print("  TIER 1 — Full WRTDS sites:")
        print(f"{'─'*65}")
        for _, r in t1.iterrows():
            print(f"  {r['site']}")
            print(f"    WQ: {r['wq_n']} obs, {r['wq_start']} to {r['wq_end']}")
            print(f"    Q:  {r['q_overlap_yrs']:.1f} yr overlap, "
                  f"{r['q_gap_pct']:.1f}% gaps, "
                  f"mean {r['q_mean_cfs']:.1f} cfs")

    # ── Print Tier 2 sites ────────────────────────────────────────────────────
    if not t2.empty:
        print(f"\n{'─'*65}")
        print("  TIER 2 — WRTDS marginal sites:")
        print(f"{'─'*65}")
        for _, r in t2.iterrows():
            print(f"  {r['site']}")
            print(f"    WQ: {r['wq_n']} obs, {r['wq_start']} to {r['wq_end']}")
            print(f"    Q:  {r['q_overlap_yrs']:.1f} yr overlap, "
                  f"{r['q_gap_pct']:.1f}% gaps")
            print(f"    Note: {r['note']}")

    # ── Print sites with no NWIS gauge ───────────────────────────────────────
    no_gauge = df[df["q_source"] == "no_nwis_gauge"]
    if not no_gauge.empty:
        print(f"\n{'─'*65}")
        print(f"  Sites with no NWIS discharge gauge ({len(no_gauge)}):")
        print(f"{'─'*65}")
        for _, r in no_gauge.iterrows():
            print(f"  {r['site']}  ({r['wq_n']} WQ obs)")

    # ── Write text report ─────────────────────────────────────────────────────
    txt_path = OUTPUT_DIR / "wrtds_feasibility.txt"
    with open(txt_path, "w") as f:
        f.write(f"WRTDS Feasibility Report\n")
        f.write(f"Generated: {datetime.now():%Y-%m-%d %H:%M}\n")
        f.write(f"Sites checked: {len(df)}\n\n")

        for tier_label, subset in [
            ("Tier 1 — Full WRTDS",         t1),
            ("Tier 2 — WRTDS marginal",     t2),
            ("Tier 3 — Seasonal Kendall+Q", t3),
            ("Tier 4 — Theil-Sen only",     t4),
            ("Tier X — Insufficient data",  tx),
        ]:
            f.write(f"\n{'='*60}\n{tier_label} ({len(subset)} sites)\n{'='*60}\n")
            for _, r in subset.iterrows():
                f.write(f"\n  {r['site']}\n")
                f.write(f"    WQ observations:  {r['wq_n']}\n")
                f.write(f"    WQ period:        {r['wq_start']} to {r['wq_end']}"
                        f" ({r['wq_span_yr']} yr)\n")
                f.write(f"    Characteristics:  {r['wq_chars']}\n")
                f.write(f"    Discharge source: {r['q_source']}\n")
                if r["q_source"] == "nwis_dv":
                    f.write(f"    Q record:         {r['q_record_start']} to "
                            f"{r['q_record_end']} ({r['q_record_yrs']} yr)\n")
                    f.write(f"    Q overlap:        {r['q_overlap_yrs']} yr, "
                            f"{r['q_gap_pct']}% gaps\n")
                    f.write(f"    Q stats:          mean {r['q_mean_cfs']} cfs, "
                            f"median {r['q_median_cfs']} cfs, "
                            f"CV {r['q_cv']}\n")
                f.write(f"    Note: {r['note']}\n")

    print(f"\n{'='*65}")
    print(f"  Full results saved to:")
    print(f"    {csv_path}")
    print(f"    {txt_path}")
    print(f"{'='*65}\n")


if __name__ == "__main__":
    main()
