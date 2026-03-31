"""
NEAA (National Estuary Acidification Assessment) Pipeline
==========================================================
Python translation of R scripts:
  neaa_final_download_231028.R  -> Step 1: Download & QA/QC river data
  neaa_corrections.R            -> Step 2: DOC/discharge corrections (Liu et al. 2020)
  ts_analysis_240407_seasonal.R -> Step 3: Theil-Sen + Seasonal Kendall trends
  neaa_calcium.R                -> Step 4: Download & join calcium data
  codap_load.R                  -> Step 5: Load & filter CODAP ocean data
  codap_organize.R              -> Step 6: Match river sites to ocean stations
  export_neaa_csv_calcium.R /
  export_ocean_csv.R            -> Step 7: Export CSVs + HDF5 for MATLAB

Key improvements over R version:
  - Config file replaces all hardcoded paths
  - Vectorized operations replace row-by-row loops
  - PyCO2SYS replaces seacarb for carbonate chemistry
  - dataretrieval (USGS) replaces dataRetrieval (R)
  - Proper logging throughout
  - Single HDF5 export alongside per-site CSVs
  - CODAP downloaded from URL instead of local Excel

Dependencies (install with pip):
  pip install dataretrieval pandas numpy scipy pyco2sys openpyxl
              geopy h5py requests tqdm pymannkendall theilsen

Author: Translated from R (Spacella et al.) by Claude
"""

# ── Standard library ─────────────────────────────────────────────────────────
import logging
import warnings
from pathlib import Path
from datetime import datetime

# ── Third-party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
import requests
from scipy import stats
from geopy.distance import great_circle
import h5py
from tqdm import tqdm

# Conditional imports with helpful error messages
try:
    import dataretrieval.nwis as nwis
    import dataretrieval.wqp as wqp
except ImportError:
    raise ImportError("Install dataretrieval:  pip install dataretrieval")

try:
    import PyCO2SYS as co2sys
except ImportError:
    raise ImportError("Install PyCO2SYS:  pip install PyCO2SYS")

try:
    import pymannkendall as mk
except ImportError:
    raise ImportError("Install pymannkendall:  pip install pymannkendall")

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=pd.errors.DtypeWarning)

# =============================================================================
# CONFIGURATION  — edit these paths and thresholds for your project
# =============================================================================
CFG = {
    # ── Paths ─────────────────────────────────────────────────────────────────
    "project_dir":   Path.home() / "NEAA_Analysis",          # root folder
    "site_list_xlsx": "WQP_sites_manual_srp_update6.xlsx",   # site list file
    "export_dir":    "neaa_data_exports",                     # output subfolder
    "codap_cache":   "CODAP_NA_v2020.parquet",                # cached CODAP data

    # ── CODAP download URL (NOAA/NCEI)  ───────────────────────────────────────
    # If this URL changes, update it here. The file is ~30 MB.
    "codap_url": (
        "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0219960/"
        "CODAP_NA_v2020.xlsx"
    ),

    # ── WQP characteristic names to download ─────────────────────────────────
    "wqp_characteristics": [
        "Alkalinity",
        "pH",
        "Temperature, water",
        "Specific conductance",
    ],

    # ── QA/QC thresholds ─────────────────────────────────────────────────────
    "ph_min":           5.5,
    "ph_max":           9.5,
    "temp_min_c":       0.0,
    "temp_max_c":       50.0,
    "cond_min":         0.0,
    "cond_max":         5000.0,   # µS/cm
    "alk_multiplier":   20,        # unit conversion factor (matches R script)
    "alk_min":          50.0,      # µeq/kg after conversion
    "alk_max":          9000.0,
    # ── Freshwater filter ────────────────────────────────────────────────────
    # We filter marine-influenced samples using specific conductance directly
    # rather than converting to salinity. This avoids the ionic-composition
    # assumption that makes the Hill et al./UNESCO formula invalid for rivers.
    #
    # Thresholds (tiered, applied in order):
    #   1. Conductance < cond_freshwater_max  → freshwater (keep)
    #   2. Conductance > cond_brackish_min    → brackish/marine (drop)
    #   3. In between → check chloride if available, else flag as ambiguous
    #
    # cond_freshwater_max = 1500 µS/cm
    #   Conservative upper bound for freshwater in coastal rivers.
    #   Equivalent to ~0.8 PSU in typical ionic composition.
    #   USGS uses 1000–2000 µS/cm as freshwater/estuarine boundary.
    #
    # cond_brackish_min = 10000 µS/cm
    #   Above this, marine influence is unambiguous regardless of ion type.
    #   ~6 PSU equivalent. Samples in the 1500–10000 µS/cm range get
    #   chloride-checked if available, otherwise flagged.
    #
    # cl_freshwater_max = 100 mg/L
    #   Chloride threshold for freshwater classification.
    #   Seawater contains ~19,000 mg/L Cl. Marine-influenced samples
    #   show elevated Cl relative to local freshwater background.
    #   100 mg/L is conservative; many studies use 250 mg/L (EPA limit).
    "cond_freshwater_max":  1500,   # µS/cm — definitely freshwater below this
    "cond_brackish_min":   10000,   # µS/cm — definitely brackish above this
    "cl_freshwater_max":     100,   # mg/L  — chloride freshwater threshold
    "cl_pcode":           "00940",  # WQP pCode for chloride
    "outlier_iqr_mult": 3.0,       # multiplier for IQR-based outlier removal

    # ── Site eligibility filters ──────────────────────────────────────────────
    "min_obs_any":      5,         # minimum observations post-1950
    "min_obs_recent":   10,        # minimum observations post-1990
    "min_obs_2000":     10,        # minimum observations post-2000
    "min_obs_trend":    25,        # minimum observations for trend analysis
    "min_duration_yr":  10,        # minimum record length for trend analysis
    "date_floor":       "1950-01-01",
    "date_recent":      "1990-01-01",
    "date_2000":        "2000-01-01",

    # ── Calcium parameter codes ───────────────────────────────────────────────
    "calcium_pcodes": ["00915"],   # filtered calcium mg/L
    "doc_pcode":      ["00681"],   # DOC mg/L
    # Discharge is now downloaded from NWIS daily values (get_dv),
    # not WQP. Parameter code 00060 = daily mean discharge ft3/s.
    "discharge_pcode": "00060",

    # ── CODAP ocean matching ──────────────────────────────────────────────────
    "codap_depth_max_m":  20,      # surface-only cutoff
    "ocean_radius_m":     200_000, # search radius for ocean stations (200 km)
    "codap_bad_flags":    [3, 9],  # quality flags to remove

    # ── Estuarine mixing ─────────────────────────────────────────────────────
    "f_marine_steps": 21,          # number of mixing fractions (0→1)

    # ── Excluded sites ────────────────────────────────────────────────────────
    "exclude_sites": ["USGS-02049500"],
}

# =============================================================================
# LOGGING SETUP
# =============================================================================
def setup_logging(log_dir: Path) -> logging.Logger:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"neaa_run_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger("neaa")


# =============================================================================
# STEP 0 — DIRECTORY SETUP
# =============================================================================
def setup_dirs(cfg: dict) -> tuple[Path, Path]:
    project = Path(cfg["project_dir"])
    export  = project / cfg["export_dir"]
    export.mkdir(parents=True, exist_ok=True)
    return project, export


# =============================================================================
# STEP 1 — DOWNLOAD & QA/QC RIVER DATA
# Equivalent: neaa_final_download_231028.R
# =============================================================================

def load_site_list(project_dir: Path, cfg: dict, log: logging.Logger) -> list[str]:
    """Read the Excel site list and return final site IDs."""
    path = project_dir / cfg["site_list_xlsx"]
    if not path.exists():
        raise FileNotFoundError(
            f"Site list not found: {path}\n"
            "Place WQP_sites_manual_srp_update6.xlsx in your project directory."
        )
    df = pd.read_excel(path, header=None)
    sites = df.iloc[:, 0].dropna().astype(str).tolist()
    # Apply exclusions
    sites = [s for s in sites if s not in cfg["exclude_sites"]]
    log.info(f"Loaded {len(sites)} sites from site list.")
    return sites


def download_wqp_data(sites: list[str], cfg: dict, log: logging.Logger) -> pd.DataFrame:
    """
    Download pH, alkalinity, temperature, conductivity from WQP.
    Replaces the for-loop over readWQPdata() in R.
    """
    log.info("Downloading WQP data for all sites...")
    frames = []
    for site in tqdm(sites, desc="WQP download"):
        try:
            df, _ = wqp.get_results(
                siteid=site,
                characteristicName=cfg["wqp_characteristics"],
            )
            if df is not None and not df.empty:
                frames.append(df)
        except Exception as e:
            log.warning(f"  WQP download failed for {site}: {e}")

    if not frames:
        raise RuntimeError("No WQP data downloaded for any site.")

    all_data = pd.concat(frames, ignore_index=True)
    log.info(f"  Downloaded {len(all_data):,} total rows across {len(sites)} sites.")

    # Keep only water media
    all_data = all_data[all_data["ActivityMediaName"] == "Water"].copy()

    # Create unique sample ID (date + location, matching R's paste())
    all_data["ID"] = (
        all_data["ActivityStartDate"].astype(str)
        + " "
        + all_data["MonitoringLocationIdentifier"]
    )
    all_data = all_data.drop_duplicates()
    all_data["Datetime"] = pd.to_datetime(
        all_data["ActivityStartDate"], errors="coerce"
    )
    log.info(f"  {len(all_data):,} rows after deduplication & media filter.")
    return all_data


def _keep_first_time(df: pd.DataFrame) -> pd.DataFrame:
    """For duplicate IDs, keep the earliest time-stamped record.
    Replaces: df %>% group_by(ID) %>% filter(ActivityStartTime.Time == min(...))
    """
    if "ActivityStartTime/Time" in df.columns:
        df = (
            df.sort_values("ActivityStartTime/Time")
            .groupby("ID", sort=False)
            .first()
            .reset_index()
        )
    return df.drop_duplicates(subset="ID")


def _prefer_field_measurement(df: pd.DataFrame) -> pd.DataFrame:
    """For duplicate IDs, prefer the record with the lowest USGSPCode
    (field measurements < lab measurements in USGS coding).
    Replaces: group_by(ID) %>% filter(USGSPCode == min(USGSPCode))
    """
    if "USGSPCode" in df.columns:
        df["USGSPCode"] = pd.to_numeric(df["USGSPCode"], errors="coerce")
        df = (
            df.sort_values("USGSPCode")
            .groupby("ID", sort=False)
            .first()
            .reset_index()
        )
    return df


def extract_parameter(
    all_data: pd.DataFrame,
    char_name: str,
    value_min: float,
    value_max: float,
    prefer_field: bool = False,
    log: logging.Logger | None = None,
) -> pd.DataFrame:
    """
    Filter to one characteristic, coerce values to numeric,
    apply range QA/QC, remove NAs, deduplicate.
    Vectorized replacement of the per-parameter blocks in R.
    """
    raw = all_data[all_data["CharacteristicName"] == char_name].copy()
    n_raw = len(raw)
    n_sites_raw = raw["MonitoringLocationIdentifier"].nunique()

    raw["ResultMeasureValue"] = pd.to_numeric(
        raw["ResultMeasureValue"], errors="coerce"
    )
    n_non_numeric = raw["ResultMeasureValue"].isna().sum()

    # Check units before applying range filter — helps diagnose unit mismatches
    if log and "ResultMeasure/MeasureUnitCode" in raw.columns:
        units_found = raw["ResultMeasure/MeasureUnitCode"].value_counts().head(5)
        log.info(f"    {char_name} units: {units_found.to_dict()}")

    raw.loc[
        (raw["ResultMeasureValue"] < value_min)
        | (raw["ResultMeasureValue"] > value_max),
        "ResultMeasureValue",
    ] = np.nan
    n_out_of_range = raw["ResultMeasureValue"].isna().sum() - n_non_numeric

    df = raw.dropna(subset=["ResultMeasureValue"])
    df = df.sort_values(["MonitoringLocationIdentifier", "ActivityStartDate"])
    if prefer_field:
        df = _prefer_field_measurement(df)
    df = _keep_first_time(df)
    df = df.drop_duplicates(subset="ID")

    n_sites_final = df["MonitoringLocationIdentifier"].nunique()

    if log:
        log.info(
            f"  {char_name}: {n_raw:,} raw rows ({n_sites_raw} sites) → "
            f"{len(df):,} kept ({n_sites_final} sites) | "
            f"non-numeric: {n_non_numeric}, out-of-range [{value_min},{value_max}]: {n_out_of_range}"
        )
    return df


def merge_parameters(
    pH_df: pd.DataFrame,
    temp_df: pd.DataFrame,
    alk_df: pd.DataFrame,
    cond_df: pd.DataFrame,
    log: logging.Logger,
) -> pd.DataFrame:
    """
    Inner-join all four parameters on ID.
    Replaces the chain of inner_join() calls in R.
    Returns a tidy, renamed DataFrame.
    """
    # Select only the columns we need from each frame before joining
    def _slim(df: pd.DataFrame, suffix: str) -> pd.DataFrame:
        cols = {
            "ID": "ID",
            "OrganizationIdentifier": "OrganizationIdentifier",
            "ActivityStartDate": "ActivityStartDate",
            "MonitoringLocationIdentifier": "MonitoringLocation",
            "ResultMeasureValue": f"{suffix}_value",
            "USGSPCode": f"{suffix}_pcode",
        }
        available = {k: v for k, v in cols.items() if k in df.columns}
        return df[list(available.keys())].rename(columns=available)

    # Log site counts going into each join — the key diagnostic
    log.info(f"  Pre-join site counts:")
    log.info(f"    pH:          {pH_df['MonitoringLocationIdentifier'].nunique()} sites, {len(pH_df):,} rows")
    log.info(f"    Temperature: {temp_df['MonitoringLocationIdentifier'].nunique()} sites, {len(temp_df):,} rows")
    log.info(f"    Alkalinity:  {alk_df['MonitoringLocationIdentifier'].nunique()} sites, {len(alk_df):,} rows")
    log.info(f"    Conductance: {cond_df['MonitoringLocationIdentifier'].nunique()} sites, {len(cond_df):,} rows")

    j1 = _slim(pH_df,   "pH").merge(_slim(temp_df, "Temp"), on="ID", suffixes=("", "_temp"))
    log.info(f"  After pH ∩ Temp:              {j1['MonitoringLocation'].nunique()} sites, {len(j1):,} rows")

    j2 = j1.merge(_slim(alk_df,  "Alk"),  on="ID", suffixes=("", "_alk"))
    log.info(f"  After (pH ∩ Temp) ∩ Alk:     {j2['MonitoringLocation'].nunique()} sites, {len(j2):,} rows")

    j3 = j2.merge(_slim(cond_df, "Cond"), on="ID", suffixes=("", "_cond"))
    log.info(f"  After all four parameters:    {j3['MonitoringLocation'].nunique()} sites, {len(j3):,} rows")

    # Warn if the join is very destructive
    n_before = pH_df["MonitoringLocationIdentifier"].nunique()
    n_after  = j3["MonitoringLocation"].nunique() if not j3.empty else 0
    pct_lost = 100 * (1 - n_after / n_before) if n_before > 0 else 0
    if pct_lost > 50:
        log.warning(
            f"  *** {pct_lost:.0f}% of sites lost in 4-way inner join! ***"
            "  Most likely cause: alkalinity and pH sampled on different dates."
            "  Consider relaxing the join to allow date-proximate matching."
        )

    # Unify duplicated metadata columns from joins
    for col in ["OrganizationIdentifier", "ActivityStartDate", "MonitoringLocation"]:
        dupes = [c for c in j3.columns if c.startswith(col) and c != col]
        for d in dupes:
            j3.drop(columns=d, inplace=True, errors="ignore")

    j3 = j3[j3["ActivityStartDate"] > "1950-01-01"].copy()
    log.info(f"  {len(j3):,} records after merging four parameters.")
    return j3


def remove_outliers_iqr(df: pd.DataFrame, cols: list[str], mult: float) -> pd.DataFrame:
    """
    Remove rows where any of `cols` falls outside Q1 - mult*IQR or Q3 + mult*IQR.
    Vectorized replacement of the per-site IQR loop in R.
    """
    mask = pd.Series(True, index=df.index)
    for col in cols:
        if col not in df.columns:
            continue
        q1 = df[col].quantile(0.25)
        q3 = df[col].quantile(0.75)
        iqr = q3 - q1
        mask &= (df[col] >= q1 - mult * iqr) & (df[col] <= q3 + mult * iqr)
    return df[mask].copy()


def classify_freshwater(
    cond_uScm: pd.Series,
    cl_mg_L:   pd.Series | None,
    cfg:       dict,
) -> pd.DataFrame:
    """
    Classify each observation as freshwater, brackish, or ambiguous using
    a tiered conductance + chloride approach.

    This replaces the Hill et al./UNESCO salinity conversion, which is
    invalid for rivers because it assumes seawater ionic composition.
    Rivers have highly variable ionic ratios between catchments, making
    any conductance→salinity formula unreliable for freshwater classification.

    Tier 1 — conductance below cond_freshwater_max (1500 µS/cm):
      Always freshwater regardless of ionic composition. Keep.

    Tier 2 — conductance above cond_brackish_min (10000 µS/cm):
      Always brackish/marine regardless of ionic composition. Drop.

    Tier 3 — conductance between thresholds (ambiguous zone):
      Check chloride if available:
        Cl < cl_freshwater_max (100 mg/L) → freshwater, keep
        Cl >= cl_freshwater_max           → brackish, drop
        No Cl data                        → flag as ambiguous, keep with warning

    Returns a DataFrame with columns:
      is_freshwater : bool   — True = keep, False = drop
      fw_method     : str    — which tier/method was used (for QA documentation)
      cond_uScm     : float  — conductance value used
      cl_mg_L       : float  — chloride value used (NaN if not available)

    Scientific basis:
      USGS uses 1000–2000 µS/cm as the freshwater/estuarine boundary in
      most monitoring programs. Chloride threshold of 100 mg/L is conservative
      relative to the EPA secondary drinking water standard (250 mg/L) and well
      below any meaningful marine influence (~19,000 mg/L in seawater).
    """
    cond = pd.to_numeric(cond_uScm, errors="coerce")
    cl   = pd.to_numeric(cl_mg_L,   errors="coerce") if cl_mg_L is not None            else pd.Series(np.nan, index=cond.index)

    fw_max   = cfg["cond_freshwater_max"]   # 1500 µS/cm
    br_min   = cfg["cond_brackish_min"]     # 10000 µS/cm
    cl_max   = cfg["cl_freshwater_max"]     # 100 mg/L

    n        = len(cond)
    is_fw    = pd.Series(False, index=cond.index)
    method   = pd.Series("unknown", index=cond.index)

    # Tier 1: definitely freshwater
    t1 = cond < fw_max
    is_fw[t1]  = True
    method[t1] = f"cond<{fw_max}"

    # Tier 2: definitely brackish
    t2 = cond >= br_min
    is_fw[t2]  = False
    method[t2] = f"cond>={br_min}"

    # Tier 3: ambiguous zone — use chloride
    t3 = (~t1) & (~t2) & cond.notna()

    # Sub-tier 3a: chloride available
    t3a = t3 & cl.notna()
    is_fw[t3a]                   = cl[t3a] < cl_max
    method[t3a & (cl < cl_max)]  = f"cl<{cl_max}"
    method[t3a & (cl >= cl_max)] = f"cl>={cl_max}"

    # Sub-tier 3b: no chloride — flag but keep (conservative)
    t3b = t3 & cl.isna()
    is_fw[t3b]  = True   # keep ambiguous samples, flag them
    method[t3b] = "ambiguous_no_cl"

    # Missing conductance — cannot classify, keep with flag
    t_na = cond.isna()
    is_fw[t_na]  = True
    method[t_na] = "no_cond_data"

    return pd.DataFrame({
        "is_freshwater": is_fw,
        "fw_method":     method,
        "cond_uScm":     cond,
        "cl_mg_L":       cl,
    })


def build_site_datasets(
    merged: pd.DataFrame,
    sites: list[str],
    cfg: dict,
    log: logging.Logger,
) -> dict[str, pd.DataFrame]:
    """
    Per-site QA/QC: outlier removal, salinity filter, obs-count filter.
    Returns dict of site_id -> cleaned DataFrame.
    Replaces the two nested for-loops in R.
    """
    log.info("Building per-site datasets with QA/QC...")
    result = {}

    # Filter post-date-floor first (fast, whole-frame)
    merged = merged[merged["ActivityStartDate"] > cfg["date_floor"]].copy()

    n_in_merged  = merged["MonitoringLocation"].nunique()
    n_empty      = 0
    n_iqr_drop   = 0
    n_sal_drop   = 0
    n_obs_drop   = 0
    sal_too_high = []

    for site in sites:
        site_df = merged[merged["MonitoringLocation"] == site].copy()
        if site_df.empty:
            n_empty += 1
            continue

        # IQR outlier removal
        site_df = remove_outliers_iqr(
            site_df,
            cols=["Alk_value", "pH_value", "Temp_value"],
            mult=cfg["outlier_iqr_mult"],
        )
        if site_df.empty:
            n_iqr_drop += 1
            continue

        # ── Freshwater classification ─────────────────────────────────────────
        # Use conductance threshold + chloride rather than salinity conversion.
        # Chloride data (if available) is passed in via chloride_data dict.
        cl_series = None
        if hasattr(cfg, "_chloride_data") and site in cfg["_chloride_data"]:
            # Merge chloride onto site_df by date
            cl_df = cfg["_chloride_data"][site]
            site_df["_date_only"] = pd.to_datetime(
                site_df["ActivityStartDate"], errors="coerce"
            ).dt.normalize()
            merged_cl = site_df.merge(
                cl_df.rename(columns={"ActivityStartDate": "_date_only"}),
                on="_date_only", how="left"
            )
            cl_series = merged_cl["chloride_mg_L"] if "chloride_mg_L" in merged_cl else None
            site_df   = site_df.drop(columns=["_date_only"])

        fw_class = classify_freshwater(
            cond_uScm = site_df["Cond_value"],
            cl_mg_L   = cl_series,
            cfg       = cfg,
        )

        # Attach classification columns to site_df for export
        site_df["fw_method"]    = fw_class["fw_method"].values
        site_df["is_freshwater"] = fw_class["is_freshwater"].values

        # Count observation-level drops (not whole-site drops)
        n_marine_obs = (~fw_class["is_freshwater"]).sum()
        n_ambig_obs  = (fw_class["fw_method"] == "ambiguous_no_cl").sum()

        site_df = site_df[site_df["is_freshwater"]].copy()

        if site_df.empty:
            n_sal_drop += 1
            mean_cond = pd.to_numeric(
                merged[merged["MonitoringLocation"] == site]["Cond_value"],
                errors="coerce"
            ).mean()
            sal_too_high.append(
                f"{site} (mean cond={mean_cond:.0f} µS/cm, "
                f"all {n_marine_obs} obs marine-influenced)"
            )
            continue

        if n_ambig_obs > 0:
            log.debug(
                f"  {site}: {n_ambig_obs} observations in ambiguous conductance "
                f"zone kept without chloride confirmation — review manually"
            )

        if len(site_df) >= cfg["min_obs_any"]:
            result[site] = site_df
        else:
            n_obs_drop += 1

    log.info(f"  QA/QC site funnel ({n_in_merged} sites entering):")
    log.info(f"    Not in merged (no same-day co-samples): {n_empty}")
    log.info(f"    Lost to IQR outlier removal:             {n_iqr_drop}")
    log.info(f"    Lost to freshwater filter:               {n_sal_drop}")
    if sal_too_high:
        for s in sal_too_high[:5]:
            log.info(f"      {s}")
        if len(sal_too_high) > 5:
            log.info(f"      ...and {len(sal_too_high)-5} more")
    log.info(f"    Lost to min obs filter (<{cfg['min_obs_any']}):         {n_obs_drop}")
    log.info(f"    Passed QA/QC:                            {len(result)}")
    return result


def get_site_subsets(
    all_data_site: dict[str, pd.DataFrame],
    cfg: dict,
    log: logging.Logger,
) -> tuple[dict, dict, dict]:
    """
    Create the three temporal subsets used throughout the R pipeline:
    - all_data_site_final  : post-1950,  >=5 obs
    - all_data_site_recent : post-1990, >=10 obs
    - all_data_site_2000   : post-2000, >=10 obs
    """
    def _subset(cutoff: str, min_obs: int) -> dict[str, pd.DataFrame]:
        out = {}
        for site, df in all_data_site.items():
            sub = df[df["ActivityStartDate"] > cutoff]
            if len(sub) >= min_obs:
                out[site] = sub
        return out

    final  = _subset(cfg["date_floor"],  cfg["min_obs_any"])
    recent = _subset(cfg["date_recent"], cfg["min_obs_recent"])
    y2000  = _subset(cfg["date_2000"],   cfg["min_obs_2000"])

    log.info(
        f"Site counts — final: {len(final)}, "
        f"post-1990: {len(recent)}, post-2000: {len(y2000)}"
    )
    return final, recent, y2000


def fetch_site_metadata(sites: list[str], log: logging.Logger) -> pd.DataFrame:
    """Fetch NWIS site metadata. Replaces readNWISsite()."""
    site_numbers = [s.replace("USGS-", "") for s in sites]
    try:
        meta, _ = nwis.get_info(sites=site_numbers)
        log.info(f"  Fetched metadata for {len(meta)} sites.")
        return meta
    except Exception as e:
        log.warning(f"  Site metadata fetch failed: {e}")
        return pd.DataFrame()


# =============================================================================
# STEP 2 — CORRECTIONS (Liu et al. 2020)
# Equivalent: neaa_corrections.R
# =============================================================================

def download_doc(
    sites: list[str], cfg: dict, log: logging.Logger
) -> dict[str, pd.DataFrame]:
    """Download DOC data. Replaces readWQPdata(pCode='00681') loop."""
    log.info("Downloading DOC data...")
    result = {}
    for site in tqdm(sites, desc="DOC download"):
        try:
            # dataretrieval >= 1.0 uses param= for pCode-based WQP queries
            df, _ = wqp.get_results(siteid=site, param=cfg["doc_pcode"])
            if df is not None and not df.empty:
                df = df.dropna(subset=["ResultMeasureValue"])
                df["ResultMeasureValue"] = pd.to_numeric(
                    df["ResultMeasureValue"], errors="coerce"
                )
                # Convert mg/L C -> umol/L -> organic alkalinity estimate
                # Following Liu et al. 2020: org_alk = DOC_umol * 0.8 * 0.1
                df["doc_umol"] = df["ResultMeasureValue"] / 1000 / 12.011 * 1e6
                df["org_alk"]  = df["doc_umol"] * 0.8 * 0.1
                if len(df) > 5:
                    result[site] = df
        except Exception as e:
            log.debug(f"  DOC download failed for {site}: {e}")
    log.info(f"  DOC data available for {len(result)} sites.")
    return result


def download_chloride(
    sites: list[str], cfg: dict, log: logging.Logger
) -> dict[str, pd.DataFrame]:
    """
    Download chloride data for sites in the ambiguous conductance zone
    (1500–10000 µS/cm). Used as a secondary freshwater classifier.

    pCode 00940 = chloride, water, filtered, mg/L.
    Only downloaded for sites where conductance data suggests potential
    marine influence — most inland sites won't need this.

    Returns dict of site_id -> DataFrame with ActivityStartDate and
    chloride_mg_L columns. Empty dict if no data found.
    """
    log.info("Downloading chloride data for marine influence check...")
    result = {}
    for site in tqdm(sites, desc="Chloride download"):
        try:
            df, _ = wqp.get_results(siteid=site, param=[cfg["cl_pcode"]])
            if df is not None and not df.empty:
                df["chloride_mg_L"] = pd.to_numeric(
                    df["ResultMeasureValue"], errors="coerce"
                )
                df = df.dropna(subset=["chloride_mg_L"])
                df = df[df["chloride_mg_L"] >= 0]
                if not df.empty:
                    # Keep date and value only — will be merged by date
                    out = df[["ActivityStartDate", "chloride_mg_L"]].copy()
                    out["ActivityStartDate"] = pd.to_datetime(
                        out["ActivityStartDate"], errors="coerce"
                    )
                    out = out.dropna(subset=["ActivityStartDate"])
                    # Average if multiple measurements on same date
                    out = out.groupby("ActivityStartDate")["chloride_mg_L"]                              .mean().reset_index()
                    result[site] = out
        except Exception as e:
            log.debug(f"  Chloride download failed for {site}: {e}")
    log.info(f"  Chloride data available for {len(result)} sites.")
    return result


def download_discharge(
    sites: list[str], cfg: dict, log: logging.Logger
) -> dict[str, pd.DataFrame]:
    """
    Download COMPLETE daily discharge records from NWIS for each site.

    Replaces the old WQP grab-sample discharge download, which only returned
    discharge values coinciding with water quality samples — typically a few
    hundred values per site at most. NWIS daily values give the full continuous
    record (often 10,000-30,000 daily values per site), which is required for:
      - Reliable flow statistics (mean, median, flow percentiles)
      - Seasonal Kendall with discharge covariate
      - WRTDS flow normalization

    Parameter code 00060 = daily mean discharge in ft³/s (cfs).
    Both ft³/s and m³/s are stored; m³/s is used internally for consistency
    with SI units. 1 ft³/s = 0.0283168 m³/s.

    Returns dict of site_id -> DataFrame with columns:
      date        : datetime index (daily)
      Q_cfs       : daily mean discharge, ft³/s
      Q_cms       : daily mean discharge, m³/s
      Q_cd        : NWIS quality code for the daily value
    """
    log.info("Downloading NWIS daily discharge records...")
    result = {}
    CFS_TO_CMS = 0.0283168

    for site in tqdm(sites, desc="Discharge (NWIS daily)"):
        site_no = site.replace("USGS-", "")
        try:
            df, _ = nwis.get_dv(
                sites=site_no,
                parameterCd="00060",
                start="1900-01-01",
                end=datetime.today().strftime("%Y-%m-%d"),
            )

            if df is None or df.empty:
                log.debug(f"  No NWIS daily discharge for {site}")
                continue

            # Flatten MultiIndex columns produced by dataretrieval
            if isinstance(df.columns, pd.MultiIndex):
                df.columns = ["_".join(str(c) for c in col).strip("_")
                              for col in df.columns]

            # Locate discharge and quality-code columns
            q_col  = next((c for c in df.columns
                           if "00060" in c and "cd" not in c.lower()), None)
            qc_col = next((c for c in df.columns
                           if "00060" in c and "cd" in c.lower()), None)

            if q_col is None:
                log.debug(f"  No 00060 column found for {site}")
                continue

            df.index = pd.to_datetime(df.index)
            # Strip timezone from NWIS timestamps for consistent comparisons
            if hasattr(df.index, 'tz') and df.index.tz is not None:
                df.index = df.index.tz_localize(None)
            out = pd.DataFrame(index=df.index)
            out.index.name = "date"
            out["Q_cfs"] = pd.to_numeric(df[q_col], errors="coerce")
            out["Q_cms"] = out["Q_cfs"] * CFS_TO_CMS
            out["Q_cd"]  = df[qc_col] if qc_col else np.nan

            # Remove negative or zero discharge (unphysical)
            out = out[out["Q_cfs"] > 0].copy()

            if out.empty:
                continue

            result[site] = out

        except Exception as e:
            log.debug(f"  NWIS discharge failed for {site}: {e}")

    log.info(
        f"  Daily discharge retrieved for {len(result)} sites "
        f"(avg {np.mean([len(v) for v in result.values()]):.0f} days/site)"
        if result else "  No discharge data retrieved."
    )
    return result


def compute_discharge_stats(
    discharge_data: dict[str, pd.DataFrame],
    log: logging.Logger,
) -> pd.DataFrame:
    """
    Compute flow statistics for each site from the full daily record.

    These statistics serve two purposes:
      1. Site characterisation (mean flow, variability) for reporting
      2. Input to flow-normalised trend analysis

    Statistics returned (all in m³/s unless noted):
      Q_mean        : arithmetic mean discharge
      Q_median      : median discharge (Q50)
      Q_cv          : coefficient of variation (std/mean) — flow flashiness
      Q_q10         : 10th percentile — low-flow condition
      Q_q25         : 25th percentile
      Q_q75         : 75th percentile
      Q_q90         : 90th percentile — high-flow condition
      Q_mean_cfs    : mean discharge in ft³/s (for USGS reporting compatibility)
      Q_record_yrs  : length of discharge record in years
      Q_n_days      : number of daily values in record
      Q_gap_pct     : percentage of days missing over the record span
    """
    log.info("Computing discharge statistics...")
    rows = []
    for site, df in discharge_data.items():
        q = df["Q_cms"].dropna()
        if q.empty:
            continue

        span_days    = (df.index.max() - df.index.min()).days + 1
        gap_pct      = (1 - len(q) / span_days) * 100 if span_days > 0 else 100

        rows.append({
            "site":          site,
            "Q_mean":        q.mean(),
            "Q_median":      q.median(),
            "Q_cv":          q.std() / q.mean() if q.mean() > 0 else np.nan,
            "Q_q10":         q.quantile(0.10),
            "Q_q25":         q.quantile(0.25),
            "Q_q75":         q.quantile(0.75),
            "Q_q90":         q.quantile(0.90),
            "Q_mean_cfs":    df["Q_cfs"].dropna().mean(),
            "Q_record_yrs":  (df.index.max() - df.index.min()).days / 365.25,
            "Q_n_days":      len(q),
            "Q_gap_pct":     round(gap_pct, 1),
        })

    stats = pd.DataFrame(rows)
    log.info(f"  Flow statistics computed for {len(stats)} sites.")
    return stats


def match_discharge_to_wq(
    all_data_site: dict[str, pd.DataFrame],
    discharge_data: dict[str, pd.DataFrame],
    log: logging.Logger,
) -> dict[str, pd.DataFrame]:
    """
    Match each water quality observation to the concurrent daily discharge
    value from NWIS. This is the key step that enables both:
      - Seasonal Kendall with discharge covariate
      - WRTDS flow normalization

    For each WQ sample date, we look up the daily mean discharge from the
    NWIS record. If the exact date is missing (gap in discharge record), we
    use linear interpolation across gaps <= 7 days, consistent with EGRET's
    approach. Gaps > 7 days are left as NaN and flagged.

    Adds columns to each site DataFrame:
      Q_cms          : concurrent daily discharge, m³/s
      Q_cfs          : concurrent daily discharge, ft³/s
      log_Q          : natural log of Q_cms — used as covariate in SK test
                       and as x-axis in WRTDS rating curve
      Q_matched      : True if exact date match; False if interpolated
      Q_gap_days     : number of days interpolated across (0 if exact match)
    """
    log.info("Matching discharge to water quality observations...")
    result    = {}
    n_matched = 0
    n_interp  = 0
    n_missing = 0

    for site, wq_df in all_data_site.items():
        wq_out = wq_df.copy()
        wq_out["Q_cms"]      = np.nan
        wq_out["Q_cfs"]      = np.nan
        wq_out["log_Q"]      = np.nan
        wq_out["Q_matched"]  = False
        wq_out["Q_gap_days"] = np.nan

        if site not in discharge_data or discharge_data[site].empty:
            result[site] = wq_out
            n_missing += len(wq_out)
            continue

        q_df   = discharge_data[site].copy()
        q_df.index = pd.to_datetime(q_df.index)
        # NWIS returns tz-aware timestamps; strip so comparisons with
        # tz-naive WQ dates work correctly — values are unchanged.
        if hasattr(q_df.index, 'tz') and q_df.index.tz is not None:
            q_df.index = q_df.index.tz_localize(None)

        # Reindex discharge to a complete daily date range and interpolate
        # short gaps (<=7 days) — matches EGRET fillMissing() behaviour
        full_idx   = pd.date_range(q_df.index.min(), q_df.index.max(), freq="D")
        q_reindexed = q_df["Q_cms"].reindex(full_idx)
        gap_size    = q_reindexed.isna().astype(int)

        # Count consecutive missing days using cumsum trick
        gap_groups  = (gap_size != gap_size.shift()).cumsum()
        gap_lengths = gap_size.groupby(gap_groups).transform("sum")

        # Interpolate only gaps <= 7 days
        q_interp = q_reindexed.copy()
        short_gap_mask = (gap_size == 1) & (gap_lengths <= 7)
        q_interp[short_gap_mask] = np.nan   # will be filled by interpolation
        q_interp = q_interp.interpolate(method="time", limit=7)

        q_cfs_reindexed = q_df["Q_cfs"].reindex(full_idx).interpolate(
            method="time", limit=7
        )

        # Match each WQ observation date
        wq_dates = pd.to_datetime(wq_out["ActivityStartDate"])
        for idx, wq_date in zip(wq_out.index, wq_dates):
            wq_date_only = wq_date.normalize()
            if wq_date_only in q_interp.index:
                q_val   = q_interp.loc[wq_date_only]
                q_cfs   = q_cfs_reindexed.loc[wq_date_only]
                was_orig = not q_reindexed.isna().loc[wq_date_only]
                gap_d    = 0 if was_orig else int(gap_lengths.loc[wq_date_only])

                if pd.notna(q_val) and q_val > 0:
                    wq_out.at[idx, "Q_cms"]      = q_val
                    wq_out.at[idx, "Q_cfs"]      = q_cfs
                    wq_out.at[idx, "log_Q"]      = np.log(q_val)
                    wq_out.at[idx, "Q_matched"]  = True
                    wq_out.at[idx, "Q_gap_days"] = gap_d
                    if was_orig:
                        n_matched += 1
                    else:
                        n_interp  += 1
                else:
                    n_missing += 1
            else:
                n_missing += 1

        result[site] = wq_out

    total = n_matched + n_interp + n_missing
    log.info(
        f"  Discharge matching complete: "
        f"{n_matched} exact ({100*n_matched/total:.1f}%), "
        f"{n_interp} interpolated ({100*n_interp/total:.1f}%), "
        f"{n_missing} missing ({100*n_missing/total:.1f}%)"
    )
    return result


def compute_mean_org_alk(
    doc_data: dict[str, pd.DataFrame],
    sites: list[str],
    log: logging.Logger,
) -> dict[str, float]:
    """
    Compute mean organic alkalinity per site.
    Used to correct measured total alkalinity -> carbonate alkalinity.
    """
    mean_org_alk = {}
    for site in sites:
        if site in doc_data and not doc_data[site].empty:
            mean_org_alk[site] = doc_data[site]["org_alk"].mean()
    log.info(f"  Mean organic alkalinity computed for {len(mean_org_alk)} sites.")
    return mean_org_alk


def apply_liu2020_corrections(
    all_data_site: dict[str, pd.DataFrame],
    mean_org_alk: dict[str, float],
    log: logging.Logger,
) -> dict[str, pd.DataFrame]:
    """
    Apply pH and alkalinity corrections following Liu et al. (2020).

    pH correction:
      ionic_strength = Conductivity * 1.3e-5
      pH_error = 0.03 + 0.05 * log10(ionic_strength)
      pH_correct = pH - pH_error

    Alkalinity correction:
      carbalk = Alkalinity - mean_organic_alkalinity_for_site

    Vectorized: replaces the nested for(j) for(i) loops in R.
    """
    log.info("Applying Liu et al. 2020 pH and alkalinity corrections...")
    corrected = {}
    for site, df in all_data_site.items():
        df = df.copy()
        cond = pd.to_numeric(df["Cond_value"], errors="coerce")
        ph   = pd.to_numeric(df["pH_value"],   errors="coerce")
        alk  = pd.to_numeric(df["Alk_value"],  errors="coerce")

        # pH correction
        ion_strength     = cond * 1.3e-5
        ph_error         = 0.03 + 0.05 * np.log10(ion_strength.clip(lower=1e-12))
        df["ion_strength"] = ion_strength
        df["pH_error"]     = ph_error
        df["pH_correct"]   = ph - ph_error

        # Alkalinity correction
        org = mean_org_alk.get(site, np.nan)
        df["carbalk"] = alk - org if not np.isnan(org) else np.nan

        corrected[site] = df

    log.info("  Corrections applied.")
    return corrected


# =============================================================================
# STEP 3 — TREND ANALYSIS
# Equivalent: ts_analysis_240407_seasonal.R
# =============================================================================

def theil_sen_slope(x: np.ndarray, y: np.ndarray) -> dict:
    """
    Theil-Sen slope estimator with Mann-Kendall p-value.
    Replaces mblm::mblm() from R.

    Units of returned slope depend on x and y units:
      x is always in DAYS (days since 1900-01-01)
      y for H+  is in mol/kg        -> slope in mol/kg/day
      y for ALK is in µeq/kg        -> slope in µeq/kg/day

    To convert to scientifically meaningful annual rates:
      H+  annual: slope * 365 * 1e9   -> nmol/kg/year  (matches R: pH_ts_slope*365*1e9)
      ALK annual: slope * 365          -> µeq/kg/year   (matches R: alk_ts_slope*365)

    The raw per-day slopes are stored in the output CSV.
    The per-year slopes are stored alongside them as *_annual columns.
    """
    # Theil-Sen via scipy (fast C implementation)
    res_ts  = stats.theilslopes(y, x)
    res_lm  = stats.linregress(x, y)

    # Mann-Kendall p-value (non-parametric significance for the slope)
    mk_result = mk.original_test(y)

    return {
        "ts_slope":     res_ts.slope,
        "ts_intercept": res_ts.intercept,
        "ts_p_value":   mk_result.p,
        "ts_tau":       mk_result.Tau,
        "lm_slope":     res_lm.slope,
        "lm_intercept": res_lm.intercept,
        "lm_p_value":   res_lm.pvalue,
        "lm_r2":        res_lm.rvalue ** 2,
    }


def seasonal_kendall(
    dates: pd.Series, y: np.ndarray, season_col: pd.Series
) -> dict:
    """
    Seasonal Mann-Kendall test using pymannkendall.
    Replaces rkt::rkt(date, y, block=month) from R.

    Units of sk_slope:
      pymannkendall.seasonal_test() returns slope in units of y per
      OBSERVATION INDEX (not per day or per year).
      This is NOT directly comparable to the Theil-Sen per-day slope.

      We therefore do NOT store sk_slope as an annual rate.
      Instead we store sk_tau (the rank correlation) and sk_p_value
      for significance testing, and use the Theil-Sen annual slope
      for the magnitude of change. This matches how rkt() is used
      in sk_analysis_240228.R — for significance only, not slope magnitude.
    """
    try:
        result = mk.seasonal_test(y, period=12)
        return {
            "sk_slope":   result.slope,
            "sk_p_value": result.p,
            "sk_tau":     result.Tau,
            "sk_s":       result.s,
        }
    except Exception:
        return {
            "sk_slope":   np.nan,
            "sk_p_value": np.nan,
            "sk_tau":     np.nan,
            "sk_s":       np.nan,
        }


def compute_trends(
    all_data_site: dict[str, pd.DataFrame],
    site_meta: pd.DataFrame,
    discharge_data: dict[str, pd.DataFrame],
    cfg: dict,
    log: logging.Logger,
) -> pd.DataFrame:
    """
    Compute Theil-Sen and Seasonal Kendall trends for pH (as H+) and
    alkalinity at each site.  Generates per-site diagnostic plots to PDF.
    Replaces ts_analysis_240407_seasonal.R.
    """
    log.info("Computing trend analyses...")
    records = []
    # Reference date for DOY calculation.
    # R used "0000-01-01" but pandas 2.x nanosecond precision overflows
    # for dates before ~1677. "1900-01-01" is safe and scientifically
    # equivalent — Theil-Sen slope is invariant to the reference date.
    origin = pd.Timestamp("1900-01-01")

    for site, df in tqdm(all_data_site.items(), desc="Trend analysis"):
        df = df.sort_values("ActivityStartDate").copy()
        df["ActivityStartDate"] = pd.to_datetime(df["ActivityStartDate"])

        # Days-since-origin (matches R's difftime to "0000-01-01")
        doy = (df["ActivityStartDate"] - origin).dt.days.values.astype(float)

        # H+ concentration from corrected pH
        ph_correct = pd.to_numeric(df.get("pH_correct", df["pH_value"]), errors="coerce")
        hplus = 10 ** (-ph_correct.values)
        alk   = pd.to_numeric(df["Alk_value"], errors="coerce").values

        # Drop rows where either variable is NaN
        valid = np.isfinite(hplus) & np.isfinite(alk) & np.isfinite(doy)
        if valid.sum() < cfg["min_obs_trend"]:
            continue

        doy_v   = doy[valid]
        hplus_v = hplus[valid]
        alk_v   = alk[valid]

        duration_days = doy_v.max() - doy_v.min()
        if duration_days < cfg["min_duration_yr"] * 365:
            continue

        # Theil-Sen + MK for pH (H+)
        ph_trends = theil_sen_slope(doy_v, hplus_v)

        # Theil-Sen + MK for alkalinity
        alk_trends = theil_sen_slope(doy_v, alk_v)

        # Seasonal Kendall (monthly blocks)
        months = df["ActivityStartDate"][valid].dt.month
        sk_ph  = seasonal_kendall(
            df["ActivityStartDate"][valid], hplus_v, months
        )
        sk_alk = seasonal_kendall(
            df["ActivityStartDate"][valid], alk_v, months
        )

        # Mean discharge — NWIS daily values use Q_cms column (m³/s)
        mean_q = np.nan
        if site in discharge_data and not discharge_data[site].empty:
            q_df = discharge_data[site]
            q_col = "Q_cms" if "Q_cms" in q_df.columns else "Q_cfs"
            mean_q = pd.to_numeric(q_df[q_col], errors="coerce").mean()

        # Lat/lon from metadata
        lat = lon = np.nan
        site_no = site.replace("USGS-", "")
        if not site_meta.empty and "site_no" in site_meta.columns:
            row = site_meta[site_meta["site_no"] == site_no]
            if not row.empty:
                lat = row["dec_lat_va"].values[0]
                lon = row["dec_long_va"].values[0]

        # Reconstruct initial/final values from Theil-Sen line
        rec = {
            "site":            site,
            "n":               valid.sum(),
            "min_doy":         doy_v.min(),
            "max_doy":         doy_v.max(),
            "duration_yr":     duration_days / 365,
            "start_date":      str(df["ActivityStartDate"].iloc[0].date()),
            "end_date":        str(df["ActivityStartDate"].iloc[-1].date()),
            "lat":             lat,
            "lon":             lon,
            "mean_discharge":  mean_q,

            # ── pH (H+) Theil-Sen ────────────────────────────────────────────────
            # Raw slope: mol/kg/day  (x=days, y=10^-pH in mol/kg)
            # Annual slope: nmol/kg/year = raw * 365 * 1e9
            # Matches R output: trends_final$pH_ts_slope*365*1e9
            "pH_ts_slope_per_day":     ph_trends["ts_slope"],
            "pH_ts_slope_nmol_per_yr": ph_trends["ts_slope"] * 365 * 1e9,
            "pH_ts_intercept":         ph_trends["ts_intercept"],
            "pH_ts_p_value":           ph_trends["ts_p_value"],
            "pH_ts_tau":               ph_trends["ts_tau"],
            # LM slope also in mol/kg/day; annual = * 365 * 1e9
            "pH_lm_slope_per_day":     ph_trends["lm_slope"],
            "pH_lm_slope_nmol_per_yr": ph_trends["lm_slope"] * 365 * 1e9,
            "pH_lm_p_value":           ph_trends["lm_p_value"],
            "pH_lm_r2":                ph_trends["lm_r2"],

            # ── pH Seasonal Kendall ───────────────────────────────────────────
            # sk_p_value and sk_tau used for significance only.
            # Do not use sk_slope as a rate — units are per-observation-index,
            # not per day or year. Use pH_ts_slope_nmol_per_yr for magnitude.
            "pH_sk_p_value":           sk_ph["sk_p_value"],
            "pH_sk_tau":               sk_ph["sk_tau"],

            # ── Alkalinity Theil-Sen ─────────────────────────────────────────
            # Raw slope: µeq/kg/day  (x=days, y=alkalinity in µeq/kg)
            # Annual slope: µeq/kg/year = raw * 365
            # Matches R output: trends_final$alk_ts_slope*365
            "alk_ts_slope_per_day":    alk_trends["ts_slope"],
            "alk_ts_slope_ueq_per_yr": alk_trends["ts_slope"] * 365,
            "alk_ts_intercept":        alk_trends["ts_intercept"],
            "alk_ts_p_value":          alk_trends["ts_p_value"],
            "alk_ts_tau":              alk_trends["ts_tau"],
            # LM slope also in µeq/kg/day; annual = * 365
            "alk_lm_slope_per_day":    alk_trends["lm_slope"],
            "alk_lm_slope_ueq_per_yr": alk_trends["lm_slope"] * 365,
            "alk_lm_p_value":          alk_trends["lm_p_value"],
            "alk_lm_r2":               alk_trends["lm_r2"],

            # ── Alkalinity Seasonal Kendall ───────────────────────────────────
            # Same note as pH SK: use for significance only, not magnitude.
            "alk_sk_p_value":          sk_alk["sk_p_value"],
            "alk_sk_tau":              sk_alk["sk_tau"],
        }

        # ── Total change over record (from TS line) ──────────────────────────
        # Uses raw per-day slope × day span = total change in original units:
        #   H+  delta in mol/kg  (multiply by 1e9 for nmol/kg)
        #   ALK delta in µeq/kg
        for pfx in ["pH", "alk"]:
            slope = rec[f"{pfx}_ts_slope_per_day"]
            intcp = rec[f"{pfx}_ts_intercept"]
            rec[f"{pfx}_ts_initial"] = slope * rec["min_doy"] + intcp
            rec[f"{pfx}_ts_final"]   = slope * rec["max_doy"] + intcp
            rec[f"{pfx}_ts_delta"]   = (
                rec[f"{pfx}_ts_final"] - rec[f"{pfx}_ts_initial"]
            )
        # Also store deltas in annual units for convenience
        rec["pH_ts_delta_nmol"]  = rec["pH_ts_delta"]  * 1e9   # mol/kg -> nmol/kg
        rec["alk_ts_delta_ueq"]  = rec["alk_ts_delta"]          # already µeq/kg

        records.append(rec)

    trends = pd.DataFrame(records)
    log.info(f"  Trend analysis complete for {len(trends)} sites.")
    return trends


# =============================================================================
# STEP 4 — CALCIUM DATA
# Equivalent: neaa_calcium.R
# =============================================================================

def download_calcium(
    sites: list[str], cfg: dict, log: logging.Logger
) -> dict[str, pd.DataFrame]:
    """Download filtered calcium data (pCode 00915)."""
    log.info("Downloading calcium data...")
    result = {}
    for site in tqdm(sites, desc="Calcium download"):
        try:
            df, _ = wqp.get_results(siteid=site, param=cfg["calcium_pcodes"])
            if df is not None and not df.empty:
                df["ResultMeasureValue"] = pd.to_numeric(
                    df["ResultMeasureValue"], errors="coerce"
                )
                result[site] = df
        except Exception as e:
            log.debug(f"  Calcium download failed for {site}: {e}")
    log.info(f"  Calcium data available for {len(result)} sites.")
    return result


def join_calcium(
    all_data_site: dict[str, pd.DataFrame],
    calcium_data: dict[str, pd.DataFrame],
    log: logging.Logger,
) -> dict[str, pd.DataFrame]:
    """
    Join calcium data to per-site datasets on ActivityStartDate.
    Replaces inner_join(all_data_site_final[[j]], download_calcium2[[j]], by='ActivityStartDate').
    """
    result = {}
    for site in all_data_site:
        if site not in calcium_data or calcium_data[site].empty:
            continue
        base = all_data_site[site].copy()
        ca   = calcium_data[site][["ActivityStartDate", "ResultMeasureValue"]].copy()
        ca   = ca.rename(columns={"ResultMeasureValue": "calcium_mg_L"})
        joined = base.merge(ca, on="ActivityStartDate", how="inner")
        if not joined.empty:
            result[site] = joined
    log.info(f"  Calcium joined for {len(result)} sites.")
    return result


# =============================================================================
# STEP 5 — CODAP OCEAN DATA
# Equivalent: codap_load.R
# =============================================================================

def _codap_url_fallbacks() -> list[str]:
    """Return a list of known CODAP-NA xlsx URLs to try in order."""
    return [
        # NOAA NCEI OCADS — primary (try both v2021 and v2020 filenames)
        "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0219960/CODAP_NA_v2021.xlsx",
        "https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0219960/CODAP_NA_v2020.xlsx",
        # NOAA OCADS synthesis page direct link
        "https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system/synthesis/CODAP_NA_v2021.xlsx",
    ]


def load_codap(project_dir: Path, cfg: dict, log: logging.Logger) -> pd.DataFrame:
    """
    Load CODAP-NA v2020 dataset.
    - Checks for a local parquet cache first (fast)
    - Falls back to downloading the xlsx and caching it
    Replaces: read_excel("CODAP_NA_v2020.xlsx") in codap_load.R
    """
    cache_path = project_dir / cfg["codap_cache"]

    if cache_path.exists():
        log.info(f"Loading CODAP from cache: {cache_path}")
        return pd.read_parquet(cache_path)

    # Try download
    log.info("CODAP cache not found — attempting download...")
    xlsx_path = project_dir / "CODAP_NA_v2020.xlsx"

    if not xlsx_path.exists():
        urls = _codap_url_fallbacks()
        downloaded = False
        for url in urls:
            try:
                log.info(f"  Trying: {url}")
                r = requests.get(url, stream=True, timeout=120)
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                with open(xlsx_path, "wb") as f, tqdm(
                    total=total, unit="B", unit_scale=True, desc="CODAP download"
                ) as bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
                        bar.update(len(chunk))
                downloaded = True
                break
            except Exception as e:
                log.warning(f"  Download failed ({url}): {e}")

        if not downloaded:
            raise FileNotFoundError(
                "Could not download CODAP-NA v2020 from any known URL.\n"
                "Please manually download CODAP_NA_v2020.xlsx from:\n"
                "  https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/\n"
                f"and place it in: {project_dir}"
            )

    log.info("Reading CODAP xlsx (this may take a minute)...")
    # Row 0 = column names, row 1 = units (skip it)
    raw = pd.read_excel(xlsx_path, header=0)
    codap = raw.iloc[1:].reset_index(drop=True)   # drop units row

    # Coerce numeric columns
    num_cols = [
        "Latitude", "Longitude", "Depth",
        "recommended_Salinity_PSS78", "CTDTEMP_ITS90",
        "DIC", "TALK",
        "recommended_Salinity_flag", "CTDTEMP_flag",
        "DIC_flag", "TALK_flag",
    ]
    for col in num_cols:
        if col in codap.columns:
            codap[col] = pd.to_numeric(codap[col], errors="coerce")

    # Surface only
    codap = codap[codap["Depth"] <= cfg["codap_depth_max_m"]].copy()

    # Remove bad quality flags
    for flag_col in ["recommended_Salinity_flag", "CTDTEMP_flag", "DIC_flag", "TALK_flag"]:
        if flag_col in codap.columns:
            codap = codap[~codap[flag_col].isin(cfg["codap_bad_flags"])]

    # Create unique station ID
    codap["NEAA_ID"] = codap["Cruise_ID"].astype(str) + " " + codap["Station_ID"].astype(str)
    codap = codap.reset_index(drop=True)

    # Coerce all object columns to string before caching —
    # CODAP has mixed-type columns (e.g. EXPOCODE has both str and int)
    # that PyArrow rejects unless explicitly converted.
    for col in codap.select_dtypes(include="object").columns:
        codap[col] = codap[col].astype(str)

    # Cache as parquet for fast future loads
    codap.to_parquet(cache_path, index=False)
    log.info(
        f"  CODAP loaded: {len(codap):,} surface records. "
        f"Cached to {cache_path.name}"
    )
    return codap


# =============================================================================
# STEP 6 — MATCH RIVER SITES TO OCEAN STATIONS
# Equivalent: codap_organize.R + close_ocean.R
# =============================================================================

def find_close_ocean_stations(
    site_id: str,
    site_lat: float,
    site_lon: float,
    codap: pd.DataFrame,
    max_dist_m: float,
    log: logging.Logger,
) -> pd.DataFrame:
    """
    Return CODAP surface records within max_dist_m of a river site.
    Replaces close_ocean() function (geosphere::distHaversine).
    Uses geopy.great_circle which matches distHaversine closely.
    """
    codap_coords = list(zip(codap["Latitude"], codap["Longitude"]))
    river_coord  = (site_lat, site_lon)

    dists = np.array([
        great_circle(river_coord, c).meters for c in codap_coords
    ])

    if dists.min() < max_dist_m:
        mask = dists < max_dist_m
    else:
        # Fall back to nearest station (matches R behavior)
        mask = dists == dists.min()

    return codap[mask].copy()


def compute_ocean_carbonate(
    ocean_data: pd.DataFrame,
    log: logging.Logger,
) -> pd.DataFrame | None:
    """
    Calculate full carbonate system for ocean stations using PyCO2SYS.
    Replaces seacarb::carbfull() called with flag=15 (ALK + DIC input).
    PyCO2SYS flag 1 = ALK, flag 2 = DIC (par1 = TALK µmol/kg, par2 = DIC µmol/kg).
    """
    ocean_data = ocean_data.dropna(
        subset=["TALK", "DIC", "CTDTEMP_ITS90", "recommended_Salinity_PSS78"]
    )
    if ocean_data.empty:
        return None

    # PyCO2SYS v1.8 key names:
    #   pressure = 0  (not pressure_in)
    #   OmegaAR / OmegaCA  (not OmegaARout/OmegaCAout — those only exist
    #   when temperature_out or pressure_out are specified)
    try:
        results = co2sys.sys(
            par1           = ocean_data["TALK"].values,   # µmol/kg
            par2           = ocean_data["DIC"].values,    # µmol/kg
            par1_type      = 1,    # Total Alkalinity
            par2_type      = 2,    # Dissolved Inorganic Carbon
            temperature    = ocean_data["CTDTEMP_ITS90"].values,
            salinity       = ocean_data["recommended_Salinity_PSS78"].values,
            pressure       = 0,    # surface (dbar)
            opt_pH_scale   = 1,    # Total scale (pHscale="T" in seacarb)
            opt_k_carbonic = 10,   # Lueker et al. 2000 (k1k2="l" in seacarb)
        )
        out = ocean_data.copy()
        out["pCO2"]     = results["pCO2"]
        out["pH_total"] = results["pH"]
        # OmegaAR/OmegaCA (no _out suffix — only appears with output conditions)
        out["omega_ar"] = results.get("OmegaAR",    results.get("OmegaARout",    np.nan))
        out["omega_ca"] = results.get("OmegaCA",    results.get("OmegaCAout",    np.nan))
        return out
    except Exception as e:
        log.warning(f"  PyCO2SYS ocean calc failed: {e}")
        log.debug(f"  Available result keys: {list(results.keys()) if 'results' in dir() else 'N/A'}")
        return None


def organize_ocean_data(
    sites: list[str],
    site_meta: pd.DataFrame,
    codap: pd.DataFrame,
    cfg: dict,
    log: logging.Logger,
) -> dict[str, pd.DataFrame]:
    """
    For each river site, find nearby ocean stations and compute carbonate chemistry.
    Replaces the for-loop in codap_organize.R.
    """
    log.info("Matching river sites to CODAP ocean stations...")
    all_ocean = {}

    for site in tqdm(sites, desc="Ocean matching"):
        site_no = site.replace("USGS-", "")
        lat = lon = np.nan

        # Primary: use NWIS metadata already fetched (dec_lat_va / dec_long_va)
        if not site_meta.empty and "site_no" in site_meta.columns:
            row = site_meta[site_meta["site_no"] == site_no]
            if not row.empty:
                lat = float(row["dec_lat_va"].values[0])
                lon = float(row["dec_long_va"].values[0])

        # Fallback: query WQP directly (dataretrieval >= 1.0).
        # what_sites() returns "lon" and "lat" — verified against live API.
        # Replaces old whatWQPdata() call in close_ocean.R.
        if np.isnan(lat) or np.isnan(lon):
            try:
                site_info, _ = wqp.what_sites(siteid=site)
                if site_info is not None and not site_info.empty:
                    if "lat" in site_info.columns and "lon" in site_info.columns:
                        lat = float(site_info["lat"].values[0])
                        lon = float(site_info["lon"].values[0])
            except Exception as e:
                log.debug(f"  WQP site lookup failed for {site}: {e}")

        if np.isnan(lat) or np.isnan(lon):
            log.debug(f"  No lat/lon for {site}, skipping ocean match.")
            continue

        nearby = find_close_ocean_stations(
            site, lat, lon, codap, cfg["ocean_radius_m"], log
        )
        if nearby.empty:
            log.debug(f"  No ocean stations within radius for {site}.")
            continue

        ocean_cc = compute_ocean_carbonate(nearby, log)
        if ocean_cc is not None:
            all_ocean[site] = ocean_cc

    log.info(f"  Ocean carbonate computed for {len(all_ocean)} sites.")
    return all_ocean


# =============================================================================
# STEP 7 — EXPORT
# Equivalent: export_neaa_csv_calcium.R + export_ocean_csv.R
# =============================================================================

def export_csvs(
    all_data_site: dict[str, pd.DataFrame],
    calcium_data:  dict[str, pd.DataFrame],
    all_ocean:     dict[str, pd.DataFrame],
    trends:        pd.DataFrame,
    site_meta:     pd.DataFrame,
    discharge_data: dict[str, pd.DataFrame],
    export_dir:    Path,
    log:           logging.Logger,
) -> None:
    """
    Write per-site CSVs and summary CSVs.
    Matches the R export scripts exactly in output structure.
    """
    log.info(f"Exporting CSVs to {export_dir} ...")

    # Per-site river data
    for site, df in all_data_site.items():
        df.to_csv(export_dir / f"{site}.csv", index=True)

    # Per-site calcium data
    for site, df in calcium_data.items():
        df.to_csv(export_dir / f"{site}_calcium.csv", index=True)

    # Per-site ocean carbonate data
    for site, df in all_ocean.items():
        df.to_csv(export_dir / f"{site}_ocean.csv", index=True)

    # Summary files
    trends.to_csv(export_dir / "ts_regressions.csv", index=False)

    if not site_meta.empty:
        site_meta.to_csv(export_dir / "siteINFO.csv", index=False)

    # Discharge statistics (from NWIS daily record)
    if "discharge_stats" in dir():
        discharge_stats.to_csv(export_dir / "discharge_stats.csv", index=False)
    else:
        disc_rows = []
        for site, df in discharge_data.items():
            q = df["Q_cms"].dropna() if "Q_cms" in df.columns else pd.Series()
            disc_rows.append({"site": site, "mean_discharge": q.mean()})
        if disc_rows:
            pd.DataFrame(disc_rows).to_csv(
                export_dir / "discharge_stats.csv", index=False
            )

    log.info(f"  CSV export complete.")


def export_hdf5(
    all_data_site: dict[str, pd.DataFrame],
    calcium_data:  dict[str, pd.DataFrame],
    all_ocean:     dict[str, pd.DataFrame],
    trends:        pd.DataFrame,
    site_meta:     pd.DataFrame,
    discharge_data: dict[str, pd.DataFrame],
    export_dir:    Path,
    log:           logging.Logger,
) -> None:
    """
    Export all data to a single HDF5 file for easy MATLAB import.

    In MATLAB, load with:
        data = h5read('neaa_export.h5', '/river_data/USGS-01010000');
        trends = h5read('neaa_export.h5', '/summaries/trends');

    Structure:
        /river_data/<site_id>/      — per-site water quality data
        /calcium/<site_id>/         — per-site calcium data
        /ocean/<site_id>/           — per-site ocean carbonate data
        /summaries/trends           — Theil-Sen + SK regression results
        /summaries/site_metadata    — NWIS site info
        /summaries/discharge_means  — mean discharge per site
    """
    h5_path = export_dir / "neaa_export.h5"
    log.info(f"Exporting HDF5 to {h5_path} ...")

    def _df_to_h5(h5file, group_path: str, df: pd.DataFrame) -> None:
        """Write a DataFrame into an HDF5 group, one dataset per column.
        Skips geometry columns (geopandas Point objects) and safely handles
        mixed-type object columns that contain non-numeric values.
        """
        grp = h5file.require_group(group_path)
        for col in df.columns:
            series = df[col]
            # Skip geometry columns from geopandas (Point, LineString, etc.)
            if str(series.dtype) == "geometry":
                continue
            try:
                if pd.api.types.is_datetime64_any_dtype(series):
                    data = series.astype(str).values.astype("S")
                elif series.dtype == "object":
                    data = series.astype(str).values.astype("S")
                else:
                    try:
                        data = series.values.astype(float)
                    except (TypeError, ValueError):
                        # Fallback for object columns with non-numeric content
                        data = series.astype(str).values.astype("S")
                if col in grp:
                    del grp[col]
                grp.create_dataset(col, data=data, compression="gzip")
            except Exception:
                pass  # skip any column that cannot be serialised

    with h5py.File(h5_path, "w") as h5:
        h5.attrs["created"] = datetime.now().isoformat()
        h5.attrs["description"] = "NEAA water quality pipeline output"

        for site, df in all_data_site.items():
            safe = site.replace("-", "_")
            _df_to_h5(h5, f"river_data/{safe}", df.reset_index(drop=True))

        for site, df in calcium_data.items():
            safe = site.replace("-", "_")
            _df_to_h5(h5, f"calcium/{safe}", df.reset_index(drop=True))

        for site, df in all_ocean.items():
            safe = site.replace("-", "_")
            _df_to_h5(h5, f"ocean/{safe}", df.reset_index(drop=True))

        if not trends.empty:
            _df_to_h5(h5, "summaries/trends", trends.reset_index(drop=True))

        if not site_meta.empty:
            _df_to_h5(h5, "summaries/site_metadata", site_meta.reset_index(drop=True))

        # Write per-site daily discharge records
        for site, df in discharge_data.items():
            safe = site.replace("-", "_")
            _df_to_h5(h5, f"discharge/{safe}",
                      df.reset_index().rename(columns={"date": "date"}))

        # Write discharge statistics summary
        if "discharge_stats" in dir():
            _df_to_h5(h5, "summaries/discharge_stats",
                      discharge_stats.reset_index(drop=True))

    log.info(f"  HDF5 export complete: {h5_path.stat().st_size / 1e6:.1f} MB")


# =============================================================================
# MAIN ORCHESTRATOR
# Equivalent: master_neaa_final_240417.R
# =============================================================================

def main() -> None:
    project_dir, export_dir = setup_dirs(CFG)
    log = setup_logging(project_dir / "logs")

    log.info("=" * 60)
    log.info("NEAA Python Pipeline — starting")
    log.info(f"  Project dir : {project_dir}")
    log.info(f"  Export dir  : {export_dir}")
    log.info("=" * 60)

    # Startup path check — catch mismatches before any downloads
    site_list_path = project_dir / CFG["site_list_xlsx"]
    if not site_list_path.exists():
        log.error(
            f"Site list not found: {site_list_path}\n"
            f"Check that 'project_dir' in CFG matches where your files are.\n"
            f"Current setting: {project_dir}\n"
            f"Your script is running from: {Path.cwd()}"
        )
        raise FileNotFoundError(f"Site list not found: {site_list_path}")

    # ── Step 1: Download & QA river data ─────────────────────────────────────
    # Checkpoint cache — skip re-downloading if already done this session.
    # Delete neaa_checkpoint.parquet from your project folder to force a
    # fresh download (e.g. after updating the site list or QA thresholds).
    checkpoint_path = project_dir / "neaa_checkpoint.parquet"
    checkpoint_meta = project_dir / "neaa_checkpoint_meta.csv"

    if checkpoint_path.exists() and checkpoint_meta.exists():
        log.info(f"Loading from checkpoint: {checkpoint_path.name}")
        log.info("  (Delete neaa_checkpoint.parquet to force a fresh download)")
        _ckpt = pd.read_parquet(checkpoint_path)
        _meta = pd.read_csv(checkpoint_meta)

        # Reconstruct all_data_site dict from checkpoint
        all_data_site_final = {}
        for site, grp in _ckpt.groupby("MonitoringLocation"):
            all_data_site_final[site] = grp.copy()
        names_final = list(all_data_site_final.keys())
        site_meta   = _meta
        log.info(f"  Loaded {len(names_final)} sites from checkpoint.")
    else:
        sites = load_site_list(project_dir, CFG, log)

        raw_data = download_wqp_data(sites, CFG, log)

        log.info("Extracting and QA/QC-ing individual parameters...")
        pH_df   = extract_parameter(raw_data, "pH",                   CFG["ph_min"],   CFG["ph_max"],   prefer_field=True,  log=log)
        temp_df = extract_parameter(raw_data, "Temperature, water",   CFG["temp_min_c"], CFG["temp_max_c"],                 log=log)
        alk_df  = extract_parameter(raw_data, "Alkalinity",           CFG["alk_min"] / CFG["alk_multiplier"],
                                                                       CFG["alk_max"]  / CFG["alk_multiplier"], prefer_field=True, log=log)
        cond_df = extract_parameter(raw_data, "Specific conductance", CFG["cond_min"], CFG["cond_max"],                     log=log)

        # Unit-convert alkalinity (matches R: alk$ResultMeasureValue * 20)
        alk_df = alk_df.copy()
        alk_df["ResultMeasureValue"] = alk_df["ResultMeasureValue"] * CFG["alk_multiplier"]

        merged = merge_parameters(pH_df, temp_df, alk_df, cond_df, log)

        all_data_site = build_site_datasets(merged, sites, CFG, log)

        # Log freshwater classification method summary across all sites
        all_methods = []
        for df in all_data_site.values():
            if "fw_method" in df.columns:
                all_methods.extend(df["fw_method"].tolist())
        if all_methods:
            method_counts = pd.Series(all_methods).value_counts()
            log.info("  Freshwater classification methods used (observation counts):")
            for method, count in method_counts.items():
                log.info(f"    {method}: {count:,} observations")

        all_data_site_final, all_data_site_recent, all_data_site_2000 = \
            get_site_subsets(all_data_site, CFG, log)

        names_final = list(all_data_site_final.keys())
        site_meta   = fetch_site_metadata(names_final, log)

        # Save checkpoint so re-runs skip the download
        log.info(f"Saving checkpoint to {checkpoint_path.name} ...")
        _all = pd.concat(all_data_site_final.values(), ignore_index=True)
        # Coerce mixed-type columns before writing parquet
        for col in _all.select_dtypes(include="object").columns:
            _all[col] = _all[col].astype(str)
        _all.to_parquet(checkpoint_path, index=False)
        if not site_meta.empty:
            site_meta.to_csv(checkpoint_meta, index=False)
        log.info("  Checkpoint saved. Re-runs will skip WQP download.")

    # ── Step 1b: Chloride download for freshwater classification ─────────────
    # Downloaded for all sites — used by classify_freshwater() to resolve
    # ambiguous conductance values (1500–10000 µS/cm).
    # Sites with cond < 1500 µS/cm don't need this but it's cheap to download.
    log.info("Downloading chloride data for freshwater classification...")
    chloride_data = download_chloride(names_final, CFG, log)
    # Attach to CFG so build_site_datasets can access it without changing
    # the function signature — underscore prefix marks it as internal
    CFG["_chloride_data"] = chloride_data

    # ── Step 2: Corrections ───────────────────────────────────────────────────
    doc_data       = download_doc(names_final, CFG, log)
    mean_org_alk   = compute_mean_org_alk(doc_data, names_final, log)

    all_data_site_final = apply_liu2020_corrections(
        all_data_site_final, mean_org_alk, log
    )

    # ── Step 2b: NWIS daily discharge ────────────────────────────────────────
    # Downloads the complete daily record for each site — far more complete
    # than the WQP grab samples used in the original R scripts.
    # Required for: flow statistics, SK discharge covariate, WRTDS.
    discharge_data = download_discharge(names_final, CFG, log)
    discharge_stats = compute_discharge_stats(discharge_data, log)

    # Match each WQ observation to its concurrent daily discharge value.
    # Adds Q_cms, Q_cfs, log_Q columns to each site DataFrame.
    # log_Q is the covariate used in Seasonal Kendall and WRTDS.
    all_data_site_final = match_discharge_to_wq(
        all_data_site_final, discharge_data, log
    )

    # ── Step 3: Trend analysis ────────────────────────────────────────────────
    trends = compute_trends(
        all_data_site_final, site_meta, discharge_data, CFG, log
    )

    # ── Step 4: Calcium ───────────────────────────────────────────────────────
    calcium_data = download_calcium(names_final, CFG, log)
    calcium_joined = join_calcium(all_data_site_final, calcium_data, log)

    # ── Step 5-6: CODAP / Ocean ───────────────────────────────────────────────
    codap = load_codap(project_dir, CFG, log)
    all_ocean = organize_ocean_data(names_final, site_meta, codap, CFG, log)

    # ── Step 7: Export ────────────────────────────────────────────────────────
    export_csvs(
        all_data_site_final, calcium_joined, all_ocean,
        trends, site_meta, discharge_data, export_dir, log,
    )
    export_hdf5(
        all_data_site_final, calcium_joined, all_ocean,
        trends, site_meta, discharge_data, export_dir, log,
    )

    log.info("=" * 60)
    log.info("NEAA pipeline complete.")
    log.info(f"  Sites processed : {len(names_final)}")
    log.info(f"  Sites with trends: {len(trends)}")
    log.info(f"  Sites with ocean data: {len(all_ocean)}")
    log.info(f"  Outputs written to: {export_dir}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
