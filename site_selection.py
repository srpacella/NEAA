"""
NEAA Systematic Site Selection
================================
Replaces the manually assembled site list with a reproducible, spatially
consistent, algorithmically derived site list.

Scientific goal: identify USGS stream gauges that are:
  1. On rivers draining to the US coast (not interior basins)
  2. In freshwater reaches above tidal influence
  3. As close to the estuary as data quality allows
  4. With sufficient water quality data for trend analysis

Method:
  Step 1  Query WQP for all coastal-state stream sites with pH + alkalinity
  Step 2  Filter to USGS-gauged freshwater stream sites only
  Step 3  Compute NHDPlus downstream distance to ocean/tidal water
  Step 4  Apply distance threshold (river-path km, not straight-line)
  Step 5  Filter to sites with adequate data for trend analysis
  Step 6  Score and rank sites by data quality + proximity to coast
  Step 7  Compare against existing manual site list

Outputs:
  site_selection_results.csv   — full scored site list
  site_selection_added.csv     — sites not in current list that qualify
  site_selection_removed.csv   — current sites that fail systematic criteria
  site_selection_report.txt    — human-readable summary

Dependencies:
    pip install dataretrieval pynhd pandas numpy requests tqdm geopy

Usage:
    conda activate neaa
    python3 site_selection.py

Reference:
  NHDPlus NLDI API: https://waterdata.usgs.gov/blog/nldi-intro/
  Hirsch & De Cicco 2015 — WRTDS site selection criteria
  Stets et al. 2014 — NEAA predecessor study site criteria
"""

import logging
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm
from geopy.distance import great_circle

warnings.filterwarnings("ignore", category=FutureWarning)

try:
    import dataretrieval.wqp as wqp
    import dataretrieval.nwis as nwis
except ImportError:
    raise ImportError("pip install dataretrieval")

# =============================================================================
# CONFIGURATION
# =============================================================================
CFG = {
    # Paths
    "project_dir":      Path.home() / "Desktop" / "NEAA",
    "existing_sites":   "WQP_sites_manual_srp_update6.xlsx",
    "output_dir":       Path.home() / "Desktop" / "NEAA" / "neaa_data_exports",

    # US coastal states — matches your original R script
    "coastal_states": [
        "ME", "NH", "MA", "RI", "CT", "NY", "NJ", "DE", "MD", "VA",
        "NC", "SC", "GA", "FL", "AL", "MS", "LA", "TX",
        "CA", "OR", "WA", "AK",
    ],

    # Site type filters — keep only freshwater non-tidal streams
    "valid_site_types": [
        "Stream",
        "River/Stream",
        "Stream: Canal",
        "River/Stream: Canal",
        "River/Stream: Regulated",
    ],
    # Explicitly exclude tidal/estuarine sites
    "exclude_site_types": [
        "Estuary",
        "Stream: Tidal stream",
        "River/Stream: Tidal stream",
        "Ocean",
        "Coastal",
        "Lake, Reservoir, Impoundment",
        "Well",
        "Spring",
        "Atmosphere",
    ],

    # Data adequacy thresholds
    "min_ph_obs":       25,    # minimum pH observations
    "min_alk_obs":      25,    # minimum alkalinity observations
    "min_record_yrs":   10,    # minimum record length in years

    # Spatial filters
    "max_river_dist_km":   500,   # maximum river-path distance to coast
    "target_river_dist_km": 200,  # preferred maximum — closer sites score higher
    "max_straight_dist_km": 300,  # straight-line pre-filter (faster than NLDI)

    # Conductance threshold for freshwater pre-filter
    # Sites with mean conductance above this get flagged for review
    "cond_review_threshold": 5000,  # µS/cm

    # NLDI navigation settings
    "nldi_downstream_dist": 9999,  # km — navigate all the way to terminal

    # Scoring weights (sum to 1.0)
    # Higher score = better site for analysis
    "score_weights": {
        "proximity":    0.30,   # closer to coast = better
        "ph_obs":       0.25,   # more pH observations = better
        "alk_obs":      0.25,   # more alkalinity observations = better
        "record_yrs":   0.20,   # longer record = better
    },

    # Score thresholds for tier assignment
    "tier1_score": 0.70,   # top tier: excellent
    "tier2_score": 0.45,   # second tier: good
    "tier3_score": 0.20,   # third tier: marginal
}

# =============================================================================
# LOGGING
# =============================================================================
def setup_logging(out_dir: Path) -> logging.Logger:
    out_dir.mkdir(parents=True, exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        handlers=[
            logging.FileHandler(out_dir / f"site_selection_{ts}.log"),
            logging.StreamHandler(),
        ],
    )
    return logging.getLogger("site_select")


# =============================================================================
# STEP 1 — QUERY WQP FOR ALL CANDIDATE SITES
# =============================================================================

def query_candidate_sites(cfg: dict, log: logging.Logger) -> pd.DataFrame:
    """
    Discover all USGS freshwater stream sites in coastal states with both
    pH and alkalinity data, using WQP as the authoritative data source.

    Design rationale:
    -----------------
    wqp.what_sites(statecode, characteristicName) returns TWO things at once:
      1. The list of sites that have that characteristic
      2. resultCount — how many observations exist at each site

    This means a single pass of 44 API calls (22 states x 2 characteristics)
    gives us site discovery AND data adequacy in one step. There is no need
    for a separate data adequacy query step.

    Queries are one state x one characteristic at a time. This is the
    pattern USGS recommends for large WQP pulls — small queries are reliable,
    large multi-state multi-characteristic queries time out silently.

    Independence:
    -------------
    This function has NO dependency on the existing manual site list.
    The manual list is used only in the final comparison step.
    """
    log.info("Discovering candidate sites via WQP (standard repository)...")
    log.info("  Strategy: 44 calls (22 states x pH + alkalinity)")
    log.info("  Each call returns site list AND observation counts — no second pass needed")

    STATE_FIPS = {
        "ME":"23","NH":"33","MA":"25","RI":"44","CT":"09","NY":"36",
        "NJ":"34","DE":"10","MD":"24","VA":"51","NC":"37","SC":"45",
        "GA":"13","FL":"12","AL":"01","MS":"28","LA":"22","TX":"48",
        "CA":"06","OR":"41","WA":"53","AK":"02",
    }

    # Collect results for each characteristic
    # key = site_id, value = dict of metadata + counts
    ph_results  = {}   # site_id -> {metadata + pH resultCount}
    alk_results = {}   # site_id -> {metadata + ALK resultCount}

    for state in tqdm(cfg["coastal_states"], desc="WQP site queries"):
        fips = STATE_FIPS.get(state, state)

        for char, store in [("pH", ph_results), ("Alkalinity", alk_results)]:
            try:
                df, _ = wqp.what_sites(
                    statecode=f"US:{fips}",
                    characteristicName=char,
                    siteType="Stream",
                )
                if df is None or df.empty:
                    continue

                # Keep only USGS sites
                df = df[
                    df["MonitoringLocationIdentifier"].str.startswith("USGS-")
                ].copy()

                # resultCount column is the number of observations at this
                # site for this characteristic — already returned, no extra call
                if "resultCount" in df.columns:
                    df["resultCount"] = pd.to_numeric(
                        df["resultCount"], errors="coerce"
                    ).fillna(0)
                else:
                    df["resultCount"] = 0

                for _, row in df.iterrows():
                    sid = row["MonitoringLocationIdentifier"]
                    store[sid] = row.to_dict()

            except Exception as e:
                log.debug(f"  WQP {state} {char}: {e}")

    # Intersect: only keep sites with BOTH pH and alkalinity
    both_ids = set(ph_results.keys()) & set(alk_results.keys())

    log.info(f"  Sites with pH data:        {len(ph_results):,}")
    log.info(f"  Sites with alkalinity:     {len(alk_results):,}")
    log.info(f"  Sites with both:           {len(both_ids):,}")

    if not both_ids:
        raise RuntimeError(
            "WQP returned 0 sites with both pH and alkalinity. "
            "This usually means the API is temporarily unavailable. "
            "Wait a few minutes and try again."
        )

    # Build unified metadata DataFrame
    rows = []
    for sid in both_ids:
        ph_row  = ph_results[sid]
        alk_row = alk_results[sid]
        row = ph_row.copy()
        # Store observation counts for each characteristic separately
        # These are used later to filter by data adequacy without any
        # additional API calls
        row["ph_result_count"]  = ph_row.get("resultCount", 0)
        row["alk_result_count"] = alk_row.get("resultCount", 0)
        rows.append(row)

    all_sites = pd.DataFrame(rows).reset_index(drop=True)

    # Standardise lat/lon — WQP returns these directly in what_sites()
    for lat_col in ["lat", "LatitudeMeasure", "dec_lat_va"]:
        if lat_col in all_sites.columns:
            all_sites["lat"] = pd.to_numeric(
                all_sites[lat_col], errors="coerce"
            )
            break

    for lon_col in ["lon", "LongitudeMeasure", "dec_long_va"]:
        if lon_col in all_sites.columns:
            all_sites["lon"] = pd.to_numeric(
                all_sites[lon_col], errors="coerce"
            )
            break

    all_sites = all_sites.drop_duplicates(
        subset="MonitoringLocationIdentifier"
    ).reset_index(drop=True)

    log.info(f"  Candidate sites with metadata: {len(all_sites):,}")
    log.info(
        f"  Observation count ranges — "
        f"pH: {all_sites['ph_result_count'].min():.0f}–"
        f"{all_sites['ph_result_count'].max():.0f}, "
        f"ALK: {all_sites['alk_result_count'].min():.0f}–"
        f"{all_sites['alk_result_count'].max():.0f}"
    )
    return all_sites


def _query_wqp_rest_api(cfg: dict, log: logging.Logger) -> list:
    """
    Fallback: query WQP Station endpoint directly via REST API.
    Used when dataretrieval what_sites() returns empty results.
    """
    STATE_FIPS = {
        "ME":"23","NH":"33","MA":"25","RI":"44","CT":"09","NY":"36",
        "NJ":"34","DE":"10","MD":"24","VA":"51","NC":"37","SC":"45",
        "GA":"13","FL":"12","AL":"01","MS":"28","LA":"22","TX":"48",
        "CA":"06","OR":"41","WA":"53","AK":"02",
    }
    base = "https://www.waterqualitydata.us/Station/search"
    frames = []
    for state in tqdm(cfg["coastal_states"], desc="REST API fallback"):
        fips = STATE_FIPS.get(state, state)
        params = {
            "statecode":   f"US:{fips}",
            "siteType":    "Stream",
            "providers":   "NWIS",
            "mimeType":    "csv",
            "zip":         "no",
        }
        try:
            r = requests.get(base, params=params, timeout=60)
            if r.status_code == 200 and len(r.content) > 100:
                from io import StringIO
                df = pd.read_csv(StringIO(r.text), low_memory=False)
                df["query_state"] = state
                frames.append(df)
        except Exception as e:
            log.debug(f"  REST fallback {state}: {e}")
    return frames


# _get_sites_with_char removed — replaced by WQP summary endpoint


# =============================================================================
# STEP 2 — FILTER TO USGS FRESHWATER STREAM SITES
# =============================================================================

def filter_site_types(
    sites: pd.DataFrame, cfg: dict, log: logging.Logger
) -> pd.DataFrame:
    """
    Keep only freshwater non-tidal stream sites.
    Exclude estuaries, tidal streams, lakes, wells, and non-USGS providers.
    """
    n_start = len(sites)
    log.info(f"  Starting with {n_start:,} USGS stream sites")

    # Filter site types — NWIS uses ST prefix for streams
    if "MonitoringLocationTypeName" in sites.columns:
        n_before  = len(sites)
        type_col  = sites["MonitoringLocationTypeName"].astype(str)
        # Keep: ST (stream), ST-CA (canal), ST-DCH (ditch), River/Stream variants
        keep_mask = (
            type_col.str.startswith("ST") |
            type_col.isin(cfg["valid_site_types"])
        )
        # Exclude explicitly tidal/estuarine
        excl_mask = type_col.isin(cfg["exclude_site_types"])
        sites = sites[keep_mask & ~excl_mask].copy()
        log.info(
            f"  Stream type filter: {len(sites):,} kept "
            f"(dropped {n_before - len(sites):,} non-stream/tidal sites)"
        )

    # Require valid coordinates
    sites["lat"] = pd.to_numeric(sites.get("lat", np.nan), errors="coerce")
    sites["lon"] = pd.to_numeric(sites.get("lon", np.nan), errors="coerce")
    n_before = len(sites)
    sites    = sites.dropna(subset=["lat", "lon"]).copy()
    log.info(f"  Valid coordinates: {len(sites):,} (dropped {n_before - len(sites):,})")

    return sites.reset_index(drop=True)


# =============================================================================
# STEP 3 — SPATIAL PRE-FILTER (straight-line distance)
# =============================================================================

# Approximate coordinates of major US coastal/tidal boundaries
# Used as a fast pre-filter before the more expensive NLDI routing
COASTAL_REFERENCE_POINTS = [
    # Atlantic coast (north to south)
    (47.0, -67.8),   # Maine
    (43.1, -70.8),   # New Hampshire
    (42.0, -70.1),   # Massachusetts Cape Cod
    (41.5, -71.3),   # Rhode Island
    (41.2, -72.9),   # Connecticut
    (40.7, -74.0),   # New York Harbor
    (39.4, -74.4),   # New Jersey
    (38.9, -75.1),   # Delaware Bay
    (38.3, -76.4),   # Chesapeake Bay mouth
    (36.9, -76.0),   # Virginia Beach
    (35.9, -75.6),   # Outer Banks NC
    (34.2, -77.8),   # Cape Fear NC
    (32.8, -79.9),   # Charleston SC
    (32.0, -81.0),   # Savannah GA
    (30.4, -81.4),   # Jacksonville FL
    (25.8, -80.2),   # Miami FL
    # Gulf coast
    (30.0, -88.0),   # Mobile Bay
    (29.9, -90.1),   # New Orleans
    (29.3, -94.8),   # Galveston TX
    (27.8, -97.4),   # Corpus Christi TX
    # Pacific coast
    (32.7, -117.2),  # San Diego
    (37.8, -122.4),  # San Francisco Bay
    (46.2, -124.0),  # Columbia River mouth
    (48.4, -122.5),  # Puget Sound
    # Alaska
    (61.2, -149.9),  # Cook Inlet
    (57.0, -135.3),  # SE Alaska
]


def straight_line_to_coast(lat: float, lon: float) -> float:
    """Return minimum straight-line distance in km to any coastal reference point."""
    return min(
        great_circle((lat, lon), ref).km
        for ref in COASTAL_REFERENCE_POINTS
    )


def spatial_prefilter(
    sites: pd.DataFrame, cfg: dict, log: logging.Logger
) -> pd.DataFrame:
    """
    Fast pre-filter using straight-line distance before the expensive NLDI call.
    Sites farther than max_straight_dist_km from coast are dropped.
    """
    log.info("Applying spatial pre-filter (straight-line distance to coast)...")
    sites = sites.copy()
    if sites.empty:
        log.warning("  No sites to filter — check upstream query step.")
        return sites
    sites["dist_coast_km_straight"] = sites.apply(
        lambda r: straight_line_to_coast(float(r["lat"]), float(r["lon"]))
        if pd.notna(r["lat"]) and pd.notna(r["lon"]) else 999.0,
        axis=1,
    )
    n_before = len(sites)
    sites = sites[
        sites["dist_coast_km_straight"] <= cfg["max_straight_dist_km"]
    ].copy()
    log.info(
        f"  Straight-line filter (<{cfg['max_straight_dist_km']} km): "
        f"{len(sites):,} sites kept (dropped {n_before - len(sites):,})"
    )
    return sites.reset_index(drop=True)


# =============================================================================
# STEP 4 — NHDPlus DOWNSTREAM DISTANCE TO OCEAN
# =============================================================================

def get_nldi_downstream_distance(
    site_id: str,
    log: logging.Logger,
) -> dict:
    """
    Use the NLDI API to navigate downstream from a USGS site and compute
    the river-path distance to the terminal (ocean/tidal) reach.

    Returns dict with:
      river_dist_km     : total downstream river-path distance to terminal
      terminal_comid    : NHDPlus ComID of the terminal reach
      n_downstream_segs : number of downstream segments
      tidal_flag        : True if a tidal reach is encountered before terminal
      nldi_error        : error message if lookup failed

    NLDI API rate limit: 3600 requests/hour. For 500+ sites this means
    running overnight or batching with delays.
    """
    base_url = "https://labs.waterdata.usgs.gov/api/nldi/linked-data"
    result = {
        "river_dist_km":     np.nan,
        "terminal_comid":    None,
        "n_downstream_segs": 0,
        "tidal_flag":        False,
        "nldi_error":        None,
    }

    # Step 1: Get the NHDPlus ComID for this site
    try:
        url = f"{base_url}/nwissite/{site_id}"
        r   = requests.get(url, timeout=30)
        if r.status_code != 200:
            result["nldi_error"] = f"site lookup HTTP {r.status_code}"
            return result
        feat = r.json()
        comid = feat["features"][0]["properties"].get("comid")
        if not comid:
            result["nldi_error"] = "no comid returned"
            return result
    except Exception as e:
        result["nldi_error"] = f"site lookup: {e}"
        return result

    # Step 2: Navigate downstream to terminal
    try:
        url  = (
            f"{base_url}/comid/{comid}/navigation/DM/flowlines"
            f"?distance=9999"
        )
        r    = requests.get(url, timeout=60)
        if r.status_code != 200:
            result["nldi_error"] = f"navigation HTTP {r.status_code}"
            return result

        data  = r.json()
        feats = data.get("features", [])
        if not feats:
            result["nldi_error"] = "no downstream features"
            return result

        # Sum lengthkm across all downstream segments
        total_km = sum(
            f["properties"].get("lengthkm", 0) for f in feats
        )
        result["river_dist_km"]     = round(total_km, 2)
        result["n_downstream_segs"] = len(feats)
        result["terminal_comid"]    = feats[-1]["properties"].get("nhdplus_comid")

        # Check for tidal flag in NHDPlus attributes
        for feat in feats:
            ftype = feat["properties"].get("ftype", "")
            if "tidal" in str(ftype).lower() or "coastal" in str(ftype).lower():
                result["tidal_flag"] = True
                break

    except Exception as e:
        result["nldi_error"] = f"navigation: {e}"

    return result


def compute_river_distances(
    sites: pd.DataFrame, cfg: dict, log: logging.Logger
) -> pd.DataFrame:
    """
    Compute river-path distance to coast for each site using NLDI.
    This is the most expensive step — ~1-2 seconds per site.
    Results are cached so re-runs are instant.
    """
    cache_path = cfg["output_dir"] / "nldi_distance_cache.csv"

    # Load cache if it exists
    cached = {}
    if cache_path.exists():
        cache_df = pd.read_csv(cache_path)
        cached   = cache_df.set_index("site")["river_dist_km"].to_dict()
        log.info(f"  Loaded {len(cached)} cached NLDI distances")

    sites    = sites.copy()
    new_rows = []

    for _, row in tqdm(sites.iterrows(), total=len(sites), desc="NLDI distances"):
        site_id = row["MonitoringLocationIdentifier"]
        usgs_id = site_id.replace("USGS-", "USGS-")  # NLDI uses full USGS- prefix

        if site_id in cached:
            sites.loc[sites["MonitoringLocationIdentifier"] == site_id,
                      "river_dist_km"] = cached[site_id]
            continue

        result = get_nldi_downstream_distance(usgs_id, log)
        dist   = result["river_dist_km"]
        cached[site_id] = dist

        sites.loc[sites["MonitoringLocationIdentifier"] == site_id,
                  "river_dist_km"]     = dist
        sites.loc[sites["MonitoringLocationIdentifier"] == site_id,
                  "tidal_flag"]        = result["tidal_flag"]
        sites.loc[sites["MonitoringLocationIdentifier"] == site_id,
                  "nldi_error"]        = result.get("nldi_error")
        sites.loc[sites["MonitoringLocationIdentifier"] == site_id,
                  "n_downstream_segs"] = result["n_downstream_segs"]

        new_rows.append({"site": site_id, "river_dist_km": dist})

    # Update cache
    if new_rows:
        new_cache = pd.DataFrame(new_rows)
        if cache_path.exists():
            old_cache = pd.read_csv(cache_path)
            new_cache = pd.concat([old_cache, new_cache]).drop_duplicates(
                subset="site", keep="last"
            )
        new_cache.to_csv(cache_path, index=False)
        log.info(f"  Cached {len(new_rows)} new NLDI distances to {cache_path.name}")

    # Apply distance filter
    n_before = len(sites)
    sites    = sites[
        sites["river_dist_km"].isna() |
        (sites["river_dist_km"] <= cfg["max_river_dist_km"])
    ].copy()
    log.info(
        f"  River distance filter (<{cfg['max_river_dist_km']} km): "
        f"{len(sites):,} sites (dropped {n_before - len(sites):,})"
    )
    return sites.reset_index(drop=True)


# =============================================================================
# STEP 5 — DATA ADEQUACY FILTER
# =============================================================================

def check_data_adequacy(
    sites: pd.DataFrame, cfg: dict, log: logging.Logger
) -> pd.DataFrame:
    """
    Filter sites by data adequacy using observation counts already retrieved
    in query_candidate_sites() — no additional API calls required.

    The ph_result_count and alk_result_count columns were populated by
    wqp.what_sites() in the discovery step. Using them here means data
    adequacy is checked in zero additional API calls rather than one per site.

    Record length is estimated from WQP start/end dates if available in
    the site metadata, or from NWIS get_info() for sites where it is not.
    NWIS calls are batched (one call per 500 sites) rather than per-site.
    """
    log.info("Applying data adequacy filters...")
    log.info("  Using observation counts from discovery step — no additional API calls")

    # Rename pre-computed columns to match expected names
    sites = sites.copy()
    if "ph_result_count" in sites.columns:
        sites["ph_obs"]  = pd.to_numeric(sites["ph_result_count"],  errors="coerce").fillna(0)
        sites["alk_obs"] = pd.to_numeric(sites["alk_result_count"], errors="coerce").fillna(0)
    else:
        # Fallback if counts weren't populated (shouldn't happen)
        sites["ph_obs"]  = 0
        sites["alk_obs"] = 0
        log.warning("  result count columns missing — all sites will pass obs filter")

    # Estimate record length from WQP site metadata fields if present
    for start_col in ["activityStartDate", "MonitoringLocationStartDate", "StartDate"]:
        if start_col in sites.columns:
            sites["record_start"] = pd.to_datetime(
                sites[start_col], errors="coerce"
            )
            break

    for end_col in ["activityEndDate", "MonitoringLocationEndDate", "EndDate"]:
        if end_col in sites.columns:
            sites["record_end"] = pd.to_datetime(
                sites[end_col], errors="coerce"
            )
            break

    if "record_start" in sites.columns and "record_end" in sites.columns:
        sites["record_yrs"] = (
            (sites["record_end"] - sites["record_start"]).dt.days / 365.25
        )
    else:
        sites["record_yrs"] = np.nan

    # For sites missing record length, batch-fetch from NWIS (one call per 500)
    missing_mask = sites["record_yrs"].isna()
    if missing_mask.sum() > 0:
        log.info(
            f"  Fetching record lengths from NWIS for "
            f"{missing_mask.sum()} sites without WQP dates..."
        )
        site_ids = sites.loc[missing_mask, "MonitoringLocationIdentifier"].tolist()
        site_nos = [s.replace("USGS-", "") for s in site_ids]

        record_map = {}
        for i in range(0, len(site_nos), 500):
            batch = site_nos[i:i+500]
            try:
                info, _ = nwis.get_info(sites=batch)
                if info is not None and not info.empty:
                    for _, row in info.iterrows():
                        try:
                            begin = pd.to_datetime(row.get("begin_date"), errors="coerce")
                            end   = pd.to_datetime(row.get("end_date"),   errors="coerce")
                            if pd.notna(begin) and pd.notna(end):
                                yrs = (end - begin).days / 365.25
                                record_map[f"USGS-{row['site_no']}"] = yrs
                        except Exception:
                            pass
            except Exception as e:
                log.debug(f"  NWIS batch record length: {e}")

        for sid, yrs in record_map.items():
            sites.loc[
                sites["MonitoringLocationIdentifier"] == sid, "record_yrs"
            ] = yrs

    # Apply filters
    n_before  = len(sites)
    meets_ph  = sites["ph_obs"]  >= cfg["min_ph_obs"]
    meets_alk = sites["alk_obs"] >= cfg["min_alk_obs"]
    meets_rec = sites["record_yrs"].isna() | (
        sites["record_yrs"] >= cfg["min_record_yrs"]
    )
    sites = sites[meets_ph & meets_alk & meets_rec].copy()

    log.info(f"  Passed adequacy filter: {len(sites):,} of {n_before:,} sites")
    log.info(
        f"    Criteria: >= {cfg['min_ph_obs']} pH obs, "
        f">= {cfg['min_alk_obs']} ALK obs, "
        f">= {cfg['min_record_yrs']} yr record"
    )
    return sites.reset_index(drop=True)


# =============================================================================
# STEP 6 — SCORE AND RANK
# =============================================================================

def score_sites(sites: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Score each site 0–1 on each criterion and compute a weighted total.

    Proximity score: 1.0 at coast, 0.0 at max_river_dist_km
    Observation scores: logarithmic — diminishing returns above ~500 obs
    Record score: 1.0 at 50+ years, 0.0 at min_record_yrs

    These scores are for ranking only, not scientific filtering.
    All sites in the output have already passed the hard filters.
    """
    sites = sites.copy()
    w     = cfg["score_weights"]

    # Proximity score — river distance to coast
    max_d = cfg["max_river_dist_km"]
    sites["score_proximity"] = sites["river_dist_km"].apply(
        lambda d: max(0.0, 1.0 - d / max_d) if pd.notna(d) else 0.5
    )

    # Observation scores — log scale, capped at 1000
    for col, score_col in [("ph_obs", "score_ph"), ("alk_obs", "score_alk")]:
        min_obs = cfg.get(f"min_{col.replace('_obs', '')}_obs", 25)
        sites[score_col] = sites[col].apply(
            lambda n: min(1.0, np.log(max(n, 1) / min_obs) / np.log(1000 / min_obs))
            if pd.notna(n) else 0.0
        )

    # Record length score
    min_yr = cfg["min_record_yrs"]
    sites["score_record"] = sites["record_yrs"].apply(
        lambda y: min(1.0, (y - min_yr) / (50 - min_yr))
        if pd.notna(y) and y >= min_yr else 0.0
    )

    # Weighted total
    sites["score_total"] = (
        w["proximity"]  * sites["score_proximity"]
        + w["ph_obs"]   * sites["score_ph"]
        + w["alk_obs"]  * sites["score_alk"]
        + w["record_yrs"] * sites["score_record"]
    )

    # Assign tiers
    sites["selection_tier"] = pd.cut(
        sites["score_total"],
        bins=[-0.01, cfg["tier3_score"], cfg["tier2_score"], cfg["tier1_score"], 1.01],
        labels=["4 — marginal", "3 — good", "2 — excellent", "1 — top"],
    )

    return sites.sort_values("score_total", ascending=False).reset_index(drop=True)


# =============================================================================
# STEP 7 — COMPARE TO EXISTING LIST
# =============================================================================

def compare_to_existing(
    scored: pd.DataFrame,
    existing_path: Path,
    cfg: dict,
    log: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compare systematic list against manual list.
    Returns (full_scored, sites_to_add, sites_to_remove).
    """
    if not existing_path.exists():
        log.warning(f"Existing site list not found: {existing_path}")
        return scored, scored, pd.DataFrame()

    existing = pd.read_excel(existing_path, header=None)
    existing_ids = set(existing.iloc[:, 0].dropna().astype(str).tolist())
    systematic_ids = set(scored["MonitoringLocationIdentifier"].tolist())

    # Sites in systematic list but not in manual list
    to_add = scored[
        ~scored["MonitoringLocationIdentifier"].isin(existing_ids)
    ].copy()
    to_add["comparison"] = "NEW — in systematic, not in manual list"

    # Sites in manual list but not in systematic list
    to_remove_ids = existing_ids - systematic_ids
    to_remove     = pd.DataFrame({
        "MonitoringLocationIdentifier": list(to_remove_ids),
        "comparison": "REVIEW — in manual list, not in systematic list",
    })

    # Sites in both
    in_both = existing_ids & systematic_ids

    log.info(f"\n  Comparison with existing manual list:")
    log.info(f"    Manual list:              {len(existing_ids)} sites")
    log.info(f"    Systematic list:          {len(systematic_ids)} sites")
    log.info(f"    In both:                  {len(in_both)} sites")
    log.info(f"    New (add to analysis):    {len(to_add)} sites")
    log.info(f"    Review (not in systematic): {len(to_remove_ids)} sites")

    # Mark in full scored table
    scored["in_manual_list"] = scored["MonitoringLocationIdentifier"].isin(
        existing_ids
    )

    return scored, to_add, to_remove


# =============================================================================
# MAIN
# =============================================================================

def main():
    CFG["output_dir"].mkdir(parents=True, exist_ok=True)
    log = setup_logging(CFG["output_dir"])

    log.info("=" * 65)
    log.info("  NEAA Systematic Site Selection")
    log.info(f"  {datetime.now():%Y-%m-%d %H:%M}")
    log.info("=" * 65)
    log.info("")
    log.info("  Filters applied in sequence:")
    log.info(f"    1. USGS stream sites with pH + alkalinity data")
    log.info(f"    2. Non-tidal freshwater stream type")
    log.info(f"    3. Within {CFG['max_straight_dist_km']} km of coast (straight-line pre-filter)")
    log.info(f"    4. Within {CFG['max_river_dist_km']} km river-path to ocean (NLDI)")
    log.info(f"    5. >= {CFG['min_ph_obs']} pH obs, >= {CFG['min_alk_obs']} ALK obs")
    log.info(f"    6. >= {CFG['min_record_yrs']} year record")
    log.info("")

    # Step 1 — query all candidates
    candidates = query_candidate_sites(CFG, log)

    # Step 2 — filter site types
    filtered = filter_site_types(candidates, CFG, log)

    # Step 3 — spatial pre-filter (fast)
    pre_filtered = spatial_prefilter(filtered, CFG, log)

    # Step 4 — Data adequacy (uses counts from Step 1 — zero extra API calls)
    # Applied BEFORE NLDI to minimise expensive river-distance queries.
    # Only sites with enough data proceed to the NLDI step.
    log.info("\nStep 4: Applying data adequacy filter...")
    adequate = check_data_adequacy(pre_filtered, CFG, log)

    # Step 5 — NLDI river distance (slow — cached after first run)
    # Now only runs for sites that passed all other filters.
    log.info(f"\nStep 5: Computing river-path distances via NLDI...")
    log.info(f"  Running for {len(adequate):,} sites (post-filter, not full candidate set)")
    log.info("  NOTE: First run takes ~2 sec/site. Results cached for re-runs.")
    with_distances = compute_river_distances(adequate, CFG, log)

    # Step 6 — score and rank
    log.info("\nStep 6: Scoring and ranking sites...")
    scored = score_sites(adequate, CFG)

    # Step 7 — compare to existing list
    log.info("\nStep 7: Comparing to existing manual site list...")
    existing_path = CFG["project_dir"] / CFG["existing_sites"]
    scored, to_add, to_remove = compare_to_existing(
        scored, existing_path, CFG, log
    )

    # ── Write outputs ─────────────────────────────────────────────────────────
    out = CFG["output_dir"]

    scored.to_csv(out / "site_selection_results.csv", index=False)
    if not to_add.empty:
        to_add.to_csv(out / "site_selection_added.csv", index=False)
    if not to_remove.empty:
        to_remove.to_csv(out / "site_selection_removed.csv", index=False)

    # Text report
    report_path = out / "site_selection_report.txt"
    with open(report_path, "w") as f:
        f.write(f"NEAA Systematic Site Selection Report\n")
        f.write(f"Generated: {datetime.now():%Y-%m-%d %H:%M}\n\n")
        f.write(f"Selection criteria:\n")
        f.write(f"  USGS freshwater non-tidal stream sites\n")
        f.write(f"  Coastal states: {', '.join(CFG['coastal_states'])}\n")
        f.write(f"  River distance to coast: <= {CFG['max_river_dist_km']} km\n")
        f.write(f"  Min pH observations: {CFG['min_ph_obs']}\n")
        f.write(f"  Min alkalinity observations: {CFG['min_alk_obs']}\n")
        f.write(f"  Min record length: {CFG['min_record_yrs']} years\n\n")
        f.write(f"Results:\n")
        f.write(f"  Total sites selected: {len(scored)}\n")
        tier_counts = scored["selection_tier"].value_counts().sort_index()
        for tier, count in tier_counts.items():
            f.write(f"  {tier}: {count} sites\n")
        f.write(f"\n  Sites to add to analysis: {len(to_add)}\n")
        f.write(f"  Manual sites to review: {len(to_remove)}\n\n")

        f.write("Top 30 sites by score:\n")
        cols = [
            "MonitoringLocationIdentifier", "MonitoringLocationName",
            "river_dist_km", "ph_obs", "alk_obs", "record_yrs",
            "score_total", "selection_tier", "in_manual_list",
        ]
        avail = [c for c in cols if c in scored.columns]
        f.write(scored[avail].head(30).to_string(index=False))
        f.write("\n")

    # ── Console summary ───────────────────────────────────────────────────────
    log.info("")
    log.info("=" * 65)
    log.info("  SITE SELECTION COMPLETE")
    log.info("=" * 65)
    log.info(f"  Total sites selected: {len(scored)}")
    for tier, count in scored["selection_tier"].value_counts().sort_index().items():
        log.info(f"    {tier}: {count}")
    log.info(f"  New sites to add:     {len(to_add)}")
    log.info(f"  Manual sites to review: {len(to_remove)}")
    log.info("")
    log.info(f"  Outputs written to: {out}")
    log.info(f"    site_selection_results.csv")
    log.info(f"    site_selection_added.csv")
    log.info(f"    site_selection_removed.csv")
    log.info(f"    site_selection_report.txt")
    log.info("=" * 65)


if __name__ == "__main__":
    main()
