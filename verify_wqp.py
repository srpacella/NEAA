"""
NEAA WQP API Verification Script
=================================
Run this BEFORE running the full pipeline to confirm that:
  1. The dataretrieval package is installed and at the right version
  2. The WQP API endpoint is reachable and returning data
  3. The column names returned by the API match what the pipeline expects
  4. NWIS site metadata download works
  5. A small real download (one site, one characteristic) works end-to-end

Usage:
    python verify_wqp.py

Expected output: a series of PASS / FAIL / WARN lines and a final summary.
Any FAIL must be fixed before running neaa_pipeline.py.
Any WARN should be reviewed — the pipeline may need a column-name update.
"""

import sys
import importlib
import traceback
from datetime import datetime

# ── Formatting helpers ────────────────────────────────────────────────────────
WIDTH = 70
PASS  = "\033[92m  PASS\033[0m"
FAIL  = "\033[91m  FAIL\033[0m"
WARN  = "\033[93m  WARN\033[0m"
INFO  = "  INFO"

results = []

def check(label: str, passed: bool, detail: str = "", warn: bool = False) -> None:
    tag = WARN if (warn and not passed) else (PASS if passed else FAIL)
    print(f"{tag}  {label}")
    if detail:
        for line in detail.splitlines():
            print(f"        {line}")
    results.append(("WARN" if (warn and not passed) else ("PASS" if passed else "FAIL"), label))

def section(title: str) -> None:
    print(f"\n{'─' * WIDTH}")
    print(f"  {title}")
    print(f"{'─' * WIDTH}")


# =============================================================================
# CHECK 1 — Python version
# =============================================================================
section("1. Python environment")

py_ok = sys.version_info >= (3, 9)
check(
    f"Python >= 3.9  (found {sys.version.split()[0]})",
    py_ok,
    "" if py_ok else "Upgrade Python. The pipeline uses match/case and union types.",
)

# =============================================================================
# CHECK 2 — Required packages
# =============================================================================
section("2. Required packages")

REQUIRED = {
    "dataretrieval": "0.9",
    "pandas":        "1.5",
    "numpy":         "1.23",
    "scipy":         "1.9",
    "PyCO2SYS":      "1.8",
    "geopy":         "2.0",
    "h5py":          "3.0",
    "requests":      "2.28",
    "tqdm":          "4.0",
    "pymannkendall": "1.4",
    "openpyxl":      "3.0",
}

pkg_versions = {}
for pkg, min_ver in REQUIRED.items():
    try:
        mod = importlib.import_module(pkg)
        ver = getattr(mod, "__version__", "unknown")
        pkg_versions[pkg] = ver

        if ver != "unknown":
            from packaging.version import Version
            ok = Version(ver) >= Version(min_ver)
        else:
            ok = True  # can't compare, assume ok

        check(f"{pkg} >= {min_ver}  (found {ver})", ok,
              f"pip install --upgrade {pkg}" if not ok else "")
    except ImportError:
        check(f"{pkg} (not installed)", False,
              f"pip install {pkg}")
        pkg_versions[pkg] = None


# =============================================================================
# CHECK 3 — dataretrieval API shape
# =============================================================================
section("3. dataretrieval package internals")

try:
    import dataretrieval.wqp as wqp
    import dataretrieval.nwis as nwis

    # Check that the functions we call actually exist
    check("wqp.get_results() exists",   hasattr(wqp,  "get_results"))
    check("nwis.get_info() exists",     hasattr(nwis, "get_info"))
    check("nwis.get_record() exists",   hasattr(nwis, "get_record"))

    # Check dataretrieval version specifically — API changed in v0.9
    dr_ver = pkg_versions.get("dataretrieval", "0")
    if dr_ver:
        from packaging.version import Version
        check(
            "dataretrieval >= 0.9 (new WQP3 API support)",
            Version(dr_ver) >= Version("0.9"),
            "Versions < 0.9 use the old WQP endpoint which may be deprecated.\n"
            "Run: pip install --upgrade dataretrieval",
            warn=True,
        )

except Exception as e:
    check("dataretrieval import", False, str(e))


# =============================================================================
# CHECK 4 — Live WQP API call (small test)
# =============================================================================
section("4. Live WQP API — small test download")
print(f"  Testing with USGS-01010000 (Aroostook River, ME) — pH only\n")

wqp_df = None
try:
    import dataretrieval.wqp as wqp
    wqp_df, meta = wqp.get_results(
        siteid="USGS-01010000",
        characteristicName="pH",
    )
    n_rows = len(wqp_df) if wqp_df is not None else 0
    check(f"WQP download returned data ({n_rows} rows)", n_rows > 0,
          "Zero rows returned — API may be down or site ID format changed.")
except Exception as e:
    check("WQP API call succeeded", False, traceback.format_exc(limit=3))


# =============================================================================
# CHECK 5 — Column names (the most likely thing to have changed)
# =============================================================================
section("5. WQP column names (API v3 vs v2 differences)")

# These are the columns the pipeline reads. If any are missing the pipeline
# will fail silently or KeyError. We check for both old and new name variants.
REQUIRED_COLS = {
    # Column the pipeline uses           Possible alternative names in newer API
    "ActivityStartDate":                ["ActivityStartDate"],
    "MonitoringLocationIdentifier":     ["MonitoringLocationIdentifier"],
    "CharacteristicName":               ["CharacteristicName"],
    "ResultMeasureValue":               ["ResultMeasureValue", "ResultMeasure/MeasureValue"],
    "ActivityMediaName":                ["ActivityMediaName"],
    "ActivityStartTime/Time":           ["ActivityStartTime/Time", "ActivityStartTime.Time"],
    "USGSPCode":                        ["USGSPCode", "ResultAnalyticalMethod/MethodIdentifier"],
    "OrganizationIdentifier":           ["OrganizationIdentifier"],
}

if wqp_df is not None and not wqp_df.empty:
    print(f"  Columns returned by API ({len(wqp_df.columns)} total):")
    for col in sorted(wqp_df.columns):
        print(f"    {col}")
    print()

    for pipeline_col, alternatives in REQUIRED_COLS.items():
        found = any(alt in wqp_df.columns for alt in alternatives)
        actual = next((alt for alt in alternatives if alt in wqp_df.columns), None)

        if found and actual == pipeline_col:
            check(f"'{pipeline_col}' present", True)
        elif found and actual != pipeline_col:
            check(
                f"'{pipeline_col}' — found as '{actual}' (renamed!)",
                False,
                f"Pipeline uses '{pipeline_col}' but API now returns '{actual}'.\n"
                f"Update the pipeline column reference.",
                warn=True,
            )
        else:
            check(
                f"'{pipeline_col}' present",
                False,
                f"Column not found. Available similar columns:\n"
                + "\n".join(f"  {c}" for c in wqp_df.columns if pipeline_col[:8].lower() in c.lower()),
            )
else:
    print("  Skipping column checks — no data downloaded.")


# =============================================================================
# CHECK 6 — Multi-characteristic download (as the pipeline does it)
# =============================================================================
section("6. Multi-characteristic download test")

multi_df = None
try:
    import dataretrieval.wqp as wqp
    multi_df, _ = wqp.get_results(
        siteid="USGS-01010000",
        characteristicName=["pH", "Alkalinity", "Temperature, water", "Specific conductance"],
    )
    n = len(multi_df) if multi_df is not None else 0
    check(f"Multi-characteristic download ({n} rows)", n > 0)

    if multi_df is not None and "CharacteristicName" in multi_df.columns:
        chars_found = multi_df["CharacteristicName"].unique()
        for c in ["pH", "Alkalinity", "Temperature, water", "Specific conductance"]:
            check(f"  Characteristic '{c}' in results", c in chars_found,
                  warn=True)
except Exception as e:
    check("Multi-characteristic download", False, traceback.format_exc(limit=3))


# =============================================================================
# CHECK 7 — NWIS site metadata
# =============================================================================
section("7. NWIS site metadata download")

try:
    import dataretrieval.nwis as nwis
    meta_df, _ = nwis.get_info(sites=["01010000"])
    check("nwis.get_info() returned data", len(meta_df) > 0)

    for col in ["site_no", "dec_lat_va", "dec_long_va", "station_nm"]:
        check(f"  Metadata column '{col}' present", col in meta_df.columns, warn=True)
except Exception as e:
    check("NWIS metadata download", False, traceback.format_exc(limit=3))


# =============================================================================
# CHECK 8 — WQP site query (whatWQPdata equivalent)
# =============================================================================
section("8. WQP site metadata query")

try:
    import dataretrieval.wqp as wqp
    # whatWQPdata() is now wqp.what_sites() in newer dataretrieval
    if hasattr(wqp, "what_sites"):
        sites_df, _ = wqp.what_sites(siteid="USGS-01010000")
        check("wqp.what_sites() works (new API)", len(sites_df) > 0)
        for col in ["lon", "lat"]:
            check(f"  Site column '{col}' present", col in sites_df.columns, warn=True)
    elif hasattr(wqp, "whatWQPdata"):
        check("wqp.whatWQPdata() exists (old API)", True,
              "Consider upgrading dataretrieval — whatWQPdata may be deprecated.", warn=True)
    else:
        check("Site query function exists", False,
              "Neither what_sites() nor whatWQPdata() found in dataretrieval.wqp")
except Exception as e:
    check("WQP site query", False, traceback.format_exc(limit=3))


# =============================================================================
# SUMMARY
# =============================================================================
section("SUMMARY")

n_pass = sum(1 for r, _ in results if r == "PASS")
n_warn = sum(1 for r, _ in results if r == "WARN")
n_fail = sum(1 for r, _ in results if r == "FAIL")

print(f"\n  {n_pass} passed  |  {n_warn} warnings  |  {n_fail} failed\n")

if n_fail == 0 and n_warn == 0:
    print("  ✓ All checks passed. Pipeline should run correctly.")
elif n_fail == 0:
    print("  ⚠ No hard failures, but review warnings above.")
    print("    Column renames or API changes may require small pipeline edits.")
else:
    print("  ✗ Fix all FAILs before running the pipeline.")
    print("\n  Failed checks:")
    for r, label in results:
        if r == "FAIL":
            print(f"    - {label}")

print(f"\n  Run timestamp: {datetime.now():%Y-%m-%d %H:%M:%S}")
print()
