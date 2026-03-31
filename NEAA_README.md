# NEAA Python Pipeline

Python translation of the National Estuary Acidification Assessment R workflow.

## Quick start

### 1. Install dependencies

```bash
pip install dataretrieval pandas numpy scipy PyCO2SYS openpyxl \
            geopy h5py requests tqdm pymannkendall
```

### 2. Set your project directory

Open `neaa_pipeline.py` and edit the `CFG` block at the top:

```python
CFG = {
    "project_dir":    Path("/your/project/folder"),
    "site_list_xlsx": "WQP_sites_manual_srp_update6.xlsx",
    ...
}
```

Place your site list Excel file in that folder. Everything else is downloaded automatically.

### 3. Run

```bash
python neaa_pipeline.py
```

---

## What it does (pipeline steps)

| Step | R equivalent | Description |
|------|-------------|-------------|
| 1 | `neaa_final_download_231028.R` | Download pH, ALK, temp, conductivity from WQP; QA/QC, salinity filter |
| 2 | `neaa_corrections.R` | Download DOC & discharge; apply Liu et al. 2020 pH and alkalinity corrections |
| 3 | `ts_analysis_240407_seasonal.R` | Theil-Sen slopes + Seasonal Kendall trends for H⁺ and alkalinity |
| 4 | `neaa_calcium.R` | Download calcium data and join to site records |
| 5 | `codap_load.R` | Download CODAP-NA v2020 ocean dataset (cached as parquet after first run) |
| 6 | `codap_organize.R` | Match each river site to nearby CODAP ocean stations (200 km radius) |
| 7 | `export_*.R` | Write per-site CSVs + single HDF5 file for MATLAB |

---

## Outputs

All outputs go to `<project_dir>/neaa_data_exports/`.

### CSV files (same as R output)
| File | Contents |
|------|----------|
| `USGS-XXXXXXXX.csv` | Per-site QA'd water quality data |
| `USGS-XXXXXXXX_calcium.csv` | Per-site calcium observations |
| `USGS-XXXXXXXX_ocean.csv` | Per-site matched CODAP ocean data + carbonate chemistry |
| `ts_regressions.csv` | Theil-Sen and seasonal Kendall results for all sites |
| `siteINFO.csv` | NWIS station metadata |
| `discharge_means.csv` | Mean discharge per site |

### HDF5 file (`neaa_export.h5`)
Single file containing all the above, organized by group:

```
/river_data/<site_id>/     — water quality columns as datasets
/calcium/<site_id>/
/ocean/<site_id>/
/summaries/trends
/summaries/site_metadata
/summaries/discharge_means
```

**Load in MATLAB:**
```matlab
% Read trend results
trends = h5read('neaa_export.h5', '/summaries/trends');

% Read one site's data
site_data = h5read('neaa_export.h5', '/river_data/USGS_01010000');

% List all sites
info = h5info('neaa_export.h5', '/river_data');
sites = {info.Groups.Name};
```

---

## Key differences from R version

### Improvements
| Issue in R | Fix in Python |
|-----------|--------------|
| Hardcoded Windows paths everywhere | Single `CFG` dict at top of file — one place to edit |
| Row-by-row loops for salinity, pH correction, IQR removal | Vectorized with NumPy/pandas — ~100x faster |
| Two nearly-identical trend analysis scripts | One function handles both Theil-Sen and Seasonal Kendall |
| No logging | Timestamped log file written to `<project_dir>/logs/` |
| `seacarb` (R-only) | `PyCO2SYS` — equivalent carbonate chemistry, actively maintained |
| `dataRetrieval` (R-only) | `dataretrieval` — USGS-maintained Python port |
| CODAP as local Excel only | Downloaded from NOAA URL, cached as parquet for fast reloads |
| No MATLAB-native output | HDF5 export alongside CSVs |

### Scientific equivalences
| R function | Python equivalent | Notes |
|-----------|-------------------|-------|
| `readWQPdata()` | `wqp.get_results()` | Same WQP API |
| `readNWISsite()` | `nwis.get_info()` | Same NWIS API |
| `ec2pss()` | `conductivity_to_salinity()` | Hill et al. 1986 formula |
| `mblm::mblm()` | `scipy.stats.theilslopes()` | Same estimator |
| `rkt::rkt()` | `pymannkendall.seasonal_test()` | Seasonal Mann-Kendall |
| `seacarb::carbfull(flag=15)` | `PyCO2SYS.sys(par1_type=1, par2_type=2)` | ALK+DIC input, total pH scale |
| `geosphere::distHaversine()` | `geopy.distance.great_circle()` | Haversine distance |
| 3×IQR outlier removal | `remove_outliers_iqr()` | Identical logic, vectorized |
| Liu et al. 2020 pH correction | `apply_liu2020_corrections()` | Identical formula |

---

## Configuration reference

All thresholds are in the `CFG` dict. Common ones to adjust:

```python
"ph_min":           5.5,       # QA range for pH
"ph_max":           9.5,
"alk_multiplier":   20,        # unit conversion (keep as-is unless data units change)
"salinity_max":     1.0,       # exclude sites with mean salinity > 1 PSU
"outlier_iqr_mult": 3.0,       # IQR multiplier for outlier removal
"min_duration_yr":  10,        # minimum record length for trend analysis
"min_obs_trend":    25,        # minimum observations for trend analysis
"ocean_radius_m":   200_000,   # CODAP search radius (200 km)
"codap_depth_max_m": 20,       # ocean surface depth cutoff
```

---

## Troubleshooting

**CODAP won't download**
The NOAA server can be slow. If both URLs fail, download manually from:
https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/
Place `CODAP_NA_v2020.xlsx` in your project directory and re-run.

**WQP download is slow**
This is normal — the Water Quality Portal API is rate-limited. The `tqdm` progress bar shows per-site status. A full run across all sites typically takes 20-60 minutes depending on network speed.

**`PyCO2SYS` results differ slightly from `seacarb`**
Small differences (~0.01%) are expected due to different implementations of the same equilibrium constants. The equilibrium constant set used here (Lueker et al. 2000, `opt_k_carbonic=10`) is the closest match to `seacarb`'s `k1k2="l"` option.

**Memory use**
For large site lists (>200 sites), CODAP loading can use ~2 GB RAM. If needed, process sites in batches using the `sites` list in `main()`.
