# Urban sunlight and wearable-derived stress: analysis scripts

This repository contains **code only** for the final study analyses. Raw wearable records, city-day stress counts, weather/source datasets, derived city-level coefficients, and model-output files are not included. They must be handled under the applicable data-use agreements and research approvals.

## Scripts

- `01_absolute_sunlight_dlnm.R`: city-level and pooled absolute sunlight DLNM.
- `02_rolling_mean_sunlight_dlnm.R`: five-day rolling-mean sunlight DLNM.
- `03_low_sunlight_days_dlnm.R`: recent low-sunlight-day count DLNM.
- `04_sunlight_anomaly_meta_regression.R`: within-city anomaly DLNM and cross-city meta-regression.
- `05_city_sunlight_constraint_clustering.R`: three-feature K-means, K=4 city typology.

The scripts retain their original statistical specifications, contrast directions, and model settings. 
## Input configuration

By default, each model looks for restricted-use inputs in `data/restricted` relative to the working directory. Set `SUNLIGHT_PROJECT_ROOT` (project root), `SUNLIGHT_DATA_DIR` (restricted input directory), and `SUNLIGHT_OUTPUT_DIR` (output directory) to override this. Individual file paths may be set using `CITY_DAY_PATH`, `WEATHER_PATH`, `SUNLIGHT_PATH`, `CITY_CLIMATE_PATH`, or `CITY_FEATURE_FILE`.

Expected files (when using defaults):

- `data/restricted/city_day_high_stress_from_personal_p90.csv`
- `data/restricted/weather_sunlight_cleaned_city_day.csv`
- `data/restricted/all_cities_daily_sunhours_cloud_aerosol_raw_totexttau_with_weather_20260720.csv`
- `data/restricted/city_centroid_latitudes_with_ERA5Land_Tmean_2021_2024.csv`
- `data/restricted/city_sunlight_constraint_features.csv`

Inputs are not included and outputs are written to `outputs/` by default. 

## Software

R packages: `data.table`, `dlnm`, `splines`, `mixmeta`; for city clustering: `cluster`, `dplyr`, `readr`. Specific R and package versions should be documented using the environment used for the final analysis.

## Publication note

Readers can inspect the analysis methods, but reproducing numerical results requires authorized access to the underlying datasets. Refer to the manuscript's Data Availability and Code Availability sections for access conditions.
