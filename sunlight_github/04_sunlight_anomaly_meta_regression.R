# Combined within-city and cross-city sunlight anomaly analysis with a shared
# pandemic-adjusted two-stage DLNM first stage
#
# Purpose
# -------
# This script uses one shared first-stage DLNM for within-city sunlight
# anomalies, then writes two second-stage result branches:
#   1. within/: intercept-only random-effects meta-analysis of the national
#      within-city anomaly-response curve.
#   2. cross/: meta-regression asking whether the anomaly-response curve differs
#      by each city's multi-year mean observed sunlight duration, adjusted for
#      city centroid latitude and ERA5-Land annual mean temperature.
#
# Stage 1: city-specific pandemic-adjusted quasibinomial DLNMs for seasonally
#          normalized sunlight anomalies, centered at 0 hours/day.
# Stage 2 within: multivariate random-effects meta-analysis of the reduced city
#                 curves with intercept only.
# Stage 2 cross:  multivariate random-effects meta-regression of the reduced city
#                 curves on each city's multi-year mean observed sunlight duration,
#                 city centroid latitude, and ERA5-Land annual mean temperature.
#
# Important interpretation
# ------------------------
# The city-level modifier is the unweighted mean of observed daily sunlight
# duration over the study period. It is not weighted by wearable sample size.
# City centroid latitude and ERA5-Land annual mean temperature are included as
# city-level covariates in the cross branch. Other city-level socioeconomic,
# demographic, climate, or healthcare covariates are not included. Daily
# meteorological controls remain in the city-specific models. Results therefore
# assess latitude- and temperature-adjusted cross-city effect modification and
# should not be described as fully adjusted causal effects.

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(splines)
  library(mixmeta)
})

options(stringsAsFactors = FALSE)
setDTthreads(max(1L, parallel::detectCores(logical = TRUE) - 1L))

# -----------------------------------------------------------------------------
# 1. Paths and analysis settings
# -----------------------------------------------------------------------------

# Input files are restricted-use and must not be committed to a public repository.
# Supply the project root or override individual locations with environment variables.
PROJECT_ROOT <- Sys.getenv("SUNLIGHT_PROJECT_ROOT", unset = ".")
DATA_DIR <- Sys.getenv("SUNLIGHT_DATA_DIR", unset = file.path(PROJECT_ROOT, "data", "restricted"))
OUTPUT_ROOT <- Sys.getenv("SUNLIGHT_OUTPUT_DIR", unset = file.path(PROJECT_ROOT, "outputs"))
CITY_DAY_PATH <- Sys.getenv("CITY_DAY_PATH", unset = file.path(DATA_DIR, "city_day_high_stress_from_personal_p90.csv"))
WEATHER_PATH <- Sys.getenv("WEATHER_PATH", unset = file.path(DATA_DIR, "weather_sunlight_cleaned_city_day.csv"))
SUNLIGHT_PATH <- Sys.getenv("SUNLIGHT_PATH", unset = file.path(DATA_DIR, "all_cities_daily_sunhours_cloud_aerosol_raw_totexttau_with_weather_20260720.csv"))
CITY_CLIMATE_PATH <- Sys.getenv("CITY_CLIMATE_PATH", unset = file.path(DATA_DIR, "city_centroid_latitudes_with_ERA5Land_Tmean_2021_2024.csv"))

STUDY_START <- as.IDate("2021-12-01")
STUDY_END <- as.IDate("2024-12-31")

# The first command-line argument sets the sunlight lag. For example:
# Rscript 04_sunlight_anomaly_meta_regression.R 28
cli_args <- commandArgs(trailingOnly = TRUE)
SUN_MAX_LAG <- if (length(cli_args) >= 1L) {
  suppressWarnings(as.integer(cli_args[1L]))
} else {
  28L
}
if (!is.finite(SUN_MAX_LAG) || SUN_MAX_LAG < 1L) {
  stop("The sunlight maximum lag must be a positive integer.")
}

# Temperature is controlled only on the same day using a nonlinear natural
# spline. No lagged temperature terms are constructed.
TEMP_MAX_LAG <- 0L

SUN_LAG_KNOTS <- c(3, 7, 14, 21)
N_SUN_LAG_KNOTS <- length(SUN_LAG_KNOTS)
SUN_VAR_KNOT_PROBS <- c(0.25, 0.50, 0.75)
TEMP_VAR_KNOT_PROBS <- c(0.25, 0.50, 0.75)

# Leave-one-year-out city-season climatology for defining sunlight anomalies.
CLIM_HALF_WINDOW <- 15L
MIN_CLIMATOLOGY_N <- 30L

TIME_DF_PER_YEAR <- 7L
TIME_DF <- 21L

MIN_MODEL_DAYS <- 365L
MIN_TOTAL_EVENTS <- 20L
MIN_TOTAL_OBSERVED <- 100L

# Quality criteria for the multi-year baseline sunlight metric. These are data
# reliability filters, not city-level confounder adjustments.
MIN_BASELINE_YEARS <- 2L
MIN_BASELINE_COVERAGE <- 0.80

CAPACITY_PROBS <- c(P10 = 0.10, P50 = 0.50, P90 = 0.90)
TARGET_ANOMALY_HOURS <- -2.768

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_ROOT <- file.path(
  OUTPUT_ROOT,
  paste0("sunlight_anomaly_lag0_", SUN_MAX_LAG, "_", timestamp)
)
WITHIN_DIR <- file.path(OUT_ROOT, "within")
CROSS_DIR <- file.path(OUT_ROOT, "cross")
dir.create(WITHIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CROSS_DIR, recursive = TRUE, showWarnings = FALSE)

# Shared data outputs are written to the within-city and cross-city result directories.
OUT_DIR <- CROSS_DIR
LOG_FILE <- file.path(OUT_ROOT, "run_log.txt")

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

log_msg <- function(...) {
  txt <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    paste(..., collapse = "")
  )
  message(txt)
  cat(txt, "\n", file = LOG_FILE, append = TRUE)
}

save_csv <- function(x, filename) {
  fwrite(as.data.table(x), file.path(OUT_DIR, filename), bom = TRUE)
}

save_csv_to <- function(x, out_dir, filename) {
  fwrite(as.data.table(x), file.path(out_dir, filename), bom = TRUE)
}

save_csv_within <- function(x, filename) {
  save_csv_to(x, WITHIN_DIR, filename)
}

save_csv_cross <- function(x, filename) {
  save_csv_to(x, CROSS_DIR, filename)
}

save_csv_both <- function(x, filename) {
  save_csv_within(x, filename)
  save_csv_cross(x, filename)
}

clean_city_code <- function(x) {
  y <- trimws(as.character(x))
  y <- sub("\\.0+$", "", y)
  y[y %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  y
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_first_character <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) NA_character_ else x[1L]
}

read_city_climate_covariates <- function(path) {
  if (!file.exists(path)) stop("Missing CITY_CLIMATE_PATH: ", path)

  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    raw <- fread(path, showProgress = FALSE)
  } else {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package 'readxl' is required to read the city climate workbook.")
    }
    raw <- as.data.table(readxl::read_excel(path, sheet = 1))
  }
  if (ncol(raw) < 2L) {
    stop("The city climate covariate table does not contain enough columns.")
  }

  code_candidates <- names(raw)[grepl("gb|code|city_code", names(raw), ignore.case = TRUE)]
  code_col <- if (length(code_candidates) > 0L) code_candidates[1L] else names(raw)[1L]

  numeric_like <- names(raw)[vapply(raw, function(x) {
    any(is.finite(suppressWarnings(as.numeric(x))))
  }, logical(1L))]
  lat_name_candidates <- numeric_like[grepl("latitude|lat", numeric_like, ignore.case = TRUE)]
  lat_candidates <- numeric_like[vapply(numeric_like, function(nm) {
    z <- suppressWarnings(as.numeric(raw[[nm]]))
    ok <- is.finite(z)
    sum(ok) > 0L &&
      max(abs(z[ok]), na.rm = TRUE) <= 90 &&
      nm != code_col
  }, logical(1L))]
  lat_col <- if (length(lat_name_candidates) > 0L) {
    lat_name_candidates[1L]
  } else if (length(lat_candidates) > 0L) {
    lat_candidates[1L]
  } else {
    stop("Could not identify a valid latitude column in CITY_CLIMATE_PATH.")
  }

  tmean_candidates <- setdiff(
    numeric_like[grepl("tmean|mean.*temp|temp.*mean|ERA5", numeric_like, ignore.case = TRUE)],
    c(code_col, lat_col)
  )
  if (length(tmean_candidates) == 0L) {
    stop("Could not identify an annual mean temperature column in CITY_CLIMATE_PATH.")
  }
  tmean_col <- tmean_candidates[1L]

  name_candidates <- setdiff(names(raw), c(code_col, lat_col, tmean_col))
  name_col <- if (length(name_candidates) > 0L) name_candidates[1L] else NA_character_

  out <- data.table(
    city_code = clean_city_code(raw[[code_col]]),
    latitude = suppressWarnings(as.numeric(raw[[lat_col]])),
    annual_tmean = suppressWarnings(as.numeric(raw[[tmean_col]]))
  )
  if (!is.na(name_col)) {
    out[, city_name_climate_table := as.character(raw[[name_col]])]
  } else {
    out[, city_name_climate_table := NA_character_]
  }

  out <- out[
    !is.na(city_code) &
      is.finite(latitude) &
      latitude >= -90 &
      latitude <= 90 &
      is.finite(annual_tmean)
  ]
  if (nrow(out) == 0L) {
    stop("No valid city climate rows were found in CITY_CLIMATE_PATH.")
  }

  out[, .(
    latitude = mean(latitude, na.rm = TRUE),
    annual_tmean = mean(annual_tmean, na.rm = TRUE),
    climate_source_rows = .N,
    city_name_climate_table = safe_first_character(city_name_climate_table)
  ), by = city_code]
}

make_pandemic_phase <- function(date) {
  phase <- fifelse(
    date <= as.IDate("2022-11-10"),
    "phase1_dynamic_zero",
    fifelse(
      date <= as.IDate("2022-12-06"),
      "phase2_policy_adjustment",
      fifelse(
        date <= as.IDate("2023-01-22"),
        "phase3_transition_peak",
        "phase4_normalization"
      )
    )
  )
  factor(
    phase,
    levels = c(
      "phase1_dynamic_zero",
      "phase2_policy_adjustment",
      "phase3_transition_peak",
      "phase4_normalization"
    )
  )
}

standardize_with <- function(x, center, scale) {
  if (!is.finite(scale) || scale <= 0) return(rep(NA_real_, length(x)))
  (x - center) / scale
}

make_sig <- function(low, high) {
  fifelse(
    !is.finite(low) | !is.finite(high),
    NA_character_,
    fifelse(
      low > 1,
      "OR > 1 significant",
      fifelse(high < 1, "OR < 1 significant", "not significant")
    )
  )
}

or_p_value <- function(or, low, high) {
  out <- rep(NA_real_, length(or))
  ok <- is.finite(or) & is.finite(low) & is.finite(high) &
    or > 0 & low > 0 & high > 0 & high > low
  se_log <- rep(NA_real_, length(or))
  se_log[ok] <- (log(high[ok]) - log(low[ok])) / (2 * qnorm(0.975))
  z <- log(or[ok]) / se_log[ok]
  out[ok] <- 2 * pnorm(abs(z), lower.tail = FALSE)
  out
}

coef_table_from_model <- function(model, model_name) {
  sm <- summary(model)
  tab <- as.data.table(sm$coefficients, keep.rownames = "term")
  setnames(tab, names(tab), make.names(names(tab)))
  tab[, model := model_name]
  setcolorder(tab, c("model", setdiff(names(tab), "model")))
  tab
}

overall_modifier_test <- function(model, modifier, model_name) {
  b <- coef(model)
  v <- vcov(model)
  idx <- grep(modifier, names(b), fixed = TRUE)
  if (length(idx) == 0L) {
    return(data.table(
      model = model_name,
      modifier = modifier,
      df = 0L,
      chisq = NA_real_,
      p_value = NA_real_
    ))
  }

  bb <- b[idx]
  vv <- v[idx, idx, drop = FALSE]
  inv <- tryCatch(solve(vv), error = function(e) MASS::ginv(vv))
  chisq <- as.numeric(t(bb) %*% inv %*% bb)
  data.table(
    model = model_name,
    modifier = modifier,
    df = length(idx),
    chisq = chisq,
    p_value = pchisq(chisq, df = length(idx), lower.tail = FALSE)
  )
}

make_circular_window_stats <- function(x, half_window) {
  n_days <- length(x)
  stopifnot(n_days == 366L)

  observed <- as.numeric(is.finite(x))
  values <- fifelse(is.finite(x), x, 0)
  values_sq <- values^2
  extend <- function(v) c(tail(v, half_window), v, head(v, half_window))
  width <- 2L * half_window + 1L
  keep <- half_window + seq_len(n_days)

  list(
    sum = frollsum(extend(values), n = width, align = "center")[keep],
    sumsq = frollsum(extend(values_sq), n = width, align = "center")[keep],
    n = frollsum(extend(observed), n = width, align = "center")[keep]
  )
}

contrast_log_or <- function(value, cen, coef_vec, vcov_mat, knots, boundary) {
  basis_values <- onebasis(
    c(value, cen),
    fun = "ns",
    knots = knots,
    Boundary.knots = boundary,
    intercept = FALSE
  )
  contrast <- as.numeric(basis_values[1L, ] - basis_values[2L, ])
  log_or <- sum(contrast * coef_vec)
  var_log_or <- as.numeric(t(contrast) %*% vcov_mat %*% contrast)
  list(log_or = log_or, var_log_or = var_log_or)
}

extract_meta_prediction <- function(x) {
  if (is.list(x) && all(c("fit", "vcov") %in% names(x))) return(x)
  if (is.list(x) && length(x) == 1L && is.list(x[[1L]]) &&
      all(c("fit", "vcov") %in% names(x[[1L]]))) {
    return(x[[1L]])
  }
  stop("Unexpected mixmeta prediction structure.")
}

predict_reduced_coefficients <- function(model, capacity_z, latitude_z, annual_tmean_z, coefficient_names) {
  raw <- predict(
    model,
    newdata = data.frame(
      baseline_actual_sunlight_z = capacity_z,
      latitude_z = latitude_z,
      annual_tmean_z = annual_tmean_z
    ),
    vcov = TRUE,
    format = "list"
  )
  pred <- extract_meta_prediction(raw)
  cf <- as.numeric(pred$fit)
  vc <- as.matrix(pred$vcov)
  names(cf) <- coefficient_names
  dimnames(vc) <- list(coefficient_names, coefficient_names)
  list(coef = cf, vcov = vc)
}

make_capacity_curve <- function(
    model,
    capacity_label,
    capacity_hours,
    capacity_z,
    latitude,
    latitude_z,
    annual_tmean,
    annual_tmean_z,
    at_grid,
    cen,
    knots,
    boundary,
    coefficient_names) {
  reduced <- predict_reduced_coefficients(
    model = model,
    capacity_z = capacity_z,
    latitude_z = latitude_z,
    annual_tmean_z = annual_tmean_z,
    coefficient_names = coefficient_names
  )
  exposure_basis <- onebasis(
    at_grid,
    fun = "ns",
    knots = knots,
    Boundary.knots = boundary,
    intercept = FALSE
  )
  pred <- crosspred(
    exposure_basis,
    coef = reduced$coef,
    vcov = reduced$vcov,
    model.link = "logit",
    at = at_grid,
    cen = cen
  )

  out <- data.table(
    capacity_label = capacity_label,
    baseline_actual_sunlight_hours = capacity_hours,
    baseline_actual_sunlight_z = capacity_z,
    latitude_for_prediction = latitude,
    latitude_z_for_prediction = latitude_z,
    annual_tmean_for_prediction = annual_tmean,
    annual_tmean_z_for_prediction = annual_tmean_z,
    sunlight_anomaly_hours = pred$predvar,
    OR = as.numeric(pred$allRRfit),
    OR_low = as.numeric(pred$allRRlow),
    OR_high = as.numeric(pred$allRRhigh)
  )
  out[, p_value := or_p_value(OR, OR_low, OR_high)]
  out[, significance := make_sig(OR_low, OR_high)]
  out
}

make_within_curve <- function(
    model,
    at_grid,
    cen,
    knots,
    boundary,
    coefficient_names) {
  cf <- as.numeric(coef(model))
  vc <- as.matrix(vcov(model))
  names(cf) <- coefficient_names
  dimnames(vc) <- list(coefficient_names, coefficient_names)

  exposure_basis <- onebasis(
    at_grid,
    fun = "ns",
    knots = knots,
    Boundary.knots = boundary,
    intercept = FALSE
  )
  pred <- crosspred(
    exposure_basis,
    coef = cf,
    vcov = vc,
    model.link = "logit",
    at = at_grid,
    cen = cen
  )

  out <- data.table(
    sunlight_anomaly_hours = pred$predvar,
    reference_anomaly_hours = cen,
    OR = as.numeric(pred$allRRfit),
    OR_low = as.numeric(pred$allRRlow),
    OR_high = as.numeric(pred$allRRhigh)
  )
  out[, p_value := or_p_value(OR, OR_low, OR_high)]
  out[, significance := make_sig(OR_low, OR_high)]
  out
}

make_lag_response_from_var_meta <- function(model, lag_values, lag_knots) {
  cf <- as.numeric(coef(model))
  vc <- as.matrix(vcov(model))

  lag_basis <- onebasis(
    lag_values,
    fun = "ns",
    knots = lag_knots,
    intercept = TRUE
  )
  if (ncol(lag_basis) != length(cf)) {
    stop(
      "Lag-basis dimension mismatch: basis columns = ",
      ncol(lag_basis),
      ", coefficients = ",
      length(cf)
    )
  }

  log_or <- as.numeric(lag_basis %*% cf)
  var_log_or <- rowSums((lag_basis %*% vc) * lag_basis)
  se_log_or <- sqrt(pmax(var_log_or, 0))
  out <- data.table(
    lag_day = lag_values,
    log_OR = log_or,
    se_log_OR = se_log_or,
    OR = exp(log_or),
    OR_low = exp(log_or - 1.96 * se_log_or),
    OR_high = exp(log_or + 1.96 * se_log_or)
  )
  out[, p_value := 2 * pnorm(abs(log_OR / se_log_OR), lower.tail = FALSE)]
  out[, significance := make_sig(OR_low, OR_high)]
  out[, p_value_BH := p.adjust(p_value, method = "BH")]
  out[, significant_BH := is.finite(p_value_BH) & p_value_BH < 0.05 &
    (OR_low > 1 | OR_high < 1)]
  out
}

make_incremental_cumulative_lag_response_from_var_meta <- function(model, lag_values, lag_knots) {
  cf <- as.numeric(coef(model))
  vc <- as.matrix(vcov(model))

  lag_basis <- onebasis(
    lag_values,
    fun = "ns",
    knots = lag_knots,
    intercept = TRUE
  )
  if (ncol(lag_basis) != length(cf)) {
    stop(
      "Lag-basis dimension mismatch: basis columns = ",
      ncol(lag_basis),
      ", coefficients = ",
      length(cf)
    )
  }

  cumulative_basis <- apply(lag_basis, 2L, cumsum)
  log_or <- as.numeric(cumulative_basis %*% cf)
  var_log_or <- rowSums((cumulative_basis %*% vc) * cumulative_basis)
  se_log_or <- sqrt(pmax(var_log_or, 0))
  out <- data.table(
    lag_day = lag_values,
    log_OR = log_or,
    se_log_OR = se_log_or,
    OR = exp(log_or),
    OR_low = exp(log_or - 1.96 * se_log_or),
    OR_high = exp(log_or + 1.96 * se_log_or)
  )
  out[, p_value := 2 * pnorm(abs(log_OR / se_log_OR), lower.tail = FALSE)]
  out[, significance := make_sig(OR_low, OR_high)]
  out[, p_value_BH := p.adjust(p_value, method = "BH")]
  out[, significant_BH := is.finite(p_value_BH) & p_value_BH < 0.05 &
    (OR_low > 1 | OR_high < 1)]
  out
}

scalar_model_indices <- function(model) {
  nm <- names(coef(model))
  int_idx <- grep("Intercept", nm, fixed = TRUE)
  capacity_idx <- grep("baseline_actual_sunlight_z", nm, fixed = TRUE)
  latitude_idx <- grep("latitude_z", nm, fixed = TRUE)
  annual_tmean_idx <- grep("annual_tmean_z", nm, fixed = TRUE)
  if (length(int_idx) != 1L || length(capacity_idx) != 1L) {
    stop("Could not identify scalar meta-regression intercept and capacity modifier.")
  }
  list(
    intercept = int_idx,
    capacity = capacity_idx,
    latitude = if (length(latitude_idx) == 1L) latitude_idx else integer(),
    annual_tmean = if (length(annual_tmean_idx) == 1L) annual_tmean_idx else integer()
  )
}

predict_scalar_log_or <- function(model, capacity_z, latitude_z = 0, annual_tmean_z = 0) {
  b <- coef(model)
  v <- vcov(model)
  idx <- scalar_model_indices(model)
  x <- rep(0, length(b))
  x[idx$intercept] <- 1
  x[idx$capacity] <- capacity_z
  if (length(idx$latitude) == 1L) x[idx$latitude] <- latitude_z
  if (length(idx$annual_tmean) == 1L) x[idx$annual_tmean] <- annual_tmean_z
  eta <- sum(x * b)
  se <- sqrt(as.numeric(t(x) %*% v %*% x))
  data.table(
    logOR = eta,
    se_logOR = se,
    OR = exp(eta),
    OR_low = exp(eta - 1.96 * se),
    OR_high = exp(eta + 1.96 * se),
    p_value = 2 * pnorm(abs(eta / se), lower.tail = FALSE)
  )
}

scalar_capacity_contrast <- function(model, z_high, z_low) {
  b <- coef(model)
  v <- vcov(model)
  idx <- scalar_model_indices(model)
  x <- rep(0, length(b))
  x[idx$capacity] <- z_high - z_low
  eta <- sum(x * b)
  se <- sqrt(as.numeric(t(x) %*% v %*% x))
  data.table(
    log_ratio_of_ORs = eta,
    se_log_ratio_of_ORs = se,
    ratio_of_ORs = exp(eta),
    ratio_of_ORs_low = exp(eta - 1.96 * se),
    ratio_of_ORs_high = exp(eta + 1.96 * se),
    p_value = 2 * pnorm(abs(eta / se), lower.tail = FALSE)
  )
}

# -----------------------------------------------------------------------------
# 3. Read cleaned city-day data
# -----------------------------------------------------------------------------

log_msg("Output root directory: ", OUT_ROOT)
log_msg("Within output directory: ", WITHIN_DIR)
log_msg("Cross output directory: ", CROSS_DIR)
log_msg("City-day outcome path: ", CITY_DAY_PATH)
log_msg("Weather controls path: ", WEATHER_PATH)
log_msg("Sunlight source path: ", SUNLIGHT_PATH, " [daily_sun_hours]")
log_msg("City climate covariate path: ", CITY_CLIMATE_PATH)
log_msg("Primary sunlight anomaly lag: 0-", SUN_MAX_LAG)
log_msg("Temperature control: same-day nonlinear natural spline; no temperature lag")
log_msg("Pandemic phase control: four national phases; phase 1 reference")
log_msg("City-level modifier: multi-year mean observed sunlight duration")
log_msg("City-level confounder adjustment: centroid latitude and ERA5-Land annual mean temperature")

if (!file.exists(CITY_DAY_PATH)) stop("Missing CITY_DAY_PATH: ", CITY_DAY_PATH)
if (!file.exists(WEATHER_PATH)) stop("Missing WEATHER_PATH: ", WEATHER_PATH)
if (!file.exists(SUNLIGHT_PATH)) stop("Missing SUNLIGHT_PATH: ", SUNLIGHT_PATH)
if (!file.exists(CITY_CLIMATE_PATH)) stop("Missing CITY_CLIMATE_PATH: ", CITY_CLIMATE_PATH)

city_day <- fread(
  CITY_DAY_PATH,
  colClasses = list(character = "city_code"),
  showProgress = FALSE
)
weather <- fread(
  WEATHER_PATH,
  colClasses = list(character = "city_code"),
  showProgress = FALSE
)
sunlight_raw <- fread(
  SUNLIGHT_PATH,
  colClasses = list(character = "CITYCODE"),
  showProgress = FALSE
)

required_outcome <- c("city_code", "date", "high_events", "n_observed")
required_weather <- c(
  "city_code", "date", "temperature",
  "precipitation", "wind_speed", "humidity"
)
if (length(setdiff(required_outcome, names(city_day))) > 0L) {
  stop("Outcome file is missing: ", paste(setdiff(required_outcome, names(city_day)), collapse = ", "))
}
if (length(setdiff(required_weather, names(weather))) > 0L) {
  stop("Weather file is missing: ", paste(setdiff(required_weather, names(weather)), collapse = ", "))
}
required_sunlight <- c("CITYCODE", "date", "daily_sun_hours")
if (length(setdiff(required_sunlight, names(sunlight_raw))) > 0L) {
  stop("Sunlight file is missing: ", paste(setdiff(required_sunlight, names(sunlight_raw)), collapse = ", "))
}

city_day[, city_code := clean_city_code(city_code)]
city_day[, date := as.IDate(date)]
city_day <- city_day[
  !is.na(city_code) & !is.na(date) &
    date >= STUDY_START & date <= STUDY_END
]
city_day <- city_day[, .(
  high_events = sum(as.numeric(high_events), na.rm = TRUE),
  n_observed = sum(as.numeric(n_observed), na.rm = TRUE)
), by = .(city_code, date)]
city_day[, high_stress_rate := high_events / n_observed]

weather[, city_code := clean_city_code(city_code)]
weather[, date := as.IDate(date)]
weather <- weather[
  !is.na(city_code) & !is.na(date) &
    date >= STUDY_START & date <= STUDY_END
]
weather <- weather[, .(
  city_name = if ("city_name" %in% names(weather)) {
    safe_first_character(city_name)
  } else {
    NA_character_
  },
  temperature = safe_mean(as.numeric(temperature)),
  precipitation = safe_mean(as.numeric(precipitation)),
  wind_speed = safe_mean(as.numeric(wind_speed)),
  humidity = safe_mean(as.numeric(humidity))
), by = .(city_code, date)]

sunlight_raw[, city_code := clean_city_code(CITYCODE)]
sunlight_raw[, date := as.IDate(date)]
sunlight <- sunlight_raw[
  !is.na(city_code) & !is.na(date) &
    date >= STUDY_START & date <= STUDY_END,
  .(sunlight_duration = safe_mean(as.numeric(daily_sun_hours))),
  by = .(city_code, date)
]
weather <- merge(
  weather,
  sunlight,
  by = c("city_code", "date"),
  all = FALSE,
  sort = FALSE
)
setorder(weather, city_code, date)

log_msg("Outcome city-day rows: ", nrow(city_day), "; cities: ", uniqueN(city_day$city_code))
log_msg("Weather/sunlight city-day rows: ", nrow(weather), "; cities: ", uniqueN(weather$city_code))

# -----------------------------------------------------------------------------
# 4. Construct the city-level baseline capacity metric
# -----------------------------------------------------------------------------

# The primary capacity metric is the equally day-weighted arithmetic mean of
# realized daily sunlight. n_observed is deliberately not used as a weight.
weather[, year := as.integer(format(date, "%Y"))]
n_study_days <- as.integer(STUDY_END - STUDY_START) + 1L

city_annual_sunlight <- weather[
  is.finite(sunlight_duration),
  .(
    valid_sunlight_days = .N,
    annual_mean_actual_sunlight = mean(sunlight_duration)
  ),
  by = .(city_code, city_name, year)
]
save_csv(city_annual_sunlight, "city_annual_mean_actual_sunlight.csv")

city_capacity <- weather[
  is.finite(sunlight_duration),
  .(
    city_name = safe_first_character(city_name),
    baseline_actual_sunlight = mean(sunlight_duration),
    baseline_actual_sunlight_sd = sd(sunlight_duration),
    valid_sunlight_days = .N,
    valid_sunlight_years = uniqueN(year),
    first_sunlight_date = min(date),
    last_sunlight_date = max(date)
  ),
  by = city_code
]
city_capacity[, expected_study_days := n_study_days]
city_capacity[, baseline_coverage := valid_sunlight_days / expected_study_days]
city_capacity[, baseline_quality_eligible :=
  valid_sunlight_years >= MIN_BASELINE_YEARS &
    baseline_coverage >= MIN_BASELINE_COVERAGE]

city_climate <- read_city_climate_covariates(CITY_CLIMATE_PATH)
city_capacity <- merge(
  city_capacity,
  city_climate,
  by = "city_code",
  all.x = TRUE,
  sort = FALSE
)
log_msg(
  "City climate covariate rows read: ",
  nrow(city_climate),
  "; matched latitude cities: ",
  sum(is.finite(city_capacity$latitude)),
  "; matched annual Tmean cities: ",
  sum(is.finite(city_capacity$annual_tmean))
)

all_capacity_mean <- mean(city_capacity$baseline_actual_sunlight, na.rm = TRUE)
all_capacity_sd <- sd(city_capacity$baseline_actual_sunlight, na.rm = TRUE)
city_capacity[, baseline_actual_sunlight_z_all_cities := standardize_with(
  baseline_actual_sunlight,
  all_capacity_mean,
  all_capacity_sd
)]
save_csv(city_capacity, "city_multi_year_mean_actual_sunlight_capacity.csv")

capacity_distribution_all <- data.table(
  statistic = c("n_cities", "mean", "sd", "min", "p10", "p25", "p50", "p75", "p90", "max"),
  value = c(
    nrow(city_capacity),
    mean(city_capacity$baseline_actual_sunlight, na.rm = TRUE),
    sd(city_capacity$baseline_actual_sunlight, na.rm = TRUE),
    min(city_capacity$baseline_actual_sunlight, na.rm = TRUE),
    quantile(city_capacity$baseline_actual_sunlight, 0.10, na.rm = TRUE),
    quantile(city_capacity$baseline_actual_sunlight, 0.25, na.rm = TRUE),
    quantile(city_capacity$baseline_actual_sunlight, 0.50, na.rm = TRUE),
    quantile(city_capacity$baseline_actual_sunlight, 0.75, na.rm = TRUE),
    quantile(city_capacity$baseline_actual_sunlight, 0.90, na.rm = TRUE),
    max(city_capacity$baseline_actual_sunlight, na.rm = TRUE)
  )
)
save_csv(capacity_distribution_all, "baseline_actual_sunlight_distribution_all_cities.csv")

latitude_distribution_all <- data.table(
  statistic = c("n_cities", "mean", "sd", "min", "p10", "p25", "p50", "p75", "p90", "max"),
  value = c(
    sum(is.finite(city_capacity$latitude)),
    mean(city_capacity$latitude, na.rm = TRUE),
    sd(city_capacity$latitude, na.rm = TRUE),
    min(city_capacity$latitude, na.rm = TRUE),
    quantile(city_capacity$latitude, 0.10, na.rm = TRUE),
    quantile(city_capacity$latitude, 0.25, na.rm = TRUE),
    quantile(city_capacity$latitude, 0.50, na.rm = TRUE),
    quantile(city_capacity$latitude, 0.75, na.rm = TRUE),
    quantile(city_capacity$latitude, 0.90, na.rm = TRUE),
    max(city_capacity$latitude, na.rm = TRUE)
  )
)
save_csv(latitude_distribution_all, "city_centroid_latitude_distribution_all_cities.csv")

annual_tmean_distribution_all <- data.table(
  statistic = c("n_cities", "mean", "sd", "min", "p10", "p25", "p50", "p75", "p90", "max"),
  value = c(
    sum(is.finite(city_capacity$annual_tmean)),
    mean(city_capacity$annual_tmean, na.rm = TRUE),
    sd(city_capacity$annual_tmean, na.rm = TRUE),
    min(city_capacity$annual_tmean, na.rm = TRUE),
    quantile(city_capacity$annual_tmean, 0.10, na.rm = TRUE),
    quantile(city_capacity$annual_tmean, 0.25, na.rm = TRUE),
    quantile(city_capacity$annual_tmean, 0.50, na.rm = TRUE),
    quantile(city_capacity$annual_tmean, 0.75, na.rm = TRUE),
    quantile(city_capacity$annual_tmean, 0.90, na.rm = TRUE),
    max(city_capacity$annual_tmean, na.rm = TRUE)
  )
)
save_csv(annual_tmean_distribution_all, "city_annual_tmean_distribution_all_cities.csv")

# -----------------------------------------------------------------------------
# 5. Construct leave-one-year-out city-season sunlight anomalies
# -----------------------------------------------------------------------------

log_msg(
  "Constructing sunlight anomaly: leave-one-year-out, +/-",
  CLIM_HALF_WINDOW,
  " calendar days"
)

# Map month-day to a leap-year calendar so all years use a common 1-366 index.
weather[, clim_doy := as.integer(format(
  as.IDate(paste0("2000-", format(date, "%m-%d"))),
  "%j"
))]

clim_daily <- weather[, .(
  sunlight_duration = safe_mean(sunlight_duration)
), by = .(city_code, year, clim_doy)]

clim_grid <- CJ(
  city_code = sort(unique(weather$city_code)),
  year = sort(unique(weather$year)),
  clim_doy = 1:366,
  unique = TRUE
)
clim_grid <- merge(
  clim_grid,
  clim_daily,
  by = c("city_code", "year", "clim_doy"),
  all.x = TRUE,
  sort = FALSE
)
setorder(clim_grid, city_code, year, clim_doy)

window_stats <- clim_grid[, {
  s <- make_circular_window_stats(sunlight_duration, CLIM_HALF_WINDOW)
  .(
    clim_doy = 1:366,
    year_window_sum = s$sum,
    year_window_sumsq = s$sumsq,
    year_window_n = s$n
  )
}, by = .(city_code, year)]

all_year_stats <- window_stats[, .(
  all_year_sum = sum(year_window_sum, na.rm = TRUE),
  all_year_sumsq = sum(year_window_sumsq, na.rm = TRUE),
  all_year_n = sum(year_window_n, na.rm = TRUE)
), by = .(city_code, clim_doy)]

window_stats <- merge(
  window_stats,
  all_year_stats,
  by = c("city_code", "clim_doy"),
  all.x = TRUE,
  sort = FALSE
)
window_stats[, `:=`(
  donor_sum = all_year_sum - year_window_sum,
  donor_sumsq = all_year_sumsq - year_window_sumsq,
  donor_n = all_year_n - year_window_n
)]

weather_anom <- merge(
  weather,
  window_stats[, .(
    city_code, year, clim_doy, donor_sum, donor_sumsq, donor_n
  )],
  by = c("city_code", "year", "clim_doy"),
  all.x = TRUE,
  sort = FALSE
)
weather_anom[, expected_sunlight := fifelse(
  donor_n >= MIN_CLIMATOLOGY_N,
  donor_sum / donor_n,
  NA_real_
)]
weather_anom[, donor_variance := fifelse(
  donor_n > 1,
  pmax((donor_sumsq - donor_sum^2 / donor_n) / (donor_n - 1), 0),
  NA_real_
)]
weather_anom[, climatology_sd := sqrt(donor_variance)]
weather_anom[, sunlight_anomaly := sunlight_duration - expected_sunlight]
weather_anom[, sunlight_anomaly_z := fifelse(
  is.finite(climatology_sd) & climatology_sd > 0,
  sunlight_anomaly / climatology_sd,
  NA_real_
)]
setorder(weather_anom, city_code, date)

save_csv(
  weather_anom[, .(
    city_code, city_name, date, year, sunlight_duration,
    expected_sunlight, sunlight_anomaly, sunlight_anomaly_z,
    donor_n, temperature, humidity, precipitation, wind_speed
  )],
  "weather_with_leave_one_year_out_sunlight_anomaly.csv"
)

# -----------------------------------------------------------------------------
# 6. Build a complete city-date panel and global basis parameters
# -----------------------------------------------------------------------------

analysis_cities <- sort(intersect(
  unique(city_day$city_code),
  unique(weather_anom$city_code)
))
all_dates <- as.IDate(seq(as.Date(STUDY_START), as.Date(STUDY_END), by = "day"))
panel <- CJ(city_code = analysis_cities, date = all_dates, unique = TRUE)
panel <- merge(panel, city_day, by = c("city_code", "date"), all.x = TRUE, sort = FALSE)
panel <- merge(
  panel,
  weather_anom[, .(
    city_code, date, city_name, sunlight_duration, expected_sunlight,
    sunlight_anomaly, sunlight_anomaly_z, donor_n,
    temperature, humidity, precipitation, wind_speed
  )],
  by = c("city_code", "date"),
  all.x = TRUE,
  sort = FALSE
)
setorder(panel, city_code, date)
panel[, failures := n_observed - high_events]
panel[, time_num := as.numeric(date - STUDY_START) + 1]
panel[, pandemic_phase := make_pandemic_phase(date)]
panel[, dow := factor(
  as.POSIXlt(as.Date(date))$wday,
  levels = 0:6,
  labels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
)]

# Full city-day analytic panel is kept in memory; it is not exported.

anomaly_global <- panel[is.finite(sunlight_anomaly), sunlight_anomaly]
temperature_global <- panel[is.finite(temperature), temperature]
if (length(anomaly_global) == 0L) stop("No finite sunlight anomaly values.")
if (length(temperature_global) == 0L) stop("No finite temperature values.")

anomaly_boundary <- range(anomaly_global, na.rm = TRUE)
temperature_boundary <- range(temperature_global, na.rm = TRUE)
anomaly_knots <- as.numeric(quantile(
  anomaly_global,
  SUN_VAR_KNOT_PROBS,
  na.rm = TRUE,
  type = 7
))
temperature_knots <- as.numeric(quantile(
  temperature_global,
  TEMP_VAR_KNOT_PROBS,
  na.rm = TRUE,
  type = 7
))
sun_lag_knots <- SUN_LAG_KNOTS
if (any(sun_lag_knots <= 0 | sun_lag_knots >= SUN_MAX_LAG)) {
  stop("SUN_LAG_KNOTS must be strictly between 0 and SUN_MAX_LAG.")
}

if (anyDuplicated(anomaly_knots)) stop("Sunlight anomaly knots are duplicated.")
if (anyDuplicated(temperature_knots)) stop("Temperature knots are duplicated.")
if (!(anomaly_boundary[1] < 0 && anomaly_boundary[2] > 0)) {
  stop("Sunlight anomaly boundary does not contain the reference value 0.")
}

anomaly_probs <- c(
  p01 = 0.01, p05 = 0.05, p10 = 0.10, p25 = 0.25,
  p50 = 0.50, p75 = 0.75, p90 = 0.90, p95 = 0.95, p99 = 0.99
)
anomaly_quantiles <- setNames(
  as.numeric(quantile(anomaly_global, anomaly_probs, na.rm = TRUE, type = 7)),
  names(anomaly_probs)
)

global_parameters <- rbindlist(list(
  data.table(
    parameter = paste0("sunlight_anomaly_exposure_knot_", seq_along(anomaly_knots)),
    value = anomaly_knots
  ),
  data.table(
    parameter = paste0("sunlight_anomaly_lag_knot_", seq_along(sun_lag_knots)),
    value = sun_lag_knots
  ),
  data.table(
    parameter = paste0("temperature_exposure_knot_", seq_along(temperature_knots)),
    value = temperature_knots
  ),
  data.table(
    parameter = c(
      "sunlight_anomaly_boundary_lower", "sunlight_anomaly_boundary_upper",
      "temperature_boundary_lower", "temperature_boundary_upper",
      "sunlight_anomaly_reference", "sunlight_max_lag",
      "temperature_max_lag", "time_df"
    ),
    value = c(
      anomaly_boundary, temperature_boundary, 0,
      SUN_MAX_LAG, TEMP_MAX_LAG, TIME_DF
    )
  )
), fill = TRUE)
save_csv_both(global_parameters, "global_dlnm_parameters.csv")
save_csv_both(
  data.table(
    percentile = names(anomaly_quantiles),
    probability = as.numeric(anomaly_probs),
    sunlight_anomaly_hours = as.numeric(anomaly_quantiles)
  ),
  "global_sunlight_anomaly_percentiles.csv"
)

log_msg("Global anomaly knots: ", paste(round(anomaly_knots, 4), collapse = ", "))
log_msg("Global anomaly lag knots: ", paste(round(sun_lag_knots, 4), collapse = ", "))
log_msg("Global temperature knots: ", paste(round(temperature_knots, 4), collapse = ", "))

# -----------------------------------------------------------------------------
# 7. Stage 1: city-specific anomaly DLNMs
# -----------------------------------------------------------------------------

run_city_model <- function(city_id) {
  dt <- copy(panel[city_code == city_id])
  setorder(dt, date)
  city_label <- safe_first_character(dt$city_name)

  base_status <- data.table(
    city_code = city_id,
    city_name = city_label,
    panel_days = nrow(dt),
    complete_days_after_lag = NA_integer_,
    total_observed_after_lag = NA_real_,
    total_high_events_after_lag = NA_real_,
    fit_status = "not_started",
    message = NA_character_
  )
  var_reduction_status <- "not_run"
  var_reduction_message <- NA_character_
  var_coef <- NULL
  var_vcov <- NULL

  tryCatch({
    cb_sun <- crossbasis(
      dt$sunlight_anomaly,
      lag = SUN_MAX_LAG,
      argvar = list(
        fun = "ns",
        knots = anomaly_knots,
        Boundary.knots = anomaly_boundary,
        intercept = FALSE
      ),
      arglag = list(
        fun = "ns",
        knots = sun_lag_knots,
        intercept = TRUE
      )
    )
    complete_after_lag <- complete.cases(cbind(
      dt[, .(
        high_stress_rate, n_observed, high_events, failures,
        temperature, humidity, precipitation, wind_speed,
        time_num, dow, pandemic_phase
      )],
      cb_sun
    )) &
      dt$n_observed > 0 & dt$high_events >= 0 & dt$failures >= 0

    n_complete <- sum(complete_after_lag)
    total_observed <- sum(dt$n_observed[complete_after_lag], na.rm = TRUE)
    total_events <- sum(dt$high_events[complete_after_lag], na.rm = TRUE)
    base_status[, `:=`(
      complete_days_after_lag = n_complete,
      total_observed_after_lag = total_observed,
      total_high_events_after_lag = total_events
    )]

    if (n_complete < MIN_MODEL_DAYS) {
      base_status[, `:=`(
        fit_status = "skipped",
        message = paste0("complete days after lag < ", MIN_MODEL_DAYS)
      )]
      return(list(status = base_status))
    }
    if (total_events < MIN_TOTAL_EVENTS) {
      base_status[, `:=`(
        fit_status = "skipped",
        message = paste0("total high events < ", MIN_TOTAL_EVENTS)
      )]
      return(list(status = base_status))
    }
    if (total_observed < MIN_TOTAL_OBSERVED) {
      base_status[, `:=`(
        fit_status = "skipped",
        message = paste0("total observed < ", MIN_TOTAL_OBSERVED)
      )]
      return(list(status = base_status))
    }

    fit <- glm(
      high_stress_rate ~
        cb_sun +
        ns(
          temperature,
          knots = temperature_knots,
          Boundary.knots = temperature_boundary
        ) +
        humidity +
        precipitation +
        wind_speed +
        ns(time_num, df = TIME_DF) +
        dow +
        pandemic_phase,
      family = quasibinomial(),
      weights = n_observed,
      data = dt,
      na.action = na.exclude,
      control = glm.control(maxit = 100)
    )

    red_sun <- crossreduce(
      cb_sun,
      fit,
      type = "overall",
      cen = TARGET_ANOMALY_HOURS
    )
    cf <- coef(red_sun)
    vc <- vcov(red_sun)
    if (any(!is.finite(cf)) || any(!is.finite(vc))) {
      stop("Non-finite reduced coefficient or covariance.")
    }

    converged <- isTRUE(fit$converged)
    if (converged) {
      var_reduction <- tryCatch(
        crossreduce(
          cb_sun,
          fit,
          type = "var",
          value = 0,
          cen = TARGET_ANOMALY_HOURS
        ),
        error = function(e) e
      )
      if (inherits(var_reduction, "error")) {
        var_reduction_status <- "failed"
        var_reduction_message <- conditionMessage(var_reduction)
      } else {
        var_coef <- coef(var_reduction)
        var_vcov <- vcov(var_reduction)
        if (any(!is.finite(var_coef)) || any(!is.finite(var_vcov))) {
          var_reduction_status <- "failed"
          var_reduction_message <- "Non-finite type='var' coefficient or covariance."
          var_coef <- NULL
          var_vcov <- NULL
        } else {
          var_reduction_status <- "ok"
        }
      }
    }
    base_status[, `:=`(
      fit_status = ifelse(converged, "ok", "not_converged"),
      model_days_used = nobs(fit),
      dispersion = summary(fit)$dispersion,
      glm_iterations = fit$iter,
      glm_converged = converged,
      var_reduction_status = var_reduction_status,
      var_reduction_message = var_reduction_message,
      message = ifelse(
        converged,
        NA_character_,
        "glm did not converge; excluded from stage 2"
      )
    )]

    list(
      status = base_status,
      coef = cf,
      vcov = vc,
      var_coef = var_coef,
      var_vcov = var_vcov
    )
  }, error = function(e) {
    base_status[, `:=`(
      fit_status = "failed",
      var_reduction_status = var_reduction_status,
      var_reduction_message = var_reduction_message,
      message = conditionMessage(e)
    )]
    list(status = base_status)
  })
}

model_cities <- sort(unique(panel$city_code))
log_msg("Fitting city-specific anomaly DLNMs for ", length(model_cities), " cities...")
stage1_results <- vector("list", length(model_cities))
names(stage1_results) <- model_cities

for (i in seq_along(model_cities)) {
  if (i == 1L || i %% 20L == 0L || i == length(model_cities)) {
    log_msg("Stage-1 progress: ", i, "/", length(model_cities))
  }
  stage1_results[[i]] <- run_city_model(model_cities[i])
}

stage1_status <- rbindlist(lapply(stage1_results, `[[`, "status"), fill = TRUE)
save_csv_both(stage1_status, "stage1_city_model_status.csv")

ok_cities_initial <- stage1_status[fit_status == "ok", city_code]
if (length(ok_cities_initial) < 5L) {
  stop("Too few successful city-specific models: ", length(ok_cities_initial))
}

coef_list_initial <- lapply(stage1_results[ok_cities_initial], `[[`, "coef")
vcov_list_initial <- lapply(stage1_results[ok_cities_initial], `[[`, "vcov")
coef_matrix_initial <- do.call(rbind, coef_list_initial)
rownames(coef_matrix_initial) <- ok_cities_initial

log_msg("Stage-1 successful models: ", length(ok_cities_initial))

ok_cities_within <- ok_cities_initial
coef_matrix_within <- coef_matrix_initial
vcov_list_within <- vcov_list_initial

var_status_within <- rbindlist(lapply(ok_cities_within, function(city_id) {
  city_result <- stage1_results[[city_id]]
  has_var <- !is.null(city_result$var_coef) && !is.null(city_result$var_vcov)
  data.table(
    city_code = city_id,
    var_reduction_status = if (has_var) "ok" else "failed_or_missing",
    var_reduction_message = city_result$status$var_reduction_message[1L]
  )
}), fill = TRUE)
save_csv_within(
  var_status_within,
  "stage1_within_A_0_vs_A_minus_2_768_var_reduction_status.csv"
)

var_ok_cities_within <- var_status_within[var_reduction_status == "ok", city_code]
if (length(var_ok_cities_within) >= 5L) {
  var_coef_list_within <- lapply(stage1_results[var_ok_cities_within], `[[`, "var_coef")
  var_vcov_list_within <- lapply(stage1_results[var_ok_cities_within], `[[`, "var_vcov")
  var_coef_matrix_within <- do.call(rbind, var_coef_list_within)
  rownames(var_coef_matrix_within) <- var_ok_cities_within
} else {
  var_coef_matrix_within <- NULL
  var_vcov_list_within <- NULL
}

# -----------------------------------------------------------------------------
# 8. Stage 2: meta-regression on multi-year mean actual sunlight
# -----------------------------------------------------------------------------

meta_data <- merge(
  data.table(city_code = ok_cities_initial),
  city_capacity,
  by = "city_code",
  all.x = TRUE,
  sort = FALSE
)
meta_data <- meta_data[match(ok_cities_initial, city_code)]
stopifnot(identical(meta_data$city_code, ok_cities_initial))

capacity_complete <- is.finite(meta_data$baseline_actual_sunlight) &
  meta_data$baseline_quality_eligible %in% TRUE &
  is.finite(meta_data$latitude) &
  is.finite(meta_data$annual_tmean)
if (sum(capacity_complete) < 5L) {
  stop("Too few cities with reliable baseline sunlight, centroid latitude, and annual mean temperature data.")
}

meta_data <- meta_data[capacity_complete]
coef_matrix <- coef_matrix_initial[capacity_complete, , drop = FALSE]
vcov_list <- vcov_list_initial[capacity_complete]
ok_cities <- ok_cities_initial[capacity_complete]
stopifnot(identical(meta_data$city_code, ok_cities))

# Standardization is based only on the cities entering stage 2.
capacity_center <- mean(meta_data$baseline_actual_sunlight)
capacity_scale <- sd(meta_data$baseline_actual_sunlight)
if (!is.finite(capacity_scale) || capacity_scale <= 0) {
  stop("The stage-2 baseline sunlight metric has zero or invalid variance.")
}
meta_data[, baseline_actual_sunlight_z := standardize_with(
  baseline_actual_sunlight,
  capacity_center,
  capacity_scale
)]

latitude_center <- mean(meta_data$latitude)
latitude_scale <- sd(meta_data$latitude)
if (!is.finite(latitude_scale) || latitude_scale <= 0) {
  stop("The stage-2 latitude metric has zero or invalid variance.")
}
meta_data[, latitude_z := standardize_with(
  latitude,
  latitude_center,
  latitude_scale
)]

annual_tmean_center <- mean(meta_data$annual_tmean)
annual_tmean_scale <- sd(meta_data$annual_tmean)
if (!is.finite(annual_tmean_scale) || annual_tmean_scale <= 0) {
  stop("The stage-2 annual mean temperature metric has zero or invalid variance.")
}
meta_data[, annual_tmean_z := standardize_with(
  annual_tmean,
  annual_tmean_center,
  annual_tmean_scale
)]

capacity_values <- as.numeric(quantile(
  meta_data$baseline_actual_sunlight,
  CAPACITY_PROBS,
  na.rm = TRUE,
  type = 7
))
capacity_levels <- data.table(
  capacity_label = names(CAPACITY_PROBS),
  probability = as.numeric(CAPACITY_PROBS),
  baseline_actual_sunlight_hours = capacity_values
)
capacity_levels[, baseline_actual_sunlight_z := standardize_with(
  baseline_actual_sunlight_hours,
  capacity_center,
  capacity_scale
)]
capacity_levels[, `:=`(
  latitude_for_prediction = latitude_center,
  latitude_z_for_prediction = 0,
  annual_tmean_for_prediction = annual_tmean_center,
  annual_tmean_z_for_prediction = 0
)]

save_csv(meta_data, "stage2_city_baseline_actual_sunlight_data.csv")
save_csv(capacity_levels, "stage2_capacity_P10_P50_P90.csv")
save_csv(data.table(
  item = c(
    "capacity_definition", "day_weighting", "wearable_sample_weighting",
    "stage2_center_hours", "stage2_scale_sd_hours",
    "latitude_adjustment", "latitude_center_degrees", "latitude_scale_sd_degrees",
    "annual_tmean_adjustment", "annual_tmean_center_C", "annual_tmean_scale_sd_C",
    "minimum_baseline_years", "minimum_baseline_coverage",
    "city_level_confounders"
  ),
  value = c(
    "multi-year mean observed daily sunlight duration",
    "each valid city-day receives equal weight",
    "none",
    capacity_center,
    capacity_scale,
    "city centroid latitude included in cross-stage meta-regression; predictions are shown at mean latitude",
    latitude_center,
    latitude_scale,
    "ERA5-Land annual mean temperature for 2021-2024 included in cross-stage meta-regression; predictions are shown at mean annual temperature",
    annual_tmean_center,
    annual_tmean_scale,
    MIN_BASELINE_YEARS,
    MIN_BASELINE_COVERAGE,
    "city centroid latitude and ERA5-Land annual mean temperature"
  )
), "baseline_capacity_definition_and_scale.csv")

shared_stage2_data <- list(
  coef_matrix_within = coef_matrix_within,
  vcov_list_within = vcov_list_within,
  ok_cities_within = ok_cities_within,
  coef_matrix_cross = coef_matrix,
  vcov_list_cross = vcov_list,
  ok_cities_cross = ok_cities,
  var_coef_matrix_within = var_coef_matrix_within,
  var_vcov_list_within = var_vcov_list_within,
  var_ok_cities_within = var_ok_cities_within,
  meta_data = meta_data,
  anomaly_knots = anomaly_knots,
  anomaly_boundary = anomaly_boundary,
  sun_lag_knots = sun_lag_knots,
  temperature_knots = temperature_knots,
  temperature_boundary = temperature_boundary,
  anomaly_cen = TARGET_ANOMALY_HOURS,
  capacity_center = capacity_center,
  capacity_scale = capacity_scale,
  latitude_center = latitude_center,
  latitude_scale = latitude_scale,
  annual_tmean_center = annual_tmean_center,
  annual_tmean_scale = annual_tmean_scale
)
saveRDS(
  shared_stage2_data,
  file.path(OUT_ROOT, "stage1_reduced_coefficients_and_stage2_data.rds")
)
saveRDS(
  shared_stage2_data,
  file.path(WITHIN_DIR, "stage1_reduced_coefficients_and_stage2_data.rds")
)
saveRDS(
  shared_stage2_data,
  file.path(CROSS_DIR, "stage1_reduced_coefficients_and_stage2_data.rds")
)

save_csv_both(
  as.data.table(coef_matrix_initial, keep.rownames = "city_code"),
  "stage1_reduced_anomaly_coefficients.csv"
)
save_csv_within(
  as.data.table(coef_matrix_within, keep.rownames = "city_code"),
  "stage1_within_reduced_anomaly_coefficients.csv"
)
save_csv_cross(
  as.data.table(coef_matrix, keep.rownames = "city_code"),
  "stage1_cross_capacity_eligible_reduced_anomaly_coefficients.csv"
)

selected_exposures <- data.table(
  exposure_label = c("P5_anomaly", "P10_anomaly", "minus_2_768h", "minus_1h", "normal_0h", "plus_1h", "plus_2h"),
  anomaly_hours = c(
    anomaly_quantiles[["p05"]],
    anomaly_quantiles[["p10"]],
    TARGET_ANOMALY_HOURS, -1, 0, 1, 2
  )
)
selected_exposures <- unique(selected_exposures[
  anomaly_hours >= anomaly_boundary[1] & anomaly_hours <= anomaly_boundary[2]
])

at_grid <- sort(unique(c(
  seq(anomaly_quantiles[["p01"]], anomaly_quantiles[["p99"]], length.out = 240L),
  selected_exposures$anomaly_hours,
  0
)))

log_msg("Fitting within intercept-only random-effects meta-analysis...")
meta0 <- mixmeta(coef_matrix_within ~ 1, S = vcov_list_within, method = "reml")

log_msg("Writing within-city pooled anomaly-response outputs...")
within_curve <- make_within_curve(
  model = meta0,
  at_grid = at_grid,
  cen = TARGET_ANOMALY_HOURS,
  knots = anomaly_knots,
  boundary = anomaly_boundary,
  coefficient_names = colnames(coef_matrix_within)
)
save_csv_within(within_curve, "within_cumulative_anomaly_OR_curve_ref_minus2_768.csv")

within_selected_or <- rbindlist(lapply(seq_len(nrow(selected_exposures)), function(j) {
  idx <- which.min(abs(within_curve$sunlight_anomaly_hours - selected_exposures$anomaly_hours[j]))
  within_curve[idx, .(
    exposure_label = selected_exposures$exposure_label[j],
    requested_anomaly_hours = selected_exposures$anomaly_hours[j],
    predicted_anomaly_hours = sunlight_anomaly_hours,
    reference_anomaly_hours,
    OR, OR_low, OR_high, p_value, significance
  )]
}))
save_csv_within(within_selected_or, "within_selected_cumulative_ORs_ref_minus2_768.csv")
save_csv_within(
  within_selected_or[exposure_label == "normal_0h"],
  "within_A_0_vs_A_minus_2_768_cumulative_OR.csv"
)

meta_var_within <- NULL
within_daily_lag_or <- data.table()
within_incremental_cumulative_lag_or <- data.table()
if (!is.null(var_coef_matrix_within) && length(var_vcov_list_within) >= 5L) {
  log_msg(
    "Fitting within type='var' lag-response meta-analysis for A=",
    0,
    " versus A=-2.768..."
  )
  meta_var_within <- mixmeta(
    var_coef_matrix_within ~ 1,
    S = var_vcov_list_within,
    method = "reml"
  )
  saveRDS(
    meta_var_within,
    file.path(WITHIN_DIR, "stage2_within_A_0_vs_A_minus_2_768_type_var_meta_model.rds")
  )
  capture.output(
    summary(meta_var_within),
    file = file.path(WITHIN_DIR, "stage2_within_A_0_vs_A_minus_2_768_type_var_meta_summary.txt")
  )

  lag_values <- 0:SUN_MAX_LAG
  within_daily_lag_or <- make_lag_response_from_var_meta(
    model = meta_var_within,
    lag_values = lag_values,
    lag_knots = sun_lag_knots
  )
  within_daily_lag_or[, `:=`(
    curve_type = "within_A_0_vs_A_minus_2_768_daily_lag_response",
    exposure_label = "A_0h",
    reference_label = "A_minus_2_768h",
    exposure_anomaly_hours = 0,
    reference_anomaly_hours = TARGET_ANOMALY_HOURS,
    meta_successful_cities = length(var_ok_cities_within)
  )]
  setcolorder(
    within_daily_lag_or,
    c(
      "curve_type", "lag_day", "exposure_label", "reference_label",
      "exposure_anomaly_hours", "reference_anomaly_hours",
      "meta_successful_cities",
      setdiff(
        names(within_daily_lag_or),
        c(
          "curve_type", "lag_day", "exposure_label", "reference_label",
          "exposure_anomaly_hours", "reference_anomaly_hours",
          "meta_successful_cities"
        )
      )
    )
  )
  save_csv_within(
    within_daily_lag_or,
    "within_A_0_vs_A_minus_2_768_daily_lag_response_OR.csv"
  )

  within_incremental_cumulative_lag_or <- make_incremental_cumulative_lag_response_from_var_meta(
    model = meta_var_within,
    lag_values = lag_values,
    lag_knots = sun_lag_knots
  )
  within_incremental_cumulative_lag_or[, `:=`(
    curve_type = "within_A_0_vs_A_minus_2_768_incremental_cumulative_lag0_k",
    exposure_label = "A_0h",
    reference_label = "A_minus_2_768h",
    exposure_anomaly_hours = 0,
    reference_anomaly_hours = TARGET_ANOMALY_HOURS,
    meta_successful_cities = length(var_ok_cities_within)
  )]
  setcolorder(
    within_incremental_cumulative_lag_or,
    c(
      "curve_type", "lag_day", "exposure_label", "reference_label",
      "exposure_anomaly_hours", "reference_anomaly_hours",
      "meta_successful_cities",
      setdiff(
        names(within_incremental_cumulative_lag_or),
        c(
          "curve_type", "lag_day", "exposure_label", "reference_label",
          "exposure_anomaly_hours", "reference_anomaly_hours",
          "meta_successful_cities"
        )
      )
    )
  )
  save_csv_within(
    within_incremental_cumulative_lag_or,
    "within_A_0_vs_A_minus_2_768_incremental_cumulative_lag_OR.csv"
  )
} else {
  log_msg("Skipping within type='var' lag-response meta-analysis: fewer than 5 valid city reductions.")
}

log_msg("Fitting latitude- and annual-Tmean-adjusted capacity multivariate meta-regression...")
meta_capacity <- mixmeta(
  coef_matrix ~ baseline_actual_sunlight_z + latitude_z + annual_tmean_z,
  S = vcov_list,
  data = meta_data,
  method = "reml"
)

saveRDS(
  list(meta0 = meta0, meta_var_within = meta_var_within, meta_capacity = meta_capacity),
  file.path(OUT_ROOT, "stage2_multivariate_meta_models.rds")
)
saveRDS(
  meta0,
  file.path(WITHIN_DIR, "stage2_within_intercept_only_meta_model.rds")
)
saveRDS(
  meta_capacity,
  file.path(CROSS_DIR, "stage2_cross_baseline_capacity_latitude_tmean_adjusted_meta_regression_model.rds")
)
capture.output(
  summary(meta0),
  file = file.path(WITHIN_DIR, "stage2_within_intercept_only_meta_summary.txt")
)
capture.output(
  summary(meta_capacity),
  file = file.path(CROSS_DIR, "stage2_cross_baseline_capacity_latitude_tmean_adjusted_meta_regression_summary.txt")
)

modifier_test <- rbindlist(list(
  overall_modifier_test(
    meta_capacity,
    "baseline_actual_sunlight_z",
    "capacity_latitude_tmean_adjusted"
  ),
  overall_modifier_test(
    meta_capacity,
    "latitude_z",
    "capacity_latitude_tmean_adjusted"
  ),
  overall_modifier_test(
    meta_capacity,
    "annual_tmean_z",
    "capacity_latitude_tmean_adjusted"
  )
), fill = TRUE)
modifier_test[, interpretation := fifelse(
  modifier == "baseline_actual_sunlight_z" & is.finite(p_value) & p_value < 0.05,
  "Evidence that the full anomaly-response curve differs by baseline sunlight capacity after latitude and annual mean temperature adjustment",
  fifelse(
    modifier == "baseline_actual_sunlight_z",
    "No statistically significant evidence that the full anomaly-response curve differs by baseline sunlight capacity after latitude and annual mean temperature adjustment",
    fifelse(
      modifier == "latitude_z" & is.finite(p_value) & p_value < 0.05,
      "Evidence that the full anomaly-response curve differs by city centroid latitude after baseline sunlight and annual mean temperature adjustment",
      fifelse(
        modifier == "latitude_z",
        "No statistically significant evidence that the full anomaly-response curve differs by city centroid latitude after baseline sunlight and annual mean temperature adjustment",
        fifelse(
          is.finite(p_value) & p_value < 0.05,
          "Evidence that the full anomaly-response curve differs by annual mean temperature after baseline sunlight and latitude adjustment",
          "No statistically significant evidence that the full anomaly-response curve differs by annual mean temperature after baseline sunlight and latitude adjustment"
        )
      )
    )
  )
)]
save_csv(modifier_test, "stage2_global_capacity_latitude_tmean_modifier_wald_tests.csv")
save_csv(
  modifier_test[modifier == "baseline_actual_sunlight_z"],
  "stage2_global_capacity_modifier_wald_test.csv"
)
save_csv(
  modifier_test[modifier == "latitude_z"],
  "stage2_global_latitude_modifier_wald_test.csv"
)
save_csv(
  modifier_test[modifier == "annual_tmean_z"],
  "stage2_global_annual_tmean_modifier_wald_test.csv"
)

stage2_coef_table <- coef_table_from_model(
  meta_capacity,
  "capacity_latitude_tmean_adjusted"
)
save_csv(stage2_coef_table, "stage2_multivariate_meta_regression_coefficients.csv")

# -----------------------------------------------------------------------------
# 9. Cross branch: capacity-specific cumulative anomaly-response curves
# -----------------------------------------------------------------------------

capacity_curves <- rbindlist(lapply(seq_len(nrow(capacity_levels)), function(i) {
  make_capacity_curve(
    model = meta_capacity,
    capacity_label = capacity_levels$capacity_label[i],
    capacity_hours = capacity_levels$baseline_actual_sunlight_hours[i],
    capacity_z = capacity_levels$baseline_actual_sunlight_z[i],
    latitude = capacity_levels$latitude_for_prediction[i],
    latitude_z = capacity_levels$latitude_z_for_prediction[i],
    annual_tmean = capacity_levels$annual_tmean_for_prediction[i],
    annual_tmean_z = capacity_levels$annual_tmean_z_for_prediction[i],
    at_grid = at_grid,
    cen = TARGET_ANOMALY_HOURS,
    knots = anomaly_knots,
    boundary = anomaly_boundary,
    coefficient_names = colnames(coef_matrix)
  )
}))
save_csv(
  capacity_curves,
  "capacity_specific_cumulative_anomaly_OR_curves_ref0.csv"
)

capacity_selected_or <- rbindlist(lapply(seq_len(nrow(capacity_levels)), function(i) {
  curve_i <- capacity_curves[capacity_label == capacity_levels$capacity_label[i]]
  rbindlist(lapply(seq_len(nrow(selected_exposures)), function(j) {
    idx <- which.min(abs(curve_i$sunlight_anomaly_hours - selected_exposures$anomaly_hours[j]))
    curve_i[idx, .(
      capacity_label,
      baseline_actual_sunlight_hours,
      baseline_actual_sunlight_z,
      latitude_for_prediction,
      latitude_z_for_prediction,
      annual_tmean_for_prediction,
      annual_tmean_z_for_prediction,
      exposure_label = selected_exposures$exposure_label[j],
      requested_anomaly_hours = selected_exposures$anomaly_hours[j],
      predicted_anomaly_hours = sunlight_anomaly_hours,
      reference_anomaly_hours = TARGET_ANOMALY_HOURS,
      OR, OR_low, OR_high, p_value, significance
    )]
  }))
}))
save_csv(
  capacity_selected_or,
  "capacity_specific_selected_cumulative_ORs_ref0.csv"
)
save_csv(
  capacity_selected_or[exposure_label == "normal_0h"],
  "cross_A_0_vs_A_minus_2_768_capacity_P10_P50_P90_OR.csv"
)

# -----------------------------------------------------------------------------
# 10. Scalar meta-regressions and ratio-of-ORs interpretation
# -----------------------------------------------------------------------------

scalar_city_effects <- rbindlist(lapply(seq_along(ok_cities), function(i) {
  rbindlist(lapply(seq_len(nrow(selected_exposures)), function(j) {
    value <- selected_exposures$anomaly_hours[j]
    if (value == TARGET_ANOMALY_HOURS) return(NULL)
    cst <- contrast_log_or(
      value = value,
      cen = TARGET_ANOMALY_HOURS,
      coef_vec = coef_matrix[i, ],
      vcov_mat = vcov_list[[i]],
      knots = anomaly_knots,
      boundary = anomaly_boundary
    )
    data.table(
      city_code = ok_cities[i],
      exposure_label = selected_exposures$exposure_label[j],
      anomaly_hours = value,
      reference_anomaly_hours = TARGET_ANOMALY_HOURS,
      logOR = cst$log_or,
      var_logOR = cst$var_log_or,
      se_logOR = sqrt(cst$var_log_or),
      OR = exp(cst$log_or),
      OR_low = exp(cst$log_or - 1.96 * sqrt(cst$var_log_or)),
      OR_high = exp(cst$log_or + 1.96 * sqrt(cst$var_log_or))
    )
  }), fill = TRUE)
}), fill = TRUE)

scalar_city_effects <- merge(
  scalar_city_effects,
  meta_data[, .(
    city_code, city_name, baseline_actual_sunlight,
    baseline_actual_sunlight_z, latitude, latitude_z,
    annual_tmean, annual_tmean_z,
    valid_sunlight_days, valid_sunlight_years, baseline_coverage
  )],
  by = "city_code",
  all.x = TRUE,
  sort = FALSE
)
save_csv(scalar_city_effects, "city_scalar_anomaly_ORs_with_baseline_capacity.csv")

scalar_models <- list()
scalar_model_coefficients <- list()
scalar_predictions <- list()
scalar_ror <- list()

for (label in unique(scalar_city_effects$exposure_label)) {
  dt <- scalar_city_effects[
    exposure_label == label &
      is.finite(logOR) & is.finite(var_logOR) & var_logOR > 0 &
      is.finite(baseline_actual_sunlight_z) &
      is.finite(latitude_z) &
      is.finite(annual_tmean_z)
  ]
  if (nrow(dt) < 5L) {
    log_msg("Skipping scalar meta-regression for ", label, ": fewer than 5 valid cities")
    next
  }

  fit <- mixmeta(
    logOR ~ baseline_actual_sunlight_z + latitude_z + annual_tmean_z,
    S = var_logOR,
    data = dt,
    method = "reml"
  )
  scalar_models[[label]] <- fit

  b <- coef(fit)
  v <- vcov(fit)
  modifier_rows <- data.table(
    modifier = c("baseline_actual_sunlight_z", "latitude_z", "annual_tmean_z"),
    original_unit = c("hours/day", "degrees latitude", "degrees C annual mean temperature"),
    scale_sd = c(capacity_scale, latitude_scale, annual_tmean_scale)
  )
  modifier_rows[, term_index := match(modifier, names(b))]
  modifier_rows <- modifier_rows[!is.na(term_index)]
  modifier_rows[, `:=`(
    exposure_label = label,
    anomaly_hours = dt$anomaly_hours[1L],
    n_cities = nrow(dt),
    log_ratio_of_ORs_per_1SD = b[term_index],
    se_per_1SD = sqrt(v[cbind(term_index, term_index)])
  )]
  modifier_rows[, `:=`(
    ratio_of_ORs_per_1SD = exp(log_ratio_of_ORs_per_1SD),
    ratio_of_ORs_per_1SD_low = exp(log_ratio_of_ORs_per_1SD - 1.96 * se_per_1SD),
    ratio_of_ORs_per_1SD_high = exp(log_ratio_of_ORs_per_1SD + 1.96 * se_per_1SD),
    p_value = 2 * pnorm(abs(log_ratio_of_ORs_per_1SD / se_per_1SD), lower.tail = FALSE),
    log_ratio_of_ORs_per_1unit = log_ratio_of_ORs_per_1SD / scale_sd,
    ratio_of_ORs_per_1unit = exp(log_ratio_of_ORs_per_1SD / scale_sd),
    city_level_confounders = "city centroid latitude and ERA5-Land annual mean temperature included in the same scalar meta-regression"
  )]
  scalar_model_coefficients[[label]] <- modifier_rows[, .(
    exposure_label, anomaly_hours, n_cities, modifier, original_unit,
    log_ratio_of_ORs_per_1SD, se_per_1SD,
    ratio_of_ORs_per_1SD, ratio_of_ORs_per_1SD_low,
    ratio_of_ORs_per_1SD_high, p_value,
    log_ratio_of_ORs_per_1unit, ratio_of_ORs_per_1unit,
    city_level_confounders
  )]

  scalar_predictions[[label]] <- rbindlist(lapply(seq_len(nrow(capacity_levels)), function(i) {
    pr <- predict_scalar_log_or(
      fit,
      capacity_levels$baseline_actual_sunlight_z[i],
      latitude_z = capacity_levels$latitude_z_for_prediction[i],
      annual_tmean_z = capacity_levels$annual_tmean_z_for_prediction[i]
    )
    cbind(
      data.table(
        exposure_label = label,
        anomaly_hours = dt$anomaly_hours[1L],
        capacity_label = capacity_levels$capacity_label[i],
        baseline_actual_sunlight_hours = capacity_levels$baseline_actual_sunlight_hours[i],
        baseline_actual_sunlight_z = capacity_levels$baseline_actual_sunlight_z[i],
        latitude_for_prediction = capacity_levels$latitude_for_prediction[i],
        latitude_z_for_prediction = capacity_levels$latitude_z_for_prediction[i],
        annual_tmean_for_prediction = capacity_levels$annual_tmean_for_prediction[i],
        annual_tmean_z_for_prediction = capacity_levels$annual_tmean_z_for_prediction[i],
        reference_anomaly_hours = TARGET_ANOMALY_HOURS
      ),
      pr
    )
  }))

  low_row <- capacity_levels[capacity_label == "P10"]
  high_row <- capacity_levels[capacity_label == "P90"]
  ror <- scalar_capacity_contrast(
    fit,
    z_high = high_row$baseline_actual_sunlight_z,
    z_low = low_row$baseline_actual_sunlight_z
  )
  scalar_ror[[label]] <- cbind(
    data.table(
      exposure_label = label,
      anomaly_hours = dt$anomaly_hours[1L],
      high_capacity_label = "P90",
      high_capacity_hours = high_row$baseline_actual_sunlight_hours,
      low_capacity_label = "P10",
      low_capacity_hours = low_row$baseline_actual_sunlight_hours,
      latitude_held_at_degrees = latitude_center,
      latitude_z_held_at = 0,
      annual_tmean_held_at_C = annual_tmean_center,
      annual_tmean_z_held_at = 0,
      comparison = "OR at P90 capacity divided by OR at P10 capacity, adjusted at mean latitude and mean annual temperature"
    ),
    ror
  )
}

saveRDS(scalar_models, file.path(OUT_DIR, "scalar_capacity_meta_regression_models.rds"))

scalar_model_coefficients_dt <- rbindlist(scalar_model_coefficients, fill = TRUE)
scalar_predictions_dt <- rbindlist(scalar_predictions, fill = TRUE)
scalar_ror_dt <- rbindlist(scalar_ror, fill = TRUE)

if (nrow(scalar_model_coefficients_dt) > 0L) {
  scalar_model_coefficients_dt[, buffering_interpretation := fifelse(
    modifier == "baseline_actual_sunlight_z" & ratio_of_ORs_per_1SD_high < 1,
    "Higher baseline actual sunlight significantly attenuates this anomaly OR",
    fifelse(
      modifier == "baseline_actual_sunlight_z" & ratio_of_ORs_per_1SD_low > 1,
      "Higher baseline actual sunlight significantly amplifies this anomaly OR",
      fifelse(
      modifier == "latitude_z" & ratio_of_ORs_per_1SD_high < 1,
        "Higher latitude significantly attenuates this anomaly OR after baseline sunlight and annual mean temperature adjustment",
        fifelse(
          modifier == "latitude_z" & ratio_of_ORs_per_1SD_low > 1,
          "Higher latitude significantly amplifies this anomaly OR after baseline sunlight and annual mean temperature adjustment",
          fifelse(
            modifier == "annual_tmean_z" & ratio_of_ORs_per_1SD_high < 1,
            "Higher annual mean temperature significantly attenuates this anomaly OR after baseline sunlight and latitude adjustment",
            fifelse(
              modifier == "annual_tmean_z" & ratio_of_ORs_per_1SD_low > 1,
              "Higher annual mean temperature significantly amplifies this anomaly OR after baseline sunlight and latitude adjustment",
              "No statistically significant modification for this contrast"
            )
          )
        )
      )
    )
  )]
}
if (nrow(scalar_ror_dt) > 0L) {
  scalar_ror_dt[, buffering_interpretation := fifelse(
    ratio_of_ORs_high < 1,
    "P90-capacity cities have a significantly smaller anomaly OR than P10-capacity cities",
    fifelse(
      ratio_of_ORs_low > 1,
      "P90-capacity cities have a significantly larger anomaly OR than P10-capacity cities",
      "No statistically significant P90-versus-P10 difference"
    )
  )]
}

save_csv(
  scalar_model_coefficients_dt,
  "scalar_capacity_latitude_tmean_adjusted_modifier_results_per_1SD_and_1unit.csv"
)
save_csv(
  scalar_model_coefficients_dt[modifier == "baseline_actual_sunlight_z"],
  "scalar_capacity_modifier_results_per_1SD_and_1hour.csv"
)
save_csv(
  scalar_model_coefficients_dt[exposure_label == "normal_0h"],
  "cross_A_0_vs_A_minus_2_768_capacity_latitude_tmean_adjusted_modifiers_per_1SD_and_1unit.csv"
)
save_csv(
  scalar_model_coefficients_dt[
    exposure_label == "normal_0h" &
      modifier == "baseline_actual_sunlight_z"
  ],
  "cross_A_0_vs_A_minus_2_768_modifier_per_1SD_and_1hour.csv"
)
save_csv(
  scalar_predictions_dt,
  "scalar_anomaly_OR_predictions_at_capacity_P10_P50_P90.csv"
)
save_csv(
  scalar_ror_dt,
  "capacity_buffering_ratio_of_ORs_P90_vs_P10.csv"
)
save_csv(
  scalar_ror_dt[exposure_label == "normal_0h"],
  "cross_A_0_vs_A_minus_2_768_capacity_P90_vs_P10_ROR.csv"
)

# -----------------------------------------------------------------------------
# 11. Diagnostics and completion
# -----------------------------------------------------------------------------
# 12. Diagnostics and completion
# -----------------------------------------------------------------------------

within_target <- within_selected_or[exposure_label == "normal_0h"]
within_diagnostics <- data.table(
  metric = c(
    "analysis_design",
    "outcome",
    "sunlight_anomaly_definition",
    "sunlight_anomaly_reference",
    "target_first_stage_contrast",
    "target_cumulative_OR",
    "target_cumulative_OR_low",
    "target_cumulative_OR_high",
    "target_cumulative_p_value",
    "sunlight_anomaly_lag",
    "temperature_lag",
    "pandemic_control",
    "stage1_time_control",
    "stage1_daily_weather_controls",
    "stage2_model",
    "within_type_var_successful_cities",
    "within_cumulative_curve_display_P1_P99",
    "within_cumulative_curve_display_target_symmetric",
    "outcome_city_day_rows",
    "weather_city_day_rows",
    "complete_panel_rows",
    "complete_panel_cities",
    "stage1_successful_cities_before_capacity_filter",
    "within_stage2_cities",
    "study_start",
    "study_end"
  ),
  value = as.character(c(
    "within-city sunlight anomaly design with intercept-only random-effects meta-analysis",
    "city-day high stress event rate; quasibinomial weighted by n_observed",
    paste0("leave-one-year-out city calendar-day climatology +/-", CLIM_HALF_WINDOW, " days"),
    "-2.768 hours/day below city-season expectation",
    "A=0 hours/day versus A=-2.768 hours/day",
    within_target$OR[1L],
    within_target$OR_low[1L],
    within_target$OR_high[1L],
    within_target$p_value[1L],
    paste0("0-", SUN_MAX_LAG),
    "none; same-day temperature only",
    "four national pandemic phases; phase 1 reference",
    paste0("natural spline of time, df=", TIME_DF, ", plus day of week"),
    "same-day nonlinear temperature spline, humidity, precipitation, wind speed",
    "multivariate random-effects meta-analysis of crossreduce(type='overall') coefficients; intercept only",
    length(var_ok_cities_within),
    paste0(anomaly_quantiles[["p01"]], " to ", anomaly_quantiles[["p99"]], " hours/day"),
    paste0(TARGET_ANOMALY_HOURS, " to +", abs(TARGET_ANOMALY_HOURS), " hours/day"),
    nrow(city_day),
    nrow(weather),
    nrow(panel),
    uniqueN(panel$city_code),
    length(ok_cities_initial),
    length(ok_cities_within),
    as.character(STUDY_START),
    as.character(STUDY_END)
  ))
)
save_csv_within(within_diagnostics, "run_diagnostics.csv")

diagnostics <- data.table(
  metric = c(
    "analysis_design",
    "outcome",
    "sunlight_anomaly_definition",
    "sunlight_anomaly_reference",
    "sunlight_anomaly_lag",
    "temperature_lag",
    "pandemic_control",
    "stage1_time_control",
    "stage1_daily_weather_controls",
    "stage2_modifier",
    "stage2_city_level_confounders",
    "stage2_prediction_latitude",
    "stage2_prediction_annual_tmean",
    "baseline_day_weighting",
    "baseline_wearable_sample_weighting",
    "outcome_city_day_rows",
    "weather_city_day_rows",
    "complete_panel_rows",
    "complete_panel_cities",
    "stage1_successful_cities_before_capacity_filter",
    "stage2_cities",
    "capacity_center_hours",
    "capacity_scale_sd_hours",
    "latitude_center_degrees",
    "latitude_scale_sd_degrees",
    "annual_tmean_center_C",
    "annual_tmean_scale_sd_C",
    "global_capacity_modifier_wald_p_value",
    "global_latitude_modifier_wald_p_value",
    "global_annual_tmean_modifier_wald_p_value",
    "study_start",
    "study_end"
  ),
  value = as.character(c(
    "cross-city baseline actual sunlight capacity effect-modification design",
    "city-day high stress event rate; quasibinomial weighted by n_observed",
    paste0("leave-one-year-out city calendar-day climatology +/-", CLIM_HALF_WINDOW, " days"),
    "0 hours/day (observed equals city-season expectation)",
    paste0("0-", SUN_MAX_LAG),
    "none; same-day temperature only",
    "four national pandemic phases; phase 1 reference",
    paste0("natural spline of time, df=", TIME_DF, ", plus day of week"),
    "same-day nonlinear temperature spline, humidity, precipitation, wind speed",
    "multi-year mean observed daily sunlight duration, standardized across stage-2 cities",
    "city centroid latitude and ERA5-Land annual mean temperature, both standardized across stage-2 cities",
    "capacity-specific curves and fitted scalar lines are predicted at mean latitude",
    "capacity-specific curves and fitted scalar lines are predicted at mean annual temperature",
    "equal weight per valid city-day",
    "none",
    nrow(city_day),
    nrow(weather),
    nrow(panel),
    uniqueN(panel$city_code),
    length(ok_cities_initial),
    length(ok_cities),
    capacity_center,
    capacity_scale,
    latitude_center,
    latitude_scale,
    annual_tmean_center,
    annual_tmean_scale,
    modifier_test[modifier == "baseline_actual_sunlight_z", p_value][1L],
    modifier_test[modifier == "latitude_z", p_value][1L],
    modifier_test[modifier == "annual_tmean_z", p_value][1L],
    as.character(STUDY_START),
    as.character(STUDY_END)
  ))
)
save_csv_cross(diagnostics, "run_diagnostics.csv")

log_msg("Stage-2 cities: ", length(ok_cities))
log_msg("Capacity mean (hours/day): ", round(capacity_center, 4))
log_msg("Capacity SD (hours/day): ", round(capacity_scale, 4))
log_msg("Latitude mean (degrees): ", round(latitude_center, 4))
log_msg("Latitude SD (degrees): ", round(latitude_scale, 4))
log_msg("Annual Tmean mean (C): ", round(annual_tmean_center, 4))
log_msg("Annual Tmean SD (C): ", round(annual_tmean_scale, 4))
log_msg(
  "Global capacity modifier Wald p-value after latitude and annual Tmean adjustment: ",
  signif(modifier_test[modifier == "baseline_actual_sunlight_z", p_value][1L], 5)
)
log_msg(
  "Global latitude modifier Wald p-value after baseline sunlight and annual Tmean adjustment: ",
  signif(modifier_test[modifier == "latitude_z", p_value][1L], 5)
)
log_msg(
  "Global annual Tmean modifier Wald p-value after baseline sunlight and latitude adjustment: ",
  signif(modifier_test[modifier == "annual_tmean_z", p_value][1L], 5)
)
log_msg("Within A=0 vs A=-2.768 cumulative OR: ", round(within_target$OR[1L], 4))
log_msg("All outputs written to root: ", OUT_ROOT)
cat(OUT_ROOT, "\n")
