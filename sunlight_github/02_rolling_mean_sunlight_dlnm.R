# Two-stage DLNM for 5-day rolling-mean sunlight, adjusted for national COVID-19 phases.
#
# Main analysis: rolling mean sunlight exposure-response.
# It fits city-specific DLNMs, reduces each city to the overall cumulative
# exposure-response association with crossreduce(type = "overall"), and pools
# those reduced coefficients in a second-stage random-effects meta-analysis.
#
# Stage 1: city-specific quasibinomial DLNMs for lag0-4 rolling-mean sunlight.
#          Temperature is controlled with a same-day natural spline. Humidity,
#          precipitation, and wind speed are same-day linear controls. A
#          four-level pandemic phase factor is included.
# Stage 2: intercept-only multivariate random-effects meta-analysis (REML).
#          No second-stage meta-regression or effect modifier is fitted.
#
# Main outputs:
#   1. National cumulative rolling-mean sunlight exposure-response curve from
#      crossreduce(type = "overall") -> meta_overall.
#   2. Exposure-response curves at every lag day.
#   3. Nationwide rolling-mean P50-versus-P10 lag-response summaries.

suppressPackageStartupMessages({
  library(data.table)
  library(dlnm)
  library(splines)
  library(mixmeta)
})

options(stringsAsFactors = FALSE)
setDTthreads(max(1L, parallel::detectCores(logical = TRUE) - 1L))

# -----------------------------------------------------------------------------
# 1. Paths and fixed analysis settings
# -----------------------------------------------------------------------------

# Input files are restricted-use and must not be committed to a public repository.
# Supply the project root or override individual locations with environment variables.
PROJECT_ROOT <- Sys.getenv("SUNLIGHT_PROJECT_ROOT", unset = ".")
DATA_DIR <- Sys.getenv("SUNLIGHT_DATA_DIR", unset = file.path(PROJECT_ROOT, "data", "restricted"))
OUTPUT_ROOT <- Sys.getenv("SUNLIGHT_OUTPUT_DIR", unset = file.path(PROJECT_ROOT, "outputs"))
CITY_DAY_PATH <- Sys.getenv("CITY_DAY_PATH", unset = file.path(DATA_DIR, "city_day_high_stress_from_personal_p90.csv"))
WEATHER_PATH <- Sys.getenv("WEATHER_PATH", unset = file.path(DATA_DIR, "weather_sunlight_cleaned_city_day.csv"))
SUNLIGHT_PATH <- Sys.getenv("SUNLIGHT_PATH", unset = file.path(DATA_DIR, "all_cities_daily_sunhours_cloud_aerosol_raw_totexttau_with_weather_20260720.csv"))

STUDY_START <- as.IDate("2021-12-01")
STUDY_END <- as.IDate("2024-12-31")

parse_numeric_env <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(default)
  values <- suppressWarnings(as.numeric(strsplit(raw, ",", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!is.finite(values))) {
    stop("Invalid numeric environment variable ", name, ": ", raw)
  }
  values
}

parse_optional_numeric_env <- function(name) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(numeric(0))
  values <- suppressWarnings(as.numeric(strsplit(raw, ",", fixed = TRUE)[[1L]]))
  if (!length(values) || any(!is.finite(values))) {
    stop("Invalid numeric environment variable ", name, ": ", raw)
  }
  values
}

parse_logical_env <- function(name, default = FALSE) {
  raw <- tolower(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(raw)) return(default)
  if (raw %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (raw %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Invalid logical environment variable ", name, ": ", raw)
}

LAG_SUN <- as.integer(parse_numeric_env("LAG_SUN", 28L)[1L])
ROLLING_MEAN_WINDOW <- as.integer(parse_numeric_env("ROLLING_MEAN_WINDOW", 5L)[1L])
if (!is.finite(ROLLING_MEAN_WINDOW) || ROLLING_MEAN_WINDOW < 1L) {
  stop("ROLLING_MEAN_WINDOW must be a positive integer.")
}
ROLLING_MEAN_LAG_END <- ROLLING_MEAN_WINDOW - 1L
EXPOSURE_COL <- paste0("mean_sunlight_w", ROLLING_MEAN_WINDOW)
EXPOSURE_LABEL <- paste0(
  "lag0-", ROLLING_MEAN_LAG_END,
  " rolling-mean sunlight"
)
SUN_VAR_KNOT_PROBS <- parse_numeric_env(
  "SUN_VAR_KNOT_PROBS",
  c(0.333, 0.667)
)
N_SUN_LAG_KNOTS <- as.integer(parse_numeric_env("N_SUN_LAG_KNOTS", 4L)[1L])
SUN_LAG_KNOTS_OVERRIDE <- parse_numeric_env(
  "SUN_LAG_KNOTS_OVERRIDE",
  c(3, 7, 14, 21)
)
format_knot_id <- function(x) {
  z <- format(x, trim = TRUE, scientific = FALSE)
  z <- gsub("\\.", "p", z)
  paste(z, collapse = "_")
}
default_knot_grid_id <- if (length(SUN_LAG_KNOTS_OVERRIDE)) {
  paste0("rolling_mean", ROLLING_MEAN_WINDOW, "_lagknots_", format_knot_id(SUN_LAG_KNOTS_OVERRIDE))
} else {
  paste0("rolling_mean", ROLLING_MEAN_WINDOW, "_logknots_k", N_SUN_LAG_KNOTS)
}
KNOT_GRID_ID <- Sys.getenv("KNOT_GRID_ID", unset = default_knot_grid_id)
TRIM_SUNLIGHT_ABOVE_FULL_P99 <- parse_logical_env("TRIM_SUNLIGHT_ABOVE_FULL_P99", FALSE)
FIT_SUPPLEMENTAL_FULL_CROSSBASIS_META <- parse_logical_env(
  "FIT_SUPPLEMENTAL_FULL_CROSSBASIS_META",
  FALSE
)
SUN_REFERENCE_MODE <- "p10"

PROVISIONAL_CEN_PROB <- 0.50
MIN_CEN_SEARCH_PROBS <- c(0.01, 0.99)
MIN_CEN_GRID_POINTS <- 5001L
LOW_LIGHT_PROB <- 0.10

TIME_DF_PER_YEAR <- 7L
TIME_DF <- TIME_DF_PER_YEAR * 3L
TEMP_NS_DF <- as.integer(parse_numeric_env("TEMP_NS_DF", 3L)[1L])

MIN_MODEL_DAYS <- 365L
MIN_TOTAL_EVENTS <- 20L
MIN_TOTAL_OBSERVED <- 100L
MIXMETA_MAXITER <- as.integer(parse_numeric_env("MIXMETA_MAXITER", 1000L)[1L])
MIXMETA_CONTROL <- mixmeta.control(
  maxiter = MIXMETA_MAXITER,
  optim = list(maxit = MIXMETA_MAXITER)
)

SELECTED_LAGS <- unique(c(0L, 1L, 3L, 5L, 7L, 10L, 14L, 21L, 28L, LAG_SUN))
SELECTED_LAGS <- SELECTED_LAGS[SELECTED_LAGS <= LAG_SUN]
LAG_LABEL <- paste0("lag0_", LAG_SUN)
LAG_TITLE <- paste0("Lag 0-", LAG_SUN)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(OUTPUT_ROOT, paste0("rolling_sunlight_", KNOT_GRID_ID, "_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOG_FILE <- file.path(OUT_DIR, "run_log.txt")

# -----------------------------------------------------------------------------
# 2. Helpers
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

safe_sum <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
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

or_p_value <- function(or, low, high) {
  out <- rep(NA_real_, length(or))
  ok <- is.finite(or) & is.finite(low) & is.finite(high) &
    or > 0 & low > 0 & high > low
  if (!any(ok)) return(out)
  se <- (log(high[ok]) - log(low[ok])) / (2 * qnorm(0.975))
  z <- log(or[ok]) / se
  out[ok] <- 2 * pnorm(abs(z), lower.tail = FALSE)
  out
}

make_significance <- function(low, high) {
  fifelse(
    !is.finite(low) | !is.finite(high),
    NA_character_,
    fifelse(low > 1 | high < 1, "significant", "not significant")
  )
}

predict_reduced_curve <- function(
    model,
    at_grid,
    cen,
    knots,
    boundary,
    curve_type,
    lag_day = NA_integer_) {
  basis_at <- onebasis(
    at_grid,
    fun = "ns",
    knots = knots,
    Boundary.knots = boundary,
    intercept = FALSE
  )
  basis_cen <- onebasis(
    cen,
    fun = "ns",
    knots = knots,
    Boundary.knots = boundary,
    intercept = FALSE
  )
  contrast <- sweep(basis_at, 2L, as.numeric(basis_cen), FUN = "-")

  beta <- coef(model)
  V <- vcov(model)
  log_or <- as.numeric(contrast %*% beta)
  var_log_or <- rowSums((contrast %*% V) * contrast)
  se_log_or <- sqrt(pmax(var_log_or, 0))

  out <- data.table(
    curve_type = curve_type,
    lag_day = as.integer(lag_day),
    sunlight_duration = at_grid,
    reference_sunlight = cen,
    log_OR = log_or,
    se_log_OR = se_log_or,
    OR = exp(log_or),
    OR_low = exp(log_or - qnorm(0.975) * se_log_or),
    OR_high = exp(log_or + qnorm(0.975) * se_log_or)
  )
  out[, p_value := or_p_value(OR, OR_low, OR_high)]
  out[, significance := make_significance(OR_low, OR_high)]
  out
}

predict_var_lag_response <- function(
    model,
    lag_days,
    lag_knots,
    max_lag,
    curve_type,
    exposure_sunlight,
    reference_sunlight,
    cumulative = FALSE) {
  basis_lag <- onebasis(
    lag_days,
    fun = "ns",
    knots = lag_knots,
    Boundary.knots = c(0, max_lag),
    intercept = TRUE
  )
  if (isTRUE(cumulative)) {
    contrast <- do.call(rbind, lapply(seq_len(nrow(basis_lag)), function(i) {
      colSums(basis_lag[seq_len(i), , drop = FALSE])
    }))
  } else {
    contrast <- basis_lag
  }

  beta <- coef(model)
  V <- vcov(model)
  if (ncol(contrast) != length(beta)) {
    stop(
      "Lag basis dimension mismatch: contrast has ", ncol(contrast),
      " columns but meta model has ", length(beta), " coefficients."
    )
  }

  log_or <- as.numeric(contrast %*% beta)
  var_log_or <- rowSums((contrast %*% V) * contrast)
  se_log_or <- sqrt(pmax(var_log_or, 0))

  out <- data.table(
    curve_type = curve_type,
    lag_day = as.integer(lag_days),
    exposure_label = "nationwide_P50",
    reference_label = "nationwide_P10",
    exposure_sunlight = exposure_sunlight,
    reference_sunlight = reference_sunlight,
    log_OR = log_or,
    se_log_OR = se_log_or,
    OR = exp(log_or),
    OR_low = exp(log_or - qnorm(0.975) * se_log_or),
    OR_high = exp(log_or + qnorm(0.975) * se_log_or)
  )
  out[, p_value := or_p_value(OR, OR_low, OR_high)]
  out[, significance := make_significance(OR_low, OR_high)]
  out
}

find_minimum_risk_cen <- function(
    model,
    knots,
    boundary,
    sun_values,
    search_probs = c(0.01, 0.99),
    grid_points = 5001L) {
  limits <- as.numeric(quantile(
    sun_values,
    search_probs,
    na.rm = TRUE,
    type = 7
  ))
  search_grid <- seq(limits[1L], limits[2L], length.out = grid_points)
  search_basis <- onebasis(
    search_grid,
    fun = "ns",
    knots = knots,
    Boundary.knots = boundary,
    intercept = FALSE
  )
  linear_predictor <- as.numeric(search_basis %*% coef(model))
  minimum_index <- which.min(linear_predictor)

  list(
    cen = search_grid[minimum_index],
    search = data.table(
      sunlight_duration = search_grid,
      relative_predicted_odds = exp(linear_predictor - min(linear_predictor)),
      is_minimum = seq_along(search_grid) == minimum_index
    ),
    summary = data.table(
      search_lower_percentile = 100 * search_probs[1L],
      search_upper_percentile = 100 * search_probs[2L],
      search_lower_hours = limits[1L],
      search_upper_hours = limits[2L],
      minimum_sunlight_hours = search_grid[minimum_index],
      empirical_percentile = 100 * mean(sun_values <= search_grid[minimum_index]),
      boundary_flag = minimum_index %in% c(1L, length(search_grid)),
      grid_index = minimum_index,
      grid_points = length(search_grid)
    )
  )
}

nearest_curve_rows <- function(curve, selected_values, selected_labels) {
  rbindlist(lapply(seq_along(selected_values), function(i) {
    idx <- which.min(abs(curve$sunlight_duration - selected_values[i]))
    curve[idx, .(
      curve_type,
      lag_day,
      exposure_label = selected_labels[i],
      requested_sunlight = selected_values[i],
      predicted_sunlight = sunlight_duration,
      reference_sunlight,
      OR,
      OR_low,
      OR_high,
      p_value,
      significance
    )]
  }))
}

predict_full_cumulative_curve <- function(
    model,
    at_grid,
    cen,
    curve_type) {
  cb_sun_template <- crossbasis(
    at_grid,
    lag = LAG_SUN,
    argvar = list(
      fun = "ns",
      knots = sun_knots,
      Boundary.knots = sun_boundary,
      intercept = FALSE
    ),
    arglag = list(
      fun = "ns",
      knots = sun_lag_knots,
      intercept = TRUE
    )
  )

  pred <- crosspred(
    cb_sun_template,
    coef = as.numeric(coef(model)),
    vcov = vcov(model),
    model.link = "logit",
    at = at_grid,
    cen = cen
  )

  out <- data.table(
    curve_type = curve_type,
    lag_day = NA_integer_,
    sunlight_duration = at_grid,
    reference_sunlight = cen,
    log_OR = as.numeric(pred$allfit),
    se_log_OR = as.numeric(pred$allse),
    OR = as.numeric(pred$allRRfit),
    OR_low = as.numeric(pred$allRRlow),
    OR_high = as.numeric(pred$allRRhigh)
  )
  out[, p_value := or_p_value(OR, OR_low, OR_high)]
  out[, significance := make_significance(OR_low, OR_high)]
  out
}

# -----------------------------------------------------------------------------
# 3. Read cached city-day outcome and weather data
# -----------------------------------------------------------------------------

log_msg("Output directory: ", OUT_DIR)
log_msg("City-day outcome: ", CITY_DAY_PATH)
log_msg("Weather controls: ", WEATHER_PATH)
log_msg("Sunlight source: ", SUNLIGHT_PATH, " [daily_sun_hours]")
log_msg("Stage 1: quasibinomial, weights=n_observed")
log_msg(
  "Rolling-mean sunlight DLNM lag 0-", LAG_SUN,
  "; exposure = ", EXPOSURE_LABEL,
  "; reference mode = ", SUN_REFERENCE_MODE
)
log_msg("Temperature: same-day natural spline control, df=", TEMP_NS_DF, "; no temperature DLNM")
log_msg("Stage 2: intercept-only multivariate random-effects meta-analysis; no meta-regression")

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
sunlight_new <- fread(
  SUNLIGHT_PATH,
  colClasses = list(character = "CITYCODE"),
  showProgress = FALSE
)

city_day[, city_code := clean_city_code(city_code)]
city_day[, date := as.IDate(date)]
city_day <- city_day[
  !is.na(city_code) & !is.na(date) &
    date >= STUDY_START & date <= STUDY_END
]
city_day <- city_day[, .(
  high_events = safe_sum(as.numeric(high_events)),
  n_observed = safe_sum(as.numeric(n_observed))
), by = .(city_code, date)]
city_day[, high_stress_rate := fifelse(
  is.finite(n_observed) & n_observed > 0,
  high_events / n_observed,
  NA_real_
)]

weather[, city_code := clean_city_code(city_code)]
weather[, date := as.IDate(date)]
weather <- weather[
  !is.na(city_code) & !is.na(date) &
    date >= STUDY_START & date <= STUDY_END
]
weather <- weather[, .(
  city_name = {
    z <- as.character(city_name)
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z)) z[1L] else NA_character_
  },
  temperature = safe_mean(as.numeric(temperature)),
  humidity = safe_mean(as.numeric(humidity)),
  precipitation = safe_mean(as.numeric(precipitation)),
  wind_speed = safe_mean(as.numeric(wind_speed))
), by = .(city_code, date)]

sunlight_new[, city_code := clean_city_code(CITYCODE)]
sunlight_new[, date := as.IDate(date)]
sunlight_new <- sunlight_new[
  !is.na(city_code) & !is.na(date) &
    date >= STUDY_START & date <= STUDY_END
]
sunlight_new <- sunlight_new[, .(
  sunlight_duration = safe_mean(as.numeric(daily_sun_hours))
), by = .(city_code, date)]

model_cities <- sort(Reduce(
  intersect,
  list(unique(city_day$city_code), unique(weather$city_code), unique(sunlight_new$city_code))
))
all_dates <- as.IDate(seq(as.Date(STUDY_START), as.Date(STUDY_END), by = "day"))
panel <- CJ(city_code = model_cities, date = all_dates, unique = TRUE)
panel <- merge(panel, city_day, by = c("city_code", "date"), all.x = TRUE, sort = FALSE)
panel <- merge(panel, weather, by = c("city_code", "date"), all.x = TRUE, sort = FALSE)
panel <- merge(panel, sunlight_new, by = c("city_code", "date"), all.x = TRUE, sort = FALSE)
setorder(panel, city_code, date)

panel[, failures := n_observed - high_events]
panel[, time_num := as.numeric(date - STUDY_START) + 1]
panel[, dow := factor(
  as.POSIXlt(as.Date(date))$wday,
  levels = 0:6,
  labels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
)]
panel[, pandemic_phase := make_pandemic_phase(date)]

full_sun_global <- panel[is.finite(sunlight_duration), sunlight_duration]
full_p99_sunlight_cutoff <- as.numeric(quantile(
  full_sun_global,
  0.99,
  na.rm = TRUE,
  type = 7
))
rows_before_p99_trim <- nrow(panel)
rows_above_full_p99 <- panel[
  is.finite(sunlight_duration) & sunlight_duration > full_p99_sunlight_cutoff,
  .N
]
values_winsorized_above_full_p99 <- 0L

if (TRIM_SUNLIGHT_ABOVE_FULL_P99) {
  values_winsorized_above_full_p99 <- rows_above_full_p99
  panel[
    is.finite(sunlight_duration) & sunlight_duration > full_p99_sunlight_cutoff,
    sunlight_duration := full_p99_sunlight_cutoff
  ]
  log_msg(
    "Applied sunlight P99 winsorization: capped ", rows_above_full_p99,
    " sunlight values above ", round(full_p99_sunlight_cutoff, 6),
    " h/day. No dates were removed."
  )
} else {
  log_msg(
    "Sunlight P99 trim not applied. Full-sample P99 = ",
    round(full_p99_sunlight_cutoff, 6), " h/day."
  )
}

save_csv(
  data.table(
    trim_applied = TRIM_SUNLIGHT_ABOVE_FULL_P99,
    full_sample_p99_sunlight_hours = full_p99_sunlight_cutoff,
    rows_before_trim = rows_before_p99_trim,
    rows_above_full_p99 = rows_above_full_p99,
    values_winsorized_above_full_p99 = values_winsorized_above_full_p99,
    rows_after_trim = nrow(panel),
    rows_removed = rows_before_p99_trim - nrow(panel),
    trim_method = if (TRIM_SUNLIGHT_ABOVE_FULL_P99) {
      "winsorized sunlight_duration to full-sample P99; no city-date rows removed"
    } else {
      "not applied"
    }
  ),
  "sunlight_full_p99_trim_summary.csv"
)

panel[, (EXPOSURE_COL) := data.table::frollmean(
  sunlight_duration,
  n = ROLLING_MEAN_WINDOW,
  align = "right",
  na.rm = FALSE
), by = city_code]

# Full city-day analytic panel is kept in memory; it is not exported.

phase_definition <- data.table(
  pandemic_phase = levels(panel$pandemic_phase),
  start_date = c("2021-12-01", "2022-11-11", "2022-12-07", "2023-01-23"),
  end_date = c("2022-11-10", "2022-12-06", "2023-01-22", "2024-12-31"),
  reference_level = c(TRUE, FALSE, FALSE, FALSE)
)
save_csv(phase_definition, "pandemic_phase_definition.csv")

phase_counts <- panel[
  !is.na(high_events) & !is.na(n_observed) & n_observed > 0,
  .(
    city_days = .N,
    cities = uniqueN(city_code),
    high_events = sum(high_events, na.rm = TRUE),
    total_observed = sum(n_observed, na.rm = TRUE),
    high_stress_rate = sum(high_events, na.rm = TRUE) / sum(n_observed, na.rm = TRUE)
  ),
  by = pandemic_phase
]
save_csv(phase_counts, "pandemic_phase_sample_counts.csv")

# -----------------------------------------------------------------------------
# 4. Global sunlight DLNM parameters
# -----------------------------------------------------------------------------

sun_global <- panel[is.finite(get(EXPOSURE_COL)), get(EXPOSURE_COL)]
sun_boundary <- range(sun_global, na.rm = TRUE)
if (SUN_REFERENCE_MODE == "p99_trim_cutoff") {
  sun_boundary[2L] <- max(sun_boundary[2L], full_p99_sunlight_cutoff, na.rm = TRUE)
}
sun_knots <- as.numeric(quantile(
  sun_global,
  SUN_VAR_KNOT_PROBS,
  na.rm = TRUE,
  type = 7
))
sun_lag_knots <- if (length(SUN_LAG_KNOTS_OVERRIDE)) {
  SUN_LAG_KNOTS_OVERRIDE
} else {
  as.numeric(dlnm::logknots(LAG_SUN, N_SUN_LAG_KNOTS))
}
sun_lag_knots <- sort(unique(as.numeric(sun_lag_knots)))
if (!length(sun_lag_knots) || any(!is.finite(sun_lag_knots))) {
  stop("Sunlight lag knots are empty or non-finite.")
}
if (any(sun_lag_knots <= 0 | sun_lag_knots >= LAG_SUN)) {
  stop("Sunlight lag knots must lie strictly between 0 and LAG_SUN.")
}
provisional_cen <- as.numeric(quantile(
  sun_global,
  PROVISIONAL_CEN_PROB,
  na.rm = TRUE,
  type = 7
))
low_light_threshold <- as.numeric(quantile(
  sun_global,
  LOW_LIGHT_PROB,
  na.rm = TRUE,
  type = 7
))

sun_probs <- c(
  p01 = 0.01,
  p05 = 0.05,
  p10 = 0.10,
  p25 = 0.25,
  p50 = 0.50,
  p75 = 0.75,
  p90 = 0.90,
  p95 = 0.95,
  p99 = 0.99
)
sun_pct <- setNames(
  as.numeric(quantile(sun_global, sun_probs, na.rm = TRUE, type = 7)),
  names(sun_probs)
)

global_parameters <- rbindlist(list(
  data.table(
    parameter = paste0("sunlight_exposure_knot_", seq_along(sun_knots)),
    value = sun_knots
  ),
  data.table(
    parameter = paste0("sunlight_lag_knot_", seq_along(sun_lag_knots)),
    value = sun_lag_knots
  ),
  data.table(
    parameter = c(
      "sunlight_boundary_lower", "sunlight_boundary_upper",
      "sunlight_provisional_reference", "sunlight_low_p10",
      "rolling_mean_window_days", "rolling_mean_window_lag_start",
      "rolling_mean_window_lag_end",
      "sunlight_lag_knot_count", "sunlight_lag_knots_from_explicit_override",
      "sunlight_max_lag", "time_df",
      "temperature_form", "temperature_ns_df", "temperature_max_lag", "temperature_exposure_knots",
      "temperature_lag_knots"
    ),
    value = c(
      sun_boundary,
      provisional_cen,
      low_light_threshold,
      ROLLING_MEAN_WINDOW,
      0,
      ROLLING_MEAN_LAG_END,
      length(sun_lag_knots),
      length(SUN_LAG_KNOTS_OVERRIDE) > 0L,
      LAG_SUN,
      TIME_DF,
      NA,
      TEMP_NS_DF,
      0,
      NA,
      NA
    )
  )
), fill = TRUE)
save_csv(global_parameters, "global_dlnm_parameters.csv")
save_csv(data.table(
  percentile = names(sun_pct),
  probability = as.numeric(sun_probs),
  rolling_mean_sunlight = as.numeric(sun_pct)
), "global_rolling_mean_sunlight_percentiles.csv")

log_msg("Knot-grid specification id: ", KNOT_GRID_ID)
log_msg("Rolling mean exposure: ", EXPOSURE_LABEL, " using column ", EXPOSURE_COL)
log_msg("Global rolling-mean sunlight exposure knots: ", paste(round(sun_knots, 4), collapse = ", "))
log_msg("Global rolling-mean sunlight lag knots: ", paste(round(sun_lag_knots, 4), collapse = ", "))
log_msg(
  "P10 rolling-mean sunlight: ", round(low_light_threshold, 6),
  "; provisional center used during coefficient reduction: ",
  round(provisional_cen, 6)
)

# -----------------------------------------------------------------------------
# 5. Stage 1: city-specific pandemic-adjusted DLNMs
# -----------------------------------------------------------------------------

run_city_model <- function(city_id) {
  dt <- copy(panel[city_code == city_id])
  setorder(dt, date)
  city_label <- dt[!is.na(city_name) & nzchar(city_name), city_name[1L]]
  if (!length(city_label)) city_label <- NA_character_

  base_status <- data.table(
    city_code = city_id,
    city_name = city_label,
    panel_days = nrow(dt),
    complete_days_after_lag = NA_integer_,
    total_observed_after_lag = NA_real_,
    total_high_events_after_lag = NA_real_,
    national_p10_sunlight = low_light_threshold,
    national_p50_sunlight = provisional_cen,
    p50_var_lag_response_status = NA_character_,
    p50_var_lag_response_message = NA_character_,
    fit_status = "not_started",
    message = NA_character_
  )

  tryCatch({
    cb_sun <- crossbasis(
      dt[[EXPOSURE_COL]],
      lag = LAG_SUN,
      argvar = list(
        fun = "ns",
        knots = sun_knots,
        Boundary.knots = sun_boundary,
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
      dt$n_observed > 0 &
      dt$high_events >= 0 &
      dt$failures >= 0

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
        ns(temperature, df = TEMP_NS_DF) +
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

    if (!isTRUE(fit$converged)) {
      stop("glm did not converge")
    }

    full_coef_names <- grep("^cb_sun", names(coef(fit)), value = TRUE)
    full_coef <- coef(fit)[full_coef_names]
    full_vcov <- vcov(fit)[full_coef_names, full_coef_names, drop = FALSE]
    if (any(!is.finite(full_coef)) || any(!is.finite(full_vcov))) {
      stop("non-finite full cross-basis coefficient or covariance")
    }

    reduced_overall <- crossreduce(
      cb_sun,
      fit,
      type = "overall",
      cen = provisional_cen
    )
    overall_coef <- coef(reduced_overall)
    overall_vcov <- vcov(reduced_overall)
    if (any(!is.finite(overall_coef)) || any(!is.finite(overall_vcov))) {
      stop("non-finite overall reduced coefficient or covariance")
    }

    daily_coef <- matrix(
      NA_real_,
      nrow = LAG_SUN + 1L,
      ncol = length(overall_coef),
      dimnames = list(paste0("lag", 0:LAG_SUN), names(overall_coef))
    )
    daily_vcov <- vector("list", LAG_SUN + 1L)

    for (lag_day in 0:LAG_SUN) {
      reduced_lag <- crossreduce(
        cb_sun,
        fit,
        type = "lag",
        value = lag_day,
        cen = provisional_cen
      )
      lag_coef <- coef(reduced_lag)
      lag_vcov <- vcov(reduced_lag)
      if (any(!is.finite(lag_coef)) || any(!is.finite(lag_vcov))) {
        stop("non-finite reduced result at lag ", lag_day)
      }
      daily_coef[lag_day + 1L, ] <- lag_coef
      daily_vcov[[lag_day + 1L]] <- lag_vcov
    }

    if (!is.finite(low_light_threshold) || !is.finite(provisional_cen) ||
      low_light_threshold >= provisional_cen) {
      stop("nationwide P10/P50 is non-finite or P10 >= P50")
    }
    reduced_var_p50 <- crossreduce(
      cb_sun,
      fit,
      type = "var",
      value = provisional_cen,
      cen = low_light_threshold
    )
    var_p50_coef <- coef(reduced_var_p50)
    var_p50_vcov <- vcov(reduced_var_p50)
    if (any(!is.finite(var_p50_coef)) || any(!is.finite(var_p50_vcov))) {
      stop("non-finite P50-vs-P10 lag-response reduced coefficient or covariance")
    }

    base_status[, `:=`(
      fit_status = "ok",
      model_days_used = nobs(fit),
      dispersion = summary(fit)$dispersion,
      glm_iterations = fit$iter,
      glm_converged = TRUE,
      p50_var_lag_response_status = "ok",
      p50_var_lag_response_message = NA_character_,
      message = NA_character_
    )]

    list(
      status = base_status,
      full_coef = full_coef,
      full_vcov = full_vcov,
      overall_coef = overall_coef,
      overall_vcov = overall_vcov,
      daily_coef = daily_coef,
      daily_vcov = daily_vcov,
      var_p50_coef = var_p50_coef,
      var_p50_vcov = var_p50_vcov
    )
  }, error = function(e) {
    base_status[, `:=`(
      fit_status = "failed",
      message = conditionMessage(e)
    )]
    list(status = base_status)
  })
}

log_msg("Fitting city-specific models for ", length(model_cities), " cities...")
stage1_results <- vector("list", length(model_cities))
names(stage1_results) <- model_cities

for (i in seq_along(model_cities)) {
  if (i == 1L || i %% 20L == 0L || i == length(model_cities)) {
    log_msg("Stage-1 progress: ", i, "/", length(model_cities))
  }
  stage1_results[[i]] <- run_city_model(model_cities[i])
}

stage1_status <- rbindlist(lapply(stage1_results, `[[`, "status"), fill = TRUE)
save_csv(stage1_status, "stage1_city_model_status.csv")

ok_cities <- stage1_status[fit_status == "ok", city_code]
if (length(ok_cities) < 5L) {
  stop("Too few successful city-specific models: ", length(ok_cities))
}

overall_coef_list <- lapply(stage1_results[ok_cities], `[[`, "overall_coef")
overall_vcov_list <- lapply(stage1_results[ok_cities], `[[`, "overall_vcov")
overall_coef_matrix <- do.call(rbind, overall_coef_list)
rownames(overall_coef_matrix) <- ok_cities

full_coef_list <- lapply(stage1_results[ok_cities], `[[`, "full_coef")
full_vcov_list <- lapply(stage1_results[ok_cities], `[[`, "full_vcov")
full_coef_matrix <- do.call(rbind, full_coef_list)
rownames(full_coef_matrix) <- ok_cities

save_csv(
  as.data.table(overall_coef_matrix, keep.rownames = "city_code"),
  "stage1_overall_reduced_coefficients.csv"
)
save_csv(
  as.data.table(full_coef_matrix, keep.rownames = "city_code"),
  "stage1_full_crossbasis_coefficients.csv"
)

log_msg("Stage-1 successful cities: ", length(ok_cities))
log_msg("Stage-1 failed or skipped cities: ", nrow(stage1_status) - length(ok_cities))

# -----------------------------------------------------------------------------
# 6. Stage 2: intercept-only random-effects meta-analysis
# -----------------------------------------------------------------------------

log_msg("Fitting cumulative intercept-only multivariate random-effects meta-analysis...")
meta_overall <- mixmeta(
  overall_coef_matrix ~ 1,
  S = overall_vcov_list,
  method = "reml",
  control = MIXMETA_CONTROL
)
capture.output(
  summary(meta_overall),
  file = file.path(OUT_DIR, "stage2_cumulative_random_effects_summary.txt")
)

meta_full <- NULL
full_valid <- apply(full_coef_matrix, 1L, function(x) all(is.finite(x))) &
  vapply(full_vcov_list, function(x) all(is.finite(x)), logical(1))

if (FIT_SUPPLEMENTAL_FULL_CROSSBASIS_META && sum(full_valid) >= 5L) {
  log_msg("Fitting supplemental full cross-basis random-effects meta-analysis...")
  meta_full <- tryCatch(
    mixmeta(
      full_coef_matrix[full_valid, , drop = FALSE] ~ 1,
      S = full_vcov_list[full_valid],
      method = "reml",
      control = MIXMETA_CONTROL
    ),
    error = function(e) {
      log_msg("WARNING: full cross-basis meta-analysis failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(meta_full)) {
    capture.output(
      summary(meta_full),
      file = file.path(OUT_DIR, "stage2_full_crossbasis_random_effects_summary.txt")
    )
  }
} else if (!FIT_SUPPLEMENTAL_FULL_CROSSBASIS_META) {
  log_msg(
    "Skipping supplemental full cross-basis meta-analysis. ",
    "The main cumulative exposure-response curve uses meta_overall from crossreduce(type='overall')."
  )
} else {
  log_msg(
    "WARNING: too few valid full cross-basis city models for supplemental meta-analysis: ",
    sum(full_valid)
  )
}

save_csv(
  data.table(
    cities_available = sum(full_valid),
    requested = FIT_SUPPLEMENTAL_FULL_CROSSBASIS_META,
    fit_status = if (!is.null(meta_full)) "ok" else "failed_or_skipped"
  ),
  "stage2_full_crossbasis_meta_status.csv"
)

minimum_cen_result <- find_minimum_risk_cen(
  meta_overall,
  sun_knots,
  sun_boundary,
  sun_global,
  search_probs = MIN_CEN_SEARCH_PROBS,
  grid_points = MIN_CEN_GRID_POINTS
)
sun_cen <- unname(sun_pct[["p10"]])
reference_label <- "P10 low-sunlight reference"
reference_detail <- paste0("P10=", round(sun_cen, 6), " hours/day")
save_csv(minimum_cen_result$search, "minimum_risk_cen_search_P1_P99.csv")
save_csv(minimum_cen_result$summary, "minimum_risk_cen_summary.csv")
save_csv(
  data.table(
    reference_mode = SUN_REFERENCE_MODE,
    reference_label = reference_label,
    reference_sunlight_hours = sun_cen,
    minimum_risk_sunlight_hours = minimum_cen_result$cen,
    full_sample_p99_sunlight_hours = full_p99_sunlight_cutoff
  ),
  "sunlight_reference_summary.csv"
)

log_msg(
  "Minimum-risk cen within P1-P99 = ", round(minimum_cen_result$cen, 6),
  " h/day (empirical percentile ",
  round(minimum_cen_result$summary$empirical_percentile, 3),
  "; boundary flag = ", minimum_cen_result$summary$boundary_flag, ")"
)
log_msg(
  "Using sunlight reference mode '", SUN_REFERENCE_MODE,
  "': ", reference_detail
)
if (isTRUE(minimum_cen_result$summary$boundary_flag)) {
  log_msg(
    "WARNING: the minimum-risk cen is at the P1-P99 search boundary; ",
    "interpret the selected reference cautiously."
  )
}

at_grid <- sort(unique(c(
  seq(sun_boundary[1L], sun_boundary[2L], length.out = 240L),
  as.numeric(sun_pct),
  sun_cen
)))

cumulative_curve <- predict_reduced_curve(
  meta_overall,
  at_grid,
  sun_cen,
  sun_knots,
  sun_boundary,
  curve_type = paste0("main_cumulative_overall_reduced_", LAG_LABEL)
)
save_csv(
  cumulative_curve,
  paste0("meta_cumulative_exposure_response_curve_", LAG_LABEL, "_vs_reference.csv")
)

daily_meta_models <- vector("list", LAG_SUN + 1L)
daily_curve_list <- vector("list", LAG_SUN + 1L)
daily_meta_status <- vector("list", LAG_SUN + 1L)

for (lag_day in 0:LAG_SUN) {
  if (lag_day == 0L || lag_day %% 5L == 0L || lag_day == LAG_SUN) {
    log_msg("Stage-2 daily lag meta-analysis: lag ", lag_day, "/", LAG_SUN)
  }

  lag_coef_matrix <- do.call(rbind, lapply(
    stage1_results[ok_cities],
    function(x) x$daily_coef[lag_day + 1L, , drop = FALSE]
  ))
  rownames(lag_coef_matrix) <- ok_cities
  lag_vcov_list <- lapply(
    stage1_results[ok_cities],
    function(x) x$daily_vcov[[lag_day + 1L]]
  )

  valid <- apply(lag_coef_matrix, 1L, function(x) all(is.finite(x))) &
    vapply(lag_vcov_list, function(x) all(is.finite(x)), logical(1))

  fit_lag <- mixmeta(
    lag_coef_matrix[valid, , drop = FALSE] ~ 1,
    S = lag_vcov_list[valid],
    method = "reml",
    control = MIXMETA_CONTROL
  )
  daily_meta_models[[lag_day + 1L]] <- fit_lag
  daily_curve_list[[lag_day + 1L]] <- predict_reduced_curve(
    fit_lag,
    at_grid,
    sun_cen,
    sun_knots,
    sun_boundary,
    curve_type = "lag_specific",
    lag_day = lag_day
  )
  daily_meta_status[[lag_day + 1L]] <- data.table(
    lag_day = lag_day,
    cities_available = sum(valid),
    fit_status = "ok"
  )
}

daily_curves <- rbindlist(daily_curve_list)
daily_meta_status <- rbindlist(daily_meta_status)
save_csv(
  daily_curves,
  paste0("meta_daily_exposure_response_curves_", LAG_LABEL, "_vs_reference.csv")
)
save_csv(daily_meta_status, "stage2_daily_lag_meta_status.csv")

p50_var_lag_response_status <- stage1_status[, .(
  city_code,
  city_name,
  fit_status,
  p50_var_lag_response_status,
  p50_var_lag_response_message,
  complete_days_after_lag,
  national_p10_sunlight,
  national_p50_sunlight
)]
save_csv(
  p50_var_lag_response_status,
  "stage1_nationwide_P50_P10_var_reduction_status.csv"
)

var_p50_coef_list <- lapply(stage1_results[ok_cities], `[[`, "var_p50_coef")
var_p50_vcov_list <- lapply(stage1_results[ok_cities], `[[`, "var_p50_vcov")
var_p50_valid <- vapply(var_p50_coef_list, function(x) {
  !is.null(x) && all(is.finite(x))
}, logical(1)) &
  vapply(var_p50_vcov_list, function(x) {
    !is.null(x) && all(is.finite(x))
  }, logical(1))

if (sum(var_p50_valid) < 5L) {
  stop("Too few valid crossreduce(type='var') P50-vs-P10 city reductions: ", sum(var_p50_valid))
}

var_p50_coef_matrix <- do.call(rbind, var_p50_coef_list[var_p50_valid])
rownames(var_p50_coef_matrix) <- ok_cities[var_p50_valid]
save_csv(
  as.data.table(var_p50_coef_matrix, keep.rownames = "city_code"),
  "stage1_nationwide_P50_vs_P10_var_reduced_coefficients.csv"
)

log_msg("Fitting nationwide P50 versus nationwide P10 lag-response multivariate meta-analysis...")
meta_var_p50 <- mixmeta(
  var_p50_coef_matrix ~ 1,
  S = var_p50_vcov_list[var_p50_valid],
  method = "reml",
  control = MIXMETA_CONTROL
)
capture.output(
  summary(meta_var_p50),
  file = file.path(OUT_DIR, "stage2_var_P50_vs_P10_meta_summary.txt")
)

p50_lag_response <- predict_var_lag_response(
  meta_var_p50,
  lag_days = 0:LAG_SUN,
  lag_knots = sun_lag_knots,
  max_lag = LAG_SUN,
  curve_type = "meta_var_nationwide_P50_vs_P10_daily_lag_response",
  exposure_sunlight = provisional_cen,
  reference_sunlight = low_light_threshold,
  cumulative = FALSE
)
p50_lag_response[, p_value_BH := p.adjust(p_value, method = "BH")]
p50_lag_response[, significant_BH := is.finite(p_value_BH) & p_value_BH < 0.05]
save_csv(
  p50_lag_response,
  "meta_nationwide_P50_vs_P10_daily_lag_response_OR.csv"
)

p50_cumulative_lag_response <- predict_var_lag_response(
  meta_var_p50,
  lag_days = 0:LAG_SUN,
  lag_knots = sun_lag_knots,
  max_lag = LAG_SUN,
  curve_type = paste0("meta_var_nationwide_P50_vs_P10_cumulative_", LAG_LABEL),
  exposure_sunlight = provisional_cen,
  reference_sunlight = low_light_threshold,
  cumulative = TRUE
)
p50_cumulative_lag_response[, p_value_BH := p.adjust(p_value, method = "BH")]
p50_cumulative_lag_response[, significant_BH := is.finite(p_value_BH) & p_value_BH < 0.05]
save_csv(
  p50_cumulative_lag_response,
  "meta_nationwide_P50_vs_P10_incremental_cumulative_lag_OR.csv"
)

saveRDS(
  list(
    meta_overall = meta_overall,
    meta_full = meta_full,
    daily_meta_models = daily_meta_models,
    meta_var_p50 = meta_var_p50,
    p50_lag_response = p50_lag_response,
    p50_cumulative_lag_response = p50_cumulative_lag_response,
    ok_cities = ok_cities,
    ok_cities_full_crossbasis = ok_cities[full_valid],
    sun_knots = sun_knots,
    sun_lag_knots = sun_lag_knots,
    sun_boundary = sun_boundary,
    sun_cen = sun_cen,
    sun_percentiles = sun_pct
  ),
  file.path(OUT_DIR, "stage2_intercept_only_meta_models.rds")
)

# -----------------------------------------------------------------------------
# 7. Selected contrasts
# -----------------------------------------------------------------------------

selected_labels <- c("P1", "P5", "P10", "P25", "P50", "P75", "P90", "P95", "P99")
selected_values <- as.numeric(sun_pct[c(
  "p01", "p05", "p10", "p25", "p50", "p75", "p90", "p95", "p99"
)])
selected_cumulative <- predict_reduced_curve(
  meta_overall,
  selected_values,
  sun_cen,
  sun_knots,
  sun_boundary,
  curve_type = paste0("selected_points_from_main_national_cumulative_curve_", LAG_LABEL)
)
selected_cumulative[, percentile_label := selected_labels]
selected_cumulative[, lag_day := LAG_SUN]
setcolorder(
  selected_cumulative,
  c("percentile_label", setdiff(names(selected_cumulative), "percentile_label"))
)
save_csv(
  selected_cumulative,
  paste0("meta_cumulative_selected_ORs_", LAG_LABEL, "_vs_reference.csv")
)

p50_from_national_cumulative_curve <- copy(
  selected_cumulative[percentile_label == "P50"]
)
p50_from_national_cumulative_curve[, `:=`(
  curve_type = paste0("P50_vs_P10_from_main_national_cumulative_curve_", LAG_LABEL),
  exposure_label = "nationwide_P50",
  reference_label = "nationwide_P10",
  exposure_rolling_mean_sunlight = unname(sun_pct[["p50"]]),
  reference_rolling_mean_sunlight = sun_cen,
  contrast_source = "main national cumulative exposure-response curve: crossreduce(type='overall') -> meta_overall"
)]
save_csv(
  p50_from_national_cumulative_curve,
  paste0("meta_P50_vs_P10_from_national_cumulative_curve_", LAG_LABEL, "_OR.csv")
)

save_csv(
  p50_lag_response,
  "meta_P50_vs_P10_daily_lag_response_OR.csv"
)

save_csv(
  p50_cumulative_lag_response,
  "meta_P50_vs_P10_incremental_cumulative_lag_OR.csv"
)

# Figure creation omitted from this public analysis script.

# -----------------------------------------------------------------------------
# 8. Diagnostics and model specification
# -----------------------------------------------------------------------------

national_contrast_successful_cities <- sum(var_p50_valid)
national_p50_vs_p10_lag28_or <- p50_from_national_cumulative_curve$OR[1L]

model_specification <- data.table(
  item = c(
    "analysis_design",
    "outcome",
    "stage1_family",
    "sunlight_exposure",
    "rolling_mean_window",
    "sunlight_data_source",
    "sunlight_lag",
    "sunlight_exposure_knots",
    "sunlight_lag_knots",
    "sunlight_reference",
    "main_cumulative_curve",
    "P50_P10_contrast_definition",
    "P50_P10_lag_response",
    "temperature_control",
    "other_weather_controls",
    "time_control",
    "pandemic_control",
    "stage2_model",
    "stage2_meta_regression"
  ),
  value = c(
    "two-stage city-specific DLNM plus random-effects meta-analysis",
    "city-day high_events / n_observed",
    "quasibinomial logit weighted by n_observed",
    paste0("X = ", EXPOSURE_LABEL, " computed from daily_sun_hours after cleaning/winsorization"),
    paste0("lag0-", ROLLING_MEAN_LAG_END, "; align='right'; na.rm=FALSE"),
    SUNLIGHT_PATH,
    paste0("0-", LAG_SUN, " days"),
    paste0(
      "natural spline with internal knots at global P",
      paste(format(round(100 * SUN_VAR_KNOT_PROBS, 1), trim = TRUE), collapse = ", P")
    ),
    paste0(
      "natural spline with lag knots at ",
      paste(round(sun_lag_knots, 4), collapse = ", "), " days"
    ),
    paste0(reference_label, " = ", round(sun_cen, 6), " hours/day of ", EXPOSURE_LABEL),
    "main curve is generated from crossreduce(type='overall') city coefficients pooled by meta_overall; full cross-basis meta is not used for the main cumulative exposure-response curve",
    "main P50-vs-P10 cumulative OR is read directly from the national cumulative exposure-response curve generated by crossreduce(type='overall') -> meta_overall; P10 and P50 are fixed nationwide pooled percentiles of the rolling-mean exposure",
    "P50-vs-P10 lag-response is generated by city crossreduce(type='var', value=nationwide rolling-mean P10, cen=nationwide rolling-mean P50), then pooled once with multivariate mixmeta; 29 separate scalar lag meta-analyses are not used",
    paste0("same-day natural spline temperature, df=", TEMP_NS_DF, "; no temperature DLNM"),
    "same-day linear humidity, precipitation, wind speed",
    paste0("natural spline of time, df=", TIME_DF, ", plus day of week"),
    "four national pandemic phases; phase 1 reference",
    "intercept-only multivariate random-effects meta-analysis, REML",
    "none"
  )
)
save_csv(model_specification, "model_specification.csv")

diagnostics <- data.table(
  metric = c(
    "city_day_outcome_rows",
    "weather_city_day_rows",
    "sunlight_city_day_rows",
    "complete_panel_rows",
    "rolling_mean_exposure",
    "panel_cities",
    "stage1_successful_cities",
    "stage1_failed_or_skipped_cities",
    "study_start",
    "study_end",
    "sunlight_lag",
    "temperature_lag",
    "temperature_form",
    "sunlight_reference",
    "main_cumulative_curve_source",
    "main_P50_P10_OR_source",
    "P50_P10_lag_response_source",
    "nationwide_P50_P10_contrast_successful_cities",
    "nationwide_P50_vs_P10_lag0_28_cumulative_OR",
    "low_light_P10",
    "pandemic_phase_reference",
    "stage2_meta_regression"
  ),
  value = as.character(c(
    nrow(city_day),
    nrow(weather),
    nrow(sunlight_new),
    nrow(panel),
    EXPOSURE_LABEL,
    uniqueN(panel$city_code),
    length(ok_cities),
    nrow(stage1_status) - length(ok_cities),
    as.character(STUDY_START),
    as.character(STUDY_END),
    paste0("0-", LAG_SUN),
    "0",
    paste0("same-day natural spline, df=", TEMP_NS_DF),
    reference_detail,
    "crossreduce(type='overall') -> meta_overall",
    "read from national cumulative exposure-response curve",
    "crossreduce(type='var') -> multivariate meta_var_p50",
    national_contrast_successful_cities,
    round(national_p50_vs_p10_lag28_or, 6),
    paste0("P10 rolling mean=", round(low_light_threshold, 6)),
    "phase1_dynamic_zero",
    "none"
  ))
)
save_csv(diagnostics, "run_diagnostics.csv")

log_msg("All outputs written to: ", OUT_DIR)
cat(OUT_DIR, "\n")

