# Single-specification two-stage DLNM using recent low-sunlight-day count exposure.
#
# Main exposure:
#   lowdays3_t = sum(I(Sun_t < city-specific raw sunlight P20),
#                    I(Sun_{t-1} < city-specific P20),
#                    I(Sun_{t-2} < city-specific P20))
#
# For this transformed exposure X_t, fit a DLNM:
#   city-day high_events / n_observed ~
#     crossbasis(X_t, lag = 0:28) +
#     ns(same-day temperature, df = 3) +
#     same-day humidity + precipitation + wind speed +
#     ns(calendar time, df = 21) + day of week + pandemic phase
#
# This script keeps the main two-stage DLNM framework used in the absolute
# sunlight analysis. Only the exposure entering crossbasis() changes from
# same-day sunlight to the recent number of days below each city's city-specific
# raw sunlight P20 threshold. The main 0-28 day cumulative P50-vs-P90 OR is read
# from the national pooled cumulative exposure-response curve:
#   city DLNM -> crossreduce(type = "overall") -> meta_overall -> pooled curve.
# P50-vs-P90 lag-response summaries are produced with city
# crossreduce(type = "var") followed by one multivariate meta-analysis.

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

parse_logical_env <- function(name, default = FALSE) {
  raw <- tolower(trimws(Sys.getenv(name, unset = "")))
  if (!nzchar(raw)) return(default)
  if (raw %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (raw %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Invalid logical environment variable ", name, ": ", raw)
}

LOW_DAYS_WINDOWS <- as.integer(parse_numeric_env(
  "LOW_DAYS_WINDOWS",
  3L
))
LOW_DAYS_WINDOWS <- unique(LOW_DAYS_WINDOWS[
  is.finite(LOW_DAYS_WINDOWS) & LOW_DAYS_WINDOWS >= 1L
])
if (!length(LOW_DAYS_WINDOWS)) stop("No valid low-sunlight-day windows specified.")

LAG_SUN <- as.integer(parse_numeric_env("LAG_SUN", 28L)[1L])
SUN_VAR_KNOT_PROBS <- parse_numeric_env(
  "SUN_VAR_KNOT_PROBS",
  c(0.333, 0.667)
)
N_SUN_LAG_KNOTS <- as.integer(parse_numeric_env("N_SUN_LAG_KNOTS", 4L)[1L])
SUN_LAG_KNOTS_OVERRIDE <- parse_numeric_env("SUN_LAG_KNOTS_OVERRIDE", c(3, 7, 14, 21))
LOW_LIGHT_THRESHOLD_PROB <- parse_numeric_env("LOW_LIGHT_THRESHOLD_PROB", 0.20)[1L]
LOW_EXPOSURE_PROB <- parse_numeric_env("LOW_EXPOSURE_PROB", 0.90)[1L]
REFERENCE_PROB <- parse_numeric_env("REFERENCE_PROB", 0.50)[1L]
prob_tag <- function(p) paste0("p", formatC(round(100 * p), width = 2L, flag = "0"))
LOW_LIGHT_THRESHOLD_TAG <- prob_tag(LOW_LIGHT_THRESHOLD_PROB)
LOW_EXPOSURE_TAG <- prob_tag(LOW_EXPOSURE_PROB)
REFERENCE_TAG <- prob_tag(REFERENCE_PROB)

TIME_DF_PER_YEAR <- 7L
TIME_DF <- TIME_DF_PER_YEAR * 3L
TEMP_NS_DF <- as.integer(parse_numeric_env("TEMP_NS_DF", 3L)[1L])

MIN_MODEL_DAYS <- 365L
MIN_TOTAL_EVENTS <- 20L
MIN_TOTAL_OBSERVED <- 100L
TRIM_SUNLIGHT_ABOVE_FULL_P99 <- parse_logical_env(
  "TRIM_SUNLIGHT_ABOVE_FULL_P99",
  FALSE
)

MIXMETA_MAXITER <- as.integer(parse_numeric_env("MIXMETA_MAXITER", 1000L)[1L])
MIXMETA_CONTROL <- mixmeta.control(
  maxiter = MIXMETA_MAXITER,
  optim = list(maxit = MIXMETA_MAXITER)
)

LAG_LABEL <- paste0("lag0_", LAG_SUN)
LAG_TITLE <- paste0("Lag 0-", LAG_SUN)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR <- file.path(OUTPUT_ROOT, paste0("low_sunlight_days_", timestamp))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "stage1_status"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "stage2_summaries"), recursive = TRUE, showWarnings = FALSE)
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
    "1_dynamic_zero_covid",
    fifelse(
      date <= as.IDate("2022-12-06"),
      "2_twenty_measures_transition",
      fifelse(
        date <= as.IDate("2023-01-22"),
        "3_rapid_transition_first_wave",
        "4_normalized_period"
      )
    )
  )
  factor(
    phase,
    levels = c(
      "1_dynamic_zero_covid",
      "2_twenty_measures_transition",
      "3_rapid_transition_first_wave",
      "4_normalized_period"
    )
  )
}

make_significance <- function(low, high) {
  fifelse(
    is.finite(low) & is.finite(high) & low > 1,
    "significantly_higher",
    fifelse(
      is.finite(low) & is.finite(high) & high < 1,
      "significantly_lower",
      "not_significant"
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

first_prediction_path <- function(x, n_expected) {
  if (is.null(x)) return(rep(NA_real_, n_expected))
  if (is.matrix(x) || is.data.frame(x)) {
    m <- as.matrix(x)
    if (nrow(m) == 1L && ncol(m) >= n_expected) {
      return(as.numeric(m[1L, seq_len(n_expected)]))
    }
    if (ncol(m) == 1L && nrow(m) >= n_expected) {
      return(as.numeric(m[seq_len(n_expected), 1L]))
    }
    if (nrow(m) >= n_expected) {
      return(as.numeric(m[seq_len(n_expected), 1L]))
    }
    if (ncol(m) >= n_expected) {
      return(as.numeric(m[1L, seq_len(n_expected)]))
    }
  }
  v <- as.numeric(x)
  if (length(v) < n_expected) {
    return(c(v, rep(NA_real_, n_expected - length(v))))
  }
  v[seq_len(n_expected)]
}

is_positive_definite <- function(mat, tol = 1e-8) {
  if (is.null(mat) || !is.matrix(mat) || any(!is.finite(mat))) {
    return(FALSE)
  }
  mat <- (mat + t(mat)) / 2
  out <- tryCatch(
    chol(mat + diag(tol, nrow(mat))),
    error = function(e) NULL
  )
  !is.null(out)
}

nearest_curve_rows <- function(curve, values, labels) {
  rbindlist(lapply(seq_along(values), function(i) {
    idx <- which.min(abs(curve$exposure_value - values[i]))
    out <- copy(curve[idx])
    out[, `:=`(
      requested_label = labels[i],
      requested_value = values[i],
      matched_value = exposure_value
    )]
    out
  }), fill = TRUE)
}

predict_reduced_curve <- function(
    model,
    at_values,
    cen,
    exposure_argvar,
    curve_type,
    spec_id = NA_character_,
    window_days = NA_integer_) {
  basis_at <- do.call(onebasis, c(list(x = at_values), exposure_argvar))
  basis_cen <- do.call(onebasis, c(list(x = cen), exposure_argvar))
  contrast <- sweep(basis_at, 2L, as.numeric(basis_cen), FUN = "-")

  beta <- coef(model)
  V <- vcov(model)
  if (ncol(contrast) != length(beta)) {
    stop(
      "Exposure basis dimension mismatch: contrast has ", ncol(contrast),
      " columns but meta model has ", length(beta), " coefficients."
    )
  }

  log_or <- as.numeric(contrast %*% beta)
  var_log_or <- rowSums((contrast %*% V) * contrast)
  se_log_or <- sqrt(pmax(var_log_or, 0))

  out <- data.table(
    spec_id = spec_id,
    window_days = window_days,
    curve_type = curve_type,
    exposure_value = at_values,
    reference_value = cen,
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
    exposure_label,
    reference_label,
    exposure_value,
    reference_value,
    cumulative = FALSE,
    spec_id = NA_character_,
    window_days = NA_integer_) {
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
    spec_id = spec_id,
    window_days = window_days,
    curve_type = curve_type,
    lag_day = as.integer(lag_days),
    exposure_label = exposure_label,
    reference_label = reference_label,
    exposure_value = exposure_value,
    reference_value = reference_value,
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

make_exposure_argvar <- function(x, knot_probs = SUN_VAR_KNOT_PROBS) {
  x <- x[is.finite(x)]
  if (length(unique(x)) < 3L) {
    return(list(
      argvar = list(fun = "lin", intercept = FALSE),
      knots = numeric(0),
      boundary = range(x, na.rm = TRUE),
      basis_type = "linear"
    ))
  }

  boundary <- range(x, na.rm = TRUE)
  knots <- as.numeric(quantile(x, knot_probs, na.rm = TRUE, type = 7))
  knots <- sort(unique(knots))
  knots <- knots[is.finite(knots) & knots > boundary[1L] & knots < boundary[2L]]

  if (!length(knots)) {
    return(list(
      argvar = list(fun = "lin", intercept = FALSE),
      knots = numeric(0),
      boundary = boundary,
      basis_type = "linear"
    ))
  }

  list(
    argvar = list(
      fun = "ns",
      knots = knots,
      Boundary.knots = boundary,
      intercept = FALSE
    ),
    knots = knots,
    boundary = boundary,
    basis_type = "natural_spline"
  )
}

make_crossbasis <- function(x, exposure_argvar, lag_knots) {
  crossbasis(
    x,
    lag = LAG_SUN,
    argvar = exposure_argvar,
    arglag = list(
      fun = "ns",
      knots = lag_knots,
      intercept = TRUE
    )
  )
}

make_prediction_template <- function(grid, exposure_argvar, lag_knots) {
  make_crossbasis(
    grid,
    exposure_argvar = exposure_argvar,
    lag_knots = lag_knots
  )
}

# -----------------------------------------------------------------------------
# 3. Read data and build complete panel
# -----------------------------------------------------------------------------

log_msg("Output directory: ", OUT_DIR)
log_msg("City-day outcome: ", CITY_DAY_PATH)
log_msg("Weather controls: ", WEATHER_PATH)
log_msg("Sunlight source: ", SUNLIGHT_PATH, " [daily_sun_hours]")
log_msg("Low-sunlight-day windows: ", paste(LOW_DAYS_WINDOWS, collapse = ", "))
log_msg("Sunlight DLNM lag: 0-", LAG_SUN)
log_msg("Temperature: same-day natural spline control, df=", TEMP_NS_DF)
log_msg("Main stage 2: crossreduce(type='overall') -> meta_overall")

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
log_msg("Sunlight path used: ", SUNLIGHT_PATH)

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

city_low_sunlight_thresholds <- panel[
  is.finite(sunlight_duration),
  .(
    city_low_sunlight_p20 = as.numeric(quantile(
      sunlight_duration,
      LOW_LIGHT_THRESHOLD_PROB,
      na.rm = TRUE,
      type = 7
    )),
    finite_sunlight_days = .N
  ),
  by = city_code
]
panel <- merge(
  panel,
  city_low_sunlight_thresholds,
  by = "city_code",
  all.x = TRUE,
  sort = FALSE
)
setorder(panel, city_code, date)

panel[, low_sun_raw_city_p20 := fifelse(
  is.finite(sunlight_duration) & is.finite(city_low_sunlight_p20),
  as.numeric(sunlight_duration < city_low_sunlight_p20),
  NA_real_
)]

city_low_counts <- panel[
  !is.na(low_sun_raw_city_p20),
  .(
    rows_classified_low_sunlight = sum(low_sun_raw_city_p20 == 1),
    rows_classified_non_low_sunlight = sum(low_sun_raw_city_p20 == 0)
  ),
  by = city_code
]
city_low_sunlight_thresholds <- merge(
  city_low_sunlight_thresholds,
  city_low_counts,
  by = "city_code",
  all.x = TRUE,
  sort = FALSE
)
city_low_sunlight_thresholds[
  is.na(rows_classified_low_sunlight),
  rows_classified_low_sunlight := 0L
]
city_low_sunlight_thresholds[
  is.na(rows_classified_non_low_sunlight),
  rows_classified_non_low_sunlight := 0L
]

save_csv(
  city_low_sunlight_thresholds[
    ,
    `:=`(
      threshold_source = "city-specific raw daily sunlight_duration",
      threshold_probability = LOW_LIGHT_THRESHOLD_PROB,
      comparison_operator = "sunlight_duration < city-specific threshold"
    )
  ],
  "low_sunlight_city_p20_threshold_summary.csv"
)

save_csv(
  data.table(
    threshold_source = "city-specific raw daily sunlight_duration",
    threshold_probability = LOW_LIGHT_THRESHOLD_PROB,
    median_city_low_sunlight_threshold_hours = median(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE),
    mean_city_low_sunlight_threshold_hours = mean(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE),
    min_city_low_sunlight_threshold_hours = min(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE),
    max_city_low_sunlight_threshold_hours = max(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE),
    comparison_operator = "sunlight_duration < city-specific threshold",
    rows_with_finite_sunlight = panel[is.finite(sunlight_duration), .N],
    rows_classified_low_sunlight = panel[low_sun_raw_city_p20 == 1, .N],
    rows_classified_non_low_sunlight = panel[low_sun_raw_city_p20 == 0, .N]
  ),
  "low_sunlight_city_p20_threshold_overall_summary.csv"
)
log_msg(
  "Low sunlight threshold: city-specific raw daily sunlight ",
  toupper(LOW_LIGHT_THRESHOLD_TAG),
  "; median threshold = ",
  round(median(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE), 6),
  " h/day."
)

for (w in LOW_DAYS_WINDOWS) {
  col <- paste0("low_sun_days_w", w)
  panel[, (col) := frollsum(
    low_sun_raw_city_p20,
    n = w,
    align = "right",
    na.rm = FALSE
  ), by = city_code]
}

panel_export_cols <- c(
  "city_code", "city_name", "date", "high_events", "n_observed",
  "high_stress_rate", "sunlight_duration", "city_low_sunlight_p20",
  "low_sun_raw_city_p20",
  paste0("low_sun_days_w", LOW_DAYS_WINDOWS),
  "temperature", "humidity", "precipitation", "wind_speed",
  "time_num", "dow", "pandemic_phase"
)
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

sun_lag_knots <- if (length(SUN_LAG_KNOTS_OVERRIDE)) {
  SUN_LAG_KNOTS_OVERRIDE
} else {
  as.numeric(dlnm::logknots(LAG_SUN, N_SUN_LAG_KNOTS))
}
log_msg("Sunlight lag knots: ", paste(round(sun_lag_knots, 6), collapse = ", "))

# -----------------------------------------------------------------------------
# 4. Fit one low-sunlight-day count DLNM specification
# -----------------------------------------------------------------------------

fit_one_window <- function(window_days) {
  spec_id <- paste0("lowdays", window_days)
  exposure_col <- paste0("low_sun_days_w", window_days)
  log_msg("Running low-sunlight-day count DLNM spec with meta_overall: ", spec_id)

  exposure_global <- panel[is.finite(get(exposure_col)), get(exposure_col)]
  if (length(unique(exposure_global)) < 3L) {
    stop("Too few unique exposure values for ", spec_id)
  }

  exposure_basis <- make_exposure_argvar(exposure_global)
  exposure_boundary <- exposure_basis$boundary
  exposure_knots <- exposure_basis$knots
  exposure_argvar <- exposure_basis$argvar
  exposure_basis_type <- exposure_basis$basis_type

  exposure_probs <- c(
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
  exposure_pct <- setNames(
    as.numeric(quantile(exposure_global, exposure_probs, na.rm = TRUE, type = 7)),
    names(exposure_probs)
  )
  high_value <- as.numeric(quantile(exposure_global, LOW_EXPOSURE_PROB, na.rm = TRUE, type = 7))
  ref_value <- as.numeric(quantile(exposure_global, REFERENCE_PROB, na.rm = TRUE, type = 7))
  if (!is.finite(high_value) || !is.finite(ref_value) || high_value <= ref_value) {
    stop("Invalid national P90/P50 low-sunlight-day contrast for ", spec_id)
  }

  exposure_knot_parameters <- if (length(exposure_knots)) {
    data.table(
      spec_id = spec_id,
      parameter = paste0("exposure_knot_", seq_along(exposure_knots)),
      value = exposure_knots
    )
  } else {
    data.table(
      spec_id = spec_id,
      parameter = "exposure_internal_knots",
      value = "none"
    )
  }

  spec_parameters <- rbindlist(list(
    exposure_knot_parameters,
    data.table(
      spec_id = spec_id,
      parameter = paste0("sunlight_lag_knot_", seq_along(sun_lag_knots)),
      value = sun_lag_knots
    ),
    data.table(
      spec_id = spec_id,
      parameter = c(
        "low_sunlight_threshold_probability",
        "median_city_low_sunlight_threshold_hours",
        "low_days_window_days",
        "low_days_window_lag_start",
        "low_days_window_lag_end",
        "dlnm_lag_start",
        "dlnm_lag_end",
        "raw_sunlight_farthest_lag_used",
        "exposure_basis_type",
        "exposure_boundary_lower",
        "exposure_boundary_upper",
        "global_exposure_P90_comparison",
        "global_exposure_p50_reference",
        "temperature_ns_df",
        "time_df",
        "main_cumulative_curve_source",
        "P90_P50_contrast_source",
        "P90_P50_lag_response_source"
      ),
      value = c(
        LOW_LIGHT_THRESHOLD_PROB,
        median(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE),
        window_days,
        0,
        window_days - 1L,
        0,
        LAG_SUN,
        LAG_SUN + window_days - 1L,
        exposure_basis_type,
        exposure_boundary[1L],
        exposure_boundary[2L],
        high_value,
        ref_value,
        TEMP_NS_DF,
        TIME_DF,
        "city crossreduce(type='overall') -> meta_overall -> national pooled cumulative curve",
        "P50-vs-P90 lag0-28 cumulative OR is read from the national pooled cumulative curve",
        "city crossreduce(type='var', value=national P90, cen=national P50) -> multivariate mixmeta"
      )
    )
  ), fill = TRUE)

  save_csv(
    spec_parameters,
    file.path("stage2_summaries", paste0("parameters_", spec_id, ".csv"))
  )
  save_csv(
    data.table(
      spec_id = spec_id,
      percentile = names(exposure_pct),
      probability = as.numeric(exposure_probs),
      exposure_value = as.numeric(exposure_pct)
    ),
    file.path("stage2_summaries", paste0("percentiles_", spec_id, ".csv"))
  )

  run_city_model <- function(city_id) {
    dt <- copy(panel[city_code == city_id])
    setorder(dt, date)
    city_label <- dt[!is.na(city_name) & nzchar(city_name), city_name[1L]]
    if (!length(city_label)) city_label <- NA_character_

    base_status <- data.table(
      spec_id = spec_id,
      window_days = window_days,
      city_code = city_id,
      city_name = city_label,
      panel_days = nrow(dt),
      complete_days_after_lag = NA_integer_,
      total_observed_after_lag = NA_real_,
      total_high_events_after_lag = NA_real_,
      fit_status = "not_started",
      message = NA_character_,
      dispersion = NA_real_,
      n_coefficients = NA_integer_,
      city_exposure_P90 = NA_real_,
      city_exposure_p50 = NA_real_,
      national_exposure_P90 = high_value,
      national_exposure_p50 = ref_value,
      overall_reduction_status = NA_character_,
      overall_reduction_message = NA_character_,
      p50_var_lag_response_status = NA_character_,
      p50_var_lag_response_message = NA_character_
    )

    tryCatch({
      cb_exposure <- make_crossbasis(
        dt[[exposure_col]],
        exposure_argvar = exposure_argvar,
        lag_knots = sun_lag_knots
      )

      complete_after_lag <- complete.cases(cbind(
        dt[, .(
          high_stress_rate, n_observed, high_events, failures,
          temperature, humidity, precipitation, wind_speed,
          time_num, dow, pandemic_phase
        )],
        cb_exposure
      )) &
        dt$n_observed > 0 & dt$high_events >= 0 & dt$failures >= 0

      city_exposure_complete <- dt[[exposure_col]][complete_after_lag]
      city_exposure_complete <- city_exposure_complete[is.finite(city_exposure_complete)]
      city_P90 <- if (length(city_exposure_complete)) {
        as.numeric(quantile(city_exposure_complete, LOW_EXPOSURE_PROB, na.rm = TRUE, type = 7))
      } else {
        NA_real_
      }
      city_p50 <- if (length(city_exposure_complete)) {
        as.numeric(quantile(city_exposure_complete, REFERENCE_PROB, na.rm = TRUE, type = 7))
      } else {
        NA_real_
      }

      n_complete <- sum(complete_after_lag)
      total_observed <- sum(dt$n_observed[complete_after_lag], na.rm = TRUE)
      total_events <- sum(dt$high_events[complete_after_lag], na.rm = TRUE)
      base_status[, `:=`(
        complete_days_after_lag = n_complete,
        total_observed_after_lag = total_observed,
        total_high_events_after_lag = total_events,
        city_exposure_P90 = city_P90,
        city_exposure_p50 = city_p50
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
          cb_exposure +
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

      full_coef_names <- grep("^cb_exposure", names(coef(fit)), value = TRUE)
      full_coef <- coef(fit)[full_coef_names]
      full_vcov <- vcov(fit)[full_coef_names, full_coef_names, drop = FALSE]
      if (any(!is.finite(full_coef)) || any(!is.finite(full_vcov))) {
        stop("non-finite full cross-basis coefficient or covariance")
      }

      reduced_overall <- crossreduce(
        cb_exposure,
        fit,
        type = "overall",
        cen = high_value
      )
      overall_coef <- coef(reduced_overall)
      overall_vcov <- vcov(reduced_overall)
      if (any(!is.finite(overall_coef)) || any(!is.finite(overall_vcov))) {
        stop("non-finite overall reduced coefficient or covariance")
      }

      reduced_var_p50 <- crossreduce(
        cb_exposure,
        fit,
        type = "var",
        value = ref_value,
        cen = high_value
      )
      var_p50_coef <- coef(reduced_var_p50)
      var_p50_vcov <- vcov(reduced_var_p50)
      if (any(!is.finite(var_p50_coef)) || any(!is.finite(var_p50_vcov))) {
        stop("non-finite comparison P50-vs-P90 lag-response reduced coefficient or covariance")
      }

      base_status[, `:=`(
        fit_status = "ok",
        message = NA_character_,
        dispersion = summary(fit)$dispersion,
        n_coefficients = length(coef(fit)),
        overall_reduction_status = "ok",
        overall_reduction_message = NA_character_,
        p50_var_lag_response_status = "ok",
        p50_var_lag_response_message = NA_character_
      )]

      list(
        status = base_status,
        full_coef = full_coef,
        full_vcov = full_vcov,
        overall_coef = overall_coef,
        overall_vcov = overall_vcov,
        var_p50_coef = var_p50_coef,
        var_p50_vcov = var_p50_vcov
      )
    }, error = function(e) {
      base_status[, `:=`(
        fit_status = "error",
        message = conditionMessage(e)
      )]
      list(status = base_status)
    })
  }

  stage1_results <- lapply(model_cities, run_city_model)
  names(stage1_results) <- model_cities

  stage1_status <- rbindlist(lapply(stage1_results, `[[`, "status"), fill = TRUE)
  save_csv(
    stage1_status,
    file.path("stage1_status", paste0("stage1_city_model_status_", spec_id, ".csv"))
  )

  ok_cities <- stage1_status[fit_status == "ok", city_code]
  log_msg(
    spec_id, " successful cities: ", length(ok_cities), " / ", length(model_cities)
  )

  city_percentiles <- stage1_status[, .(
    spec_id,
    window_days,
    low_days_window_label = paste0("lowdays", window_days),
    low_days_window_lag = paste0("lag0-", window_days - 1L),
    city_code,
    city_name,
    complete_days_after_lag,
    total_observed_after_lag,
    total_high_events_after_lag,
    city_exposure_P90,
    city_exposure_p50,
    national_exposure_P90,
    national_exposure_p50,
    overall_reduction_status,
    p50_var_lag_response_status,
    p50_var_lag_response_message
  )]
  save_csv(
    city_percentiles,
    file.path("stage2_summaries", paste0("stage1_city_and_national_P90_P50_", spec_id, ".csv"))
  )

  if (length(ok_cities) < 5L) {
    warning("Too few successful cities for ", spec_id)
    return(list(
      spec_id = spec_id,
      status = stage1_status,
      ok_cities = ok_cities,
      meta_status = data.table(
        spec_id = spec_id,
        window_days = window_days,
        stage1_successful_cities = length(ok_cities),
        meta_eligible_cities = length(ok_cities),
        converged = FALSE,
        message = "too few successful cities"
      )
    ))
  }

  full_coef_mat <- do.call(rbind, lapply(ok_cities, function(city_id) {
    stage1_results[[city_id]]$full_coef
  }))
  rownames(full_coef_mat) <- ok_cities
  full_vcov_list <- lapply(ok_cities, function(city_id) {
    stage1_results[[city_id]]$full_vcov
  })
  names(full_vcov_list) <- ok_cities

  coef_long <- rbindlist(lapply(ok_cities, function(city_id) {
    data.table(
      spec_id = spec_id,
      window_days = window_days,
      city_code = city_id,
      coefficient = names(stage1_results[[city_id]]$full_coef),
      beta = as.numeric(stage1_results[[city_id]]$full_coef)
    )
  }))
  save_csv(
    coef_long,
    file.path("stage2_summaries", paste0("stage1_full_crossbasis_coefficients_", spec_id, ".csv"))
  )

  overall_coef_mat <- do.call(rbind, lapply(ok_cities, function(city_id) {
    stage1_results[[city_id]]$overall_coef
  }))
  rownames(overall_coef_mat) <- ok_cities
  overall_vcov_list <- lapply(ok_cities, function(city_id) {
    stage1_results[[city_id]]$overall_vcov
  })
  names(overall_vcov_list) <- ok_cities
  save_csv(
    as.data.table(overall_coef_mat, keep.rownames = "city_code"),
    file.path("stage2_summaries", paste0("stage1_overall_reduced_coefficients_", spec_id, ".csv"))
  )

  log_msg("Fitting ", spec_id, " cumulative meta_overall...")
  meta_overall <- mixmeta(
    overall_coef_mat ~ 1,
    S = overall_vcov_list,
    method = "reml",
    control = MIXMETA_CONTROL
  )
  capture.output(
    summary(meta_overall),
    file = file.path(
      OUT_DIR,
      "stage2_summaries",
      paste0("stage2_overall_reduced_meta_summary_", spec_id, ".txt")
    )
  )

  meta_status <- data.table(
    spec_id = spec_id,
    window_days = window_days,
    stage1_successful_cities = length(ok_cities),
    meta_eligible_cities = length(ok_cities),
    meta_source = "crossreduce(type='overall') -> meta_overall",
    converged = isTRUE(meta_overall$converged),
    method = meta_overall$method,
    logLik = as.numeric(logLik(meta_overall)),
    AIC = AIC(meta_overall),
    BIC = BIC(meta_overall),
    message = NA_character_
  )

  grid_range <- range(exposure_global, na.rm = TRUE)
  at_grid <- seq(grid_range[1L], grid_range[2L], length.out = 250L)
  fixed_values <- seq(
    ceiling(exposure_boundary[1L]),
    floor(exposure_boundary[2L]),
    by = 1
  )
  fixed_values <- fixed_values[
    fixed_values >= exposure_boundary[1L] &
      fixed_values <= exposure_boundary[2L]
  ]
  selected_pred_values <- sort(unique(c(
    fixed_values,
    as.numeric(exposure_pct)
  )))

  cb_template <- make_prediction_template(
    at_grid,
    exposure_argvar = exposure_argvar,
    lag_knots = sun_lag_knots
  )
  cumulative_curve <- predict_reduced_curve(
    meta_overall,
    at_grid,
    high_value,
    exposure_argvar,
    curve_type = paste0("main_cumulative_overall_reduced_", spec_id, "_", LAG_LABEL),
    spec_id = spec_id,
    window_days = window_days
  )
  cumulative_curve[, `:=`(
    exposure_p90 = high_value,
    low_days_window_label = paste0("lowdays", window_days),
    dlnm_lag_label = LAG_LABEL
  )]

  selected_dt <- rbindlist(list(
    data.table(
      value_type = "percentile",
      requested_label = c("P1", "P5", "P10", "P25", "P50", "P75", "P90", "P95", "P99"),
      requested_value = as.numeric(exposure_pct[c(
        "p01", "p05", "p10", "p25", "p50", "p75", "p90", "p95", "p99"
      )])
    ),
    data.table(
      value_type = "fixed_low_sun_days",
      requested_label = paste0(fixed_values, " days"),
      requested_value = fixed_values
    )
  ), fill = TRUE)
  selected_dt <- unique(selected_dt, by = c("value_type", "requested_label", "requested_value"))
  setorder(selected_dt, requested_value, value_type, requested_label)
  selected_pred_values <- sort(unique(selected_dt$requested_value[is.finite(selected_dt$requested_value)]))

  selected_prediction <- predict_reduced_curve(
    meta_overall,
    selected_pred_values,
    high_value,
    exposure_argvar,
    curve_type = paste0("selected_points_from_main_national_cumulative_curve_", spec_id, "_", LAG_LABEL),
    spec_id = spec_id,
    window_days = window_days
  )
  setnames(selected_prediction, "exposure_value", "requested_value")
  selected_cumulative <- merge(
    selected_dt,
    selected_prediction[, .(
      requested_value, curve_type, log_OR, se_log_OR,
      OR, OR_low, OR_high, p_value, significance
    )],
    by = "requested_value",
    all.x = TRUE,
    sort = FALSE,
    allow.cartesian = TRUE
  )
  selected_cumulative[, `:=`(
    spec_id = spec_id,
    window_days = window_days,
    reference_value = high_value
  )]
  setcolorder(
    selected_cumulative,
    c(
      "spec_id", "window_days", "value_type", "requested_label",
      "requested_value", "reference_value", "log_OR", "se_log_OR",
      "OR", "OR_low", "OR_high", "p_value", "significance", "curve_type"
    )
  )
  selected_cumulative[, comparison := sprintf(
    "%s %.6f low-sunlight days versus national P90 %.6f low-sunlight days",
    requested_label,
    requested_value,
    high_value
  )]

  p50_from_national_cumulative_curve <- copy(
    selected_cumulative[value_type == "percentile" & requested_label == "P50"]
  )
  p50_from_national_cumulative_curve[, `:=`(
    curve_type = paste0("P50_vs_P90_from_main_national_cumulative_curve_", spec_id, "_", LAG_LABEL),
    lag_day = LAG_SUN,
    exposure_label = "national_P50",
    reference_label = "national_P90",
    exposure_value = ref_value,
    reference_value = high_value,
    contrast_source = "main national cumulative exposure-response curve: crossreduce(type='overall') -> meta_overall"
  )]

  save_csv(
    p50_from_national_cumulative_curve,
    file.path("stage2_summaries", paste0("meta_P50_vs_P90_from_national_cumulative_curve_", spec_id, "_", LAG_LABEL, "_OR.csv"))
  )

  p50_var_lag_response_status <- stage1_status[, .(
    spec_id,
    window_days,
    city_code,
    city_name,
    fit_status,
    p50_var_lag_response_status,
    p50_var_lag_response_message,
    complete_days_after_lag,
    national_exposure_P90,
    national_exposure_p50
  )]
  save_csv(
    p50_var_lag_response_status,
    file.path("stage2_summaries", paste0("stage1_national_P50_P90_var_reduction_status_", spec_id, ".csv"))
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
    stop("Too few valid crossreduce(type='var') P50-vs-P90 city reductions for ", spec_id)
  }

  var_p50_coef_mat <- do.call(rbind, var_p50_coef_list[var_p50_valid])
  rownames(var_p50_coef_mat) <- ok_cities[var_p50_valid]
  save_csv(
    as.data.table(var_p50_coef_mat, keep.rownames = "city_code"),
    file.path("stage2_summaries", paste0("stage1_national_P50_vs_P90_var_reduced_coefficients_", spec_id, ".csv"))
  )

  log_msg("Fitting ", spec_id, " national P50 versus P90 comparison lag-response meta_var...")
  meta_var_p50 <- mixmeta(
    var_p50_coef_mat ~ 1,
    S = var_p50_vcov_list[var_p50_valid],
    method = "reml",
    control = MIXMETA_CONTROL
  )
  capture.output(
    summary(meta_var_p50),
    file = file.path(OUT_DIR, "stage2_summaries", paste0("stage2_var_P50_vs_P90_meta_summary_", spec_id, ".txt"))
  )

  daily_lag <- predict_var_lag_response(
    meta_var_p50,
    lag_days = 0:LAG_SUN,
    lag_knots = sun_lag_knots,
    max_lag = LAG_SUN,
    curve_type = paste0("meta_var_national_P50_vs_P90_daily_lag_response_", spec_id),
    exposure_label = "national_P50",
    reference_label = "national_P90",
    exposure_value = ref_value,
    reference_value = high_value,
    cumulative = FALSE,
    spec_id = spec_id,
    window_days = window_days
  )
  daily_lag[, p_value_BH := p.adjust(p_value, method = "BH")]
  daily_lag[, significant_BH := is.finite(p_value_BH) & p_value_BH < 0.05]

  cumulative_lag <- predict_var_lag_response(
    meta_var_p50,
    lag_days = 0:LAG_SUN,
    lag_knots = sun_lag_knots,
    max_lag = LAG_SUN,
    curve_type = paste0("meta_var_national_P50_vs_P90_cumulative_", spec_id, "_", LAG_LABEL),
    exposure_label = "national_P50",
    reference_label = "national_P90",
    exposure_value = ref_value,
    reference_value = high_value,
    cumulative = TRUE,
    spec_id = spec_id,
    window_days = window_days
  )
  cumulative_lag[, p_value_BH := p.adjust(p_value, method = "BH")]
  cumulative_lag[, significant_BH := is.finite(p_value_BH) & p_value_BH < 0.05]

  save_csv(
    cumulative_curve,
    file.path("stage2_summaries", paste0("meta_cumulative_exposure_response_curve_", spec_id, ".csv"))
  )
  save_csv(
    selected_cumulative,
    file.path("stage2_summaries", paste0("meta_cumulative_selected_ORs_", spec_id, ".csv"))
  )
  save_csv(
    daily_lag,
    file.path("stage2_summaries", paste0("meta_national_P50_vs_P90_daily_lag_OR_", spec_id, ".csv"))
  )
  save_csv(
    cumulative_lag,
    file.path("stage2_summaries", paste0("meta_national_P50_vs_P90_incremental_cumulative_lag_OR_", spec_id, ".csv"))
  )

  # Figure generation omitted.

  full_high <- p50_from_national_cumulative_curve[1L]
  onset <- cumulative_lag[
    significant_BH == TRUE & OR_high < 1,
    .SD[1L]
  ]
  if (!nrow(onset)) {
    onset <- data.table(
      lag_day = NA_integer_,
      OR = NA_real_,
      OR_low = NA_real_,
      OR_high = NA_real_,
      p_value = NA_real_,
      p_value_BH = NA_real_
    )
  }

  summary <- data.table(
    spec_id = spec_id,
    window_days = window_days,
    low_days_window_lag = paste0("lag0-", window_days - 1L),
    dlnm_lag = LAG_LABEL,
    raw_sunlight_farthest_lag_used = LAG_SUN + window_days - 1L,
    stage1_successful_cities = length(ok_cities),
    meta_eligible_cities = length(ok_cities),
    p50_var_lag_response_cities = sum(var_p50_valid),
    low_sunlight_definition = paste0("raw daily_sun_hours < city-specific ", toupper(LOW_LIGHT_THRESHOLD_TAG)),
    low_sunlight_threshold_probability = LOW_LIGHT_THRESHOLD_PROB,
    median_city_low_sunlight_threshold_hours = median(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE),
    exposure_p90 = high_value,
    exposure_p50_reference = ref_value,
    P50_vs_p90_lag0_28_cumulative_OR_source = "national pooled cumulative exposure-response curve from meta_overall",
    p50_vs_p90_lag0_28_cumulative_OR = full_high$OR,
    p50_vs_p90_lag0_28_cumulative_OR_low = full_high$OR_low,
    p50_vs_p90_lag0_28_cumulative_OR_high = full_high$OR_high,
    p50_vs_p90_lag0_28_cumulative_p_value = full_high$p_value,
    p50_vs_p90_lag0_28_cumulative_significance = full_high$significance,
    first_BH_significant_higher_cumulative_lag = onset$lag_day,
    first_BH_significant_higher_cumulative_OR = onset$OR,
    first_BH_significant_higher_cumulative_OR_low = onset$OR_low,
    first_BH_significant_higher_cumulative_OR_high = onset$OR_high
  )

  list(
    spec_id = spec_id,
    window_days = window_days,
    status = stage1_status,
    meta_status = meta_status,
    summary = summary,
    cumulative_curve = cumulative_curve,
    selected_cumulative = selected_cumulative,
    daily_lag = daily_lag,
    cumulative_lag = cumulative_lag,
    meta_overall = meta_overall,
    meta_var_p50 = meta_var_p50
  )
}

# -----------------------------------------------------------------------------
# 5. Run all requested specifications
# -----------------------------------------------------------------------------

results <- lapply(LOW_DAYS_WINDOWS, fit_one_window)
names(results) <- paste0("lowdays", LOW_DAYS_WINDOWS)

stage1_all <- rbindlist(lapply(results, `[[`, "status"), fill = TRUE)
meta_status_all <- rbindlist(lapply(results, `[[`, "meta_status"), fill = TRUE)
summary_all <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
curve_all <- rbindlist(lapply(results, `[[`, "cumulative_curve"), fill = TRUE)
selected_all <- rbindlist(lapply(results, `[[`, "selected_cumulative"), fill = TRUE)
daily_all <- rbindlist(lapply(results, `[[`, "daily_lag"), fill = TRUE)
cumlag_all <- rbindlist(lapply(results, `[[`, "cumulative_lag"), fill = TRUE)

save_csv(stage1_all, "stage1_city_model_status_all_windows.csv")
save_csv(meta_status_all, "stage2_overall_reduced_meta_status_all_windows.csv")
save_csv(summary_all, "summary_national_P50_vs_P90_low_sun_days_windows.csv")
save_csv(curve_all, "meta_cumulative_exposure_response_curves_all_windows.csv")
save_csv(selected_all, "meta_cumulative_selected_ORs_all_windows.csv")
save_csv(daily_all, "meta_national_P50_vs_P90_daily_lag_OR_all_windows.csv")
save_csv(cumlag_all, "meta_national_P50_vs_P90_incremental_cumulative_lag_OR_all_windows.csv")

model_specification <- data.table(
  item = c(
    "analysis_unit",
    "outcome",
    "raw_sunlight_source",
    "raw_sunlight_variable",
    "exposure_main",
    "low_sunlight_definition",
    "low_days_window_alignment",
    "sunlight_dlnm",
    "exposure_spline",
    "lag_spline",
    "comparison",
    "temperature_control",
    "other_weather_controls",
    "time_control",
    "weekday_control",
    "pandemic_phase_control",
    "stage1_family",
    "stage1_weights",
    "stage2_model",
    "minimum_city_days",
    "minimum_city_events",
    "minimum_city_observed"
  ),
  value = c(
    "city-date",
    "city-day high-stress event rate",
    SUNLIGHT_PATH,
    "daily_sun_hours renamed to sunlight_duration before defining low-sunlight days",
    paste0(
      "lowdays windows: ",
      paste(LOW_DAYS_WINDOWS, collapse = ", "),
      "; each exposure is the count of low-sunlight days over current day plus previous K-1 days"
    ),
    paste0(
      "low sunlight is raw daily_sun_hours < city-specific ",
      toupper(LOW_LIGHT_THRESHOLD_TAG),
      "; median city threshold = ",
      round(median(city_low_sunlight_thresholds$city_low_sunlight_p20, na.rm = TRUE), 6), " h/day"
    ),
    "right-aligned frollsum: X_t = sum(LowSun_t ... LowSun_{t-k+1})",
    paste0("crossbasis of recent low-sunlight-day count, lag 0-", LAG_SUN),
    paste0(
      "natural spline with internal knots at global low-sunlight-day count quantiles ",
      paste(SUN_VAR_KNOT_PROBS, collapse = ", ")
    ),
    paste0("natural spline with lag knots at ", paste(round(sun_lag_knots, 6), collapse = ", "), " days"),
    "national P90 low-sunlight-day count versus national P50 reference, read from the pooled meta_overall cumulative exposure-response curve",
    paste0("same-day natural spline temperature, df=", TEMP_NS_DF),
    "same-day linear humidity, precipitation, and wind speed",
    paste0("natural spline of time, df=", TIME_DF),
    "day-of-week factor, Sunday reference",
    "four-level pandemic phase factor, dynamic zero-COVID reference",
    "quasibinomial GLM",
    "n_observed",
    "main: crossreduce(type='overall') city coefficients pooled by intercept-only multivariate random-effects meta-analysis, REML; lag-response: crossreduce(type='var') pooled by multivariate meta",
    MIN_MODEL_DAYS,
    MIN_TOTAL_EVENTS,
    MIN_TOTAL_OBSERVED
  )
)
save_csv(model_specification, "model_specification.csv")

# Figure generation omitted.

saveRDS(
  list(
    results = results,
    stage1_all = stage1_all,
    meta_status_all = meta_status_all,
    summary_all = summary_all,
    curve_all = curve_all,
    selected_all = selected_all,
    daily_all = daily_all,
    cumlag_all = cumlag_all,
    model_specification = model_specification
  ),
  file.path(OUT_DIR, "low_sun_days_dlnm_results.rds")
)

log_msg("Finished all low-sunlight-day count DLNM specifications.")
log_msg("Output directory: ", OUT_DIR)
message("Output directory: ", OUT_DIR)



