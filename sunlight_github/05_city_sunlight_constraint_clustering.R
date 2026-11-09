# Three-feature descriptive typology of city sunlight constraints.
# Clustering variables: total_loss_per_day, cloud_share, aerosol_noncloud_share.
# K-means, K = 4, z-scored features and Euclidean distances.

suppressPackageStartupMessages({
  library(cluster)
  library(dplyr)
  library(readr)
})

options(encoding = "UTF-8")
set.seed(20260928)

PROJECT_ROOT <- Sys.getenv("SUNLIGHT_PROJECT_ROOT", unset = ".")
DATA_DIR <- Sys.getenv("SUNLIGHT_DATA_DIR", unset = file.path(PROJECT_ROOT, "data", "restricted"))
OUTPUT_ROOT <- Sys.getenv("SUNLIGHT_OUTPUT_DIR", unset = file.path(PROJECT_ROOT, "outputs"))
feature_file <- Sys.getenv("CITY_FEATURE_FILE", unset = file.path(DATA_DIR, "city_sunlight_constraint_features.csv"))
out_dir <- file.path(OUTPUT_ROOT, "city_sunlight_constraint_clustering")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_out_csv <- function(x, filename) {
  readr::write_excel_csv(x, file.path(out_dir, filename), na = "")
}

zscore <- function(v) {
  s <- sd(v, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(v)))
  as.numeric((v - mean(v, na.rm = TRUE)) / s)
}

message("Reading city-level sunlight-constraint feature data...")
city_features <- read_csv(
  feature_file,
  locale = locale(encoding = "UTF-8"),
  col_types = cols(CITYCODE = col_character(), CITYNAME = col_character(), .default = col_guess()),
  show_col_types = FALSE
) %>% arrange(CITYCODE)

stopifnot(nrow(city_features) == 286L)
stopifnot(!anyDuplicated(city_features$CITYCODE))
feature_cols <- c("total_loss_per_day", "cloud_share", "aerosol_noncloud_share")
if (anyNA(city_features[, feature_cols])) {
  stop("The three clustering variables contain missing values.")
}
feature_matrix <- city_features %>% select(all_of(feature_cols)) %>% as.data.frame()
rownames(feature_matrix) <- paste(city_features$CITYCODE, city_features$CITYNAME, sep = "_")
scaled_x <- scale(feature_matrix)
dist_x <- dist(scaled_x, method = "euclidean")

calc_ch <- function(x, cl) {
  cl <- as.factor(cl)
  n <- nrow(x)
  k <- nlevels(cl)
  overall <- colMeans(x)
  wss <- 0
  bss <- 0
  for (lv in levels(cl)) {
    xi <- x[cl == lv, , drop = FALSE]
    ci <- colMeans(xi)
    wss <- wss + sum(rowSums((sweep(xi, 2, ci, "-"))^2))
    bss <- bss + nrow(xi) * sum((ci - overall)^2)
  }
  (bss / (k - 1)) / (wss / (n - k))
}

calc_db <- function(x, cl) {
  cl <- as.factor(cl)
  lev <- levels(cl)
  centers <- do.call(rbind, lapply(lev, function(lv) {
    colMeans(x[cl == lv, , drop = FALSE])
  }))
  scatter <- sapply(seq_along(lev), function(i) {
    xi <- x[cl == lev[i], , drop = FALSE]
    mean(sqrt(rowSums((sweep(xi, 2, centers[i, ], "-"))^2)))
  })
  center_dist <- as.matrix(dist(centers))
  diag(center_dist) <- Inf
  db_i <- sapply(seq_along(lev), function(i) {
    max((scatter[i] + scatter) / center_dist[i, ], na.rm = TRUE)
  })
  mean(db_i)
}

calc_sil <- function(cl) {
  mean(silhouette(as.integer(as.factor(cl)), dist_x)[, "sil_width"])
}

score_solution <- function(method, k, cl) {
  tibble(
    method = method,
    k = k,
    silhouette = calc_sil(cl),
    calinski_harabasz = calc_ch(as.matrix(scaled_x), cl),
    davies_bouldin = calc_db(as.matrix(scaled_x), cl),
    min_cluster_size = min(table(cl)),
    max_cluster_size = max(table(cl))
  )
}


# Keep the K = 4 K-means model used for the city typology.
k4_fit <- kmeans(scaled_x, centers = 4, nstart = 300, iter.max = 300)
k4_clusters <- k4_fit$cluster
k4_metrics <- score_solution("kmeans", 4, k4_clusters)
write_out_csv(k4_metrics, "clustering_model_metrics_kmeans_k4.csv")

label_low_en <- "Low modifiable-loss burden"
label_aerosol_en <- "High-cloud aerosol-elevated"
label_building_en <- "Building-high"
label_compound_en <- "Compound high-loss"
label_mixed_en <- "Cloud-building mixed"

label_low_zh <- "\u4f4e\u53ef\u5e72\u9884\u8d1f\u62c5\u578b"
label_aerosol_zh <- "\u9ad8\u4e91\u635f\u5931-\u6c14\u6eb6\u80f6\u76f8\u5bf9\u504f\u9ad8\u578b"
label_building_zh <- "\u5efa\u7b51\u906e\u6321\u9ad8\u8d21\u732e\u578b"
label_compound_zh <- "\u590d\u5408\u9ad8\u635f\u5931\u578b"
label_mixed_zh <- "\u4e91-\u5efa\u7b51\u6df7\u5408\u4e2d\u7b49\u7ea6\u675f\u578b"

make_labels <- function(profile) {
  profile <- profile %>%
    mutate(
      building_noncloud_share = 1 - aerosol_noncloud_share,
      label_en = NA_character_,
      label_zh = NA_character_
    )
  assigned <- rep(FALSE, nrow(profile))

  assign_one <- function(idx, en, zh) {
    if (length(idx) == 1 && !is.na(idx) && !assigned[idx]) {
      profile$label_en[idx] <<- en
      profile$label_zh[idx] <<- zh
      assigned[idx] <<- TRUE
    }
  }

  low_idx <- which.min(
    zscore(profile$total_loss_per_day) + zscore(profile$controllable_loss_per_day)
  )
  assign_one(low_idx, label_low_en, label_low_zh)

  aerosol_score <- zscore(profile$aerosol_noncloud_share) +
    zscore(profile$aerosol_loss_per_day)
  aerosol_score[assigned] <- -Inf
  assign_one(which.max(aerosol_score), label_aerosol_en, label_aerosol_zh)

  building_score <- zscore(profile$building_noncloud_share) +
    zscore(profile$building_loss_per_day)
  building_score[assigned] <- -Inf
  assign_one(which.max(building_score), label_building_en, label_building_zh)

  compound_score <- zscore(profile$total_loss_per_day) +
    zscore(profile$cloud_loss_per_day) +
    zscore(profile$controllable_loss_per_day)
  compound_score[assigned] <- -Inf
  assign_one(which.max(compound_score), label_compound_en, label_compound_zh)

  remaining <- which(!assigned)
  if (length(remaining) > 0) {
    mixed_score <- zscore(profile$building_loss_per_day) +
      zscore(profile$cloud_loss_per_day)
    remaining <- remaining[order(mixed_score[remaining], decreasing = TRUE)]
    for (idx in remaining) assign_one(idx, label_mixed_en, label_mixed_zh)
  }
  profile
}

profile_solution <- function(cl, solution_name) {
  city_features %>%
    mutate(cluster = as.integer(as.factor(cl))) %>%
    group_by(cluster) %>%
    summarise(
      n_cities = n(),
      total_loss_per_day = mean(total_loss_per_day),
      cloud_loss_per_day = mean(cloud_loss_per_day),
      aerosol_loss_per_day = mean(aerosol_loss_per_day),
      building_loss_per_day = mean(building_loss_per_day),
      controllable_loss_per_day = mean(controllable_loss_per_day),
      cloud_share = mean(cloud_share),
      aerosol_share = mean(aerosol_share),
      building_share = mean(building_share),
      aerosol_noncloud_share = mean(aerosol_noncloud_share),
      building_noncloud_share = mean(building_noncloud_share),
      .groups = "drop"
    ) %>%
    make_labels() %>%
    mutate(solution = solution_name, .before = 1)
}


k4_profile <- profile_solution(k4_clusters, "kmeans_k4")
write_out_csv(k4_profile, "cluster_profiles_kmeans_k4.csv")

k4_results <- city_features %>%
  mutate(kmeans_k4_cluster = as.integer(as.factor(k4_clusters))) %>%
  left_join(
    k4_profile %>% select(cluster, kmeans_k4_label_en = label_en, kmeans_k4_label_zh = label_zh),
    by = c("kmeans_k4_cluster" = "cluster")
  )

k4_results <- k4_results %>% select(
  CITYCODE, CITYNAME, cloud_pinyin,
  total_loss_per_day, cloud_loss_per_day, aerosol_loss_per_day,
  building_loss_per_day, controllable_loss_per_day,
  cloud_share, aerosol_share, building_share,
  aerosol_noncloud_share, building_noncloud_share,
  kmeans_k4_cluster, kmeans_k4_label_en, kmeans_k4_label_zh
)
write_out_csv(k4_results, "city_three_feature_clusters_kmeans_k4.csv")

scaled_tbl <- as_tibble(scaled_x) %>%
  rename_with(~ paste0(.x, "_z")) %>%
  mutate(CITYCODE = city_features$CITYCODE, .before = 1)
three_feature_input <- city_features %>%
  select(
    CITYCODE, CITYNAME, cloud_pinyin,
    total_loss_per_day, cloud_share, aerosol_noncloud_share,
    cloud_loss_per_day, aerosol_loss_per_day, building_loss_per_day,
    controllable_loss_per_day, aerosol_share, building_share
  ) %>% left_join(scaled_tbl, by = "CITYCODE")
write_out_csv(three_feature_input, "city_three_feature_inputs.csv")

province_lookup <- tibble(
  province_code = c(
    "11", "12", "13", "14", "15", "21", "22", "23", "31", "32", "33",
    "34", "35", "36", "37", "41", "42", "43", "44", "45", "46", "50",
    "51", "52", "53", "54", "61", "62", "63", "64", "65"
  ),
  province = c(
    "\u5317\u4eac", "\u5929\u6d25", "\u6cb3\u5317", "\u5c71\u897f", "\u5185\u8499\u53e4",
    "\u8fbd\u5b81", "\u5409\u6797", "\u9ed1\u9f99\u6c5f", "\u4e0a\u6d77", "\u6c5f\u82cf",
    "\u6d59\u6c5f", "\u5b89\u5fbd", "\u798f\u5efa", "\u6c5f\u897f", "\u5c71\u4e1c",
    "\u6cb3\u5357", "\u6e56\u5317", "\u6e56\u5357", "\u5e7f\u4e1c", "\u5e7f\u897f",
    "\u6d77\u5357", "\u91cd\u5e86", "\u56db\u5ddd", "\u8d35\u5dde", "\u4e91\u5357",
    "\u897f\u85cf", "\u9655\u897f", "\u7518\u8083", "\u9752\u6d77", "\u5b81\u590f", "\u65b0\u7586"
  ),
  region = c(
    rep("\u534e\u5317", 5), rep("\u4e1c\u5317", 3), rep("\u534e\u4e1c", 7),
    rep("\u534e\u4e2d", 3), rep("\u534e\u5357", 3), rep("\u897f\u5357", 5),
    rep("\u897f\u5317", 5)
  )
)


regional_city_k4 <- k4_results %>%
  mutate(province_code = substr(CITYCODE, 4, 5)) %>%
  left_join(province_lookup, by = "province_code")
stopifnot(!anyNA(regional_city_k4$region))
regional_composition_k4 <- regional_city_k4 %>%
  count(kmeans_k4_label_zh, region, name = "n_cities") %>%
  group_by(kmeans_k4_label_zh) %>%
  mutate(proportion_within_cluster = n_cities / sum(n_cities)) %>%
  ungroup() %>% arrange(kmeans_k4_label_zh, desc(n_cities))
province_composition_k4 <- regional_city_k4 %>%
  count(kmeans_k4_label_zh, province, region, name = "n_cities") %>%
  group_by(kmeans_k4_label_zh) %>%
  mutate(proportion_within_cluster = n_cities / sum(n_cities)) %>%
  ungroup() %>% arrange(kmeans_k4_label_zh, desc(n_cities))
write_out_csv(regional_composition_k4, "regional_composition_kmeans_k4.csv")
write_out_csv(province_composition_k4, "province_composition_kmeans_k4.csv")

summary_file <- file.path(out_dir, "city_typology_run_summary.txt")
sink(summary_file)
cat("Three-feature city sunlight-constraint typology (K-means K=4)\n")
cat("Cities:", nrow(city_features), "\n")
cat("Variables after z-standardization:\n")
print(feature_cols)
cat("\nClustering quality:\n")
print(k4_metrics)
cat("\nCluster profiles:\n")
print(k4_profile %>% arrange(cluster))
cat("\nRegional composition:\n")
print(regional_composition_k4)
sink()
message("Results written to: ", out_dir)
