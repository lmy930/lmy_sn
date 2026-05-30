# ============================================================
# 实际数据窗口筛选：数学严格真子集 + IDW 空间同步性
# 单变量 PM10/PM2.5 vs 双变量 PM10+PM2.5
#
# 主要目标：
#   1. 筛选 6 个站点 PM10 和 PM2.5 均无缺失的时间窗口；
#   2. 用 IDW 构造目标站点的空间加权参考序列；
#   3. 结合 6 站点同步性和目标站点-IDW参考序列一致性筛选正常背景窗口；
#   4. 对目标站点做：
#        单变量均值变点检测：single_pollutant，可设为 "PM25" 或 "PM10"；
#        双变量均值变点检测：PM10 + PM25；
#   5. 筛选数学意义上的严格真子集：
#        单变量变点集 ⊂ 双变量变点集；
#   6. 输出两类图：
#        图1：单变量检测 vs 双变量检测；
#        图2：目标站点 vs IDW加权参考值 + 6站点原始序列。
#
# 重要设置：
#   - SNSeg 参数固定：confidence = 0.95, grid_size_scale = 0.05；
#   - 检测使用全局标准化数据；
#   - 分段最短长度 min_segment_length 已由 12 改为 24；
#   - 单变量变点集可以为空集；若双变量变点集非空，则空集是真子集。
# ============================================================

library(tidyverse)
library(lubridate)
library(ggplot2)
library(SNSeg)

# ----------------------------
# 用户需要确认或修改的参数
# ----------------------------
data_dir <- "D:/03/2025/Beijing"
target_station <- "Dongsi"

# 单变量检测使用哪个变量：
#   "PM25" 表示 PM2.5；
#   "PM10" 表示 PM10。
single_pollutant <- "PM25"

# 只考虑 2020 年之后
start_date_filter <- as.POSIXct("2020-01-01 00:00:00", tz = "Asia/Shanghai")

# 输出文件夹
out_dir <- file.path(data_dir, paste0("screen_exact_subset_idw_", single_pollutant, "_vs_dual"))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# SNSeg 参数
sn_confidence <- 0.95
sn_grid_size_scale <- 0.05

# 滑动窗口设置
window_hours_vec <- c(72, 96, 120, 168, 240)
step_hours <- 12

# 每种窗口长度下，保留空间同步性最高的若干窗口进入 SNSeg 检测
max_spatial_candidates_per_length <- 300

# 分段最短长度：由 12 改为 24，保证阶段均值更稳定
min_segment_length <- 24

# 双变量变点至少 1 个；不设置变点数量上限
min_dual_cp <- 1

# IDW 幂参数
idw_power <- 2

# 空间同步性筛选阈值
min_target_cor <- 0.40
min_direction_agreement <- 0.40
max_deviation_rate <- 0.15

# IDW参考序列一致性阈值
min_idw_cor <- 0.45
min_idw_direction_agreement <- 0.40

# ------------------------------------------------------------
# 站点经纬度：请根据你的数据来源进一步核对。
# 如果你有更准确经纬度，直接修改这里即可。
# ------------------------------------------------------------
station_coords <- tibble::tribble(
  ~station,          ~lat,    ~lon,
  "Aotizhongxin",    39.982,  116.397,
  "Dongsi",          39.929,  116.417,
  "Guanyuan",        39.929,  116.339,
  "Nongzhanguan",   39.933,  116.461,
  "Tiantan",         39.886,  116.407,
  "Wanshouxigong",  39.878,  116.352
)


# ============================================================
# 1. 读取并整理数据
# ============================================================

read_one_station <- function(file) {
  station_name <- basename(file) |>
    str_remove("^Beijing_") |>
    str_remove("\\.csv$")

  dat <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)

  need_vars <- c("year", "month", "day", "hour", "PM2.5", "PM10")
  miss_vars <- setdiff(need_vars, names(dat))
  if (length(miss_vars) > 0) {
    stop(paste0("文件 ", basename(file), " 缺少变量：", paste(miss_vars, collapse = ", ")))
  }

  dat |>
    transmute(
      station = station_name,
      datetime = make_datetime(year, month, day, hour, tz = "Asia/Shanghai"),
      PM25 = as.numeric(`PM2.5`),
      PM10 = as.numeric(PM10)
    ) |>
    filter(datetime >= start_date_filter)
}

files <- list.files(data_dir, pattern = "^Beijing_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("未找到 Beijing_*.csv 文件，请检查 data_dir。")

raw_long <- purrr::map_dfr(files, read_one_station)

if (!(target_station %in% unique(raw_long$station))) {
  stop(paste0("未找到目标站点：", target_station, "。请检查文件名。"))
}

wide_dat <- raw_long |>
  pivot_wider(
    id_cols = datetime,
    names_from = station,
    values_from = c(PM25, PM10),
    names_glue = "{.value}_{station}"
  ) |>
  arrange(datetime)

# 补齐逐小时时间序列。若原始数据时间断裂，会出现 NA，后续会被“无缺失”规则排除。
full_time <- tibble(
  datetime = seq(
    min(wide_dat$datetime, na.rm = TRUE),
    max(wide_dat$datetime, na.rm = TRUE),
    by = "hour"
  )
)

wide_dat <- full_time |>
  left_join(wide_dat, by = "datetime") |>
  arrange(datetime)

missing_coord_stations <- setdiff(unique(raw_long$station), station_coords$station)
if (length(missing_coord_stations) > 0) {
  stop(paste0("以下站点缺少经纬度，请补充 station_coords：", paste(missing_coord_stations, collapse = ", ")))
}

cat("读取完成。\n")
cat("站点：", paste(sort(unique(raw_long$station)), collapse = ", "), "\n")
cat("数据范围：", format(min(wide_dat$datetime)), " 至 ", format(max(wide_dat$datetime)), "\n")
cat("单变量检测变量：", single_pollutant, "\n")


# ============================================================
# 2. 全局标准化与 IDW 加权参考序列
# ============================================================

get_station_cols <- function(dat, pollutant) {
  grep(paste0("^", pollutant, "_"), names(dat), value = TRUE)
}

get_target_col <- function(pollutant, station = target_station) {
  paste0(pollutant, "_", station)
}

# 全局标准化：对 2020 年之后所有站点数据统一计算每个污染物的均值和标准差
global_pm10_values <- as.numeric(as.matrix(wide_dat[, get_station_cols(wide_dat, "PM10"), drop = FALSE]))
global_pm25_values <- as.numeric(as.matrix(wide_dat[, get_station_cols(wide_dat, "PM25"), drop = FALSE]))

global_scale_info <- tibble(
  pollutant = c("PM10", "PM25"),
  mean = c(mean(global_pm10_values, na.rm = TRUE), mean(global_pm25_values, na.rm = TRUE)),
  sd = c(sd(global_pm10_values, na.rm = TRUE), sd(global_pm25_values, na.rm = TRUE))
)

write.csv(
  global_scale_info,
  file.path(out_dir, "00_global_scale_info.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

add_global_z_cols <- function(dat) {
  pm10_mean <- global_scale_info$mean[global_scale_info$pollutant == "PM10"]
  pm10_sd <- global_scale_info$sd[global_scale_info$pollutant == "PM10"]
  pm25_mean <- global_scale_info$mean[global_scale_info$pollutant == "PM25"]
  pm25_sd <- global_scale_info$sd[global_scale_info$pollutant == "PM25"]

  out <- dat

  for (col in get_station_cols(dat, "PM10")) {
    out[[paste0(col, "_z")]] <- (out[[col]] - pm10_mean) / pm10_sd
  }
  for (col in get_station_cols(dat, "PM25")) {
    out[[paste0(col, "_z")]] <- (out[[col]] - pm25_mean) / pm25_sd
  }

  out
}

wide_dat <- add_global_z_cols(wide_dat)

# Haversine 距离，单位 km
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R <- 6371
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * R * asin(sqrt(a))
}

compute_idw_weights <- function(coords, target_station, power = 2) {
  target_coord <- coords |> filter(station == target_station)
  other_coords <- coords |> filter(station != target_station)

  dist <- haversine_km(
    lat1 = target_coord$lat,
    lon1 = target_coord$lon,
    lat2 = other_coords$lat,
    lon2 = other_coords$lon
  )

  w_raw <- 1 / (dist^power)
  w <- w_raw / sum(w_raw)

  tibble(
    station = other_coords$station,
    distance_km = dist,
    weight = w
  ) |>
    arrange(desc(weight))
}

idw_weights <- compute_idw_weights(station_coords, target_station, idw_power)

write.csv(
  idw_weights,
  file.path(out_dir, "00_idw_weights.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("IDW 权重：\n")
print(idw_weights)

add_idw_reference_cols <- function(dat, weights) {
  out <- dat

  for (pollutant in c("PM10", "PM25")) {
    cols <- paste0(pollutant, "_", weights$station)
    z_cols <- paste0(cols, "_z")

    out[[paste0(pollutant, "_", target_station, "_IDW")]] <-
      as.numeric(as.matrix(out[, cols, drop = FALSE]) %*% weights$weight)

    out[[paste0(pollutant, "_", target_station, "_IDW_z")]] <-
      as.numeric(as.matrix(out[, z_cols, drop = FALSE]) %*% weights$weight)
  }

  out
}

wide_dat <- add_idw_reference_cols(wide_dat, idw_weights)


# ============================================================
# 3. 空间同步性工具函数
# ============================================================

safe_cor <- function(x, y) {
  suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
}

has_complete_required_data <- function(dat_window) {
  pm10_cols <- get_station_cols(dat_window, "PM10")
  pm25_cols <- get_station_cols(dat_window, "PM25")
  raw_cols <- c(pm10_cols, pm25_cols)
  raw_cols <- raw_cols[!str_detect(raw_cols, "_z$|_IDW$|_IDW_z$")]
  all(!is.na(as.matrix(dat_window[, raw_cols, drop = FALSE])))
}

mean_pairwise_cor_among_stations <- function(dat_window, pollutant) {
  cols <- get_station_cols(dat_window, pollutant)
  cols <- cols[!str_detect(cols, "_z$|_IDW$|_IDW_z$")]
  if (length(cols) < 2) return(NA_real_)
  cm <- suppressWarnings(cor(dat_window[, cols, drop = FALSE], use = "pairwise.complete.obs"))
  mean(cm[upper.tri(cm)], na.rm = TRUE)
}

target_other_cor <- function(dat_window, pollutant, station = target_station) {
  target_col <- get_target_col(pollutant, station)
  other_cols <- setdiff(get_station_cols(dat_window, pollutant), target_col)
  other_cols <- other_cols[!str_detect(other_cols, "_z$|_IDW$|_IDW_z$")]
  if (!(target_col %in% names(dat_window)) || length(other_cols) == 0) return(NA_real_)
  vals <- sapply(other_cols, function(col) safe_cor(dat_window[[target_col]], dat_window[[col]]))
  mean(vals, na.rm = TRUE)
}

target_idw_cor <- function(dat_window, pollutant, station = target_station) {
  target_col <- get_target_col(pollutant, station)
  idw_col <- paste0(pollutant, "_", station, "_IDW")
  safe_cor(dat_window[[target_col]], dat_window[[idw_col]])
}

direction_agreement_two_series <- function(x, y, eps = 0) {
  dx <- diff(x)
  dy <- diff(y)
  ok <- !is.na(dx) & !is.na(dy) & abs(dx) > eps & abs(dy) > eps
  if (sum(ok) == 0) return(NA_real_)
  mean(sign(dx[ok]) == sign(dy[ok]))
}

target_other_direction_agreement <- function(dat_window,
                                             pollutant,
                                             station = target_station,
                                             min_agree_other = 3,
                                             eps = 0) {
  target_col <- get_target_col(pollutant, station)
  other_cols <- setdiff(get_station_cols(dat_window, pollutant), target_col)
  other_cols <- other_cols[!str_detect(other_cols, "_z$|_IDW$|_IDW_z$")]
  if (!(target_col %in% names(dat_window)) || length(other_cols) < min_agree_other) return(NA_real_)

  d0 <- diff(dat_window[[target_col]])
  dm <- apply(dat_window[, other_cols, drop = FALSE], 2, diff)
  if (is.null(dim(dm))) dm <- matrix(dm, ncol = 1)

  agree <- rep(NA, length(d0))

  for (i in seq_along(d0)) {
    if (is.na(d0[i]) || abs(d0[i]) <= eps) next
    ds <- dm[i, ]
    if (sum(!is.na(ds)) < min_agree_other) next

    if (d0[i] > eps) {
      agree[i] <- sum(ds > eps, na.rm = TRUE) >= min_agree_other
    } else {
      agree[i] <- sum(ds < -eps, na.rm = TRUE) >= min_agree_other
    }
  }

  mean(agree, na.rm = TRUE)
}

target_idw_direction_agreement <- function(dat_window, pollutant, station = target_station) {
  target_col <- get_target_col(pollutant, station)
  idw_col <- paste0(pollutant, "_", station, "_IDW")
  direction_agreement_two_series(dat_window[[target_col]], dat_window[[idw_col]])
}

target_local_deviation_rate <- function(dat_window,
                                        pollutant,
                                        station = target_station,
                                        z_cut = 3.5) {
  target_col <- get_target_col(pollutant, station)
  other_cols <- setdiff(get_station_cols(dat_window, pollutant), target_col)
  other_cols <- other_cols[!str_detect(other_cols, "_z$|_IDW$|_IDW_z$")]
  if (!(target_col %in% names(dat_window)) || length(other_cols) < 3) return(NA_real_)

  x0 <- dat_window[[target_col]]
  xm <- as.matrix(dat_window[, other_cols, drop = FALSE])

  flag <- rep(FALSE, length(x0))

  for (i in seq_along(x0)) {
    others <- xm[i, ]
    med <- median(others, na.rm = TRUE)
    mad_val <- mad(others, constant = 1.4826, na.rm = TRUE)

    if (is.na(mad_val) || mad_val == 0) {
      flag[i] <- FALSE
    } else {
      flag[i] <- abs((x0[i] - med) / mad_val) > z_cut
    }
  }

  mean(flag, na.rm = TRUE)
}

calc_spatial_metrics <- function(dat_window, station = target_station) {
  pm10_cols <- get_station_cols(dat_window, "PM10")
  pm25_cols <- get_station_cols(dat_window, "PM25")
  pm10_cols <- pm10_cols[!str_detect(pm10_cols, "_z$|_IDW$|_IDW_z$")]
  pm25_cols <- pm25_cols[!str_detect(pm25_cols, "_z$|_IDW$|_IDW_z$")]
  vals <- as.matrix(dat_window[, c(pm10_cols, pm25_cols), drop = FALSE])

  missing_rate <- mean(is.na(vals))
  negative_rate <- mean(vals < 0, na.rm = TRUE)

  pm10_target_cor <- target_other_cor(dat_window, "PM10", station)
  pm25_target_cor <- target_other_cor(dat_window, "PM25", station)
  pm10_all_cor <- mean_pairwise_cor_among_stations(dat_window, "PM10")
  pm25_all_cor <- mean_pairwise_cor_among_stations(dat_window, "PM25")
  pm10_dir <- target_other_direction_agreement(dat_window, "PM10", station)
  pm25_dir <- target_other_direction_agreement(dat_window, "PM25", station)
  pm10_dev <- target_local_deviation_rate(dat_window, "PM10", station)
  pm25_dev <- target_local_deviation_rate(dat_window, "PM25", station)
  deviation_rate <- mean(c(pm10_dev, pm25_dev), na.rm = TRUE)

  pm10_idw_cor <- target_idw_cor(dat_window, "PM10", station)
  pm25_idw_cor <- target_idw_cor(dat_window, "PM25", station)
  pm10_idw_dir <- target_idw_direction_agreement(dat_window, "PM10", station)
  pm25_idw_dir <- target_idw_direction_agreement(dat_window, "PM25", station)

  spatial_score <-
    0.15 * pm10_target_cor +
    0.15 * pm25_target_cor +
    0.10 * pm10_all_cor +
    0.10 * pm25_all_cor +
    0.10 * pm10_dir +
    0.10 * pm25_dir +
    0.15 * pm10_idw_cor +
    0.15 * pm25_idw_cor +
    0.05 * pm10_idw_dir +
    0.05 * pm25_idw_dir -
    0.60 * deviation_rate -
    1.00 * negative_rate

  tibble(
    missing_rate = missing_rate,
    negative_rate = negative_rate,
    pm10_target_other_cor = pm10_target_cor,
    pm25_target_other_cor = pm25_target_cor,
    pm10_all_station_cor = pm10_all_cor,
    pm25_all_station_cor = pm25_all_cor,
    pm10_direction_agreement = pm10_dir,
    pm25_direction_agreement = pm25_dir,
    pm10_deviation_rate = pm10_dev,
    pm25_deviation_rate = pm25_dev,
    deviation_rate = deviation_rate,
    pm10_target_idw_cor = pm10_idw_cor,
    pm25_target_idw_cor = pm25_idw_cor,
    pm10_target_idw_direction_agreement = pm10_idw_dir,
    pm25_target_idw_direction_agreement = pm25_idw_dir,
    spatial_score = spatial_score
  )
}


# ============================================================
# 4. 严格真子集与变点强度函数
# ============================================================

is_exact_subset <- function(A, B) {
  A <- sort(unique(A))
  B <- sort(unique(B))
  if (length(A) == 0) return(TRUE)
  all(A %in% B)
}

is_exact_proper_subset <- function(A, B) {
  A <- sort(unique(A))
  B <- sort(unique(B))
  is_exact_subset(A, B) && length(setdiff(B, A)) >= 1
}

extra_cp_count_exact <- function(big_set, small_set) {
  length(setdiff(sort(unique(big_set)), sort(unique(small_set))))
}

segment_lengths_from_cp <- function(cp_index, n) {
  cuts <- c(0, sort(unique(cp_index)), n)
  diff(cuts)
}

cp_spacing_ok <- function(cp_index, n, min_seg_len = min_segment_length) {
  if (length(cp_index) == 0) return(FALSE)
  min(segment_lengths_from_cp(cp_index, n)) >= min_seg_len
}

mean_shift_strength <- function(dat_window, cols, cp_index, min_seg_len = min_segment_length) {
  if (length(cp_index) == 0) return(NA_real_)

  x <- dat_window |> select(all_of(cols)) |> as.matrix()
  n <- nrow(x)
  cps <- sort(unique(cp_index))
  cuts <- c(0, cps, n)
  seg_lengths <- diff(cuts)

  if (any(seg_lengths < min_seg_len)) return(NA_real_)

  means <- vector("list", length(cuts) - 1)
  for (j in seq_len(length(cuts) - 1)) {
    idx <- (cuts[j] + 1):cuts[j + 1]
    means[[j]] <- colMeans(x[idx, , drop = FALSE])
  }

  jumps <- numeric(length(means) - 1)
  for (j in seq_len(length(means) - 1)) {
    jumps[j] <- sqrt(sum((means[[j + 1]] - means[[j]])^2))
  }

  mean(jumps, na.rm = TRUE)
}


# ============================================================
# 5. SNSeg 检测函数
# ============================================================

extract_cp_from_snseg <- function(res) {
  if (is.null(res)) return(integer(0))

  possible_names <- c("est_cp", "cp", "cpts", "changepoints")
  for (nm in possible_names) {
    if (nm %in% names(res)) {
      cp <- res[[nm]]
      if (length(cp) == 0 || all(is.na(cp))) return(integer(0))
      return(as.integer(cp))
    }
  }

  integer(0)
}

run_snseg_uni_mean <- function(dat_window, col_z) {
  x <- dat_window[[col_z]]
  if (anyNA(x)) stop("run_snseg_uni_mean 输入含缺失值，请先过滤。")
  if (length(x) < 50) {
    return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  }

  res <- tryCatch(
    SNSeg_Uni(
      ts = x,
      paras_to_test = "mean",
      confidence = sn_confidence,
      grid_size_scale = sn_grid_size_scale,
      grid_size = NULL,
      plot_SN = FALSE,
      est_cp_loc = TRUE
    ),
    error = function(e) {
      message("SNSeg_Uni 出错：", e$message)
      NULL
    }
  )

  if (is.null(res)) {
    return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  }

  cp <- extract_cp_from_snseg(res)
  list(ok = TRUE, n_cp = length(cp), cp_index = cp, cp_datetime = dat_window$datetime[cp])
}

run_snseg_multi_mean <- function(dat_window, cols_z) {
  x <- dat_window |> select(all_of(cols_z)) |> as.matrix()
  if (anyNA(x)) stop("run_snseg_multi_mean 输入含缺失值，请先过滤。")
  if (nrow(x) < 50) {
    return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  }

  res <- tryCatch(
    SNSeg_Multi(
      ts = x,
      paras_to_test = "mean",
      confidence = sn_confidence,
      grid_size_scale = sn_grid_size_scale,
      grid_size = NULL,
      plot_SN = FALSE,
      est_cp_loc = TRUE
    ),
    error = function(e) {
      message("SNSeg_Multi 出错：", e$message)
      NULL
    }
  )

  if (is.null(res)) {
    return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  }

  cp <- extract_cp_from_snseg(res)
  list(ok = TRUE, n_cp = length(cp), cp_index = cp, cp_datetime = dat_window$datetime[cp])
}


# ============================================================
# 6. 单窗口检测与筛选
# ============================================================

evaluate_window <- function(dat_window, station = target_station) {
  single_col_z <- paste0(single_pollutant, "_", station, "_z")
  pm10_col_z <- paste0("PM10_", station, "_z")
  pm25_col_z <- paste0("PM25_", station, "_z")

  res_single <- run_snseg_uni_mean(dat_window, single_col_z)
  res_dual <- run_snseg_multi_mean(dat_window, c(pm10_col_z, pm25_col_z))

  subset_flag <- is_exact_subset(
    A = res_single$cp_index,
    B = res_dual$cp_index
  )

  proper_subset_flag <- is_exact_proper_subset(
    A = res_single$cp_index,
    B = res_dual$cp_index
  )

  dual_extra_n <- extra_cp_count_exact(
    big_set = res_dual$cp_index,
    small_set = res_single$cp_index
  )

  dual_spacing <- cp_spacing_ok(res_dual$cp_index, nrow(dat_window))

  dual_strength <- mean_shift_strength(
    dat_window,
    cols = c(pm10_col_z, pm25_col_z),
    cp_index = res_dual$cp_index
  )

  single_strength <- mean_shift_strength(
    dat_window,
    cols = c(single_col_z),
    cp_index = res_single$cp_index
  )

  pass_cp <- !is.na(res_dual$n_cp) &&
    res_dual$n_cp >= min_dual_cp &&
    proper_subset_flag &&
    dual_spacing &&
    !is.na(dual_strength)

  tibble(
    single_pollutant = single_pollutant,
    n_cp_single = res_single$n_cp,
    cp_index_single = paste(res_single$cp_index, collapse = ";"),
    cp_datetime_single = paste(format(res_single$cp_datetime, "%Y-%m-%d %H:%M"), collapse = ";"),

    n_cp_dual = res_dual$n_cp,
    cp_index_dual = paste(res_dual$cp_index, collapse = ";"),
    cp_datetime_dual = paste(format(res_dual$cp_datetime, "%Y-%m-%d %H:%M"), collapse = ";"),

    exact_single_subset_of_dual = subset_flag,
    exact_single_proper_subset_of_dual = proper_subset_flag,
    dual_extra_cp_count = dual_extra_n,
    dual_spacing_ok = dual_spacing,
    single_strength = single_strength,
    dual_strength = dual_strength,
    pass_cp = pass_cp
  )
}


# ============================================================
# 7. 筛选候选窗口
# ============================================================

pre_screen_spatial_windows <- function(wide_dat) {
  result_list <- list()
  idx <- 1
  n_total <- nrow(wide_dat)

  for (L in window_hours_vec) {
    cat("\n空间同步性初筛：窗口长度 =", L, "小时\n")

    starts <- seq(1, n_total - L + 1, by = step_hours)

    for (s in starts) {
      e <- s + L - 1
      dat_window <- wide_dat[s:e, ]

      if (!has_complete_required_data(dat_window)) next

      metrics <- calc_spatial_metrics(dat_window)

      pass <- with(
        metrics,
        missing_rate == 0 &&
          negative_rate <= 0 &&
          pm10_target_other_cor >= min_target_cor &&
          pm25_target_other_cor >= min_target_cor &&
          pm10_direction_agreement >= min_direction_agreement &&
          pm25_direction_agreement >= min_direction_agreement &&
          pm10_target_idw_cor >= min_idw_cor &&
          pm25_target_idw_cor >= min_idw_cor &&
          pm10_target_idw_direction_agreement >= min_idw_direction_agreement &&
          pm25_target_idw_direction_agreement >= min_idw_direction_agreement &&
          deviation_rate <= max_deviation_rate
      )

      result_list[[idx]] <- metrics |>
        mutate(
          window_hours = L,
          start_index = s,
          end_index = e,
          start_time = dat_window$datetime[1],
          end_time = dat_window$datetime[nrow(dat_window)],
          pass_spatial = pass
        )

      idx <- idx + 1
    }
  }

  if (length(result_list) == 0) return(tibble())

  bind_rows(result_list) |>
    arrange(desc(pass_spatial), desc(spatial_score))
}

find_candidate_windows <- function(wide_dat) {
  cat("\n============================================================\n")
  cat("开始筛选数学意义上的严格真子集窗口：单变量变点集 ⊂ 双变量变点集\n")
  cat("单变量：", single_pollutant, "\n")
  cat("============================================================\n")

  spatial_res <- pre_screen_spatial_windows(wide_dat)

  write.csv(
    spatial_res,
    file.path(out_dir, "01_spatial_screen_complete_windows.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  if (nrow(spatial_res) == 0) {
    warning("没有无缺失窗口。")
    return(list(all = tibble(), candidates = tibble()))
  }

  spatial_keep <- spatial_res |>
    filter(pass_spatial) |>
    group_by(window_hours) |>
    slice_max(order_by = spatial_score, n = max_spatial_candidates_per_length, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(spatial_score))

  write.csv(
    spatial_keep,
    file.path(out_dir, "02_spatial_kept_for_snseg.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  cat("进入 SNSeg 检测的窗口数：", nrow(spatial_keep), "\n")

  if (nrow(spatial_keep) == 0) {
    warning("没有窗口通过空间同步性初筛。可适当放宽空间同步性阈值。")
    return(list(all = tibble(), candidates = tibble()))
  }

  result_list <- list()

  for (i in seq_len(nrow(spatial_keep))) {
    s <- spatial_keep$start_index[i]
    e <- spatial_keep$end_index[i]
    dat_window <- wide_dat[s:e, ]

    det <- evaluate_window(dat_window)

    result_list[[i]] <- spatial_keep[i, ] |>
      bind_cols(det) |>
      mutate(
        step_hours = step_hours,
        min_segment_length = min_segment_length,
        confidence = sn_confidence,
        grid_size_scale = sn_grid_size_scale,
        idw_power = idw_power,
        global_standardization = TRUE,
        final_score =
          0.35 * spatial_score +
          0.35 * pmin(ifelse(is.na(dual_strength), 0, dual_strength) / 1.0, 1) +
          0.20 * pmin(dual_extra_cp_count / 3, 1) +
          0.10 * ifelse(exact_single_proper_subset_of_dual, 1, 0)
      )

    if (i %% 20 == 0) cat("已完成 SNSeg 检测：", i, "/", nrow(spatial_keep), "\n")
  }

  all_res <- bind_rows(result_list) |>
    arrange(desc(pass_cp), desc(final_score), desc(spatial_score))

  candidates <- all_res |>
    filter(pass_cp) |>
    arrange(desc(final_score), desc(dual_extra_cp_count), desc(spatial_score))

  write.csv(
    all_res,
    file.path(out_dir, "03_all_windows_detection_results.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    candidates,
    file.path(out_dir, "04_exact_strict_subset_candidates.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  list(all = all_res, candidates = candidates)
}


# ============================================================
# 8. 绘图函数
# ============================================================

parse_cp_times <- function(cp_text) {
  if (is.na(cp_text) || cp_text == "") return(as.POSIXct(character(0)))
  ymd_hm(str_split(cp_text, ";")[[1]], tz = "Asia/Shanghai")
}

plot_two_panel_detection <- function(candidate_row,
                                     wide_dat,
                                     station = target_station,
                                     prefix = "candidate") {
  dat_window <- wide_dat[candidate_row$start_index:candidate_row$end_index, ]

  single_col_raw <- paste0(single_pollutant, "_", station)
  pm10_col_z <- paste0("PM10_", station, "_z")
  pm25_col_z <- paste0("PM25_", station, "_z")

  cp_single <- parse_cp_times(candidate_row$cp_datetime_single)
  cp_dual <- parse_cp_times(candidate_row$cp_datetime_dual)

  single_label <- ifelse(single_pollutant == "PM25", "PM2.5", "PM10")

  upper_dat <- dat_window |>
    transmute(
      datetime = datetime,
      value = .data[[single_col_raw]],
      panel = paste0("① 使用 ", single_label, " 做单变量均值变点检测")
    )

  lower_dat <- tibble(
    datetime = rep(dat_window$datetime, 2),
    value = c(dat_window[[pm10_col_z]], dat_window[[pm25_col_z]]),
    variable = rep(c("PM10（全局标准化）", "PM2.5（全局标准化）"), each = nrow(dat_window)),
    panel = "② 使用 PM10 和 PM2.5 做双变量均值变点检测"
  )

  p <- ggplot() +
    geom_line(data = upper_dat, aes(x = datetime, y = value), linewidth = 0.70) +
    geom_line(data = lower_dat, aes(x = datetime, y = value, linetype = variable), linewidth = 0.70) +
    geom_vline(
      data = tibble(x = cp_single, panel = paste0("① 使用 ", single_label, " 做单变量均值变点检测")),
      aes(xintercept = as.numeric(x)),
      linetype = "dashed",
      linewidth = 0.75
    ) +
    geom_vline(
      data = tibble(x = cp_dual, panel = "② 使用 PM10 和 PM2.5 做双变量均值变点检测"),
      aes(xintercept = as.numeric(x)),
      linetype = "solid",
      linewidth = 0.80
    ) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(
      title = paste0(
        single_label, " 单变量与 PM10-PM2.5 双变量检测对比：",
        format(dat_window$datetime[1], "%Y-%m-%d %H:%M"),
        " 至 ",
        format(dat_window$datetime[nrow(dat_window)], "%Y-%m-%d %H:%M")
      ),
      subtitle = "上图虚线为单变量变点；下图实线为双变量变点。筛选条件：单变量变点集为双变量变点集的数学真子集。",
      x = "时间",
      y = "观测值 / 全局标准化观测值",
      linetype = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  print(p)

  fname <- paste0(
    prefix, "_detection_two_panel_",
    single_pollutant, "_",
    format(dat_window$datetime[1], "%Y%m%d%H"),
    "_",
    format(dat_window$datetime[nrow(dat_window)], "%Y%m%d%H"),
    ".png"
  )

  ggsave(file.path(out_dir, fname), p, width = 12, height = 7, dpi = 300)
  invisible(p)
}

plot_idw_spatial_overview <- function(candidate_row,
                                      wide_dat,
                                      station = target_station,
                                      prefix = "candidate") {
  dat_window <- wide_dat[candidate_row$start_index:candidate_row$end_index, ]

  pm10_target <- paste0("PM10_", station)
  pm25_target <- paste0("PM25_", station)
  pm10_idw <- paste0("PM10_", station, "_IDW")
  pm25_idw <- paste0("PM25_", station, "_IDW")

  idw_dat <- bind_rows(
    dat_window |> transmute(datetime, value = .data[[pm10_target]], series = "目标站点", panel = "① PM10：目标站点与IDW加权参考值"),
    dat_window |> transmute(datetime, value = .data[[pm10_idw]], series = "IDW加权参考值", panel = "① PM10：目标站点与IDW加权参考值"),
    dat_window |> transmute(datetime, value = .data[[pm25_target]], series = "目标站点", panel = "② PM2.5：目标站点与IDW加权参考值"),
    dat_window |> transmute(datetime, value = .data[[pm25_idw]], series = "IDW加权参考值", panel = "② PM2.5：目标站点与IDW加权参考值")
  )

  pm10_cols <- get_station_cols(dat_window, "PM10")
  pm25_cols <- get_station_cols(dat_window, "PM25")
  pm10_cols <- pm10_cols[!str_detect(pm10_cols, "_z$|_IDW$|_IDW_z$")]
  pm25_cols <- pm25_cols[!str_detect(pm25_cols, "_z$|_IDW$|_IDW_z$")]

  raw_dat <- dat_window |>
    select(datetime, all_of(pm10_cols), all_of(pm25_cols)) |>
    pivot_longer(cols = -datetime, names_to = "name", values_to = "value") |>
    mutate(
      pollutant = if_else(str_detect(name, "^PM10_"), "PM10", "PM2.5"),
      station_name = name |>
        str_remove("^PM10_") |>
        str_remove("^PM25_"),
      series = station_name,
      panel = if_else(
        pollutant == "PM10",
        "③ PM10：6个站点原始浓度变化",
        "④ PM2.5：6个站点原始浓度变化"
      )
    )

  p <- ggplot() +
    geom_line(
      data = idw_dat,
      aes(x = datetime, y = value, linetype = series),
      linewidth = 0.75
    ) +
    geom_line(
      data = raw_dat,
      aes(x = datetime, y = value, group = series),
      linewidth = 0.45,
      alpha = 0.60
    ) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(
      title = paste0(
        "目标站点、IDW加权参考值与6站点原始序列：",
        format(dat_window$datetime[1], "%Y-%m-%d %H:%M"),
        " 至 ",
        format(dat_window$datetime[nrow(dat_window)], "%Y-%m-%d %H:%M")
      ),
      subtitle = "前两行为目标站点与临近站点IDW加权参考值；后两行为6个站点原始浓度变化趋势。",
      x = "时间",
      y = "污染物浓度",
      linetype = NULL
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )

  print(p)

  fname <- paste0(
    prefix, "_idw_spatial_overview_",
    format(dat_window$datetime[1], "%Y%m%d%H"),
    "_",
    format(dat_window$datetime[nrow(dat_window)], "%Y%m%d%H"),
    ".png"
  )

  ggsave(file.path(out_dir, fname), p, width = 12, height = 10, dpi = 300)
  invisible(p)
}


# ============================================================
# 9. 执行筛选与作图
# ============================================================

screen_res <- find_candidate_windows(wide_dat)
all_results <- screen_res$all
target_candidates <- screen_res$candidates

cat("\n============================================================\n")
cat("筛选完成。\n")
cat("输出文件夹：", out_dir, "\n")
cat("满足数学严格真子集条件的候选窗口数：", nrow(target_candidates), "\n")
cat("============================================================\n")

cat("\n候选窗口前 10 个：\n")
print(head(target_candidates, 10))

if (nrow(target_candidates) > 0) {
  for (i in seq_len(min(10, nrow(target_candidates)))) {
    plot_two_panel_detection(
      candidate_row = target_candidates[i, ],
      wide_dat = wide_dat,
      prefix = paste0("exact_subset_top_", i)
    )

    plot_idw_spatial_overview(
      candidate_row = target_candidates[i, ],
      wide_dat = wide_dat,
      prefix = paste0("exact_subset_top_", i)
    )
  }
}

cat("\n输出文件说明：\n")
cat("00_global_scale_info.csv：全局标准化均值和标准差。\n")
cat("00_idw_weights.csv：IDW距离和权重。\n")
cat("01_spatial_screen_complete_windows.csv：所有无缺失滑动窗口的空间同步性指标。\n")
cat("02_spatial_kept_for_snseg.csv：进入 SNSeg 检测的空间同步窗口。\n")
cat("03_all_windows_detection_results.csv：所有检测窗口的结果。\n")
cat("04_exact_strict_subset_candidates.csv：满足数学严格真子集条件的候选窗口。\n")
cat("*_detection_two_panel_*.png：单变量与双变量检测结果图。\n")
cat("*_idw_spatial_overview_*.png：目标站点、IDW参考值和6站点原始序列图。\n")
cat("============================================================\n")
