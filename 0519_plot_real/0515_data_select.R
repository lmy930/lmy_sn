# ============================================================
# 实际数据候选时间段筛选：6站点 PM10/PM2.5，目标站点 Dongsi
# 目的：
#   1) 筛选无缺失值、空间同步性较好的正常时间段；
#   2) 在候选时间段中比较：
#        单变量方法检测变量1（PM10）；
#        单变量方法检测变量2（PM2.5）；
#        双变量方法检测二维联合序列（PM10, PM2.5）；
#   3) 找到双变量检测相对单变量具有补充价值的时间段；
#   4) 输出一张三行子图：
#        第1行：变量1 + 变量1单变量检测变点；
#        第2行：变量2 + 变量2单变量检测变点；
#        第3行：变量1/变量2标准化序列 + 双变量检测变点。
#
# SNSeg 参数：confidence = 0.95, grid_size_scale = 0.05
# ============================================================

library(tidyverse)
library(lubridate)
library(ggplot2)
library(SNSeg)

# ----------------------------
# 0. 用户参数
# ----------------------------
data_dir <- "D:/03/2025/Beijing"
target_station <- "Dongsi"
start_date_filter <- as.POSIXct("2020-01-01 00:00:00", tz = "Asia/Shanghai")

var1_name <- "PM10"
var2_name <- "PM25"   # 读入时把 PM2.5 重命名为 PM25

out_dir <- file.path(data_dir, "actual_data_window_screening_v2")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

sn_confidence <- 0.95
sn_grid_size_scale <- 0.05

window_hours_vec <- c(72, 96, 120, 168)
step_hours <- 12
max_spatial_candidates_per_length <- 150
cp_match_tolerance <- 6
min_segment_length <- 12

# ----------------------------
# 1. 读取与整理数据
# ----------------------------
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

# 补齐小时时间轴；断裂处会产生 NA，后续被“无缺失”规则排除
full_time <- tibble(datetime = seq(min(wide_dat$datetime, na.rm = TRUE),
                                   max(wide_dat$datetime, na.rm = TRUE),
                                   by = "hour"))
wide_dat <- full_time |>
  left_join(wide_dat, by = "datetime") |>
  arrange(datetime)

cat("读取完成。站点：", paste(sort(unique(raw_long$station)), collapse = ", "), "\n")
cat("数据范围：", format(min(wide_dat$datetime)), " 至 ", format(max(wide_dat$datetime)), "\n")

# ----------------------------
# 2. 空间同步性工具函数
# ----------------------------
get_station_cols <- function(dat, pollutant) grep(paste0("^", pollutant, "_"), names(dat), value = TRUE)
get_target_col <- function(pollutant, station = target_station) paste0(pollutant, "_", station)
safe_cor <- function(x, y) suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))

has_complete_required_data <- function(dat_window) {
  all_cols <- c(get_station_cols(dat_window, "PM10"), get_station_cols(dat_window, "PM25"))
  all(!is.na(as.matrix(dat_window[, all_cols, drop = FALSE])))
}

mean_pairwise_cor_among_stations <- function(dat_window, pollutant) {
  cols <- get_station_cols(dat_window, pollutant)
  if (length(cols) < 2) return(NA_real_)
  cm <- suppressWarnings(cor(dat_window[, cols, drop = FALSE], use = "pairwise.complete.obs"))
  mean(cm[upper.tri(cm)], na.rm = TRUE)
}

target_other_cor <- function(dat_window, pollutant, station = target_station) {
  target_col <- get_target_col(pollutant, station)
  other_cols <- setdiff(get_station_cols(dat_window, pollutant), target_col)
  if (!(target_col %in% names(dat_window)) || length(other_cols) == 0) return(NA_real_)
  mean(sapply(other_cols, function(col) safe_cor(dat_window[[target_col]], dat_window[[col]])), na.rm = TRUE)
}

target_other_direction_agreement <- function(dat_window, pollutant, station = target_station,
                                             min_agree_other = 3, eps = 0) {
  target_col <- get_target_col(pollutant, station)
  other_cols <- setdiff(get_station_cols(dat_window, pollutant), target_col)
  if (!(target_col %in% names(dat_window)) || length(other_cols) < min_agree_other) return(NA_real_)
  d0 <- diff(dat_window[[target_col]])
  dm <- apply(dat_window[, other_cols, drop = FALSE], 2, diff)
  if (is.null(dim(dm))) dm <- matrix(dm, ncol = 1)
  agree <- rep(NA, length(d0))
  for (i in seq_along(d0)) {
    if (is.na(d0[i]) || abs(d0[i]) <= eps) next
    ds <- dm[i, ]
    if (sum(!is.na(ds)) < min_agree_other) next
    agree[i] <- if (d0[i] > eps) sum(ds > eps, na.rm = TRUE) >= min_agree_other else sum(ds < -eps, na.rm = TRUE) >= min_agree_other
  }
  mean(agree, na.rm = TRUE)
}

target_local_deviation_rate <- function(dat_window, pollutant, station = target_station, z_cut = 3.5) {
  target_col <- get_target_col(pollutant, station)
  other_cols <- setdiff(get_station_cols(dat_window, pollutant), target_col)
  if (!(target_col %in% names(dat_window)) || length(other_cols) < 3) return(NA_real_)
  x0 <- dat_window[[target_col]]
  xm <- as.matrix(dat_window[, other_cols, drop = FALSE])
  flag <- rep(FALSE, length(x0))
  for (i in seq_along(x0)) {
    others <- xm[i, ]
    med <- median(others, na.rm = TRUE)
    mad_val <- mad(others, constant = 1.4826, na.rm = TRUE)
    flag[i] <- if (is.na(mad_val) || mad_val == 0) FALSE else abs((x0[i] - med) / mad_val) > z_cut
  }
  mean(flag, na.rm = TRUE)
}

calc_spatial_metrics <- function(dat_window, station = target_station) {
  all_cols <- c(get_station_cols(dat_window, "PM10"), get_station_cols(dat_window, "PM25"))
  vals <- as.matrix(dat_window[, all_cols, drop = FALSE])
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
  spatial_score <- 0.20 * pm10_target_cor + 0.20 * pm25_target_cor +
    0.15 * pm10_all_cor + 0.15 * pm25_all_cor +
    0.15 * pm10_dir + 0.15 * pm25_dir -
    0.80 * missing_rate - 0.60 * deviation_rate - 1.00 * negative_rate
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
    spatial_score = spatial_score
  )
}

# ----------------------------
# 3. SNSeg 检测函数
# ----------------------------
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

run_snseg_multi_mean <- function(dat_window, cols) {
  x <- dat_window |> select(all_of(cols)) |> as.matrix()
  if (anyNA(x)) stop("run_snseg_multi_mean 输入含缺失值，请先过滤。")
  if (nrow(x) < 50) return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  x <- scale(x)
  res <- tryCatch(
    SNSeg_Multi(ts = x, paras_to_test = "mean", confidence = sn_confidence,
                grid_size_scale = sn_grid_size_scale, grid_size = NULL,
                plot_SN = FALSE, est_cp_loc = TRUE),
    error = function(e) { message("SNSeg_Multi 出错：", e$message); NULL }
  )
  if (is.null(res)) return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  cp <- extract_cp_from_snseg(res)
  list(ok = TRUE, n_cp = length(cp), cp_index = cp, cp_datetime = dat_window$datetime[cp])
}

run_snseg_uni_mean <- function(dat_window, col) {
  x <- dat_window[[col]]
  if (anyNA(x)) stop("run_snseg_uni_mean 输入含缺失值，请先过滤。")
  if (length(x) < 50) return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  x <- as.numeric(scale(x))
  if (!exists("SNSeg_Uni", mode = "function")) {
    stop("当前 SNSeg 包中未找到 SNSeg_Uni()。请检查本地 SNSeg 版本中的单变量函数名称。")
  }
  res <- tryCatch(
    SNSeg_Uni(ts = x, paras_to_test = "mean", confidence = sn_confidence,
              grid_size_scale = sn_grid_size_scale, grid_size = NULL,
              plot_SN = FALSE, est_cp_loc = TRUE),
    error = function(e) { message("SNSeg_Uni 出错：", e$message); NULL }
  )
  if (is.null(res)) return(list(ok = FALSE, n_cp = NA_integer_, cp_index = integer(0), cp_datetime = as.POSIXct(character(0))))
  cp <- extract_cp_from_snseg(res)
  list(ok = TRUE, n_cp = length(cp), cp_index = cp, cp_datetime = dat_window$datetime[cp])
}

# ----------------------------
# 4. 变点比较与评分
# ----------------------------
min_distance_to_set <- function(x, set) if (length(set) == 0) Inf else min(abs(x - set))

compare_multi_with_single <- function(cp_multi, cp_var1, cp_var2, tolerance = cp_match_tolerance) {
  if (length(cp_multi) == 0) {
    return(tibble(multi_cp_count = 0, multi_only_count = 0,
                  captured_by_var1_count = 0, captured_by_var2_count = 0,
                  captured_by_any_single_count = 0, captured_by_both_single_count = 0,
                  multi_only_ratio = NA_real_, captured_by_var1_ratio = NA_real_,
                  captured_by_var2_ratio = NA_real_, captured_by_any_single_ratio = NA_real_,
                  captured_by_both_single_ratio = NA_real_))
  }
  cap1 <- sapply(cp_multi, function(cp) min_distance_to_set(cp, cp_var1) <= tolerance)
  cap2 <- sapply(cp_multi, function(cp) min_distance_to_set(cp, cp_var2) <= tolerance)
  cap_any <- cap1 | cap2
  cap_both <- cap1 & cap2
  tibble(
    multi_cp_count = length(cp_multi),
    multi_only_count = sum(!cap_any),
    captured_by_var1_count = sum(cap1),
    captured_by_var2_count = sum(cap2),
    captured_by_any_single_count = sum(cap_any),
    captured_by_both_single_count = sum(cap_both),
    multi_only_ratio = sum(!cap_any) / length(cp_multi),
    captured_by_var1_ratio = sum(cap1) / length(cp_multi),
    captured_by_var2_ratio = sum(cap2) / length(cp_multi),
    captured_by_any_single_ratio = sum(cap_any) / length(cp_multi),
    captured_by_both_single_ratio = sum(cap_both) / length(cp_multi)
  )
}

segment_lengths_from_cp <- function(cp_index, n) diff(c(0, sort(unique(cp_index)), n))
cp_spacing_score <- function(cp_index, n, min_seg_len = min_segment_length) {
  if (length(cp_index) == 0) return(NA_real_)
  pmin(min(segment_lengths_from_cp(cp_index, n)) / min_seg_len, 1)
}

mean_shift_strength <- function(dat_window, cols, cp_index, min_seg_len = min_segment_length) {
  if (length(cp_index) == 0) return(NA_real_)
  x <- dat_window |> select(all_of(cols)) |> as.matrix() |> scale()
  n <- nrow(x)
  cuts <- c(0, sort(unique(cp_index)), n)
  if (any(diff(cuts) < min_seg_len)) return(NA_real_)
  means <- lapply(seq_len(length(cuts) - 1), function(j) colMeans(x[(cuts[j] + 1):cuts[j + 1], , drop = FALSE]))
  jumps <- sapply(seq_len(length(means) - 1), function(j) sqrt(sum((means[[j + 1]] - means[[j]])^2)))
  mean(jumps, na.rm = TRUE)
}

evaluate_window_cp <- function(dat_window, station = target_station) {
  var1_col <- paste0(var1_name, "_", station)
  var2_col <- paste0(var2_name, "_", station)
  res_var1 <- run_snseg_uni_mean(dat_window, var1_col)
  res_var2 <- run_snseg_uni_mean(dat_window, var2_col)
  res_multi <- run_snseg_multi_mean(dat_window, c(var1_col, var2_col))
  cmp <- compare_multi_with_single(res_multi$cp_index, res_var1$cp_index, res_var2$cp_index)
  strength_multi <- mean_shift_strength(dat_window, c(var1_col, var2_col), res_multi$cp_index)
  strength_var1 <- mean_shift_strength(dat_window, c(var1_col), res_var1$cp_index)
  strength_var2 <- mean_shift_strength(dat_window, c(var2_col), res_var2$cp_index)
  spacing <- cp_spacing_score(res_multi$cp_index, nrow(dat_window))
  clean_structure_score <- ifelse(is.na(strength_multi) || length(res_multi$cp_index) == 0, NA_real_,
                                  0.70 * pmin(strength_multi / 1.0, 1) + 0.30 * ifelse(is.na(spacing), 0, spacing))
  dual_value_score <- ifelse(is.na(clean_structure_score), NA_real_,
                             0.30 * clean_structure_score +
                               0.25 * ifelse(is.na(cmp$multi_only_ratio), 0, cmp$multi_only_ratio) +
                               0.25 * (1 - ifelse(is.na(cmp$captured_by_var1_ratio), 0, cmp$captured_by_var1_ratio)) +
                               0.10 * (1 - ifelse(is.na(cmp$captured_by_both_single_ratio), 0, cmp$captured_by_both_single_ratio)) +
                               0.10 * ifelse(!is.na(strength_multi) && strength_multi < 1.0, 1, 0))
  tibble(
    n_cp_var1 = res_var1$n_cp,
    cp_index_var1 = paste(res_var1$cp_index, collapse = ";"),
    cp_datetime_var1 = paste(format(res_var1$cp_datetime, "%Y-%m-%d %H:%M"), collapse = ";"),
    n_cp_var2 = res_var2$n_cp,
    cp_index_var2 = paste(res_var2$cp_index, collapse = ";"),
    cp_datetime_var2 = paste(format(res_var2$cp_datetime, "%Y-%m-%d %H:%M"), collapse = ";"),
    n_cp_multi = res_multi$n_cp,
    cp_index_multi = paste(res_multi$cp_index, collapse = ";"),
    cp_datetime_multi = paste(format(res_multi$cp_datetime, "%Y-%m-%d %H:%M"), collapse = ";"),
    strength_var1 = strength_var1,
    strength_var2 = strength_var2,
    strength_multi = strength_multi,
    spacing_score_multi = spacing,
    clean_structure_score = clean_structure_score,
    dual_value_score = dual_value_score
  ) |> bind_cols(cmp)
}

# ----------------------------
# 5. 候选窗口筛选
# ----------------------------
pre_screen_spatial_windows <- function(wide_dat,
                                       max_negative_rate = 0,
                                       min_target_cor = 0.50,
                                       min_direction = 0.45,
                                       max_deviation_rate = 0.12) {
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
      pass <- with(metrics,
                   missing_rate == 0 &&
                     negative_rate <= max_negative_rate &&
                     pm10_target_other_cor >= min_target_cor &&
                     pm25_target_other_cor >= min_target_cor &&
                     pm10_direction_agreement >= min_direction &&
                     pm25_direction_agreement >= min_direction &&
                     deviation_rate <= max_deviation_rate)
      result_list[[idx]] <- metrics |>
        mutate(window_hours = L, start_index = s, end_index = e,
               start_time = dat_window$datetime[1], end_time = dat_window$datetime[nrow(dat_window)],
               pass_spatial = pass)
      idx <- idx + 1
    }
  }
  if (length(result_list) == 0) return(tibble())
  bind_rows(result_list) |> arrange(desc(pass_spatial), desc(spatial_score))
}

find_candidate_windows <- function(wide_dat, min_multi_cp = 1, max_multi_cp = 8, min_clean_score = 0.25) {
  cat("\n============================================================\n开始筛选候选时间段\n============================================================\n")
  spatial_res <- pre_screen_spatial_windows(wide_dat)
  write.csv(spatial_res, file.path(out_dir, "01_spatial_screen_all_complete_windows.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  spatial_keep <- spatial_res |>
    filter(pass_spatial) |>
    group_by(window_hours) |>
    slice_max(order_by = spatial_score, n = max_spatial_candidates_per_length, with_ties = FALSE) |>
    ungroup() |>
    arrange(desc(spatial_score))
  write.csv(spatial_keep, file.path(out_dir, "02_spatial_screen_kept_for_snseg.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat("进入 SNSeg 检测的窗口数：", nrow(spatial_keep), "\n")
  if (nrow(spatial_keep) == 0) return(list(all = tibble(), good = tibble(), type_A = tibble(), type_B = tibble(), type_C = tibble()))
  res_list <- list()
  for (i in seq_len(nrow(spatial_keep))) {
    dat_window <- wide_dat[spatial_keep$start_index[i]:spatial_keep$end_index[i], ]
    cp_res <- evaluate_window_cp(dat_window)
    pass_cp <- with(cp_res, !is.na(n_cp_multi) && n_cp_multi >= min_multi_cp && n_cp_multi <= max_multi_cp &&
                      !is.na(clean_structure_score) && clean_structure_score >= min_clean_score)
    res_list[[i]] <- spatial_keep[i, ] |> bind_cols(cp_res) |> mutate(pass_cp = pass_cp)
    if (i %% 20 == 0) cat("已完成 SNSeg 检测：", i, "/", nrow(spatial_keep), "\n")
  }
  all_res <- bind_rows(res_list) |>
    mutate(final_score = 0.35 * spatial_score +
             0.30 * ifelse(is.na(dual_value_score), 0, dual_value_score) +
             0.25 * ifelse(is.na(clean_structure_score), 0, clean_structure_score) +
             0.10 * pmin(ifelse(is.na(strength_multi), 0, strength_multi) / 1.0, 1)) |>
    arrange(desc(pass_cp), desc(final_score), desc(dual_value_score), desc(spatial_score))
  good <- all_res |> filter(pass_cp) |> arrange(desc(final_score), desc(dual_value_score), desc(spatial_score))
  type_A <- good |> filter(n_cp_multi >= 1, captured_by_var1_ratio < 1) |>
    arrange(desc(final_score), desc(1 - captured_by_var1_ratio), desc(dual_value_score))
  type_B <- good |> filter(n_cp_multi >= 1, captured_by_var1_ratio < 1, captured_by_var2_ratio < 1) |>
    arrange(desc(final_score), desc(multi_only_ratio), desc(dual_value_score))
  type_C <- good |>
    mutate(single_cp_count_gap = abs(ifelse(is.na(n_cp_var1), 0, n_cp_var1) - ifelse(is.na(n_cp_var2), 0, n_cp_var2))) |>
    filter(n_cp_multi >= 1, single_cp_count_gap >= 1, clean_structure_score >= min_clean_score) |>
    arrange(desc(final_score), desc(single_cp_count_gap), desc(dual_value_score))
  write.csv(all_res, file.path(out_dir, "03_all_candidate_windows_with_cp.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(good, file.path(out_dir, "04_good_candidate_windows.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(type_A, file.path(out_dir, "05_type_A_var1_misses_multi_cps.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(type_B, file.path(out_dir, "06_type_B_weak_joint_signal_candidates.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(type_C, file.path(out_dir, "07_type_C_single_results_inconsistent.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  list(all = all_res, good = good, type_A = type_A, type_B = type_B, type_C = type_C)
}

# ----------------------------
# 6. 可视化：三行子图
# ----------------------------
parse_cp_times <- function(cp_text) {
  if (is.na(cp_text) || cp_text == "") return(as.POSIXct(character(0)))
  ymd_hm(str_split(cp_text, ";")[[1]], tz = "Asia/Shanghai")
}

plot_three_panel_detection <- function(candidate_row, wide_dat, station = target_station, prefix = "candidate") {
  dat_window <- wide_dat[candidate_row$start_index:candidate_row$end_index, ]
  var1_col <- paste0(var1_name, "_", station)
  var2_col <- paste0(var2_name, "_", station)
  cp_var1 <- parse_cp_times(candidate_row$cp_datetime_var1)
  cp_var2 <- parse_cp_times(candidate_row$cp_datetime_var2)
  cp_multi <- parse_cp_times(candidate_row$cp_datetime_multi)
  p1_dat <- dat_window |> transmute(datetime = datetime, value = .data[[var1_col]], panel = "① 单变量方法检测变量1")
  p2_dat <- dat_window |> transmute(datetime = datetime, value = .data[[var2_col]], panel = "② 单变量方法检测变量2")
  x1 <- as.numeric(scale(dat_window[[var1_col]]))
  x2 <- as.numeric(scale(dat_window[[var2_col]]))
  p3_dat <- tibble(datetime = rep(dat_window$datetime, 2), value = c(x1, x2),
                   variable = rep(c("变量1（标准化）", "变量2（标准化）"), each = nrow(dat_window)),
                   panel = "③ 双变量方法检测二维联合序列")
  p12_dat <- bind_rows(p1_dat |> mutate(variable = "变量1"), p2_dat |> mutate(variable = "变量2"))
  p <- ggplot() +
    geom_line(data = p12_dat, aes(x = datetime, y = value), linewidth = 0.65) +
    geom_line(data = p3_dat, aes(x = datetime, y = value, linetype = variable), linewidth = 0.65) +
    geom_vline(data = tibble(x = cp_var1, panel = "① 单变量方法检测变量1"), aes(xintercept = as.numeric(x)), linetype = "dashed", linewidth = 0.7) +
    geom_vline(data = tibble(x = cp_var2, panel = "② 单变量方法检测变量2"), aes(xintercept = as.numeric(x)), linetype = "dashed", linewidth = 0.7) +
    geom_vline(data = tibble(x = cp_multi, panel = "③ 双变量方法检测二维联合序列"), aes(xintercept = as.numeric(x)), linetype = "solid", linewidth = 0.75) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(title = paste0("单变量与双变量变点检测结果对比：", format(dat_window$datetime[1], "%Y-%m-%d %H:%M"), " 至 ", format(dat_window$datetime[nrow(dat_window)], "%Y-%m-%d %H:%M")),
         subtitle = "第1行和第2行虚线为对应单变量检测变点；第3行实线为双变量检测变点。",
         x = "时间", y = "观测值 / 标准化观测值", linetype = NULL) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
  print(p)
  fname <- paste0(prefix, "_three_panel_", format(dat_window$datetime[1], "%Y%m%d%H"), "_", format(dat_window$datetime[nrow(dat_window)], "%Y%m%d%H"), ".png")
  ggsave(file.path(out_dir, fname), p, width = 12, height = 8, dpi = 300)
  invisible(p)
}

# ----------------------------
# 7. 执行筛选与画图
# ----------------------------
screen_res <- find_candidate_windows(wide_dat)
all_candidates <- screen_res$all
good_candidates <- screen_res$good
type_A <- screen_res$type_A
type_B <- screen_res$type_B
type_C <- screen_res$type_C

cat("\n============================================================\n")
cat("筛选完成。输出文件夹：", out_dir, "\n")
cat("通过筛选的候选窗口数：", nrow(good_candidates), "\n")
cat("类型 A 窗口数：", nrow(type_A), "\n")
cat("类型 B 窗口数：", nrow(type_B), "\n")
cat("类型 C 窗口数：", nrow(type_C), "\n")
cat("============================================================\n")

cat("\n总体候选前10个：\n"); print(head(good_candidates, 10))
cat("\n类型A前10个：双变量变点未被变量1单变量完全捕捉。\n"); print(head(type_A, 10))
cat("\n类型B前10个：两个单变量对双变量变点捕捉都不充分。\n"); print(head(type_B, 10))
cat("\n类型C前10个：两个单变量结果不一致。\n"); print(head(type_C, 10))

if (nrow(type_A) > 0) {
  for (i in seq_len(min(5, nrow(type_A)))) plot_three_panel_detection(type_A[i, ], wide_dat, prefix = paste0("type_A_top_", i))
}
if (nrow(type_B) > 0) {
  for (i in seq_len(min(5, nrow(type_B)))) plot_three_panel_detection(type_B[i, ], wide_dat, prefix = paste0("type_B_top_", i))
}
if (nrow(good_candidates) > 0) {
  for (i in seq_len(min(5, nrow(good_candidates)))) plot_three_panel_detection(good_candidates[i, ], wide_dat, prefix = paste0("overall_top_", i))
}

cat("\n结果文件：\n")
cat("01_spatial_screen_all_complete_windows.csv：所有无缺失滑动窗口的空间同步性指标。\n")
cat("02_spatial_screen_kept_for_snseg.csv：进入SNSeg检测的空间同步窗口。\n")
cat("03_all_candidate_windows_with_cp.csv：所有进入检测窗口的变点结果。\n")
cat("04_good_candidate_windows.csv：空间同步性好且双变量变点较清晰的窗口。\n")
cat("05_type_A_var1_misses_multi_cps.csv：类型A窗口。\n")
cat("06_type_B_weak_joint_signal_candidates.csv：类型B窗口。\n")
cat("07_type_C_single_results_inconsistent.csv：类型C窗口。\n")
cat("图片：*_three_panel_*.png，三行子图。\n")
cat("============================================================\n")
