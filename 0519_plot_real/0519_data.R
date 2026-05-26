# ============================================================
# 实际数据双变量异常识别：SNSeg_Multi + 反距离加权空间差值 + CBB Hotelling 检验
# 锚点式完整流程版：L = 10；先用原始数据确定第2、第3个候选变点位置，
# 情形2和情形3均在这些“最开始的候选变点位置”左侧局部窗口内施加人工异常，
# 然后对处理后的数据重新执行第一阶段变点检测与第二阶段异常识别。
#
# 三种情形：
#   1）原始情形：不做人工处理，直接进行异常识别；
#   2）单人工异常变点：在原始第2个候选变点左侧，对目标站点东四施加双变量小幅变化；
#   3）双人工异常变点：在原始第2、第3个候选变点左侧，对目标站点东四分别施加小幅变化。
#
# 说明：
#   - 第一阶段：对目标站点 PM10 与 PM2.5 二维序列做 Z-score 标准化后，
#               使用 SNSeg_Multi 进行均值变点检测。
#   - 第二阶段：使用原始浓度序列构造“目标站点 - 临近站点反距离加权值”的空间差值序列，
#               对第一阶段候选变点逐个进行正则化 Hotelling 型统计量 + CBB 检验。
#   - 分类规则：p_value < alpha 判为“异常”，否则判为“正常”。
#   - 绘图：正常变点用绿色竖线，异常变点用红色竖线。
# ============================================================

rm(list = ls())
options(scipen = 999, digits = 4)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
  library(SNSeg)
  library(Rcpp)
  library(RcppArmadillo)
})

set.seed(42)

# ============================================================
# 0. 用户参数设置
# ============================================================

data_dir <- "D:/03/2025/Beijing"
target_station <- "Dongsi"

# 这里改成你筛选好的“正常基础时间段”
start_time <- as.POSIXct("2024-12-07 00:00:00", tz = "Asia/Shanghai")
end_time   <- as.POSIXct("2024-12-16 23:00:00", tz = "Asia/Shanghai")

# 北京 6 个站点
selected_stations <- c(
  "Dongsi",
  "Tiantan",
  "Guanyuan",
  "Wanshouxigong",
  "Nongzhanguan",
  "Aotizhongxin"
)

# 第一阶段 SNSeg 参数：与 lmy_0515 中的设置保持一致
sn_confidence <- 0.95
sn_grid_size_scale <- 0.05

# 第二阶段异常识别参数：以 0518 代码为主
# 最小修改：只将 L 从 40 改为 10
L <- 10
M <- 1000
lambda_frac <- 0.15
alpha <- 0.05

choose_cbb_block_size <- function(L) {
  max(5L, min(as.integer(L), as.integer(round(sqrt(2 * L)))))
}
block_size <- choose_cbb_block_size(L)

# 候选变点与人工变点匹配时，仅用于结果表标注，不影响分类
cp_match_tolerance <- 6

# 是否显示曲线图例
show_series_legend <- TRUE

# 绘图参数
font_family <- "SimHei"     # 可改为 "Microsoft YaHei"
base_size <- 14
title_size <- 18
axis_title_size <- 14
axis_text_size <- 12
legend_text_size <- 12
x_date_breaks <- "1 day"

# 输出文件夹
out_dir <- file.path(data_dir, "actual_stage2_anomaly_identification_L5_anchor_cp2_cp3_rerun")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ============================================================
# 1. 站点经纬度与反距离权重
# ============================================================

# 经纬度来自 Beijing Multi-Site Air-Quality 常用站点信息。
# lon 为经度，lat 为纬度。
station_coord <- tibble::tribble(
  ~station,         ~station_zh,   ~lon,     ~lat,
  "Dongsi",         "东四",        116.417,  39.929,
  "Tiantan",        "天坛",        116.407,  39.886,
  "Guanyuan",       "官园",        116.339,  39.929,
  "Wanshouxigong",  "万寿西宫",    116.352,  39.878,
  "Nongzhanguan",   "农展馆",      116.461,  39.937,
  "Aotizhongxin",   "奥体中心",    116.397,  39.982,
  "Changping",      "昌平",        116.231,  40.217,
  "Dingling",       "定陵",        116.220,  40.292,
  "Gucheng",        "古城",        116.184,  39.914,
  "Huairou",        "怀柔",        116.628,  40.328,
  "Shunyi",         "顺义",        116.655,  40.127,
  "Wanliu",         "万柳",        116.287,  39.987
)

target_station_zh <- station_coord$station_zh[match(target_station, station_coord$station)]
if (is.na(target_station_zh)) target_station_zh <- target_station

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371
  to_rad <- pi / 180
  dlon <- (lon2 - lon1) * to_rad
  dlat <- (lat2 - lat1) * to_rad
  a <- sin(dlat / 2)^2 +
    cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * R * atan2(sqrt(a), sqrt(1 - a))
}

calc_idw_weights <- function(station_coord, target_station, selected_stations) {
  target_row <- station_coord |> dplyr::filter(station == target_station)
  if (nrow(target_row) != 1) stop("目标站点经纬度缺失或不唯一。")

  ref_df <- station_coord |>
    dplyr::filter(station %in% setdiff(selected_stations, target_station)) |>
    dplyr::mutate(
      distance_km = haversine_km(
        lon1 = target_row$lon,
        lat1 = target_row$lat,
        lon2 = lon,
        lat2 = lat
      ),
      inv_distance = 1 / pmax(distance_km, 1e-6),
      weight = inv_distance / sum(inv_distance)
    )

  if (nrow(ref_df) == 0) stop("没有可用参考站点。")

  ref_df |>
    dplyr::select(station, station_zh, lon, lat, distance_km, weight) |>
    dplyr::arrange(desc(weight))
}

idw_weights <- calc_idw_weights(
  station_coord = station_coord,
  target_station = target_station,
  selected_stations = selected_stations
)

write.csv(
  idw_weights,
  file.path(out_dir, "01_IDW_weights.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n反距离加权权重：\n")
print(idw_weights)


# ============================================================
# 2. 读取与整理实际数据
# ============================================================

extract_station_name <- function(file) {
  raw_name <- basename(file) |>
    stringr::str_remove("\\.csv$") |>
    stringr::str_remove("^Beijing_")

  # 优先精确匹配；若文件名带日期后缀，则用包含关系识别站点名。
  if (raw_name %in% station_coord$station) return(raw_name)

  hits <- station_coord$station[
    stringr::str_detect(raw_name, paste0("(^|_)", station_coord$station, "($|_)"))
  ]

  if (length(hits) >= 1) return(hits[1])
  raw_name
}

read_one_station <- function(file) {
  station_name <- extract_station_name(file)

  dat <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  need_vars <- c("year", "month", "day", "hour", "PM2.5", "PM10")
  miss_vars <- setdiff(need_vars, names(dat))

  if (length(miss_vars) > 0) {
    stop(
      paste0(
        "文件 ", basename(file), " 缺少变量：",
        paste(miss_vars, collapse = ", ")
      )
    )
  }

  dat |>
    dplyr::transmute(
      station = station_name,
      datetime = lubridate::make_datetime(
        year = year,
        month = month,
        day = day,
        hour = hour,
        tz = "Asia/Shanghai"
      ),
      PM25 = as.numeric(`PM2.5`),
      PM10 = as.numeric(PM10)
    )
}

files <- list.files(
  data_dir,
  pattern = "^Beijing_.*\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("未找到 Beijing_*.csv 文件，请检查 data_dir。")
}

raw_long <- purrr::map_dfr(files, read_one_station) |>
  dplyr::filter(station %in% selected_stations)

miss_station <- setdiff(selected_stations, unique(raw_long$station))
if (length(miss_station) > 0) {
  stop(
    paste0(
      "以下站点没有读取到数据：",
      paste(miss_station, collapse = ", "),
      "。请检查文件名或 selected_stations。"
    )
  )
}

wide_dat_all <- raw_long |>
  dplyr::filter(datetime >= start_time, datetime <= end_time) |>
  tidyr::pivot_wider(
    id_cols = datetime,
    names_from = station,
    values_from = c(PM10, PM25),
    names_glue = "{.value}_{station}"
  ) |>
  dplyr::arrange(datetime)

if (nrow(wide_dat_all) == 0) {
  stop("当前时间段没有数据，请检查 start_time 和 end_time。")
}

# 补齐小时序列；若有缺失小时或缺失浓度，后面直接停止，避免影响变点检测与异常识别。
full_time <- tibble::tibble(
  datetime = seq(
    from = min(wide_dat_all$datetime, na.rm = TRUE),
    to   = max(wide_dat_all$datetime, na.rm = TRUE),
    by   = "hour"
  )
)

wide_dat <- full_time |>
  dplyr::left_join(wide_dat_all, by = "datetime") |>
  dplyr::arrange(datetime)

required_cols <- as.vector(outer(c("PM10", "PM25"), selected_stations, paste, sep = "_"))

miss_cols <- setdiff(required_cols, names(wide_dat))
if (length(miss_cols) > 0) {
  stop(paste0("宽格式数据缺少列：", paste(miss_cols, collapse = ", ")))
}

if (anyNA(wide_dat[, required_cols, drop = FALSE])) {
  miss_rows <- wide_dat |>
    dplyr::filter(if_any(dplyr::all_of(required_cols), is.na)) |>
    dplyr::select(datetime, dplyr::all_of(required_cols))

  print(head(miss_rows, 20))
  stop("当前时间段存在缺失小时或污染物缺失值，请更换时间段或先处理缺失。")
}

n_time <- nrow(wide_dat)
if (n_time < 2 * L + 20) {
  stop("当前时间段过短，无法进行局部窗口异常识别。请扩大时间范围或减小 L。")
}

cat("\n数据读取完成：", n_time, "个小时。\n")
cat("时间范围：", format(min(wide_dat$datetime)), "至", format(max(wide_dat$datetime)), "\n")


# ============================================================
# 3. Rcpp：CBB + 正则化 Hotelling 型 H 统计量
# ============================================================

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

static inline mat estimate_reg_cov(const mat& Z_centered, double lambda_frac) {
  mat S = cov(Z_centered);
  int d = S.n_cols;

  double trS = trace(S);
  double lambda = lambda_frac * trS / std::max(1, d);

  if (!std::isfinite(lambda) || lambda <= 1e-10) {
    lambda = 1e-6;
  }

  S += lambda * eye<mat>(d, d);
  return S;
}

static inline double calc_hotelling_stat(const mat& Z, int L, const mat& Sinv) {
  rowvec mean_before = mean(Z.rows(0, L - 1), 0);
  rowvec mean_after  = mean(Z.rows(L, 2 * L - 1), 0);

  rowvec diff = mean_after - mean_before;
  mat val = diff * Sinv * diff.t();

  return (static_cast<double>(L) / 2.0) * val(0, 0);
}

// [[Rcpp::export]]
List block_boot_cbb_hotelling_cpp(
    const arma::mat& Z,
    int L,
    int M = 500,
    int block_size = 10,
    double lambda_frac = 0.05
) {
  int n = Z.n_rows;
  int d = Z.n_cols;

  if (n != 2 * L) stop("Z.n_rows must be exactly 2 * L.");
  if (L < 5 || d < 1 || M < 10) stop("Invalid inputs.");

  rowvec mean_before = mean(Z.rows(0, L - 1), 0);
  rowvec mean_after  = mean(Z.rows(L, 2 * L - 1), 0);

  mat Z_tilde = Z;
  for (int t = 0; t < L; ++t) {
    Z_tilde.row(t) -= mean_before;
  }
  for (int t = L; t < n; ++t) {
    Z_tilde.row(t) -= mean_after;
  }

  mat Sigma_reg = estimate_reg_cov(Z_tilde, lambda_frac);

  mat Sinv;
  bool ok = inv_sympd(Sinv, Sigma_reg);
  if (!ok) {
    Sinv = pinv(Sigma_reg);
  }

  double H_obs = calc_hotelling_stat(Z, L, Sinv);

  int b = std::max(1, block_size);
  b = std::min(b, n);

  int exceed_count = 0;
  vec H_boot(M, fill::zeros);

  for (int m = 0; m < M; ++m) {
    mat Z_star(n, d, fill::zeros);
    int t = 0;

    while (t < n) {
      int start_idx = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(n))));

      for (int i = 0; i < b && t < n; ++i) {
        int idx = (start_idx + i) % n;
        Z_star.row(t) = Z_tilde.row(idx);
        t++;
      }
    }

    double H_star = calc_hotelling_stat(Z_star, L, Sinv);
    H_boot(m) = H_star;

    if (H_star >= H_obs) {
      exceed_count++;
    }
  }

  double p_value = static_cast<double>(exceed_count + 1) / static_cast<double>(M + 1);

  vec H_sorted = sort(H_boot);
  int idx95 = std::min(M - 1, std::max(0, static_cast<int>(std::ceil(0.95 * M)) - 1));

  return List::create(
    _["p_value"] = p_value,
    _["H_obs"] = H_obs,
    _["H_boot_q95"] = H_sorted(idx95)
  );
}
'

Rcpp::sourceCpp(code = cpp_code)


# ============================================================
# 4. 第一阶段：SNSeg_Multi 双变量均值变点检测
# ============================================================

safe_scale_matrix <- function(x) {
  if (anyNA(x)) stop("输入矩阵含缺失值，不能标准化。")

  sds <- apply(x, 2, sd)
  if (any(!is.finite(sds) | sds <= 1e-12)) {
    bad <- colnames(x)[which(!is.finite(sds) | sds <= 1e-12)]
    stop(paste0("以下变量标准差为 0 或不可用：", paste(bad, collapse = ", ")))
  }

  scale(x)
}

extract_cp_from_snseg <- function(res, n) {
  if (is.null(res)) return(integer(0))

  possible_names <- c(
    "est_cp", "cp", "cpts", "changepoints",
    "change_points", "cp_index"
  )

  for (nm in possible_names) {
    if (nm %in% names(res)) {
      cp <- res[[nm]]
      if (length(cp) == 0 || all(is.na(cp))) return(integer(0))

      cp <- as.integer(round(as.numeric(unlist(cp))))
      cp <- cp[!is.na(cp)]
      cp <- cp[cp > 0 & cp < n]
      return(sort(unique(cp)))
    }
  }

  message("未从 SNSeg 返回对象中识别到变点字段。返回对象字段名为：")
  print(names(res))
  integer(0)
}

get_target_matrix <- function(wide_dat, target_station) {
  x <- as.matrix(wide_dat[, c(
    paste0("PM10_", target_station),
    paste0("PM25_", target_station)
  )])
  colnames(x) <- c("PM10", "PM2.5")
  x
}

detect_candidate_cps <- function(
    wide_dat,
    target_station,
    confidence = 0.95,
    grid_size_scale = 0.05
) {
  x_raw <- get_target_matrix(wide_dat, target_station)
  x_std <- safe_scale_matrix(x_raw)

  res <- tryCatch(
    SNSeg_Multi(
      ts = x_std,
      paras_to_test = "mean",
      confidence = confidence,
      grid_size_scale = grid_size_scale,
      grid_size = NULL,
      plot_SN = FALSE,
      est_cp_loc = TRUE
    ),
    error = function(e) {
      message("SNSeg_Multi 出错：", e$message)
      NULL
    }
  )

  cp_index <- extract_cp_from_snseg(res, nrow(wide_dat))

  list(
    result = res,
    cp_index = cp_index,
    cp_datetime = wide_dat$datetime[cp_index]
  )
}


# ============================================================
# 5. 第二阶段：反距离加权空间差值 + Hotelling CBB 异常识别
# ============================================================

get_station_matrix_from_wide <- function(wide_dat, station) {
  x <- as.matrix(wide_dat[, c(
    paste0("PM10_", station),
    paste0("PM25_", station)
  )])
  colnames(x) <- c("PM10", "PM2.5")
  x
}

calc_weighted_reference <- function(wide_dat, target_station, idw_weights) {
  ref_stations <- idw_weights$station
  weights <- idw_weights$weight

  ref_mat <- matrix(
    0,
    nrow = nrow(wide_dat),
    ncol = 2
  )
  colnames(ref_mat) <- c("PM10", "PM2.5")

  for (j in seq_along(ref_stations)) {
    ref_mat <- ref_mat + weights[j] * get_station_matrix_from_wide(wide_dat, ref_stations[j])
  }

  ref_mat
}

build_residual_window_raw <- function(
    wide_dat,
    tau,
    L,
    target_station,
    idw_weights
) {
  n_local <- nrow(wide_dat)

  s_idx <- tau - L + 1
  e_idx <- tau + L

  if (s_idx < 1 || e_idx > n_local) {
    return(NULL)
  }

  win_dat <- wide_dat[s_idx:e_idx, , drop = FALSE]

  target_mat <- get_station_matrix_from_wide(win_dat, target_station)
  ref_mat <- calc_weighted_reference(
    wide_dat = win_dat,
    target_station = target_station,
    idw_weights = idw_weights
  )

  Z_raw <- target_mat - ref_mat
  colnames(Z_raw) <- c("PM10差值", "PM2.5差值")

  list(
    Z_raw = Z_raw,
    s_idx = s_idx,
    e_idx = e_idx
  )
}

test_one_candidate <- function(
    wide_dat,
    tau,
    L,
    M,
    block_size,
    lambda_frac,
    target_station,
    idw_weights
) {
  feat <- build_residual_window_raw(
    wide_dat = wide_dat,
    tau = tau,
    L = L,
    target_station = target_station,
    idw_weights = idw_weights
  )

  if (is.null(feat)) return(NULL)

  out <- block_boot_cbb_hotelling_cpp(
    Z = feat$Z_raw,
    L = L,
    M = M,
    block_size = block_size,
    lambda_frac = lambda_frac
  )

  data.frame(
    tau = tau,
    datetime = wide_dat$datetime[tau],
    H_obs = as.numeric(out$H_obs),
    H_boot_q95 = as.numeric(out$H_boot_q95),
    p_value = as.numeric(out$p_value),
    stringsAsFactors = FALSE
  )
}

classify_candidates <- function(
    wide_dat,
    candidate_cps,
    L,
    M,
    block_size,
    lambda_frac,
    alpha,
    target_station,
    idw_weights
) {
  if (length(candidate_cps) == 0) {
    return(tibble::tibble(
      tau = integer(0),
      datetime = as.POSIXct(character(0), tz = "Asia/Shanghai"),
      H_obs = numeric(0),
      H_boot_q95 = numeric(0),
      p_value = numeric(0),
      pred_label = character(0)
    ))
  }

  res_list <- lapply(candidate_cps, function(tau) {
    test_one_candidate(
      wide_dat = wide_dat,
      tau = tau,
      L = L,
      M = M,
      block_size = block_size,
      lambda_frac = lambda_frac,
      target_station = target_station,
      idw_weights = idw_weights
    )
  })

  res_df <- dplyr::bind_rows(res_list)

  if (nrow(res_df) == 0) {
    return(tibble::tibble(
      tau = integer(0),
      datetime = as.POSIXct(character(0), tz = "Asia/Shanghai"),
      H_obs = numeric(0),
      H_boot_q95 = numeric(0),
      p_value = numeric(0),
      pred_label = character(0)
    ))
  }

  res_df |>
    dplyr::mutate(
      pred_label = ifelse(p_value < alpha, "异常", "正常")
    ) |>
    dplyr::arrange(tau)
}


# ============================================================
# 6. 人工均值变点设置与施加
# ============================================================

clamp_cp_index <- function(idx, n, L) {
  idx <- as.integer(round(idx))
  idx <- max(L + 2L, idx)
  idx <- min(n - L - 2L, idx)
  idx
}

nearest_time_index <- function(datetime_vec, target_time) {
  which.min(abs(as.numeric(datetime_vec - target_time)))
}

calc_auto_shifts <- function(wide_dat, target_station) {
  pm10 <- wide_dat[[paste0("PM10_", target_station)]]
  pm25 <- wide_dat[[paste0("PM25_", target_station)]]

  pm10_shift <- 0.30 * sd(pm10, na.rm = TRUE)
  pm25_shift <- 0.30 * sd(pm25, na.rm = TRUE)

  # 防止变化过小或过大。你可以根据图形效果继续调小或调大。
  pm10_shift <- round(max(2.0, min(8.0, pm10_shift)), 2)
  pm25_shift <- round(max(1.5, min(6.0, pm25_shift)), 2)

  list(
    both_vars = c(PM10 = pm10_shift, PM25 = -pm25_shift),
    pm10_only = c(PM10 = -pm10_shift, PM25 = 0)
  )
}

auto_shifts <- calc_auto_shifts(wide_dat, target_station)

cat("\n自动设置的人工均值变化幅度：\n")
print(auto_shifts)

# 你也可以直接手动指定，例如：
 auto_shifts$both_vars <- c(PM10 = -6, PM25 = -4)
 auto_shifts$pm10_only <- c(PM10 = -10, PM25 = 0)

# 先在原始数据上运行一次第一阶段检测，用于确定人工变化锚点。
# 这里优先选择原始第一阶段检测到的第 2 个和第 3 个均值变点，避免使用较不稳定的第 4 个变点。
base_stage1_res_for_anchor <- detect_candidate_cps(
  wide_dat = wide_dat,
  target_station = target_station,
  confidence = sn_confidence,
  grid_size_scale = sn_grid_size_scale
)

base_candidate_cps <- base_stage1_res_for_anchor$cp_index
base_candidate_times <- base_stage1_res_for_anchor$cp_datetime

cat("\n原始数据第一阶段候选均值变点：\n")
if (length(base_candidate_cps) == 0) {
  cat("  未检测到候选变点。后续人工锚点将使用比例位置兜底。\n")
} else {
  print(data.frame(
    序号 = seq_along(base_candidate_cps),
    位置 = base_candidate_cps,
    时间 = format(base_candidate_times, "%Y-%m-%d %H:%M")
  ))
}

# 若你想手动指定人工变点时间，可在这里填写 POSIXct 向量；保持 NULL 则按候选变点序号自动选择。
# 当前默认：情形2使用第 2 个候选变点；情形3使用第 2 个和第 3 个候选变点。
manual_cp_time_scenario2 <- NULL
manual_cp_time_scenario3 <- NULL

manual_cp_order_scenario2 <- 2L
manual_cp_order_scenario3 <- c(2L, 3L)

choose_anchor_cp <- function(candidate_cps, order_no, fallback_ratio, n_time, L) {
  if (length(candidate_cps) >= order_no && !is.na(candidate_cps[order_no])) {
    return(clamp_cp_index(candidate_cps[order_no], n_time, L))
  }

  warning(
    paste0(
      "原始第一阶段候选变点数量不足，无法选择第 ", order_no,
      " 个候选变点；改用比例位置 ", fallback_ratio, " 作为兜底。"
    )
  )

  clamp_cp_index(round(fallback_ratio * n_time), n_time, L)
}

if (is.null(manual_cp_time_scenario2)) {
  cp_s2 <- choose_anchor_cp(
    candidate_cps = base_candidate_cps,
    order_no = manual_cp_order_scenario2,
    fallback_ratio = 0.50,
    n_time = n_time,
    L = L
  )
} else {
  cp_s2 <- clamp_cp_index(nearest_time_index(wide_dat$datetime, manual_cp_time_scenario2[1]), n_time, L)
}

if (is.null(manual_cp_time_scenario3)) {
  cp_s3_1 <- choose_anchor_cp(
    candidate_cps = base_candidate_cps,
    order_no = manual_cp_order_scenario3[1],
    fallback_ratio = 0.35,
    n_time = n_time,
    L = L
  )
  cp_s3_2 <- choose_anchor_cp(
    candidate_cps = base_candidate_cps,
    order_no = manual_cp_order_scenario3[2],
    fallback_ratio = 0.65,
    n_time = n_time,
    L = L
  )
} else {
  cp_s3_1 <- clamp_cp_index(nearest_time_index(wide_dat$datetime, manual_cp_time_scenario3[1]), n_time, L)
  cp_s3_2 <- clamp_cp_index(nearest_time_index(wide_dat$datetime, manual_cp_time_scenario3[2]), n_time, L)
}

if (cp_s3_2 <= cp_s3_1 + 2 * L) {
  warning("情形3两个人工变点距离较近，可能影响局部窗口检验。建议调整 manual_cp_time_scenario3 或 manual_cp_order_scenario3。")
}

# 保存原始第一阶段候选变点和本次人工变化锚点，便于结果核对。
base_anchor_table <- tibble::tibble(
  anchor_role = c("情形2锚点", "情形3第一个锚点", "情形3第二个锚点"),
  cp_index = c(cp_s2, cp_s3_1, cp_s3_2),
  cp_time = format(wide_dat$datetime[c(cp_s2, cp_s3_1, cp_s3_2)], "%Y-%m-%d %H:%M:%S"),
  source = c("原始第2个候选变点", "原始第2个候选变点", "原始第3个候选变点")
)

write.csv(
  base_anchor_table,
  file.path(out_dir, "02_artificial_anchor_points.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# 人工均值变化施加方式（只影响情形2和情形3；情形1不调用该函数）：
#   本版按你的要求，将人工变化施加在变点左侧，即第二阶段局部检验窗口的左半段。
#   对候选变点 tau，第二阶段窗口为 [tau - L + 1, tau + L]；
#   人工变化只施加在左半窗 [tau - L + 1, tau]。
#   这样在局部窗口内仍表现为一个均值变点，并且不会影响该变点右侧后续序列。
injection_pre_hours <- L

apply_mean_shift_to_wide <- function(wide_dat, cp_index, shift_vec, scope, target_station, selected_stations) {
  out <- wide_dat

  if (!all(names(shift_vec) %in% c("PM10", "PM25"))) {
    stop("shift_vec 必须命名为 PM10 和 PM25。")
  }

  affected_stations <- if (scope == "all") {
    selected_stations
  } else if (scope == "target") {
    target_station
  } else {
    stop("scope 只能为 'all' 或 'target'。")
  }

  # 关键修改：只在局部检验窗口左半段施加人工变化。
  # 原始情形 shift_specs 为空，不会调用本函数，因此情形1完全不受影响。
  start_idx <- max(1L, cp_index - injection_pre_hours + 1L)
  end_idx <- cp_index
  rows <- if (start_idx <= end_idx) start_idx:end_idx else integer(0)

  if (length(rows) == 0) return(out)

  for (st in affected_stations) {
    out[[paste0("PM10_", st)]][rows] <- out[[paste0("PM10_", st)]][rows] + unname(shift_vec["PM10"])
    out[[paste0("PM25_", st)]][rows] <- out[[paste0("PM25_", st)]][rows] + unname(shift_vec["PM25"])
  }

  out
}

apply_scenario_shifts <- function(wide_dat, shift_specs) {
  out <- wide_dat

  if (length(shift_specs) == 0) {
    return(out)
  }

  # 按变点顺序施加。
  shift_specs <- shift_specs[order(vapply(shift_specs, function(x) x$cp_index, numeric(1)))]

  for (sp in shift_specs) {
    out <- apply_mean_shift_to_wide(
      wide_dat = out,
      cp_index = sp$cp_index,
      shift_vec = sp$shift_vec,
      scope = sp$scope,
      target_station = target_station,
      selected_stations = selected_stations
    )
  }

  out
}

scenario_list <- list(
  original = list(
    scenario_id = "S1_original",
    scenario_name = "情形1：原始时间段（不做人工处理）",
    shift_specs = list()
  ),

  one_artificial_cp = list(
    scenario_id = "S2_one_artificial_cp",
    scenario_name = "情形2：原始第2个候选变点左侧施加目标站点双变量异常",
    shift_specs = list(
      list(
        cp_index = cp_s2,
        cp_time = wide_dat$datetime[cp_s2],
        shift_vec = auto_shifts$both_vars,
        scope = "target",
        true_label = "异常",
        change_desc = "原始第2个候选变点左侧：目标站点 PM10 与 PM2.5 同时发生小幅均值变化"
      )
    )
  ),

  two_artificial_cp = list(
    scenario_id = "S3_two_artificial_cp",
    scenario_name = "情形3：原始第2和第3个候选变点左侧均施加目标站点异常",
    shift_specs = list(
      list(
        cp_index = cp_s3_1,
        cp_time = wide_dat$datetime[cp_s3_1],
        shift_vec = auto_shifts$both_vars,
        scope = "target",
        true_label = "异常",
        change_desc = "原始第2个候选变点左侧：目标站点 PM10 与 PM2.5 同时发生小幅均值变化"
      ),
      list(
        cp_index = cp_s3_2,
        cp_time = wide_dat$datetime[cp_s3_2],
        shift_vec = auto_shifts$pm10_only,
        scope = "target",
        true_label = "异常",
        change_desc = "原始第3个候选变点左侧：目标站点仅 PM10 发生小幅均值变化"
      )
    )
  )
)

cat("\n人工变化施加方式：局部检验窗口左半段；持续小时数 = ", injection_pre_hours, "\n", sep = "")
cat("\n人工变点设置：\n")
for (sc in scenario_list) {
  cat("\n", sc$scenario_name, "\n", sep = "")
  if (length(sc$shift_specs) == 0) {
    cat("  无人工处理。\n")
  } else {
    for (sp in sc$shift_specs) {
      inj_start <- max(1L, sp$cp_index - injection_pre_hours + 1L)
      inj_end <- sp$cp_index
      cat(
        "  位置 =", sp$cp_index,
        "；时间 =", format(sp$cp_time, "%Y-%m-%d %H:%M"),
        "；范围 =", sp$scope,
        "；真实类别 =", sp$true_label,
        "；变化 = (PM10:", sp$shift_vec["PM10"], ", PM2.5:", sp$shift_vec["PM25"], ")",
        "；左侧局部施加区间 =", format(wide_dat$datetime[inj_start], "%Y-%m-%d %H:%M"),
        "至", format(wide_dat$datetime[inj_end], "%Y-%m-%d %H:%M"), "\n"
      )
    }
  }
}


# ============================================================
# 7. 结果标注、绘图与保存
# ============================================================

annotate_candidate_truth <- function(candidate_df, shift_specs, tolerance = 6) {
  if (nrow(candidate_df) == 0) {
    return(candidate_df |>
             dplyr::mutate(
               nearest_true_cp = integer(0),
               nearest_true_time = as.POSIXct(character(0), tz = "Asia/Shanghai"),
               true_label = character(0),
               change_desc = character(0),
               abs_error = numeric(0)
             ))
  }

  if (length(shift_specs) == 0) {
    return(candidate_df |>
             dplyr::mutate(
               nearest_true_cp = NA_integer_,
               nearest_true_time = as.POSIXct(NA, tz = "Asia/Shanghai"),
               true_label = "未预设",
               change_desc = "原始数据中的候选变点",
               abs_error = NA_real_
             ))
  }

  true_cps <- vapply(shift_specs, function(x) x$cp_index, numeric(1))

  candidate_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      nearest_idx = which.min(abs(tau - true_cps)),
      nearest_true_cp = as.integer(true_cps[nearest_idx]),
      nearest_true_time = wide_dat$datetime[nearest_true_cp],
      abs_error = abs(tau - nearest_true_cp),
      true_label = ifelse(
        abs_error <= tolerance,
        shift_specs[[nearest_idx]]$true_label,
        "未匹配"
      ),
      change_desc = ifelse(
        abs_error <= tolerance,
        shift_specs[[nearest_idx]]$change_desc,
        "未匹配人工变点"
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-nearest_idx)
}

make_plot_data <- function(wide_dat_scenario, target_station, idw_weights) {
  target_mat <- get_station_matrix_from_wide(wide_dat_scenario, target_station)
  ref_mat <- calc_weighted_reference(
    wide_dat = wide_dat_scenario,
    target_station = target_station,
    idw_weights = idw_weights
  )

  df_target <- tibble::tibble(
    datetime = wide_dat_scenario$datetime,
    `PM10-目标站点` = target_mat[, "PM10"],
    `PM2.5-目标站点` = target_mat[, "PM2.5"]
  ) |>
    tidyr::pivot_longer(
      cols = -datetime,
      names_to = "序列",
      values_to = "取值"
    )

  df_ref <- tibble::tibble(
    datetime = wide_dat_scenario$datetime,
    `PM10-临近站点加权值` = ref_mat[, "PM10"],
    `PM2.5-临近站点加权值` = ref_mat[, "PM2.5"]
  ) |>
    tidyr::pivot_longer(
      cols = -datetime,
      names_to = "序列",
      values_to = "取值"
    )

  list(df_target = df_target, df_ref = df_ref)
}

plot_scenario_result <- function(
    wide_dat_scenario,
    scenario,
    candidate_df,
    idw_weights,
    save_path
) {
  plot_data <- make_plot_data(
    wide_dat_scenario = wide_dat_scenario,
    target_station = target_station,
    idw_weights = idw_weights
  )

  color_values <- c(
    "PM10-目标站点" = "#B2182B",
    "PM10-临近站点加权值" = "#EF8A62",
    "PM2.5-目标站点" = "#2166AC",
    "PM2.5-临近站点加权值" = "#67A9CF"
  )

  linetype_values <- c(
    "PM10-目标站点" = "solid",
    "PM10-临近站点加权值" = "longdash",
    "PM2.5-目标站点" = "solid",
    "PM2.5-临近站点加权值" = "longdash"
  )

  linewidth_values <- c(
    "PM10-目标站点" = 0.88,
    "PM10-临近站点加权值" = 0.70,
    "PM2.5-目标站点" = 0.88,
    "PM2.5-临近站点加权值" = 0.70
  )

  legend_order <- names(color_values)
  legend_position_use <- ifelse(show_series_legend, "bottom", "none")

  p <- ggplot() +
    geom_line(
      data = plot_data$df_ref,
      aes(x = datetime, y = 取值, color = 序列, linetype = 序列, linewidth = 序列),
      alpha = 0.78
    ) +
    geom_line(
      data = plot_data$df_target,
      aes(x = datetime, y = 取值, color = 序列, linetype = 序列, linewidth = 序列),
      alpha = 0.96
    ) +
    scale_color_manual(values = color_values, breaks = legend_order) +
    scale_linetype_manual(values = linetype_values, breaks = legend_order) +
    scale_linewidth_manual(values = linewidth_values, breaks = legend_order) +
    scale_x_datetime(
      date_labels = "%m-%d",
      date_breaks = x_date_breaks
    ) +
    labs(
      title = scenario$scenario_name,
      subtitle = paste0(
        "目标站点：", target_station_zh,
        "；绿色竖线表示正常，红色竖线表示异常"
      ),
      x = "时间",
      y = "浓度",
      color = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    theme_minimal(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = title_size),
      plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey35"),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      legend.position = legend_position_use,
      legend.box = "vertical",
      legend.text = element_text(size = legend_text_size),
      legend.key.width = grid::unit(2.4, "cm"),
      panel.border = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
      plot.margin = margin(10, 16, 10, 10)
    )

  # 画第二阶段判定结果。正常为绿色，异常为红色。
  if (!is.null(candidate_df) && nrow(candidate_df) > 0) {
    for (i in seq_len(nrow(candidate_df))) {
      cp_color <- ifelse(candidate_df$pred_label[i] == "异常", "#D62728", "#2CA02C")
      p <- p +
        geom_vline(
          xintercept = as.numeric(candidate_df$datetime[i]),
          color = cp_color,
          linetype = "solid",
          linewidth = 1.05,
          alpha = 0.95
        )
    }
  }

  ggsave(
    filename = save_path,
    plot = p,
    width = 12.5,
    height = 7.2,
    dpi = 300,
    bg = "white"
  )

  message("图已保存：", normalizePath(save_path, winslash = "/", mustWork = FALSE))
  invisible(p)
}


# ============================================================
# 8. 主流程：三种情形逐一检测、识别、绘图
# ============================================================

all_result_list <- list()
all_plot_list <- list()

for (sc in scenario_list) {
  cat("\n============================================================\n")
  cat("正在处理：", sc$scenario_name, "\n")
  cat("============================================================\n")

  wide_scenario <- apply_scenario_shifts(
    wide_dat = wide_dat,
    shift_specs = sc$shift_specs
  )

  stage1_res <- detect_candidate_cps(
    wide_dat = wide_scenario,
    target_station = target_station,
    confidence = sn_confidence,
    grid_size_scale = sn_grid_size_scale
  )

  # 保存该情形重新运行第一阶段后的候选变点。
  stage1_save <- tibble::tibble(
    scenario_id = sc$scenario_id,
    scenario_name = sc$scenario_name,
    tau = stage1_res$cp_index,
    datetime = format(stage1_res$cp_datetime, "%Y-%m-%d %H:%M:%S")
  )
  if (nrow(stage1_save) == 0) {
    stage1_save <- tibble::tibble(
      scenario_id = sc$scenario_id,
      scenario_name = sc$scenario_name,
      tau = NA_integer_,
      datetime = NA_character_
    )
  }
  write.csv(
    stage1_save,
    file.path(out_dir, paste0(sc$scenario_id, "_stage1_candidates.csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  candidate_df <- classify_candidates(
    wide_dat = wide_scenario,
    candidate_cps = stage1_res$cp_index,
    L = L,
    M = M,
    block_size = block_size,
    lambda_frac = lambda_frac,
    alpha = alpha,
    target_station = target_station,
    idw_weights = idw_weights
  )

  candidate_df <- annotate_candidate_truth(
    candidate_df = candidate_df,
    shift_specs = sc$shift_specs,
    tolerance = cp_match_tolerance
  )

  if (nrow(candidate_df) == 0) {
    cat("未检测到可进入第二阶段的候选变点。\n")
  } else {
    cat("候选变点与第二阶段分类结果：\n")
    print(candidate_df)
  }

  result_save <- candidate_df |>
    dplyr::mutate(
      scenario_id = sc$scenario_id,
      scenario_name = sc$scenario_name,
      datetime = format(datetime, "%Y-%m-%d %H:%M:%S"),
      nearest_true_time = ifelse(
        is.na(nearest_true_time),
        NA_character_,
        format(nearest_true_time, "%Y-%m-%d %H:%M:%S")
      )
    ) |>
    dplyr::select(
      scenario_id,
      scenario_name,
      tau,
      datetime,
      H_obs,
      H_boot_q95,
      p_value,
      pred_label,
      nearest_true_cp,
      nearest_true_time,
      abs_error,
      true_label,
      change_desc
    )

  if (nrow(result_save) == 0) {
    result_save <- tibble::tibble(
      scenario_id = sc$scenario_id,
      scenario_name = sc$scenario_name,
      tau = NA_integer_,
      datetime = NA_character_,
      H_obs = NA_real_,
      H_boot_q95 = NA_real_,
      p_value = NA_real_,
      pred_label = NA_character_,
      nearest_true_cp = NA_integer_,
      nearest_true_time = NA_character_,
      abs_error = NA_real_,
      true_label = NA_character_,
      change_desc = "未检测到候选变点"
    )
  }

  all_result_list[[sc$scenario_id]] <- result_save

  result_csv <- file.path(out_dir, paste0(sc$scenario_id, "_classification_result.csv"))
  write.csv(result_save, result_csv, row.names = FALSE, fileEncoding = "UTF-8")

  plot_path <- file.path(out_dir, paste0(sc$scenario_id, "_anomaly_identification.png"))
  p <- plot_scenario_result(
    wide_dat_scenario = wide_scenario,
    scenario = sc,
    candidate_df = candidate_df,
    idw_weights = idw_weights,
    save_path = plot_path
  )

  all_plot_list[[sc$scenario_id]] <- p
}

all_results <- dplyr::bind_rows(all_result_list)

write.csv(
  all_results,
  file.path(out_dir, "00_all_scenarios_classification_result.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n============================================================\n")
cat("全部完成。输出文件夹：\n")
cat(normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("主要输出：\n")
cat("  01_IDW_weights.csv：反距离加权权重\n")
cat("  02_artificial_anchor_points.csv：原始候选变点锚点设置\n")
cat("  00_all_scenarios_classification_result.csv：三种情形汇总结果\n")
cat("  S1/S2/S3_*_stage1_candidates.csv：各情形重新运行第一阶段后的候选变点\n")
cat("  S1/S2/S3_*_classification_result.csv：各情形分类结果\n")
cat("  S1/S2/S3_*_anomaly_identification.png：各情形异常识别图\n")
cat("============================================================\n")
