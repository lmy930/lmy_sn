# ============================================================
# 实际数据双变量异常识别：SNSeg_Multi + 反距离加权空间差值 + CBB Hotelling 检验
# 锚点式完整流程版
#
# 【修复记录】
#   v1：修复随机数流问题
#       对每个 tau 用 set.seed(global_seed + tau)，保证三种情形可比。
#   v2：修复异方差问题
#       Z(t) 改为对数差值：log(目标+1) - log(IDW参考+1)
#       从根本上消除浓度水平对差值方差的乘性影响。
#   v3：人工突变与对数差值尺度匹配
#       突变方式从绝对加减改为乘性比例，保证在对数尺度下
#       突变强度均匀，不随浓度水平变化。
#       绘图仍使用原始浓度序列，不受影响。
#   v4：人工突变改回加法（原始浓度单位）
#       模拟阶段只需在原始浓度上施加固定偏移即可验证算法，
#       物理含义清晰（传感器系统性偏低/偏高 δ μg/m³）。
#       算法内部仍使用对数差值 Z(t)，检验框架不变。
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

# ============================================================
# 0. 用户参数设置
# ============================================================

data_dir       <- "D:/03/2025/Beijing"
target_station <- "Dongsi"

start_time <- as.POSIXct("2024-12-07 00:00:00", tz = "Asia/Shanghai")
end_time   <- as.POSIXct("2024-12-16 23:00:00", tz = "Asia/Shanghai")

selected_stations <- c(
  "Dongsi", "Tiantan", "Guanyuan",
  "Wanshouxigong", "Nongzhanguan", "Aotizhongxin"
)

# 第一阶段 SNSeg 参数
sn_confidence      <- 0.95
sn_grid_size_scale <- 0.05

# 第二阶段参数
L           <- 10
M           <- 1000
lambda_frac <- 0.05
alpha       <- 0.05
global_seed <- 42

# 【v2】对数差值偏移常数，防止浓度为0时log报错
# 对高浓度数据（PM10/PM2.5通常>10）影响可忽略
log_offset <- 0

# 【v4】加法突变设置（原始浓度单位，μg/m³）
# δ < 0：目标站点系统性偏低；δ > 0：目标站点系统性偏高
# 建议设为该污染物正常均值的 20%–40%，太小不易检测，太大不像真实异常
# 两种人工场景的偏移量，主流程会用这里的值
shift_delta_both <- c(PM10 = -15, PM25 = -10)  # 情形2：PM10和PM2.5同时偏低
shift_delta_pm10 <- c(PM10 = -20, PM25 =   0)  # 情形3第二个变点：仅PM10偏低

choose_cbb_block_size <- function(L) {
  max(5L, min(as.integer(L), as.integer(round(sqrt(2 * L)))))
}
block_size <- choose_cbb_block_size(L)

cp_match_tolerance <- 6
show_series_legend <- TRUE

font_family      <- "SimHei"
base_size        <- 14
title_size       <- 18
axis_title_size  <- 14
axis_text_size   <- 12
legend_text_size <- 12
x_date_breaks    <- "1 day"

out_dir <- file.path(data_dir, "actual_stage2_anomaly_identification_v4")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# ============================================================
# 1. 站点经纬度与反距离权重
# ============================================================

station_coord <- tibble::tribble(
  ~station,         ~station_zh,  ~lon,     ~lat,
  "Dongsi",         "东四",       116.417,  39.929,
  "Tiantan",        "天坛",       116.407,  39.886,
  "Guanyuan",       "官园",       116.339,  39.929,
  "Wanshouxigong",  "万寿西宫",   116.352,  39.878,
  "Nongzhanguan",   "农展馆",     116.461,  39.937,
  "Aotizhongxin",   "奥体中心",   116.397,  39.982,
  "Changping",      "昌平",       116.231,  40.217,
  "Dingling",       "定陵",       116.220,  40.292,
  "Gucheng",        "古城",       116.184,  39.914,
  "Huairou",        "怀柔",       116.628,  40.328,
  "Shunyi",         "顺义",       116.655,  40.127,
  "Wanliu",         "万柳",       116.287,  39.987
)

target_station_zh <- station_coord$station_zh[
  match(target_station, station_coord$station)]
if (is.na(target_station_zh)) target_station_zh <- target_station

haversine_km <- function(lon1, lat1, lon2, lat2) {
  R <- 6371; to_rad <- pi / 180
  dlon <- (lon2 - lon1) * to_rad
  dlat <- (lat2 - lat1) * to_rad
  a <- sin(dlat/2)^2 + cos(lat1*to_rad) * cos(lat2*to_rad) * sin(dlon/2)^2
  2 * R * atan2(sqrt(a), sqrt(1 - a))
}

calc_idw_weights <- function(station_coord, target_station, selected_stations) {
  target_row <- station_coord |> dplyr::filter(station == target_station)
  if (nrow(target_row) != 1) stop("目标站点经纬度缺失或不唯一。")
  ref_df <- station_coord |>
    dplyr::filter(station %in% setdiff(selected_stations, target_station)) |>
    dplyr::mutate(
      distance_km  = haversine_km(target_row$lon, target_row$lat, lon, lat),
      inv_distance = 1 / pmax(distance_km, 1e-6),
      weight       = inv_distance / sum(inv_distance)
    )
  if (nrow(ref_df) == 0) stop("没有可用参考站点。")
  ref_df |>
    dplyr::select(station, station_zh, lon, lat, distance_km, weight) |>
    dplyr::arrange(desc(weight))
}

idw_weights <- calc_idw_weights(station_coord, target_station, selected_stations)
write.csv(idw_weights, file.path(out_dir, "01_IDW_weights.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("\n反距离加权权重：\n"); print(idw_weights)


# ============================================================
# 2. 读取与整理实际数据
# ============================================================

extract_station_name <- function(file) {
  raw_name <- basename(file) |>
    stringr::str_remove("\\.csv$") |>
    stringr::str_remove("^Beijing_")
  if (raw_name %in% station_coord$station) return(raw_name)
  hits <- station_coord$station[
    stringr::str_detect(raw_name,
                        paste0("(^|_)", station_coord$station, "($|_)"))]
  if (length(hits) >= 1) return(hits[1])
  raw_name
}

read_one_station <- function(file) {
  sn  <- extract_station_name(file)
  dat <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  need_vars <- c("year", "month", "day", "hour", "PM2.5", "PM10")
  miss_vars <- setdiff(need_vars, names(dat))
  if (length(miss_vars) > 0)
    stop(paste0("文件 ", basename(file), " 缺少变量：",
                paste(miss_vars, collapse = ", ")))
  dat |> dplyr::transmute(
    station  = sn,
    datetime = lubridate::make_datetime(year, month, day, hour,
                                        tz = "Asia/Shanghai"),
    PM25 = as.numeric(`PM2.5`),
    PM10 = as.numeric(PM10)
  )
}

files <- list.files(data_dir, pattern = "^Beijing_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("未找到 Beijing_*.csv 文件，请检查 data_dir。")

raw_long <- purrr::map_dfr(files, read_one_station) |>
  dplyr::filter(station %in% selected_stations)

miss_station <- setdiff(selected_stations, unique(raw_long$station))
if (length(miss_station) > 0)
  stop(paste0("以下站点没有读取到数据：", paste(miss_station, collapse = ", ")))

wide_dat_all <- raw_long |>
  dplyr::filter(datetime >= start_time, datetime <= end_time) |>
  tidyr::pivot_wider(
    id_cols     = datetime,
    names_from  = station,
    values_from = c(PM10, PM25),
    names_glue  = "{.value}_{station}"
  ) |>
  dplyr::arrange(datetime)

if (nrow(wide_dat_all) == 0)
  stop("当前时间段没有数据，请检查 start_time 和 end_time。")

full_time <- tibble::tibble(
  datetime = seq(min(wide_dat_all$datetime, na.rm = TRUE),
                 max(wide_dat_all$datetime, na.rm = TRUE),
                 by = "hour")
)

wide_dat <- full_time |>
  dplyr::left_join(wide_dat_all, by = "datetime") |>
  dplyr::arrange(datetime)

required_cols <- as.vector(outer(c("PM10","PM25"), selected_stations,
                                  paste, sep = "_"))
miss_cols <- setdiff(required_cols, names(wide_dat))
if (length(miss_cols) > 0)
  stop(paste0("宽格式数据缺少列：", paste(miss_cols, collapse = ", ")))

if (anyNA(wide_dat[, required_cols, drop = FALSE])) {
  miss_rows <- wide_dat |>
    dplyr::filter(if_any(dplyr::all_of(required_cols), is.na)) |>
    dplyr::select(datetime, dplyr::all_of(required_cols))
  print(head(miss_rows, 20))
  stop("当前时间段存在缺失小时或污染物缺失值，请更换时间段或先处理缺失。")
}

n_time <- nrow(wide_dat)
if (n_time < 2 * L + 20)
  stop("当前时间段过短，无法进行局部窗口异常识别。")

cat("\n数据读取完成：", n_time, "个小时。\n")
cat("时间范围：", format(min(wide_dat$datetime)), "至",
    format(max(wide_dat$datetime)), "\n")


# ============================================================
# 3. Rcpp：CBB + 正则化 Hotelling 型 H 统计量
# ============================================================

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; using namespace arma;

static inline mat estimate_reg_cov(const mat& Zc, double lf) {
  mat S = cov(Zc); int d = S.n_cols;
  double trS = trace(S);
  double lam = lf * trS / std::max(1, d);
  if (!std::isfinite(lam) || lam <= 1e-10) lam = 1e-6;
  S += lam * eye<mat>(d, d);
  return S;
}

static inline double hotelling(const mat& Z, int L, const mat& Si) {
  rowvec mb = mean(Z.rows(0, L-1), 0);
  rowvec ma = mean(Z.rows(L, 2*L-1), 0);
  rowvec dv = ma - mb;
  mat v = dv * Si * dv.t();
  return (static_cast<double>(L) / 2.0) * v(0, 0);
}

// [[Rcpp::export]]
List block_boot_cbb_hotelling_cpp(
    const arma::mat& Z, int L, int M = 500,
    int block_size = 10, double lambda_frac = 0.05) {
  int n = Z.n_rows, d = Z.n_cols;
  if (n != 2 * L) stop("Z.n_rows must be exactly 2 * L.");
  if (L < 5 || d < 1 || M < 10) stop("Invalid inputs.");

  rowvec mb = mean(Z.rows(0, L-1), 0);
  rowvec ma = mean(Z.rows(L, 2*L-1), 0);
  mat Zt = Z;
  for (int t = 0; t < L; ++t) Zt.row(t) -= mb;
  for (int t = L; t < n; ++t) Zt.row(t) -= ma;

  mat Sr = estimate_reg_cov(Zt, lambda_frac);
  mat Si; bool ok = inv_sympd(Si, Sr);
  if (!ok) Si = pinv(Sr);

  double Ho = hotelling(Z, L, Si);
  int b = std::max(1, std::min(block_size, n));
  int ex = 0;
  vec Hb(M, fill::zeros);

  for (int m = 0; m < M; ++m) {
    mat Zs(n, d, fill::zeros); int t = 0;
    while (t < n) {
      int si = static_cast<int>(std::floor(R::runif(0.0, (double)n)));
      for (int i = 0; i < b && t < n; ++i) {
        Zs.row(t) = Zt.row((si + i) % n); t++;
      }
    }
    double Hs = hotelling(Zs, L, Si);
    Hb(m) = Hs;
    if (Hs >= Ho) ex++;
  }

  double pv = static_cast<double>(ex + 1) / static_cast<double>(M + 1);
  vec Hs = sort(Hb);
  int i95 = std::min(M-1, std::max(0, (int)std::ceil(0.95*M) - 1));
  return List::create(
    _["p_value"]    = pv,
    _["H_obs"]      = Ho,
    _["H_boot_q95"] = Hs(i95)
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
    stop(paste0("以下变量标准差为0或不可用：", paste(bad, collapse = ", ")))
  }
  scale(x)
}

extract_cp_from_snseg <- function(res, n) {
  if (is.null(res)) return(integer(0))
  for (nm in c("est_cp","cp","cpts","changepoints","change_points","cp_index")) {
    if (nm %in% names(res)) {
      cp <- as.integer(round(as.numeric(unlist(res[[nm]]))))
      return(sort(unique(cp[!is.na(cp) & cp > 0 & cp < n])))
    }
  }
  message("未从SNSeg返回对象中识别到变点字段，字段名：")
  print(names(res))
  integer(0)
}

get_target_matrix <- function(wide_dat, target_station) {
  x <- as.matrix(wide_dat[, c(paste0("PM10_", target_station),
                               paste0("PM25_", target_station))])
  colnames(x) <- c("PM10", "PM2.5")
  x
}

detect_candidate_cps <- function(wide_dat, target_station,
                                 confidence = 0.95,
                                 grid_size_scale = 0.05) {
  x_std <- safe_scale_matrix(get_target_matrix(wide_dat, target_station))
  res <- tryCatch(
    SNSeg_Multi(
      ts              = x_std,
      paras_to_test   = "mean",
      confidence      = confidence,
      grid_size_scale = grid_size_scale,
      grid_size       = NULL,
      plot_SN         = FALSE,
      est_cp_loc      = TRUE
    ),
    error = function(e) { message("SNSeg_Multi出错：", e$message); NULL }
  )
  cp_index <- extract_cp_from_snseg(res, nrow(wide_dat))
  list(result      = res,
       cp_index    = cp_index,
       cp_datetime = wide_dat$datetime[cp_index])
}


# ============================================================
# 5. 第二阶段辅助函数
# ============================================================

get_station_matrix_from_wide <- function(wide_dat, station) {
  x <- as.matrix(wide_dat[, c(paste0("PM10_", station),
                               paste0("PM25_", station))])
  colnames(x) <- c("PM10", "PM2.5")
  x
}

calc_weighted_reference <- function(wide_dat, target_station, idw_weights) {
  ref_mat <- matrix(0, nrow = nrow(wide_dat), ncol = 2)
  colnames(ref_mat) <- c("PM10", "PM2.5")
  for (j in seq_along(idw_weights$station))
    ref_mat <- ref_mat +
      idw_weights$weight[j] *
      get_station_matrix_from_wide(wide_dat, idw_weights$station[j])
  ref_mat
}

# ============================================================
# 【v2】build_residual_window_raw：Z改为对数差值
#
# Z(t) = log(目标(t) + log_offset) - log(IDW参考(t) + log_offset)
#
# 物理含义：目标站点相对区域背景的对数比偏离。
# 高浓度时期各站点比值若稳定，对数差值就稳定，
# 从根本上消除浓度水平对差值方差的乘性影响。
# ============================================================

build_residual_window_raw <- function(wide_dat, tau, L,
                                      target_station, idw_weights) {
  s_idx <- tau - L + 1
  e_idx <- tau + L
  if (s_idx < 1 || e_idx > nrow(wide_dat)) return(NULL)

  win_dat    <- wide_dat[s_idx:e_idx, , drop = FALSE]
  target_mat <- get_station_matrix_from_wide(win_dat, target_station)
  ref_mat    <- calc_weighted_reference(win_dat, target_station, idw_weights)

  # 对数差值
  Z_log <- log(target_mat + log_offset) - log(ref_mat + log_offset)
  colnames(Z_log) <- c("log比值_PM10", "log比值_PM2.5")

  list(Z_raw = Z_log, s_idx = s_idx, e_idx = e_idx)
}

# ============================================================
# 【v1】test_one_candidate：每个tau用固定局部种子
# ============================================================

test_one_candidate <- function(wide_dat, tau, L, M, block_size, lambda_frac,
                               target_station, idw_weights,
                               global_seed = 42) {
  feat <- build_residual_window_raw(wide_dat, tau, L,
                                    target_station, idw_weights)
  if (is.null(feat)) return(NULL)
  set.seed(global_seed + tau)
  out <- block_boot_cbb_hotelling_cpp(
    Z = feat$Z_raw, L = L, M = M,
    block_size = block_size, lambda_frac = lambda_frac
  )
  data.frame(
    tau        = tau,
    datetime   = wide_dat$datetime[tau],
    H_obs      = as.numeric(out$H_obs),
    H_boot_q95 = as.numeric(out$H_boot_q95),
    p_value    = as.numeric(out$p_value),
    stringsAsFactors = FALSE
  )
}

classify_candidates <- function(wide_dat, candidate_cps, L, M, block_size,
                                lambda_frac, alpha, target_station,
                                idw_weights, global_seed = 42) {
  empty <- tibble::tibble(
    tau        = integer(0),
    datetime   = as.POSIXct(character(0), tz = "Asia/Shanghai"),
    H_obs      = numeric(0),
    H_boot_q95 = numeric(0),
    p_value    = numeric(0),
    pred_label = character(0)
  )
  if (length(candidate_cps) == 0) return(empty)

  res_df <- dplyr::bind_rows(lapply(candidate_cps, function(tau)
    test_one_candidate(wide_dat, tau, L, M, block_size, lambda_frac,
                       target_station, idw_weights, global_seed)))

  if (nrow(res_df) == 0) return(empty)

  res_df |>
    dplyr::mutate(pred_label = ifelse(p_value < alpha, "异常", "正常")) |>
    dplyr::arrange(tau)
}


# ============================================================
# 6. 人工均值变点设置与施加
# ============================================================

clamp_cp_index <- function(idx, n, L) {
  as.integer(max(L + 2L, min(n - L - 2L, round(idx))))
}

nearest_time_index <- function(datetime_vec, target_time) {
  which.min(abs(as.numeric(datetime_vec - target_time)))
}

# 用原始数据确定锚点
base_stage1        <- detect_candidate_cps(wide_dat, target_station,
                                           sn_confidence, sn_grid_size_scale)
base_candidate_cps <- base_stage1$cp_index
base_candidate_times <- base_stage1$cp_datetime

cat("\n原始数据第一阶段候选均值变点：\n")
if (length(base_candidate_cps) == 0) {
  cat("  未检测到候选变点。\n")
} else {
  print(data.frame(
    序号 = seq_along(base_candidate_cps),
    位置 = base_candidate_cps,
    时间 = format(base_candidate_times, "%Y-%m-%d %H:%M")
  ))
}

manual_cp_order_scenario2 <- 3L
manual_cp_order_scenario3 <- c(3L, 2L)

choose_anchor_cp <- function(cps, k, fb, n, L) {
  if (length(cps) >= k && !is.na(cps[k]))
    return(clamp_cp_index(cps[k], n, L))
  warning(paste0("候选变点数不足，用比例位置 ", fb, " 兜底"))
  clamp_cp_index(round(fb * n), n, L)
}

cp_s2   <- choose_anchor_cp(base_candidate_cps,
                             manual_cp_order_scenario2, 0.50, n_time, L)
cp_s3_1 <- choose_anchor_cp(base_candidate_cps,
                             manual_cp_order_scenario3[1], 0.35, n_time, L)
cp_s3_2 <- choose_anchor_cp(base_candidate_cps,
                             manual_cp_order_scenario3[2], 0.65, n_time, L)

if (cp_s3_2 <= cp_s3_1 + 2 * L)
  warning("情形3两个人工变点距离较近，可能影响局部窗口检验。")

base_anchor_table <- tibble::tibble(
  anchor_role = c("情形2锚点", "情形3第一个锚点", "情形3第二个锚点"),
  cp_index    = c(cp_s2, cp_s3_1, cp_s3_2),
  cp_time     = format(wide_dat$datetime[c(cp_s2, cp_s3_1, cp_s3_2)],
                       "%Y-%m-%d %H:%M:%S"),
  source      = c("原始第2个候选变点", "原始第2个候选变点", "原始第3个候选变点")
)
write.csv(base_anchor_table,
          file.path(out_dir, "02_artificial_anchor_points.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

injection_pre_hours <- L

# ============================================================
# 【v4】apply_mean_shift_to_wide：改回加法突变（原始浓度单位）
#
# 直接在原始浓度上加减固定值 delta（μg/m³），物理含义清晰：
# 模拟传感器系统性偏低/偏高固定浓度。
# 算法内部的对数差值 Z(t) 和 Hotelling 检验框架不受影响。
# ============================================================

apply_mean_shift_to_wide <- function(wide_dat, cp_index, shift_vec, scope,
                                     target_station, selected_stations) {
  out <- wide_dat
  if (!all(names(shift_vec) %in% c("PM10", "PM25")))
    stop("shift_vec 必须命名为 PM10 和 PM25。")

  stations <- if (scope == "all") selected_stations else
    if (scope == "target") target_station else
      stop("scope 只能为 'all' 或 'target'。")

  rows <- seq(max(1L, cp_index - injection_pre_hours + 1L), cp_index)
  if (length(rows) == 0) return(out)

  for (st in stations) {
    # 【v4关键修改】加法突变，在原始浓度上直接加减固定偏移量（单位：μg/m³）
    out[[paste0("PM10_", st)]][rows] <-
      out[[paste0("PM10_", st)]][rows] + unname(shift_vec["PM10"])
    out[[paste0("PM25_", st)]][rows] <-
      out[[paste0("PM25_", st)]][rows] + unname(shift_vec["PM25"])
  }
  # 防止出现负浓度值
  for (st in stations) {
    out[[paste0("PM10_", st)]][rows] <-
      pmax(0, out[[paste0("PM10_", st)]][rows])
    out[[paste0("PM25_", st)]][rows] <-
      pmax(0, out[[paste0("PM25_", st)]][rows])
  }
  out
}

apply_scenario_shifts <- function(wide_dat, shift_specs) {
  if (length(shift_specs) == 0) return(wide_dat)
  shift_specs <- shift_specs[order(
    vapply(shift_specs, function(x) x$cp_index, numeric(1)))]
  for (sp in shift_specs)
    wide_dat <- apply_mean_shift_to_wide(
      wide_dat, sp$cp_index, sp$shift_vec,
      sp$scope, target_station, selected_stations)
  wide_dat
}

scenario_list <- list(
  original = list(
    scenario_id   = "S1_original",
    scenario_name = "东四站多均值正常变点",
    shift_specs   = list()
  ),
  one_artificial_cp = list(
    scenario_id   = "S2_one_artificial_cp",
    scenario_name = "东四站单均值异常变点",
    shift_specs   = list(
      list(
        cp_index    = cp_s2,
        cp_time     = wide_dat$datetime[cp_s2],
        shift_vec   = shift_delta_both,   # 加法偏移（μg/m³）
        scope       = "target",
        true_label  = "异常",
        change_desc = paste0(
          "原始第2个候选变点左侧：目标站点PM10+", shift_delta_both["PM10"],
          "，PM2.5+", shift_delta_both["PM25"])
      )
    )
  ),
  two_artificial_cp = list(
    scenario_id   = "S3_two_artificial_cp",
    scenario_name = "东四站多均值异常变点",
    shift_specs   = list(
      list(
        cp_index    = cp_s3_1,
        cp_time     = wide_dat$datetime[cp_s3_1],
        shift_vec   = shift_delta_both,   # 加法偏移（μg/m³）
        scope       = "target",
        true_label  = "异常",
        change_desc = paste0(
          "原始第2个候选变点左侧：目标站点PM10+", shift_delta_both["PM10"],
          "，PM2.5+", shift_delta_both["PM25"])
      ),
      list(
        cp_index    = cp_s3_2,
        cp_time     = wide_dat$datetime[cp_s3_2],
        shift_vec   = shift_delta_pm10,   # 加法偏移（μg/m³）
        scope       = "target",
        true_label  = "异常",
        change_desc = paste0(
          "原始第3个候选变点左侧：目标站点仅PM10+", shift_delta_pm10["PM10"])
      )
    )
  )
)

cat("\n人工变化：局部检验窗口左半段；持续小时数 =", injection_pre_hours, "\n")
cat("\n人工变点设置：\n")
for (sc in scenario_list) {
  cat("\n", sc$scenario_name, "\n", sep = "")
  if (length(sc$shift_specs) == 0) { cat("  无人工处理。\n"); next }
  for (sp in sc$shift_specs) {
    s <- max(1L, sp$cp_index - injection_pre_hours + 1L)
    cat("  位置 =", sp$cp_index,
        "；时间 =", format(sp$cp_time, "%Y-%m-%d %H:%M"),
        "；范围 =", sp$scope,
        "；变化 = (PM10+", sp$shift_vec["PM10"],
        ", PM2.5+", sp$shift_vec["PM25"], ")",
        "；区间 =", format(wide_dat$datetime[s], "%Y-%m-%d %H:%M"),
        "至", format(wide_dat$datetime[sp$cp_index], "%Y-%m-%d %H:%M"), "\n")
  }
}


# ============================================================
# 7. 结果标注与绘图
# ============================================================

annotate_candidate_truth <- function(candidate_df, shift_specs,
                                     tolerance = 6) {
  if (nrow(candidate_df) == 0)
    return(candidate_df |> dplyr::mutate(
      nearest_true_cp   = integer(0),
      nearest_true_time = as.POSIXct(character(0), tz = "Asia/Shanghai"),
      true_label        = character(0),
      change_desc       = character(0),
      abs_error         = numeric(0)
    ))

  if (length(shift_specs) == 0)
    return(candidate_df |> dplyr::mutate(
      nearest_true_cp   = NA_integer_,
      nearest_true_time = as.POSIXct(NA, tz = "Asia/Shanghai"),
      true_label        = "未预设",
      change_desc       = "原始数据中的候选变点",
      abs_error         = NA_real_
    ))

  true_cps <- vapply(shift_specs, function(x) x$cp_index, numeric(1))

  candidate_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      nearest_idx       = which.min(abs(tau - true_cps)),
      nearest_true_cp   = as.integer(true_cps[nearest_idx]),
      nearest_true_time = wide_dat$datetime[nearest_true_cp],
      abs_error         = abs(tau - nearest_true_cp),
      true_label  = ifelse(abs_error <= tolerance,
                           shift_specs[[nearest_idx]]$true_label, "未匹配"),
      change_desc = ifelse(abs_error <= tolerance,
                           shift_specs[[nearest_idx]]$change_desc, "未匹配人工变点")
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-nearest_idx)
}

# 绘图使用原始浓度（不受对数变换影响）
make_plot_data <- function(wide_dat_scenario, target_station, idw_weights) {
  target_mat <- get_station_matrix_from_wide(wide_dat_scenario, target_station)
  ref_mat    <- calc_weighted_reference(wide_dat_scenario,
                                        target_station, idw_weights)
  df_target <- tibble::tibble(
    datetime          = wide_dat_scenario$datetime,
    `PM10-目标站点`   = target_mat[, "PM10"],
    `PM2.5-目标站点`  = target_mat[, "PM2.5"]
  ) |> tidyr::pivot_longer(-datetime, names_to = "序列", values_to = "取值")

  df_ref <- tibble::tibble(
    datetime               = wide_dat_scenario$datetime,
    `PM10-临近站点加权值`  = ref_mat[, "PM10"],
    `PM2.5-临近站点加权值` = ref_mat[, "PM2.5"]
  ) |> tidyr::pivot_longer(-datetime, names_to = "序列", values_to = "取值")

  list(df_target = df_target, df_ref = df_ref)
}

plot_scenario_result <- function(wide_dat_scenario, scenario, candidate_df,
                                 idw_weights, save_path) {
  pd <- make_plot_data(wide_dat_scenario, target_station, idw_weights)

  color_values <- c(
    "PM10-目标站点"        = "#B2182B",
    "PM10-临近站点加权值"  = "#EF8A62",
    "PM2.5-目标站点"       = "#2166AC",
    "PM2.5-临近站点加权值" = "#67A9CF"
  )
  linetype_values <- c(
    "PM10-目标站点"        = "solid",
    "PM10-临近站点加权值"  = "longdash",
    "PM2.5-目标站点"       = "solid",
    "PM2.5-临近站点加权值" = "longdash"
  )
  linewidth_values <- c(
    "PM10-目标站点"        = 0.88,
    "PM10-临近站点加权值"  = 0.70,
    "PM2.5-目标站点"       = 0.88,
    "PM2.5-临近站点加权值" = 0.70
  )
  lo <- names(color_values)

  p <- ggplot() +
    geom_line(data  = pd$df_ref,
              aes(datetime, 取值, color = 序列,
                  linetype = 序列, linewidth = 序列),
              alpha = 0.78) +
    geom_line(data  = pd$df_target,
              aes(datetime, 取值, color = 序列,
                  linetype = 序列, linewidth = 序列),
              alpha = 0.96) +
    scale_color_manual(values = color_values,      breaks = lo) +
    scale_linetype_manual(values = linetype_values, breaks = lo) +
    scale_linewidth_manual(values = linewidth_values, breaks = lo) +
    scale_x_datetime(date_labels = "%m-%d", date_breaks = x_date_breaks) +
    labs(
      title    = scenario$scenario_name,
    #   subtitle = paste0("目标站点：", target_station_zh,
    #                     "；绿色竖线表示正常，红色竖线表示异常"),
    x = "时间", y = "浓度",
      color = NULL, linetype = NULL, linewidth = NULL
    ) +
    theme_minimal(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
  
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      # axis.text.x = element_text(size = 13),
      # axis.text.y = element_text(size = 13),
  
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 14),
      legend.key.width = grid::unit(2.6, "cm"),
  
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray88"),
      
    #   plot.title       = element_text(face = "bold", hjust = 0.5,
    #                                   size = title_size),
    #   plot.subtitle    = element_text(hjust = 0.5, size = 11,
    #                                   color = "grey35"),
    #   axis.title       = element_text(size = axis_title_size),
    #   axis.text        = element_text(size = axis_text_size),
    #   legend.position  = ifelse(show_series_legend, "bottom", "none"),
    #   legend.box       = "vertical",
    #   legend.text      = element_text(size = legend_text_size),
    #   legend.key.width = grid::unit(2.4, "cm"),
      panel.border     = element_blank(),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
    #   panel.grid.minor = element_blank(),
    #   panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
      plot.margin      = margin(10, 16, 10, 10)
    )

  if (!is.null(candidate_df) && nrow(candidate_df) > 0)
    for (i in seq_len(nrow(candidate_df)))
      p <- p + geom_vline(
        xintercept = as.numeric(candidate_df$datetime[i]),
        color      = ifelse(candidate_df$pred_label[i] == "异常",
                            "#D62728", "#2CA02C"),
        linetype   = "solid", linewidth = 1.05, alpha = 0.95)

  ggsave(save_path, p, width = 12.5, height = 7.2, dpi = 300, bg = "white")
  message("图已保存：",
          normalizePath(save_path, winslash = "/", mustWork = FALSE))
  invisible(p)
}


# ============================================================
# 8. 主流程：三种情形逐一检测、识别、绘图
# ============================================================

all_result_list <- list()
all_plot_list   <- list()

for (sc in scenario_list) {
  cat("\n============================================================\n")
  cat("正在处理：", sc$scenario_name, "\n")
  cat("============================================================\n")

  wide_scenario <- apply_scenario_shifts(wide_dat, sc$shift_specs)

  stage1_res <- detect_candidate_cps(wide_scenario, target_station,
                                     sn_confidence, sn_grid_size_scale)

  stage1_save <- if (length(stage1_res$cp_index) == 0)
    tibble::tibble(scenario_id   = sc$scenario_id,
                   scenario_name = sc$scenario_name,
                   tau           = NA_integer_,
                   datetime      = NA_character_)
  else
    tibble::tibble(scenario_id   = sc$scenario_id,
                   scenario_name = sc$scenario_name,
                   tau           = stage1_res$cp_index,
                   datetime      = format(stage1_res$cp_datetime,
                                          "%Y-%m-%d %H:%M:%S"))

  write.csv(stage1_save,
            file.path(out_dir,
                      paste0(sc$scenario_id, "_stage1_candidates.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")

  candidate_df <- classify_candidates(
    wide_dat       = wide_scenario,
    candidate_cps  = stage1_res$cp_index,
    L              = L, M = M,
    block_size     = block_size,
    lambda_frac    = lambda_frac,
    alpha          = alpha,
    target_station = target_station,
    idw_weights    = idw_weights,
    global_seed    = global_seed
  )

  candidate_df <- annotate_candidate_truth(candidate_df, sc$shift_specs,
                                           cp_match_tolerance)

  if (nrow(candidate_df) == 0) {
    cat("未检测到可进入第二阶段的候选变点。\n")
  } else {
    cat("候选变点与第二阶段分类结果：\n")
    print(candidate_df)
  }

  result_save <- if (nrow(candidate_df) == 0)
    tibble::tibble(
      scenario_id       = sc$scenario_id,
      scenario_name     = sc$scenario_name,
      tau               = NA_integer_,
      datetime          = NA_character_,
      H_obs             = NA_real_,
      H_boot_q95        = NA_real_,
      p_value           = NA_real_,
      pred_label        = NA_character_,
      nearest_true_cp   = NA_integer_,
      nearest_true_time = NA_character_,
      abs_error         = NA_real_,
      true_label        = NA_character_,
      change_desc       = "未检测到候选变点"
    )
  else
    candidate_df |>
      dplyr::mutate(
        scenario_id       = sc$scenario_id,
        scenario_name     = sc$scenario_name,
        datetime          = format(datetime, "%Y-%m-%d %H:%M:%S"),
        nearest_true_time = ifelse(
          is.na(nearest_true_time), NA_character_,
          format(nearest_true_time, "%Y-%m-%d %H:%M:%S"))
      ) |>
      dplyr::select(scenario_id, scenario_name, tau, datetime,
                    H_obs, H_boot_q95, p_value, pred_label,
                    nearest_true_cp, nearest_true_time,
                    abs_error, true_label, change_desc)

  all_result_list[[sc$scenario_id]] <- result_save

  write.csv(result_save,
            file.path(out_dir,
                      paste0(sc$scenario_id, "_classification_result.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")

  p <- plot_scenario_result(
    wide_dat_scenario = wide_scenario,
    scenario          = sc,
    candidate_df      = candidate_df,
    idw_weights       = idw_weights,
    save_path         = file.path(
      out_dir, paste0(sc$scenario_id, "_anomaly_identification.png"))
  )
  all_plot_list[[sc$scenario_id]] <- p
}

write.csv(
  dplyr::bind_rows(all_result_list),
  file.path(out_dir, "00_all_scenarios_classification_result.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

cat("\n============================================================\n")
cat("全部完成。输出文件夹：\n")
cat(normalizePath(out_dir, winslash = "/", mustWork = FALSE), "\n")
cat("主要输出：\n")
cat("  01_IDW_weights.csv：反距离加权权重\n")
cat("  02_artificial_anchor_points.csv：原始候选变点锚点设置\n")
cat("  00_all_scenarios_classification_result.csv：三种情形汇总结果\n")
cat("  S1/S2/S3_*_stage1_candidates.csv：各情形第一阶段候选变点\n")
cat("  S1/S2/S3_*_classification_result.csv：各情形分类结果\n")
cat("  S1/S2/S3_*_anomaly_identification.png：各情形异常识别图\n")
cat("============================================================\n")