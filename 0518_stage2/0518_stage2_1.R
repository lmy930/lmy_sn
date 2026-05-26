# ============================================================
# 双变量单变点异常识别模拟：CBB + 正则化 Hotelling 空间差值检验 + 并行版
# 仅保留情形1和情形2
#
# 修改要点：
#   1. 第一阶段：仅使用目标站点自身 Z-score 标准化序列进行双变量均值变点检测；
#   2. 第二阶段：使用原始序列构造目标站点与邻近站点的空间差值序列；
#   3. 第二阶段：在候选变点局部窗口内基于空间差值序列构造正则化 Hotelling 型 H 统计量；
#   4. Bootstrap：使用 CBB 循环块自助法，默认块长 b = round(sqrt(2L))；
#   5. 每个候选变点单独检验，不合并候选点，不做 p 值校正；
#   6. 第一阶段指标：变点检测率、H 距离、ARI；
#   7. 第二阶段指标：分类准确率；
#   8. 并行运行 Monte Carlo 模拟；
#   9. 仅设置一个全局随机数种子：3047；
#  10. 保留可视化功能：在 Monte Carlo 计算过程中同步更新最佳图样本，不额外重复模拟。
# ============================================================

rm(list = ls())
options(scipen = 999, digits = 3)
suppressPackageStartupMessages({
  library(SNSeg)
  library(Rcpp)
  library(RcppArmadillo)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(parallel)
})

set.seed(3047)

get_default_cores <- function(max_cores = 8) {
  cores <- parallel::detectCores(logical = TRUE)
  if (is.na(cores) || !is.finite(cores) || cores < 2) {
    return(1)
  }
  max(1, min(max_cores, cores - 1))
}

choose_cbb_block_size <- function(L) {
  # 局部窗口长度为 2L。默认块长取 round(sqrt(2L))。
  # 当 L = 40 时，窗口长度为 80，默认块长为 9。
  max(5L, min(as.integer(L), as.integer(round(sqrt(2 * L)))))
}

# ============================================================
# 0. 参数设置
# ============================================================

N_sim <- 500
n <- 500
cp_true <- 250
loc_tol <- 20      # 仅用于最佳效果图筛选：第一阶段是否定位到真实变点附近

L <- 40
M <- 500
block_size <- choose_cbb_block_size(L)
lambda_frac <- 0.05
alpha <- 0.05

confidence <- 0.95
grid_size_scale <- 0.05

target_station <- "东四"
output_dir <- "single_cp_anomaly_Hotelling_CBB_case12_inline_best_figures"

# 并行参数：Windows / macOS / Linux 均可用 PSOCK cluster
use_parallel <- TRUE
n_cores <- get_default_cores()

# ----------------------------
# VAR(1) 模型参数：保持原样
# ----------------------------

A <- matrix(
  c(
    0.30, 0.10,
    0.12, 0.25
  ),
  nrow = 2,
  byrow = TRUE
)

sd_eps <- c(5.0, 5.0)
rho_eps <- 0.65

Sigma <- matrix(
  c(
    sd_eps[1]^2, rho_eps * prod(sd_eps),
    rho_eps * prod(sd_eps), sd_eps[2]^2
  ),
  nrow = 2,
  byrow = TRUE
)

mu0 <- c(25, 10)

delta_case1 <- c(0, -7.5)
delta_case2 <- c(2.5, -2.0)

# ============================================================
# 1. Rcpp：CBB 循环块自助法 + 正则化 Hotelling 型 H 统计量
# ============================================================

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

static inline mat estimate_reg_cov(const mat& Z_centered, double lambda_frac) {
  mat S = cov(Z_centered);
  int d = S.n_cols;

  if (d < 1) stop("Invalid dimension.");

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

  // 构造原假设下的局部均值对齐序列，用于协方差估计和 CBB 重抽样。
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
# 2. 基础工具函数
# ============================================================

rmvnorm2 <- function(Sigma) {
  as.numeric(t(chol(Sigma)) %*% rnorm(2))
}

safe_div <- function(a, b) {
  ifelse(b == 0, NA_real_, a / b)
}

mean_if_available <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

standardize_matrix <- function(x) {
  center <- colMeans(x)
  scale_val <- apply(x, 2, sd)
  scale_val[!is.finite(scale_val) | scale_val < 1e-8] <- 1

  out <- sweep(sweep(x, 2, center, "-"), 2, scale_val, "/")
  colnames(out) <- colnames(x)
  out
}

# ============================================================
# 3. 第一阶段评价函数：检测率、ARI、H 距离
# ============================================================

get_segment_labels <- function(n, cps) {
  labels <- rep(1L, n)

  cps <- sort(unique(as.integer(cps)))
  cps <- cps[cps > 0 & cps < n]

  if (length(cps) == 0) return(labels)

  for (i in seq_along(cps)) {
    labels[(cps[i] + 1):n] <- i + 1L
  }

  labels
}

ari_score <- function(label_true, label_est) {
  tab <- table(label_true, label_est)
  comb2 <- function(x) x * (x - 1) / 2

  n_total <- sum(tab)
  total <- comb2(n_total)

  if (total == 0) return(NA_real_)

  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))

  expected <- sum_i * sum_j / total
  max_index <- 0.5 * (sum_i + sum_j)

  if (max_index == expected) return(1)

  (sum_ij - expected) / (max_index - expected)
}

hausdorff_distance <- function(detected_cps, true_cps) {
  detected_cps <- sort(unique(as.numeric(detected_cps)))
  true_cps <- sort(unique(as.numeric(true_cps)))

  if (length(detected_cps) == 0 || length(true_cps) == 0) {
    return(NA_real_)
  }

  dist_mat <- abs(outer(detected_cps, true_cps, "-"))

  max(
    max(apply(dist_mat, 1, min)),
    max(apply(dist_mat, 2, min))
  )
}

calc_stage1_ari_hausdorff <- function(n, true_cp, detected_cps) {
  detected_cps <- detected_cps[detected_cps > 0 & detected_cps < n]

  if (length(detected_cps) == 0) {
    return(list(ARI = NA_real_, H = NA_real_))
  }

  true_labels <- get_segment_labels(n, true_cp)
  est_labels <- get_segment_labels(n, detected_cps)

  list(
    ARI = ari_score(true_labels, est_labels),
    H = hausdorff_distance(detected_cps, true_cp)
  )
}

# ============================================================
# 4. 站点信息与基础数据生成
# ============================================================

make_station_info <- function() {
  data.frame(
    station = c("东四", "天坛", "官园", "万寿西宫", "奥体中心", "农展馆"),

    mean_var1 = c(
      mu0[1],
      mu0[1] - 1.0,
      mu0[1] + 1.5,
      mu0[1] - 2.0,
      mu0[1] + 2.0,
      mu0[1] - 0.5
    ),

    mean_var2 = c(
      mu0[2],
      mu0[2] - 0.8,
      mu0[2] + 1.0,
      mu0[2] - 1.2,
      mu0[2] + 1.4,
      mu0[2] - 0.4
    ),

    stringsAsFactors = FALSE
  )
}

simulate_station_series <- function(n, base_mean, A, Sigma, burn = 100) {
  total_n <- n + burn
  y <- matrix(0, total_n, 2)

  for (t in 2:total_n) {
    y[t, ] <- as.numeric(A %*% y[t - 1, ] + rmvnorm2(Sigma))
  }

  y <- y[(burn + 1):total_n, , drop = FALSE]
  x <- matrix(rep(base_mean, each = n), nrow = n, ncol = 2) + y
  colnames(x) <- c("变量1", "变量2")
  x
}

generate_base_station_data <- function(n, station_info) {
  data_list <- vector("list", nrow(station_info))
  names(data_list) <- station_info$station

  for (i in seq_len(nrow(station_info))) {
    base_mean <- c(station_info$mean_var1[i], station_info$mean_var2[i])

    data_list[[i]] <- simulate_station_series(
      n = n,
      base_mean = base_mean,
      A = A,
      Sigma = Sigma
    )
  }

  data_list
}

apply_single_cp_change <- function(
    base_data,
    station_info,
    cp,
    target_station = "东四",
    change_type = c("normal", "anomaly"),
    shift_vec
) {
  change_type <- match.arg(change_type)

  data_list <- lapply(base_data, function(x) {
    out <- x
    colnames(out) <- colnames(x)
    out
  })

  target_idx <- which(station_info$station == target_station)

  if (length(target_idx) != 1) {
    stop("target_station must be unique in station_info$station.")
  }

  affected_stations <- if (change_type == "normal") {
    seq_along(data_list)
  } else {
    target_idx
  }

  for (i in affected_stations) {
    data_list[[i]][(cp + 1):nrow(data_list[[i]]), ] <-
      sweep(
        data_list[[i]][(cp + 1):nrow(data_list[[i]]), , drop = FALSE],
        2,
        shift_vec,
        "+"
      )
  }

  data_list
}

# ============================================================
# 5. 第一阶段：目标站点自身标准化 + SNSeg_Multi
# ============================================================

detect_candidate_cps <- function(
    target_data,
    confidence = 0.95,
    grid_size_scale = 0.05
) {
  out <- tryCatch(
    SNSeg_Multi(
      target_data,
      paras_to_test = "mean",
      confidence = confidence,
      grid_size_scale = grid_size_scale,
      plot_SN = FALSE,
      est_cp_loc = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(out) || is.null(out$est_cp) || length(out$est_cp) == 0) {
    return(integer(0))
  }

  cps <- out$est_cp
  if (is.list(cps)) cps <- unlist(cps)

  cps <- sort(unique(as.integer(round(cps))))
  cps <- cps[is.finite(cps)]
  cps <- cps[cps > 0 & cps < nrow(target_data)]

  cps
}

# ============================================================
# 6. 第二阶段：原始空间差值 + 正则化 Hotelling 型 H 统计量 + CBB 检验
# ============================================================

weighted_neighbors_equal <- function(
    data_list,
    station_info,
    target_station = "东四"
) {
  target_idx <- which(station_info$station == target_station)

  if (length(target_idx) != 1) {
    stop("target_station must be unique in station_info$station.")
  }

  ref_idx <- setdiff(seq_along(data_list), target_idx)

  ref_mat <- matrix(
    0,
    nrow = nrow(data_list[[target_idx]]),
    ncol = ncol(data_list[[target_idx]])
  )

  w <- 1 / length(ref_idx)

  for (j in ref_idx) {
    ref_mat <- ref_mat + w * data_list[[j]]
  }

  colnames(ref_mat) <- colnames(data_list[[target_idx]])
  ref_mat
}

build_residual_window_raw <- function(
    data_list_raw,
    station_info,
    tau,
    L = 50,
    target_station = "东四"
) {
  target_idx <- which(station_info$station == target_station)
  n_local <- nrow(data_list_raw[[target_idx]])

  s_idx <- tau - L + 1
  e_idx <- tau + L

  if (s_idx < 1 || e_idx > n_local) {
    return(NULL)
  }

  win_list <- lapply(data_list_raw, function(x) {
    x[s_idx:e_idx, , drop = FALSE]
  })

  ref_mat <- weighted_neighbors_equal(
    data_list = win_list,
    station_info = station_info,
    target_station = target_station
  )

  Z_raw <- win_list[[target_idx]] - ref_mat

  list(
    Z_raw = Z_raw,
    s_idx = s_idx,
    e_idx = e_idx
  )
}

test_one_candidate <- function(
    data_list_raw,
    station_info,
    tau,
    L = 50,
    M = 500,
    block_size = choose_cbb_block_size(L),
    lambda_frac = 0.05,
    target_station = "东四"
) {
  feat <- build_residual_window_raw(
    data_list_raw = data_list_raw,
    station_info = station_info,
    tau = tau,
    L = L,
    target_station = target_station
  )

  if (is.null(feat)) {
    return(NULL)
  }

  out <- block_boot_cbb_hotelling_cpp(
    Z = feat$Z_raw,
    L = L,
    M = M,
    block_size = block_size,
    lambda_frac = lambda_frac
  )

  data.frame(
    tau = tau,
    H_obs = as.numeric(out$H_obs),
    H_boot_q95 = as.numeric(out$H_boot_q95),
    p_value = as.numeric(out$p_value),
    stringsAsFactors = FALSE
  )
}

classify_candidates <- function(
    data_list_raw,
    station_info,
    candidate_cps,
    L = 50,
    M = 500,
    block_size = choose_cbb_block_size(L),
    lambda_frac = 0.05,
    alpha = 0.05,
    target_station = "东四"
) {
  if (length(candidate_cps) == 0) {
    return(data.frame())
  }

  res_list <- lapply(candidate_cps, function(tau) {
    test_one_candidate(
      data_list_raw = data_list_raw,
      station_info = station_info,
      tau = tau,
      L = L,
      M = M,
      block_size = block_size,
      lambda_frac = lambda_frac,
      target_station = target_station
    )
  })

  res_df <- do.call(rbind, res_list)

  if (is.null(res_df) || nrow(res_df) == 0) {
    return(data.frame())
  }

  res_df$pred_label <- ifelse(res_df$p_value < alpha, "异常", "正常")
  res_df <- res_df[order(res_df$tau), , drop = FALSE]
  rownames(res_df) <- NULL

  res_df
}

# ============================================================
# 7. 单次模拟结果提取
# ============================================================

extract_stage_results <- function(
    candidate_cps,
    candidate_df,
    cp,
    true_label,
    n,
    loc_tol = 20
) {
  stage1_detected <- length(candidate_cps) > 0

  if (stage1_detected) {
    stage1_eval <- calc_stage1_ari_hausdorff(
      n = n,
      true_cp = cp,
      detected_cps = candidate_cps
    )

    stage1_ari <- stage1_eval$ARI
    stage1_h <- stage1_eval$H
    stage1_nearest_tau <- candidate_cps[which.min(abs(candidate_cps - cp))]
    stage1_abs_error <- abs(stage1_nearest_tau - cp)
    stage1_loc_success <- any(abs(candidate_cps - cp) <= loc_tol)
  } else {
    stage1_ari <- NA_real_
    stage1_h <- NA_real_
    stage1_nearest_tau <- NA_integer_
    stage1_abs_error <- NA_real_
    stage1_loc_success <- FALSE
  }

  if (!is.null(candidate_df) && nrow(candidate_df) > 0) {
    near_df <- candidate_df[order(abs(candidate_df$tau - cp), candidate_df$p_value), , drop = FALSE]

    stage2_available <- TRUE
    pred_label <- near_df$pred_label[1]
    pred_tau <- near_df$tau[1]
    pred_abs_error <- abs(pred_tau - cp)
    pred_p <- near_df$p_value[1]
    pred_H <- near_df$H_obs[1]
    stage2_correct <- pred_label == true_label
  } else {
    stage2_available <- FALSE
    pred_label <- "未分类"
    pred_tau <- NA_integer_
    pred_abs_error <- NA_real_
    pred_p <- NA_real_
    pred_H <- NA_real_
    stage2_correct <- FALSE
  }

  data.frame(
    true_label = true_label,
    candidate_num = length(candidate_cps),

    stage1_detected = stage1_detected,
    stage1_loc_success = stage1_loc_success,
    stage1_nearest_tau = stage1_nearest_tau,
    stage1_abs_error = stage1_abs_error,
    stage1_ari = stage1_ari,
    stage1_h = stage1_h,

    stage2_available = stage2_available,
    pred_label = pred_label,
    pred_tau = pred_tau,
    pred_abs_error = pred_abs_error,
    p_value = pred_p,
    H_obs = pred_H,
    stage2_correct = stage2_correct,

    stringsAsFactors = FALSE
  )
}

run_one_simulation <- function(
    base_data,
    scenario,
    station_info,
    target_station = "东四",
    n = 500,
    cp = 250,
    L = 50,
    M = 500,
    block_size = choose_cbb_block_size(L),
    lambda_frac = 0.05,
    alpha = 0.05,
    confidence = 0.95,
    grid_size_scale = 0.05,
    loc_tol = 20,
    return_run_object = FALSE
) {
  data_raw <- apply_single_cp_change(
    base_data = base_data,
    station_info = station_info,
    cp = cp,
    target_station = target_station,
    change_type = scenario$change_type,
    shift_vec = scenario$shift_vec
  )

  target_idx <- which(station_info$station == target_station)

  # 第一阶段：只使用目标站点自身 Z-score 标准化后的双变量序列
  target_data_stage1 <- standardize_matrix(data_raw[[target_idx]])

  candidate_cps <- detect_candidate_cps(
    target_data = target_data_stage1,
    confidence = confidence,
    grid_size_scale = grid_size_scale
  )

  # 第二阶段：使用原始序列构造空间差值序列，并对每个候选变点单独进行 Hotelling 型检验
  candidate_df <- classify_candidates(
    data_list_raw = data_raw,
    station_info = station_info,
    candidate_cps = candidate_cps,
    L = L,
    M = M,
    block_size = block_size,
    lambda_frac = lambda_frac,
    alpha = alpha,
    target_station = target_station
  )

  true_label <- ifelse(scenario$change_type == "anomaly", "异常", "正常")

  result_row <- extract_stage_results(
    candidate_cps = candidate_cps,
    candidate_df = candidate_df,
    cp = cp,
    true_label = true_label,
    n = n,
    loc_tol = loc_tol
  )

  result_row$case_id <- scenario$case_id
  result_row$case_name <- scenario$case_name
  result_row$scenario_id <- scenario$scenario_id
  result_row$scenario_name <- scenario$scenario_name
  result_row$change_form <- scenario$change_form
  result_row$candidate_signature <- paste(candidate_cps, collapse = ",")

  result_row <- result_row[, c(
    "case_id",
    "case_name",
    "scenario_id",
    "scenario_name",
    "change_form",
    "true_label",
    "candidate_num",
    "candidate_signature",
    "stage1_detected",
    "stage1_loc_success",
    "stage1_nearest_tau",
    "stage1_abs_error",
    "stage1_ari",
    "stage1_h",
    "stage2_available",
    "pred_label",
    "pred_tau",
    "pred_abs_error",
    "p_value",
    "H_obs",
    "stage2_correct"
  )]

  if (!isTRUE(return_run_object)) {
    return(result_row)
  }

  list(
    data_raw = data_raw,
    candidate_cps = candidate_cps,
    candidate_df = candidate_df,
    result = result_row,
    scenario = scenario,
    cp = cp
  )
}

# ============================================================
# 8. 可视化：保存每个场景的最佳效果图
# ============================================================

score_plot_run <- function(run_obj, cp, loc_tol = 20) {
  result_row <- run_obj$result

  stage1_ok <- isTRUE(result_row$stage1_loc_success)
  stage2_ok <- isTRUE(result_row$stage2_correct)
  perfect_ok <- stage1_ok && stage2_ok

  score <- 0

  # 最高优先级：必须满足“第一阶段定位正确 + 第二阶段分类正确”。
  if (perfect_ok) {
    score <- score + 10000000
  }

  if (isTRUE(result_row$stage1_detected)) {
    score <- score + 100000
  } else {
    score <- score - 100000
  }

  if (stage1_ok) {
    score <- score + 200000
  } else {
    score <- score - 200000
  }

  if (isTRUE(result_row$stage2_available)) {
    score <- score + 50000
  }

  if (stage2_ok) {
    score <- score + 200000
  } else {
    score <- score - 200000
  }

  # 效果图优先选择定位误差更小、候选点更少的结果，避免图中过多无关候选线。
  if (!is.na(result_row$pred_abs_error)) {
    score <- score - 5000 * result_row$pred_abs_error
  }

  if (!is.na(result_row$candidate_num)) {
    score <- score - 2000 * max(0, result_row$candidate_num - 1)
  }

  # 分类置信度：异常场景 p 值越小越好；正常场景 p 值越大越好。
  if (!is.na(result_row$p_value)) {
    if (result_row$true_label == "异常") {
      score <- score + (1 - result_row$p_value) * 1000
    } else {
      score <- score + result_row$p_value * 1000
    }
  }

  score
}

plot_best_case <- function(
    run_obj,
    station_info,
    target_station = "东四",
    save_path
) {
  data_raw <- run_obj$data_raw
  candidate_df <- run_obj$candidate_df
  scenario <- run_obj$scenario
  cp <- run_obj$cp

  target_idx <- which(station_info$station == target_station)
  target <- data_raw[[target_idx]]

  ref <- weighted_neighbors_equal(
    data_list = data_raw,
    station_info = station_info,
    target_station = target_station
  )

  df_target <- data.frame(
    时间 = seq_len(nrow(target)),
    `变量1-目标站点` = target[, "变量1"],
    `变量2-目标站点` = target[, "变量2"],
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = -时间,
      names_to = "序列",
      values_to = "取值"
    )

  df_ref <- data.frame(
    时间 = seq_len(nrow(ref)),
    `变量1-临近站点加权值` = ref[, "变量1"],
    `变量2-临近站点加权值` = ref[, "变量2"],
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = -时间,
      names_to = "序列",
      values_to = "取值"
    )

  true_label <- ifelse(scenario$change_type == "anomaly", "异常", "正常")
  true_color <- ifelse(true_label == "异常", "#D62728", "#2CA02C")

  color_values <- c(
    "变量1-目标站点" = "#B2182B",
    "变量1-临近站点加权值" = "#EF8A62",
    "变量2-目标站点" = "#2166AC",
    "变量2-临近站点加权值" = "#67A9CF"
  )

  linetype_values <- c(
    "变量1-目标站点" = "solid",
    "变量1-临近站点加权值" = "longdash",
    "变量2-目标站点" = "solid",
    "变量2-临近站点加权值" = "longdash"
  )

  linewidth_values <- c(
    "变量1-目标站点" = 0.95,
    "变量1-临近站点加权值" = 0.72,
    "变量2-目标站点" = 0.95,
    "变量2-临近站点加权值" = 0.72
  )

  legend_order <- names(color_values)

  p <- ggplot() +
    geom_line(
      data = df_ref,
      aes(x = 时间, y = 取值, color = 序列, linetype = 序列, linewidth = 序列),
      alpha = 0.78
    ) +
    geom_line(
      data = df_target,
      aes(x = 时间, y = 取值, color = 序列, linetype = 序列, linewidth = 序列),
      alpha = 0.96
    ) +
    geom_vline(
      xintercept = cp,
      color = true_color,
      linetype = "dashed",
      linewidth = 1.10
    ) +
    scale_color_manual(values = color_values, breaks = legend_order) +
    scale_linetype_manual(values = linetype_values, breaks = legend_order) +
    scale_linewidth_manual(values = linewidth_values, breaks = legend_order) +
    labs(
      title = scenario$scenario_name,
      # subtitle = paste0(
      #   "目标站点：", target_station,
      #   "；虚线竖线为真实变点；实线竖线为第二阶段判定结果；绿色表示正常，红色表示异常"
      # ),
      x = "时间",
      y = "观测值",
      color = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    theme_minimal(base_size = 14) +
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
      panel.grid.major = element_line(color = "gray88")
    )

  if (!is.null(candidate_df) && nrow(candidate_df) > 0) {
    for (i in seq_len(nrow(candidate_df))) {
      pred_color <- ifelse(candidate_df$pred_label[i] == "异常", "#D62728", "#2CA02C")

      p <- p +
        geom_vline(
          xintercept = candidate_df$tau[i],
          color = pred_color,
          linetype = "solid",
          linewidth = 1.00
        )
    }
  }

  ggsave(
    filename = save_path,
    plot = p,
    width = 12.5,
    height = 7.2,
    dpi = 200,
    bg = "white"
  )

  message("图已保存：", normalizePath(save_path, winslash = "/", mustWork = FALSE))
  invisible(p)
}

# ============================================================
# 9. 场景设置与评价指标
# ============================================================

make_case_settings <- function(active_case_ids = c("case1", "case2")) {
  all_cases <- list(
    case1 = list(
      case_id = "case1",
      case_name = "（仅变量2有均值变化）",
      change_form = "仅变量2下降",
      shift_vec = delta_case1
    ),
    case2 = list(
      case_id = "case2",
      case_name = "（变量1和变量2同步地微小变化）",
      change_form = "变量1小幅上升、变量2小幅下降",
      shift_vec = delta_case2
    )
  )

  if (!all(active_case_ids %in% names(all_cases))) {
    stop("active_case_ids 只能包含 case1 和 case2。")
  }

  all_cases[active_case_ids]
}

build_scenarios_for_case <- function(case_obj) {
  list(
    list(
      scenario_id = paste0(case_obj$case_id, "_normal"),
      scenario_name = paste0("正常均值单变点",case_obj$case_name),
      case_id = case_obj$case_id,
      case_name = case_obj$case_name,
      change_form = case_obj$change_form,
      change_type = "normal",
      shift_vec = case_obj$shift_vec
    ),
    list(
      scenario_id = paste0(case_obj$case_id, "_anomaly"),
      scenario_name = paste0("异常均值单变点", case_obj$case_name),
      case_id = case_obj$case_id,
      case_name = case_obj$case_name,
      change_form = case_obj$change_form,
      change_type = "anomaly",
      shift_vec = case_obj$shift_vec
    )
  )
}

calc_metrics_table <- function(result_df) {
  result_df %>%
    group_by(case_id, case_name, scenario_name, true_label, change_form) %>%
    summarise(
      模拟次数 = n(),

      第一阶段变点检测率 = round(mean(stage1_detected), 3),
      第一阶段平均H距离 = round(mean_if_available(stage1_h), 2),
      第一阶段平均ARI = round(mean_if_available(stage1_ari), 3),

      第二阶段分类样本数 = sum(stage2_available),
      第二阶段分类准确率 = round(
        safe_div(
          sum(stage2_available & stage2_correct),
          sum(stage2_available)
        ),
        3
      ),

      .groups = "drop"
    ) %>%
    select(
      情形编号 = case_id,
      情形名称 = case_name,
      场景名称 = scenario_name,
      真实类别 = true_label,
      变化形式 = change_form,
      模拟次数,
      第一阶段变点检测率,
      第一阶段平均H距离,
      第一阶段平均ARI,
      第二阶段分类样本数,
      第二阶段分类准确率
    ) %>%
    arrange(情形编号, 真实类别)
}

# ============================================================
# 10. Monte Carlo 主函数
# ============================================================

split_iterations <- function(iter_vec, n_chunks) {
  n_chunks <- max(1L, min(as.integer(n_chunks), length(iter_vec)))
  split(iter_vec, cut(seq_along(iter_vec), breaks = n_chunks, labels = FALSE))
}

run_iterations_chunk_for_case <- function(
    iter_vec,
    case_obj,
    scenarios,
    station_info,
    n,
    cp,
    L,
    M,
    block_size,
    lambda_frac,
    alpha,
    confidence,
    grid_size_scale,
    loc_tol,
    target_station
) {
  # 每个并行任务内部连续运行一组模拟，并同步保留该组内每个场景的最佳图样本。
  result_rows <- vector("list", length(iter_vec) * length(scenarios))
  row_id <- 0L

  scenario_ids <- vapply(scenarios, function(x) x$scenario_id, character(1))
  best_runs <- setNames(vector("list", length(scenario_ids)), scenario_ids)
  best_scores <- setNames(as.list(rep(-Inf, length(scenario_ids))), scenario_ids)

  for (iter in iter_vec) {
    # 同一情形、同一轮模拟中，正常/异常共享同一套基础数据。
    base_data <- generate_base_station_data(
      n = n,
      station_info = station_info
    )

    for (scenario in scenarios) {
      run_obj <- run_one_simulation(
        base_data = base_data,
        scenario = scenario,
        station_info = station_info,
        target_station = target_station,
        n = n,
        cp = cp,
        L = L,
        M = M,
        block_size = block_size,
        lambda_frac = lambda_frac,
        alpha = alpha,
        confidence = confidence,
        grid_size_scale = grid_size_scale,
        loc_tol = loc_tol,
        return_run_object = TRUE
      )

      run_obj$result$iter <- iter
      row_id <- row_id + 1L
      result_rows[[row_id]] <- run_obj$result

      key <- scenario$scenario_id
      score <- score_plot_run(run_obj, cp = cp, loc_tol = loc_tol)

      if (is.null(best_runs[[key]]) || score > best_scores[[key]]) {
        best_scores[[key]] <- score
        best_runs[[key]] <- run_obj
      }
    }
  }

  result_rows <- result_rows[seq_len(row_id)]

  list(
    results = do.call(rbind, result_rows),
    best_runs = best_runs,
    best_scores = best_scores
  )
}

setup_parallel_cluster <- function(n_cores) {
  cl <- parallel::makeCluster(n_cores)

  objects_to_export <- c(
    "cpp_code",
    "A", "Sigma", "mu0",
    "delta_case1", "delta_case2",
    "get_default_cores", "choose_cbb_block_size",
    "rmvnorm2", "safe_div", "mean_if_available",
    "standardize_matrix",
    "get_segment_labels", "ari_score", "hausdorff_distance",
    "calc_stage1_ari_hausdorff",
    "make_station_info", "simulate_station_series",
    "generate_base_station_data", "apply_single_cp_change",
    "detect_candidate_cps",
    "weighted_neighbors_equal", "build_residual_window_raw",
    "test_one_candidate", "classify_candidates",
    "extract_stage_results", "run_one_simulation",
    "score_plot_run",
    "make_case_settings", "build_scenarios_for_case",
    "calc_metrics_table", "split_iterations", "run_iterations_chunk_for_case"
  )

  parallel::clusterExport(
    cl,
    varlist = objects_to_export,
    envir = .GlobalEnv
  )

  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(SNSeg)
      library(Rcpp)
      library(RcppArmadillo)
    })
    Rcpp::sourceCpp(code = cpp_code)
    NULL
  })

  parallel::clusterSetRNGStream(cl, iseed = 3047)

  cl
}

run_monte_carlo_single_cp_revised <- function(
    N_sim = 500,
    n = 500,
    cp = 250,
    L = 40,
    M = 500,
    block_size = choose_cbb_block_size(L),
    lambda_frac = 0.05,
    alpha = 0.05,
    confidence = 0.95,
    grid_size_scale = 0.05,
    target_station = "东四",
    active_case_ids = c("case1", "case2"),
    output_dir = "single_cp_anomaly_Hotelling_CBB_case12_inline_best_figures",
    save_plots = TRUE,
    loc_tol = 20,
    use_parallel = TRUE,
    n_cores = get_default_cores()
) {
  set.seed(3047)

  station_info <- make_station_info()
  case_settings <- make_case_settings(active_case_ids)

  cat("CBB 参数：M =", M, "；block_size =", block_size, "；局部窗口长度 =", 2 * L, "\n")

  all_results <- list()
  all_best_runs <- list()
  all_best_scores <- list()

  if (use_parallel && n_cores > 1) {
    cat("启用并行计算，核心数：", n_cores, "\n")
    cl <- setup_parallel_cluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    cl <- NULL
    cat("使用串行计算。\n")
  }

  for (case_idx in seq_along(case_settings)) {
    case_obj <- case_settings[[case_idx]]
    scenarios <- build_scenarios_for_case(case_obj)
    scenario_ids <- vapply(scenarios, function(x) x$scenario_id, character(1))

    cat("\n====================================================\n")
    cat("正在运行：", case_obj$case_name, "\n")
    cat("====================================================\n")

    if (!is.null(cl)) {
      iter_chunks <- split_iterations(seq_len(N_sim), n_chunks = n_cores * 4L)

      chunk_outputs <- parallel::parLapplyLB(
        cl,
        X = iter_chunks,
        fun = run_iterations_chunk_for_case,
        case_obj = case_obj,
        scenarios = scenarios,
        station_info = station_info,
        n = n,
        cp = cp,
        L = L,
        M = M,
        block_size = block_size,
        lambda_frac = lambda_frac,
        alpha = alpha,
        confidence = confidence,
        grid_size_scale = grid_size_scale,
        loc_tol = loc_tol,
        target_station = target_station
      )
    } else {
      chunk_outputs <- list(
        run_iterations_chunk_for_case(
          iter_vec = seq_len(N_sim),
          case_obj = case_obj,
          scenarios = scenarios,
          station_info = station_info,
          n = n,
          cp = cp,
          L = L,
          M = M,
          block_size = block_size,
          lambda_frac = lambda_frac,
          alpha = alpha,
          confidence = confidence,
          grid_size_scale = grid_size_scale,
          loc_tol = loc_tol,
          target_station = target_station
        )
      )
    }

    all_results[[case_obj$case_id]] <- do.call(
      rbind,
      lapply(chunk_outputs, function(x) x$results)
    )

    # 从各个并行任务的局部最佳结果中选出全局最佳图样本。
    for (sid in scenario_ids) {
      best_score <- -Inf
      best_run <- NULL

      for (chunk_out in chunk_outputs) {
        score <- chunk_out$best_scores[[sid]]
        if (!is.null(score) && is.finite(score) && score > best_score) {
          best_score <- score
          best_run <- chunk_out$best_runs[[sid]]
        }
      }

      all_best_scores[[sid]] <- best_score
      all_best_runs[[sid]] <- best_run
    }
  }

  result_df <- bind_rows(all_results)
  metrics_table <- calc_metrics_table(result_df)

  cat("\n\n============================================================\n")
  cat("表1  单变点异常识别性能评价结果：情形1和情形2，CBB + 正则化 Hotelling 空间差值检验\n")
  cat("说明：第一阶段使用目标站点 Z-score 标准化序列；第二阶段使用原始空间差值序列构造正则化 Hotelling 型 H 统计量。\n")
  cat("      第二阶段对每个候选变点单独检验，不合并候选点，不做 p 值校正；CBB 默认块长为 round(sqrt(2L))。\n")
  cat("      最佳效果图样本在 Monte Carlo 主循环中同步筛选，不额外重复模拟。\n")
  cat("============================================================\n\n")

  print(metrics_table, n = Inf, width = Inf)

  if (isTRUE(save_plots)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    cat("\n正在保存每个场景的最佳效果图...\n")

    for (key in names(all_best_runs)) {
      if (is.null(all_best_runs[[key]])) next

      rr <- all_best_runs[[key]]$result
      if (!(isTRUE(rr$stage1_loc_success) && isTRUE(rr$stage2_correct))) {
        message("警告：", key, " 在本次 Monte Carlo 模拟中未找到‘定位正确且分类正确’的样本，已保存退而求其次的最佳结果。")
      }

      save_path <- file.path(
        output_dir,
        paste0("Best_", key, ".png")
      )

      plot_best_case(
        run_obj = all_best_runs[[key]],
        station_info = station_info,
        target_station = target_station,
        save_path = save_path
      )
    }

    cat("\n图片已保存到文件夹：\n")
    cat(normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")
  }

  invisible(list(
    raw_results = result_df,
    metrics_table = metrics_table,
    station_info = station_info,
    best_plot_runs = all_best_runs,
    best_plot_scores = all_best_scores
  ))
}

# ============================================================
# 11. 运行
# ============================================================

mc_results <- run_monte_carlo_single_cp_revised(
  N_sim = 1000,
  n = n,
  cp = cp_true,
  L = L,
  M = M,
  block_size = block_size,
  lambda_frac = lambda_frac,
  alpha = alpha,
  confidence = confidence,
  grid_size_scale = grid_size_scale,
  target_station = target_station,
  active_case_ids = c("case1", "case2"),
  output_dir = output_dir,
  save_plots = TRUE,
  loc_tol = loc_tol,
  use_parallel = use_parallel,
  n_cores = n_cores
)
