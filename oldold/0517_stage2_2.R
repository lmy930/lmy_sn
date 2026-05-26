# ============================================================
# 双变量双变点异常识别模拟：修正版（仅模拟双变点情况1）
# 目标：
#   第一阶段：仅基于目标站点自身标准化序列进行双变量均值变点检测
#   第二阶段：基于各站点分别标准化后的 IDW 空间残差进行正常/异常识别
#
# 关键修正：
#   1. 第一阶段不使用所有站点合并标准化；
#   2. 同一轮模拟中，四类双变点场景共享同一套基础数据；
#   3. 四类场景中目标站点的两个变点完全相同，差别只在临近站点是否同步变化；
#   4. 第一阶段检测次数定义为 length(candidate_cps) > 0；
#   5. 第二阶段才引入 IDW 和临近站点信息。
#
# 双变点情况1：
#   第1个变点：仅变量2发生突变，delta1_case1 = c(0, -7.5)
#   第2个变点：两个变量均发生微小变化，delta2_case1 = c(-2.5, 2.0)
#
# 四类场景：
#   1. 正常-正常：两个变点均为正常变点，所有站点同步变化
#   2. 异常-异常：两个变点均为异常变点，仅目标站点变化
#   3. 正常-异常：第一个变点正常，第二个变点异常
#   4. 异常-正常：第一个变点异常，第二个变点正常
#
# 输出：
#   1. 每个场景保存一张最佳效果图；
#   2. 控制台打印评价指标表；
#   3. 不保存 CSV 文件。
# ============================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(SNSeg)
  library(Rcpp)
  library(RcppArmadillo)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

# ============================================================
# 0. 参数设置
# ============================================================

set.seed(3047)

B <- 1000
n <- 600
cps_true <- c(200, 400)

tol_best <- 2      # 仅用于最佳图选择
loc_tol <- 20      # 用于准确定位次数和额外异常误报判断，不用于定义是否检测到变点

grid_size_scale <- 0.05
confidence <- 0.95

# ----------------------------
# VAR(1) 模型参数
# ----------------------------

A <- matrix(
  c(
    0.30, 0.10,
    0.12, 0.25
  ),
  nrow = 2,
  byrow = TRUE
)

phi_ar <- diag(A)

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

# 仅使用双变点情况1的参数
delta1_case1 <- c(0, -7.5)
delta2_case1 <- c(-2.5, 2.0)

cols <- c(
  "变量1" = "#D55E00",
  "变量2" = "#0072B2"
)

# 第二阶段参数
L <- 50
M <- 500
block_size <- 15
boot_method <- "SBB"      # 可选："SBB", "MBB", "CBB"
lambda_frac <- 0.05
alpha <- 0.05

merge_close <- TRUE
merge_gap <- 30

target_station <- "东四"
output_dir <- "double_cp_case1_anomaly_figures_fixed"


# ============================================================
# 1. Rcpp：正则化 Hotelling 统计量 + 块抽样检验
# ============================================================

cpp_code <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

static inline mat estimate_reg_cov(const mat& Z_tilde, double lambda_frac) {
  mat S = cov(Z_tilde);
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

  rowvec delta = mean_after - mean_before;
  mat val = delta * Sinv * delta.t();

  return (static_cast<double>(L) / 2.0) * val(0, 0);
}

// [[Rcpp::export]]
List block_boot_hotelling_cpp(
    mat Z,
    int L,
    int M = 500,
    double expected_block_size = 15.0,
    double lambda_frac = 0.05,
    std::string boot_method = "SBB"
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

  for (int t = L; t < 2 * L; ++t) {
    Z_tilde.row(t) -= mean_after;
  }

  mat Sigma_reg = estimate_reg_cov(Z_tilde, lambda_frac);

  mat Sinv;
  bool ok = inv_sympd(Sinv, Sigma_reg);

  if (!ok) {
    Sinv = pinv(Sigma_reg);
  }

  double S_obs = calc_hotelling_stat(Z, L, Sinv);

  int b = std::max(1, static_cast<int>(std::round(expected_block_size)));
  b = std::min(b, n);

  double p_geom = 1.0 / std::max(2.0, expected_block_size);

  int exceed_count = 0;
  vec S_boot(M, fill::zeros);

  for (int m = 0; m < M; ++m) {
    mat Z_star(n, d, fill::zeros);

    int t = 0;

    while (t < n) {
      int block_length;
      int start_idx;

      if (boot_method == "SBB") {
        block_length = static_cast<int>(R::rgeom(p_geom)) + 1;
        start_idx = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(n))));
      } else if (boot_method == "MBB") {
        block_length = b;
        int max_start = std::max(0, n - b);
        start_idx = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(max_start + 1))));
      } else if (boot_method == "CBB") {
        block_length = b;
        start_idx = static_cast<int>(std::floor(R::runif(0.0, static_cast<double>(n))));
      } else {
        stop("boot_method must be one of: SBB, MBB, CBB.");
      }

      for (int i = 0; i < block_length && t < n; ++i) {
        int idx;

        if (boot_method == "MBB") {
          idx = start_idx + i;
          if (idx >= n) break;
        } else {
          idx = (start_idx + i) % n;
        }

        Z_star.row(t) = Z_tilde.row(idx);
        t++;
      }
    }

    double S_star = calc_hotelling_stat(Z_star, L, Sinv);
    S_boot(m) = S_star;

    if (S_star >= S_obs) {
      exceed_count++;
    }
  }

  double p_value = static_cast<double>(exceed_count + 1) / static_cast<double>(M + 1);

  vec S_sorted = sort(S_boot);
  int idx95 = std::min(M - 1, std::max(0, static_cast<int>(std::floor(0.95 * M)) - 1));

  return List::create(
    _["p_value"] = p_value,
    _["S_obs"] = S_obs,
    _["S_boot_q95"] = S_sorted(idx95)
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


haversine_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  R <- 6371.0

  dlon <- (lon2 - lon1) * rad
  dlat <- (lat2 - lat1) * rad

  a <- sin(dlat / 2)^2 +
    cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2

  R * 2 * atan2(sqrt(a), sqrt(1 - a))
}


safe_div <- function(a, b) {
  ifelse(b == 0, NA_real_, a / b)
}


mean_if_available <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(NA_real_)
  }

  mean(x)
}


standardize_matrix <- function(x) {
  center <- colMeans(x)
  scale_val <- apply(x, 2, sd)
  scale_val[scale_val < 1e-8] <- 1

  out <- sweep(sweep(x, 2, center, "-"), 2, scale_val, "/")
  colnames(out) <- colnames(x)

  out
}


standardize_each_station <- function(data_list) {
  lapply(data_list, standardize_matrix)
}


weighted_neighbors_idw <- function(data_list, station_info, target_station = "东四", power = 2) {
  target_idx <- which(station_info$站点 == target_station)

  if (length(target_idx) != 1) {
    stop("target_station 在 station_info$站点 中必须唯一存在。")
  }

  ref_idx <- setdiff(seq_along(data_list), target_idx)

  lon_t <- station_info$经度[target_idx]
  lat_t <- station_info$纬度[target_idx]

  dists <- sapply(ref_idx, function(i) {
    haversine_dist(
      lon_t, lat_t,
      station_info$经度[i], station_info$纬度[i]
    )
  })

  dists[dists < 1e-6] <- 1e-6

  weights <- (1 / dists^power) / sum(1 / dists^power)

  ref_mat <- matrix(
    0,
    nrow = nrow(data_list[[target_idx]]),
    ncol = ncol(data_list[[target_idx]])
  )

  for (j in seq_along(ref_idx)) {
    ref_mat <- ref_mat + data_list[[ref_idx[j]]] * weights[j]
  }

  colnames(ref_mat) <- colnames(data_list[[target_idx]])

  ref_mat
}


# ============================================================
# 3. 第一阶段评价函数：检测次数、ARI、H距离
# ============================================================

get_segment_labels <- function(n, cps) {
  labels <- rep(1L, n)

  if (length(cps) == 0) {
    return(labels)
  }

  cps <- sort(unique(as.integer(cps)))
  cps <- cps[cps > 0 & cps < n]

  if (length(cps) == 0) {
    return(labels)
  }

  for (i in seq_along(cps)) {
    labels[(cps[i] + 1):n] <- i + 1L
  }

  labels
}


ari_score <- function(label_true, label_est) {
  tab <- table(label_true, label_est)

  comb2 <- function(x) x * (x - 1) / 2

  n_total <- sum(tab)
  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total <- comb2(n_total)

  expected <- sum_i * sum_j / total
  max_index <- 0.5 * (sum_i + sum_j)

  if (max_index == expected) {
    return(1)
  }

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


calc_stage1_ari_hausdorff <- function(n, true_cps, detected_cps) {
  detected_cps <- detected_cps[detected_cps > 0 & detected_cps < n]

  if (length(detected_cps) == 0) {
    return(list(ARI = NA_real_, H = NA_real_))
  }

  true_labels <- get_segment_labels(n, true_cps)
  est_labels <- get_segment_labels(n, detected_cps)

  list(
    ARI = ari_score(true_labels, est_labels),
    H = hausdorff_distance(detected_cps, true_cps)
  )
}


match_detected_to_true <- function(true_cps, detected_cps, loc_tol = 20) {
  true_cps <- as.integer(true_cps)
  detected_cps <- sort(unique(as.integer(detected_cps)))

  out <- data.frame(
    true_cp = true_cps,
    matched = FALSE,
    matched_tau = NA_integer_,
    abs_error = NA_real_
  )

  if (length(detected_cps) == 0) {
    return(out)
  }

  used <- rep(FALSE, length(detected_cps))

  pair_df <- expand.grid(
    true_id = seq_along(true_cps),
    cand_id = seq_along(detected_cps)
  )

  pair_df$dist <- abs(true_cps[pair_df$true_id] - detected_cps[pair_df$cand_id])
  pair_df <- pair_df %>% arrange(dist)

  for (i in seq_len(nrow(pair_df))) {
    true_id <- pair_df$true_id[i]
    cand_id <- pair_df$cand_id[i]

    if (!out$matched[true_id] && !used[cand_id] && pair_df$dist[i] <= loc_tol) {
      out$matched[true_id] <- TRUE
      out$matched_tau[true_id] <- detected_cps[cand_id]
      out$abs_error[true_id] <- pair_df$dist[i]
      used[cand_id] <- TRUE
    }
  }

  out
}


# ============================================================
# 4. 站点信息与基础数据生成
# ============================================================

make_station_info <- function() {
  data.frame(
    站点 = c("东四", "天坛", "官园", "万寿西宫", "奥体中心", "农展馆"),
    经度 = c(116.417, 116.407, 116.339, 116.352, 116.397, 116.461),
    纬度 = c(39.929, 39.886, 39.929, 39.878, 39.982, 39.937),

    变量1均值 = c(
      mu0[1],
      mu0[1] - 1.0,
      mu0[1] + 1.5,
      mu0[1] - 2.0,
      mu0[1] + 2.0,
      mu0[1] - 0.5
    ),

    变量2均值 = c(
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


simulate_station_series <- function(
    n,
    base_mean,
    A,
    Sigma,
    burn = 100
) {
  total_n <- n + burn
  y <- matrix(0, total_n, 2)

  for (t in 2:total_n) {
    y[t, ] <- A %*% y[t - 1, ] + rmvnorm2(Sigma)
  }

  y <- y[(burn + 1):total_n, ]

  x <- matrix(rep(base_mean, each = n), ncol = 2) + y
  colnames(x) <- c("变量1", "变量2")

  x
}


generate_base_station_data <- function(
    n,
    station_info
) {
  data_list <- vector("list", nrow(station_info))
  names(data_list) <- station_info$站点

  for (i in seq_len(nrow(station_info))) {
    base_mean <- c(station_info$变量1均值[i], station_info$变量2均值[i])

    data_list[[i]] <- simulate_station_series(
      n = n,
      base_mean = base_mean,
      A = A,
      Sigma = Sigma
    )
  }

  data_list
}


apply_double_cp_change <- function(
    base_data,
    station_info,
    cps,
    target_station = "东四",
    change_types = c("normal", "normal"),
    shift_vecs = list(delta1_case1, delta2_case1)
) {
  if (length(cps) != length(change_types) || length(cps) != length(shift_vecs)) {
    stop("cps、change_types 和 shift_vecs 的长度必须一致。")
  }

  data_list <- lapply(base_data, function(x) {
    out <- x
    colnames(out) <- colnames(x)
    out
  })

  target_idx <- which(station_info$站点 == target_station)

  if (length(target_idx) != 1) {
    stop("target_station 在 station_info$站点 中必须唯一存在。")
  }

  for (k in seq_along(cps)) {
    change_type <- match.arg(change_types[k], choices = c("normal", "anomaly"))

    affected_stations <- if (change_type == "normal") {
      seq_along(data_list)
    } else {
      target_idx
    }

    for (i in affected_stations) {
      data_list[[i]][(cps[k] + 1):nrow(data_list[[i]]), ] <-
        sweep(
          data_list[[i]][(cps[k] + 1):nrow(data_list[[i]]), , drop = FALSE],
          2,
          shift_vecs[[k]],
          "+"
        )
    }
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

  cps <- sort(unique(as.integer(round(out$est_cp))))
  cps <- cps[is.finite(cps)]

  cps
}


# ============================================================
# 6. 第二阶段：各站点分别标准化 + IDW 残差 + Hotelling 块抽样
# ============================================================

build_residual_window <- function(
    data_list_std,
    station_info,
    tau,
    L = 50,
    target_station = "东四",
    power = 2
) {
  target_idx <- which(station_info$站点 == target_station)

  n_local <- nrow(data_list_std[[target_idx]])

  s_idx <- tau - L + 1
  e_idx <- tau + L

  if (s_idx < 1 || e_idx > n_local) {
    return(NULL)
  }

  win_list <- lapply(data_list_std, function(x) {
    x[s_idx:e_idx, , drop = FALSE]
  })

  ref_mat <- weighted_neighbors_idw(
    data_list = win_list,
    station_info = station_info,
    target_station = target_station,
    power = power
  )

  Z <- win_list[[target_idx]] - ref_mat

  list(
    Z = Z,
    s_idx = s_idx,
    e_idx = e_idx
  )
}


test_one_candidate <- function(
    data_list_std,
    station_info,
    tau,
    L = 50,
    M = 500,
    block_size = 15,
    boot_method = "SBB",
    lambda_frac = 0.05,
    target_station = "东四",
    power = 2
) {
  feat <- build_residual_window(
    data_list_std = data_list_std,
    station_info = station_info,
    tau = tau,
    L = L,
    target_station = target_station,
    power = power
  )

  if (is.null(feat)) {
    return(NULL)
  }

  out <- block_boot_hotelling_cpp(
    Z = feat$Z,
    L = L,
    M = M,
    expected_block_size = block_size,
    lambda_frac = lambda_frac,
    boot_method = boot_method
  )

  data.frame(
    tau = tau,
    S_obs = as.numeric(out$S_obs),
    S_boot_q95 = as.numeric(out$S_boot_q95),
    p_value = as.numeric(out$p_value)
  )
}


merge_close_candidates <- function(candidate_df, merge_gap = 30) {
  if (is.null(candidate_df) || nrow(candidate_df) == 0) {
    return(candidate_df)
  }

  candidate_df <- candidate_df %>%
    arrange(tau)

  group_id <- cumsum(c(TRUE, diff(candidate_df$tau) > merge_gap))
  candidate_df$group_id <- group_id

  candidate_df %>%
    group_by(group_id) %>%
    slice_min(order_by = p_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-group_id) %>%
    arrange(tau)
}


classify_candidates <- function(
    data_list_std,
    station_info,
    candidate_cps,
    L = 50,
    M = 500,
    block_size = 15,
    boot_method = "SBB",
    lambda_frac = 0.05,
    alpha = 0.05,
    p_adjust_method = "BH",
    merge_close = TRUE,
    merge_gap = 30,
    target_station = "东四",
    power = 2
) {
  if (length(candidate_cps) == 0) {
    return(data.frame())
  }

  res_list <- lapply(candidate_cps, function(tau) {
    test_one_candidate(
      data_list_std = data_list_std,
      station_info = station_info,
      tau = tau,
      L = L,
      M = M,
      block_size = block_size,
      boot_method = boot_method,
      lambda_frac = lambda_frac,
      target_station = target_station,
      power = power
    )
  })

  res_df <- bind_rows(res_list)

  if (nrow(res_df) == 0) {
    return(data.frame())
  }

  if (merge_close) {
    res_df <- merge_close_candidates(res_df, merge_gap = merge_gap)
  }

  res_df %>%
    mutate(
      p_adj = p.adjust(p_value, method = p_adjust_method),
      pred_label = ifelse(p_adj < alpha, "异常", "正常"),
      line_color = ifelse(pred_label == "异常", "#d62728", "#2ca02c")
    ) %>%
    arrange(tau)
}


# ============================================================
# 7. 单次模拟结果提取
# ============================================================

extract_stage_results <- function(
    candidate_cps,
    candidate_df,
    cps,
    true_labels,
    n,
    loc_tol = 20
) {
  stage1_detected <- length(candidate_cps) > 0

  if (stage1_detected) {
    match_df <- match_detected_to_true(
      true_cps = cps,
      detected_cps = candidate_cps,
      loc_tol = loc_tol
    )

    stage1_loc_hit_num <- sum(match_df$matched)
    stage1_complete_loc_success <- all(match_df$matched)
    stage1_abs_error <- mean_if_available(match_df$abs_error[match_df$matched])

    stage1_eval <- calc_stage1_ari_hausdorff(
      n = n,
      true_cps = cps,
      detected_cps = candidate_cps
    )

    stage1_ari <- stage1_eval$ARI
    stage1_h <- stage1_eval$H
  } else {
    match_df <- data.frame(
      true_cp = cps,
      matched = FALSE,
      matched_tau = NA_integer_,
      abs_error = NA_real_
    )

    stage1_loc_hit_num <- 0
    stage1_complete_loc_success <- FALSE
    stage1_abs_error <- NA_real_
    stage1_ari <- NA_real_
    stage1_h <- NA_real_
  }

  stage2_available_num <- 0
  stage2_correct_num <- 0

  TP <- 0
  FP <- 0
  TN <- 0
  FN <- 0

  p_values <- c()
  p_adj_values <- c()

  if (stage1_detected && !is.null(candidate_df) && nrow(candidate_df) > 0) {
    for (k in seq_along(cps)) {
      if (!match_df$matched[k]) {
        next
      }

      tau_k <- match_df$matched_tau[k]

      near_df <- candidate_df %>%
        mutate(dist_to_tau = abs(tau - tau_k)) %>%
        arrange(dist_to_tau, p_adj)

      if (nrow(near_df) > 0) {
        pred_label <- near_df$pred_label[1]
        true_label <- true_labels[k]

        stage2_available_num <- stage2_available_num + 1
        p_values <- c(p_values, near_df$p_value[1])
        p_adj_values <- c(p_adj_values, near_df$p_adj[1])

        if (pred_label == true_label) {
          stage2_correct_num <- stage2_correct_num + 1
        }

        if (true_label == "异常" && pred_label == "异常") {
          TP <- TP + 1
        } else if (true_label == "正常" && pred_label == "异常") {
          FP <- FP + 1
        } else if (true_label == "正常" && pred_label == "正常") {
          TN <- TN + 1
        } else if (true_label == "异常" && pred_label == "正常") {
          FN <- FN + 1
        }
      }
    }
  }

  stage2_available <- stage2_available_num > 0
  stage2_complete_available <- stage2_available_num == length(cps)
  stage2_correct <- stage2_correct_num == stage2_available_num && stage2_available_num > 0
  stage2_complete_correct <- stage2_correct_num == length(cps)

  extra_anomaly_fp <- if (!is.null(candidate_df) && nrow(candidate_df) > 0) {
    candidate_df %>%
      filter(
        pred_label == "异常",
        sapply(tau, function(x) min(abs(x - cps)) > loc_tol)
      ) %>%
      nrow()
  } else {
    0
  }

  data.frame(
    true_label = paste(true_labels, collapse = "-"),
    candidate_num = length(candidate_cps),

    stage1_detected = stage1_detected,
    stage1_loc_hit_num = stage1_loc_hit_num,
    stage1_complete_loc_success = stage1_complete_loc_success,
    stage1_abs_error = stage1_abs_error,
    stage1_ari = stage1_ari,
    stage1_h = stage1_h,

    stage2_available = stage2_available,
    stage2_complete_available = stage2_complete_available,
    stage2_available_num = stage2_available_num,
    stage2_correct = stage2_correct,
    stage2_complete_correct = stage2_complete_correct,
    stage2_correct_num = stage2_correct_num,

    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,

    p_value = mean_if_available(p_values),
    p_adj = mean_if_available(p_adj_values),
    extra_anomaly_fp = extra_anomaly_fp
  )
}


score_best_run <- function(result_row) {
  score <- 0

  if (isTRUE(result_row$stage1_detected)) {
    score <- score + 100000
  } else {
    score <- score - 100000
  }

  score <- score + 30000 * result_row$stage1_loc_hit_num

  if (isTRUE(result_row$stage1_complete_loc_success)) {
    score <- score + 50000
  }

  if (!is.na(result_row$stage1_abs_error)) {
    score <- score - 1000 * result_row$stage1_abs_error
  } else {
    score <- score - 50000
  }

  if (isTRUE(result_row$stage2_complete_available)) {
    score <- score + 30000
  } else if (isTRUE(result_row$stage2_available)) {
    score <- score + 10000
  }

  score <- score + 20000 * result_row$stage2_correct_num

  if (isTRUE(result_row$stage2_complete_correct)) {
    score <- score + 50000
  }

  score <- score - 5000 * result_row$extra_anomaly_fp

  if (!is.na(result_row$p_adj)) {
    score <- score + (1 - result_row$p_adj) * 100
  }

  score
}


run_one_simulation <- function(
    base_data,
    scenario,
    station_info,
    target_station = "东四",
    n = 600,
    cps = c(200, 400),
    L = 50,
    M = 500,
    block_size = 15,
    boot_method = "SBB",
    lambda_frac = 0.05,
    alpha = 0.05,
    loc_tol = 20,
    grid_size_scale = 0.05,
    merge_close = TRUE,
    merge_gap = 30
) {
  data_raw <- apply_double_cp_change(
    base_data = base_data,
    station_info = station_info,
    cps = cps,
    target_station = target_station,
    change_types = scenario$change_types,
    shift_vecs = scenario$shift_vecs
  )

  target_idx <- which(station_info$站点 == target_station)

  # 第一阶段：只使用目标站点自身标准化后的双变量序列
  target_data_stage1 <- standardize_matrix(data_raw[[target_idx]])

  candidate_cps <- detect_candidate_cps(
    target_data = target_data_stage1,
    confidence = confidence,
    grid_size_scale = grid_size_scale
  )

  # 第二阶段：各站点分别标准化后构造 IDW 空间残差
  data_std <- standardize_each_station(data_raw)

  candidate_df <- classify_candidates(
    data_list_std = data_std,
    station_info = station_info,
    candidate_cps = candidate_cps,
    L = L,
    M = M,
    block_size = block_size,
    boot_method = boot_method,
    lambda_frac = lambda_frac,
    alpha = alpha,
    p_adjust_method = "BH",
    merge_close = merge_close,
    merge_gap = merge_gap,
    target_station = target_station,
    power = 2
  )

  true_labels <- ifelse(scenario$change_types == "anomaly", "异常", "正常")

  result_row <- extract_stage_results(
    candidate_cps = candidate_cps,
    candidate_df = candidate_df,
    cps = cps,
    true_labels = true_labels,
    n = n,
    loc_tol = loc_tol
  )

  result_row$scenario_id <- scenario$scenario_id
  result_row$scenario_name <- scenario$scenario_name
  result_row$change_form <- scenario$change_form

  result_row <- result_row %>%
    select(
      scenario_id,
      scenario_name,
      change_form,
      everything()
    )

  result_row$best_score <- score_best_run(result_row)

  list(
    data_raw = data_raw,
    data_std = data_std,
    candidate_cps = candidate_cps,
    candidate_df = candidate_df,
    result = result_row,
    scenario = scenario,
    cps = cps
  )
}


# ============================================================
# 8. 绘图：原始数据
# ============================================================

plot_best_case <- function(
    run_obj,
    station_info,
    target_station = "东四",
    save_path
) {
  data_raw <- run_obj$data_raw
  candidate_df <- run_obj$candidate_df
  scenario <- run_obj$scenario
  cps <- run_obj$cps
  true_labels <- ifelse(scenario$change_types == "anomaly", "异常", "正常")

  target_idx <- which(station_info$站点 == target_station)

  target <- data_raw[[target_idx]]

  ref <- weighted_neighbors_idw(
    data_list = data_raw,
    station_info = station_info,
    target_station = target_station,
    power = 2
  )

  df_ref <- data.frame(
    时间 = seq_len(nrow(target)),
    `变量1-临近站点加权参考` = ref[, "变量1"],
    `变量2-临近站点加权参考` = ref[, "变量2"],
    check.names = FALSE
  ) %>%
    pivot_longer(
      cols = -时间,
      names_to = "序列",
      values_to = "取值"
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

  color_values <- c(
    "变量1-目标站点" = "#B2182B",
    "变量1-临近站点加权参考" = "#EF8A62",
    "变量2-目标站点" = "#2166AC",
    "变量2-临近站点加权参考" = "#67A9CF"
  )

  linetype_values <- c(
    "变量1-目标站点" = "solid",
    "变量1-临近站点加权参考" = "longdash",
    "变量2-目标站点" = "solid",
    "变量2-临近站点加权参考" = "longdash"
  )

  linewidth_values <- c(
    "变量1-目标站点" = 0.95,
    "变量1-临近站点加权参考" = 0.72,
    "变量2-目标站点" = 0.95,
    "变量2-临近站点加权参考" = 0.72
  )

  legend_order <- c(
    "变量1-目标站点",
    "变量1-临近站点加权参考",
    "变量2-目标站点",
    "变量2-临近站点加权参考"
  )

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
    )

  for (k in seq_along(cps)) {
    true_color <- ifelse(true_labels[k] == "异常", "#D62728", "#2CA02C")

    p <- p +
      geom_vline(
        xintercept = cps[k],
        color = true_color,
        linetype = "dashed",
        linewidth = 1.15
      )
  }

  p <- p +
    scale_color_manual(values = color_values, breaks = legend_order) +
    scale_linetype_manual(values = linetype_values, breaks = legend_order) +
    scale_linewidth_manual(values = linewidth_values, breaks = legend_order) +
    labs(
      title = scenario$scenario_name,
      subtitle = paste0(
        "目标站点：", target_station,
        "；虚线竖线为真实变点，实线竖线为算法判定结果；绿色表示正常变点，红色表示异常变点"
      ),
      x = "时间",
      y = "变量取值",
      color = NULL,
      linetype = NULL,
      linewidth = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = unit(2.6, "cm"),
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
          linewidth = 1.05
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

make_double_cp_case1_scenarios <- function() {
  list(
    list(
      scenario_id = "case1_normal_normal",
      scenario_name = "正常双变点：第一个变点正常，第二个变点正常",
      change_form = "正常-正常",
      change_types = c("normal", "normal"),
      shift_vecs = list(delta1_case1, delta2_case1)
    ),
    list(
      scenario_id = "case1_anomaly_anomaly",
      scenario_name = "异常双变点：第一个变点异常，第二个变点异常",
      change_form = "异常-异常",
      change_types = c("anomaly", "anomaly"),
      shift_vecs = list(delta1_case1, delta2_case1)
    ),
    list(
      scenario_id = "case1_normal_anomaly",
      scenario_name = "混合双变点：第一个变点正常，第二个变点异常",
      change_form = "正常-异常",
      change_types = c("normal", "anomaly"),
      shift_vecs = list(delta1_case1, delta2_case1)
    ),
    list(
      scenario_id = "case1_anomaly_normal",
      scenario_name = "混合双变点：第一个变点异常，第二个变点正常",
      change_form = "异常-正常",
      change_types = c("anomaly", "normal"),
      shift_vecs = list(delta1_case1, delta2_case1)
    )
  )
}


calc_metrics_table <- function(result_df) {
  result_df %>%
    group_by(scenario_id, scenario_name, true_label, change_form) %>%
    summarise(
      模拟次数 = n(),

      第一阶段检测次数 = sum(stage1_detected),
      第一阶段检测率 = round(mean(stage1_detected), 3),
      第一阶段平均检出变点数 = round(mean(stage1_loc_hit_num), 3),
      第一阶段完整定位次数 = sum(stage1_complete_loc_success),
      第一阶段完整定位率 = round(mean(stage1_complete_loc_success), 3),
      第一阶段平均定位误差 = round(mean_if_available(stage1_abs_error[stage1_detected]), 2),
      第一阶段平均ARI = round(mean_if_available(stage1_ari[stage1_detected]), 3),
      第一阶段平均H距离 = round(mean_if_available(stage1_h[stage1_detected]), 2),

      第二阶段分类样本数 = sum(stage2_available_num[stage1_detected]),
      第二阶段正确分类次数 = sum(stage2_correct_num[stage1_detected]),
      第二阶段完整分类次数 = sum(stage2_complete_correct[stage1_detected]),
      第二阶段完整分类率 = round(mean(stage2_complete_correct[stage1_detected]), 3),

      第二阶段准确率 = round(
        safe_div(
          sum(TP[stage1_detected] + TN[stage1_detected]),
          sum(TP[stage1_detected] + FP[stage1_detected] + TN[stage1_detected] + FN[stage1_detected])
        ),
        3
      ),
      第二阶段精确率 = round(
        safe_div(
          sum(TP[stage1_detected]),
          sum(TP[stage1_detected] + FP[stage1_detected])
        ),
        3
      ),
      第二阶段召回率 = round(
        safe_div(
          sum(TP[stage1_detected]),
          sum(TP[stage1_detected] + FN[stage1_detected])
        ),
        3
      ),
      第二阶段F1分数 = round(
        safe_div(
          2 *
            safe_div(sum(TP[stage1_detected]), sum(TP[stage1_detected] + FP[stage1_detected])) *
            safe_div(sum(TP[stage1_detected]), sum(TP[stage1_detected] + FN[stage1_detected])),
          safe_div(sum(TP[stage1_detected]), sum(TP[stage1_detected] + FP[stage1_detected])) +
            safe_div(sum(TP[stage1_detected]), sum(TP[stage1_detected] + FN[stage1_detected]))
        ),
        3
      ),
      平均p值 = round(
        mean_if_available(p_value[stage1_detected]),
        4
      ),
      平均校正p值 = round(
        mean_if_available(p_adj[stage1_detected]),
        4
      ),
      额外异常误报平均数 = round(mean(extra_anomaly_fp[stage1_detected], na.rm = TRUE), 3),

      .groups = "drop"
    ) %>%
    select(
      场景编号 = scenario_id,
      场景名称 = scenario_name,
      真实类别组合 = true_label,
      变化形式 = change_form,
      everything()
    ) %>%
    arrange(场景编号)
}


check_stage1_consistency <- function(result_df) {
  result_df %>%
    group_by(iter) %>%
    summarise(
      场景数 = n(),
      第一阶段候选结果种类数 = n_distinct(candidate_signature),
      第一阶段是否一致 = 第一阶段候选结果种类数 == 1,
      .groups = "drop"
    ) %>%
    summarise(
      检查轮数 = n(),
      第一阶段一致次数 = sum(第一阶段是否一致),
      第一阶段一致比例 = round(mean(第一阶段是否一致), 3),
      .groups = "drop"
    )
}


# ============================================================
# 10. Monte Carlo 主函数
# ============================================================

run_monte_carlo_double_cp_case1_fixed <- function(
    N_sim = 1000,
    n = 600,
    cps = c(200, 400),
    L = 50,
    M = 500,
    block_size = 15,
    boot_method = "SBB",
    lambda_frac = 0.05,
    alpha = 0.05,
    loc_tol = 20,
    tol_best = 2,
    grid_size_scale = 0.05,
    merge_close = TRUE,
    merge_gap = 30,
    target_station = "东四",
    output_dir = "double_cp_case1_anomaly_figures_fixed"
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  station_info <- make_station_info()
  scenarios <- make_double_cp_case1_scenarios()

  all_result_rows <- list()
  all_run_objects <- list()

  cat("\n====================================================\n")
  cat("正在运行：双变点情况1异常识别模拟\n")
  cat("====================================================\n")

  pb <- txtProgressBar(min = 0, max = N_sim, style = 3)

  for (iter in seq_len(N_sim)) {
    set.seed(3047 + iter)

    # 同一轮模拟中，四类场景共享同一套基础数据
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
        cps = cps,
        L = L,
        M = M,
        block_size = block_size,
        boot_method = boot_method,
        lambda_frac = lambda_frac,
        alpha = alpha,
        loc_tol = loc_tol,
        grid_size_scale = grid_size_scale,
        merge_close = merge_close,
        merge_gap = merge_gap
      )

      run_obj$result$iter <- iter
      run_obj$result$candidate_signature <- paste(run_obj$candidate_cps, collapse = ",")

      key <- scenario$scenario_id

      if (is.null(all_result_rows[[key]])) {
        all_result_rows[[key]] <- list()
      }

      if (is.null(all_run_objects[[key]])) {
        all_run_objects[[key]] <- list()
      }

      all_result_rows[[key]][[iter]] <- run_obj$result
      all_run_objects[[key]][[iter]] <- run_obj
    }

    setTxtProgressBar(pb, iter)
  }

  close(pb)

  result_df <- bind_rows(lapply(all_result_rows, bind_rows))

  # 保存每个场景最佳图
  best_runs <- list()

  scenario_keys <- names(all_result_rows)

  for (key in scenario_keys) {
    scenario_df <- bind_rows(all_result_rows[[key]])

    scenario_df <- scenario_df %>%
      mutate(
        best_score_final = best_score +
          ifelse(stage1_complete_loc_success, 50000, 0) +
          ifelse(!is.na(stage1_abs_error) & stage1_abs_error <= tol_best, 20000, 0)
      )

    best_iter <- scenario_df$iter[which.max(scenario_df$best_score_final)]
    best_run <- all_run_objects[[key]][[best_iter]]
    best_runs[[key]] <- best_run

    save_path <- file.path(
      output_dir,
      paste0("Best_", key, ".png")
    )

    plot_best_case(
      run_obj = best_run,
      station_info = station_info,
      target_station = target_station,
      save_path = save_path
    )
  }

  metrics_table <- calc_metrics_table(result_df)
  consistency_table <- check_stage1_consistency(result_df)

  cat("\n\n============================================================\n")
  cat("表1  修正后：双变点情况1异常识别性能评价结果\n")
  cat("说明：第一阶段检测次数按 length(candidate_cps) > 0 统计；ARI 和 H 距离仅在第一阶段检测到变点的模拟中计算。\n")
  cat("      第一阶段仅使用目标站点自身标准化数据；第二阶段才使用各站点分别标准化后的 IDW 空间残差。\n")
  cat("============================================================\n\n")

  print(metrics_table, n = Inf, width = Inf)

  cat("\n\n============================================================\n")
  cat("表2  四类双变点场景第一阶段结果一致性检查\n")
  cat("说明：同一轮模拟中四类场景共享同一套基础数据，目标站点变化相同，因此第一阶段应一致。\n")
  cat("============================================================\n\n")

  print(consistency_table, n = Inf, width = Inf)

  cat("\n图片已保存到文件夹：\n")
  cat(normalizePath(output_dir, winslash = "/", mustWork = FALSE), "\n")

  invisible(list(
    raw_results = result_df,
    metrics_table = metrics_table,
    consistency_table = consistency_table,
    best_runs = best_runs,
    station_info = station_info
  ))
}


# ============================================================
# 11. 运行
# ============================================================

mc_results <- run_monte_carlo_double_cp_case1_fixed(
  N_sim = 500,
  n = n,
  cps = cps_true,
  L = L,
  M = M,
  block_size = block_size,
  boot_method = boot_method,
  lambda_frac = lambda_frac,
  alpha = alpha,
  loc_tol = loc_tol,
  tol_best = tol_best,
  grid_size_scale = grid_size_scale,
  merge_close = merge_close,
  merge_gap = merge_gap,
  target_station = target_station,
  output_dir = output_dir
)
