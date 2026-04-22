library(mclust)
library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(zoo)

# ==========================================
# 1. Rcpp 加速的 SBB 统计量模块
# ==========================================
cpp_code <- "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

static inline double calc_stat(const mat& Z, int L) {
  rowvec mean_before = mean(Z.rows(0, L - 1), 0);
  rowvec mean_after = mean(Z.rows(L, 2 * L - 1), 0);
  return accu(square(mean_after - mean_before));
}

// [[Rcpp::export]]
List sbb_test_cpp(mat Z, int L, int M = 500, double expected_block_size = 5.0) {
  int n = Z.n_rows;
  int d = Z.n_cols;

  if (n != 2 * L) stop(\"Z.n_rows must be exactly 2 * L.\");
  if (L < 2 || d < 1 || M < 10) stop(\"Invalid inputs.\");

  mat Z_norm = Z;
  rowvec g_mean = mean(Z, 0);
  for (int t = 0; t < n; ++t) Z_norm.row(t) -= g_mean;
  
  for (int j = 0; j < d; ++j) {
    double col_sd = stddev(Z_norm.col(j));
    if (col_sd > 1e-8) Z_norm.col(j) /= col_sd;
  }
  
  double S_obs = calc_stat(Z_norm, L);

  mat Z_tilde = Z_norm;
  rowvec m_left = mean(Z_norm.rows(0, L - 1), 0);
  rowvec m_right = mean(Z_norm.rows(L, 2 * L - 1), 0);
  
  for (int t = 0; t < L; ++t) Z_tilde.row(t) -= m_left;
  for (int t = L; t < 2 * L; ++t) Z_tilde.row(t) -= m_right;

  int exceed_count = 0;
  vec S_boot(M, fill::zeros);
  double p_geom = 1.0 / std::max(2.0, expected_block_size);

  for (int m = 0; m < M; ++m) {
    mat Z_star(n, d, fill::zeros);
    int t = 0;
    while (t < n) {
      int block_length = R::rgeom(p_geom) + 1; 
      int start_idx = std::floor(R::runif(0.0, n)); 
      for (int i = 0; i < block_length && t < n; ++i) {
        int idx = (start_idx + i) % n; 
        Z_star.row(t) = Z_tilde.row(idx);
        t++;
      }
    }
    double S_star = calc_stat(Z_star, L);
    S_boot(m) = S_star;
    if (S_star >= S_obs) exceed_count++;
  }

  double p_value = static_cast<double>(exceed_count + 1) / static_cast<double>(M + 1);
  vec S_sorted = sort(S_boot);
  int idx95 = std::min(M - 1, std::max(0, static_cast<int>(std::floor(0.95 * M)) - 1));
  
  return List::create(_[\"p_value\"] = p_value, _[\"S_obs\"] = S_obs, _[\"S_boot_q95\"] = S_sorted(idx95));
}
"
Rcpp::sourceCpp(code = cpp_code)

# ==========================================
# 2. 基础工具模块
# ==========================================
haversine_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180; R <- 6371.0
  dlon <- (lon2 - lon1) * rad; dlat <- (lat2 - lat1) * rad
  a <- sin(dlat/2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon/2)^2
  return(R * 2 * atan2(sqrt(a), sqrt(1 - a)))
}

find_and_reorder_stations <- function(coords) {
  K <- nrow(coords)
  dist_mat <- matrix(0, K, K)
  for (i in 1:K) for (j in 1:K) dist_mat[i, j] <- haversine_dist(coords[i,1], coords[i,2], coords[j,1], coords[j,2])
  center_idx <- which.min(rowSums(dist_mat))
  cat(sprintf("=> 选定中心站点 ID: %d\n", center_idx))
  return(coords[c(center_idx, setdiff(1:K, center_idx)), ])
}

weighted_neighbors_idw <- function(data_list, coords, target_idx = 1, power = 2) {
  K <- length(data_list)
  lon_t <- coords[target_idx, 1]; lat_t <- coords[target_idx, 2]
  ref_idx <- setdiff(1:K, target_idx)
  
  dists <- sapply(ref_idx, function(i) haversine_dist(lon_t, lat_t, coords[i, 1], coords[i, 2]))
  dists[dists < 1e-6] <- 1e-6 
  w <- (1 / (dists ^ power)) / sum(1 / (dists ^ power))
  
  ref_mat <- matrix(0, nrow = nrow(data_list[[1]]), ncol = ncol(data_list[[1]]))
  for (j in seq_along(ref_idx)) ref_mat <- ref_mat + data_list[[ref_idx[j]]] * w[j]
  return(ref_mat)
}

# ==========================================
# 3. 实际数据加载与异常注入模块
# ==========================================
load_real_data <- function(dir_path) {
  station_files <- c(
    "Beijing_Aotizhongxin.csv", "Beijing_Dongsi.csv", "Beijing_Guanyuan.csv", 
    "Beijing_Nongzhanguan.csv", "Beijing_Tiantan.csv", "Beijing_Wanshouxigong.csv"
  )
  data_list <- list()
  for (i in seq_along(station_files)) {
    file_path <- file.path(dir_path, station_files[i])
    df <- read.csv(file_path)
    mat <- as.matrix(df[, c("PM2.5", "PM10")])
    # 填补实际数据中的缺失值
    mat[,1] <- na.approx(mat[,1], rule=2)
    mat[,2] <- na.approx(mat[,2], rule=2)
    data_list[[i]] <- mat
  }
  return(data_list)
}

generate_spatial_data_real_1cp <- function(base_data_list, cp_loc, shift_size, condition) {
  data_list <- base_data_list
  K <- length(data_list)
  n <- nrow(data_list[[1]])
  
  for (i in seq_len(K)) {
    if ((condition == "Anomaly_Both" && i == 1) || (condition == "Normal") ||
        (condition == "Anomaly_Var1" && i == 1) || (condition == "Anomaly_Var2")) {
      data_list[[i]][cp_loc:n, 1] <- data_list[[i]][cp_loc:n, 1] + shift_size
    }
    
    if ((condition == "Anomaly_Both" && i == 1) || (condition == "Normal") ||
        (condition == "Anomaly_Var1") || (condition == "Anomaly_Var2" && i == 1)) {
      data_list[[i]][cp_loc:n, 2] <- data_list[[i]][cp_loc:n, 2] + shift_size
    }
  }
  return(data_list)
}

# ==========================================
# 4. 分析与严格单变点逻辑
# ==========================================
build_feature_matrix_mean <- function(data_list, tau, L, coords) {
  n <- nrow(data_list[[1]]); p <- ncol(data_list[[1]])
  s_idx <- max(1, tau - L + 1); e_idx <- min(n, tau + L)
  actual_L <- tau - s_idx + 1
  if (actual_L < 4 || (e_idx - s_idx + 1) != 2 * actual_L) return(NULL) 
  
  win_list <- lapply(data_list, function(mat) mat[s_idx:e_idx, , drop = FALSE])
  ref_mat <- weighted_neighbors_idw(win_list, coords, 1, 2)
  return(list(Z = win_list[[1]] - ref_mat, actual_L = actual_L))
}

test_anomaly_sbb <- function(data_list, tau, L, M, block_size, coords) {
  feat <- build_feature_matrix_mean(data_list, tau, L, coords)
  if (is.null(feat)) return(NULL)
  out <- sbb_test_cpp(feat$Z, feat$actual_L, M, block_size)
  return(c(out, list(tau = tau)))
}

scan_anomaly_statistics <- function(data_list, tau_grid, L, M, block_size, coords) {
  rows <- lapply(tau_grid, function(tau) {
    res <- test_anomaly_sbb(data_list, tau, L, M, block_size, coords)
    if (!is.null(res)) data.frame(tau=tau, S_obs=res$S_obs, S_boot_q95=res$S_boot_q95, p_value=res$p_value) else NULL
  })
  do.call(rbind, rows)
}

evaluate_metrics <- function(n, cp_true, detected_cps, tol = 12, condition) {
  is_anomaly_cond <- condition != "Normal"
  
  get_labels <- function(n, cps) {
    lbls <- rep(1, n)
    if(length(cps) == 0) return(lbls)
    cps <- sort(cps[cps > 0 & cps < n])
    for(i in seq_along(cps)) lbls[(cps[i]+1):n] <- i + 1
    return(lbls)
  }
  
  true_labels <- get_labels(n, if(is_anomaly_cond) cp_true else c())
  est_labels <- get_labels(n, detected_cps)
  ARI <- adjustedRandIndex(true_labels, est_labels)
  
  if (is_anomaly_cond) {
    dist_mat <- abs(outer(detected_cps, cp_true, "-"))
    if (length(detected_cps) > 0) {
      min_dists <- apply(dist_mat, 2, min)
      TP <- sum(min_dists <= tol)
      FP <- max(0, length(detected_cps) - TP)
      MAE <- min(min_dists)
      H_dist <- max(max(min_dists), max(apply(dist_mat, 1, min)))
    } else {
      TP <- 0; FP <- 0; MAE <- NA; H_dist <- n
    }
    Precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
    Recall <- TP / 1
    F1 <- ifelse(Precision + Recall == 0, 0, 2 * Precision * Recall / (Precision + Recall))
    return(list(TP=TP, FP=FP, MAE=MAE, Hausdorff=H_dist, ARI=ARI, F1=F1, Precision=Precision, Recall=Recall))
  } else {
    FP <- length(detected_cps)
    return(list(TP=0, FP=FP, MAE=NA, Hausdorff=NA, ARI=ARI, F1=NA, Precision=0, Recall=0))
  }
}

run_case_logic <- function(condition, base_data_list, coords, cp_true, shift_size, L=12, M=1000, block_size=4) {
  n <- nrow(base_data_list[[1]])
  
  data_list <- generate_spatial_data_real_1cp(base_data_list, cp_true, shift_size, condition)
  sn_res <- SNSeg_Multi(data_list[[1]], paras_to_test = "mean", confidence = 0.95, grid_size_scale = 0.05, plot_SN = FALSE)
  
  # === 核心修改点：严格限定只能挑出 1 个变点 ===
  candidates <- list()
  for (ecp in sn_res$est_cp) {
    stage2 <- test_anomaly_sbb(data_list, tau = ecp, L = L, M = M, block_size = block_size, coords = coords)
    if (!is.null(stage2) && stage2$p_value < 0.05) {
      candidates[[length(candidates) + 1]] <- stage2
    }
  }
  
  final_detected_cps <- c()
  # 如果有多个候选点通过检验，仅选择 S_obs（统计量差值）绝对最大的那一个
  if (length(candidates) > 0) {
    s_obs_vals <- sapply(candidates, function(x) x$S_obs)
    best_idx <- which.max(s_obs_vals)
    final_detected_cps <- candidates[[best_idx]]$tau
  }
  # ============================================

  metrics <- evaluate_metrics(n, cp_true, final_detected_cps, tol = L, condition)
  
  if (condition != "Normal") {
    score <- if (metrics$TP > 0) (1000 - metrics$MAE - metrics$FP * 200) else (-5000 - metrics$FP * 200)
  } else {
    score <- -metrics$FP * 100
  }
  
  return(list(data_list=data_list, cps=final_detected_cps, metrics=metrics, score=score))
}

# ==========================================
# 5. 美化版绘图函数
# ==========================================
generate_best_plots <- function(best_run, condition, cp_true, L, M, block_size, coords, prefix, n, p) {
  data_list <- best_run$data_list
  detected_cps <- best_run$cps
  
  tau_grid <- seq(max(L + 2, 5), min(n - L - 2, n - 5), by = 1)
  scan_df <- scan_anomaly_statistics(data_list, tau_grid, L, max(100, floor(M/2)), block_size, coords)
  
  target <- data_list[[1]]; ref <- weighted_neighbors_idw(data_list, coords, 1, 2)
  var_names <- c("PM2.5", "PM10")
  
  # 图1：时间序列美化
  png(filename = sprintf("%s_%s_TimeSeries.png", prefix, condition), width = 1200, height = 1000, res = 150)
  par(mfrow=c(p + 1, 1), mar=c(4, 4.5, 2.5, 1), cex.main=1.2, cex.lab=1.1) 
  
  for(j in 1:p) {
    # 核心美化点：type="o" (点+线)，pch=20 (实心小圆), lwd=2 (线条加粗)
    plot(1:n, target[,j], type="o", pch=20, cex=0.8, lwd=2, col="#d7301f", 
         main=sprintf("%s Time Series (%s)", var_names[j], condition), 
         ylab=sprintf("%s (ug/m3)", var_names[j]), xlab="Time Index", ylim=range(c(target[,j], ref[,j])))
    grid(col = "gray85", lty = 1, lwd = 1) # 添加浅色网格背景
    lines(1:n, ref[,j], type="o", pch=20, cex=0.8, col="#2171b5", lty=2, lwd=2)
    
    if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=3)
    
    is_anomaly_var <- FALSE
    if (condition == "Anomaly_Both") is_anomaly_var <- TRUE
    else if (condition == "Anomaly_Var1" && j == 1) is_anomaly_var <- TRUE
    else if (condition == "Anomaly_Var2" && j == 2) is_anomaly_var <- TRUE
    
    cp_col <- if(is_anomaly_var) "red" else "green"
    if(length(detected_cps) > 0) abline(v=detected_cps, col=cp_col, lty=1, lwd=3)
  }
  
  # 联合视图美化
  plot(1:n, target[,1], type="o", pch=20, cex=0.7, lwd=2, col="#e41a1c", 
       main=sprintf("Combined View: PM2.5 & PM10 (%s)", condition), 
       ylab="Concentration", xlab="Time Index", ylim=range(c(target, ref)))
  grid(col = "gray85", lty = 1, lwd = 1)
  lines(1:n, ref[,1], type="o", pch=20, cex=0.7, col="#e41a1c", lty=2, lwd=2)
  lines(1:n, target[,2], type="o", pch=18, cex=0.7, col="#377eb8", lty=1, lwd=2)
  lines(1:n, ref[,2], type="o", pch=18, cex=0.7, col="#377eb8", lty=2, lwd=2)
  
  if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=3)
  if(length(detected_cps) > 0) abline(v=detected_cps, col="black", lty=1, lwd=3) 
  
  # 图例加白底框以防遮挡线条
  legend("topright", legend=c("PM2.5 Target", "PM2.5 Ref", "PM10 Target", "PM10 Ref"),
         col=c("#e41a1c", "#e41a1c", "#377eb8", "#377eb8"),
         pch=c(20, 20, 18, 18), lty=c(1, 2, 1, 2), lwd=2, cex=0.9, bg="white", box.col="gray")
  
  dev.off()
  
  # 图2：统计量美化
  png(filename = sprintf("%s_%s_Statistics.png", prefix, condition), width = 1200, height = 800, res = 150)
  par(mfrow=c(2, 1), mar=c(4, 4.5, 2.5, 1), cex.main=1.2)
  if(!is.null(scan_df)) {
    plot(scan_df$tau, scan_df$S_obs, type="o", pch=20, cex=0.6, lwd=2, col="#08519c", 
         main="SBB Statistic Trend", ylab="S_obs", xlab="Time Index")
    grid(col = "gray85", lty = 1)
    lines(scan_df$tau, scan_df$S_boot_q95, type="l", col="#cb181d", lty=2, lwd=2)
    if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=3)
    if(length(detected_cps) > 0) abline(v=detected_cps, col="#2ca25f", lty=1, lwd=3)
    
    plot(scan_df$tau, scan_df$p_value, type="o", pch=20, cex=0.6, lwd=2, col="#54278f", 
         main="P-value Trend", ylab="P-value", xlab="Time Index", ylim=c(0,1))
    grid(col = "gray85", lty = 1)
    abline(h=0.05, col="#f16913", lty=2, lwd=2)
    if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=3)
    if(length(detected_cps) > 0) abline(v=detected_cps, col="#2ca25f", lty=1, lwd=3)
  }
  dev.off()
}

# ==========================================
# 6. 主程序运行
# ==========================================
run_real_data_eval <- function(data_dir, N_sim = 5, raw_coords) {
  base_data_list <- load_real_data(data_dir)
  n_len <- nrow(base_data_list[[1]])
  
  cp_true <- floor(n_len / 2) 
  # 核心修改点：大幅度增加真实数据的注入偏移量，盖过真实数据的本底波动噪音
  shift_size <- 100 
  
  conditions <- c("Anomaly_Both", "Anomaly_Var1", "Anomaly_Var2", "Normal")
  coords <- find_and_reorder_stations(raw_coords)
  summary_df <- data.frame()
  
  for (cond in conditions) {
    cat(sprintf("\n=== 正在运行: 实际数据严格单变点 | 场景: %s ===\n", toupper(cond)))
    
    best_score <- -Inf
    best_run <- NULL
    
    all_F1 <- c(); all_MAE <- c(); all_Hausdorff <- c(); all_ARI <- c()
    correct_count <- 0; total_fp <- 0
    
    pb <- txtProgressBar(min = 0, max = N_sim, style = 3)
    for (iter in 1:N_sim) {
      res <- run_case_logic(cond, base_data_list, coords, cp_true=cp_true, shift_size=shift_size, L=12, M=500, block_size=4)
      
      is_anomaly_cond <- cond != "Normal"
      
      if(is_anomaly_cond) {
        if(res$metrics$TP > 0) correct_count <- correct_count + 1
        all_F1 <- c(all_F1, res$metrics$F1)
        if(!is.na(res$metrics$MAE)) all_MAE <- c(all_MAE, res$metrics$MAE)
        all_Hausdorff <- c(all_Hausdorff, res$metrics$Hausdorff)
      } else {
        if(res$metrics$FP == 0) correct_count <- correct_count + 1
      }
      total_fp <- total_fp + res$metrics$FP
      all_ARI <- c(all_ARI, res$metrics$ARI)
      
      if (res$score > best_score) {
        best_score <- res$score
        best_run <- res
      }
      setTxtProgressBar(pb, iter)
    }
    close(pb)
    
    generate_best_plots(best_run, cond, cp_true=cp_true, L=12, M=500, block_size=4, coords=coords, prefix="BestRun_RealData_1CP", n=n_len, p=2)
    
    summary_df <- rbind(summary_df, data.frame(
      Condition = cond,
      Accuracy = round(correct_count / N_sim, 3),
      Avg_FP = round(total_fp / N_sim, 3),
      F1_Score = ifelse(cond != "Normal", round(mean(all_F1, na.rm=TRUE), 3), NA),
      Avg_MAE = ifelse(cond != "Normal", round(mean(all_MAE, na.rm=TRUE), 2), NA),
      Hausdorff = ifelse(cond != "Normal", round(mean(all_Hausdorff, na.rm=TRUE), 2), NA),
      ARI = round(mean(all_ARI, na.rm=TRUE), 3)
    ))
  }
  
  print(summary_df)
  return(summary_df)
}

raw_station_coords <- matrix(c(
  116.397, 39.982, 116.417, 39.929, 116.339, 39.929,
  116.461, 39.937, 116.407, 39.886, 116.352, 39.878 
), ncol = 2, byrow = TRUE)

set.seed(2026)
mc_results <- run_real_data_eval(data_dir = "2025_new", N_sim = 5, raw_coords = raw_station_coords)