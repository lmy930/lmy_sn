library(mclust)
library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)

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
List sbb_test_cpp(mat Z, int L, int M = 500, double expected_block_size = 15.0) {
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
# 1. 基础工具模块
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

# 修改点：加入 condition 参数以支持 4 种情况
generate_spatial_data_ar1_1cp <- function(n=1000, p=2, K, cp_loc=500, shift_size=2.5, condition="Anomaly_Both") {
  Sigma_eps <- outer(1:p, 1:p, function(i, j) 0.45 ^ abs(i - j))
  A <- if (p == 2) matrix(c(0.50, 0.08, 0.05, 0.45), p, p, byrow = TRUE) else diag(0.5, p)
  
  data_list <- vector("list", K)
  for (i in seq_len(K)) {
    X <- matrix(0, nrow = n, ncol = p)
    for (t in 2:n) X[t, ] <- as.numeric(A %*% X[t - 1, ]) + as.numeric(rmvnorm(1, mean = rep(0, p), sigma = Sigma_eps))
    data_list[[i]] <- X
  }
  
  # 分别对每个变量注入偏移，以此体现“异常”与“正常（一致变化）”的情况
  for (i in seq_len(K)) {
    # 变量 1 的偏移：双变量异常时仅中心站变，正常时全变，Var1异常时仅中心站变，Var2异常时全变(正常一致)
    if ((condition == "Anomaly_Both" && i == 1) ||
        (condition == "Normal") ||
        (condition == "Anomaly_Var1" && i == 1) ||
        (condition == "Anomaly_Var2")) {
      data_list[[i]][cp_loc:n, 1] <- data_list[[i]][cp_loc:n, 1] + shift_size
    }
    
    # 变量 2 的偏移：双变量异常时仅中心站变，正常时全变，Var1异常时全变(正常一致)，Var2异常时仅中心站变
    if ((condition == "Anomaly_Both" && i == 1) ||
        (condition == "Normal") ||
        (condition == "Anomaly_Var1") ||
        (condition == "Anomaly_Var2" && i == 1)) {
      data_list[[i]][cp_loc:n, 2] <- data_list[[i]][cp_loc:n, 2] + shift_size
    }
  }
  return(data_list)
}

build_feature_matrix_mean <- function(data_list, tau, L, coords) {
  n <- nrow(data_list[[1]]); p <- ncol(data_list[[1]])
  s_idx <- max(1, tau - L + 1); e_idx <- min(n, tau + L)
  actual_L <- tau - s_idx + 1
  if (actual_L < 20 || (e_idx - s_idx + 1) != 2 * actual_L) return(NULL)
  
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


evaluate_metrics <- function(n, cp_true, detected_cps, tol = 50, condition) {
  is_anomaly_cond <- condition != "Normal" # 修改点：判断是否含有异常
  
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
      FP <- length(detected_cps) - TP
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

run_case_logic <- function(condition, coords, n=1000, p=2, cp_true=500, L=50, M=1000, block_size=15) {
  data_list <- generate_spatial_data_ar1_1cp(n, p, nrow(coords), cp_true, 2.5, condition)
  sn_res <- SNSeg_Multi(data_list[[1]], paras_to_test = "mean", confidence = 0.95, grid_size_scale = 0.05, plot_SN = FALSE)
  
  final_detected_cps <- c()
  for (ecp in sn_res$est_cp) {
    stage2 <- test_anomaly_sbb(data_list, tau = ecp, L = L, M = M, block_size = block_size, coords = coords)
    if (!is.null(stage2) && stage2$p_value < 0.05) {
      if (!any(abs(final_detected_cps - ecp) <= 30)) final_detected_cps <- c(final_detected_cps, ecp)
    }
  }
  
  metrics <- evaluate_metrics(n, cp_true, final_detected_cps, tol = 50, condition)
  
  if (condition != "Normal") {
    score <- if (metrics$TP > 0) (1000 - metrics$MAE - metrics$FP * 200) else (-5000 - metrics$FP * 200)
  } else {
    score <- -metrics$FP * 100
  }
  
  return(list(data_list=data_list, cps=final_detected_cps, metrics=metrics, score=score))
}


generate_best_plots <- function(best_run, condition, cp_true, L, M, block_size, coords, prefix, n, p) {
  data_list <- best_run$data_list
  detected_cps <- best_run$cps
  
  tau_grid <- seq(max(L + 5, 30), min(n - L - 5, n - 30), by = 10)
  scan_df <- scan_anomaly_statistics(data_list, tau_grid, L, max(200, floor(M/2)), block_size, coords)
  
  target <- data_list[[1]]; ref <- weighted_neighbors_idw(data_list, coords, 1, 2)
  
  # 修改点：扩展高度，并添加第三张合并图
  png(filename = sprintf("%s_%s_TimeSeries.png", prefix, condition), width = 1200, height = 1000, res = 120)
  par(mfrow=c(p + 1, 1), mar=c(4, 4, 2, 1)) # 修改点：增加一行用于联合图
  
  # 图1：逐一绘制时间序列 (包括正常/异常的竖线颜色控制)
  for(j in 1:p) {
    plot(1:n, target[,j], type="l", col="#d7301f", main=sprintf("Var %d Time Series (%s)", j, condition), ylab="Value", ylim=range(c(target[,j], ref[,j])))
    lines(1:n, ref[,j], col="#2171b5", lty=2)
    
    if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
    
    # 动态判断当前变量(j)是否被注入了异常(目标站与临近站不一致)
    is_anomaly_var <- FALSE
    if (condition == "Anomaly_Both") is_anomaly_var <- TRUE
    else if (condition == "Anomaly_Var1" && j == 1) is_anomaly_var <- TRUE
    else if (condition == "Anomaly_Var2" && j == 2) is_anomaly_var <- TRUE
    
    # 核心修改：异常变化变点用红色竖线，正常变化用绿色竖线
    cp_col <- if(is_anomaly_var) "red" else "green"
    
    if(length(detected_cps) > 0) abline(v=detected_cps, col=cp_col, lty=1, lwd=2)
  }
  
  # 图2：新增联合绘图区块
  plot(1:n, target[,1], type="l", col="#e41a1c", main=sprintf("Combined View: Both Variables (%s)", condition), ylab="Value", ylim=range(c(target, ref)))
  lines(1:n, ref[,1], col="#e41a1c", lty=2, lwd=1.5)
  lines(1:n, target[,2], col="#377eb8", lty=1)
  lines(1:n, ref[,2], col="#377eb8", lty=2, lwd=1.5)
  
  if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
  if(length(detected_cps) > 0) abline(v=detected_cps, col="black", lty=1, lwd=2) # 联合图用黑色凸显变点
  
  legend("topleft", legend=c("Var1 Target", "Var1 Ref", "Var2 Target", "Var2 Ref"),
         col=c("#e41a1c", "#e41a1c", "#377eb8", "#377eb8"),
         lty=c(1, 2, 1, 2), cex=0.8, horiz=TRUE, bty="n", lwd=c(1,1.5,1,1.5))
  
  dev.off()
  
  
  png(filename = sprintf("%s_%s_Statistics.png", prefix, condition), width = 1200, height = 800, res = 120)
  par(mfrow=c(2, 1), mar=c(4, 4, 2, 1))
  if(!is.null(scan_df)) {
    plot(scan_df$tau, scan_df$S_obs, type="l", col="blue", main="SBB Statistic Trend", ylab="S_obs", lwd=1.5)
    lines(scan_df$tau, scan_df$S_boot_q95, col="red", lty=2)
    if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
    if(length(detected_cps) > 0) abline(v=detected_cps, col="#2ca25f", lty=1, lwd=2)
    
    plot(scan_df$tau, scan_df$p_value, type="l", col="purple", main="P-value Trend", ylab="P-value", ylim=c(0,1), lwd=1.5)
    abline(h=0.05, col="orange", lty=2, lwd=2)
    if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
    if(length(detected_cps) > 0) abline(v=detected_cps, col="#2ca25f", lty=1, lwd=2)
  }
  dev.off()
  cat(sprintf("\n=> 最佳运行图表已分离保存: [%s] 和 [%s]\n", 
              sprintf("%s_%s_TimeSeries.png", prefix, condition), 
              sprintf("%s_%s_Statistics.png", prefix, condition)))
}

run_monte_carlo_1cp <- function(N_sim = 50, raw_coords) {
  # 修改点：将 condition 数组扩展至 4 种情况
  conditions <- c("Anomaly_Both", "Anomaly_Var1", "Anomaly_Var2", "Normal")
  coords <- find_and_reorder_stations(raw_coords)
  summary_df <- data.frame()
  
  for (cond in conditions) {
    cat(sprintf("\n=== 正在运行: AR(1) Mean | 场景: %s ===\n", toupper(cond)))
    
    best_score <- -Inf
    best_run <- NULL
    
    all_F1 <- c(); all_MAE <- c(); all_Hausdorff <- c(); all_ARI <- c()
    correct_count <- 0; total_fp <- 0
    
    pb <- txtProgressBar(min = 0, max = N_sim, style = 3)
    for (iter in 1:N_sim) {
      res <- run_case_logic(cond, coords, n=1000, p=2, cp_true=500, L=50, M=1000, block_size=15)
      
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
    
    generate_best_plots(best_run, cond, cp_true=500, L=50, M=1000, block_size=15, coords=coords, prefix="BestRun_AR1", n=1000, p=2)
    
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
  116.41,39.92, 116.397,39.982, 116.339,39.929,
  116.461,39.937, 116.407,39.866, 116.352,39.878 
), ncol = 2, byrow = TRUE)

set.seed(2026)
mc_results <- run_monte_carlo_1cp(N_sim = 5, raw_coords = raw_station_coords)