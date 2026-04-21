library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)

# ==========================================
# 0. 动态编译修复后的 C++ DWB 核心逻辑
# (使用全局中心化，完美修复假阳性Bug)
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

static inline vec gen_dwb_weights(int n, double bandwidth) {
  double bw = std::max(1.0, bandwidth);
  double a = std::exp(-1.0 / bw);
  double s = std::sqrt(std::max(1e-12, 1.0 - a * a));
  vec xi(n, fill::zeros);
  xi(0) = R::rnorm(0.0, 1.0);
  for (int t = 1; t < n; ++t) {
    xi(t) = a * xi(t - 1) + s * R::rnorm(0.0, 1.0);
  }
  xi -= mean(xi);
  double sd_xi = stddev(xi);
  if (sd_xi > 1e-12) xi /= sd_xi;
  return xi;
}

// [[Rcpp::export]]
List dwb_test_cpp(mat Z, int L, int M = 500, double bandwidth = 10.0) {
  int n = Z.n_rows;
  int d = Z.n_cols;

  if (n != 2 * L) stop(\"Z.n_rows must be exactly 2 * L.\");
  if (L < 2 || d < 1 || M < 10) stop(\"Invalid inputs.\");

  // 【核心修复1】：在Null Hypothesis下，必须使用全局中心化（Global Centering）！
  // 这样才能保证 Bootstrap 样本方差不会被低估，彻底消除假阳性。
  mat Z_centered = Z;
  rowvec g_mean = mean(Z, 0);
  for (int t = 0; t < n; ++t) {
    Z_centered.row(t) -= g_mean;
  }

  mat Z_norm = Z_centered;
  for (int j = 0; j < d; ++j) {
    double col_sd = stddev(Z_centered.col(j));
    if (col_sd > 1e-8) {
      Z_norm.col(j) /= col_sd;
    }
  }

  // 观测统计量
  double S_obs = calc_stat(Z_norm, L);

  // 构造 Bootstrap 残差
  mat Z_tilde = Z_norm;
  int exceed_count = 0;
  vec S_boot(M, fill::zeros);

  for (int m = 0; m < M; ++m) {
    vec xi = gen_dwb_weights(n, bandwidth);
    mat Z_star = Z_tilde;
    for (int t = 0; t < n; ++t) {
      Z_star.row(t) *= xi(t);
    }
    double S_star = calc_stat(Z_star, L);
    S_boot(m) = S_star;
    if (S_star >= S_obs) {
      exceed_count++;
    }
  }

  double p_value = static_cast<double>(exceed_count + 1) / static_cast<double>(M + 1);

  vec S_sorted = sort(S_boot);
  int idx95 = std::min(M - 1, std::max(0, static_cast<int>(std::floor(0.95 * M)) - 1));
  double S_q95 = S_sorted(idx95);

  return List::create(
    _[\"p_value\"] = p_value,
    _[\"S_obs\"] = S_obs,
    _[\"S_boot_q95\"] = S_q95
  );
}
"
# 编译执行 C++
Rcpp::sourceCpp(code = cpp_code)

# ==========================================
# 1. 坐标处理与距离加权 (IDW) 模块
# ==========================================
vech <- function(mat) mat[lower.tri(mat, diag = TRUE)]

haversine_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180; R <- 6371.0
  dlon <- (lon2 - lon1) * rad; dlat <- (lat2 - lat1) * rad
  a <- sin(dlat/2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  return(R * c)
}

find_and_reorder_stations <- function(coords) {
  K <- nrow(coords)
  dist_mat <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      dist_mat[i, j] <- haversine_dist(coords[i,1], coords[i,2], coords[j,1], coords[j,2])
    }
  }
  sum_dists <- rowSums(dist_mat)
  center_idx <- which.min(sum_dists)
  
  cat(sprintf("=> 选定第 %d 号站点为中心目标站点 (距离总和: %.2f km)\n", center_idx, sum_dists[center_idx]))
  
  new_order <- c(center_idx, setdiff(1:K, center_idx))
  return(coords[new_order, ])
}

# 仅用于绘图直观展示
weighted_neighbors_idw <- function(data_list, coords, target_idx = 1, power = 2) {
  K <- length(data_list)
  lon_target <- coords[target_idx, 1]; lat_target <- coords[target_idx, 2]
  ref_indices <- setdiff(1:K, target_idx)
  
  distances <- sapply(ref_indices, function(i) {
    haversine_dist(lon_target, lat_target, coords[i, 1], coords[i, 2])
  })
  distances[distances < 1e-6] <- 1e-6 
  norm_weights <- (1 / (distances ^ power)) / sum(1 / (distances ^ power))
  
  ref_matrix <- matrix(0, nrow = nrow(data_list[[1]]), ncol = ncol(data_list[[1]]))
  for (j in seq_along(ref_indices)) {
    ref_matrix <- ref_matrix + data_list[[ref_indices[j]]] * norm_weights[j]
  }
  return(ref_matrix)
}

# ==========================================
# 2. 4个变点的数据生成模块
# ==========================================
apply_covariance_shift_segment <- function(X, start_idx, end_idx, cov_scale) {
  n <- nrow(X); p <- ncol(X)
  if (start_idx > end_idx || start_idx > n) return(X)
  end_idx <- min(end_idx, n)
  
  if (start_idx > 1) {
    idx_pre <- 1:(start_idx - 1)
    mu_pre <- colMeans(X[idx_pre, , drop = FALSE])
  } else {
    mu_pre <- colMeans(X[start_idx:end_idx, , drop = FALSE])
  }
  
  B <- diag(rep(1, p))
  diag(B) <- seq(cov_scale, max(1.05, cov_scale * 0.7), length.out = p)
  if (p >= 2) B[1, 2] <- 0.25; if (p >= 3) B[2, 3] <- -0.15
  
  X_seg <- X[start_idx:end_idx, , drop = FALSE]
  X_centered <- sweep(X_seg, 2, mu_pre, "-")
  X_seg_new <- X_centered %*% t(B)
  X_seg_new <- sweep(X_seg_new, 2, mu_pre, "+")
  
  mean_gap <- colMeans(X_seg_new) - mu_pre
  X_seg_new <- sweep(X_seg_new, 2, mean_gap, "-")
  
  X[start_idx:end_idx, ] <- X_seg_new
  return(X)
}

generate_spatial_data_4cp <- function(n = 1000, p = 2, K_stations = 6, 
                                      cp_locs = c(400, 800, 1200, 1600), 
                                      shift_sizes = c(2.5, -3.0), 
                                      anomaly_only = TRUE, scenario = "mean") {
  Sigma_eps <- outer(1:p, 1:p, function(i, j) 0.45 ^ abs(i - j))
  if (p == 2) {
    A <- matrix(c(0.50, 0.08, 0.05, 0.45), p, p, byrow = TRUE)
  } else {
    A <- diag(0.5, p)
  }
  
  data_list <- vector("list", K_stations)
  for (i in seq_len(K_stations)) {
    X <- matrix(0, nrow = n, ncol = p)
    for (t in 2:n) {
      X[t, ] <- as.numeric(A %*% X[t - 1, ]) + as.numeric(rmvnorm(1, mean = rep(0, p), sigma = Sigma_eps))
    }
    data_list[[i]] <- X
  }
  
  stations_to_change <- if (anomaly_only) 1 else seq_len(K_stations)
  cp1 <- cp_locs[1]; cp2 <- cp_locs[2]
  cp3 <- cp_locs[3]; cp4 <- cp_locs[4]
  
  for (i in stations_to_change) {
    X <- data_list[[i]]
    if (scenario == "mean") {
      X[cp1:(cp2-1), ] <- sweep(X[cp1:(cp2-1), , drop=FALSE], 2, rep(shift_sizes[1], p), "+")
      X[cp3:(cp4-1), ] <- sweep(X[cp3:(cp4-1), , drop=FALSE], 2, rep(shift_sizes[2], p), "+")
    } else if (scenario == "cov") {
      X <- apply_covariance_shift_segment(X, cp1, cp2-1, cov_scale = shift_sizes[1])
      X <- apply_covariance_shift_segment(X, cp3, cp4-1, cov_scale = shift_sizes[2])
    } else if (scenario == "multiparam") {
      X[cp1:(cp2-1), 1] <- X[cp1:(cp2-1), 1] * shift_sizes[1] + shift_sizes[1]
      X[cp3:(cp4-1), 1] <- X[cp3:(cp4-1), 1] * shift_sizes[2] + shift_sizes[2]
    }
    data_list[[i]] <- X
  }
  return(data_list)
}

# ==========================================
# 3. 特征提取与变点检验逻辑（【核心修复2】：完美非线性抵消法）
# ==========================================
build_feature_matrix <- function(data_list, tau, L, scenario = "mean", coords = NULL) {
  K <- length(data_list)
  p <- ncol(data_list[[1]])
  n <- nrow(data_list[[1]])
  start_idx <- max(1, tau - L + 1)
  end_idx <- min(n, tau + L)
  actual_L <- tau - start_idx + 1
  if (actual_L < 20 || (end_idx - start_idx + 1) != 2 * actual_L) return(NULL)
  
  window_list <- lapply(data_list, function(mat) mat[start_idx:end_idx, , drop = FALSE])
  len <- end_idx - start_idx + 1
  
  lon_target <- coords[1, 1]; lat_target <- coords[1, 2]
  ref_indices <- 2:K
  distances <- sapply(ref_indices, function(i) haversine_dist(lon_target, lat_target, coords[i, 1], coords[i, 2]))
  distances[distances < 1e-6] <- 1e-6 
  norm_weights <- (1 / (distances ^ 2)) / sum(1 / (distances ^ 2))
  
  Z <- NULL
  
  # 数学绝杀：直接提取平方特征后再做加权减法，完美抵消任何全局参数的跳跃期望！
  if (scenario == "mean") {
    ref_data <- matrix(0, nrow = len, ncol = p)
    for (j in seq_along(ref_indices)) {
      ref_data <- ref_data + window_list[[ref_indices[j]]] * norm_weights[j]
    }
    Z <- window_list[[1]] - ref_data
    
  } else if (scenario == "cov") {
    d_vech <- p * (p + 1) / 2
    Z <- matrix(0, nrow = len, ncol = d_vech)
    
    for (t in seq_len(len)) {
      v_t <- window_list[[1]][t, ]
      v_target <- vech(v_t %*% t(v_t))  # 目标站点的协方差特征
      
      v_ref <- rep(0, d_vech)
      for (j in seq_along(ref_indices)) {
        v_r <- window_list[[ref_indices[j]]][t, ]
        v_ref <- v_ref + vech(v_r %*% t(v_r)) * norm_weights[j]  # 参考站点特征的加权和
      }
      Z[t, ] <- v_target - v_ref
    }
    
  } else if (scenario == "multiparam") {
    Z <- matrix(0, nrow = len, ncol = 2)
    
    for (t in seq_len(len)) {
      val_t <- window_list[[1]][t, 1]
      v_target <- c(val_t, val_t^2)  # 目标站点的一阶、二阶特征
      
      v_ref <- c(0, 0)
      for (j in seq_along(ref_indices)) {
        val_r <- window_list[[ref_indices[j]]][t, 1]
        v_ref <- v_ref + c(val_r, val_r^2) * norm_weights[j]  # 参考站点特征的加权和
      }
      Z[t, ] <- v_target - v_ref
    }
  }
  
  return(list(Z = Z, actual_L = actual_L, start_idx = start_idx, end_idx = end_idx))
}

test_anomaly_dwb <- function(data_list, tau, L, M, bandwidth = 10, scenario = "mean", coords = NULL) {
  feat <- build_feature_matrix(data_list, tau, L, scenario, coords)
  if (is.null(feat)) return(NULL)
  out <- dwb_test_cpp(feat$Z, feat$actual_L, M, bandwidth)
  return(c(out, list(actual_L = feat$actual_L, tau = tau)))
}

detect_cp_stage1 <- function(target_mat, scenario) {
  if (scenario == "mean") {
    SNSeg_Multi(target_mat, paras_to_test = "mean", confidence = 0.99, grid_size_scale = 0.05, plot_SN = FALSE)
  } else if (scenario == "cov") {
    SNSeg_Multi(target_mat, paras_to_test = "covariance", confidence = 0.95, grid_size_scale = 0.05, plot_SN = FALSE)
  } else if (scenario == "multiparam") {
    SNSeg_Uni(target_mat[, 1], paras_to_test = c("mean", "variance"), confidence = 0.99, grid_size_scale = 0.05, plot_SN = FALSE)
  }
}

scan_anomaly_statistics <- function(data_list, tau_grid, L, M, bandwidth, scenario, coords) {
  rows <- lapply(tau_grid, function(tau) {
    res <- test_anomaly_dwb(data_list, tau, L, M, bandwidth, scenario, coords)
    if (!is.null(res)) data.frame(tau=tau, S_obs=res$S_obs, S_boot_q95=res$S_boot_q95, p_value=res$p_value) else NULL
  })
  do.call(rbind, rows)
}

# ==========================================
# 4. 单次运行与绘图函数
# ==========================================
run_case_4cp <- function(scenario = "mean", condition = "Anomaly", coords, n = 2000, p = 2, 
                         cp_true = c(400, 800, 1200, 1600), shift_sizes = c(2.5, -2.5), 
                         alpha = 0.05, L = 50, M = 1200, bandwidth = 10,  # 调小bandwidth适应局部窗口
                         save_plot = FALSE, plot_prefix = "plot") {
  
  if (scenario == "cov") shift_sizes <- c(1.6, 2.2)
  
  anomaly_only <- (condition == "Anomaly")
  data_list <- generate_spatial_data_4cp(n, p, nrow(coords), cp_true, shift_sizes, anomaly_only, scenario)

  # 第1阶段：预检
  sn_res <- detect_cp_stage1(data_list[[1]], scenario)
  est_cps_raw <- sn_res$est_cp
  
  final_detected_cps <- c()
  p_values_list <- c()
  
  # 第2阶段：严格使用重构后的 DWB 空间检验
  if (length(est_cps_raw) > 0) {
    for (ecp in est_cps_raw) {
      stage2 <- test_anomaly_dwb(data_list, tau = ecp, L = L, M = M, bandwidth = bandwidth, scenario = scenario, coords = coords)
      if (!is.null(stage2)) {
        p_val <- as.numeric(stage2$p_value)
        if (p_val < alpha) {
          if (!any(abs(final_detected_cps - ecp) <= 30)) {
            final_detected_cps <- c(final_detected_cps, ecp)
            p_values_list <- c(p_values_list, p_val)
          }
        }
      }
    }
  }
  
  # 准确率统计
  cp1_found <- FALSE; cp2_found <- FALSE; cp3_found <- FALSE; cp4_found <- FALSE
  if (length(final_detected_cps) > 0) {
    cp1_found <- any(abs(final_detected_cps - cp_true[1]) <= 50)
    cp2_found <- any(abs(final_detected_cps - cp_true[2]) <= 50)
    cp3_found <- any(abs(final_detected_cps - cp_true[3]) <= 50)
    cp4_found <- any(abs(final_detected_cps - cp_true[4]) <= 50)
  }
  
  is_correct <- FALSE
  if (condition == "Anomaly") {
    if (cp1_found && cp2_found && cp3_found && cp4_found) is_correct <- TRUE
  } else {
    if (length(final_detected_cps) == 0) is_correct <- TRUE
  }
  
  # ----------------- 画图模块 -----------------
  if (save_plot) {
    tau_grid <- seq(max(L + 5, 30), min(n - L - 5, n - 30), by = 10)
    scan_df <- scan_anomaly_statistics(data_list, tau_grid, L=L, M=max(200, floor(M/2)), bandwidth=bandwidth, scenario=scenario, coords=coords)
    
    pdf_filename <- sprintf("%s_%s_%s_4CP.pdf", plot_prefix, scenario, condition)
    pdf(pdf_filename, width = 12, height = 8)
    
    target <- data_list[[1]]
    ref <- weighted_neighbors_idw(data_list, coords, 1, 2)
    
    plot_p <- ifelse(scenario == "multiparam", 1, p)
    
    par(mfrow=c(plot_p, 1))
    for(j in 1:plot_p) {
      y <- c(target[,j], ref[,j])
      plot(1:n, target[,j], type="l", col="#d7301f", main=sprintf("Var %d (%s-%s)", j, scenario, condition), ylab="Val", ylim=range(y))
      lines(1:n, ref[,j], col="#2171b5", lty=2)
      abline(v=cp_true, col="gray40", lty=3, lwd=2)
      if(length(final_detected_cps) > 0) abline(v=final_detected_cps, col="#2ca25f", lty=1, lwd=2)
    }
    
    par(mfrow=c(2,1))
    if (!is.null(scan_df) && nrow(scan_df) > 0) {
      plot(scan_df$tau, scan_df$S_obs, type="l", col="blue", main="Statistic Trend (4 CPs)", ylab="S_obs", lwd=1.5)
      lines(scan_df$tau, scan_df$S_boot_q95, col="red", lty=2)
      abline(v=cp_true, col="gray40", lty=3, lwd=2)
      if(length(final_detected_cps) > 0) abline(v=final_detected_cps, col="#2ca25f", lty=1, lwd=2)
      
      plot(scan_df$tau, scan_df$p_value, type="l", col="purple", main="P-value Trend", ylab="P-value", ylim=c(0,1), lwd=1.5)
      abline(h=alpha, col="orange", lty=2, lwd=2)
      abline(v=cp_true, col="gray40", lty=3, lwd=2)
      if(length(final_detected_cps) > 0) abline(v=final_detected_cps, col="#2ca25f", lty=1, lwd=2)
    }
    dev.off()
    cat(sprintf("=> 已保存过程图至: %s\n", pdf_filename))
  }
  
  return(list(
    scenario = scenario, condition = condition,
    cp1_found = cp1_found, cp2_found = cp2_found,
    cp3_found = cp3_found, cp4_found = cp4_found,
    false_positives = sum(!sapply(final_detected_cps, function(x) any(abs(x - cp_true) <= 50))),
    is_correct = is_correct
  ))
}

# ==========================================
# 5. 蒙特卡洛主循环器
# ==========================================
run_monte_carlo_4cp <- function(N_sim = 50, raw_coords) {
  scenarios <- c("mean", "cov", "multiparam")
  conditions <- c("Anomaly", "Normal")
  
  cat("================= 准备阶段 =================\n")
  coords <- find_and_reorder_stations(raw_coords)
  cat("============================================\n\n")
  
  results_summary <- data.frame()
  
  for (scen in scenarios) {
    for (cond in conditions) {
      cat(sprintf("Running 4-CP Simulation: [%s] | Condition:[%s] ... ", toupper(scen), toupper(cond)))
      
      correct_count <- 0
      cp1_rec <- 0; cp2_rec <- 0; cp3_rec <- 0; cp4_rec <- 0
      total_fp <- 0
      
      pb <- txtProgressBar(min = 0, max = N_sim, style = 3)
      for (iter in 1:N_sim) {
        do_save_plot <- (iter == 1) 
        
        res <- run_case_4cp(scenario = scen, condition = cond, coords = coords, 
                            save_plot = do_save_plot, plot_prefix = "MC_4CP")
        
        if (res$is_correct) correct_count <- correct_count + 1
        if (res$cp1_found) cp1_rec <- cp1_rec + 1
        if (res$cp2_found) cp2_rec <- cp2_rec + 1
        if (res$cp3_found) cp3_rec <- cp3_rec + 1
        if (res$cp4_found) cp4_rec <- cp4_rec + 1
        total_fp <- total_fp + res$false_positives
        
        setTxtProgressBar(pb, iter)
      }
      close(pb)
      
      acc <- correct_count / N_sim
      r1 <- cp1_rec / N_sim; r2 <- cp2_rec / N_sim
      r3 <- cp3_rec / N_sim; r4 <- cp4_rec / N_sim
      fpr <- total_fp / N_sim 
      
      results_summary <- rbind(results_summary, data.frame(
        Scenario = scen, Condition = cond, Total_Sim = N_sim,
        Strict_Accuracy = acc, CP1_Recall = r1, CP2_Recall = r2, 
        CP3_Recall = r3, CP4_Recall = r4, Avg_False_Positives = fpr
      ))
    }
  }
  
  cat("\n=================== FINAL REPORT (4 CPs) ===================\n")
  print(results_summary)
  cat("============================================================\n")
  
  return(results_summary)
}

# ==========================================
# 6. 执行示例代码
# ==========================================
raw_station_coords <- matrix(c(
  116.41,39.92,
  116.397,39.982,
  116.339,39.929,
  116.461,39.937,
  116.407,39.866,
  116.352,39.878 
), ncol = 2, byrow = TRUE)

set.seed(42)

# 执行验证 (测试运行 1 次)
mc_results <- run_monte_carlo_4cp(N_sim = 1, raw_coords = raw_station_coords)