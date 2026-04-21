library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)

# 编译 C++ 文件（确保 0415.cpp 与 R 脚本在同一目录下）
Rcpp::sourceCpp("0415.cpp")

# ==========================================
# 1. 坐标处理与距离加权 (IDW) 模块
# ==========================================
vech <- function(mat) mat[lower.tri(mat, diag = TRUE)]

# 计算两点经纬度之间的球面距离 (单位: km)
haversine_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180; R <- 6371.0
  dlon <- (lon2 - lon1) * rad; dlat <- (lat2 - lat1) * rad
  a <- sin(dlat/2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  return(R * c)
}

# 需求1：寻找空间中心站点，并重排坐标矩阵（目标站点设为第1行）
find_and_reorder_stations <- function(coords) {
  K <- nrow(coords)
  if (K != 6) warning("输入站点数量不是6个，但代码仍会继续寻找中心点。")
  
  dist_mat <- matrix(0, K, K)
  for (i in 1:K) {
    for (j in 1:K) {
      dist_mat[i, j] <- haversine_dist(coords[i,1], coords[i,2], coords[j,1], coords[j,2])
    }
  }
  # 计算每个站点到其他站点的距离总和，最小的即为中心
  sum_dists <- rowSums(dist_mat)
  center_idx <- which.min(sum_dists)
  
  cat(sprintf("=> 选定第 %d 号站点为中心目标站点 (距离总和: %.2f km)\n", center_idx, sum_dists[center_idx]))
  
  # 重排：把中心站点放在第1行，其余顺延
  new_order <- c(center_idx, setdiff(1:K, center_idx))
  return(coords[new_order, ])
}

# 反距离加权提取参考数据
weighted_neighbors_idw <- function(data_list, coords, target_idx = 1, power = 2) {
  K <- length(data_list)
  lon_target <- coords[target_idx, 1]; lat_target <- coords[target_idx, 2]
  ref_indices <- setdiff(1:K, target_idx)
  
  distances <- sapply(ref_indices, function(i) {
    haversine_dist(lon_target, lat_target, coords[i, 1], coords[i, 2])
  })
  distances[distances < 1e-6] <- 1e-6 # 防止除零
  norm_weights <- (1 / (distances ^ power)) / sum(1 / (distances ^ power))
  
  ref_matrix <- matrix(0, nrow = nrow(data_list[[1]]), ncol = ncol(data_list[[1]]))
  for (j in seq_along(ref_indices)) {
    ref_matrix <- ref_matrix + data_list[[ref_indices[j]]] * norm_weights[j]
  }
  return(ref_matrix)
}

# ==========================================
# 2. 多变点时间序列生成模块
# ==========================================
# 局部分段协方差偏移函数
apply_covariance_shift_segment <- function(X, start_idx, end_idx, cov_scale) {
  n <- nrow(X); p <- ncol(X)
  if (start_idx > end_idx || start_idx > n) return(X)
  end_idx <- min(end_idx, n)
  
  idx_pre <- 1:(start_idx - 1)
  mu_pre <- colMeans(X[idx_pre, , drop = FALSE])
  
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

# 生成包含2个变点的空间数据
generate_spatial_data_2cp <- function(n = 600, p = 4, K_stations = 6, 
                                      cp_locs = c(200, 400), 
                                      shift_sizes = c(2, -2), # 两个变点对应的偏移量
                                      anomaly_only = TRUE, scenario = "mean") {
  Sigma_eps <- outer(1:p, 1:p, function(i, j) 0.45 ^ abs(i - j))
  A <- matrix(c(0.50, 0.08, 0.00, 0.00, 0.05, 0.45, 0.07, 0.00, 0.00, 0.06, 0.42, 0.08, 0.00, 0.00, 0.06, 0.40), p, p, byrow = TRUE)
  
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
  
  for (i in stations_to_change) {
    X <- data_list[[i]]
    if (scenario == "mean") {
      X[cp1:(cp2-1), ] <- sweep(X[cp1:(cp2-1), , drop=FALSE], 2, rep(shift_sizes[1], p), "+")
      X[cp2:n, ] <- sweep(X[cp2:n, , drop=FALSE], 2, rep(shift_sizes[2], p), "+")
    } else if (scenario == "cov") {
      # 协方差的 shift_sizes 对应缩放系数，例如 c(1.6, 2.2)
      X <- apply_covariance_shift_segment(X, cp1, cp2-1, cov_scale = shift_sizes[1])
      X <- apply_covariance_shift_segment(X, cp2, n, cov_scale = shift_sizes[2])
    } else if (scenario == "multiparam") {
      X[cp1:(cp2-1), 1] <- X[cp1:(cp2-1), 1] * shift_sizes[1] + shift_sizes[1]
      X[cp2:n, 1] <- X[cp2:n, 1] * shift_sizes[2] + shift_sizes[2]
    }
    data_list[[i]] <- X
  }
  return(data_list)
}

# ==========================================
# 3. 特征提取与多变点检验逻辑
# ==========================================
# (此段与单变点逻辑相同，依然计算局部窗口特征)
build_feature_matrix <- function(data_list, tau, L = 100, scenario = "mean", coords = NULL) {
  K <- length(data_list); p <- ncol(data_list[[1]]); n <- nrow(data_list[[1]])
  start_idx <- max(1, tau - L + 1); end_idx <- min(n, tau + L)
  actual_L <- tau - start_idx + 1
  if (actual_L < 20 || (end_idx - start_idx + 1) != 2 * actual_L) return(NULL)
  
  window_list <- lapply(data_list, function(mat) mat[start_idx:end_idx, , drop = FALSE])
  target_data <- window_list[[1]]
  ref_data <- weighted_neighbors_idw(window_list, coords, 1, 2)
  
  local_center <- function(mat, tau_loc) {
    out <- mat
    out[1:tau_loc, ] <- sweep(out[1:tau_loc, , drop=F], 2, colMeans(mat[1:tau_loc, , drop=F]), "-")
    out[(tau_loc+1):nrow(mat), ] <- sweep(out[(tau_loc+1):nrow(mat), , drop=F], 2, colMeans(mat[(tau_loc+1):nrow(mat), , drop=F]), "-")
    return(out)
  }
  
  Z <- NULL
  if (scenario == "mean") {
    Z <- target_data - ref_data
  } else if (scenario == "cov") {
    d_vech <- p*(p+1)/2
    Z <- matrix(0, nrow = nrow(target_data), ncol = d_vech)
    tc <- local_center(target_data, actual_L); rc <- lapply(window_list[-1], local_center, tau_loc=actual_L)
    for (t in seq_len(nrow(tc))) {
      v_ref <- rep(0, d_vech)
      for (i in seq_along(rc)) v_ref <- v_ref + vech(tcrossprod(rc[[i]][t, ])) / length(rc)
      Z[t, ] <- vech(tcrossprod(tc[t, ])) - v_ref
    }
  } else if (scenario == "multiparam") {
    tc <- local_center(target_data, actual_L); rc <- lapply(window_list[-1], local_center, tau_loc=actual_L)
    r_var <- rep(0, nrow(target_data))
    for (i in seq_along(rc)) r_var <- r_var + (rc[[i]][, 1] ^ 2) / length(rc)
    Z <- cbind(target_data[, 1] - ref_data[, 1], tc[, 1] ^ 2 - r_var)
  }
  return(list(Z = Z, actual_L = actual_L, start_idx = start_idx, end_idx = end_idx))
}

test_anomaly_dwb <- function(data_list, tau, L = 100, M = 500, bandwidth = 10, scenario = "mean", coords = NULL) {
  feat <- build_feature_matrix(data_list, tau, L, scenario, coords)
  if (is.null(feat)) return(NULL)
  out <- dwb_test_cpp(feat$Z, feat$actual_L, M, bandwidth)
  return(c(out, list(actual_L = feat$actual_L, tau = tau)))
}

detect_cp_stage1 <- function(target_mat, scenario) {
  # 注意此处放宽了搜索条件，以捕捉多个变点
  if (scenario == "mean") {
    SNSeg_Multi(target_mat, paras_to_test = "mean", confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
  } else if (scenario == "cov") {
    SNSeg_Multi(target_mat, paras_to_test = "covariance", confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
  } else {
    SNSeg_Multi(target_mat, paras_to_test = "mean", confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
  }
}

scan_anomaly_statistics <- function(data_list, tau_grid, L = 100, M = 200, bandwidth = 10, scenario = "mean", coords = NULL) {
  rows <- lapply(tau_grid, function(tau) {
    res <- test_anomaly_dwb(data_list, tau, L, M, bandwidth, scenario, coords)
    if (!is.null(res)) data.frame(tau=tau, S_obs=res$S_obs, S_boot_q95=res$S_boot_q95, p_value=res$p_value) else NULL
  })
  do.call(rbind, rows)
}

# ==========================================
# 4. 单次运行与绘图函数 (多变点专属逻辑)
# ==========================================
run_case_2cp <- function(scenario = "mean", condition = "Anomaly", coords, n = 600, p = 4, 
                         cp_true = c(200, 400), shift_sizes = c(2.5, -2.5), 
                         alpha = 0.05, L = 100, M = 500, bandwidth = 10,
                         save_plot = FALSE, plot_prefix = "plot") {
  
  if (scenario == "cov") shift_sizes <- c(1.6, 2.2) # 若为协方差，自动替换为缩放系数
  
  anomaly_only <- (condition == "Anomaly")
  data_list <- generate_spatial_data_2cp(n, p, nrow(coords), cp_true, shift_sizes, anomaly_only, scenario)

  # 第一阶段：SNSeg 预检出多个点
  sn_res <- detect_cp_stage1(data_list[[1]], scenario)
  est_cps_raw <- sn_res$est_cp
  
  final_detected_cps <- c()
  p_values_list <- c()
  
  # 第二阶段：遍历每个找出来的嫌疑点进行 DWB 检验
  if (length(est_cps_raw) > 0) {
    for (ecp in est_cps_raw) {
      stage2 <- test_anomaly_dwb(data_list, tau = ecp, L = L, M = M, bandwidth = bandwidth, scenario = scenario, coords = coords)
      if (!is.null(stage2)) {
        p_val <- as.numeric(stage2$p_value)
        if (p_val < alpha) {
          final_detected_cps <- c(final_detected_cps, ecp)
          p_values_list <- c(p_values_list, p_val)
        }
      }
    }
  }
  
  # 多变点正确率判定
  cp1_found <- FALSE; cp2_found <- FALSE
  if (length(final_detected_cps) > 0) {
    cp1_found <- any(abs(final_detected_cps - cp_true[1]) <= 50)
    cp2_found <- any(abs(final_detected_cps - cp_true[2]) <= 50)
  }
  
  is_correct <- FALSE
  if (condition == "Anomaly") {
    # 异常场景下：要求两个真实的异常变点都被检出且显著
    if (cp1_found && cp2_found) is_correct <- TRUE
  } else {
    # 正常场景下：全局漂移，不应有局部异常通过显著性检验，即 final 集合应为空
    if (length(final_detected_cps) == 0) is_correct <- TRUE
  }
  
  # ----------------- 画图模块 -----------------
  if (save_plot) {
    tau_grid <- seq(max(L + 5, 30), min(n - L - 5, n - 30), by = 10)
    scan_df <- scan_anomaly_statistics(data_list, tau_grid, L=L, M=max(200, floor(M/2)), bandwidth=bandwidth, scenario=scenario, coords=coords)
    
    pdf_filename <- sprintf("%s_%s_%s_2CP.pdf", plot_prefix, scenario, condition)
    pdf(pdf_filename, width = 10, height = 8)
    
    target <- data_list[[1]]; ref <- weighted_neighbors_idw(data_list, coords, 1, 2)
    par(mfrow=c(2,2))
    for(j in 1:min(p, 4)) {
      y <- c(target[,j], ref[,j])
      plot(1:n, target[,j], type="l", col="#d7301f", main=sprintf("Var %d (%s-%s)", j, scenario, condition), ylab="Val", ylim=range(y))
      lines(1:n, ref[,j], col="#2171b5", lty=2)
      abline(v=cp_true, col="gray40", lty=3, lwd=2)
      if(length(final_detected_cps) > 0) abline(v=final_detected_cps, col="#2ca25f", lty=1, lwd=2)
    }
    
    par(mfrow=c(2,1))
    if (!is.null(scan_df) && nrow(scan_df) > 0) {
      plot(scan_df$tau, scan_df$S_obs, type="l", col="blue", main="Statistic Trend (Multi-CP)", ylab="S_obs", lwd=1.5)
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
    false_positives = sum(!sapply(final_detected_cps, function(x) any(abs(x - cp_true) <= 50))),
    is_correct = is_correct
  ))
}

# ==========================================
# 5. 蒙特卡洛主循环器 (适配多变点)
# ==========================================
run_monte_carlo_2cp <- function(N_sim = 50, raw_coords) {
  scenarios <- c("mean", "cov", "multiparam")
  conditions <- c("Anomaly", "Normal")
  
  # 在最外层进行一次目标点寻址与坐标重排
  cat("================= 准备阶段 =================\n")
  coords <- find_and_reorder_stations(raw_coords)
  cat("============================================\n\n")
  
  results_summary <- data.frame()
  
  for (scen in scenarios) {
    for (cond in conditions) {
      cat(sprintf("Running Multi-CP Simulation: [%s] | Condition: [%s] ... ", toupper(scen), toupper(cond)))
      
      correct_count <- 0
      cp1_recall <- 0; cp2_recall <- 0
      total_false_positives <- 0
      
      pb <- txtProgressBar(min = 0, max = N_sim, style = 3)
      for (iter in 1:N_sim) {
        do_save_plot <- (iter == 1) # 仅保存第一次的图像
        
        res <- run_case_2cp(scenario = scen, condition = cond, coords = coords, 
                            save_plot = do_save_plot, plot_prefix = "MC_2CP")
        
        if (res$is_correct) correct_count <- correct_count + 1
        if (res$cp1_found) cp1_recall <- cp1_recall + 1
        if (res$cp2_found) cp2_recall <- cp2_recall + 1
        total_false_positives <- total_false_positives + res$false_positives
        
        setTxtProgressBar(pb, iter)
      }
      close(pb)
      
      acc <- correct_count / N_sim
      r1 <- cp1_recall / N_sim; r2 <- cp2_recall / N_sim
      fpr <- total_false_positives / N_sim # 平均每次模拟产生的误报个数
      
      results_summary <- rbind(results_summary, data.frame(
        Scenario = scen, Condition = cond, Total_Sim = N_sim,
        Strict_Accuracy = acc, CP1_Recall = r1, CP2_Recall = r2, Avg_False_Positives = fpr
      ))
    }
  }
  
  cat("\n=================== FINAL REPORT (2 CPs) ===================\n")
  print(results_summary)
  cat("============================================================\n")
  
  return(results_summary)
}

# ==========================================
# 6. 执行示例代码
# ==========================================
# 随机生成 6 个站点的乱序经纬度坐标
raw_station_coords <- matrix(c(
  120.0, 30.0, # 站点 A
  120.1, 30.2, # 站点 B
  119.8, 29.9, # 站点 C
  120.2, 29.8, # 站点 D
  119.9, 30.3, # 站点 E
  120.5, 30.0  # 站点 F
), ncol = 2, byrow = TRUE)

set.seed(2024)

# 执行蒙特卡洛模拟 (为了节约演示时间设为 N=10，实测可调整为 N=100)
# 此函数会自动寻找 raw_station_coords 中最中心的一个作为目标站！
mc_results <- run_monte_carlo_2cp(N_sim = 5, raw_coords = raw_station_coords)