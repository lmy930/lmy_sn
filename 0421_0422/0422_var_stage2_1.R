library(mclust)
library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)

# ==========================================
# 1. C++ 模块 (无需修改：通过对数映射后，继续使用均值跳跃检验)
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
# 3. 数据生成模块（修改点：方差突变，保持均值和相关性）
# ==========================================
generate_spatial_data_ar1_1cp_var <- function(n=1000, p=2, K, cp_loc=500, var_scale=2.5, condition="Anomaly_Both") {
  # 基础噪声协方差矩阵和AR1转移矩阵
  Sigma_eps <- outer(1:p, 1:p, function(i, j) 0.45 ^ abs(i - j))
  A <- if (p == 2) matrix(c(0.50, 0.08, 0.05, 0.45), p, p, byrow = TRUE) else diag(0.5, p)
  
  data_list <- vector("list", K)
  for (i in seq_len(K)) {
    X <- matrix(0, nrow = n, ncol = p)
    for (t in 2:n) {
      # 每一步生成独立的白噪声
      noise <- as.numeric(rmvnorm(1, mean = rep(0, p), sigma = Sigma_eps))
      scale_vec <- c(1, 1) # 默认不放缩
      
      if (t >= cp_loc) {
        # Var1 的方差放缩逻辑
        if ((condition == "Anomaly_Both" && i == 1) ||
            (condition == "Normal") ||
            (condition == "Anomaly_Var1" && i == 1) ||
            (condition == "Anomaly_Var2")) {
          scale_vec[1] <- var_scale
        }
        
        # Var2 的方差放缩逻辑
        if ((condition == "Anomaly_Both" && i == 1) ||
            (condition == "Normal") ||
            (condition == "Anomaly_Var1") ||
            (condition == "Anomaly_Var2" && i == 1)) {
          scale_vec[2] <- var_scale
        }
      }
      # 将放缩比例乘在当前步的噪声上，累加进自回归过程
      X[t, ] <- as.numeric(A %*% X[t - 1, ]) + noise * scale_vec
    }
    data_list[[i]] <- X
  }
  return(data_list)
}

# ==========================================
# 4. 检测算法与特征工程模块（修改点：对数方差映射）
# ==========================================
#对数方差映射
# build_feature_matrix_var <- function(data_list, tau, L, coords) {
#   n <- nrow(data_list[[1]]); p <- ncol(data_list[[1]])
#   s_idx <- max(1, tau - L + 1); e_idx <- min(n, tau + L)
#   actual_L <- tau - s_idx + 1
#   if (actual_L < 20 || (e_idx - s_idx + 1) != 2 * actual_L) return(NULL)
  
#   win_list <- lapply(data_list, function(mat) mat[s_idx:e_idx, , drop = FALSE])
#   target_mat <- win_list[[1]]
#   ref_mat <- weighted_neighbors_idw(win_list, coords, 1, 2)
  
#   # a. 局部强制中心化（只测波动，剥离局部小漂移）
#   target_centered <- scale(target_mat, center = TRUE, scale = FALSE)
#   ref_centered <- scale(ref_mat, center = TRUE, scale = FALSE)
  
#   # b. 对数平方映射（将方差突变转换为均值突变）
#   epsilon <- 1e-6
#   Y_target <- log(target_centered^2 + epsilon)
#   Y_ref <- log(ref_centered^2 + epsilon)
  
#   # c. 空间差值：在对数域相减，完美抵消气候背景同等倍数的方差扩大
#   Z_mapped <- Y_target - Y_ref
  
#   return(list(Z = Z_mapped, actual_L = actual_L))
# }

# 绝对值方差映射
build_feature_matrix_var <- function(data_list, tau, L, coords) {
  n <- nrow(data_list[[1]])
  s_idx <- max(1, tau - L + 1); e_idx <- min(n, tau + L)
  actual_L <- tau - s_idx + 1
  if (actual_L < 20 || (e_idx - s_idx + 1) != 2 * actual_L) return(NULL)
  
  win_list <- lapply(data_list, function(mat) mat[s_idx:e_idx, , drop = FALSE])
  target_mat <- win_list[[1]]
  ref_mat <- weighted_neighbors_idw(win_list, coords, 1, 2)
  
  # 1. 局部去均值 (只留波动)
  target_centered <- scale(target_mat, center = TRUE, scale = FALSE)
  ref_centered <- scale(ref_mat, center = TRUE, scale = FALSE)
  
  # 2. 核心修改：使用绝对值映射 (稳健、无极值)
  Y_target <- abs(target_centered)
  Y_ref <- abs(ref_centered)
  
  # 3. 空间差分
  Z <- Y_target - Y_ref
  
  return(list(Z = Z, actual_L = actual_L))
}


test_anomaly_sbb_var <- function(data_list, tau, L, M, block_size, coords) {
  feat <- build_feature_matrix_var(data_list, tau, L, coords)
  if (is.null(feat)) return(NULL)
  out <- sbb_test_cpp(feat$Z, feat$actual_L, M, block_size)
  return(c(out, list(tau = tau)))
}

scan_anomaly_statistics_var <- function(data_list, tau_grid, L, M, block_size, coords) {
  rows <- lapply(tau_grid, function(tau) {
    res <- test_anomaly_sbb_var(data_list, tau, L, M, block_size, coords)
    if (!is.null(res)) data.frame(tau=tau, S_obs=res$S_obs, S_boot_q95=res$S_boot_q95, p_value=res$p_value) else NULL
  })
  do.call(rbind, rows)
}

evaluate_metrics <- function(n, cp_true, detected_cps, tol = 50, condition) {
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

run_case_logic_var <- function(condition, coords, n=1000, p=2, cp_true=500, L=50, M=1000, block_size=15) {
  data_list <- generate_spatial_data_ar1_1cp_var(n, p, nrow(coords), cp_true, 2.5, condition)
  
  # 修改点：Stage 1 的粗筛必须使用 "covariance" 才能捕获方差/协方差的突变
  sn_res <- SNSeg_Multi(data_list[[1]], paras_to_test = "covariance", confidence = 0.95, grid_size_scale = 0.05, plot_SN = FALSE)
  
  final_detected_cps <- c()
  for (ecp in sn_res$est_cp) {
    # Stage 2 调用更新后的对数映射检验函数
    stage2 <- test_anomaly_sbb_var(data_list, tau = ecp, L = L, M = M, block_size = block_size, coords = coords)
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


# ==========================================
# 5. 可视化与模拟模块 (完整保留)
# ==========================================
# generate_best_plots_var <- function(best_run, condition, cp_true, L, M, block_size, coords, prefix, n, p) {
#   data_list <- best_run$data_list
#   detected_cps <- best_run$cps
  
#   # 生成扫描统计量绘图数据，调用方差版扫描函数
#   tau_grid <- seq(max(L + 5, 30), min(n - L - 5, n - 30), by = 10)
#   scan_df <- scan_anomaly_statistics_var(data_list, tau_grid, L, max(200, floor(M/2)), block_size, coords)
  
#   target <- data_list[[1]]; ref <- weighted_neighbors_idw(data_list, coords, 1, 2)
  
#   png(filename = sprintf("%s_%s_Variance_TimeSeries.png", prefix, condition), width = 1200, height = 1000, res = 120)
#   par(mfrow=c(p + 1, 1), mar=c(4, 4, 2, 1)) 
  
#   for(j in 1:p) {
#     # 绘制原始时间序列，方差突变表现为包络线（波动）变宽
#     plot(1:n, target[,j], type="l", col="#d7301f", main=sprintf("Var %d Time Series (Variance | %s)", j, condition), ylab="Value", ylim=range(c(target[,j], ref[,j])))
#     lines(1:n, ref[,j], col="#2171b5", lty=2)
    
#     if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
    
#     is_anomaly_var <- FALSE
#     if (condition == "Anomaly_Both") is_anomaly_var <- TRUE
#     else if (condition == "Anomaly_Var1" && j == 1) is_anomaly_var <- TRUE
#     else if (condition == "Anomaly_Var2" && j == 2) is_anomaly_var <- TRUE
    
#     cp_col <- if(is_anomaly_var) "red" else "green"
#     if(length(detected_cps) > 0) abline(v=detected_cps, col=cp_col, lty=1, lwd=2)
#   }
  
#   plot(1:n, target[,1], type="l", col="#e41a1c", main=sprintf("Combined View: Both Variables (%s)", condition), ylab="Value", ylim=range(c(target, ref)))
#   lines(1:n, ref[,1], col="#e41a1c", lty=2, lwd=1.5)
#   lines(1:n, target[,2], col="#377eb8", lty=1)
#   lines(1:n, ref[,2], col="#377eb8", lty=2, lwd=1.5)
  
#   if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
#   if(length(detected_cps) > 0) abline(v=detected_cps, col="black", lty=1, lwd=2)
  
#   legend("topleft", legend=c("Var1 Target", "Var1 Ref", "Var2 Target", "Var2 Ref"),
#          col=c("#e41a1c", "#e41a1c", "#377eb8", "#377eb8"),
#          lty=c(1, 2, 1, 2), cex=0.8, horiz=TRUE, bty="n", lwd=c(1,1.5,1,1.5))
  
#   dev.off()
  
  
#   png(filename = sprintf("%s_%s_Variance_Statistics.png", prefix, condition), width = 1200, height = 800, res = 120)
#   par(mfrow=c(2, 1), mar=c(4, 4, 2, 1))
#   if(!is.null(scan_df)) {
#     plot(scan_df$tau, scan_df$S_obs, type="l", col="blue", main="Log-Variance SBB Statistic Trend", ylab="S_obs", lwd=1.5)
#     lines(scan_df$tau, scan_df$S_boot_q95, col="red", lty=2)
#     if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
#     if(length(detected_cps) > 0) abline(v=detected_cps, col="#2ca25f", lty=1, lwd=2)
    
#     plot(scan_df$tau, scan_df$p_value, type="l", col="purple", main="P-value Trend", ylab="P-value", ylim=c(0,1), lwd=1.5)
#     abline(h=0.05, col="orange", lty=2, lwd=2)
#     if(condition != "Normal") abline(v=cp_true, col="gray40", lty=3, lwd=2)
#     if(length(detected_cps) > 0) abline(v=detected_cps, col="#2ca25f", lty=1, lwd=2)
#   }
#   dev.off()
#   cat(sprintf("\n=> 最佳运行图表已分离保存: [%s] 和 [%s]\n", 
#               sprintf("%s_%s_Variance_TimeSeries.png", prefix, condition), 
#               sprintf("%s_%s_Variance_Statistics.png", prefix, condition)))
# }


generate_best_plots_var <- function(best_run, condition, cp_true, L, M, block_size, coords, prefix, n, p) {
  require(ggplot2)
  require(tidyr)
  require(dplyr)
  require(patchwork) 
  
  data_list <- best_run$data_list
  detected_cps <- best_run$cps
  
  target <- data_list[[1]]
  ref <- weighted_neighbors_idw(data_list, coords, 1, 2)
  
  # ==========================================
  # 1. 数据转换
  # ==========================================
  df_time <- data.frame(
    Time = 1:n,
    Var1_Target = target[, 1],
    Var1_Ref    = ref[, 1],
    Var2_Target = target[, 2],
    Var2_Ref    = ref[, 2]
  )
  
  df_long <- df_time %>%
    pivot_longer(cols = -Time, names_to = c("Variable", "Type"), names_sep = "_", values_to = "Value") %>%
    mutate(Group = paste(Variable, Type, sep = " - "))
  
  # ==========================================
  # 2. 颜色与线型强对比映射
  # ==========================================
  # 彻底区分目标站和临近站的颜色
  custom_colors <- c(
    "Var1 - Target" = "#e31a1c", # 鲜红色 (Var1 目标)
    "Var1 - Ref"    = "#06ef58", # 亮橙色 (Var1 临近)
    "Var2 - Target" = "#1aa0c1", # 深蓝色 (Var2 目标)
    "Var2 - Ref"    = "#420df0"  # 浅绿色 (Var2 临近)
  )
  custom_linetypes <- c(
    "Var1 - Target" = "solid", "Var1 - Ref" = "dashed",
    "Var2 - Target" = "solid", "Var2 - Ref" = "dashed"
  )
  custom_linewidths <- c(
    "Var1 - Target" = 1.0, "Var1 - Ref" = 0.7,
    "Var2 - Target" = 1.0, "Var2 - Ref" = 0.7
  )
  
  # ==========================================
  # 3. 封装时序图绘图函数 (新增 cp_color 参数控制竖线颜色)
  # ==========================================
  plot_ts <- function(data, title, cp_color) {
    p_base <- ggplot(data, aes(x = Time, y = Value, color = Group, linetype = Group)) +
      geom_line(aes(linewidth = Group), alpha = 0.85) +
      scale_color_manual(values = custom_colors) +
      scale_linetype_manual(values = custom_linetypes) +
      scale_linewidth_manual(values = custom_linewidths) +
      # 真实突变点 (灰色虚线)
      {if(condition != "Normal") geom_vline(xintercept = cp_true, color = "gray40", linetype = "dotted", linewidth = 1.2)} +
      # 检测到的突变点 (根据传入的颜色：红/绿)
      {if(length(detected_cps) > 0) geom_vline(xintercept = detected_cps, color = cp_color, linetype = "solid", linewidth = 1.2)} +
      labs(title = title, x = NULL, y = "Value") +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.width = unit(2.5, "cm"), # 拉长图例，方便看清虚线
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "gray90")
      )
    return(p_base)
  }
  
  # ==========================================
  # 4. 判定竖线颜色与生成图表
  # ==========================================
  # Var 1 的突变判定：如果场景是 Both 或 Var1 异常，则标红，否则标绿
  is_anomaly_var1 <- condition %in% c("Anomaly_Both", "Anomaly_Var1")
  cp_col1 <- ifelse(is_anomaly_var1, "red", "green")
  
  # Var 2 的突变判定
  is_anomaly_var2 <- condition %in% c("Anomaly_Both", "Anomaly_Var2")
  cp_col2 <- ifelse(is_anomaly_var2, "red", "green")
  
  # 组合图和统计图突变判定：只要不是全正常(Normal)，统一用红色标示异常
  cp_col_combined <- ifelse(condition != "Normal", "red", "green")
  
  # 生成三张子图
  p1 <- plot_ts(df_long %>% filter(Variable == "Var1"), sprintf("Var 1 Time Series (Variance | %s)", condition), cp_col1)
  p2 <- plot_ts(df_long %>% filter(Variable == "Var2"), sprintf("Var 2 Time Series (Variance | %s)", condition), cp_col2)
  p3 <- plot_ts(df_long, sprintf("Combined View: Both Variables (%s)", condition), cp_col_combined) +
        labs(x = "Time") +
        theme(legend.position = "bottom")
  
  # 使用 patchwork 拼图
  combined_plot <- p1 / p2 / p3
  
  ts_filename <- sprintf("%s_%s_Variance_TimeSeries.png", prefix, condition)
  ggsave(ts_filename, plot = combined_plot, width = 12, height = 10, dpi = 150, bg = "white")
  
  # ==========================================
  # 5. 生成统计量图 (同样适配红色/绿色竖线)
  # ==========================================
  tau_grid <- seq(max(L + 5, 30), min(n - L - 5, n - 30), by = 10)
  scan_df <- scan_anomaly_statistics_var(data_list, tau_grid, L, max(200, floor(M/2)), block_size, coords)
  
  stat_filename <- sprintf("%s_%s_Variance_Statistics.png", prefix, condition)
  
  if(!is.null(scan_df)) {
    p_stat <- ggplot(scan_df, aes(x = tau)) +
      geom_line(aes(y = S_obs, color = "S_obs"), linewidth = 1.2) +
      geom_line(aes(y = S_boot_q95, color = "95% Threshold"), linetype = "dashed", linewidth = 1) +
      scale_color_manual(values = c("S_obs" = "#3182bd", "95% Threshold" = "#de2d26")) +
      {if(condition != "Normal") geom_vline(xintercept = cp_true, color = "gray40", linetype = "dotted", linewidth = 1.2)} +
      {if(length(detected_cps) > 0) geom_vline(xintercept = detected_cps, color = cp_col_combined, linetype = "solid", linewidth = 1.2)} +
      labs(title = "Log-Variance SBB Statistic Trend", x = NULL, y = "S_obs") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.title = element_blank(), legend.position = "right")
    
    p_pval <- ggplot(scan_df, aes(x = tau)) +
      geom_line(aes(y = p_value), color = "#756bb1", linewidth = 1.2) +
      geom_hline(yintercept = 0.05, color = "#ff7f00", linetype = "dashed", linewidth = 1) +
      {if(condition != "Normal") geom_vline(xintercept = cp_true, color = "gray40", linetype = "dotted", linewidth = 1.2)} +
      {if(length(detected_cps) > 0) geom_vline(xintercept = detected_cps, color = cp_col_combined, linetype = "solid", linewidth = 1.2)} +
      labs(title = "P-value Trend", x = "Time (tau)", y = "P-value") +
      coord_cartesian(ylim = c(0, 1)) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    stat_plot <- p_stat / p_pval
    ggsave(stat_filename, plot = stat_plot, width = 12, height = 8, dpi = 150, bg = "white")
  }
  
  cat(sprintf("\n=> 最佳运行图表已生成:\n   [%s]\n   [%s]\n", ts_filename, stat_filename))
}

run_monte_carlo_1cp_var <- function(N_sim = 50, raw_coords) {
  conditions <- c("Anomaly_Both", "Anomaly_Var1", "Anomaly_Var2", "Normal")
  coords <- find_and_reorder_stations(raw_coords)
  summary_df <- data.frame()
  
  for (cond in conditions) {
    cat(sprintf("\n=== 正在运行: AR(1) VARIANCE | 场景: %s ===\n", toupper(cond)))
    
    best_score <- -Inf
    best_run <- NULL
    
    all_F1 <- c(); all_MAE <- c(); all_Hausdorff <- c(); all_ARI <- c()
    correct_count <- 0; total_fp <- 0
    
    pb <- txtProgressBar(min = 0, max = N_sim, style = 3)
    for (iter in 1:N_sim) {
      res <- run_case_logic_var(cond, coords, n=1000, p=2, cp_true=500, L=50, M=1000, block_size=15)
      
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
    
    generate_best_plots_var(best_run, cond, cp_true=500, L=50, M=1000, block_size=15, coords=coords, prefix="BestRun_AR1", n=1000, p=2)
    
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

# 运行入口
raw_station_coords <- matrix(c(
  116.41,39.92, 116.397,39.982, 116.339,39.929,
  116.461,39.937, 116.407,39.866, 116.352,39.878 
), ncol = 2, byrow = TRUE)

# N_sim设置为2用于快速测试，正式跑建议设置为50-100
set.seed(2026)
mc_results_var <- run_monte_carlo_1cp_var(N_sim = 5, raw_coords = raw_station_coords)