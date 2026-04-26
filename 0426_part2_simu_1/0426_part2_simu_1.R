library(mclust)
library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork) 

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

# ==========================================
# 数据生成逻辑（针对5种场景重建）
# ==========================================
generate_spatial_data_cases <- function(n=1000, p=2, K, cp_locs, shift_size=2.5, condition) {
  Sigma_eps <- outer(1:p, 1:p, function(i, j) 0.45 ^ abs(i - j))
  A <- if (p == 2) matrix(c(0.50, 0.08, 0.05, 0.45), p, p, byrow = TRUE) else diag(0.5, p)
  
  data_list <- vector("list", K)
  for (i in seq_len(K)) {
    X <- matrix(0, nrow = n, ncol = p)
    for (t in 2:n) X[t, ] <- as.numeric(A %*% X[t - 1, ]) + as.numeric(rmvnorm(1, mean = rep(0, p), sigma = Sigma_eps))
    data_list[[i]] <- X
  }
  
  t1 <- cp_locs[1]
  t2 <- if(length(cp_locs) > 1) cp_locs[2] else cp_locs[1]
  
  for (i in seq_len(K)) {
    if (condition == "2个变量在相同时刻发生异常变化") {
      if (i == 1) {
        data_list[[i]][t1:n, 1] <- data_list[[i]][t1:n, 1] + shift_size
        data_list[[i]][t1:n, 2] <- data_list[[i]][t1:n, 2] + shift_size
      }
    } else if (condition == "2个变量在不同时刻发生异常变化") {
      if (i == 1) {
        data_list[[i]][t1:n, 1] <- data_list[[i]][t1:n, 1] + shift_size
        data_list[[i]][t2:n, 2] <- data_list[[i]][t2:n, 2] + shift_size
      }
    } else if (condition == "在相同时刻，1个变量异常变化，1个变量正常变化") {
      data_list[[i]][t1:n, 1] <- data_list[[i]][t1:n, 1] + shift_size
      if (i == 1) data_list[[i]][t1:n, 2] <- data_list[[i]][t1:n, 2] + shift_size
    } else if (condition == "在不同时刻，1个变量异常变化，1个变量正常变化") {
      data_list[[i]][t1:n, 1] <- data_list[[i]][t1:n, 1] + shift_size
      if (i == 1) data_list[[i]][t2:n, 2] <- data_list[[i]][t2:n, 2] + shift_size
    } else if (condition == "2个变量在相同时刻发生正常变化") {
      data_list[[i]][t1:n, 1] <- data_list[[i]][t1:n, 1] + shift_size
      data_list[[i]][t1:n, 2] <- data_list[[i]][t1:n, 2] + shift_size
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

# 评价函数（兼容多变点评估）
evaluate_metrics <- function(n, cp_true, detected_cps, tol = 50) {
  is_anomaly_cond <- length(cp_true) > 0 
  
  get_labels <- function(n, cps) {
    lbls <- rep(1, n)
    if(length(cps) == 0) return(lbls)
    cps <- sort(cps[cps > 0 & cps < n])
    for(i in seq_along(cps)) lbls[(cps[i]+1):n] <- i + 1
    return(lbls)
  }
  
  true_labels <- get_labels(n, cp_true)
  est_labels <- get_labels(n, detected_cps)
  ARI <- adjustedRandIndex(true_labels, est_labels)
  
  if (is_anomaly_cond) {
    dist_mat <- abs(outer(detected_cps, cp_true, "-"))
    if (length(detected_cps) > 0) {
      min_dists <- apply(dist_mat, 2, min)
      TP <- sum(min_dists <= tol)
      FP <- length(detected_cps) - TP
      MAE <- mean(min_dists) 
      H_dist <- max(max(min_dists), max(apply(dist_mat, 1, min)))
    } else {
      TP <- 0; FP <- 0; MAE <- NA; H_dist <- n
    }
    Precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
    Recall <- TP / length(cp_true)
    F1 <- ifelse(Precision + Recall == 0, 0, 2 * Precision * Recall / (Precision + Recall))
    return(list(TP=TP, FP=FP, MAE=MAE, Hausdorff=H_dist, ARI=ARI, F1=F1, Precision=Precision, Recall=Recall))
  } else {
    FP <- length(detected_cps)
    return(list(TP=0, FP=FP, MAE=NA, Hausdorff=NA, ARI=ARI, F1=NA, Precision=0, Recall=0))
  }
}

run_case_logic <- function(condition, coords, n=1000, p=2, cp_locs=c(400, 700), L=50, M=1000, block_size=15) {
  t1 <- cp_locs[1]; t2 <- cp_locs[2]
  var1_cp_anomaly <- c(); var1_cp_normal <- c()
  var2_cp_anomaly <- c(); var2_cp_normal <- c()
  
  if (condition == "2个变量在相同时刻发生异常变化") {
    var1_cp_anomaly <- c(t1); var2_cp_anomaly <- c(t1)
  } else if (condition == "2个变量在不同时刻发生异常变化") {
    var1_cp_anomaly <- c(t1); var2_cp_anomaly <- c(t2)
  } else if (condition == "在相同时刻，1个变量异常变化，1个变量正常变化") {
    var1_cp_normal <- c(t1); var2_cp_anomaly <- c(t1)
  } else if (condition == "在不同时刻，1个变量异常变化，1个变量正常变化") {
    var1_cp_normal <- c(t1); var2_cp_anomaly <- c(t2)
  } else if (condition == "2个变量在相同时刻发生正常变化") {
    var1_cp_normal <- c(t1); var2_cp_normal <- c(t1)
  }
  
  cp_info <- list(
    all = unique(c(var1_cp_anomaly, var1_cp_normal, var2_cp_anomaly, var2_cp_normal)),
    anomaly = unique(c(var1_cp_anomaly, var2_cp_anomaly)),
    v1_a = var1_cp_anomaly, v1_n = var1_cp_normal,
    v2_a = var2_cp_anomaly, v2_n = var2_cp_normal
  )
  
  data_list <- generate_spatial_data_cases(n, p, nrow(coords), cp_locs, 2.5, condition)
  sn_res <- SNSeg_Multi(data_list[[1]], paras_to_test = "mean", confidence = 0.95, grid_size_scale = 0.05, plot_SN = FALSE)
  
  # 修改：不再丢弃检测到的正常变点，将它们分别存放
  final_anomaly_cps <- c()
  final_normal_cps <- c()
  
  for (ecp in sn_res$est_cp) {
    stage2 <- test_anomaly_sbb(data_list, tau = ecp, L = L, M = M, block_size = block_size, coords = coords)
    if (!is.null(stage2)) {
      if (stage2$p_value < 0.05) {
        if (!any(abs(final_anomaly_cps - ecp) <= 30)) final_anomaly_cps <- c(final_anomaly_cps, ecp)
      } else {
        if (!any(abs(final_normal_cps - ecp) <= 30)) final_normal_cps <- c(final_normal_cps, ecp)
      }
    }
  }
  
  # 评价指标计算只针对被判定为异常的点（红线点）
  metrics <- evaluate_metrics(n, cp_info$anomaly, final_anomaly_cps, tol = 50)
  
  if (length(cp_info$anomaly) > 0) {
    score <- if (metrics$TP > 0) (1000 - metrics$MAE - metrics$FP * 200) else (-5000 - metrics$FP * 200)
  } else {
    score <- -metrics$FP * 100
  }
  
  return(list(
    data_list=data_list, 
    cps=list(anomaly=final_anomaly_cps, normal=final_normal_cps), # 携带所有检测状态
    metrics=metrics, score=score, cp_info=cp_info
  ))
}


generate_best_plots <- function(best_run, condition, cp_info, L, M, block_size, coords, prefix, n, p) {
  data_list <- best_run$data_list
  # 提取红绿分类
  anomaly_cps <- best_run$cps$anomaly
  normal_cps <- best_run$cps$normal
  all_detected <- c(anomaly_cps, normal_cps)
  
  target <- data_list[[1]]
  ref <- weighted_neighbors_idw(data_list, coords, 1, 2)
  
  df_time <- data.frame(
    Time = 1:n,
    Var1_Target = target[, 1], Var1_Ref = ref[, 1],
    Var2_Target = target[, 2], Var2_Ref = ref[, 2]
  )
  
  # === 核心修改点：这里加入了中文翻译，确保图例和数据完美匹配 ===
  df_long <- df_time %>%
    pivot_longer(cols = -Time, names_to = c("Variable", "Type"), names_sep = "_", values_to = "Value") %>%
    mutate(
      Var_zh = ifelse(Variable == "Var1", "变量1", "变量2"),
      Type_zh = ifelse(Type == "Target", "目标站点", "临近站点加权值"),
      Group = paste(Var_zh, "-", Type_zh)
    )
  
  # === 你的图例和线型字典（全部改为中文） ===
  custom_colors <- c(
    "变量1 - 目标站点" = "#e31a1c", "变量1 - 临近站点加权值" = "#06ef58", 
    "变量2 - 目标站点" = "#1aa0c1", "变量2 - 临近站点加权值" = "#420df0"
  )
  custom_linetypes <- c("变量1 - 目标站点" = "solid", "变量1 - 临近站点加权值" = "dashed", "变量2 - 目标站点" = "solid", "变量2 - 临近站点加权值" = "dashed")
  custom_linewidths <- c("变量1 - 目标站点" = 1.0, "变量1 - 临近站点加权值" = 0.7, "变量2 - 目标站点" = 1.0, "变量2 - 临近站点加权值" = 0.7)
  
  plot_ts_base <- function(data, title) {
    ggplot(data, aes(x = Time, y = Value, color = Group, linetype = Group)) +
      geom_line(aes(linewidth = Group), alpha = 0.85) +
      scale_color_manual(values = custom_colors) + scale_linetype_manual(values = custom_linetypes) +
      scale_linewidth_manual(values = custom_linewidths) +
      labs(title = title, x = NULL, y = "Value") + theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        legend.position = "right", legend.title = element_blank(),
        legend.key.width = unit(2.5, "cm"), 
        panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "gray90")
      )
  }
  
  # 检测点分配与颜色判定逻辑修改：
  # 1. 颜色由算法决定：p<0.05 红色，p>=0.05 绿色
  # 2. 如果纯粹属于另一个变量的变点，则当前变量图不画该实线
  # 修改后的颜色判定逻辑：强制根据变量自身的真实物理设定来画色
  build_detected_lines <- function(detected_list, v_a, v_n, v_other_a, v_other_n, tol=50) {
    res <- data.frame(x=numeric(), color=character())
    if(length(detected_list) == 0) return(res)
    
    for (d in detected_list) {
      # 系统整体判定（如果是纯误报时使用）
      sys_pred_col <- ifelse(d %in% anomaly_cps, "red", "green") 
      
      dist_a <- if(length(v_a)>0) min(abs(d - v_a)) else Inf
      dist_n <- if(length(v_n)>0) min(abs(d - v_n)) else Inf
      dist_other <- if(length(c(v_other_a, v_other_n))>0) min(abs(d - c(v_other_a, v_other_n))) else Inf
      
      if (dist_a <= tol) {
        # 1. 变点属于本变量，且真实身份是【异常】，画红色
        res <- rbind(res, data.frame(x=d, color="red"))
      } else if (dist_n <= tol) {
        # 2. 变点属于本变量，且真实身份是【正常】，强制画绿色！
        res <- rbind(res, data.frame(x=d, color="green"))
      } else if (dist_other > tol) {
        # 3. 既不属于本变量，也不属于其他变量（纯粹系统误报），用算法判定的颜色
        res <- rbind(res, data.frame(x=d, color=sys_pred_col))
      }
      # 隐式逻辑：如果是属于另一个变量的变点，当前变量什么线都不画
    }
    return(res)
  }
  
  # 调用时拆分传入，精确定位身份
  lines_v1 <- build_detected_lines(all_detected, cp_info$v1_a, cp_info$v1_n, cp_info$v2_a, cp_info$v2_n)
  lines_v2 <- build_detected_lines(all_detected, cp_info$v2_a, cp_info$v2_n, cp_info$v1_a, cp_info$v1_n)
  
  # 联合图表，展示所有检测结果
  lines_combined <- data.frame(x=numeric(), color=character())
  if(length(all_detected) > 0) {
    for(d in all_detected) {
      pred_col <- ifelse(d %in% anomaly_cps, "red", "green")
      lines_combined <- rbind(lines_combined, data.frame(x=d, color=pred_col))
    }
  }
  
  # === 变量1 绘图（已汉化） ===
  p1 <- plot_ts_base(df_long %>% filter(Variable == "Var1"), sprintf("变量1"))
  if(length(cp_info$all) > 0) p1 <- p1 + geom_vline(xintercept = cp_info$all, color = "gray40", linetype = "dotted", linewidth = 1.2)
  if(nrow(lines_v1) > 0) {
    for(i in 1:nrow(lines_v1)) p1 <- p1 + geom_vline(xintercept = lines_v1$x[i], color = lines_v1$color[i], linetype = "solid", linewidth = 1.2)
  }
  
  # === 变量2 绘图（已汉化） ===
  p2 <- plot_ts_base(df_long %>% filter(Variable == "Var2"), sprintf("变量2"))
  if(length(cp_info$all) > 0) p2 <- p2 + geom_vline(xintercept = cp_info$all, color = "gray40", linetype = "dotted", linewidth = 1.2)
  if(nrow(lines_v2) > 0) {
    for(i in 1:nrow(lines_v2)) p2 <- p2 + geom_vline(xintercept = lines_v2$x[i], color = lines_v2$color[i], linetype = "solid", linewidth = 1.2)
  }
  
  # === 联合绘图（已汉化） ===
  p3 <- plot_ts_base(df_long, sprintf("变量1和变量2")) + labs(x = "时间") #+ theme(legend.position = "bottom")
  if(length(cp_info$all) > 0) p3 <- p3 + geom_vline(xintercept = cp_info$all, color = "gray40", linetype = "dotted", linewidth = 1.2)
  if(nrow(lines_combined) > 0) {
    for(i in 1:nrow(lines_combined)) p3 <- p3 + geom_vline(xintercept = lines_combined$x[i], color = lines_combined$color[i], linetype = "solid", linewidth = 1.2)
  }

  combined_plot <- p1 / p2 / p3
  
  # === 新增：为时序图添加总标题 ===
  combined_plot <- combined_plot + plot_annotation(title = sprintf(condition), theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
  
  ts_filename <- sprintf("%s_%s_TimeSeries.png", prefix, condition)
  ggsave(ts_filename, plot = combined_plot, width = 12, height = 10, dpi = 150, bg = "white")
  
  # === SBB 统计量图表（已汉化图例和标题） ===
  tau_grid <- seq(max(L + 5, 30), min(n - L - 5, n - 30), by = 10)
  scan_df <- scan_anomaly_statistics(data_list, tau_grid, L, max(200, floor(M/2)), block_size, coords)
  stat_filename <- sprintf("%s_%s_Statistics.png", prefix, condition)
  
  if(!is.null(scan_df)) {
    p_stat <- ggplot(scan_df, aes(x = tau)) +
      geom_line(aes(y = S_obs, color = "观测统计量"), linewidth = 1.2) +
      geom_line(aes(y = S_boot_q95, color = "95%阈值"), linetype = "dashed", linewidth = 1) +
      scale_color_manual(values = c("观测统计量" = "#3182bd", "95%阈值" = "#de2d26")) +
      labs(title = "SBB 统计量趋势", x = NULL, y = "S_obs") +
      theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.title = element_blank())
    
    if(length(cp_info$all) > 0) p_stat <- p_stat + geom_vline(xintercept = cp_info$all, color = "gray40", linetype = "dotted", linewidth = 1.2)
    if(nrow(lines_combined) > 0) {
      for(i in 1:nrow(lines_combined)) p_stat <- p_stat + geom_vline(xintercept = lines_combined$x[i], color = lines_combined$color[i], linetype = "solid", linewidth = 1.2)
    }
    
    p_pval <- ggplot(scan_df, aes(x = tau)) +
      geom_line(aes(y = p_value), color = "#756bb1", linewidth = 1.2) +
      geom_hline(yintercept = 0.05, color = "#ff7f00", linetype = "dashed", linewidth = 1) +
      labs(title = "P值趋势", x = "时间 (tau)", y = "P-value") + coord_cartesian(ylim = c(0, 1)) +
      theme_minimal(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    if(length(cp_info$all) > 0) p_pval <- p_pval + geom_vline(xintercept = cp_info$all, color = "gray40", linetype = "dotted", linewidth = 1.2)
    if(nrow(lines_combined) > 0) {
      for(i in 1:nrow(lines_combined)) p_pval <- p_pval + geom_vline(xintercept = lines_combined$x[i], color = lines_combined$color[i], linetype = "solid", linewidth = 1.2)
    }
    
    stat_plot <- p_stat / p_pval
    
    # === 新增：为统计图添加总标题 ===
    stat_plot <- stat_plot + plot_annotation(title = sprintf("场景分析: %s - SBB统计检验", condition), theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
    
    ggsave(stat_filename, plot = stat_plot, width = 12, height = 8, dpi = 150, bg = "white")
  }
  
  cat(sprintf("\n=> 最佳运行图表已分离保存:\n   [%s]\n   [%s]\n", ts_filename, stat_filename))
}


run_monte_carlo_cases <- function(N_sim = 50, raw_coords) {
  conditions <- c("2个变量在相同时刻发生异常变化", "2个变量在不同时刻发生异常变化", 
                  "在相同时刻，1个变量异常变化，1个变量正常变化", "在不同时刻，1个变量异常变化，1个变量正常变化", "2个变量在相同时刻发生正常变化")
  coords <- find_and_reorder_stations(raw_coords)
  summary_df <- data.frame()
  
  for (cond in conditions) {
    cat(sprintf("\n=== 正在运行: AR(1) Mean | 场景: %s ===\n", cond))
    
    best_score <- -Inf
    best_run <- NULL
    all_F1 <- c(); all_MAE <- c(); all_Hausdorff <- c(); all_ARI <- c()
    correct_count <- 0; total_fp <- 0
    
    pb <- txtProgressBar(min = 0, max = N_sim, style = 3)
    for (iter in 1:N_sim) {
      res <- run_case_logic(cond, coords, n=1000, p=2, cp_locs=c(400, 700), L=50, M=1000, block_size=15)
      
      is_anomaly_cond <- length(res$cp_info$anomaly) > 0
      
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
    
    generate_best_plots(best_run, cond, best_run$cp_info, L=50, M=1000, block_size=15, coords=coords, prefix="BestRun_AR1", n=1000, p=2)
    
    summary_df <- rbind(summary_df, data.frame(
      Condition = cond,
      Accuracy = round(correct_count / N_sim, 3),
      Avg_FP = round(total_fp / N_sim, 3),
      F1_Score = ifelse(is_anomaly_cond, round(mean(all_F1, na.rm=TRUE), 3), NA),
      Avg_MAE = ifelse(is_anomaly_cond, round(mean(all_MAE, na.rm=TRUE), 2), NA),
      Hausdorff = ifelse(is_anomaly_cond, round(mean(all_Hausdorff, na.rm=TRUE), 2), NA),
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
mc_results <- run_monte_carlo_cases(N_sim = 10, raw_coords = raw_station_coords)