library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(readr)
library(zoo)

# ==========================================
# 0. 动态编译 C++ DWB 核心逻辑 (严格保持原版，逻辑完全一致)
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
  double S_obs = calc_stat(Z_norm, L);
  int exceed_count = 0;
  vec S_boot(M, fill::zeros);
  for (int m = 0; m < M; ++m) {
    vec xi = gen_dwb_weights(n, bandwidth);
    mat Z_star = Z_norm;
    for (int t = 0; t < n; ++t) {
      Z_star.row(t) *= xi(t);
    }
    double S_star = calc_stat(Z_star, L);
    S_boot(m) = S_star;
    if (S_star >= S_obs) exceed_count++;
  }
  double p_value = static_cast<double>(exceed_count + 1) / static_cast<double>(M + 1);
  return List::create(_[\"p_value\"] = p_value, _[\"S_obs\"] = S_obs);
}
"
Rcpp::sourceCpp(code = cpp_code)

# ==========================================
# 1. 基础工具逻辑 (vech & haversine)
# ==========================================
vech <- function(mat) mat[lower.tri(mat, diag = TRUE)]
haversine_dist <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180; R <- 6371.0
  dlon <- (lon2 - lon1) * rad; dlat <- (lat2 - lat1) * rad
  a <- sin(dlat/2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon/2)^2
  return(R * 2 * atan2(sqrt(a), sqrt(1 - a)))
}

# ==========================================
# 2. 注入逻辑 (逻辑与数值模拟代码 apply_covariance_shift_segment 完全一致)
# ==========================================
apply_shift_rigorous <- function(X, start_idx, end_idx, scenario) {
  p <- ncol(X); n <- nrow(X)
  mu_pre <- colMeans(X[1:(start_idx - 1), , drop = FALSE])
  
  if (scenario == "mean") {
    shift_size <- 60
    X[start_idx:end_idx, ] <- sweep(X[start_idx:end_idx, , drop=FALSE], 2, rep(shift_size, p), "+")
    
  } else if (scenario == "cov") {
    cov_scale <- 2.5
    B <- diag(rep(1, p))
    diag(B) <- seq(cov_scale, max(1.05, cov_scale * 0.7), length.out = p)
    if (p >= 2) B[1, 2] <- 0.25 
    
    X_seg <- X[start_idx:end_idx, , drop = FALSE]
    X_centered <- sweep(X_seg, 2, mu_pre, "-")
    X_seg_new <- X_centered %*% t(B)
    X_seg_new <- sweep(X_seg_new, 2, mu_pre, "+")
    
    # 严格执行消除均值跳变逻辑
    mean_gap <- colMeans(X_seg_new) - mu_pre
    X_seg_new <- sweep(X_seg_new, 2, mean_gap, "-")
    X[start_idx:end_idx, ] <- X_seg_new
    
  } else if (scenario == "multiparam") {
    shift_size <- 2.2
    # X * a + a 逻辑
    X[start_idx:end_idx, 1] <- X[start_idx:end_idx, 1] * shift_size + shift_size
  }
  return(X)
}

# ==========================================
# 3. 特征矩阵 Z 构造 (逻辑与数值模拟代码 build_feature_matrix 完全一致)
# ==========================================
build_Z_rigorous <- function(data_list, tau, L, scenario, coords) {
  K <- length(data_list); p <- ncol(data_list[[1]]); n <- nrow(data_list[[1]])
  start_idx <- max(1, tau - L + 1); end_idx <- min(n, tau + L)
  actual_L <- tau - start_idx + 1
  if (actual_L < 10 || (end_idx - start_idx + 1) != 2 * actual_L) return(NULL)
  
  window_list <- lapply(data_list, function(mat) mat[start_idx:end_idx, , drop = FALSE])
  len <- 2 * actual_L
  
  # IDW 权重计算
  ref_indices <- 2:K
  dists <- sapply(ref_indices, function(i) haversine_dist(coords[1,1], coords[1,2], coords[i,1], coords[i,2]))
  norm_weights <- (1 / (dists ^ 2)) / sum(1 / (dists ^ 2))
  
  if (scenario == "mean") {
    ref_data <- matrix(0, nrow = len, ncol = p)
    for (j in seq_along(ref_indices)) ref_data <- ref_data + window_list[[ref_indices[j]]] * norm_weights[j]
    Z <- window_list[[1]] - ref_data
    
  } else if (scenario == "cov") {
    d_vech <- p * (p + 1) / 2; Z <- matrix(0, nrow = len, ncol = d_vech)
    for (t in seq_len(len)) {
      v_target <- vech(window_list[[1]][t, ] %*% t(window_list[[1]][t, ]))
      v_ref <- rep(0, d_vech)
      for (j in seq_along(ref_indices)) {
        v_ref <- v_ref + vech(window_list[[ref_indices[j]]][t, ] %*% t(window_list[[ref_indices[j]]][t, ])) * norm_weights[j]
      }
      Z[t, ] <- v_target - v_ref
    }
    
  } else if (scenario == "multiparam") {
    Z <- matrix(0, nrow = len, ncol = 2)
    for (t in seq_len(len)) {
      v_target <- c(window_list[[1]][t, 1], window_list[[1]][t, 1]^2)
      v_ref <- c(0, 0)
      for (j in seq_along(ref_indices)) {
        v_ref <- v_ref + c(window_list[[ref_indices[j]]][t, 1], window_list[[ref_indices[j]]][t, 1]^2) * norm_weights[j]
      }
      Z[t, ] <- v_target - v_ref
    }
  }
  return(list(Z = Z, actual_L = actual_L))
}

# ==========================================
# 4. 实际数据分析流程 (完全对齐要求)
# ==========================================
run_real_analysis <- function(scenario = "mean", condition = "Anomaly", L = 12) {
  coords <- matrix(c(116.40,39.98, 116.42,39.93, 116.34,39.93, 116.46,39.94, 116.41,39.88, 116.35,39.88), ncol=2, byrow=T)
  files <- c("Beijing_Aotizhongxin.csv", "Beijing_Dongsi.csv", "Beijing_Guanyuan.csv", 
             "Beijing_Nongzhanguan.csv", "Beijing_Tiantan.csv", "Beijing_Wanshouxigong.csv")
  
  data_list <- lapply(files, function(f) {
    df <- read_csv(file.path("2025_new", f), show_col_types = FALSE)
    mat <- as.matrix(df[, c("PM2.5", "PM10")])
    for(j in 1:2) mat[,j] <- na.approx(mat[,j], na.rm = FALSE)
    return(mat)
  })
  
  n <- nrow(data_list[[1]])
  cp_ref <- c(25, 45) # 始终显示的灰色虚线变点
  
  # 数据处理
  work_data <- data_list
  if(condition == "Anomaly") {
    work_data[[1]] <- apply_shift_rigorous(work_data[[1]], 25, 45, scenario)
  }
  
  # 第一阶段 SNSeg
  target_mat <- work_data[[1]]
  if (scenario == "mean") {
    sn_res <- SNSeg_Multi(target_mat, paras_to_test = "mean", confidence = 0.99, grid_size_scale = 0.15)
  } else if (scenario == "cov") {
    sn_res <- SNSeg_Multi(target_mat, paras_to_test = "covariance", confidence = 0.99, grid_size_scale = 0.1)
  } else {
    sn_res <- SNSeg_Uni(target_mat[,1], paras_to_test = c("mean", "variance"), confidence = 0.99, grid_size_scale = 0.1)
  }
  raw_cps <- sn_res$est_cp
  
  # 第二阶段 DWB 严格校验
  final_cps <- c()
  alpha <- if(condition == "Normal") 0.01 else 0.05
  if (length(raw_cps) > 0) {
    for (cp in raw_cps) {
      feat <- build_Z_rigorous(work_data, cp, L, scenario, coords)
      if (!is.null(feat)) {
        test_out <- dwb_test_cpp(feat$Z, feat$actual_L, M = 1600, bandwidth = 1) # 带宽保持10
        if (test_out$p_value < alpha) final_cps <- c(final_cps, cp)
      }
    }
  }
  
  # --- 绘图与保存 ---
  pdf_name <- sprintf("Final_Analysis_%s_%s.pdf", toupper(scenario), condition)
  pdf(pdf_name, width = 11, height = if(scenario == "multiparam") 6 else 9)
  
  plot_p <- if(scenario == "multiparam") 1 else 2
  par(mfrow = c(plot_p, 1), mar = c(4, 4.5, 3, 1))
  
  # IDW 空间趋势用于绘图
  dists <- sapply(2:6, function(i) haversine_dist(coords[1,1], coords[1,2], coords[i,1], coords[i,2]))
  w <- (1/dists^2)/sum(1/dists^2); v_names <- c("PM2.5", "PM10")
  
  for (j in 1:plot_p) {
    t_s <- work_data[[1]][, j]; r_s <- matrix(0, n, 1)
    for (k in 1:5) r_s <- r_s + work_data[[k+1]][, j] * w[k]
    
    plot(1:n, t_s, type="l", col="#d7301f", lwd=2, ylim=range(c(t_s, r_s)),
         main=sprintf("%s (%s - %s)", v_names[j], toupper(scenario), condition), ylab="Val")
    lines(1:n, r_s, col="#2171b5", lty=2)
    
    # 注入窗口参考线 (始终显示灰色虚线)
    abline(v = cp_ref, col = "gray45", lty = 3, lwd = 1.5)
    
    # 检测到的异常变点 (仅在P显著时显示绿色竖线)
    if (length(final_cps) > 0) abline(v = final_cps, col = "#2ca25f", lwd = 2.5)
    
    if (j == 1) legend("topright", bty="n", cex=0.8, lty=c(1,2,3,1), lwd=c(2,1.2,1.5,2.5),
                       legend=c("Target Site","Spatial Ref","Injected Window","Detected CP"),
                       col=c("#d7301f","#2171b5","gray45","#2ca25f"))
  }
  dev.off()
  cat(sprintf("[SUCCESS] Scenario: %s, Condition: %s. Saved to %s\n", scenario, condition, pdf_name))
}

# ==========================================
# 5. 执行测试
# ==========================================
set.seed(42)
for (sc in c("mean", "cov", "multiparam")) {
  run_real_analysis(sc, "Anomaly")
  run_real_analysis(sc, "Normal")
}