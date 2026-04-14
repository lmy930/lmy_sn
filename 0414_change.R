######在0414new.R的基础上更改数据生成方式，调用R包中的MAR函数
library(SNSeg)
library(Rcpp)
library(RcppArmadillo)
library(dplyr)


cpp_code <- "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 计算多维特征差的 L2 范数平方
double calc_stat(const mat& Z, int L) {
  rowvec mean_before = mean(Z.rows(0, L - 1), 0);
  rowvec mean_after  = mean(Z.rows(L, 2 * L - 1), 0);
  return accu(square(mean_after - mean_before));
}

// [[Rcpp::export]]
double boot_test_cpp(mat Z, int L, int block_size, int M) {
  int n = 2 * L;
  int d = Z.n_cols;
  
  // 特征标准化 (Feature Standardization): 消除不同参数的量纲差异
  mat Z_norm = Z;
  for(int j = 0; j < d; ++j) {
    double col_sd = stddev(Z.col(j));
    if(col_sd > 1e-6) {
      Z_norm.col(j) /= col_sd;
    }
  }
  
  double S_obs = calc_stat(Z_norm, L);
  
  // 强加原假设 (Null-adjustment)
  rowvec mean_before = mean(Z_norm.rows(0, L - 1), 0);
  rowvec mean_after  = mean(Z_norm.rows(L, n - 1), 0);
  rowvec diff = mean_after - mean_before;
  
  mat Z_tilde = Z_norm;
  for(int t = L; t < n; ++t) {
    Z_tilde.row(t) -= diff;
  }
  
  // MBB 抽样
  int n_blocks = std::ceil((double)n / block_size);
  int exceed_count = 0;
  
  for(int m = 0; m < M; ++m) {
    mat Z_star(n, d, fill::zeros);
    int current_idx = 0;
    
    for(int b = 0; b < n_blocks; ++b) {
      int start_idx = R::runif(0, n - block_size + 1); 
      for(int k = 0; k < block_size; ++k) {
        if (current_idx < n) {
          Z_star.row(current_idx) = Z_tilde.row(start_idx + k);
          current_idx++;
        }
      }
    }
    
    double S_star = calc_stat(Z_star, L);
    if(S_star > S_obs) exceed_count++;
  }
  
  return (double)exceed_count / M;
}
"

Rcpp::sourceCpp(code = cpp_code)


generate_spatial_data_pkg <- function(n, p, K_stations, cp_loc, shift_size, anomaly_only=TRUE, scenario="mean") {
  data_list <- list()
  
  
  common_factor <- MAR(n, reptime = p, rho = 0.5) 
  
  for(i in 1:K_stations) {
    station_noise <- MAR(n, reptime = p, rho = 0.5)
    # 结合公共因子和独立噪声，构建空间相关性
    data_list[[i]] <- 0.7 * common_factor + 0.3 * station_noise
  }
  
  # 植入变点特征
  for(t in cp_loc:n) {
    if(scenario == "mean") {
      data_list[[1]][t, ] <- data_list[[1]][t, ] + shift_size
      if(!anomaly_only) {
        for(i in 2:K_stations) data_list[[i]][t, ] <- data_list[[i]][t, ] + shift_size
      }
    } else if (scenario == "cov") {
      # 放大协方差 (乘以系数)
      data_list[[1]][t, ] <- data_list[[1]][t, ] * shift_size
      if(!anomaly_only) {
        for(i in 2:K_stations) data_list[[i]][t, ] <- data_list[[i]][t, ] * shift_size
      }
    } else if (scenario == "multiparam") {
      # 均值和方差同时改变
      data_list[[1]][t, ] <- data_list[[1]][t, ] * shift_size + shift_size
      if(!anomaly_only) {
        for(i in 2:K_stations) data_list[[i]][t, ] <- data_list[[i]][t, ] * shift_size + shift_size
      }
    }
  }
  return(data_list)
}


vech <- function(mat) { mat[lower.tri(mat, diag = TRUE)] }

test_anomaly <- function(data_list, tau, L=100, block_size=10, M=500, scenario="mean") {
  K <- length(data_list)
  p <- ncol(data_list[[1]])
  
  start_idx <- max(1, tau - L + 1)
  end_idx <- min(nrow(data_list[[1]]), tau + L)
  actual_L <- tau - start_idx + 1
  if(actual_L < 20) return(NA)
  
  target_data <- data_list[[1]][start_idx:end_idx, , drop=FALSE]
  ref_data_list <- lapply(data_list, function(mat) mat[start_idx:end_idx, , drop=FALSE])
  
  # 局部中心化
  local_center <- function(mat, tau_local) {
    mat_centered <- mat
    mu_before <- colMeans(mat[1:tau_local, , drop=FALSE])
    mu_after  <- colMeans(mat[(tau_local+1):nrow(mat), , drop=FALSE])
    for(t in 1:tau_local) mat_centered[t, ] <- mat_centered[t, ] - mu_before
    for(t in (tau_local+1):nrow(mat)) mat_centered[t, ] <- mat_centered[t, ] - mu_after
    return(mat_centered)
  }
  
  target_centered <- local_center(target_data, actual_L)
  ref_centered_list <- lapply(ref_data_list, local_center, tau_local=actual_L)
  
  # 场景特征映射
  if (scenario == "mean") {
    ref_data_avg <- matrix(0, nrow=nrow(target_data), ncol=p)
    for(i in 2:K) ref_data_avg <- ref_data_avg + ref_data_list[[i]] / (K - 1)
    Z <- target_data - ref_data_avg
    
  } else if (scenario == "cov") {
    d_vech <- p * (p + 1) / 2
    Z <- matrix(0, nrow=nrow(target_centered), ncol=d_vech)
    for(t in 1:nrow(target_centered)) {
      v_target <- vech(tcrossprod(target_centered[t, ]))
      v_ref <- rep(0, d_vech)
      for(i in 2:K) v_ref <- v_ref + vech(tcrossprod(ref_centered_list[[i]][t, ])) / (K - 1)
      Z[t, ] <- v_target - v_ref
    }
    
  } else if (scenario == "multiparam") {
    Z <- matrix(0, nrow=nrow(target_data), ncol=2)
    ref_data_avg <- matrix(0, nrow=nrow(target_data), ncol=p)
    for(i in 2:K) ref_data_avg <- ref_data_avg + ref_data_list[[i]] / (K - 1)
    Z[, 1] <- target_data[, 1] - ref_data_avg[, 1]
    
    ref_var_avg <- rep(0, nrow(target_data))
    for(i in 2:K) ref_var_avg <- ref_var_avg + (ref_centered_list[[i]][, 1]^2) / (K - 1)
    Z[, 2] <- (target_centered[, 1]^2) - ref_var_avg
  }
  
  # 调用 C++
  return(boot_test_cpp(Z, actual_L, block_size, M))
}


run_mc_simulations <- function(N_sim = 5, n = 500, cp_true = 250) {
  
  cat(sprintf("开始蒙特卡洛模拟，模拟次数: %d 次...\n", N_sim))
  alpha <- 0.05
  scenarios <- c("mean", "cov", "multiparam")
  conditions <- c("Anomaly", "Normal")
  
  all_results <- list()
  
  for(scen in scenarios) {
    p <- ifelse(scen == "multiparam", 1, 3) 
    
    for(cond in conditions) {
      cat(sprintf("当前测试: 场景 [%s] - 状态 [%s]\n", scen, cond))
      is_anomaly_truth <- (cond == "Anomaly")
      
      # 记录单次模拟状态
      sim_metrics <- data.frame()
      
      for(sim in 1:N_sim) {
        # 1. 数据生成
        data_list <- generate_spatial_data_pkg(n, p, K_stations=4, cp_loc=cp_true, 
                                               shift_size=2.0, anomaly_only=is_anomaly_truth, scenario=scen)
        
        # 2. 阶段一：SNCP 变点检测
        invisible(capture.output({
          if(scen == "mean") {
            sn_res <- SNSeg_Multi(data_list[[1]], paras_to_test = "mean", confidence = 0.9, grid_size_scale = 0.05, plot_SN = FALSE)
          } else if (scen == "cov") {
            sn_res <- SNSeg_Multi(data_list[[1]], paras_to_test = "covariance", confidence = 0.9, grid_size_scale = 0.05, plot_SN = FALSE)
          } else if (scen == "multiparam") {
            sn_res <- SNSeg_Uni(data_list[[1]][,1], paras_to_test = c("mean", "variance"), confidence = 0.9, grid_size_scale = 0.05, plot_SN = FALSE)
          }
        }))
        
        est_cps <- sn_res$est_cp
        if(length(est_cps) == 0) {
          best_cp <- NA; detected <- FALSE; h_dist <- NA
        } else {
          h_dist <- min(abs(est_cps - cp_true))
          detected <- (h_dist <= 50) # 容差阈值设定为 50
          best_cp <- est_cps[which.min(abs(est_cps - cp_true))]
        }
        
        # 3. 阶段二：异常检验
        p_val <- NA
        pred_anomaly <- NA
        if(detected) {
          p_val <- test_anomaly(data_list, tau=best_cp, L=100, block_size=10, M=300, scenario=scen)
          if(!is.na(p_val)) {
            pred_anomaly <- (p_val < alpha)
          }
        }
        
        sim_metrics <- rbind(sim_metrics, data.frame(
          Detected = detected,
          HDist = h_dist,
          P_Val = p_val,
          Pred_Anomaly = pred_anomaly
        ))
      }
      
   
      # 计算该组合下的汇总评价指标
  
      stage1_hit_rate <- mean(sim_metrics$Detected, na.rm = TRUE) * 100
      avg_hdist <- mean(sim_metrics$HDist[sim_metrics$Detected], na.rm = TRUE)
      
      # 阶段二条件准确率 (仅在 Stage 1 成功检测时计算)
      valid_stage2 <- sim_metrics[!is.na(sim_metrics$Pred_Anomaly), ]
      if(nrow(valid_stage2) > 0) {
        if(cond == "Anomaly") {
          stage2_acc <- mean(valid_stage2$Pred_Anomaly == TRUE) * 100
        } else {
          stage2_acc <- mean(valid_stage2$Pred_Anomaly == FALSE) * 100
        }
      } else {
        stage2_acc <- NA
      }
      
      # 端到端总准确率 (Pipeline)
      pipeline_acc <- stage1_hit_rate * (stage2_acc / 100)
      
      all_results[[paste(scen, cond, sep="_")]] <- data.frame(
        Scenario = toupper(scen),
        GroundTruth = cond,
        Stage1_HitRate_Pct = round(stage1_hit_rate, 2),
        Hausdorff_Err = round(avg_hdist, 2),
        Stage2_Accuracy_Pct = round(stage2_acc, 2),
        Pipeline_Accuracy_Pct = round(pipeline_acc, 2)
      )
    }
  }
  
  final_df <- do.call(rbind, all_results)
  rownames(final_df) <- NULL
  return(final_df)
}

final_summary <- run_mc_simulations(N_sim = 20)

cat("\n================================================================\n")
cat("               多重复杂时空异常变点检测算法 - 蒙特卡洛评估报告          \n")
cat("================================================================\n")
print(final_summary)