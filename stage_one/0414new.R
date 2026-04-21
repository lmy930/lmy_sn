######在0414_new.R的基础上添加特征标准化，消除量纲差异
library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)


Rcpp::sourceCpp("0414new.cpp")


generate_spatial_data <- function(n, p, K_stations, cp_loc, shift_size, anomaly_only=TRUE, scenario="mean") {
  data_list <- list()
  for(i in 1:K_stations) {
    data_list[[i]] <- matrix(0, nrow=n, ncol=p)
  }
  
  # 基础 VAR(1) 过程 + 噪声
  for(t in 2:n) {
    for(i in 1:K_stations) {
      for(dim in 1:p) {
        data_list[[i]][t, dim] <- 0.5 * data_list[[i]][t-1, dim] + rnorm(1, 0, 1)
      }
    }
  }
  
  # 植入变点
  for(t in cp_loc:n) {
    if(scenario == "mean") {
      # 均值平移
      data_list[[1]][t, ] <- data_list[[1]][t, ] + shift_size
      if(!anomaly_only) {
        for(i in 2:K_stations) data_list[[i]][t, ] <- data_list[[i]][t, ] + shift_size
      }
    } else if (scenario == "cov") {
      # 协方差/方差放大
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

# 半向量化函数 (vech) 用于协方差特征映射
vech <- function(mat) {
  mat[lower.tri(mat, diag = TRUE)]
}



test_anomaly <- function(data_list, tau, L=150, block_size=10, M=500, scenario="mean") {
  
  K <- length(data_list)
  p <- ncol(data_list[[1]])
  
  start_idx <- max(1, tau - L + 1)
  end_idx <- min(nrow(data_list[[1]]), tau + L)
  actual_L <- tau - start_idx + 1
  
  if(actual_L < 20) return(NA)
  
  target_data <- data_list[[1]][start_idx:end_idx, , drop=FALSE]
  ref_data_list <- lapply(data_list, function(mat) mat[start_idx:end_idx, , drop=FALSE])
  
  # 分段中心化辅助函数
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
  
  Z <- NULL
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
      for(i in 2:K) {
        v_ref <- v_ref + vech(tcrossprod(ref_centered_list[[i]][t, ])) / (K - 1)
      }
      Z[t, ] <- v_target - v_ref
    }
    
  } else if (scenario == "multiparam") {
    Z <- matrix(0, nrow=nrow(target_data), ncol=2)
    ref_data_avg <- matrix(0, nrow=nrow(target_data), ncol=p)
    for(i in 2:K) ref_data_avg <- ref_data_avg + ref_data_list[[i]] / (K - 1)
    
    Z[, 1] <- target_data[, 1] - ref_data_avg[, 1]
    
    ref_var_avg <- rep(0, nrow(target_data))
    for(i in 2:K) {
      ref_var_avg <- ref_var_avg + (ref_centered_list[[i]][, 1]^2) / (K - 1)
    }
    Z[, 2] <- (target_centered[, 1]^2) - ref_var_avg
  }
  
  # 调用 C++ 做 Block Bootstrap 检验
  p_val <- boot_test_cpp(Z, actual_L, block_size, M)
  return(p_val)
}


# ---------------------------------------------------------
# 3. 绘图可视化函数 (高度还原师兄论文风格)
# ---------------------------------------------------------
plot_spatial_cp <- function(data_list, est_cp, final_pred, title_str, scenario) {
  n <- nrow(data_list[[1]])
  
  
  target_ts <- data_list[[1]][, 1]
  
  # 辅助站点加权平均
  K <- length(data_list)
  ref_ts <- rep(0, n)
  for(i in 2:K) {
    ref_ts <- ref_ts + data_list[[i]][, 1] / (K - 1)
  }
  
  y_min <- min(c(target_ts, ref_ts), na.rm = TRUE)
  y_max <- max(c(target_ts, ref_ts), na.rm = TRUE)
  

  plot(1:n, target_ts, type="l", col="red", lwd=1.2,
       ylim=c(y_min, y_max*1.1), xlab="Time", ylab="Value (Dim 1)",
       main=title_str, cex.main=1.1)
  lines(1:n, ref_ts, type="l", col="blue", lty=2, lwd=1.2)
  
  
  if(!is.na(est_cp) && est_cp > 0 && est_cp <= n) {
    if(final_pred == "Anomaly") {
      
      abline(v = est_cp, col = "red", lwd = 2, lty = 1)
      text(est_cp, y_max*1.05, expression(hat(tau) ~ "(Anomaly)"), col="red", pos=4, cex=1)
    } else if (final_pred == "Normal") {
      
      abline(v = est_cp, col = "green3", lwd = 2, lty = 1)
      text(est_cp, y_max*1.05, expression(hat(tau) ~ "(Normal)"), col="green3", pos=4, cex=1)
    }
  }
  
  legend("topleft", legend=c("Target Station 1", "Weighted Neighbors"),
         col=c("red", "blue"), lty=c(1, 2), lwd=1.5, cex=0.9, bty="n")
}



run_simulation_and_plot <- function() {
  n <- 500
  cp_true <- 250
  alpha <- 0.05
  set.seed(42)
  
  scenarios <- c("mean", "cov", "multiparam")
  conditions <- c("Anomaly", "Normal")
  
  results <- data.frame()
  par(mfrow=c(3, 2), mar=c(4, 4, 3, 1))
  
  for(scen in scenarios) {
    p <- ifelse(scen == "multiparam", 1, 3) 
    
    for(cond in conditions) {
      is_anomaly_truth <- (cond == "Anomaly")
      
      current_data <- generate_spatial_data(n, p, K_stations=4, cp_loc=cp_true, 
                                            shift_size=2, anomaly_only=is_anomaly_truth, scenario=scen)
      
      # 第一阶段: SNCP 变点检测
      if(scen == "mean") {
        sn_res <- SNSeg_Multi(current_data[[1]], paras_to_test = "mean", confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
      } else if (scen == "cov") {
        sn_res <- SNSeg_Multi(current_data[[1]], paras_to_test = "covariance", confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
      } else if (scen == "multiparam") {
        sn_res <- SNSeg_Uni(current_data[[1]][,1], paras_to_test = c("mean", "variance"), confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
      }
      
      est_cps <- sn_res$est_cp
      if(length(est_cps) == 0) {
        best_cp <- NA; detected <- FALSE
      } else {
        detected <- (min(abs(est_cps - cp_true)) <= 50)
        best_cp <- est_cps[which.min(abs(est_cps - cp_true))]
      }
      
      # 第二阶段: Bootstrap 异常分析
      is_anomaly_detected <- FALSE
      p_value <- NA
      final_pred <- "No CP"
      
      if(detected) {
        p_value <- test_anomaly(current_data, tau=best_cp, L=150, block_size=10, M=500, scenario=scen)
        if(!is.na(p_value) && p_value < alpha) {
          is_anomaly_detected <- TRUE
          final_pred <- "Anomaly"
        } else {
          final_pred <- "Normal"
        }
      }
      
      title_str <- sprintf("Scenario: %s | Truth: %s\nPred: %s (P=%.3f)", 
                           toupper(scen), cond, final_pred, p_value)
      plot_spatial_cp(current_data, best_cp, final_pred, title_str, scen)
      
      results <- rbind(results, data.frame(
        Scenario = scen,
        GroundTruth = cond,
        Stage1_Detected = detected,
        Est_CP = best_cp,
        Stage2_Pvalue = p_value,
        Final_Prediction = final_pred
      ))
    }
  }
  par(mfrow=c(1,1))
  print(results)
  return(results)
}

sim_results <- run_simulation_and_plot()