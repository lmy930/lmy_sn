
library(SNSeg)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)

Rcpp::sourceCpp("0414.cpp")


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


test_anomaly <- function(data_list, tau, L=40, block_size=5, M=500, scenario="mean") {
  K <- length(data_list)
  p <- ncol(data_list[[1]])
  
  start_idx <- max(1, tau - L + 1)
  end_idx <- min(nrow(data_list[[1]]), tau + L)
  actual_L <- tau - start_idx + 1
  
  if(actual_L < 10) return(NA)
  
  target_data <- data_list[[1]][start_idx:end_idx, , drop=FALSE]
  ref_data <- matrix(0, nrow=(end_idx - start_idx + 1), ncol=p)
  for(i in 2:K) {
    ref_data <- ref_data + data_list[[i]][start_idx:end_idx, , drop=FALSE] / (K - 1)
  }
  
  # 特征映射：将各种场景统一为多元向量的检验
  if (scenario == "mean") {
    Z <- target_data - ref_data
  } else if (scenario == "cov") {
    d_vech <- p * (p + 1) / 2
    Z <- matrix(0, nrow=nrow(target_data), ncol=d_vech)
    for(t in 1:nrow(target_data)) {
      v_target <- vech(tcrossprod(target_data[t, ]))
      v_ref <- rep(0, d_vech)
      for(i in 2:K) {
        v_ref <- v_ref + vech(tcrossprod(data_list[[i]][start_idx+t-1, ])) / (K - 1)
      }
      Z[t, ] <- v_target - v_ref
    }
  } else if (scenario == "multiparam") {
    # 均值+方差映射为 (X, X^2)
    Z <- matrix(0, nrow=nrow(target_data), ncol=2)
    Z[, 1] <- target_data[, 1] - ref_data[, 1]
    Z[, 2] <- (target_data[, 1]^2) - (ref_data[, 1]^2)
  }
  
  # 调用 C++ 做 Block Bootstrap 检验
  p_val <- boot_test_cpp(Z, actual_L, block_size, M)
  return(p_val)
}


plot_spatial_cp <- function(data_list, est_cp, final_pred, title_str, scenario) {
  n <- nrow(data_list[[1]])
  
  #无论哪种场景，只画出第一个维度的数据
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
  
  # 设置画图排版: 3行2列，共张图
  par(mfrow=c(3, 2), mar=c(4, 4, 3, 1))
  
  for(scen in scenarios) {
    p <- ifelse(scen == "multiparam", 1, 3) # 多参数为单变量，其他为3变量
    
    for(cond in conditions) {
      is_anomaly_truth <- (cond == "Anomaly")
      
      # 数据生成
      current_data <- generate_spatial_data(n, p, K_stations=4, cp_loc=cp_true, 
                                            shift_size=2, anomaly_only=is_anomaly_truth, scenario=scen)
      
      # --- 第一阶段: SNCP 变点检测 ---
      if(scen == "mean") {
        sn_res <- SNSeg_Multi(current_data[[1]], paras_to_test = "mean", confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
      } else if (scen == "cov") {
        sn_res <- SNSeg_Multi(current_data[[1]], paras_to_test = "covariance", confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
      } else if (scen == "multiparam") {
        sn_res <- SNSeg_Uni(current_data[[1]][,1], paras_to_test = c("mean", "variance"), confidence = 0.9, grid_size_scale = 0.05, grid_size = 25, plot_SN = FALSE)
      }
      
      est_cps <- sn_res$est_cp
      if(length(est_cps) == 0) {
        best_cp <- NA
        detected <- FALSE
        hausdorff_dist <- NA
      } else {
        hausdorff_dist <- min(abs(est_cps - cp_true))
        detected <- (hausdorff_dist <= 50)
        best_cp <- est_cps[which.min(abs(est_cps - cp_true))]
      }
      
      # --- 第二阶段: Bootstrap 异常分析---
      is_anomaly_detected <- FALSE
      p_value <- NA
      final_pred <- "No CP"
      
      if(detected) {
        p_value <- test_anomaly(current_data, tau=best_cp, L=40, block_size=5, M=500, scenario=scen)
        if(!is.na(p_value) && p_value < alpha) {
          is_anomaly_detected <- TRUE
          final_pred <- "Anomaly"
        } else {
          final_pred <- "Normal"
        }
      }
      
      # --- 第三阶段: 绘图 ---
      title_str <- sprintf("Scenario: %s | Truth: %s\nPred: %s (P=%.3f)", 
                           toupper(scen), cond, final_pred, p_value)
      plot_spatial_cp(current_data, best_cp, final_pred, title_str, scen)
      
      # 记录结果
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
  
  print("==== 模拟实验评价指标输出 ====")
  print(results)
  return(results)
}


sim_results <- run_simulation_and_plot()