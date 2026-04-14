
library(SNSeg)
library(strucchange)
library(mclust) # 新增：用于计算 ARI (Adjusted Rand Index)
library(Rcpp)
Rcpp::sourceCpp('D:/lmy/sn_stat.cpp')


cp_to_labels <- function(cp, n) {
  
  cp <- unique(sort(cp[cp > 0 & cp < n]))
  cp_padded <- c(0, cp, n)
  labels <- rep(1:(length(cp_padded) - 1), times = diff(cp_padded))
  return(labels)
}


eval_metrics <- function(est_cp, true_cp, n) {

  if(any(is.na(est_cp))) est_cp <- numeric(0)
  
  m_hat <- length(est_cp)
  m_o <- length(true_cp)
  delta_m <- m_hat - m_o # 变点数量估计误差
  
  # 计算距离类指标 (dH, d1, d2)
  if (m_hat == 0 && m_o == 0) {
    d1 <- 0; d2 <- 0; dH <- 0
  } else if (m_hat == 0) {
    d1 <- 0; d2 <- n; dH <- n # 严重欠估计
  } else if (m_o == 0) {
    d1 <- n; d2 <- 0; dH <- n # 严重过估计
  } else {
   
    d1 <- max(sapply(est_cp, function(x) min(abs(x - true_cp))))
    
    d2 <- max(sapply(true_cp, function(x) min(abs(x - est_cp))))
  
    dH <- max(d1, d2)
  }
  
  # 计算 ARI 
  labels_est <- cp_to_labels(est_cp, n)
  labels_true <- cp_to_labels(true_cp, n)
  ari <- mclust::adjustedRandIndex(labels_est, labels_true)
  

  return(c(dH = dH, d1 = d1, d2 = d2, delta_m = delta_m, ARI = ari))
}


SN_BS_fast <- function(Y, offset=0, cv_threshold=31.87, min_size=30) {
  n <- length(Y)
  if (n < 2 * min_size) return(numeric(0))
  T_stats <- calc_global_sn_stat_cpp(Y)
  search_idx <- min_size:(n - min_size)
  max_k <- search_idx[which.max(T_stats[search_idx])]
  if (T_stats[max_k] > cv_threshold) {
    cp <- offset + max_k
    return(sort(c(cp, 
                  SN_BS_fast(Y[1:max_k], offset, cv_threshold, min_size),
                  SN_BS_fast(Y[(max_k+1):n], cp, cv_threshold, min_size))))
  } else return(numeric(0))
}


CUSUM_Multi_BS <- function(Y, offset=0, threshold=1.5, min_size=30) {
  n <- nrow(Y)
  if (n < 2 * min_size) return(numeric(0))
  stats <- numeric(n)
  for (k in min_size:(n - min_size)) {
    mean1 <- colMeans(Y[1:k, , drop = FALSE])
    mean2 <- colMeans(Y[(k+1):n, , drop = FALSE])
    stats[k] <- sqrt(k * (n - k) / n) * sqrt(sum((mean1 - mean2)^2))
  }
  max_k <- which.max(stats)
  if (stats[max_k] > threshold) {
    cp <- offset + max_k
    return(sort(c(cp, 
                  CUSUM_Multi_BS(Y[1:max_k, , drop=FALSE], offset, threshold, min_size),
                  CUSUM_Multi_BS(Y[(max_k+1):n, , drop=FALSE], cp, threshold, min_size))))
  } else return(numeric(0))
}


gen_data <- function(type, n, rho, d=1) {
  Y_raw <- SNSeg::MAR(n, reptime=d, rho=rho)
  if(d == 1) Y <- as.numeric(Y_raw) else Y <- Y_raw
  shift <- if(d == 1) 1 else 1/sqrt(d)
  
  if (type == "M1") {
    if(d==1) { Y[101:200]<-Y[101:200]+2; Y[301:400]<-Y[301:400]+2; Y[501:600]<-Y[501:600]+2
    } else { Y[101:200,]<-Y[101:200,]+2*shift; Y[301:400,]<-Y[301:400,]+2*shift; Y[501:600,]<-Y[501:600,]+2*shift }
    return(list(Y=Y, cpts=c(100, 200, 300, 400, 500)))
  } else if (type == "M2") {
    s <- 3*shift
    if(d==1) { Y[1:75]<-Y[1:75]-s; Y[376:425]<-Y[376:425]+s; Y[526:575]<-Y[526:575]-s
    } else { Y[1:75,]<-Y[1:75,]-s; Y[376:425,]<-Y[376:425,]+s; Y[526:575,]<-Y[526:575,]-s }
    return(list(Y=Y, cpts=c(75, 375, 425, 525, 575)))
  } else if (type == "M3") {
    s <- 0.4*shift
    if(d==1) { Y[1:1000]<-Y[1:1000]+s; Y[1501:2000]<-Y[1501:2000]+s
    } else { Y[1:1000,]<-Y[1:1000,]+s; Y[1501:2000,]<-Y[1501:2000,]+s }
    return(list(Y=Y, cpts=c(1000, 1500)))
  } else if (type == "M4") {
    Y[1:1000] <- Y[1:1000] + 0.8; Y[1501:2000] <- Y[1501:2000] + 0.8
    return(list(Y=Y, cpts=c(1000, 1500)))
  } else if (type == "M5") {
    Y[1001:1500] <- Y[1001:1500] + 0.8; Y[1501:2000] <- Y[1501:2000] + 1.6
    return(list(Y=Y, cpts=c(1000, 1500)))
  } else if (type == "MV") {
    Y[513:n] <- 1 + 1.5 * Y[513:n]
    return(list(Y=Y, cpts=c(512)))
  }
}


summarize_metrics <- function(metrics_matrix, model_name, method_name) {
  avg_res <- colMeans(metrics_matrix, na.rm = TRUE)
  return(data.frame(
    Model = model_name, Method = method_name,
    dH = avg_res[1], d1 = avg_res[2], d2 = avg_res[3], 
    Delta_m = avg_res[4], ARI = avg_res[5]
  ))
}


run_experiment_univariate <- function(iter = 100) {
  models <- list(
    list(name="M1", n=600, rho=0.2), list(name="M2", n=1000, rho=0.5),
    list(name="M3", n=2000, rho=-0.7), list(name="M4", n=2000, rho=0.7),
    list(name="M5", n=2000, rho=0.7)
  )
  final_res <- data.frame()
  
  for (m in models) {
    cat("正在模拟单变量模型:", m$name, "\n")
    m_sncp <- matrix(0, nrow=iter, ncol=5)
    m_cusum <- matrix(0, nrow=iter, ncol=5)
    m_snbs <- matrix(0, nrow=iter, ncol=5)
    
    for (i in 1:iter) {
      data <- gen_data(m$name, m$n, m$rho)
      
      # SNCP
      fit_sncp <- SNSeg_Uni(data$Y, paras_to_test="mean", grid_size_scale=0.05)
      m_sncp[i,] <- eval_metrics(fit_sncp$est_cp, data$cpts, m$n)
      
      # CUSUM+BP
      bp <- tryCatch(breakpoints(data$Y ~ 1, h=0.05), error=function(e) NULL)
      bp_pts <- if(!is.null(bp)) bp$breakpoints else numeric(0)
      m_cusum[i,] <- eval_metrics(bp_pts, data$cpts, m$n)
      
      # SN+BS
      est_snbs <- SN_BS_fast(data$Y)
      m_snbs[i,] <- eval_metrics(est_snbs, data$cpts, m$n)
    }
    
    final_res <- rbind(final_res, summarize_metrics(m_sncp, m$name, "SNCP"))
    final_res <- rbind(final_res, summarize_metrics(m_cusum, m$name, "CUSUM_BP"))
    final_res <- rbind(final_res, summarize_metrics(m_snbs, m$name, "SN_BS"))
  }
  return(final_res)
}


run_experiment_multivariate <- function(iter = 100) {
  models <- list(
    list(name="M1", n=600, rho=0.2), list(name="M2", n=1000, rho=0.5), list(name="M3", n=2000, rho=-0.7)
  )
  final_res <- data.frame()
  d <- 5
  
  for (m in models) {
    cat("正在模拟多变量模型:", m$name, "\n")
    m_sncp <- matrix(0, nrow=iter, ncol=5)
    m_cusum <- matrix(0, nrow=iter, ncol=5)
    
    for (i in 1:iter) {
      data <- gen_data(m$name, m$n, m$rho, d=d)
      
      fit_sncp <- SNSeg_Multi(data$Y, paras_to_test="mean", grid_size_scale=0.05)
      m_sncp[i,] <- eval_metrics(fit_sncp$est_cp, data$cpts, m$n)
      
      est_cusum <- CUSUM_Multi_BS(data$Y, threshold=2.0)
      m_cusum[i,] <- eval_metrics(est_cusum, data$cpts, m$n)
    }
    final_res <- rbind(final_res, summarize_metrics(m_sncp, m$name, "SNCP_Multi"))
    final_res <- rbind(final_res, summarize_metrics(m_cusum, m$name, "CUSUM_Multi"))
  }
  return(final_res)
}


run_experiment_multiparam <- function(iter = 100) {
  rhos <- c(0.5, -0.5)
  final_res <- data.frame()
  n <- 1024
  
  for (r in rhos) {
    cat("正在模拟多参数模型 rho =", r, "\n")
    m_sncp <- matrix(0, nrow=iter, ncol=5)
    
    for (i in 1:iter) {
      data <- gen_data("MV", n, r)
      fit <- SNSeg_Uni(data$Y, paras_to_test=c("mean", "variance"), grid_size_scale=0.05)
      m_sncp[i,] <- eval_metrics(fit$est_cp, data$cpts, n)
    }
    mod_name <- paste0("MV(rho=", r, ")")
    final_res <- rbind(final_res, summarize_metrics(m_sncp, mod_name, "SNCP_Mean_Var"))
  }
  return(final_res)
}


set.seed(7)
it <- 100 

cat("--- 实验一：单变量多变点检测有效性 (M1-M5) ---\n")
result_uni <- run_experiment_univariate(iter = it)
print(result_uni)

cat("\n--- 实验二：多变量多变点检测可行性 (d=5, M1-M3) ---\n")
result_multi <- run_experiment_multivariate(iter = it)
print(result_multi)

cat("\n--- 实验三：单变量多参数变点检测可行性 (MV1-MV2) ---\n")
result_param <- run_experiment_multiparam(iter = it)
print(result_param)