
library(SNSeg)
library(mclust)
library(changepoint)
library(ecp)
library(Rcpp)
library(RcppArmadillo)
library(MASS) # mvrnorm


Rcpp::sourceCpp("D:/lmy/sn_stat.cpp") 


eval_metrics <- function(est_cp, true_cp, n) {
  est_cp <- as.numeric(est_cp)
  est_cp <- est_cp[est_cp > 5 & est_cp < (n-5)]
  true_cp <- true_cp[true_cp > 0 & true_cp < n]
  
  m_hat <- length(est_cp); m_o <- length(true_cp)
  
  if (m_hat == 0 && m_o == 0) {
    d1 <- 0; d2 <- 0; dH <- 0
  } else if (m_hat == 0 || m_o == 0) {
    d1 <- if(m_hat > 0) 1 else 0
    d2 <- if(m_o > 0) 1 else 0
    dH <- 1
  } else {
    d1 <- max(sapply(est_cp, function(x) min(abs(x - true_cp)))) / n
    d2 <- max(sapply(true_cp, function(x) min(abs(x - est_cp)))) / n
    dH <- max(d1, d2)
  }
  
  cp_to_lab <- function(p, total) {
    pts <- unique(sort(p[p > 0 & p < total]))
    rep(1:(length(pts)+1), times = diff(c(0, pts, total)))
  }
  ari <- adjustedRandIndex(cp_to_lab(est_cp, n), cp_to_lab(true_cp, n))
  
  return(c(m_hat = m_hat, Delta_m = m_hat - m_o, d1 = d1*100, d2 = d2*100, dH = dH*100, ARI = ari))
}


gen_data <- function(type, d=1) {
 
  if(type %in% c("M1", "M2", "M3", "M4", "M5")) {
    if(type == "M1") { n=600; rho=0.2; cpts=c(100,200,300,400,500); val=2/sqrt(d) }
    if(type == "M2") { n=1000; rho=0.5; cpts=c(75,375,425,525,575) }
    if(type == "M3") { n=2000; rho=-0.7; cpts=c(1000,1500); val=0.4/sqrt(d) }
    if(type == "M4") { n=2000; rho=0.7; cpts=c(1000,1500) }
    if(type == "M5") { n=2000; rho=0.7; cpts=c(1000,1500) }
    
    Y <- matrix(0, n, d)
    for(j in 1:d) Y[,j] <- as.numeric(arima.sim(n = n, list(ar = rho), sd = 1))
    
    mu <- matrix(0, n, d)
    if(type == "M1") {
      mu[101:200,] <- val; mu[301:400,] <- val; mu[501:600,] <- val
    } else if(type == "M2") {
      val <- 3/sqrt(d)
      mu[1:75,] <- -val; mu[376:425,] <- val; mu[526:575,] <- -val
    } else if(type == "M3") {
      mu[1:1000,] <- val; mu[1501:2000,] <- val
    } else if(type == "M4") {
      mu[1:1000,] <- 0.8; mu[1501:2000,] <- 0.8
    } else if(type == "M5") {
      mu[1001:1500,] <- 0.8; mu[1501:2000,] <- 1.6
    }
    Y <- Y + mu
    if(d==1) Y <- as.numeric(Y)
    return(list(Y=Y, cpts=cpts))
  }
  

  if(type %in% c("C1", "C2", "C3", "C4")) {
    n = 1000; cpts = c(333, 667)
    Y <- matrix(0, n, d)
    
    Sigma_05 <- matrix(0.5, d, d); diag(Sigma_05) <- 1
    Sigma_02 <- matrix(0.2, d, d); diag(Sigma_02) <- 1
    L0 <- matrix(c(1,1,0,0, 0,0,1,1), nrow=4, ncol=2)
    
    for(t in 1:n) {
      if(type %in% c("C1", "C2")) {
        Ft <- mvrnorm(1, rep(0,2), diag(2))
        et <- mvrnorm(1, rep(0,d), diag(d))
        
        if(type == "C1") {
          if(t <= 333 || t > 667) Y[t,] <- L0 %*% Ft + et
          else Y[t,] <- sqrt(3) * (L0 %*% Ft) + et
        } else {
          if(t <= 333) Y[t,] <- L0 %*% Ft + et
          else if(t <= 667) Y[t,] <- sqrt(3) * (L0 %*% Ft) + et
          else Y[t,] <- 3 * (L0 %*% Ft) + et
        }
      } else { # C3, C4 (VAR 模型)
        if(t == 1) prev <- rep(0, d) else prev <- Y[t-1,]
        if(type == "C3") {
          if(t <= 333 || t > 667) et <- mvrnorm(1, rep(0,d), Sigma_05) * sqrt(2)
          else et <- mvrnorm(1, rep(0,d), Sigma_02)
        } else {
          if(t <= 333) et <- mvrnorm(1, rep(0,d), Sigma_02)
          else if(t <= 667) et <- mvrnorm(1, rep(0,d), Sigma_05) * sqrt(2)
          else et <- mvrnorm(1, rep(0,d), Sigma_05) * 2
        }
        Y[t,] <- 0.3 * prev + et
      }
    }
    return(list(Y=Y, cpts=cpts))
  }
  

  if(type %in% c("MV1", "MV2")) {
    n = 1024; cpts = c(512)
    rho = if(type == "MV1") 0.5 else -0.5
    X <- as.numeric(arima.sim(n = n, list(ar = rho), sd = 1))
    Y <- X
    Y[513:1024] <- 1 + 1.5 * X[513:1024]
    return(list(Y=Y, cpts=cpts))
  }
}


# (A) SN+BS 失效方法
SN_BS_Failing <- function(Y, threshold=141.9) {
  n <- length(Y); cpts <- c()
  recursive_bs <- function(st, ed, offset) {
    if(ed - st < 60) return()
    stats <- calc_global_sn_univariate(Y[st:ed])
    if(max(stats, na.rm=T) > threshold) {
      idx <- which.max(stats)
      cpts <<- c(cpts, idx + offset)
      recursive_bs(1, idx, offset)
      recursive_bs(idx + 1, ed - st + 1, offset + idx)
    }
  }
  recursive_bs(1, n, 0)
  return(sort(cpts))
}

# (B) 多元 CUSUM + BS (用于多维均值)
CUSUM_Multi_BS <- function(Y, threshold=2.5) {
  n <- nrow(Y); cpts <- c()
  recursive_bs <- function(st, ed, offset) {
    if(ed - st < 60) return()
    stats <- calc_cusum_multi(as.matrix(Y[st:ed,]))
    if(max(stats, na.rm=T) > threshold) {
      idx <- which.max(stats)
      cpts <<- c(cpts, idx + offset)
      recursive_bs(1, idx, offset)
      recursive_bs(idx + 1, ed - st + 1, offset + idx)
    }
  }
  recursive_bs(1, n, 0)
  return(sort(cpts))
}

# (C) AHHR + BS (用于协方差阵检测)
# 逻辑：取 Y_t Y_t^T 的上三角拉直为新向量 Z_t，然后跑 CUSUM
AHHR_BS <- function(Y, threshold=3.5) {
  n <- nrow(Y); d <- ncol(Y)
  Z <- matrix(0, n, d*(d+1)/2)
  for(i in 1:n) {
    mat <- Y[i,] %o% Y[i,]
    Z[i,] <- mat[upper.tri(mat, diag=TRUE)]
  }
  return(CUSUM_Multi_BS(Z, threshold))
}


run_experiment <- function(iter = 5) {
  
  cat("===============================================================\n")
  cat("[任务1] 单变量多变点: SNCP vs SN+BS (验证失效) \n")
  cat("模型: M4, M5 (非单调与单调变化)\n")
  
  res_t1 <- list()
  for(m in c("M4", "M5")) {
    tmp <- replicate(iter, {
      d <- gen_data(m)
      cp_sncp <- SNSeg_Uni(d$Y, paras_to_test="mean")$est_cp
      cp_snbs <- SN_BS_Failing(d$Y, 141.9)
      c(eval_metrics(cp_sncp, d$cpts, 2000), eval_metrics(cp_snbs, d$cpts, 2000))
    })
    res_t1[[m]] <- rowMeans(tmp)
  }
  print(do.call(rbind, res_t1))
  
  # -------------------------------------------------------------
  cat("\n===============================================================\n")
  cat("[任务2] 单变量多变点: SNCP vs CUSUM+BS\n")
  cat("模型: M1, M2, M3\n")
  
  res_t2 <- list()
  for(m in c("M1", "M2", "M3")) {
    n <- ifelse(m=="M1", 600, ifelse(m=="M2", 1000, 2000))
    tmp <- replicate(iter, {
      d <- gen_data(m)
      cp_sncp <- SNSeg_Uni(d$Y, paras_to_test="mean")$est_cp
      fit_cusum <- cpt.mean(d$Y, method="BinSeg", Q=10)
      cp_cusum <- cpts(fit_cusum)
      c(eval_metrics(cp_sncp, d$cpts, n), eval_metrics(cp_cusum, d$cpts, n))
    })
    res_t2[[m]] <- rowMeans(tmp)
  }
  print(do.call(rbind, res_t2))
  
  # -------------------------------------------------------------
  cat("\n===============================================================\n")
  cat("[任务3] 多变量多变点检测 (d=4)\n")
  cat("(A) 均值检测: SNCP vs CUSUM_Multi | 模型: M1, M2, M3\n")
  cat("(B) 协方差检测: SNCP vs AHHR | 模型: C1, C2, C3, C4\n")
  
  res_t3a <- list()
  for(m in c("M1", "M2", "M3")) {
    n <- ifelse(m=="M1", 600, ifelse(m=="M2", 1000, 2000))
    tmp <- replicate(iter, {
      d <- gen_data(m, d=4)
      cp_sncp <- SNSeg_Multi(d$Y, paras_to_test="mean")$est_cp
      cp_cusum <- CUSUM_Multi_BS(d$Y)
      c(eval_metrics(cp_sncp, d$cpts, n), eval_metrics(cp_cusum, d$cpts, n))
    })
    res_t3a[[m]] <- rowMeans(tmp)
  }
  print("--- 均值检测结果 ---")
  print(do.call(rbind, res_t3a))
  
  res_t3b <- list()
  for(m in c("C1", "C2", "C3", "C4")) {
    tmp <- replicate(iter, {
      d <- gen_data(m, d=4)
      cp_sncp <- SNSeg_Multi(d$Y, paras_to_test="covariance")$est_cp
      cp_ahhr <- AHHR_BS(d$Y)
      c(eval_metrics(cp_sncp, d$cpts, 1000), eval_metrics(cp_ahhr, d$cpts, 1000))
    })
    res_t3b[[m]] <- rowMeans(tmp)
  }
  print("--- 协方差检测结果 ---")
  print(do.call(rbind, res_t3b))
  
  # -------------------------------------------------------------
  cat("\n===============================================================\n")
  cat("[任务4] 单变量多参数(均值+方差): SNCP vs ECF\n")
  cat("模型: MV1, MV2\n")
  
  res_t4 <- list()
  for(m in c("MV1", "MV2")) {
    tmp <- replicate(iter, {
      d <- gen_data(m)
      cp_sncp <- SNSeg_Uni(d$Y, paras_to_test=c("mean", "variance"))$est_cp
      fit_ecp <- e.divisive(as.matrix(d$Y), sig.lvl=0.05, min.size=30)
      cp_ecp <- fit_ecp$estimates[c(-1, -length(fit_ecp$estimates))]
      c(eval_metrics(cp_sncp, d$cpts, 1024), eval_metrics(cp_ecp, d$cpts, 1024))
    })
    res_t4[[m]] <- rowMeans(tmp)
  }
  print(do.call(rbind, res_t4))
  
  cat("===============================================================\n")
}


set.seed(7)
run_experiment(iter = 100)