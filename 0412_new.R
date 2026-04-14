
library(SNSeg)
library(mclust) 
library(ecp)         
library(changepoint) 
library(Rcpp)


Rcpp::sourceCpp('D:/lmy/sn_stat.cpp')


eval_metrics <- function(est_cp, true_cp, n) {
  if(any(is.na(est_cp))) est_cp <- numeric(0)
  m_hat <- length(est_cp)
  m_o <- length(true_cp)
  

  if (m_hat == 0 && m_o == 0) {
    d1 <- 0; d2 <- 0; dH <- 0
  } else if (m_hat == 0) { 
    d1 <- 0; d2 <- n; dH <- n
  } else if (m_o == 0) {  
    d1 <- n; d2 <- 0; dH <- n
  } else {
   
    d1 <- max(sapply(est_cp, function(x) min(abs(x - true_cp))))
   
    d2 <- max(sapply(true_cp, function(x) min(abs(x - est_cp))))
    dH <- max(d1, d2)
  }
  
 
  cp_to_lab <- function(p, total) {
    pts <- unique(sort(p[p > 0 & p < total]))
    rep(1:(length(pts)+1), times = diff(c(0, pts, total)))
  }
  ari <- mclust::adjustedRandIndex(cp_to_lab(est_cp, n), cp_to_lab(true_cp, n))
  
  
  return(c(d1 = d1, d2 = d2, dH = dH, delta_m = m_hat - m_o, ARI = ari))
}

summarize <- function(mat, m_name, meth_name) {
  avg <- colMeans(mat, na.rm=T)
  return(data.frame(Model=m_name, Method=meth_name, 
                    d1=avg[1], d2=avg[2], dH=avg[3], Delta_m=avg[4], ARI=avg[5]))
}


# 1. 严谨版 SN+BS 
SN_BS_Rigorous <- function(Y, offset=0, cv_threshold=31.87, min_size=30) {
  n <- length(Y)
  if (n < 2 * min_size) return(numeric(0))
  
  
  T_stats <- calc_global_sn_stat_cpp(Y) 
  search_idx <- min_size:(n - min_size)
  
  max_k <- search_idx[which.max(T_stats[search_idx])]
  if (T_stats[max_k] > cv_threshold) {
    cp <- offset + max_k
    # 严格的二分递归
    return(sort(c(cp, 
                  SN_BS_Rigorous(Y[1:max_k], offset, cv_threshold, min_size),
                  SN_BS_Rigorous(Y[(max_k+1):n], cp, cv_threshold, min_size))))
  } else return(numeric(0))
}

# 2. 多变量均值基准：多元 CUSUM + BS
CUSUM_Multi_BS <- function(Y, threshold=1.5, min_size=30) {
  n <- nrow(Y); d <- ncol(Y)
  find_cp <- function(data, start) {
    nn <- nrow(data)
    if (nn < 2 * min_size) return(NULL)
    stats <- sapply(min_size:(nn - min_size), function(k) {
      sqrt(k * (nn - k) / nn) * sqrt(sum((colMeans(data[1:k,,drop=F]) - colMeans(data[(k+1):nn,,drop=F]))^2))
    })
    if (max(stats) > threshold) {
      best_k <- which.max(stats) + min_size - 1
      return(best_k + start)
    }
    return(NULL)
  }
  
  all_cp <- NULL
  recursive_bs <- function(data, start) {
    cp <- find_cp(data, start)
    if (!is.null(cp)) {
      all_cp <<- c(all_cp, cp)
      recursive_bs(data[1:(cp-start),,drop=F], start)
      recursive_bs(data[(cp-start+1):nrow(data),,drop=F], cp)
    }
  }
  recursive_bs(Y, 0)
  return(sort(all_cp))
}

# 3. 多变量协方差矩阵基准：AHHR 
AHHR_BS <- function(Y, threshold=2.0) {
  n <- nrow(Y); d <- ncol(Y)
  Z <- matrix(0, n, d*(d+1)/2)
  for(i in 1:n) {
    op <- Y[i,] %o% Y[i,]
    Z[i,] <- op[upper.tri(op, diag=TRUE)]
  }
  return(CUSUM_Multi_BS(Z, threshold=threshold))
}

# =====================================================================
# 模块三：数据生成过程 (DGP)
# =====================================================================

gen_data_v2 <- function(type, n, rho, d=1) {
  Y <- SNSeg::MAR(n, reptime=d, rho=rho)
  if(d == 1) Y <- as.numeric(Y)
  
  if (type == "M4") { 
    Y[1:(n/2)] <- Y[1:(n/2)] + 0.8; Y[(n*0.75+1):n] <- Y[(n*0.75+1):n] + 0.8
    return(list(Y=Y, cpts=c(n/2, n*0.75)))
  } else if (type == "M5") { 
    Y[(n/2+1):(n*0.75)] <- Y[(n/2+1):(n*0.75)] + 0.8; Y[(n*0.75+1):n] <- Y[(n*0.75+1):n] + 1.6
    return(list(Y=Y, cpts=c(n/2, n*0.75)))
  } else if (type == "M1") { 
    Y[101:200] <- Y[101:200] + 2; Y[301:400] <- Y[301:400] + 2; Y[501:600] <- Y[501:600] + 2
    return(list(Y=Y, cpts=c(100, 200, 300, 400, 500)))
  } else if (type == "MV_Multi_Mean") { 
    Y[(n/2+1):n, ] <- Y[(n/2+1):n, ] + 1.5
    return(list(Y=Y, cpts=c(n/2)))
  } else if (type == "MV_Multi_Cov") {  
    Y[(n/2+1):n, ] <- Y[(n/2+1):n, ] * 2.0 
    return(list(Y=Y, cpts=c(n/2)))
  } else if (type == "MV_Single") { 
    Y[(n/2+1):n] <- 1 + 1.5 * Y[(n/2+1):n]
    return(list(Y=Y, cpts=c(n/2)))
  }
}


exp_1 <- function(iter) {
  res <- data.frame()
  for (m in c("M4", "M5")) {
    m_sncp <- m_snbs <- matrix(0, iter, 5)
    for (i in 1:iter) {
      data <- gen_data_v2(m, n=2000, rho=0.7)
      f <- SNSeg::SNSeg_Uni(data$Y, paras_to_test="mean")
      m_sncp[i,] <- eval_metrics(f$est_cp, data$cpts, 2000)
      m_snbs[i,] <- eval_metrics(SN_BS_Rigorous(data$Y), data$cpts, 2000)
    }
    res <- rbind(res, summarize(m_sncp, m, "SNCP"), summarize(m_snbs, m, "SN+BS (Rigorous)"))
  }
  return(res)
}

exp_2 <- function(iter) {
  res <- data.frame()
  m_sncp <- m_cusum <- matrix(0, iter, 5)
  for (i in 1:iter) {
    data <- gen_data_v2("M1", n=600, rho=0.2)
    f_sncp <- SNSeg::SNSeg_Uni(data$Y, paras_to_test="mean")
    m_sncp[i,] <- eval_metrics(f_sncp$est_cp, data$cpts, 600)
    
    f_cusum <- changepoint::cpt.mean(data$Y, method="BinSeg", Q=10)
    m_cusum[i,] <- eval_metrics(changepoint::cpts(f_cusum), data$cpts, 600)
  }
  rbind(summarize(m_sncp, "M1(Mean)", "SNCP"), summarize(m_cusum, "M1(Mean)", "CUSUM+BS"))
}

exp_3 <- function(iter) {
  n <- 1000; res <- data.frame()
  d_m <- 5; d_c <- 4 
  m_sncp_m <- m_cusum_m <- m_sncp_c <- m_ahhr_c <- matrix(0, iter, 5)
  
  for (i in 1:iter) {
    dm <- gen_data_v2("MV_Multi_Mean", n, rho=0.5, d=d_m)
    fm <- SNSeg::SNSeg_Multi(dm$Y, paras_to_test="mean")
    m_sncp_m[i,] <- eval_metrics(fm$est_cp, dm$cpts, n)
    m_cusum_m[i,] <- eval_metrics(CUSUM_Multi_BS(dm$Y), dm$cpts, n)
    
    dc <- gen_data_v2("MV_Multi_Cov", n, rho=0.5, d=d_c)
    fc <- SNSeg::SNSeg_Multi(dc$Y, paras_to_test="covariance")
    m_sncp_c[i,] <- eval_metrics(fc$est_cp, dc$cpts, n)
    m_ahhr_c[i,] <- eval_metrics(AHHR_BS(dc$Y), dc$cpts, n)
  }
  res <- rbind(summarize(m_sncp_m, paste0("Multi_Mean(d=",d_m,")"), "SNCP_Multi"), 
               summarize(m_cusum_m, paste0("Multi_Mean(d=",d_m,")"), "CUSUM_Multi"),
               summarize(m_sncp_c, paste0("Multi_Cov(d=",d_c,")"), "SNCP_Multi"), 
               summarize(m_ahhr_c, paste0("Multi_Cov(d=",d_c,")"), "AHHR_BS"))
  return(res)
}

exp_4 <- function(iter) {
  n <- 1000; m_sncp <- m_ecf <- matrix(0, iter, 5)
  for (i in 1:iter) {
    data <- gen_data_v2("MV_Single", n, rho=0.5)
    f_s <- SNSeg::SNSeg_Uni(data$Y, paras_to_test=c("mean", "variance"))
    m_sncp[i,] <- eval_metrics(f_s$est_cp, data$cpts, n)
    
    f_e <- ecp::e.divisive(as.matrix(data$Y), sig.lvl=0.05, min.size=30)
    m_ecf[i,] <- eval_metrics(f_e$estimates[c(-1, -length(f_e$estimates))], data$cpts, n)
  }
  rbind(summarize(m_sncp, "Single_Mean+Var", "SNCP_Uni"), 
        summarize(m_ecf, "Single_Mean+Var", "ECF(e.divisive)"))
}


set.seed(7)
it <- 100

cat("\n==========================================================\n")
cat("--- [结论一]：单变量多变点场景下 SN+BS 失效性验证 ---\n")
print(exp_1(it))

cat("\n--- [结论二]：SNCP 均值检测性能对比 CUSUM+BS ---\n")
print(exp_2(it))

cat("\n--- [结论三]：多变量均值(d=5)与协方差(d=4)检测性能对比 ---\n")
print(exp_3(it))

cat("\n--- [结论四]：单变量多参数(均值+方差)性能对比 ECF 方法 ---\n")
print(exp_4(it))
cat("==========================================================\n")