######单变点

library(SNSeg)
library(ggplot2)
library(patchwork)

n <- 500
cp_true <- 250
set.seed(2026) 

v1_pre          <- rnorm(cp_true, mean = 0, sd = 1)
v1_post_normal  <- rnorm(n - cp_true, mean = 0, sd = 1)
v1_post_change  <- rnorm(n - cp_true, mean = 0, sd = 4) # 信号强度 1 -> 4

v2_pre          <- rnorm(cp_true, mean = 0, sd = 1)
v2_post_normal  <- rnorm(n - cp_true, mean = 0, sd = 1)
v2_post_change  <- rnorm(n - cp_true, mean = 0, sd = 4)

s1_v1 <- c(v1_pre, v1_post_change)
s1_v2 <- c(v2_pre, v2_post_normal)

s2_v1 <- c(v1_pre, v1_post_change)
s2_v2 <- c(v2_pre, v2_post_change)

s3_v1 <- c(v1_pre, v1_post_normal)
s3_v2 <- c(v2_pre, v2_post_change)

scenarios <- list(
  S1 = data.frame(Time=1:n, Var1=s1_v1, Var2=s1_v2),
  S2 = data.frame(Time=1:n, Var1=s2_v1, Var2=s2_v2),
  S3 = data.frame(Time=1:n, Var1=s3_v1, Var2=s3_v2)
)


# 噪声过滤函数 (Noise Handling)
filter_noisy_cps <- function(cps, threshold = 30) {
  if (is.null(cps) || length(cps) == 0) return(NULL)
  if (length(cps) == 1) return(cps)
  
  cps <- sort(cps)
  # 将距离小于 threshold 的连续变点归为同一聚类簇
  group <- cumsum(c(1, diff(cps) >= threshold))
  # 取每个簇的中位数作为最终的去噪变点
  robust_cps <- round(tapply(cps, group, median))
  return(unname(robust_cps))
}


process_and_save_perfect <- function(data, filename, title_prefix) {
  
  local_max <- max(abs(c(data$Var1, data$Var2))) * 1.05
  
  res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "variance")
  cp_uni <- filter_noisy_cps(res_uni$est_cp)
  
  res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "covariance")
  cp_multi <- filter_noisy_cps(res_multi$est_cp)
  
 
  create_sub_plot <- function(df, detected_cp, sub_title, show_var2 = FALSE) {
    p <- ggplot(df, aes(x = Time)) +
      geom_line(aes(y = Var1, color = "Var1"), alpha = 0.7) +
      theme_minimal() +
      labs(subtitle = sub_title, y = "观测值", x = "时间") +
      scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "orange")) +
      coord_cartesian(ylim = c(-local_max, local_max)) +
     
      theme(legend.title = element_blank()) +
      geom_vline(xintercept = cp_true, color = "grey", linetype = "dashed", size = 1.2)
    
    if(show_var2) {
      p <- p + geom_line(aes(y = Var2, color = "Var2"), alpha = 0.7)
    }
    
    if(!is.null(detected_cp)) {
      p <- p + geom_vline(xintercept = detected_cp, color = "green", linetype = "solid", size = 0.8)
    }
    
    return(p)
  }
  
  p_left <- create_sub_plot(data, cp_uni, "Var1")
  p_right <- create_sub_plot(data, cp_multi, "Var1和Var2", show_var2 = TRUE)
  
  
  combined_plot <- (p_left | p_right) +
    plot_annotation(title = title_prefix, 
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.spacing.x = unit(1.5, "cm"))
  
  ggsave(filename, combined_plot, width = 12, height = 5, dpi = 300)
}


process_and_save_perfect(scenarios$S1, "Variance_S1_1.png", "仅Var1变化")
process_and_save_perfect(scenarios$S2, "Variance_S2_1.png", "Var1和Var2均发生变化")
process_and_save_perfect(scenarios$S3, "Variance_S3_1.png", "仅Var2变化")





######两个变点
library(SNSeg)
library(ggplot2)
library(patchwork)
n <- 900
cp_true <- c(300, 600) 
set.seed(2026)

vol_low  <- rnorm(300, mean = 0, sd = 1)
vol_high <- rnorm(300, mean = 0, sd = 4)
vol_mid  <- rnorm(300, mean = 0, sd = 1) 

s1_v1 <- c(vol_low, vol_high, vol_mid)
s1_v2 <- rnorm(n, mean = 0, sd = 1) 

s2_v1 <- c(vol_low, vol_high, vol_mid)
s2_v2 <- c(rnorm(300, 0, 1), rnorm(300, 0, 4), rnorm(300, 0, 1))

s3_v1 <- rnorm(n, mean = 0, sd = 1) 
s3_v2 <- c(vol_low, vol_high, vol_mid)

scenarios <- list(
  S1 = data.frame(Time=1:n, Var1=s1_v1, Var2=s1_v2),
  S2 = data.frame(Time=1:n, Var1=s2_v1, Var2=s2_v2),
  S3 = data.frame(Time=1:n, Var1=s3_v1, Var2=s3_v2)
)


# 噪声过滤函数 
filter_noisy_cps <- function(cps, threshold = 40) {
  if (is.null(cps) || length(cps) == 0) return(NULL)
  if (length(cps) == 1) return(cps)
  
  cps <- sort(cps)
 
  group <- cumsum(c(1, diff(cps) >= threshold))
  robust_cps <- round(tapply(cps, group, median))
  return(unname(robust_cps))
}


process_and_save_final <- function(data, filename, title_prefix) {
  
  local_max <- max(abs(c(data$Var1, data$Var2))) * 1.05
  

  res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "variance")
  cp_uni <- filter_noisy_cps(res_uni$est_cp)
  
 
  res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "covariance")
  cp_multi <- filter_noisy_cps(res_multi$est_cp)
  

  create_sub_plot <- function(df, detected_cps, sub_title, show_var2 = FALSE) {
    p <- ggplot(df, aes(x = Time)) +
      geom_line(aes(y = Var1, color = "Var1"), alpha = 0.6) +
      theme_minimal() +
      labs(subtitle = sub_title, y = "观测值", x = "时间") +
      scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "orange")) +
      coord_cartesian(ylim = c(-local_max, local_max)) +
      theme(legend.title = element_blank()) +
      geom_vline(xintercept = cp_true, color = "grey50", linetype = "dashed", size = 1)
    
    if(show_var2) {
      p <- p + geom_line(aes(y = Var2, color = "Var2"), alpha = 0.6)
    }
    
  
    if(!is.null(detected_cps)) {
      p <- p + geom_vline(xintercept = detected_cps, color = "green", linetype = "solid", size = 0.8)
    }
    
    return(p)
  }
  
  p_left <- create_sub_plot(data, cp_uni, "Var1")
  p_right <- create_sub_plot(data, cp_multi, "Var1和Var2", show_var2 = TRUE)
  

  combined_plot <- (p_left | p_right) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = title_prefix, 
      #subtitle = "真实变点 (灰色虚线) vs 算法识别 (绿色实线)",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) &
    theme(legend.position = "bottom")
  
  ggsave(filename, combined_plot, width = 12, height = 5, dpi = 300)
}


process_and_save_final(scenarios$S1, "Variance_S1_2.png", "仅Var1变化")
process_and_save_final(scenarios$S2, "Variance_S2_2.png", "Var1和Var2均发生变化")
process_and_save_final(scenarios$S3, "Variance_S3_2.png", "仅Var2变化")


