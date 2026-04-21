######单变点
library(SNSeg)
library(ggplot2)
library(patchwork)
library(MASS)

n <- 1000
cp_true <- 500
set.seed(2026) 


gen_original_1000 <- function(n, cp, rho_pre, rho_post) {
  sigma_pre <- matrix(c(1, rho_pre, rho_pre, 1), nrow = 2)
  data_pre  <- mvrnorm(cp, mu = c(0, 0), Sigma = sigma_pre)
  
  sigma_post <- matrix(c(1, rho_post, rho_post, 1), nrow = 2)
  data_post  <- mvrnorm(n - cp, mu = c(0, 0), Sigma = sigma_post)
  
  full_matrix <- rbind(data_pre, data_post)
  return(data.frame(Time = 1:n, Var1 = full_matrix[,1], Var2 = full_matrix[,2]))
}


scenarios <- list(
  S1 = gen_original_1000(n, cp_true, 0, 0.9),   
  S2 = gen_original_1000(n, cp_true, 0.8, -0.8), 
  S3 = gen_original_1000(n, cp_true, 0, 0.7)     
)


process_and_save_split <- function(data, scene_id, title_prefix, rho_info) {
  
 
  res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "variance")
  cp_uni <- if(length(res_uni$est_cp) > 0) res_uni$est_cp else NULL
  
  res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "covariance")
  cp_multi <- if(length(res_multi$est_cp) > 0) res_multi$est_cp else NULL

 
  local_max <- max(abs(c(data$Var1, data$Var2))) * 1.1

 
  create_ts_plot <- function(df, detected_cp, sub_title, show_var2 = FALSE) {
    p <- ggplot(df, aes(x = Time)) +
      geom_line(aes(y = Var1, color = "Var1"), alpha = 0.5, linewidth = 0.3) +
      theme_minimal() +
      labs(subtitle = sub_title, y = "观测值", x = "时间") +
      scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "orange")) +
      coord_cartesian(ylim = c(-local_max, local_max)) +
      geom_vline(xintercept = cp_true, color = "grey40", linetype = "dashed", linewidth = 1)
    
    if(show_var2) {
      p <- p + geom_line(aes(y = Var2, color = "Var2"), alpha = 0.5, linewidth = 0.3)
    }
    if(!is.null(detected_cp)) {
      p <- p + geom_vline(xintercept = detected_cp, color = "forestgreen", linetype = "solid", linewidth = 1)
    }
    return(p + theme(legend.position = "none"))
  }

  p_ts_left  <- create_ts_plot(data, cp_uni, "Var1")
  p_ts_right <- create_ts_plot(data, cp_multi, "Var1和Var2", show_var2 = TRUE)

 
  combined_ts <- (p_ts_left | p_ts_right) + 
    plot_annotation(
      title = paste(title_prefix),
      #subtitle = paste("相关系数跳变：", rho_info, " | 目标：左图静默，右图报点"),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  ggsave(paste0(scene_id, ".png"), combined_ts, width = 12, height = 5, dpi = 300)


  data_pre  <- data[data$Time <= cp_true, ]
  data_post <- data[data$Time > cp_true, ]
  
  create_scatter <- function(df, sub_title, color_val) {
    ggplot(df, aes(x = Var1, y = Var2)) +
      geom_point(alpha = 0.3, color = color_val, size = 0.8) +
      geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE, linetype = "dotted") +
      theme_minimal() +
      labs(subtitle = sub_title, x = "Var1", y = "Var2") +
      coord_fixed(xlim = c(-4, 4), ylim = c(-4, 4))
  }

  p_sc_left  <- create_scatter(data_pre,"变点前", "grey40")
  p_sc_right <- create_scatter(data_post,"变点后", "steelblue")

  combined_scat <- (p_sc_left | p_sc_right) + 
    plot_annotation(
      title = paste(title_prefix),
      #subtitle = "左图（独立）与右图（相关）的分布差异是检测的统计学依据",
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  ggsave(paste0(scene_id, "_scatter_1.png"), combined_scat, width = 12, height = 6, dpi = 300)
}


process_and_save_split(scenarios$S1, "S1",  "0 -> 0.9")
process_and_save_split(scenarios$S2, "S2",  "0.8 -> -0.8")
process_and_save_split(scenarios$S3, "S3",  "0 -> 0.7")







######多变点
library(SNSeg)
library(ggplot2)
library(patchwork)
library(MASS)


n <- 1500
cp_true <- c(500, 1000) 
set.seed(2026) 


gen_original_two_cp <- function(n, cp, rho1, rho2, rho3) {
  
  sigma1 <- matrix(c(1, rho1, rho1, 1), nrow = 2)
  data1  <- mvrnorm(cp[1], mu = c(0, 0), Sigma = sigma1)
  
  
  sigma2 <- matrix(c(1, rho2, rho2, 1), nrow = 2)
  data2  <- mvrnorm(cp[2] - cp[1], mu = c(0, 0), Sigma = sigma2)
  
 
  sigma3 <- matrix(c(1, rho3, rho3, 1), nrow = 2)
  data3  <- mvrnorm(n - cp[2], mu = c(0, 0), Sigma = sigma3)
  
  full_matrix <- rbind(data1, data2, data3)
  return(data.frame(Time = 1:n, Var1 = full_matrix[,1], Var2 = full_matrix[,2]))
}


scenarios <- list(
  S1 = gen_original_two_cp(n, cp_true, 0, 0.9, 0),    
  S2 = gen_original_two_cp(n, cp_true, 0.8, -0.8, 0.8), 
  S3 = gen_original_two_cp(n, cp_true, 0, 0.7, 0)     
)


process_and_save_split_two_cp <- function(data, scene_id, title_prefix) {
  

  res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "variance")
  cp_uni <- if(length(res_uni$est_cp) > 0) res_uni$est_cp else NULL
  
  res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "covariance")
  cp_multi <- if(length(res_multi$est_cp) > 0) res_multi$est_cp else NULL

  local_max <- max(abs(c(data$Var1, data$Var2))) * 1.1


  create_ts_plot <- function(df, detected_cp, sub_title, show_var2 = FALSE) {
    p <- ggplot(df, aes(x = Time)) +
      geom_line(aes(y = Var1, color = "Var1"), alpha = 0.5, linewidth = 0.3) +
      theme_minimal() +
      labs(subtitle = sub_title, y = "观测值", x = "时间") +
      scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "orange")) +
      coord_cartesian(ylim = c(-local_max, local_max)) +
      
      geom_vline(xintercept = cp_true, color = "grey40", linetype = "dashed", linewidth = 1)
    
    if(show_var2) {
      p <- p + geom_line(aes(y = Var2, color = "Var2"), alpha = 0.5, linewidth = 0.3)
    }
    
    if(!is.null(detected_cp)) {
      p <- p + geom_vline(xintercept = detected_cp, color = "forestgreen", linetype = "solid", linewidth = 1)
    }
    return(p + theme(legend.position = "none"))
  }

  p_ts_left  <- create_ts_plot(data, cp_uni, "Var1")
  p_ts_right <- create_ts_plot(data, cp_multi, "Var1和Var2", show_var2 = TRUE)

  combined_ts <- (p_ts_left | p_ts_right) + 
    plot_annotation(
      title = paste(title_prefix),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )
  
  ggsave(paste0(scene_id, "_2.png"), combined_ts, width = 12, height = 5, dpi = 300)

 
  data_seg1 <- data[data$Time <= cp_true[1], ]
  data_seg2 <- data[data$Time > cp_true[1] & data$Time <= cp_true[2], ]
  data_seg3 <- data[data$Time > cp_true[2], ]
  
  create_scatter <- function(df, sub_title, color_val) {
    ggplot(df, aes(x = Var1, y = Var2)) +
      geom_point(alpha = 0.3, color = color_val, size = 0.8) +
      geom_smooth(method = "lm", formula = y ~ x, color = "red", se = FALSE, linetype = "dotted") +
      theme_minimal() +
      labs(subtitle = sub_title, x = "Var1", y = "Var2") +
      coord_fixed(xlim = c(-4, 4), ylim = c(-4, 4))
  }

  p_sc_1 <- create_scatter(data_seg1, "正常_1 ", "grey40")
  p_sc_2 <- create_scatter(data_seg2, "变点", "steelblue")
  p_sc_3 <- create_scatter(data_seg3, "正常_2", "darkorange")


  combined_scat <- (p_sc_1 | p_sc_2 | p_sc_3) + 
    plot_annotation(
      title = paste(title_prefix),
      theme = theme(plot.title = element_text(size = 14, face = "bold"))
    )

  ggsave(paste0(scene_id, "_scatter_2.png"), combined_scat, width = 15, height = 5, dpi = 300)
}


process_and_save_split_two_cp(scenarios$S1, "S1", "0 -> 0.9 -> 0")
process_and_save_split_two_cp(scenarios$S2, "S2", "0.8 -> -0.8 -> 0.8")
process_and_save_split_two_cp(scenarios$S3, "S3", "0 -> 0.7 -> 0")
