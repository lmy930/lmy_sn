# ######单独保存
# library(SNSeg)
# library(ggplot2)

# n <- 300
# cp_true <- 150
# set.seed(2026)


# s1_v1 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))
# s1_v2 <- rnorm(n, 0, 1)
# s2_v1 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))
# s2_v2 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))
# s3_v1 <- rnorm(n, 0, 1)
# s3_v2 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))

# scenarios <- list(
#   S1 = data.frame(Time=1:n, Var1=s1_v1, Var2=s1_v2),
#   S2 = data.frame(Time=1:n, Var1=s2_v1, Var2=s2_v2),
#   S3 = data.frame(Time=1:n, Var1=s3_v1, Var2=s3_v2)
# )

# 逻辑
# process_and_save_separate <- function(data, prefix_name, main_title) {
  
#   res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "mean")
#   cp_uni <- if(length(res_uni$est_cp) > 0) res_uni$est_cp[1] else NULL
  
#   res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "mean")
#   cp_multi <- if(length(res_multi$est_cp) > 0) res_multi$est_cp[1] else NULL
  

#   base_plot <- function(df, detected_cp, title_text, show_var2) {
#     p <- ggplot(df, aes(x = Time)) +
#       geom_line(aes(y = Var1, color = "Var1"), alpha = 0.7) +
#       theme_minimal() +
#       labs(title = title_text, y = "观测值", x = "时间") +
#       scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "orange")) +
#       theme(legend.title = element_blank(), 
#             legend.position = if(show_var2) "top" else "none",
#             plot.title = element_text(face = "bold", size = 12),legend.spacing.x = unit(1.5, "cm")) +
#       geom_vline(xintercept = cp_true, color = "grey", linetype = "dashed", size = 1.2)
    
#     if(show_var2) p <- p + geom_line(aes(y = Var2, color = "Var2"), alpha = 0.7)
#     if(!is.null(detected_cp)) p <- p + geom_vline(xintercept = detected_cp, color = "green", linetype = "solid", size = 0.8)
    
#     return(p)
#   }
  
  
#   p_uni <- base_plot(data, cp_uni, paste0(main_title, " - Var1"), FALSE)
#   ggsave(paste0(prefix_name, "_Univariate.png"), p_uni, width = 8, height = 5, dpi = 300)
  

#   p_multi <- base_plot(data, cp_multi, paste0(main_title, " - Var1和Var2"), TRUE)
#   ggsave(paste0(prefix_name, "_Multivariate.png"), p_multi, width = 8, height = 5.5, dpi = 300)
# }


# process_and_save_separate(scenarios$S1, "S1", "仅Var1变化")
# process_and_save_separate(scenarios$S2, "S2", "Var1和Var2均变化")
# process_and_save_separate(scenarios$S3, "S3", "仅Var2变化")




######单个变点
library(SNSeg)
library(ggplot2)
library(patchwork)


n <- 300
cp_true <- 150
set.seed(2026)

s1_v1 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))
s1_v2 <- rnorm(n, 0, 1)

s2_v1 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))
s2_v2 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))

s3_v1 <- rnorm(n, 0, 1)
s3_v2 <- c(rnorm(cp_true, 0, 1), rnorm(n - cp_true, 2, 1))

scenarios <- list(
  S1 = data.frame(Time=1:n, Var1=s1_v1, Var2=s1_v2),
  S2 = data.frame(Time=1:n, Var1=s2_v1, Var2=s2_v2),
  S3 = data.frame(Time=1:n, Var1=s3_v1, Var2=s3_v2)
)


process_and_save_side_by_side <- function(data, filename, title_prefix) {
  
  
  res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "mean")
  cp_uni <- if(length(res_uni$est_cp) > 0) res_uni$est_cp[1] else NULL
  
  res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "mean")
  cp_multi <- if(length(res_multi$est_cp) > 0) res_multi$est_cp[1] else NULL
  
  
  create_sub_plot <- function(df, detected_cp, sub_title, show_var2 = FALSE) {
    p <- ggplot(df, aes(x = Time)) +
      geom_line(aes(y = Var1, color = "Var1"), alpha = 0.7) +
      theme_minimal() +
      labs(subtitle = sub_title, y = "观测值", x = "时间") +
      scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "orange")) +
      
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
    plot_layout(guides = "collect") +
    plot_annotation(title = title_prefix, 
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) &
    theme(legend.position = "bottom") 
  
  ggsave(filename, combined_plot, width = 12, height = 5, dpi = 300)
}


process_and_save_side_by_side(scenarios$S1, "Scenario_1_mean_1.png", "仅Var1发生变化")
process_and_save_side_by_side(scenarios$S2, "Scenario_2_mean_1.png", "Var1和Var2均发生变化")
process_and_save_side_by_side(scenarios$S3, "Scenario_3_mean_1.png", "仅Var2发生变化")



######双变点
library(SNSeg)
library(ggplot2)
library(patchwork)

n <- 300
cp_true <- c(100, 200) 
set.seed(2026)

# 场景1：仅Var1发生两次变点
s1_v1 <- c(rnorm(100, 0, 1), rnorm(100, 2, 1), rnorm(100, -1, 1))
s1_v2 <- rnorm(n, 0, 1)

# 场景2：Var1和Var2同步发生两次变点
s2_v1 <- c(rnorm(100, 0, 1), rnorm(100, 2, 1), rnorm(100, -1, 1))
s2_v2 <- c(rnorm(100, 0, 1), rnorm(100, 2, 1), rnorm(100, -1, 1))

# 场景3：仅Var2发生两次变点
s3_v1 <- rnorm(n, 0, 1)
s3_v2 <- c(rnorm(100, 0, 1), rnorm(100, 2, 1), rnorm(100, -1, 1))

scenarios <- list(
  S1 = data.frame(Time=1:n, Var1=s1_v1, Var2=s1_v2),
  S2 = data.frame(Time=1:n, Var1=s2_v1, Var2=s2_v2),
  S3 = data.frame(Time=1:n, Var1=s3_v1, Var2=s3_v2)
)

process_and_save_side_by_side <- function(data, filename, title_prefix) {
  
  res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "mean")
  cp_uni <- if(length(res_uni$est_cp) > 0) res_uni$est_cp else NULL
  
  res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "mean")
  cp_multi <- if(length(res_multi$est_cp) > 0) res_multi$est_cp else NULL
  
  create_sub_plot <- function(df, detected_cp, sub_title, show_var2 = FALSE) {
    p <- ggplot(df, aes(x = Time)) +
      geom_line(aes(y = Var1, color = "Var1"), alpha = 0.7) +
      theme_minimal() +
      labs(subtitle = sub_title, y = "观测值", x = "时间") +
      scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "orange")) +
      
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
    plot_layout(guides = "collect") + 
    plot_annotation(title = title_prefix, 
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) &
    theme(legend.position = "bottom")
  
  ggsave(filename, combined_plot, width = 12, height = 5, dpi = 300)
}

process_and_save_side_by_side(scenarios$S1, "Scenario_1_mean_2.png", "仅Var1发生变化")
process_and_save_side_by_side(scenarios$S2, "Scenario_2_mean_2.png", "Var1和Var2均发生变化")
process_and_save_side_by_side(scenarios$S3, "Scenario_3_mean_2.png", "仅Var2发生变化")
