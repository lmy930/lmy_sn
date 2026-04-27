library(SNSeg)
library(ggplot2)
library(patchwork)
library(mclust)  

n <- 800            
margin_tol <- 10      
num_sim <- 5         


means_v1_0 <- c(0)
means_v1_1 <- c(0, 2.5)         # Var1 有1个变点时：0 -> 2.5
means_v1_2 <- c(0, 2.5, 5.0)    # Var1 有2个变点时：0 -> 2.5 -> 5.0

means_v2_0 <- c(0)
means_v2_1 <- c(0, 4.0)         # Var2 有1个变点时：0 -> 4.0
means_v2_2 <- c(0, 4.0, -1.0)   # Var2 有2个变点时：0 -> 4.0 -> -1.0

set.seed(2026)


generate_ar1_data <- function(n, cps = NULL, means = c(0), ar_coef = 0.3, sd = 0.5) {
  series <- as.numeric(arima.sim(n = n, model = list(ar = ar_coef), sd = sd))
  
  if (is.null(cps) || length(cps) == 0) {
    series <- series + means[1]
  } else {
    cps <- sort(cps)
    start_idx <- 1
    for (i in seq_along(cps)) {
      series[start_idx:cps[i]] <- series[start_idx:cps[i]] + means[i]
      start_idx <- cps[i] + 1
    }
    series[start_idx:n] <- series[start_idx:n] + means[length(cps) + 1]
  }
  return(series)
}


get_segment_labels <- function(cps, n) {
  labels <- rep(1, n)
  if (length(cps) > 0) {
    cps <- sort(cps)
    for (i in seq_along(cps)) {
      if (cps[i] >= 1 && cps[i] < n) {
        labels[(cps[i] + 1):n] <- i + 1
      }
    }
  }
  return(labels)
}


calc_metrics <- function(true_cps, est_cps, n, margin = 15) {
  # F1-score
  if (length(est_cps) == 0) { precision <- 0 } else {
    tp_p <- sum(sapply(est_cps, function(e) any(abs(true_cps - e) <= margin)))
    precision <- tp_p / length(est_cps)
  }
  if (length(true_cps) == 0) { recall <- ifelse(length(est_cps) == 0, 1, 0) } else {
    tp_r <- sum(sapply(true_cps, function(t) any(abs(est_cps - t) <= margin)))
    recall <- tp_r / length(true_cps)
  }
  f1 <- ifelse((precision + recall) == 0, 0, 2 * precision * recall / (precision + recall))
  
  # ARI
  true_labels <- get_segment_labels(true_cps, n)
  est_labels <- get_segment_labels(est_cps, n)
  ari <- adjustedRandIndex(true_labels, est_labels)
  
  # Hausdorff Distance
  if (length(true_cps) == 0 && length(est_cps) == 0) { hd <- 0 } 
  else if (length(est_cps) == 0 || length(true_cps) == 0) { hd <- n } 
  else {
    dist_T_to_E <- max(sapply(true_cps, function(t) min(abs(est_cps - t))))
    dist_E_to_T <- max(sapply(est_cps, function(e) min(abs(true_cps - e))))
    hd <- max(dist_T_to_E, dist_E_to_T)
  }
  return(c(F1 = f1, ARI = ari, HD = hd))
}


global_true_cps <- list(
  S1 = c(250, 550),              # Var1 两点
  S2 = c(250, 550),              # Var2 两点
  S3 = c(250, 550),              # 均两点，时刻完全相同
  S4 = c(200, 300, 500, 600),    # 均两点，时刻完全不同
  S5 = c(250, 500, 600),         # 均两点，仅 250 相同
  
  # 新增的 4 种情况
  S6 = c(250, 400, 600),         # Var1 一点(400), Var2 两点(250, 600) -> 全不同
  S7 = c(250, 600),              # Var1 一点(250), Var2 两点(250, 600) -> Var1与Var2之一相同
  S8 = c(250, 600),              # Var1 两点(250, 600), Var2 一点(250) -> Var2与Var1之一相同
  S9 = c(250, 400, 600)          # Var1 两点(250, 600), Var2 一点(400) -> 全不同
)

results_list <- list()
all_sim_data <- vector("list", num_sim)

for (sim in 1:num_sim) {
     if(sim %% 5 == 0) cat(sprintf("  进行到第 %d 次模拟...\n", sim))
  
  # 生成当次模拟的 9 种场景数据
  S1 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, c(250, 550), means_v1_2), Var2 = generate_ar1_data(n, NULL, means_v2_0))
  S2 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, NULL, means_v1_0),       Var2 = generate_ar1_data(n, c(250, 550), means_v2_2))
  S3 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, c(250, 550), means_v1_2), Var2 = generate_ar1_data(n, c(250, 550), means_v2_2))
  S4 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, c(200, 600), means_v1_2), Var2 = generate_ar1_data(n, c(300, 500), means_v2_2))
  S5 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, c(250, 600), means_v1_2), Var2 = generate_ar1_data(n, c(250, 500), means_v2_2))
  
  S6 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, 400, means_v1_1),         Var2 = generate_ar1_data(n, c(250, 600), means_v2_2))
  S7 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, 250, means_v1_1),         Var2 = generate_ar1_data(n, c(250, 600), means_v2_2))
  S8 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, c(250, 600), means_v1_2), Var2 = generate_ar1_data(n, 250, means_v2_1))
  S9 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, c(250, 600), means_v1_2), Var2 = generate_ar1_data(n, 400, means_v2_1))
  
  scenarios_data <- list(S1=S1, S2=S2, S3=S3, S4=S4, S5=S5, S6=S6, S7=S7, S8=S8, S9=S9)
  all_sim_data[[sim]] <- scenarios_data 
  
  for (scene in names(scenarios_data)) {
    data <- scenarios_data[[scene]]
    t_cps <- global_true_cps[[scene]]
    
   
    res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "mean")
    res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "mean")
    
    est_uni <- if(length(res_uni$est_cp) > 0) res_uni$est_cp else numeric(0)
    est_multi <- if(length(res_multi$est_cp) > 0) res_multi$est_cp else numeric(0)
    
    metric_uni <- calc_metrics(t_cps, est_uni, n, margin_tol)
    metric_multi <- calc_metrics(t_cps, est_multi, n, margin_tol)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Simulation = sim, Scenario = scene, Method = "Univariate (Var1)",
      F1 = metric_uni["F1"], ARI = metric_uni["ARI"], HD = metric_uni["HD"]
    )
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Simulation = sim, Scenario = scene, Method = "Multivariate",
      F1 = metric_multi["F1"], ARI = metric_multi["ARI"], HD = metric_multi["HD"]
    )
  }
}

results <- do.call(rbind, results_list)
rownames(results) <- NULL

mean_results <- aggregate(cbind(F1, ARI, HD) ~ Scenario + Method, data = results, FUN = mean)
mean_results[, 3:5] <- round(mean_results[, 3:5], 3)
#cat("\n=== 混合变点综合模拟结果均值统计表 ===\n")
print(mean_results)

process_and_save_side_by_side <- function(data, filename, title_prefix, true_cps) {
  
  res_uni <- SNSeg_Uni(ts = data$Var1, paras_to_test = "mean")
  cp_uni <- if(length(res_uni$est_cp) > 0) res_uni$est_cp else NULL
  
  res_multi <- SNSeg_Multi(ts = cbind(data$Var1, data$Var2), paras_to_test = "mean")
  cp_multi <- if(length(res_multi$est_cp) > 0) res_multi$est_cp else NULL
  
  create_sub_plot <- function(df, detected_cp, sub_title, show_var2 = FALSE) {
    p <- ggplot(df, aes(x = Time)) +
      geom_line(aes(y = Var1, color = "Var1"), alpha = 0.85, linewidth = 0.6) +
      theme_bw() +
      labs(subtitle = sub_title, y = "观测值", x = "时间") +
      scale_color_manual(values = c("Var1" = "steelblue", "Var2" = "#E69F00")) +
      theme(
        legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank()
      )
    
    if(show_var2) { 
      p <- p + geom_line(aes(y = Var2, color = "Var2"), alpha = 0.85, linewidth = 0.6) 
    }
    
    if(!is.null(detected_cp)) {
      p <- p + geom_vline(xintercept = detected_cp, color = "#009E73", linetype = "solid", linewidth = 1.0, alpha = 0.9)
    }
    
    if(length(true_cps) > 0) {
      p <- p + geom_vline(xintercept = true_cps, color = "grey50", linetype = "dashed", linewidth = 1.2, alpha = 0.7)
    }
    
    return(p)
  }
  
  p_left <- create_sub_plot(data, cp_uni, "仅检测 Var1", show_var2 = FALSE)
  p_right <- create_sub_plot(data, cp_multi, "检测 Var1 和 Var2", show_var2 = TRUE)
  
  combined_plot <- (p_left | p_right) +
    plot_layout(guides = "collect") +
    plot_annotation(title = title_prefix, 
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))) &
    theme(legend.position = "bottom") 
  
  ggsave(filename, combined_plot, width = 12, height = 5, dpi = 300)
}


titles <- c(
  S1 = "仅 Var1 有2变点",
  S2 = "仅 Var2 有2变点",
  S3 = "Var1 和 Var2 均有2变点，且时刻均相同",
  S4 = "Var1 和 Var2 均有2变点，且时刻均不同",
  S5 = "Var1 和 Var2 均有2变点，且1个变点时刻相同",
  S6 = "Var1 1变点, Var2 2变点，且时刻不同",
  S7 = "Var1 1变点, Var2 2变点，且1个变点时刻相同",
  S8 = "Var1 2变点, Var1 1变点，且1个变点时刻相同",
  S9 = "Var1 2变点, Var1 1变点，且时刻不同"
)

for (scene in names(global_true_cps)) {
  sub_res <- results[results$Scenario == scene & results$Method == "Multivariate", ]
  
  # 寻找最佳性能：优先 F1 最高，其次 HD 最小，最后 ARI 最高
  ordered_res <- sub_res[order(-sub_res$F1, sub_res$HD, -sub_res$ARI), ]
  best_sim_id <- ordered_res$Simulation[1]

  best_data <- all_sim_data[[best_sim_id]][[scene]]
  filename <- sprintf("Scenario_MixCP_%s_best.png", scene)
  
  process_and_save_side_by_side(best_data, filename, titles[scene], global_true_cps[[scene]])
}
