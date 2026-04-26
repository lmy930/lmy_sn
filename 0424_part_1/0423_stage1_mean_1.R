library(SNSeg)
library(ggplot2)
library(patchwork)
library(mclust)  

n <- 500              
shift_v1 <- 2.5       
shift_v2 <- 4.0       
margin_tol <- 15      
num_sim <- 50        

set.seed(2026)

generate_ar1_data <- function(n, cp = NULL, mean_before = 0, mean_after = 3, ar_coef = 0.3, sd = 0.5) {
  series <- as.numeric(arima.sim(n = n, model = list(ar = ar_coef), sd = sd))
  if (!is.null(cp)) {
    series[1:cp] <- series[1:cp] + mean_before
    series[(cp + 1):n] <- series[(cp + 1):n] + mean_after
  } else {
    series <- series + mean_before
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
 #计算 F1-score
  if (length(est_cps) == 0) {
    precision <- 0
  } else {
    tp_p <- sum(sapply(est_cps, function(e) any(abs(true_cps - e) <= margin)))
    precision <- tp_p / length(est_cps)
  }
  
  if (length(true_cps) == 0) {
    recall <- ifelse(length(est_cps) == 0, 1, 0) 
  } else {
    tp_r <- sum(sapply(true_cps, function(t) any(abs(est_cps - t) <= margin)))
    recall <- tp_r / length(true_cps)
  }
  
  f1 <- ifelse((precision + recall) == 0, 0, 2 * precision * recall / (precision + recall))
  
  #计算 ARI
  true_labels <- get_segment_labels(true_cps, n)
  est_labels <- get_segment_labels(est_cps, n)
  ari <- adjustedRandIndex(true_labels, est_labels)
  
  # 计算 Hausdorff Distance
  if (length(true_cps) == 0 && length(est_cps) == 0) {
    hd <- 0
  } else if (length(est_cps) == 0 || length(true_cps) == 0) {
    hd <- n  
  } else {
    dist_T_to_E <- max(sapply(true_cps, function(t) min(abs(est_cps - t))))
    dist_E_to_T <- max(sapply(est_cps, function(e) min(abs(true_cps - e))))
    hd <- max(dist_T_to_E, dist_E_to_T)
  }
  
  return(c(F1 = f1, ARI = ari, HD = hd))
}


global_true_cps <- list(
  S1 = 250,           # 仅 Var1 变
  S2 = 250,           # Var1, Var2 同时变
  S3 = 250,           # 仅 Var2 变
  S4 = c(150, 350)    # 均变，但时刻不同
)

results_list <- list()
all_sim_data <- vector("list", num_sim)

for (sim in 1:num_sim) {
  if(sim %% 10 == 0) #cat(sprintf("  进行到第 %d 次模拟...\n", sim))
  
  # 生成当次模拟的四种场景数据
  S1 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, 250, mean_after=shift_v1), Var2 = generate_ar1_data(n, NULL))
  S2 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, 250, mean_after=shift_v1), Var2 = generate_ar1_data(n, 250, mean_after=shift_v2))
  S3 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, NULL),                     Var2 = generate_ar1_data(n, 250, mean_after=shift_v2))
  S4 <- data.frame(Time=1:n, Var1 = generate_ar1_data(n, 150, mean_after=shift_v1), Var2 = generate_ar1_data(n, 350, mean_after=shift_v2))
  
  scenarios_data <- list(S1=S1, S2=S2, S3=S3, S4=S4)
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
#cat("\n=== 模拟结果均值统计表 ===\n")
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
      p <- p + geom_vline(xintercept = detected_cp, color = "#20b78e", linetype = "solid", linewidth = 1.0, alpha = 0.9)
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
                    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))) &
    theme(legend.position = "bottom") 
  
  ggsave(filename, combined_plot, width = 12, height = 5, dpi = 300)
}



for (scene in names(global_true_cps)) {
 
  sub_res <- results[results$Scenario == scene & results$Method == "Multivariate", ]
  
  ordered_res <- sub_res[order(-sub_res$F1, sub_res$HD, -sub_res$ARI), ]
  
  best_sim_id <- ordered_res$Simulation[1]
  
  best_data <- all_sim_data[[best_sim_id]][[scene]]
  
  titles <- c(S1 = "仅 Var1 发生变化", S2 = "Var1 和 Var2 发生变化，且时刻相同",
              S3 = "仅 Var2 发生变化", S4 = "Var1 和 Var2 发生变化，且时刻不同")
  filename <- sprintf("Scenario_%s_best.png", scene)
  
  process_and_save_side_by_side(best_data, filename, titles[scene], global_true_cps[[scene]])
}
