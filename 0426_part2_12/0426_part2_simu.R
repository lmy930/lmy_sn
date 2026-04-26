######介绍正常/异常的变点
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

generate_demo_data <- function(n = 600, cps = c(300), shift = 4.5, condition = "Single_Normal") {
  set.seed(2026) 
  
  common_trend <- as.numeric(arima.sim(model = list(ar = 0.6), n = n))
  
  base_target <- 0
  base_nb1    <- 8   
  base_nb2    <- -8  
  
  target <- base_target + common_trend + rnorm(n, mean = 0, sd = 0.6)
  nb1    <- base_nb1    + common_trend + rnorm(n, mean = 0, sd = 0.6)
  nb2    <- base_nb2    + common_trend + rnorm(n, mean = 0, sd = 0.6)
  
  if (condition == "Single_Normal") {
    target[cps[1]:n] <- target[cps[1]:n] + shift
    nb1[cps[1]:n]    <- nb1[cps[1]:n] + shift
    nb2[cps[1]:n]    <- nb2[cps[1]:n] + shift
    
  } else if (condition == "Single_Anomaly") {
    target[cps[1]:n] <- target[cps[1]:n] + shift
    
  } else if (condition == "Double_Normal") {
    target[cps[1]:n] <- target[cps[1]:n] + shift
    nb1[cps[1]:n]    <- nb1[cps[1]:n] + shift
    nb2[cps[1]:n]    <- nb2[cps[1]:n] + shift
    
    target[cps[2]:n] <- target[cps[2]:n] - shift * 1.2
    nb1[cps[2]:n]    <- nb1[cps[2]:n] - shift * 1.2
    nb2[cps[2]:n]    <- nb2[cps[2]:n] - shift * 1.2
    
  } else if (condition == "Double_Anomaly") {
    target[cps[1]:n] <- target[cps[1]:n] + shift
    target[cps[2]:n] <- target[cps[2]:n] - shift * 1.2
  }
  
  data.frame(
    Time = 1:n,
    `站点1` = target,
    `站点2` = nb1,
    `站点3` = nb2,
    check.names = FALSE
  )
}


n_len <- 600
cp1 <- 200
cp2 <- 400

df_single_norm <- generate_demo_data(n = n_len, cps = c(cp1), condition = "Single_Normal")
df_single_anom <- generate_demo_data(n = n_len, cps = c(cp1), condition = "Single_Anomaly")
df_double_norm <- generate_demo_data(n = n_len, cps = c(cp1, cp2), condition = "Double_Normal")
df_double_anom <- generate_demo_data(n = n_len, cps = c(cp1, cp2), condition = "Double_Anomaly")


plot_scenario <- function(df, title, cp_locs, cp_color) {
  df_long <- df %>% pivot_longer(cols = -Time, names_to = "Station", values_to = "Value")
  df_long$Station <- factor(df_long$Station, levels = c("站点1", "站点2", "站点3"))
  
  custom_colors <- c("站点1" = "#e31a1c",      
                     "站点2" = "#1f78b4", 
                     "站点3" = "#33a02c") 
  custom_linetypes <- c("站点1" = "solid", 
                        "站点2" = "dashed", 
                        "站点3" = "dashed")
  custom_linewidths <- c("站点1" = 1.2, 
                         "站点2" = 0.8, 
                         "站点3" = 0.8)
  
  p <- ggplot(df_long, aes(x = Time, y = Value, color = Station, linetype = Station)) +
    geom_line(aes(linewidth = Station), alpha = 0.85) +
    scale_color_manual(values = custom_colors) +
    scale_linetype_manual(values = custom_linetypes) +
    scale_linewidth_manual(values = custom_linewidths) +
    scale_y_continuous(limits = c(-12, 16), breaks = seq(-10, 15, by = 5)) +
    labs(title = title, x = "时间", y = "观测值") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_blank(),
      legend.key.width = unit(2.5, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray85", linetype = "solid")
    )
  
  for(loc in cp_locs) {
    p <- p + geom_vline(xintercept = loc, color = cp_color, linetype = "dotted", linewidth = 1.2)
  }
  
  return(p)
}

p1 <- plot_scenario(df_single_norm, "正常均值单变点", cp_locs = c(cp1), cp_color = "#33a02c")
p2 <- plot_scenario(df_single_anom, "异常均值单变点", cp_locs = c(cp1), cp_color = "#e31a1c")

p3 <- plot_scenario(df_double_norm, "正常均值双变点", cp_locs = c(cp1, cp2), cp_color = "#33a02c")
p4 <- plot_scenario(df_double_anom, "异常均值双变点", cp_locs = c(cp1, cp2), cp_color = "#e31a1c")


plot_single <- (p1 | p2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

plot_double <- (p3 | p4) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# print(plot_single)
# print(plot_double)

 ggsave("Single_ChangePoint_Comparison.png", plot = plot_single, width = 14, height = 6, dpi = 300, bg = "white")
 ggsave("Double_ChangePoint_Comparison.png", plot = plot_double, width = 14, height = 6, dpi = 300, bg = "white")