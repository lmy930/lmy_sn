###### 均值变点示意图：单维 / 二维；单变点 / 双变点
library(SNSeg)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

set.seed(2026)

# ===============================
# 1. 数据生成函数
# ===============================

# ---- 单维：一个序列，均值变点 ----
generate_uni_data <- function(n = 600, cps = c(300), shift = 3.5,
                              type = c("single", "double")) {
  type <- match.arg(type)
  
  x <- as.numeric(arima.sim(model = list(ar = 0.5), n = n)) +
    rnorm(n, mean = 0, sd = 0.5)
  
  if (type == "single") {
    x[cps[1]:n] <- x[cps[1]:n] + shift
  }
  
  if (type == "double") {
    x[cps[1]:n] <- x[cps[1]:n] + shift
    x[cps[2]:n] <- x[cps[2]:n] - shift * 1.2
  }
  
  data.frame(
    Time = 1:n,
    Value = x
  )
}


# ---- 二维：两个序列，均值变点 ----
# 设置不同基线，保证两个均值不重合
generate_bi_data <- function(n = 600, cps = c(300), shift = c(3.5, 3.0),
                             type = c("single", "double")) {
  type <- match.arg(type)
  
  common_trend <- as.numeric(arima.sim(model = list(ar = 0.5), n = n))
  
  base1 <- 2.5
  base2 <- -2.5
  
  x1 <- base1 + common_trend + rnorm(n, mean = 0, sd = 0.5)
  x2 <- base2 + common_trend + rnorm(n, mean = 0, sd = 0.5)
  
  if (type == "single") {
    x1[cps[1]:n] <- x1[cps[1]:n] + shift[1]
    x2[cps[1]:n] <- x2[cps[1]:n] + shift[2]
  }
  
  if (type == "double") {
    x1[cps[1]:n] <- x1[cps[1]:n] + shift[1]
    x2[cps[1]:n] <- x2[cps[1]:n] + shift[2]
    
    x1[cps[2]:n] <- x1[cps[2]:n] - shift[1] * 1.2
    x2[cps[2]:n] <- x2[cps[2]:n] - shift[2] * 1.2
  }
  
  data.frame(
    Time = 1:n,
    `变量1` = x1,
    `变量2` = x2,
    check.names = FALSE
  )
}


# ===============================
# 2. SNSeg 检测函数
# ===============================

detect_uni_cp <- function(df) {
  fit <- SNSeg_Uni(
    ts = df$Value,
    paras_to_test = "mean",
    confidence = 0.95,
    grid_size_scale = 0.05,
    plot_SN = FALSE,
    est_cp_loc = TRUE
  )
  as.numeric(fit$est_cp)
}

detect_bi_cp <- function(df) {
  x_mat <- as.matrix(df[, c("变量1", "变量2")])
  
  fit <- SNSeg_Multi(
    ts = x_mat,
    paras_to_test = "mean",
    confidence = 0.95,
    grid_size_scale = 0.05,
    plot_SN = FALSE,
    est_cp_loc = TRUE
  )
  as.numeric(fit$est_cp)
}


# ===============================
# 3. 绘图函数
# ===============================

# ---- 单维绘图：蓝色 ----
plot_uni_scenario <- function(df, title, cp_locs) {
  p <- ggplot(df, aes(x = Time, y = Value)) +
    geom_line(color = "#e31a1c", linetype = "dashed",linewidth = 1.1,alpha = 0.9) +
    labs(title = title, x = "时间", y = "观测值") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray85", linetype = "solid")
    )
  
  for (loc in cp_locs) {
    p <- p + geom_vline(
      xintercept = loc,
      color = "#33a02c",
      linetype = "dotted",
      linewidth = 1.1
    )
  }
  
  p
}


# ---- 二维绘图：蓝色 + 绿色 ----
plot_bi_scenario <- function(df, title, cp_locs) {
  df_long <- df %>%
    pivot_longer(cols = -Time, names_to = "Variable", values_to = "Value")
  
  df_long$Variable <- factor(df_long$Variable, levels = c("变量1", "变量2"))
  
  custom_colors <- c("变量1" = "#e31a1c",
                     "变量2" = "#1f78b4")
  
  custom_linetypes <- c("变量1" = "dashed",
                        "变量2" = "dashed")
  
  custom_linewidths <- c("变量1" = 1.1,
                         "变量2" = 1.0)
  
  p <- ggplot(df_long, aes(x = Time, y = Value,
                           color = Variable, linetype = Variable)) +
    geom_line(aes(linewidth = Variable), alpha = 0.9) +
    scale_color_manual(values = custom_colors) +
    scale_linetype_manual(values = custom_linetypes) +
    scale_linewidth_manual(values = custom_linewidths) +
    labs(title = title, x = "时间", y = "观测值") +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.title = element_blank(),
      legend.key.width = unit(2.5, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray85", linetype = "solid")
    )
  
  for (loc in cp_locs) {
    p <- p + geom_vline(
      xintercept = loc,
      color = "#33a02c",
      linetype = "dashed",
      linewidth = 1.1
    )
  }
  
  p
}


# ===============================
# 4. 生成数据
# ===============================

n_len <- 600
cp1 <- 200
cp2 <- 400

# 单维
df_uni_single <- generate_uni_data(
  n = n_len,
  cps = c(cp1),
  type = "single"
)

df_uni_double <- generate_uni_data(
  n = n_len,
  cps = c(cp1, cp2),
  type = "double"
)

# 二维
df_bi_single <- generate_bi_data(
  n = n_len,
  cps = c(cp1),
  type = "single"
)

df_bi_double <- generate_bi_data(
  n = n_len,
  cps = c(cp1, cp2),
  type = "double"
)


# ===============================
# 5. SNSeg 检测
# ===============================

est_uni_single <- detect_uni_cp(df_uni_single)
est_uni_double <- detect_uni_cp(df_uni_double)

est_bi_single <- detect_bi_cp(df_bi_single)
est_bi_double <- detect_bi_cp(df_bi_double)

cat("单维单变点估计结果：", est_uni_single, "\n")
cat("单维双变点估计结果：", est_uni_double, "\n")
cat("二维单变点估计结果：", est_bi_single, "\n")
cat("二维双变点估计结果：", est_bi_double, "\n")


# ===============================
# 6. 绘图
# ===============================

p1 <- plot_uni_scenario(
  df_uni_single,
  title = "单维均值单变点",
  cp_locs = c(cp1)
)

p2 <- plot_uni_scenario(
  df_uni_double,
  title = "单维均值双变点",
  cp_locs = c(cp1, cp2)
)

p3 <- plot_bi_scenario(
  df_bi_single,
  title = "二维均值单变点",
  cp_locs = c(cp1)
)

p4 <- plot_bi_scenario(
  df_bi_double,
  title = "二维均值双变点",
  cp_locs = c(cp1, cp2)
)

plot_uni <- p1 | p2

plot_bi <- (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


# ===============================
# 7. 展示与保存
# ===============================

print(plot_uni)
print(plot_bi)

ggsave(
  "Mean_ChangePoint_Univariate.png",
  plot = plot_uni,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "Mean_ChangePoint_Bivariate.png",
  plot = plot_bi,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)