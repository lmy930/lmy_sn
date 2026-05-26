# ============================================================
# 单变点情形：单变量 AR(1) vs 双变量 VAR(1) 均值变点检测
# 情形1：仅变量2发生一次均值突变
# 情形2：变量1和变量2均发生微小均值突变
#
# 单变量方法：原始单变量序列 + SNSeg_Uni()
# 双变量方法：逐变量标准化后的双变量序列 + SNSeg_Multi()
# 绘图：仍使用原始数据
# 评价指标：检测比例、H距离、ARI
# ============================================================

library(SNSeg)
library(ggplot2)
library(patchwork)

# ----------------------------
# 1. 基本参数
# ----------------------------
set.seed(3047)

B <- 50                  # 模拟次数；正式实验可改为100或200
n <- 500                 # 序列长度
cp_true <- 250           # 真实变点位置
tol <- 2                 # 仅用于挑选最优图，不影响检测结果
grid_size_scale <- 0.05
confidence <- 0.95
# ----------------------------
# 2. VAR(1) 模型参数
# ----------------------------
# VAR(1) 系数矩阵
# 对角项与单变量 AR(1) 的参数对应，用于控制变量
A <- matrix(c(0.30, 0.10,
              0.12, 0.25), nrow = 2, byrow = TRUE)
phi_ar <- diag(A)
# 创新项协方差矩阵
# 变量1和变量2正相关
sd_eps <- c(5.0, 5.0)
rho_eps <- 0.65
Sigma <- matrix(c(sd_eps[1]^2, rho_eps * prod(sd_eps),
                  rho_eps * prod(sd_eps), sd_eps[2]^2),
                nrow = 2, byrow = TRUE)
# 初始均值
mu0 <- c(25, 10)
# 情形1：仅变量2发生一次下降型均值突变
delta_case1 <- c(0, -7.5)
# 情形2：两个变量均发生微小均值突变
delta_case2 <- c(2.5, -2.0)
# 颜色
cols <- c("变量1" = "#B2182B", "变量2" = "#2166AC")


# ----------------------------
# 3. 数据生成函数：VAR(1)
# X_t = mu_t + Y_t
# Y_t = A Y_{t-1} + eps_t
# ----------------------------

rmvnorm2 <- function(Sigma) {
  as.numeric(t(chol(Sigma)) %*% rnorm(2))
}

simulate_var1 <- function(n, cp, mu0, delta, A, Sigma, burn = 100) {
  total_n <- n + burn
  y <- matrix(0, total_n, 2)
  
  for (t in 2:total_n) {
    y[t, ] <- A %*% y[t - 1, ] + rmvnorm2(Sigma)
  }
  
  y <- y[(burn + 1):total_n, ]
  
  mu_mat <- matrix(rep(mu0, each = n), ncol = 2)
  mu_mat[(cp + 1):n, ] <- sweep(mu_mat[(cp + 1):n, ], 2, delta, "+")
  
  x <- mu_mat + y
  colnames(x) <- c("变量1", "变量2")
  
  x
}


# ----------------------------
# 4. SNSeg 检测函数
# ----------------------------

get_cp <- function(res) {
  cp <- res$est_cp
  cp <- cp[!is.na(cp)]
  as.integer(cp)
}

detect_uni <- function(x) {
  res <- tryCatch(
    SNSeg_Uni(
      x,
      paras_to_test = "mean",
      confidence = confidence,
      grid_size_scale = grid_size_scale,
      plot_SN = FALSE,
      est_cp_loc = TRUE
      # 如需固定 grid_size，可添加：
      # , grid_size = grid_size
    ),
    error = function(e) NULL
  )
  
  if (is.null(res)) integer(0) else get_cp(res)
}

detect_multi <- function(xmat) {
  res <- tryCatch(
    SNSeg_Multi(
      xmat,
      paras_to_test = "mean",
      confidence = confidence,
      grid_size_scale = grid_size_scale,
      plot_SN = FALSE,
      est_cp_loc = TRUE
      # 如需固定 grid_size，可添加：
      # , grid_size = grid_size
    ),
    error = function(e) NULL
  )
  
  if (is.null(res)) integer(0) else get_cp(res)
}


# ----------------------------
# 5. 评价指标函数
# ----------------------------

# 是否检测到变点：检测到记1，否则记0
detect_indicator <- function(cp_hat) {
  as.integer(length(cp_hat) > 0)
}

# H距离，即 Hausdorff 距离
# 如果没有检测到变点，则返回 NA
hausdorff_cp <- function(est_cp, true_cp) {
  if (length(est_cp) == 0 || length(true_cp) == 0) return(NA_real_)
  
  est_cp <- sort(unique(est_cp))
  true_cp <- sort(unique(true_cp))
  
  d1 <- max(sapply(est_cp, function(x) min(abs(x - true_cp))))
  d2 <- max(sapply(true_cp, function(x) min(abs(x - est_cp))))
  
  max(d1, d2)
}

# 根据变点位置生成分段标签
# 例如 n=500, cp=250，则 1:250 为第1段，251:500 为第2段
cp_to_label <- function(cp, n) {
  cp <- sort(unique(cp))
  cp <- cp[cp > 0 & cp < n]
  
  label <- integer(n)
  
  if (length(cp) == 0) {
    label[] <- 1
    return(label)
  }
  
  start <- 1
  
  for (k in seq_along(cp)) {
    label[start:cp[k]] <- k
    start <- cp[k] + 1
  }
  
  label[start:n] <- length(cp) + 1
  
  label
}

# 手动计算 ARI，避免额外依赖 mclust 包
ari_score <- function(label_true, label_est) {
  tab <- table(label_true, label_est)
  
  comb2 <- function(x) x * (x - 1) / 2
  
  n_total <- sum(tab)
  sum_ij <- sum(comb2(tab))
  sum_i <- sum(comb2(rowSums(tab)))
  sum_j <- sum(comb2(colSums(tab)))
  total <- comb2(n_total)
  
  expected <- sum_i * sum_j / total
  max_index <- 0.5 * (sum_i + sum_j)
  
  if (max_index == expected) return(1)
  
  (sum_ij - expected) / (max_index - expected)
}

# 根据估计变点和真实变点计算 ARI
# 如果没有检测到变点，则返回 NA
ari_cp <- function(est_cp, true_cp, n) {
  if (length(est_cp) == 0) return(NA_real_)
  
  label_true <- cp_to_label(true_cp, n)
  label_est <- cp_to_label(est_cp, n)
  
  ari_score(label_true, label_est)
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# 情形1选图评分：
# 希望单变量1不检出，双变量检出，且H距离小、ARI大
score_case1 <- function(cp_uni1, cp_multi) {
  uni1_no <- length(cp_uni1) == 0
  multi_yes <- length(cp_multi) > 0
  
  h <- hausdorff_cp(cp_multi, cp_true)
  if (is.na(h)) h <- 999
  
  ari <- ari_cp(cp_multi, cp_true, n)
  if (is.na(ari)) ari <- 0
  
  1000 * uni1_no +
    1000 * multi_yes +
    100 * ari -
    h -
    20 * max(0, length(cp_multi) - 1)
}

# 情形2选图评分：
# 希望两个单变量均不检出，双变量检出，且H距离小、ARI大
score_case2 <- function(cp_uni1, cp_uni2, cp_multi) {
  uni1_no <- length(cp_uni1) == 0
  uni2_no <- length(cp_uni2) == 0
  multi_yes <- length(cp_multi) > 0
  
  h <- hausdorff_cp(cp_multi, cp_true)
  if (is.na(h)) h <- 999
  
  ari <- ari_cp(cp_multi, cp_true, n)
  if (is.na(ari)) ari <- 0
  
  1000 * uni1_no +
    1000 * uni2_no +
    1000 * multi_yes +
    100 * ari -
    h -
    20 * max(0, length(cp_multi) - 1)
}


# ----------------------------
# 6. 绘图函数
# ----------------------------

make_df <- function(xmat, vars) {
  data.frame(
    时间 = rep(seq_len(nrow(xmat)), times = length(vars)),
    数值 = as.vector(xmat[, vars, drop = FALSE]),
    变量 = rep(vars, each = nrow(xmat))
  )
}

plot_result <- function(xmat, vars, cp_hat, title_text) {
  df <- make_df(xmat, vars)
  
  p <- ggplot(df, aes(x = 时间, y = 数值, color = 变量)) +
    geom_line(linewidth = 0.55, alpha = 0.9) +
    geom_vline(
      xintercept = cp_true,
      linetype = "dashed",
      color = "grey35",
      linewidth = 0.7
    ) +
    scale_color_manual(values = cols[vars]) +
    labs(
      title = title_text,
      x = "时间",
      y = "观测值",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  if (length(cp_hat) > 0) {
    p <- p +
      geom_vline(
        xintercept = cp_hat,
        color = "#0dfe91",
        linewidth = 0.8,
        alpha = 0.85
      )
  }
  
  p
}


# ============================================================
# 7. 情形1：仅变量2发生一次均值突变
# ============================================================

case1_list <- vector("list", B)

for (b in 1:B) {
  set.seed(2026 + b)
  
  x <- simulate_var1(
    n = n,
    cp = cp_true,
    mu0 = mu0,
    delta = delta_case1,
    A = A,
    Sigma = Sigma
  )
  
  # 单变量检测：使用原始变量1
  cp_uni1 <- detect_uni(x[, "变量1"])
  
  # 双变量检测：先逐变量标准化，再进行多变量检测
  x_std <- scale(x)
  cp_multi <- detect_multi(x_std)
  
  case1_list[[b]] <- list(
    x = x,
    cp_uni1 = cp_uni1,
    cp_multi = cp_multi,
    
    det_uni1 = detect_indicator(cp_uni1),
    det_multi = detect_indicator(cp_multi),
    
    h_uni1 = hausdorff_cp(cp_uni1, cp_true),
    h_multi = hausdorff_cp(cp_multi, cp_true),
    
    ari_uni1 = ari_cp(cp_uni1, cp_true, n),
    ari_multi = ari_cp(cp_multi, cp_true, n),
    
    score = score_case1(cp_uni1, cp_multi)
  )
}

# 情形1：选择最符合要求的一次模拟并绘图
best1_id <- which.max(sapply(case1_list, function(z) z$score))
best1 <- case1_list[[best1_id]]

p1_left <- plot_result(
  best1$x,
  vars = "变量1",
  cp_hat = best1$cp_uni1,
  title_text = NULL
)

p1_right <- plot_result(
  best1$x,
  vars = c("变量1", "变量2"),
  cp_hat = best1$cp_multi,
  title_text = NULL
)

fig_case1 <- p1_left + p1_right + plot_layout(ncol = 2)+
plot_annotation(
    title = "仅变量2有均值变化",
   theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16,
        margin = margin(t = 15, b = -5)
      )
    )
  ) &
  theme(
    legend.text = element_text(size = 12),
    panel.border = element_blank()
  )

ggsave(
  filename = "情形1_最佳结果.png",
  plot = fig_case1,
  width = 12,
  height = 4.8,
  dpi = 320
)


# ============================================================
# 8. 情形2：两个变量均发生微小均值突变
# ============================================================

case2_list <- vector("list", B)

for (b in 1:B) {
  set.seed(3026 + b)
  
  x <- simulate_var1(
    n = n,
    cp = cp_true,
    mu0 = mu0,
    delta = delta_case2,
    A = A,
    Sigma = Sigma
  )
  
  # 单变量检测：分别使用原始变量1和原始变量2
  cp_uni1 <- detect_uni(x[, "变量1"])
  cp_uni2 <- detect_uni(x[, "变量2"])
  
  # 双变量检测：先逐变量标准化，再进行多变量检测
  x_std <- scale(x)
  cp_multi <- detect_multi(x_std)
  
  case2_list[[b]] <- list(
    x = x,
    cp_uni1 = cp_uni1,
    cp_uni2 = cp_uni2,
    cp_multi = cp_multi,
    
    det_uni1 = detect_indicator(cp_uni1),
    det_uni2 = detect_indicator(cp_uni2),
    det_multi = detect_indicator(cp_multi),
    
    h_uni1 = hausdorff_cp(cp_uni1, cp_true),
    h_uni2 = hausdorff_cp(cp_uni2, cp_true),
    h_multi = hausdorff_cp(cp_multi, cp_true),
    
    ari_uni1 = ari_cp(cp_uni1, cp_true, n),
    ari_uni2 = ari_cp(cp_uni2, cp_true, n),
    ari_multi = ari_cp(cp_multi, cp_true, n),
    
    score = score_case2(cp_uni1, cp_uni2, cp_multi)
  )
}

# 情形2：选择最符合要求的一次模拟并绘图
best2_id <- which.max(sapply(case2_list, function(z) z$score))
best2 <- case2_list[[best2_id]]

p2_left <- plot_result(
  best2$x,
  vars = "变量1",
  cp_hat = best2$cp_uni1,
  title_text = NULL
)

p2_right <- plot_result(
  best2$x,
  vars = "变量2",
  cp_hat = best2$cp_uni2,
  title_text = NULL
)

fig_case2_uni <- p2_left + p2_right + plot_layout(ncol = 2)+
 plot_annotation(
    title = "变量1和变量2同步地有微小均值变化",
    theme = theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16,
        margin = margin(t = 15, b = -5)
      )
    )
  ) &
  theme(
    legend.text = element_text(size = 12),
    panel.border = element_blank()
  )

p2_multi <- plot_result(
  best2$x,
  vars = c("变量1", "变量2"),
  cp_hat = best2$cp_multi,
  title_text = "变量1和变量2同步地有微小均值变化"
)+
theme(
    panel.border = element_blank()
  )

ggsave(
  filename = "情形2_单变量最佳结果.png",
  plot = fig_case2_uni,
  width = 12,
  height = 4.8,
  dpi = 320
)

ggsave(
  filename = "情形2_双变量最佳结果.png",
  plot = p2_multi,
  width = 7.2,
  height = 4.8,
  dpi = 320
)


# ============================================================
# 9. 指标表格整理与打印
# ============================================================

result_case1 <- data.frame(
  实验情形 = c("情形1：仅变量2发生均值突变",
           "情形1：仅变量2发生均值突变"),
  方法 = c("单变量方法（变量1）",
         "双变量方法（变量1+变量2）"),
  检测到变点的次数 = c(
    sum(sapply(case1_list, function(z) z$det_uni1)),
    sum(sapply(case1_list, function(z) z$det_multi))
  ),
  模拟次数 = c(B, B),
  检测到变点的比例 = c(
    mean(sapply(case1_list, function(z) z$det_uni1)),
    mean(sapply(case1_list, function(z) z$det_multi))
  ),
  平均H距离 = c(
    safe_mean(sapply(case1_list, function(z) z$h_uni1)),
    safe_mean(sapply(case1_list, function(z) z$h_multi))
  ),
  平均ARI = c(
    safe_mean(sapply(case1_list, function(z) z$ari_uni1)),
    safe_mean(sapply(case1_list, function(z) z$ari_multi))
  )
)

result_case2 <- data.frame(
  实验情形 = c("情形2：两个变量均发生微小均值突变",
           "情形2：两个变量均发生微小均值突变",
           "情形2：两个变量均发生微小均值突变"),
  方法 = c("单变量方法（变量1）",
         "单变量方法（变量2）",
         "双变量方法（变量1+变量2）"),
  检测到变点的次数 = c(
    sum(sapply(case2_list, function(z) z$det_uni1)),
    sum(sapply(case2_list, function(z) z$det_uni2)),
    sum(sapply(case2_list, function(z) z$det_multi))
  ),
  模拟次数 = c(B, B, B),
  检测到变点的比例 = c(
    mean(sapply(case2_list, function(z) z$det_uni1)),
    mean(sapply(case2_list, function(z) z$det_uni2)),
    mean(sapply(case2_list, function(z) z$det_multi))
  ),
  平均H距离 = c(
    safe_mean(sapply(case2_list, function(z) z$h_uni1)),
    safe_mean(sapply(case2_list, function(z) z$h_uni2)),
    safe_mean(sapply(case2_list, function(z) z$h_multi))
  ),
  平均ARI = c(
    safe_mean(sapply(case2_list, function(z) z$ari_uni1)),
    safe_mean(sapply(case2_list, function(z) z$ari_uni2)),
    safe_mean(sapply(case2_list, function(z) z$ari_multi))
  )
)

result_table <- rbind(result_case1, result_case2)

result_table$检测到变点的比例 <- round(result_table$检测到变点的比例, 4)
result_table$平均H距离 <- round(result_table$平均H距离, 4)
result_table$平均ARI <- round(result_table$平均ARI, 4)

cat("\n==================== 数值模拟评价指标汇总表 ====================\n")
print(result_table, row.names = FALSE)

# 如需保存指标表格，可取消注释
 write.csv(result_table, "单变点数值模拟评价指标汇总表_情形1和情形2.csv",
           row.names = FALSE, fileEncoding = "UTF-8")
