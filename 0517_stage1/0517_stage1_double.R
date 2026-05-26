# ============================================================
# 双变点情形：单变量 AR(1) vs 双变量 VAR(1) 均值变点检测
#
# 情况1：
#   第1个变点：仅变量2发生突变（对应单变点情形1）
#   第2个变点：两个变量均发生微小变化（对应单变点情形2）
#   单变量方法：变量1
#   双变量方法：二维联合序列
#
# 情况2：
#   第1个变点：仅变量1发生突变（对应单变点情形3）
#   第2个变点：变量1和变量2同步发生突变（对应单变点情形4）
#   单变量方法：变量1
#   双变量方法：二维联合序列
#
# 绘图：
#   每种情况均保存一张图：左图为单变量1结果，右图为双变量结果
#
# 评价指标：
#   1) 检测到变点的次数
#   2) 检测到变点的比例
#   3) 平均 H 距离
#   4) 平均 ARI
#
# 说明：
#   双变量检测前对联合序列逐变量标准化：scale(x)
#   H距离和ARI：只要该次模拟检测到至少1个变点，就计算
# ============================================================

library(SNSeg)
library(ggplot2)
library(patchwork)

# ----------------------------
# 1. 基本参数
# ----------------------------
set.seed(3047)

B <- 50                     # 模拟次数；后续正式实验可改大
n <- 600                  # 序列长度
cp_true <- c(200, 400)    # 两个真实变点位置
tol <- 2                # 仅用于挑选最优图时判断命中，不影响检测本身

grid_size_scale <- 0.05
confidence <- 0.95

# 如需固定 grid_size，可取消下一行注释，
# 并在 detect_uni / detect_multi 中加入 grid_size = grid_size
# grid_size <- 20


# ----------------------------
# 2. VAR(1) 模型参数
# ----------------------------

# 对角项与单变量 AR(1) 对应，用于控制变量
A <- matrix(c(0.30, 0.10,
              0.12, 0.25), nrow = 2, byrow = TRUE)

phi_ar <- diag(A)

sd_eps <- c(5.0, 5.0)
rho_eps <- 0.65

Sigma <- matrix(c(sd_eps[1]^2, rho_eps * prod(sd_eps),
                  rho_eps * prod(sd_eps), sd_eps[2]^2),
                nrow = 2, byrow = TRUE)

mu0 <- c(25, 10)

# ----------------------------
# 3. 两种双变点情况的均值突变设置
# ----------------------------

# 情况1：
# 第1个变点：仅变量2突变（对应原情形1）
# 第2个变点：两个变量均发生微小变化（对应原情形2）
#
# 这里采用“累积变化”的方式：
# 第1段：mu0
# 第2段：mu0 + delta1_case1
# 第3段：mu0 + delta1_case1 + delta2_case1
delta1_case1 <- c(0, -7.5)
delta2_case1 <- c(-2.5, 2.0)

# 情况2：
# 第1个变点：仅变量1突变（对应原情形3）
# 第2个变点：变量1和变量2同步突变（对应原情形4）
delta1_case2 <- c(-7.5, 0)
delta2_case2 <- c(7.5, -6.5)

# 配色
cols <- c("变量1" = "#B2182B", "变量2" = "#2166AC")


# ----------------------------
# 4. 数据生成函数：双变点 VAR(1)
# X_t = mu_t + Y_t
# Y_t = A Y_{t-1} + eps_t
# ----------------------------

rmvnorm2 <- function(Sigma) {
  as.numeric(t(chol(Sigma)) %*% rnorm(2))
}

simulate_var1_2cp <- function(n, cp, mu0, delta1, delta2, A, Sigma, burn = 100) {
  cp <- sort(cp)
  cp1 <- cp[1]
  cp2 <- cp[2]

  total_n <- n + burn
  y <- matrix(0, total_n, 2)

  for (t in 2:total_n) {
    y[t, ] <- A %*% y[t - 1, ] + rmvnorm2(Sigma)
  }

  y <- y[(burn + 1):total_n, ]

  mu_mat <- matrix(rep(mu0, each = n), ncol = 2)

  # 第2段
  mu_mat[(cp1 + 1):n, ] <- sweep(mu_mat[(cp1 + 1):n, ], 2, delta1, "+")

  # 第3段（累积变化）
  mu_mat[(cp2 + 1):n, ] <- sweep(mu_mat[(cp2 + 1):n, ], 2, delta2, "+")

  x <- mu_mat + y
  colnames(x) <- c("变量1", "变量2")

  x
}


# ----------------------------
# 5. SNSeg 检测函数
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
# 6. 评价指标函数
# ----------------------------

# 是否检测到至少1个变点
detect_indicator <- function(cp_hat) {
  as.integer(length(cp_hat) > 0)
}

# H距离（Hausdorff 距离）
# 针对“估计变点集合”与“真实变点集合”计算
hausdorff_cp <- function(est_cp, true_cp) {
  if (length(est_cp) == 0 || length(true_cp) == 0) return(NA_real_)

  est_cp <- sort(unique(est_cp))
  true_cp <- sort(unique(true_cp))

  d1 <- max(sapply(est_cp, function(x) min(abs(x - true_cp))))
  d2 <- max(sapply(true_cp, function(x) min(abs(x - est_cp))))

  max(d1, d2)
}

# 根据变点位置生成分段标签
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

# 计算 ARI
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

# 根据估计变点与真实变点计算 ARI
ari_cp <- function(est_cp, true_cp, n) {
  if (length(est_cp) == 0) return(NA_real_)

  label_true <- cp_to_label(true_cp, n)
  label_est <- cp_to_label(est_cp, n)

  ari_score(label_true, label_est)
}

# 统计估计变点中命中了多少个真实变点
count_hit_cp <- function(est_cp, true_cp, tol) {
  if (length(est_cp) == 0) return(0L)

  sum(sapply(true_cp, function(c0) any(abs(est_cp - c0) <= tol)))
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# ----------------------------
# 7. 挑选最优图的评分函数
# ----------------------------

# 情况1理想结果：
# 单变量1：尽量不检出（因为第1个变点不在变量1上，第2个变点是微小变化）
# 双变量：尽量检出两个变点
score_case1 <- function(cp_uni1, cp_multi) {
  uni1_no <- length(cp_uni1) == 0
  multi_hit_n <- count_hit_cp(cp_multi, cp_true, tol)

  h_multi <- hausdorff_cp(cp_multi, cp_true)
  if (is.na(h_multi)) h_multi <- 999

  ari_multi <- ari_cp(cp_multi, cp_true, n)
  if (is.na(ari_multi)) ari_multi <- 0

  1500 * uni1_no +
    1200 * multi_hit_n +
    100 * ari_multi -
    h_multi -
    25 * max(0, length(cp_multi) - length(cp_true)) -
    10 * length(cp_uni1)
}

# 情况2理想结果：
# 单变量1：尽量检出两个变点
# 双变量：尽量检出两个变点
score_case2 <- function(cp_uni1, cp_multi) {
  uni1_hit_n <- count_hit_cp(cp_uni1, cp_true, tol)
  multi_hit_n <- count_hit_cp(cp_multi, cp_true, tol)

  h_uni1 <- hausdorff_cp(cp_uni1, cp_true)
  h_multi <- hausdorff_cp(cp_multi, cp_true)
  if (is.na(h_uni1)) h_uni1 <- 999
  if (is.na(h_multi)) h_multi <- 999

  ari_uni1 <- ari_cp(cp_uni1, cp_true, n)
  ari_multi <- ari_cp(cp_multi, cp_true, n)
  if (is.na(ari_uni1)) ari_uni1 <- 0
  if (is.na(ari_multi)) ari_multi <- 0

  1200 * uni1_hit_n +
    1200 * multi_hit_n +
    100 * ari_uni1 +
    100 * ari_multi -
    h_uni1 -
    h_multi -
    25 * max(0, length(cp_uni1) - length(cp_true)) -
    25 * max(0, length(cp_multi) - length(cp_true))
}


# ----------------------------
# 8. 绘图函数
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
        alpha = 0.9
      )
  }

  p
}


# ============================================================
# 9. 双变点情况1
# ============================================================

case1_list <- vector("list", B)

for (b in 1:B) {
  set.seed(3026 + b)

  x <- simulate_var1_2cp(
    n = n,
    cp = cp_true,
    mu0 = mu0,
    delta1 = delta1_case1,
    delta2 = delta2_case1,
    A = A,
    Sigma = Sigma
  )

  # 单变量方法：变量1
  cp_uni1 <- detect_uni(x[, "变量1"])

  # 双变量方法：二维联合序列（标准化后检测）
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
    title = "变量2有均值变化、变量1和变量2同步地有微小均值变化",
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
  filename = "双变点情况1_最佳结果.png",
  plot = fig_case1,
  width = 12,
  height = 4.8,
  dpi = 320
)


# ============================================================
# 10. 双变点情况2
# ============================================================

case2_list <- vector("list", B)

for (b in 1:B) {
  set.seed(4026 + b)

  x <- simulate_var1_2cp(
    n = n,
    cp = cp_true,
    mu0 = mu0,
    delta1 = delta1_case2,
    delta2 = delta2_case2,
    A = A,
    Sigma = Sigma
  )

  # 单变量方法：变量1
  cp_uni1 <- detect_uni(x[, "变量1"])

  # 双变量方法：二维联合序列（标准化后检测）
  x_std <- scale(x)
  cp_multi <- detect_multi(x_std)

  case2_list[[b]] <- list(
    x = x,
    cp_uni1 = cp_uni1,
    cp_multi = cp_multi,

    det_uni1 = detect_indicator(cp_uni1),
    det_multi = detect_indicator(cp_multi),

    h_uni1 = hausdorff_cp(cp_uni1, cp_true),
    h_multi = hausdorff_cp(cp_multi, cp_true),

    ari_uni1 = ari_cp(cp_uni1, cp_true, n),
    ari_multi = ari_cp(cp_multi, cp_true, n),

    score = score_case2(cp_uni1, cp_multi)
  )
}

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
  vars = c("变量1", "变量2"),
  cp_hat = best2$cp_multi,
  title_text = NULL
)

fig_case2 <- p2_left + p2_right + plot_layout(ncol = 2)+
plot_annotation(
    title = "变量1有均值变化、变量1和变量2同步地有均值变化",
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
  filename = "双变点情况2_最佳结果.png",
  plot = fig_case2,
  width = 12,
  height = 4.8,
  dpi = 320
)


# ============================================================
# 11. 指标表格整理与打印
# ============================================================

result_case1 <- data.frame(
  实验情形 = c("双变点情况1：第1个变点仅变量2突变，第2个变点两个变量均微小变化",
           "双变点情况1：第1个变点仅变量2突变，第2个变点两个变量均微小变化"),
  方法 = c("单变量方法（变量1）", "双变量方法（变量1+变量2）"),
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
  实验情形 = c("双变点情况2：第1个变点仅变量1突变，第2个变点两个变量同步突变",
           "双变点情况2：第1个变点仅变量1突变，第2个变点两个变量同步突变"),
  方法 = c("单变量方法（变量1）", "双变量方法（变量1+变量2）"),
  检测到变点的次数 = c(
    sum(sapply(case2_list, function(z) z$det_uni1)),
    sum(sapply(case2_list, function(z) z$det_multi))
  ),
  模拟次数 = c(B, B),
  检测到变点的比例 = c(
    mean(sapply(case2_list, function(z) z$det_uni1)),
    mean(sapply(case2_list, function(z) z$det_multi))
  ),
  平均H距离 = c(
    safe_mean(sapply(case2_list, function(z) z$h_uni1)),
    safe_mean(sapply(case2_list, function(z) z$h_multi))
  ),
  平均ARI = c(
    safe_mean(sapply(case2_list, function(z) z$ari_uni1)),
    safe_mean(sapply(case2_list, function(z) z$ari_multi))
  )
)

result_table <- rbind(result_case1, result_case2)

result_table$检测到变点的比例 <- round(result_table$检测到变点的比例, 4)
result_table$平均H距离 <- round(result_table$平均H距离, 4)
result_table$平均ARI <- round(result_table$平均ARI, 4)

cat("\n==================== 双变点数值模拟评价指标汇总表 ====================\n")
print(result_table, row.names = FALSE)

# 如需保存指标表格，可取消注释
# write.csv(result_table, "双变点数值模拟评价指标汇总表.csv",
#           row.names = FALSE, fileEncoding = "UTF-8")

