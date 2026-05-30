# ============================================================
# 实际数据均值变点检测可视化
# 上图：PM2.5 单变量均值变点检测结果
# 下图：PM10 与 PM2.5 联合序列多变量均值变点检测结果
# 说明：
#   1）检测时仍使用标准化数据，与前面筛选代码保持一致；
#   2）绘图时：上图用 PM2.5 原始浓度；下图用 PM10 和 PM2.5 原始浓度。
# ============================================================

library(tidyverse)
library(lubridate)
library(ggplot2)
library(SNSeg)
library(patchwork)

# ----------------------------
# 1. 用户参数设置
# ----------------------------

data_dir <- "D:/03/2025/Beijing"
target_station <- "Dongsi"

# 修改为你筛选好的时间段
start_time <- as.POSIXct("2024-12-07 00:00:00", tz = "Asia/Shanghai")
end_time   <- as.POSIXct("2024-12-16 23:00:00", tz = "Asia/Shanghai")

# SNSeg 参数：与前面筛选代码保持一致
sn_confidence <- 0.95
sn_grid_size_scale <- 0.05

# 是否显示图例：TRUE 显示，FALSE 不显示
show_legend <- FALSE

# 字体与字号设置
font_family <- "SimHei"

base_size          <- 14
title_size         <- 17
subtitle_size      <- 12
axis_title_size    <- 13
axis_text_size     <- 11
legend_title_size  <- 12
legend_text_size   <- 11

# 横坐标时间间隔
x_date_breaks <- "1 day"

# 输出文件夹
out_dir <- file.path(data_dir, "actual_data_cp_visualization")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 站点中文名
station_name_zh_map <- c(
  "Dongsi"        = "东四",
  "Tiantan"       = "天坛",
  "Guanyuan"      = "官园",
  "Wanshouxigong" = "万寿西宫",
  "Nongzhanguan"  = "农展馆",
  "Aotizhongxin"  = "奥体中心"
)

target_station_zh <- ifelse(
  target_station %in% names(station_name_zh_map),
  station_name_zh_map[[target_station]],
  target_station
)

output_image <- file.path(
  out_dir,
  paste0(
    target_station, "_PM25_uni_PM10_PM25_multi_rawplot_",
    format(start_time, "%Y%m%d%H"), "_",
    format(end_time, "%Y%m%d%H"),
    ".png"
  )
)

output_csv <- file.path(
  out_dir,
  paste0(target_station, "_change_points_result.csv")
)


# ----------------------------
# 2. 读取目标站点数据
# ----------------------------

read_one_station <- function(file) {
  
  station_name <- basename(file) |>
    stringr::str_remove("^Beijing_") |>
    stringr::str_remove("\\.csv$")
  
  dat <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  need_vars <- c("year", "month", "day", "hour", "PM2.5", "PM10")
  miss_vars <- setdiff(need_vars, names(dat))
  
  if (length(miss_vars) > 0) {
    stop(
      paste0(
        "文件 ", basename(file),
        " 缺少变量：",
        paste(miss_vars, collapse = ", ")
      )
    )
  }
  
  dat |>
    dplyr::transmute(
      station = station_name,
      datetime = lubridate::make_datetime(
        year = year,
        month = month,
        day = day,
        hour = hour,
        tz = "Asia/Shanghai"
      ),
      PM25 = as.numeric(`PM2.5`),
      PM10 = as.numeric(PM10)
    )
}

files <- list.files(
  data_dir,
  pattern = "^Beijing_.*\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("未找到 Beijing_*.csv 文件，请检查 data_dir。")
}

raw_long <- purrr::map_dfr(files, read_one_station)

if (!(target_station %in% unique(raw_long$station))) {
  stop(
    paste0(
      "未找到目标站点：", target_station,
      "。当前读取到的站点包括：",
      paste(sort(unique(raw_long$station)), collapse = ", ")
    )
  )
}

dat_window_raw <- raw_long |>
  dplyr::filter(
    station == target_station,
    datetime >= start_time,
    datetime <= end_time
  ) |>
  dplyr::arrange(datetime)

if (nrow(dat_window_raw) == 0) {
  stop("筛选后的时间段没有数据，请检查 start_time 和 end_time。")
}

# 补齐小时序列
full_time <- tibble(
  datetime = seq(
    from = min(dat_window_raw$datetime),
    to   = max(dat_window_raw$datetime),
    by   = "hour"
  )
)

dat_window <- full_time |>
  dplyr::left_join(dat_window_raw, by = "datetime") |>
  dplyr::mutate(station = target_station) |>
  dplyr::arrange(datetime)

if (anyNA(dat_window[, c("PM10", "PM25")])) {
  miss_time <- dat_window |>
    dplyr::filter(is.na(PM10) | is.na(PM25)) |>
    dplyr::select(datetime, PM10, PM25)
  
  print(head(miss_time, 20))
  stop("当前时间段 PM10 或 PM2.5 含有缺失值，或存在缺失小时，请更换时间段。")
}

if (nrow(dat_window) < 50) {
  stop("当前时间段样本量小于 50，不建议进行 SNSeg 变点检测。")
}


# ----------------------------
# 3. 标准化与变点提取函数
# ----------------------------

safe_scale_vector <- function(x, var_name = "变量") {
  
  if (anyNA(x)) {
    stop(paste0(var_name, " 含有缺失值。"))
  }
  
  if (sd(x) == 0) {
    stop(paste0(var_name, " 标准差为 0，无法标准化。"))
  }
  
  as.numeric(scale(x))
}

safe_scale_matrix <- function(x) {
  
  if (anyNA(x)) {
    stop("多变量输入矩阵含有缺失值。")
  }
  
  sds <- apply(x, 2, sd)
  
  if (any(sds == 0)) {
    stop(
      paste0(
        "以下变量标准差为 0，无法标准化：",
        paste(colnames(x)[sds == 0], collapse = ", ")
      )
    )
  }
  
  scale(x)
}

extract_cp_from_snseg <- function(res, n) {
  
  if (is.null(res)) return(integer(0))
  
  possible_names <- c(
    "est_cp", "cp", "cpts", "changepoints",
    "change_points", "cp_index"
  )
  
  for (nm in possible_names) {
    if (nm %in% names(res)) {
      cp <- res[[nm]]
      
      if (length(cp) == 0 || all(is.na(cp))) {
        return(integer(0))
      }
      
      cp <- as.integer(round(as.numeric(unlist(cp))))
      cp <- cp[!is.na(cp)]
      cp <- cp[cp >= 1 & cp <= n]
      cp <- sort(unique(cp))
      
      return(cp)
    }
  }
  
  message("未从 SNSeg 返回对象中识别到变点字段。返回对象字段名为：")
  print(names(res))
  
  integer(0)
}


# ----------------------------
# 4. 均值变点检测函数
# ----------------------------

run_pm25_uni_detection <- function(dat) {
  
  x_pm25 <- safe_scale_vector(dat$PM25, "PM2.5")
  
  res <- SNSeg_Uni(
    ts = x_pm25,
    paras_to_test = "mean",
    confidence = sn_confidence,
    grid_size_scale = sn_grid_size_scale,
    grid_size = NULL,
    plot_SN = FALSE,
    est_cp_loc = TRUE
  )
  
  cp_index <- extract_cp_from_snseg(res, nrow(dat))
  
  list(
    result = res,
    cp_index = cp_index,
    cp_datetime = dat$datetime[cp_index],
    error = NA_character_
  )
}

run_pm10_pm25_multi_detection <- function(dat) {
  
  x_multi <- as.matrix(dat[, c("PM10", "PM25")])
  colnames(x_multi) <- c("PM10", "PM2.5")
  
  # 检测时仍使用标准化数据
  x_multi_std <- safe_scale_matrix(x_multi)
  
  res <- SNSeg_Multi(
    ts = x_multi_std,
    paras_to_test = "mean",
    confidence = sn_confidence,
    grid_size_scale = sn_grid_size_scale,
    grid_size = NULL,
    plot_SN = FALSE,
    est_cp_loc = TRUE
  )
  
  cp_index <- extract_cp_from_snseg(res, nrow(dat))
  
  list(
    result = res,
    cp_index = cp_index,
    cp_datetime = dat$datetime[cp_index],
    error = NA_character_
  )
}

empty_detection_result <- function(error_msg) {
  list(
    result = NULL,
    cp_index = integer(0),
    cp_datetime = as.POSIXct(character(0), tz = "Asia/Shanghai"),
    error = error_msg
  )
}


# ----------------------------
# 5. 执行变点检测
# ----------------------------

res_uni_pm25 <- tryCatch(
  run_pm25_uni_detection(dat_window),
  error = function(e) {
    message("PM2.5 单变量检测失败，真实错误为：", e$message)
    empty_detection_result(e$message)
  }
)

res_multi <- tryCatch(
  run_pm10_pm25_multi_detection(dat_window),
  error = function(e) {
    message("PM10-PM2.5 多变量检测失败，真实错误为：", e$message)
    empty_detection_result(e$message)
  }
)

cat("\nPM2.5 单变量均值变点检测结果：\n")
if (!is.na(res_uni_pm25$error)) {
  cat("检测失败：", res_uni_pm25$error, "\n")
} else if (length(res_uni_pm25$cp_index) == 0) {
  cat("未检测到变点。\n")
} else {
  print(data.frame(
    变点序号 = seq_along(res_uni_pm25$cp_index),
    变点位置 = res_uni_pm25$cp_index,
    变点时间 = format(res_uni_pm25$cp_datetime, "%Y-%m-%d %H:%M")
  ))
}

cat("\nPM10 与 PM2.5 多变量均值变点检测结果：\n")
if (!is.na(res_multi$error)) {
  cat("检测失败：", res_multi$error, "\n")
} else if (length(res_multi$cp_index) == 0) {
  cat("未检测到变点。\n")
} else {
  print(data.frame(
    变点序号 = seq_along(res_multi$cp_index),
    变点位置 = res_multi$cp_index,
    变点时间 = format(res_multi$cp_datetime, "%Y-%m-%d %H:%M")
  ))
}


# ----------------------------
# 6. 整理绘图数据
# ----------------------------

# 上图：PM2.5 原始浓度
plot_uni_dat <- dat_window |>
  dplyr::transmute(
    datetime = datetime,
    value = PM25,
    variable = "PM2.5"
  )

# 下图：PM10 和 PM2.5 原始浓度
plot_multi_dat <- dat_window |>
  dplyr::transmute(
    datetime = datetime,
    PM10 = PM10,
    `PM2.5` = PM25
  ) |>
  tidyr::pivot_longer(
    cols = c(PM10, `PM2.5`),
    names_to = "variable",
    values_to = "value"
  )

make_cp_dat <- function(cp_datetime) {
  
  if (length(cp_datetime) == 0) {
    return(tibble(
      cp_datetime = as.POSIXct(character(0), tz = "Asia/Shanghai"),
      cp_x = numeric(0)
    ))
  }
  
  tibble(
    cp_datetime = cp_datetime,
    cp_x = as.numeric(cp_datetime)
  )
}

cp_uni_dat <- make_cp_dat(res_uni_pm25$cp_datetime)
cp_multi_dat <- make_cp_dat(res_multi$cp_datetime)


# ----------------------------
# 7. 颜色与主题设置
# ----------------------------

pollutant_colors <- c(
  "PM10"  = "#B2182B",
  "PM2.5" = "#2166AC"
)

legend_position_use <- ifelse(show_legend, "right", "none")

theme_cp <- theme_minimal(
  base_family = font_family,
  base_size = base_size
) +
  theme(
    plot.title = element_text(
      size = title_size,
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = subtitle_size,
      hjust = 0.5,
      color = "grey35"
    ),
    axis.title = element_text(size = axis_title_size),
    axis.text = element_text(size = axis_text_size),
    
    legend.position = legend_position_use,
    legend.title = element_text(size = legend_title_size, face = "bold"),
    legend.text = element_text(size = legend_text_size),
    
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    axis.line = element_blank(),
    
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(8, 16, 8, 8)
  )


# ----------------------------
# 8. 上图：PM2.5 单变量均值变点检测
# ----------------------------

p_uni <- ggplot() +
  geom_line(
    data = plot_uni_dat,
    aes(x = datetime, y = value, color = variable),
    linewidth = 0.85,
    alpha = 0.95,
    na.rm = TRUE
  ) +
  geom_vline(
    data = cp_uni_dat,
    aes(xintercept = cp_x),
    color = "#4D4D4D",
    linetype = "dashed",
    linewidth = 0.80,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = pollutant_colors,
    name = "变量"
  ) +
  scale_x_datetime(
    date_labels = "%m-%d",
    date_breaks = x_date_breaks
  ) +
  labs(
    title = paste0(target_station_zh, "站 PM2.5 均值变点检测结果"),
    # subtitle = paste0(
    #   "时间段：",
    #   format(start_time, "%Y-%m-%d %H:%M"), " 至 ",
    #   format(end_time, "%Y-%m-%d %H:%M"),
    #   "；灰色虚线为单变量均值变点"
    # ),
    x = NULL,
    y = "浓度"
  ) +
  theme_cp


# ----------------------------
# 9. 下图：PM10 与 PM2.5 多变量均值变点检测
# ----------------------------

p_multi <- ggplot() +
  geom_line(
    data = plot_multi_dat,
    aes(x = datetime, y = value, color = variable),
    linewidth = 0.85,
    alpha = 0.95,
    na.rm = TRUE
  ) +
  geom_vline(
    data = cp_multi_dat,
    aes(xintercept = cp_x),
    color = "#4D4D4D",
    linetype = "solid",
    linewidth = 0.85,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = pollutant_colors,
    name = "变量"
  ) +
  scale_x_datetime(
    date_labels = "%m-%d",
    date_breaks = x_date_breaks
  ) +
  labs(
    title = paste0(target_station_zh, "站 PM10 与 PM2.5 均值变点检测结果"),
    # subtitle = "灰色实线为多变量均值变点",
    x = "时间",
    y = "浓度"
  ) +
  theme_cp


# ----------------------------
# 10. 合并上下两图并保存
# ----------------------------

p_final <- p_uni / p_multi +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = legend_position_use
  )

print(p_final)

ggsave(
  filename = output_image,
  plot = p_final,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

message("图片已保存：", output_image)


# ----------------------------
# 11. 保存变点结果
# ----------------------------

make_cp_result <- function(method_name, res) {
  
  if (!is.na(res$error)) {
    return(tibble(
      方法 = method_name,
      变点位置 = NA_integer_,
      变点时间 = NA_character_,
      备注 = paste0("检测失败：", res$error)
    ))
  }
  
  if (length(res$cp_index) == 0) {
    return(tibble(
      方法 = method_name,
      变点位置 = NA_integer_,
      变点时间 = NA_character_,
      备注 = "未检测到变点"
    ))
  }
  
  tibble(
    方法 = method_name,
    变点位置 = res$cp_index,
    变点时间 = format(res$cp_datetime, "%Y-%m-%d %H:%M:%S"),
    备注 = "检测成功"
  )
}

cp_result <- dplyr::bind_rows(
  make_cp_result("PM2.5单变量均值变点检测", res_uni_pm25),
  make_cp_result("PM10-PM2.5多变量均值变点检测", res_multi)
)

write.csv(
  cp_result,
  output_csv,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("变点结果已保存：", output_csv)