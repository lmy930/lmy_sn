######可视化站点
library(sf)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# ================== 1. 配置区域 ==================
TARGET_CITY <- "Beijing"
ADCODE <- "110000"
CSV_FILE <- "D:/file/sha/station/station.xlsx"
OUTPUT_IMAGE <- paste0(TARGET_CITY, "_custom_stations.png")

# 选择需要绘制的站点
SELECTED_STATIONS <- c("Dongsi", "Tiantan", "Guanyuan", "Wanshouxigong", "Nongzhanguan", "Aotizhongxin")

# ================== 2. 数据处理 ==================
# 读取 Excel 文件
df <- read_excel(CSV_FILE)

# 过滤城市并去除站名两端的空格
city_df <- df %>%
  mutate(station_name = trimws(station_name)) %>%
  filter(city == TARGET_CITY)

# 如果指定了 SELECTED_STATIONS，则进一步过滤
if (length(SELECTED_STATIONS) > 0) {
  city_df <- city_df %>% filter(station_name %in% SELECTED_STATIONS)
}

if (nrow(city_df) == 0) {
  stop("错误: 筛选后的站点列表为空，请检查数据。")
}

# 将数据框转换为 sf 空间对象 (设置 WGS84 坐标系 EPSG:4326)
stations_sf <- st_as_sf(city_df, coords = c("lon", "lat"), crs = 4326)

# ================== 3. 获取地图边界 ==================
geojson_url <- paste0("https://geo.datav.aliyun.com/areas_v3/bound/", ADCODE, "_full.json")
message("正在从阿里云获取 ", TARGET_CITY, " 的地图底图...")
map_bounds <- st_read(geojson_url, quiet = TRUE)

# ================== 4. ggplot2 绘制高颜值地图 ==================
# 如果地图有中文，请在后面加上 family = "SimHei" (Windows) 或 "STHeiti" (Mac)
font_family <- "sans" 

p <- ggplot() +
  # A. 绘制底层区划多边形
  geom_sf(data = map_bounds, 
          fill = "#f8f9fa",       # 柔和的浅灰背景色
          color = "#adb5bd",      # 边框颜色
          linewidth = 0.6) +
  
  # B. 绘制站点 (使用 21 号形状：可同时修改填充色和边框色)
  geom_sf(data = stations_sf, 
          aes(fill = station_name), 
          shape = 21, 
          color = "white",        # 白色的点描边，增加质感
          size = 4.5,             # 点的大小
          stroke = 0.8) +         # 描边粗细
  
  # C. 添加智能文本标签 (带有引线且防重叠)
  geom_text_repel(data = city_df, 
                  aes(x = lon, y = lat, label = station_name),
                  family = font_family,
                  size = 3.5,
                  color = "#343a40",
                  fontface = "bold",
                  bg.color = "white",     # 文字外发光效果（需 ggplot2 版本较新）
                  bg.r = 0.15,
                  box.padding = 0.6,
                  point.padding = 0.6,
                  segment.color = "gray50") +
  
  # D. 颜色映射 (使用经典的 Set1 色板)
  scale_fill_brewer(palette = "Set1", name = "Stations") +
  
  # E. 标题与排版
  labs(#title = paste("Station Locations:", TARGET_CITY),
       #subtitle = "Spatial distribution of selected observation stations",
       x = NULL, y = NULL) +
  
  # F. 极简美化主题
  theme_minimal(base_family = font_family) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, color = "gray50", hjust = 0.5, margin = margin(b = 20)),
    panel.grid = element_blank(),         # 隐藏经纬度网格线
    axis.text = element_blank(),          # 隐藏经纬度刻度文本
    axis.ticks = element_blank(),
    legend.position = "right",            # 图例放在右侧
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA), # 确保输出背景纯白
    panel.background = element_rect(fill = "white", color = NA)
  )

# ================== 5. 保存和展示 ==================
print(p)

# 保存高质量图片
ggsave(OUTPUT_IMAGE, plot = p, width = 10, height = 8, dpi = 300)
message("成功生成图片: ", OUTPUT_IMAGE)






######可视化污染物
library(tidyverse)  # 包含 readr, dplyr, ggplot2, purrr 等
library(lubridate)  # 专门用于处理时间日期的神器
library(viridis)    # 高级科学配色包

# ================== 1. 配置区域 ==================
INPUT_FOLDER <- "D:/03/2025/Beijing" 
OUTPUT_FOLDER <- "2025_new"            

# 注意：使用 lubridate::ymd_hm 解析不带秒的时间字符串
START_TIME <- ymd_hm("2023-02-18 12:00")
END_TIME <- ymd_hm("2023-02-21 12:00")

POLLUTANTS <- c("PM10", "PM2.5")

# 设置字体族 (如果图表中有中文，Windows建议用 "SimHei"，Mac建议用 "STHeiti")
# 如果只有英文，使用默认的 "sans" 即可
FONT_FAMILY <- "sans" 
# =================================================

# --- 核心功能 1：数据提取 ---
extract_data <- function(input_path, output_path, start_time, end_time) {
  # 创建输出目录
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    message("[*] 已创建输出文件夹: ", output_path)
  }
  
  # 获取所有 csv 文件路径
  csv_files <- list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    stop("[!] 错误：在 ", input_path, " 中未找到CSV文件。")
  }
  
  message(sprintf("[*] 开始提取数据：共找到 %d 个文件...", length(csv_files)))
  
  for (file_path in csv_files) {
    tryCatch({
      # 1. 读取数据 (使用 read_csv 比 base R 的 read.csv 更快)
      df <- read_csv(file_path, show_col_types = FALSE)
      
      # 2. 检查必要列
      req_cols <- c("year", "month", "day")
      if (!all(req_cols %in% names(df))) next
      
      # 3. 处理 hour 列 (缺失补 0)
      if (!"hour" %in% names(df)) {
        df$hour <- 0
      } else {
        df$hour <- replace_na(as.numeric(df$hour), 0)
      }
      
      # 4. 生成时间序列并过滤
      filtered_df <- df %>%
        mutate(datetime = make_datetime(year, month, day, hour)) %>%
        filter(!is.na(datetime),
               datetime >= start_time,
               datetime <= end_time)
      
      # 5. 保存结果 (为了后续画图方便，这里保留了 datetime 列)
      if (nrow(filtered_df) > 0) {
        file_name <- basename(file_path)
        save_path <- file.path(output_path, file_name)
        
        write_csv(filtered_df, save_path, na = "NA")
        message(sprintf("    - 已提取: %s (%d 行)", file_name, nrow(filtered_df)))
      }
      
    }, error = function(e) {
      message(sprintf("[!] 处理 %s 时出错: %s", file_path, e$message))
    })
  }
  
  message("[*] 数据提取阶段完成。\n")
  return(TRUE)
}

# --- 核心功能 2：高颜值可视化 ---
visualize_data <- function(data_dir) {
  message("[*] 开始生成可视化图表，读取路径: ", data_dir, "...")
  
  csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    stop("[!] 错误：输出文件夹中没有可绘图的数据。")
  }
  
  # ggplot2 喜欢“长数据”。这里我们将所有站点的 CSV 读取并合并为一个大数据框
  all_data <- map_dfr(csv_files, function(file) {
    station_name <- tools::file_path_sans_ext(basename(file))
    df <- read_csv(file, show_col_types = FALSE)
    # 增加一列记录站点名称
    df %>% mutate(station = station_name)
  })
  
  # 对指定的污染物循环绘图
  for (pollutant in POLLUTANTS) {
    # 检查该污染物是否存在于数据中
    if (!pollutant %in% names(all_data)) next
    
    # 过滤掉该污染物完全缺失的数据点
    plot_data <- all_data %>% filter(!is.na(.data[[pollutant]]))
    if (nrow(plot_data) == 0) next
    
    # 绘图逻辑
    p <- ggplot(plot_data, aes(x = datetime, y = .data[[pollutant]], color = station)) +
      # 绘制折线图，设置线宽和透明度防遮挡
      geom_line(linewidth = 1.2, alpha = 0.8) +
      # 使用 viridis 的 turbo 调色板，色彩分辨度极高，适合多站点
      scale_color_viridis_d(option = "turbo", name = "Station") +
      # 格式化 X 轴时间刻度
      scale_x_datetime(date_labels = "%m-%d %H:%M", date_breaks = "12 hours") +
      # 标题与轴标签
      labs(
        #title = paste("Concentration Trends of", pollutant),
        #subtitle = paste("Observation Period:", format(START_TIME, "%Y-%m-%d %H:%M"), 
        #                 "to", format(END_TIME, "%Y-%m-%d %H:%M")),
        x = "Date & Time",
        y = paste(pollutant, "(μg/m³)")
      ) +
      # 使用极简主题为底版
      theme_minimal(base_size = 14, base_family = FONT_FAMILY) +
      # 精细排版微调
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 8)),
        plot.subtitle = element_text(color = "gray40", size = 12, hjust = 0.5, margin = margin(b = 20)),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
        panel.grid.minor = element_blank(),                       # 移除次级网格线
        panel.grid.major = element_line(color = "gray90"),        # 柔化主网格线
        legend.position = "right",                                # 图例放右侧
        legend.title = element_text(face = "bold"),
        plot.background = element_rect(fill = "white", color = NA), # 确保输出背景纯白
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    # 保存图片 (高分辨率 300dpi)
    output_img <- paste0(pollutant, "_trends_merged.png")
    ggsave(output_img, plot = p, width = 14, height = 7, dpi = 300)
    message("[*] 图表已保存: ", output_img)
  }
}

# ================== 执行 ==================
success <- extract_data(INPUT_FOLDER, OUTPUT_FOLDER, START_TIME, END_TIME)
if (success) {
  visualize_data(OUTPUT_FOLDER)
  message("\n[√] 任务全部完成！请查看生成的文件夹和图片。")
}