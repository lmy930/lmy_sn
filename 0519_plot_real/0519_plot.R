# ============================================================
# 北京市6个空气质量监测站点 PM2.5 和 PM10 小时数据可视化
# ============================================================

library(tidyverse)
library(lubridate)

# -----------------------------
# 1. 基本参数设置
# -----------------------------

data_dir <- "D:/03/2025/Beijing/"

# 后续只需要改这里控制绘图时间段
start_time <- ymd_hms("2024-12-07 00:00:00")
end_time   <- ymd_hms("2024-12-16 23:00:00")

# 字体设置
font_family <- "Microsoft YaHei"

# 字体大小
base_size   <- 14
title_size  <- 18
axis_size   <- 13
legend_size <- 12

# 横坐标时间间隔
x_date_breaks <- "3 days"

# 输出文件夹
out_dir <- "D:/lmy_sn"
#file.path(data_dir, "figures")
#dir.create(out_dir, showWarnings = FALSE)


# -----------------------------
# 2. 站点中文名称设置
# -----------------------------

station_name_map <- c(
  "Aotizhongxin"   = "奥体中心",
  "Changping"     = "昌平",
  "Dingling"      = "定陵",
  "Dongsi"        = "东四",
  "Guanyuan"      = "官园",
  "Gucheng"       = "古城",
  "Huairou"       = "怀柔",
  "Nongzhanguan"  = "农展馆",
  "Shunyi"        = "顺义",
  "Tiantan"       = "天坛",
  "Wanliu"        = "万柳",
  "Wanshouxigong" = "万寿西宫"
)


# -----------------------------
# 3. 读取 CSV 文件
# -----------------------------

csv_files <- list.files(
  path = data_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)

if (length(csv_files) == 0) {
  stop("指定文件夹下没有找到 CSV 文件，请检查路径是否正确。")
}

read_one_station <- function(file_path) {
  
  station_raw <- tools::file_path_sans_ext(basename(file_path))
  
  station_key <- station_raw
  
  for (nm in names(station_name_map)) {
    if (str_detect(station_raw, regex(nm, ignore_case = TRUE))) {
      station_key <- nm
      break
    }
  }
  
  station_cn <- ifelse(
    station_key %in% names(station_name_map),
    station_name_map[station_key],
    station_raw
  )
  
  df <- read_csv(
    file_path,
    show_col_types = FALSE,
    locale = locale(encoding = "UTF-8")
  )
  
  if ("datetime" %in% names(df)) {
    df <- df %>%
      mutate(datetime = ymd_hms(datetime, quiet = TRUE))
  } else {
    df <- df %>%
      mutate(
        datetime = make_datetime(
          year = year,
          month = month,
          day = day,
          hour = hour
        )
      )
  }
  
  df %>%
    mutate(站点 = station_cn) %>%
    select(站点, datetime, PM2.5, PM10)
}

air_data <- map_dfr(csv_files, read_one_station)


# -----------------------------
# 4. 数据筛选
# -----------------------------

plot_data <- air_data %>%
  mutate(
    PM2.5 = as.numeric(PM2.5),
    PM10  = as.numeric(PM10)
  ) %>%
  filter(
    datetime >= start_time,
    datetime <= end_time
  )

if (nrow(plot_data) == 0) {
  stop("筛选后的数据为空，请检查 start_time 和 end_time 是否在数据时间范围内。")
}


# -----------------------------
# 5. 设置站点颜色
# -----------------------------

station_levels <- sort(unique(plot_data$站点))
plot_data$站点 <- factor(plot_data$站点, levels = station_levels)

# 6个站点使用区分度较高、视觉比较柔和的颜色
station_colors_6 <- c(
  "#0072B2",  # 蓝色
  "#D55E00",  # 橙红色
  "#009E73",  # 绿色
  "#CC79A7",  # 紫粉色
  "#E69F00",  # 金黄色
  "#56B4E9"   # 浅蓝色
)

# 如果刚好是6个站点，直接使用这6种颜色；
# 如果不是6个，则自动生成对应数量的颜色，避免报错。
if (length(station_levels) == 6) {
  station_colors <- setNames(station_colors_6, station_levels)
} else {
  station_colors <- setNames(
    colorRampPalette(station_colors_6)(length(station_levels)),
    station_levels
  )
}


# -----------------------------
# 6. 通用主题：彻底去掉四周黑色边框
# -----------------------------

theme_air <- theme_minimal(
  base_family = font_family,
  base_size = base_size
) +
  theme(
    plot.title = element_text(
      size = title_size,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(size = axis_size),
    axis.text = element_text(size = axis_size - 1),
    
    legend.position = "right",
    legend.title = element_text(size = legend_size),
    legend.text = element_text(size = legend_size),
    
    # 关键：去掉绘图区四周的长方形黑色边框
    panel.border = element_blank(),
    
    # 关键：去掉整张图背景边框
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    
    # 不画坐标轴黑线，避免看起来像边框
    axis.line = element_blank(),
    
    # 网格线保留浅灰色，更适合时间序列图
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(10, 20, 10, 10)
  )


# -----------------------------
# 7. 绘制 PM2.5 图
# -----------------------------
#\n%H:%M"
p_pm25 <- ggplot(
  plot_data,
  aes(x = datetime, y = PM2.5, color = 站点)
) +
  geom_line(linewidth = 0.75, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(values = station_colors) +
  scale_x_datetime(
    date_labels = "%m-%d",
    date_breaks = x_date_breaks
  ) +
  labs(
    title = "北京市空气质量监测站点的PM2.5浓度",
    x = "时间",
    y = "PM2.5浓度",
    color = "监测站点"
  ) +
  theme_air

print(p_pm25)


# -----------------------------
# 8. 绘制 PM10 图
# -----------------------------
#\n%H:%M
p_pm10 <- ggplot(
  plot_data,
  aes(x = datetime, y = PM10, color = 站点)
) +
  geom_line(linewidth = 0.75, alpha = 0.9, na.rm = TRUE) +
  scale_color_manual(values = station_colors) +
  scale_x_datetime(
    date_labels = "%m-%d",
    date_breaks = x_date_breaks
  ) +
  labs(
    title = "北京市空气质量监测站点的PM10浓度",
    x = "时间",
    y = "PM10浓度",
    color = "监测站点"
  ) +
  theme_air

print(p_pm10)


# -----------------------------
# 9. 保存图片
# -----------------------------

ggsave(
  filename = file.path(out_dir, "北京市六站点_PM25_小时浓度变化.png"),
  plot = p_pm25,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(out_dir, "北京市六站点_PM10_小时浓度变化.png"),
  plot = p_pm10,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "white"
)
























library(sf)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

TARGET_CITY <- "Beijing"
ADCODE <- "110000"
CSV_FILE <- "D:/file/sha/station/station.xlsx"
OUTPUT_IMAGE <- paste0(TARGET_CITY, "_custom_stations.png")

SELECTED_STATIONS <- c( "Tiantan", "Guanyuan", "Wanshouxigong", "Nongzhanguan", "Aotizhongxin",
                       "Shunyi")

# ============================================================
# 字体与字号设置：后续主要修改这里即可
# ============================================================

# 图片标题
plot_title <- "北京市12个空气质量监测站点分布图"

# 中文字体名称
# Windows 常用："SimHei"、"SimSun"、"Microsoft YaHei"
title_font_family  <- "SimHei"             # 标题字体
legend_font_family <- "SimHei"             # 图例字体
label_font_family  <- "SimHei"             # 图中站点名称字体
base_font_family   <- "SimHei"             # 整体基础字体

# 字号设置
title_size        <- 19
legend_title_size <- 14
legend_text_size  <- 14
label_size        <- 4

df <- read_excel(CSV_FILE)

city_df <- df %>%
  mutate(station_name = trimws(station_name)) %>%
  filter(city == TARGET_CITY)

if (length(SELECTED_STATIONS) > 0) {
  city_df <- city_df %>% filter(station_name %in% SELECTED_STATIONS)
}

city_df <- city_df %>%
  mutate(station_name_zh = recode(station_name,
    "Dongsi" = "东四",
    "Tiantan" = "天坛",
    "Guanyuan" = "官园",
    "Wanshouxigong" = "万寿西宫",
    "Nongzhanguan" = "农展馆",
    "Aotizhongxin" = "奥体中心",
    "Changping" = "昌平",
    #"Dingling" = "定陵",
    "Gucheng" = "古城",
    "Huairou" = "怀柔",
    "Shunyi" = "顺义",
    "Wanliu" = "万柳"
  ))

stations_sf <- st_as_sf(city_df, coords = c("lon", "lat"), crs = 4326)

geojson_url <- paste0("https://geo.datav.aliyun.com/areas_v3/bound/", ADCODE, "_full.json")
message("正在从阿里云获取 ", TARGET_CITY, " 的地图底图...")
map_bounds <- st_read(geojson_url, quiet = TRUE)

p <- ggplot() +
  # 绘制底层区划多边形
  geom_sf(
    data = map_bounds,
    fill = "#f8f9fa",
    color = "#adb5bd",
    linewidth = 0.6
  ) +
  
  # 绘制站点
  geom_sf(
    data = stations_sf,
    aes(fill = station_name_zh),
    shape = 21,
    color = "white",
    size = 4.5,
    stroke = 0.8
  ) +
  
  # 绘制站点中文标签
  geom_text_repel(
    data = city_df,
    aes(x = lon, y = lat, label = station_name_zh),
    family = label_font_family,
    size = label_size,
    color = "#343a40",
    fontface = "bold",
    bg.color = "white",
    bg.r = 0.15,
    box.padding = 0.6,
    point.padding = 0.6,
    segment.color = "gray50"
  ) +
  
  scale_fill_brewer(
    palette = "Set1",
    name = "监测站点"
  ) +
  
  labs(
    title = plot_title,
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal(base_family = base_font_family) +
  theme(
    plot.title = element_text(
      family = title_font_family,
      size = title_size,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    
    plot.subtitle = element_text(
      size = 12,
      color = "gray50",
      hjust = 0.5,
      margin = margin(b = 20)
    ),
    
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    
    legend.position = "right",
    legend.title = element_text(
      family = legend_font_family,
      size = legend_title_size,
      face = "bold"
    ),
    legend.text = element_text(
      family = legend_font_family,
      size = legend_text_size
    ),
    
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p)

ggsave(
  OUTPUT_IMAGE,
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

message("成功生成图片: ", OUTPUT_IMAGE)












# 1. 加载必要的包
library(leaflet)
library(dplyr)

# 2. 准备你的数据 (直接复制你原有的数据)
station_coord <- tibble::tribble(
  ~station, ~station_zh, ~lon, ~lat,
  "Dongsi", "东四", 116.417, 39.929,
  "Tiantan", "天坛", 116.407, 39.886,
  "Guanyuan", "官园", 116.339, 39.929,
  "Wanshouxigong", "万寿西宫", 116.352, 39.878,
  "Nongzhanguan", "农展馆", 116.461, 39.937,
  "Aotizhongxin", "奥体中心", 116.397, 39.982,
  "Changping", "昌平", 116.231, 40.217,
  "Dingling", "定陵", 116.220, 40.292,
  "Gucheng", "古城", 116.184, 39.914,
  "Huairou", "怀柔", 116.628, 40.328,
  "Shunyi", "顺义", 116.655, 40.127,
  "Wanliu", "万柳", 116.287, 39.987
)

# 3. 为了直观，我们把你特别关注的"东四"站标记为红色，其他站标记为蓝色
station_coord <- station_coord %>%
  mutate(point_color = ifelse(station == "Dongsi", "red", "#3388ff"))

# 4. 生成可视化交互地图
map <- leaflet(data = station_coord) %>%
  # 添加一个干净漂亮的底图 (CartoDB底图，颜色较淡，适合突出数据点)
  addProviderTiles(providers$CartoDB.Positron) %>%
  # 在地图上添加圆点标记
  addCircleMarkers(
    lng = ~lon,               # 经度
    lat = ~lat,               # 纬度
    color = ~point_color,     # 边框颜色
    fillColor = ~point_color, # 填充颜色
    fillOpacity = 0.8,        # 透明度
    radius = 6,               # 圆点大小
    # 鼠标放上去时显示的标签（中文名）
    label = ~station_zh,      
    # 标签的样式设置（加大字号，更清楚）
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", "font-size" = "14px", padding = "3px 8px")
    )
  )

# 5. 显示地图
map



# 1. 加载必要的包
library(leaflet)
library(dplyr)

# 2. 准备成都市空气质量国控点数据
# 这里收录了成都最经典的 8 个环境空气质量国控点
station_coord <- tibble::tribble(
  ~station,         ~station_zh,  ~lon,    ~lat,
  "Jinquanlianghe", "金泉两河",   103.973, 30.7236,
  "Shilidian",      "十里店",     104.142, 30.6764,
  "Sanwayao",       "三瓦窑",     104.059, 30.5767,
  "Shahepu",        "沙河铺",     104.112, 30.6306,
  #"Caotangsi",      "草堂寺",     104.026, 30.6544,
  #"Liangjiaxiang",  "梁家巷",     104.074, 30.685,
  "Junpingjie",      "君平街",     104.043, 30.6556,
  #"Lingyansi",      "灵岩寺(对照点)", 103.62, 31.0201, # 位于都江堰的远郊清洁对照点
  "Dashixilu",      "大石西路",   104.022, 30.6558,
  #"Longquanyiqu",     "龙泉驿区",     104.273, 30.5589,
  "Jinbolu",        "金博路",     104.011, 30.6931,
  #"Linjianglu",      "临江路",     104.846, 30.6994,
  #"Huayang",        "华阳",       104.058, 30.5225,
  #"Chechengdongqilu", "车城东七路", 104.239, 30.5319,
  #"Qingbaijiangjishifenyuan", "青白江技师分院", 104.385, 30.7453,
  #"Kexuecheng",      "科学城",     104.079, 30.4072
)

# 3. 假设你特别关注的站点是"三瓦窑" (Sanwayao)，将其标记为红色，其他标记为蓝色
target_station_name <- "Sanwayao" 

station_coord <- station_coord %>%
  mutate(point_color = ifelse(station == target_station_name, "red", "#3388ff"))

# 4. 生成可视化交互地图
map <- leaflet(data = station_coord) %>%
  # 添加一个干净漂亮的底图 (CartoDB底图，颜色较淡，不抢戏)
  addProviderTiles(providers$CartoDB.Positron) %>%
  # 在地图上添加圆点标记
  addCircleMarkers(
    lng = ~lon,               # 经度
    lat = ~lat,               # 纬度
    color = ~point_color,     # 边框颜色
    fillColor = ~point_color, # 填充颜色
    fillOpacity = 0.8,        # 透明度
    radius = 6,               # 圆点大小
    # 鼠标放上去时显示的标签（中文名）
    label = ~station_zh,      
    # 标签的样式设置（加大字号，更清楚）
    labelOptions = labelOptions(
      style = list("font-weight" = "normal", "font-size" = "14px", padding = "3px 8px")
    )
  )

# 5. 显示地图
map