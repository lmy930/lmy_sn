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

SELECTED_STATIONS <- c("Dongsi", "Tiantan", "Guanyuan", "Wanshouxigong", "Nongzhanguan", "Aotizhongxin")

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
    "Aotizhongxin" = "奥体中心"
  ))

stations_sf <- st_as_sf(city_df, coords = c("lon", "lat"), crs = 4326)


geojson_url <- paste0("https://geo.datav.aliyun.com/areas_v3/bound/", ADCODE, "_full.json")
message("正在从阿里云获取 ", TARGET_CITY, " 的地图底图...")
map_bounds <- st_read(geojson_url, quiet = TRUE)

font_family <- "SimHei" 

p <- ggplot() +
  # 绘制底层区划多边形
  geom_sf(data = map_bounds, 
          fill = "#f8f9fa",       #
          color = "#adb5bd",      
          linewidth = 0.6) +
  
  #  绘制站点 
  geom_sf(data = stations_sf, 
          aes(fill = station_name_zh),  
          shape = 21, 
          color = "white",        
          size = 4.5,             
          stroke = 0.8) +         
  
  
  geom_text_repel(data = city_df, 
                  aes(x = lon, y = lat, label = station_name_zh), 
                  family = font_family,
                  size = 3.5,
                  color = "#343a40",
                  fontface = "bold",
                  bg.color = "white",    
                  bg.r = 0.15,
                  box.padding = 0.6,
                  point.padding = 0.6,
                  segment.color = "gray50") +
  

  scale_fill_brewer(palette = "Set1", name = "监测站点") + # 修改点：图例标题中文化

  labs(title = paste("北京市6个空气质量监测站点的分布图"),
       #subtitle = "Spatial distribution of selected observation stations",
       x = NULL, y = NULL) +

  theme_minimal(base_family = font_family) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 12, color = "gray50", hjust = 0.5, margin = margin(b = 20)),
    panel.grid = element_blank(),         
    axis.text = element_blank(),          
    axis.ticks = element_blank(),
    legend.position = "right",           
    legend.title = element_text(face = "bold"),
    plot.background = element_rect(fill = "white", color = NA), 
    panel.background = element_rect(fill = "white", color = NA)
  )


print(p)

ggsave(OUTPUT_IMAGE, plot = p, width = 10, height = 8, dpi = 300)
message("成功生成图片: ", OUTPUT_IMAGE)









library(tidyverse)  
library(lubridate)  
library(viridis)    

INPUT_FOLDER <- "D:/03/2025/Beijing" 
OUTPUT_FOLDER <- "2025_new"            


START_TIME <- ymd_hm("2023-02-18 12:00")
END_TIME <- ymd_hm("2023-02-21 12:00")

POLLUTANTS <- c("PM10", "PM2.5")

FONT_FAMILY <- "SimHei"  

extract_data <- function(input_path, output_path, start_time, end_time) {
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
    #message("[*] 已创建输出文件夹: ", output_path)
  }
  
  csv_files <- list.files(input_path, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("[!] 错误：未找到CSV文件。")
  
  #message(sprintf("[*] 开始提取数据：共找到 %d 个文件...", length(csv_files)))
  
  for (file_path in csv_files) {
    tryCatch({
      df <- read_csv(file_path, show_col_types = FALSE)
      req_cols <- c("year", "month", "day")
      if (!all(req_cols %in% names(df))) next
      
      if (!"hour" %in% names(df)) {
        df$hour <- 0
      } else {
        df$hour <- replace_na(as.numeric(df$hour), 0)
      }
      
      filtered_df <- df %>%
        mutate(datetime = make_datetime(year, month, day, hour)) %>%
        filter(!is.na(datetime), datetime >= start_time, datetime <= end_time)
      
      if (nrow(filtered_df) > 0) {
        file_name <- basename(file_path)
        save_path <- file.path(output_path, file_name)
        write_csv(filtered_df, save_path, na = "NA")
        #message(sprintf("    - 已提取: %s (%d 行)", file_name, nrow(filtered_df)))
      }
    }, error = function(e) {
      #message(sprintf("[!] 处理 %s 时出错: %s", file_path, e$message))
    })
  }
  #message("[*] 数据提取阶段完成。\n")
  return(TRUE)
}


visualize_data <- function(data_dir) {
  #message("[*] 开始生成可视化图表...")
  
  csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  
  all_data <- map_dfr(csv_files, function(file) {
    station_name <- tools::file_path_sans_ext(basename(file))
    df <- read_csv(file, show_col_types = FALSE)
    df %>% mutate(station = station_name)
  })
  
  all_data <- all_data %>%
    mutate(station_zh = case_when(
      str_detect(tolower(station), "dongsi") ~ "东四",
      str_detect(tolower(station), "tiantan") ~ "天坛",
      str_detect(tolower(station), "guanyuan") ~ "官园",
      str_detect(tolower(station), "wanshouxigong") ~ "万寿西宫",
      str_detect(tolower(station), "nongzhanguan") ~ "农展馆",
      str_detect(tolower(station), "aotizhongxin") ~ "奥体中心",
      TRUE ~ station 
    ))
  
  for (pollutant in POLLUTANTS) {
    if (!pollutant %in% names(all_data)) next
    
    plot_data <- all_data %>% filter(!is.na(.data[[pollutant]]))
    if (nrow(plot_data) == 0) next
    
    p <- ggplot(plot_data, aes(x = datetime, y = .data[[pollutant]], color = station_zh)) +
      geom_line(linewidth = 1.2, alpha = 0.8) +
      scale_color_viridis_d(option = "turbo", name = "监测站点") +
      scale_x_datetime(date_labels = "%m-%d %H:%M", date_breaks = "12 hours") +
      labs(
        x = "时间",
        y = paste(pollutant, "浓度 (μg/m³)")
      ) +
      theme_minimal(base_size = 14, base_family = FONT_FAMILY) +
      theme(
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 8)),
        plot.subtitle = element_text(color = "gray40", size = 12, hjust = 0.5, margin = margin(b = 20)),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 15), face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), face = "bold"),
        panel.grid.minor = element_blank(),                       
        panel.grid.major = element_line(color = "gray90"),        
        legend.position = "right",                                
        legend.title = element_text(face = "bold", size = 16), 
        legend.text = element_text(size = 14),
        plot.background = element_rect(fill = "white", color = NA), 
        panel.background = element_rect(fill = "white", color = NA)
      )
    
    output_img <- paste0(pollutant, "_trends_merged.png")
    ggsave(output_img, plot = p, width = 14, height = 7, dpi = 300)
    #message("[*] 图表已保存: ", output_img)
  }
}

# ================== 执行 ==================
success <- extract_data(INPUT_FOLDER, OUTPUT_FOLDER, START_TIME, END_TIME)
if (success) {
  visualize_data(OUTPUT_FOLDER)
 #message("\n[√] 任务全部完成！请查看生成的文件夹和图片。")
}