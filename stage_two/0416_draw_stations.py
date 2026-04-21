# import pandas as pd
# import geopandas as gpd
# import matplotlib.pyplot as plt
# import requests
# import os

# # ================== 配置区域 ==================
# # 1. 设置目标城市（需匹配 CSV 中的 city 列）
# # 推荐尝试: "Beijing" 或 "Shanghai"
# TARGET_CITY = "Beijing" 

# # 2. 设置对应的 ADCODE (用于从阿里云获取地图边界)
# # 北京: 110000 | 上海: 310000 | 浙江: 330000 | 宁波: 330200
# ADCODE = "110000" 

# # 3. 文件路径

# CSV_FILE = r"D:/file/sha/station/station.xlsx"
# OUTPUT_IMAGE = f"{TARGET_CITY}_station_map.png"
# # =============================================

# def draw_static_map():
#     # --- 1. 读取站点数据 ---
#     df = pd.read_excel(CSV_FILE)
#     city_df = df[df['city'] == TARGET_CITY].copy()
    
#     if city_df.empty:
#         print(f"错误: 在 CSV 中找不到城市 '{TARGET_CITY}'。请检查拼写。")
#         return

#     # --- 2. 获取地图 JSON (GeoJSON) ---
#     # 我们使用 DataV.GeoAtlas 提供的 API
#     # '_full' 代表包含下级区县边界
#     geojson_url = f"https://geo.datav.aliyun.com/areas_v3/bound/{ADCODE}_full.json"
    
#     print(f"正在从网络加载 {TARGET_CITY} 的地图边界...")
#     try:
#         # 直接使用 geopandas 读取远程 URL
#         gdf = gpd.read_file(geojson_url)
#     except Exception as e:
#         print(f"地图加载失败，请检查网络或 ADCODE: {e}")
#         return

#     # --- 3. 开始绘图 ---
#     # 设置中文字体（如果是中文名建议设置，英文名则无需）
#     plt.rcParams['font.sans-serif'] = ['SimHei'] 
#     plt.rcParams['axes.unicode_minus'] = False

#     fig, ax = plt.subplots(figsize=(12, 10))

#     # A. 绘制底图（行政区划背景）
#     gdf.plot(
#         ax=ax, 
#         color='#f7f7f7',      # 地图填充色
#         edgecolor='#333333',  # 边界颜色
#         linewidth=0.8,
#         alpha=0.8
#     )

#     # B. 绘制站点（每个站点一个颜色）
#     # 获取唯一的站点列表
#     stations = city_df['station_name'].unique()
#     # 使用自带的颜色映射方案（如 tab20 或 set3）
#     colormap = plt.cm.get_cmap('tab20', len(stations))

#     for i, name in enumerate(stations):
#         data = city_df[city_df['station_name'] == name]
#         ax.scatter(
#             data['lon'], data['lat'], 
#             label=name, 
#             color=colormap(i), 
#             s=80,                # 点的大小
#             edgecolors='black',  # 点的描边，更醒目
#             linewidths=0.5,
#             zorder=5             # 确保点在地图层之上
#         )

#     # --- 4. 美化与细节 ---
#     ax.set_title(f"{TARGET_CITY} 观测站点地理分布图", fontsize=18, pad=20)
    
#     # 隐藏经纬度坐标轴（让图看起来更专业）
#     ax.set_axis_off() 

#     # 动态调整图例：如果站点很多，图例会自动换列显示
#     ncol = 1 if len(stations) < 15 else 2
#     ax.legend(
#         title="站点名称", 
#         bbox_to_anchor=(1.05, 1), 
#         loc='upper left', 
#         fontsize=9,
#         ncol=ncol,
#         frameon=True
#     )

#     # --- 5. 保存图片 ---
#     plt.tight_layout()
#     plt.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches='tight')
#     print(f"成功导出图片：{OUTPUT_IMAGE}")
#     plt.show()

# if __name__ == "__main__":
#     draw_static_map()





######to draw station map and look the distance between stations
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import requests
import os

TARGET_CITY = "Beijing"   # 选项: "Beijing" 或 "Shanghai"
ADCODE = "110000"        # 北京: 110000 | 上海: 310000


# SELECTED_STATIONS = ["Dongsi", "Tiantan", "Guanyuan", "Wanshouxigong", "Nongzhanguan", "Aotizhongxin", 
# "Haidianwanliu", "Gucheng", "Haidianbeianhe", "Shijingshan", "Fengtaihuayuan", "Yungang", "Daxing", 
# "Tongzhou", "Shunyi", "Changping", "Mentougou", "Fangshan", "Pinggu", "Huairou", "Miyun", "Yanqing"]

SELECTED_STATIONS = ["Dongsi", "Tiantan", "Guanyuan", "Wanshouxigong", "Nongzhanguan", "Aotizhongxin"]

# SELECTED_STATIONS = [
#     # 丰台区 (Fengtai)
#     "Fengtaihuayuan",
#       "Nansanhuan", 
#       "Yunzhiguan",
#     # 石景山区 (Shijingshan)
#     "Gucheng",
#     # 门头沟区 (Mentougou)
#     "Longquan",
#     # 房山区 (Fangshan)
#     "Liangxiang",
#     # 大兴区 (Daxing)
#     "Huangcunzhen",
#     # 通州区 (Tongzhou)
#     "Tongzhouxincheng",
#     # 顺义区 (Shunyi)
#     "Shunyi",
#     # 昌平区 (Changping)
#     "Changping",
#     # 延庆区 (Yanqing)
#     "Yanqing",
#     # 怀柔区 (Huairou)
#     "Huairou",
#     # 密云区 (Miyun)
#     "Miyun",
#     # 平谷区 (Pinggu)
#     "Pinggu"
# ]



CSV_FILE = r"D:/file/sha/station/station.xlsx"
OUTPUT_IMAGE = f"{TARGET_CITY}_custom_stations.png"

def draw_custom_station_map():
    if not os.path.exists(CSV_FILE):
        print(f"错误: 找不到文件 {CSV_FILE}")
        return
        
    df = pd.read_excel(CSV_FILE)
  
    df['station_name'] = df['station_name'].str.strip()
    
  
    all_city_data = df[df['city'] == TARGET_CITY]
    
    if all_city_data.empty:
        print(f"错误: 在文件中未找到城市 '{TARGET_CITY}'。")
        return

    print(f"--- 城市 {TARGET_CITY} 的所有可用站点名如下 ---")
    print(list(all_city_data['station_name'].unique()))
    print("-" * 40)


    if SELECTED_STATIONS:
        city_df = all_city_data[all_city_data['station_name'].isin(SELECTED_STATIONS)].copy()
  
        missing = set(SELECTED_STATIONS) - set(city_df['station_name'])
        if missing:
            print(f"注意: 以下站点在文件中未找到，请检查拼写: {missing}")
    else:
        city_df = all_city_data.copy()

    if city_df.empty:
        print("错误: 筛选后的站点列表为空，请检查 SELECTED_STATIONS 中的拼写。")
        return

    geojson_url = f"https://geo.datav.aliyun.com/areas_v3/bound/{ADCODE}_full.json"
    print(f"正在从网络获取 {TARGET_CITY} 的地图底图...")
    try:
        gdf = gpd.read_file(geojson_url)
    except Exception as e:
        print(f"地图下载失败: {e}")
        return

    fig, ax = plt.subplots(figsize=(14, 10))

    gdf.plot(ax=ax, color='#f8f9fa', edgecolor='#495057', linewidth=0.7)


    stations_to_draw = city_df['station_name'].unique()

    cmap = plt.cm.get_cmap('Set1', len(stations_to_draw))

    for i, name in enumerate(stations_to_draw):
        data = city_df[city_df['station_name'] == name]
        ax.scatter(
            data['lon'], data['lat'], 
            label=name, 
            color=cmap(i), 
            s=130,               
            edgecolors='black',  
            linewidths=0.8,
            zorder=10            
        )


    ax.set_title(f"Station Locations: {TARGET_CITY}", fontsize=18, pad=20)
    ax.set_axis_off() 

   
    ax.legend(
        title="Stations", 
        bbox_to_anchor=(1.02, 1), 
        loc='upper left', 
        fontsize=10,
        frameon=True,
        ncol=1 if len(stations_to_draw) < 15 else 2
    )

   
    plt.tight_layout()
    plt.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches='tight')
    print(f"成功生成图片: {OUTPUT_IMAGE}")
    plt.show()

if __name__ == "__main__":
    draw_custom_station_map()