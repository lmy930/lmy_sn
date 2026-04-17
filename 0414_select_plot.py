import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import warnings

warnings.filterwarnings('ignore')

INPUT_FOLDER = r"D:/03/2025/Beijing"  
OUTPUT_FOLDER = "2025_new"            


START_TIME = "2023-02-18 12:00"
END_TIME = "2023-02-21 12:00"


POLLUTANTS = ['PM10', 'PM2.5']        
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans'] 
plt.rcParams['axes.unicode_minus'] = False
# =================================================

def extract_data(input_path, output_path, start, end):
  
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        print(f"[*] 已创建输出文件夹: {output_path}")

    start_dt = pd.to_datetime(start)
    end_dt = pd.to_datetime(end)
    
    csv_files = glob.glob(os.path.join(input_path, "*.csv"))
    if not csv_files:
        print(f"[!] 错误：在 {input_path} 中未找到CSV文件。")
        return False

    print(f"[*] 开始提取数据：共找到 {len(csv_files)} 个文件...")

    for file_path in csv_files:
        try:
            df = pd.read_csv(file_path)
           
            if not {'year', 'month', 'day'}.issubset(df.columns):
                continue

            
            if 'hour' not in df.columns:
                df['hour'] = 0
            else:
                df['hour'] = df['hour'].fillna(0).astype(int)

            
            df['datetime'] = pd.to_datetime(df[['year', 'month', 'day', 'hour']], errors='coerce')
            
            if df['datetime'].isna().all():
                continue

            
            mask = (df['datetime'] >= start_dt) & (df['datetime'] <= end_dt)
            filtered_df = df[mask].drop(columns=['datetime'])

            if not filtered_df.empty:
                file_name = os.path.basename(file_path)
                save_path = os.path.join(output_path, file_name)
                filtered_df.to_csv(save_path, index=False, na_rep='NA')
                print(f"    - 已提取: {file_name} ({len(filtered_df)} 行)")
            
        except Exception as e:
            print(f"[!] 处理 {file_path} 时出错: {e}")
    
    print("[*] 数据提取阶段完成。\n")
    return True

def visualize_data(data_dir):

    
    print(f"[*] 开始生成可视化图表，读取路径: {data_dir}...")
    
    csv_files = [f for f in os.listdir(data_dir) if f.endswith('.csv')]
    if not csv_files:
        print("[!] 错误：输出文件夹中没有可绘图的数据。")
        return

    csv_files.sort()
    all_stations_data = {}

    for csv_file in csv_files:
        station_name = csv_file.replace('.csv', '')
        df = pd.read_csv(os.path.join(data_dir, csv_file))
        df['datetime'] = pd.to_datetime(df[['year', 'month', 'day', 'hour']])
        df = df.sort_values('datetime')
        all_stations_data[station_name] = df

  
    for pollutant in POLLUTANTS:
        plt.figure(figsize=(14, 7))
        
        has_data = False
        for station_name, df in all_stations_data.items():
            if pollutant in df.columns:
                plt.plot(df['datetime'], df[pollutant], 
                        label=station_name, linewidth=1.5, alpha=0.8)
                has_data = True
        
        if not has_data:
            plt.close()
            continue

        plt.title(f'{pollutant} 浓度变化趋势 ({START_TIME} 至 {END_TIME})', fontsize=16, fontweight='bold')
        plt.xlabel('时间', fontsize=12)
        plt.ylabel(f'{pollutant} 浓度 (μg/m³)', fontsize=12)
        plt.legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=9)
        plt.grid(True, alpha=0.3)
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        output_img = f'{pollutant}_trends_merged.png'
        plt.savefig(output_img, dpi=300, bbox_inches='tight')
        print(f"[*] 图表已保存: {output_img}")
        plt.show()

if __name__ == "__main__":
   
    success = extract_data(INPUT_FOLDER, OUTPUT_FOLDER, START_TIME, END_TIME)
    
    if success:
        visualize_data(OUTPUT_FOLDER)
        print("\n[√] 任务全部完成！请查看生成的文件夹和图片。")