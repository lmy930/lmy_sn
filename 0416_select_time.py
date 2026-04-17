######code to select expected time period,based on correlation and PCA metrics
import os
import glob
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import warnings

warnings.filterwarnings('ignore')

def calculate_metrics(df_window, columns):
    data = df_window[columns]
    
  
    corr_matrix = data.corr()
    upper_triangle = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
    mean_corr = upper_triangle.mean().mean()
    
    if data.std().min() > 1e-6:
        pca = PCA(n_components=1)
       
        scaled_data = (data - data.mean()) / data.std()
        pca.fit(scaled_data)
        pc1_ratio = pca.explained_variance_ratio_[0]
    else:
        pc1_ratio = 0.0
        
    return mean_corr, pc1_ratio

def find_golden_periods(folder_path, window_hours=48, step_hours=12):
 
    print(f"正在读取 {folder_path} 下的 CSV 文件...")
    all_files = glob.glob(os.path.join(folder_path, "*.csv"))
    
    if len(all_files) != 6:
        print(f"警告：找到了 {len(all_files)} 个文件，期望数量是 6 个。")
    
   
    station_dfs = []
    for file in all_files:
        
        station_name = os.path.basename(file).split('.')[0]
        
        df = pd.read_csv(file)
        
      
        df['datetime'] = pd.to_datetime(df[['year', 'month', 'day', 'hour']])
        
        
        df_sub = df[['datetime', 'PM2.5', 'PM10']].set_index('datetime')
        df_sub.columns = [f'PM25_{station_name}', f'PM10_{station_name}']
        
        station_dfs.append(df_sub)
    
   
    print("正在对齐时间序列...")
    merged_df = pd.concat(station_dfs, axis=1)
    
   
    pm10_cols = [col for col in merged_df.columns if col.startswith('PM10_')]
    pm25_cols = [col for col in merged_df.columns if col.startswith('PM25_')]
    
    results = []
    total_hours = len(merged_df)
    
    print(f"开始滑动窗口扫描 (窗口={window_hours}小时, 步长={step_hours}小时)...")
    

    for start_idx in range(0, total_hours - window_hours + 1, step_hours):
        end_idx = start_idx + window_hours
        window_df = merged_df.iloc[start_idx:end_idx]
        
        
        na_ratio = window_df.isna().sum().sum() / (window_df.shape[0] * window_df.shape[1])
        if na_ratio > 0.05:
            continue
            
   
        window_clean = window_df.interpolate(method='linear').bfill().ffill()
        
   
        pm10_corr, pm10_pc1 = calculate_metrics(window_clean, pm10_cols)
        pm25_corr, pm25_pc1 = calculate_metrics(window_clean, pm25_cols)
        
        results.append({
            'Start_Time': window_clean.index[0],
            'End_Time': window_clean.index[-1],
            'NA_Ratio': na_ratio,
            'PM10_Corr': pm10_corr,
            'PM25_Corr': pm25_corr,
            'PM10_PC1': pm10_pc1,
            'PM25_PC1': pm25_pc1,
      
            'Total_Score': pm10_corr + pm25_corr + pm10_pc1 + pm25_pc1
        })
        
    results_df = pd.DataFrame(results)
    

    filtered_df = results_df[
        (results_df['PM10_Corr'] > 0.75) & 
        (results_df['PM25_Corr'] > 0.75) &
        (results_df['PM10_PC1'] > 0.70) & 
        (results_df['PM25_PC1'] > 0.70)
    ]
    

    final_top = filtered_df.sort_values(by='Total_Score', ascending=False).head(10)
    
    print("\n" + "="*80)
    print("✨ 扫描完成！为你找到最完美的黄金 48 小时测试时段：")
    print("="*80)
    

    for i, row in final_top.iterrows():
        print(f"起止时间: {row['Start_Time']} 到 {row['End_Time']}")
        print(f"  -> PM10  表现: 相关系数 = {row['PM10_Corr']:.3f}, PC1占比 = {row['PM10_PC1']:.3f}")
        print(f"  -> PM2.5 表现: 相关系数 = {row['PM25_Corr']:.3f}, PC1占比 = {row['PM25_PC1']:.3f}")
        print("-" * 60)

    return final_top


if __name__ == "__main__":
   
    YOUR_FOLDER_PATH = "D:/03/2025/Beijing" 
    
    best_periods = find_golden_periods(YOUR_FOLDER_PATH, window_hours=240, step_hours=6)



