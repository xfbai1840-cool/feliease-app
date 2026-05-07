import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io

# ==========================================
# 0. 页面配置
# ==========================================
st.set_page_config(page_title="FeliEase 3.0 | 双荧光重复分析", layout="wide")

st.title("🧪 FeliEase 3.0 高通量双荧光筛选系统")
st.markdown("""
**专为 Discovery 实验室双荧光筛选定制版：**
- **双重归一化**：$最终移码效率 = (药物 F / 溶剂 F 均值) / (药物 R / 溶剂 R 均值)$。
- **智能毒性排雷**：低于阳性对照 (M/AA列) 的孔位会被自动标记为毒性，**且不参与全局 Z-score 基线计算**。
- **多样本合并**：支持 3 次重复实验自动合并与变异系数 (CV) 评估。
""")

# 侧边栏参数
with st.sidebar:
    st.header("⚙️ 分析参数设置")
    PLATE_ROWS = st.number_input("板内数据行数 (如 A-H 共8行)", value=8)
    st.divider()
    Z_THRESHOLD = st.slider("Hit 判定阈值 (Z-score < ?)", -5.0, 0.0, -2.0, 0.1)
    st.info("注：Z-score < -2 代表化合物对核糖体移码效率有极其显著的抑制作用。")

# ==========================================
# 1. 核心解析函数 (Discovery 实验室特定版式)
# ==========================================
def parse_raw_csv(file):
    """
    左侧 F：B列(1)为溶剂, M列(12)为毒性, C-L列(2:11)为药物
    右侧 R：P列(15)为溶剂, AA列(26)为毒性, Q-Z列(16:25)为药物
    """
    df = pd.read_csv(file, header=None)
    results = []
    
    i = 0
    while i < len(df):
        val_a = str(df.iloc[i, 0]).strip()
        
        # 寻找每板起始标志 (如 ZY-3K-1#-F)
        if val_a and val_a.lower() != 'nan':
            if i + int(PLATE_ROWS) <= len(df):
                plate_id = val_a.replace("-F", "").replace("-R", "").strip()
                
                # --- 1. 计算溶剂对照平均值 (B列 和 P列) ---
                ctrl_f_vals = pd.to_numeric(df.iloc[i : i+int(PLATE_ROWS), 1], errors='coerce').dropna()
                avg_ctrl_f = ctrl_f_vals.mean() if not ctrl_f_vals.empty else 1.0
                
                ctrl_r_vals = pd.to_numeric(df.iloc[i : i+int(PLATE_ROWS), 15], errors='coerce').dropna()
                avg_ctrl_r = ctrl_r_vals.mean() if not ctrl_r_vals.empty else 1.0
                
                # --- 2. 计算毒性底线平均值 (M1~M5 和 AA1~AA5) ---
                tox_f_vals = pd.to_numeric(df.iloc[i : i+5, 12], errors='coerce').dropna()
                min_tox_f = tox_f_vals.mean() if not tox_f_vals.empty else 0
                
                tox_r_vals = pd.to_numeric(df.iloc[i : i+5, 26], errors='coerce').dropna()
                min_tox_r = tox_r_vals.mean() if not tox_r_vals.empty else 0
                
                # --- 3. 提取药物数据并计算 (C~L 和 Q~Z) ---
                for r in range(int(PLATE_ROWS)):
                    for c_idx in range(10): # 每行10个药
                        col_f = 2 + c_idx
                        col_r = 16 + c_idx
                        
                        raw_f = pd.to_numeric(df.iloc[i+r, col_f], errors='coerce')
                        raw_r = pd.to_numeric(df.iloc[i+r, col_r], errors='coerce')
                        
                        if pd.notna(raw_f) and pd.notna(raw_r):
                            # 计算相对荧光素酶信号
                            rel_f = raw_f / avg_ctrl_f
                            rel_r = raw_r / avg_ctrl_r
                            
                            # 核糖体移码效率
                            ratio = rel_f / rel_r if rel_r > 0 else np.nan
                            
                            # 毒性判定
                            is_toxic = "Yes" if (raw_f < min_tox_f or raw_r < min_tox_r) else "No"
                            
                            # 🎯 生成药物编号：
                            # 根据您的描述，Excel的 C1 (r=0, c_idx=0) 对应 编号B1
                            # chr(66) 是 'B'。如果您其实想让它对应 A1，把下面的 66 改成 65 即可。
                            drug_id = f"{chr(66+r)}{c_idx+1}" 
                            
                            results.append({
                                "Plate_ID": plate_id,
                                "药物编号": drug_id,
                                "Raw_F": raw_f,
                                "Raw_R": raw_r,
                                "Rel_F": rel_f,  # 记录相对信号以备核查
                                "Rel_R": rel_r,
                                "Ratio": ratio,
                                "Toxicity": is_toxic
                            })
                i += int(PLATE_ROWS)
            else:
                break
        else:
            i += 1
            
    return pd.DataFrame(results)

# ==========================================
# 2. 主程序逻辑
# ==========================================
col_u1, col_u2 = st.columns([2, 1])
with col_u1:
    uploaded_files = st.file_uploader("1️⃣ 上传重复实验文件 (支持上传3个 CSV 自动合并)", type="csv", accept_multiple_files=True)
with col_u2:
    lib_file = st.file_uploader("2️⃣ 可选：上传化合物库信息表", type=["xlsx", "csv"])

if uploaded_files:
    all_reps = []
    for idx, f in enumerate(uploaded_files):
        rep_df = parse_raw_csv(f)
        all_reps.append(rep_df)
    
    if all_reps:
        combined_raw = pd.concat(all_reps)
        
        # 聚合重复实验数据，使用“药物编号”代替“Well”
        final_df = combined_raw.groupby(['Plate_ID', '药物编号']).agg({
            'Ratio': ['mean', 'std'],
            'Rel_F': 'mean',
            'Rel_R': 'mean',
            'Raw_F': 'mean',
            'Raw_R': 'mean',
            # 只要在任何一次重复中判定为有毒，整体就标记为有毒
            'Toxicity': lambda x: 'Yes' if 'Yes' in x.values else 'No'
        }).reset_index()
        
        final_df.columns = ['Plate_ID', '药物编号', 'Avg_Ratio', 'Std_Ratio', 'Mean_Rel_F', 'Mean_Rel_R', 'Mean_Raw_F', 'Mean_Raw_R', 'Toxicity']
        final_df['CV_%'] = (final_df['Std_Ratio'] / final_df['Avg_Ratio']) * 100
        
        # -----------------------------------------------------
        # 核心改进：剔除毒性孔后，计算全局基线均值和标准差
        # -----------------------------------------------------
        healthy_df = final_df[final_df['Toxicity'] == 'No']
        
        if healthy_df.empty:
            st.error("⚠️ 警告：所有孔位均被判定为有毒！无法计算基线。请检查数据或阳性对照是否异常。")
        else:
            mu = healthy_df['Avg_Ratio'].mean()
            sigma = healthy_df['Avg_Ratio'].std()
            
            # 对所有孔计算 Z-score
            final_df['Z_score'] = (final_df['Avg_Ratio'] - mu) / sigma
            
            # Hit 判定分类逻辑
            def determine_hit(row):
                if row['Toxicity'] == 'Yes':
                    return 'Toxic (Excluded)'
                elif row['Z_score'] < Z_THRESHOLD:
                    return 'Yes (Hit)'
                else:
                    return 'No'
                    
            final_df['Result'] = final_df.apply(determine_hit, axis=1)

            # 合并化合物库
            if lib_file:
                ldf = pd.read_excel(lib_file) if lib_file.name.endswith('xlsx') else pd.read_csv(lib_file)
                # 统一列名以确保匹配，将库中的坐标列与我们的“药物编号”对齐
                ldf.rename(columns={'Physical_Plate': 'Plate_ID', 'Physical_Well': '药物编号', 'Well': '药物编号'}, inplace=True)
                final_df = pd.merge(final_df, ldf, on=['Plate_ID', '药物编号'], how='left')

            # ==========================
            # 界面展示
            # ==========================
            st.divider()
            total_drugs = len(final_df)
            toxic_drugs = len(final_df[final_df['Toxicity'] == 'Yes'])
            hit_drugs = len(final_df[final_df['Result'] == 'Yes (Hit)'])
            
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("总筛药孔数", total_drugs)
            m2.metric("被剔除的毒性假阳性", toxic_drugs)
            m3.metric("健康基线均值 (Mean Ratio)", f"{mu:.3f}")
            m4.metric("有效抑制剂 (Hits)", hit_drugs)

            # 绘图展示
            st.subheader("📊 筛选散点图 (基线与抑制效果)")
            fig, ax = plt.subplots(figsize=(12, 5))
            
            # 分类画点
            hits = final_df[final_df['Result'] == 'Yes (Hit)']
            normal = final_df[final_df['Result'] == 'No']
            toxic = final_df[final_df['Toxicity'] == 'Yes']
            
            ax.scatter(normal.index, normal['Z_score'], c='#A0AEC0', s=15, alpha=0.5, label='Normal')
            ax.scatter(toxic.index, toxic['Z_score'], c='#ED8936', s=15, marker='x', alpha=0.8, label='Toxic (Discarded)')
            ax.scatter(hits.index, hits['Z_score'], c='#E53E3E', s=35, label='Hits')
            
            ax.axhline(Z_THRESHOLD, color='#3182CE', linestyle='--', label=f'Hit Threshold (Z={Z_THRESHOLD})')
            ax.axhline(0, color='black', linewidth=0.8, alpha=0.5)
            
            ax.set_ylabel("Z-score")
            ax.set_xlabel("Compound Sequence Index")
            ax.legend(loc='upper right')
            ax.grid(True, linestyle=':', alpha=0.4)
            st.pyplot(fig)

            # 结果表格
            st.subheader("🚩 强效抑制剂清单 (仅展示健康 Hits)")
            # 优化展示列的顺序，优先展示核心指标
            display_cols = ['Plate_ID', '药物编号', 'Z_score', 'Avg_Ratio', 'CV_%', 'Mean_Rel_F', 'Mean_Rel_R', 'Result']
            # 如果有化合物库的名字，也加进去
            if 'Product Name' in final_df.columns: display_cols.insert(2, 'Product Name')
            
            hits_display = final_df[final_df['Result'] == 'Yes (Hit)'].sort_values('Z_score')
            st.dataframe(hits_display[[c for c in display_cols if c in hits_display.columns]])

            # 下载
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                # 按照类别排序，先看 Hit，再看 Normal，最后是有毒的
                export_df = final_df.sort_values(by=['Result', 'Z_score'], ascending=[True, True])
                export_df.to_excel(writer, index=False, sheet_name='All_Screen_Results')
                
            st.download_button(
                label="📥 下载带有毒性标记的完整报告 (Excel)", 
                data=output.getvalue(), 
                file_name="Discovery_Dual_Luciferase_Report.xlsx",
                type="primary"
            )
