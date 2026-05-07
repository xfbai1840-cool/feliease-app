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
该版本专为**重复实验**设计：
1. **自动归一化**：计算 $Ratio = Firefly / Renilla$。
2. **多样本合并**：上传 3 个重复文件，自动按孔位计算均值与偏差。
3. **Z-score 统计**：基于整板数据计算 $Z = (x - \mu) / \sigma$，科学评估抑制强度。
""")

# 侧边栏参数
with st.sidebar:
    st.header("⚙️ 分析参数设置")
    PLATE_ROWS = st.number_input("板内数据行数 (如 B-M 共12行则填12)", value=8)
    st.divider()
    Z_THRESHOLD = st.slider("Hit 判定阈值 (Z-score < ?)", -5.0, 0.0, -2.0, 0.1)
    st.info("注：Z-score < -2 通常代表具有 95% 以上的统计学显著抑制。")

# ==========================================
# 1. 核心解析函数 (适配原始 96 孔板 CSV)
# ==========================================
def parse_raw_csv(file):
    """解析原始导出的 CSV，提取 ID, F 和 R"""
    df = pd.read_csv(file, header=None)
    results = []
    
    # 自动找板逻辑
    i = 0
    plate_count = 1
    while i < len(df):
        val = str(df.iloc[i, 0]).strip()
        if val and val.lower() != 'nan':  # A列有板ID
            if i + int(PLATE_ROWS) <= len(df):
                plate_id = val
                # 提取 F (B列=1) 和 R (M列=12) 以及中间的药物孔 (C-L列 = 2:12)
                # 简化逻辑：这里假设用户上传的是已经简单清洗或标准的 96 孔数据分布
                # 我们遍历 C-L 列 (药物孔)
                for r in range(int(PLATE_ROWS)):
                    for c in range(2, 12): # C列到L列
                        f_val = pd.to_numeric(df.iloc[i+r, 1], errors='coerce')  # B列作为单行内参或对照
                        r_val = pd.to_numeric(df.iloc[i+r, 12], errors='coerce') # M列作为海神内参
                        drug_f = pd.to_numeric(df.iloc[i+r, c], errors='coerce')
                        
                        if not np.isnan(drug_f) and r_val > 0:
                            well_id = f"{chr(65+r)}{c+1:02d}"
                            results.append({
                                "Plate_ID": plate_id,
                                "Well": well_id,
                                "F": drug_f,
                                "R": r_val,
                                "Ratio": drug_f / r_val
                            })
                plate_count += 1
                i += int(PLATE_ROWS)
            else: break
        else: i += 1
    return pd.DataFrame(results)

# ==========================================
# 2. 主程序逻辑
# ==========================================
col_u1, col_u2 = st.columns([2, 1])
with col_u1:
    uploaded_files = st.file_uploader("1️⃣ 上传重复实验文件 (支持多个 CSV)", type="csv", accept_multiple_files=True)
with col_u2:
    lib_file = st.file_uploader("2️⃣ 可选：上传化合物库信息", type=["xlsx", "csv"])

if uploaded_files:
    all_reps = []
    for idx, f in enumerate(uploaded_files):
        rep_df = parse_raw_csv(f)
        rep_df['Rep'] = f"Rep_{idx+1}"
        all_reps.append(rep_df)
    
    if all_reps:
        combined_raw = pd.concat(all_reps)
        
        # 按板号和孔位进行聚合计算
        final_df = combined_raw.groupby(['Plate_ID', 'Well']).agg({
            'Ratio': ['mean', 'std', 'count'],
            'F': 'mean',
            'R': 'mean'
        }).reset_index()
        
        # 展平列名
        final_df.columns = ['Plate_ID', 'Well', 'Avg_Ratio', 'Std_Ratio', 'Count', 'Mean_F', 'Mean_R']
        
        # 计算变异系数 CV (衡量重复性)
        final_df['CV_%'] = (final_df['Std_Ratio'] / final_df['Avg_Ratio']) * 100
        
        # 计算全局 Z-score
        mu = final_df['Avg_Ratio'].mean()
        sigma = final_df['Avg_Ratio'].std()
        final_df['Z_score'] = (final_df['Avg_Ratio'] - mu) / sigma
        
        # 判定 Hit
        final_df['Is_Hit'] = final_df['Z_score'].apply(lambda x: "Yes" if x < Z_THRESHOLD else "No")

        # 关联库信息
        if lib_file:
            ldf = pd.read_excel(lib_file) if lib_file.name.endswith('xlsx') else pd.read_csv(lib_file)
            final_df = pd.merge(final_df, ldf, on=['Plate_ID', 'Well'], how='left')

        # --- 展示结果 ---
        st.divider()
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("总分析孔数", len(final_df))
        m2.metric("平均比值", f"{mu:.3f}")
        m3.metric("筛选标准差", f"{sigma:.3f}")
        m4.metric("检测到 Hits", len(final_df[final_df['Is_Hit']=='Yes']))

        # 绘图
        st.subheader("📊 筛选散点图 (Z-score 分布)")
        fig, ax = plt.subplots(figsize=(12, 5))
        hit_data = final_df[final_df['Is_Hit'] == 'Yes']
        norm_data = final_df[final_df['Is_Hit'] == 'No']
        
        ax.scatter(range(len(norm_data)), norm_data['Z_score'], c='lightgray', s=15, alpha=0.6, label='Normal')
        ax.scatter(range(len(hit_data)), hit_data['Z_score'], c='red', s=30, label='Hits')
        ax.axhline(Z_THRESHOLD, color='blue', linestyle='--', label=f'Threshold ({Z_THRESHOLD})')
        ax.set_ylabel("Z-score")
        ax.set_xlabel("Compound Index")
        ax.legend()
        st.pyplot(fig)

        # 结果表
        st.subheader("🚩 强效抑制剂清单 (Hits)")
        hits_display = final_df[final_df['Is_Hit'] == 'Yes'].sort_values('Z_score')
        st.dataframe(hits_display)

        # 下载按钮
        output = io.BytesIO()
        with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
            final_df.to_excel(writer, index=False, sheet_name='All_Results')
        st.download_button("📥 下载完整分析报告 (Excel)", output.getvalue(), "FeliEase_3.0_Report.xlsx")
