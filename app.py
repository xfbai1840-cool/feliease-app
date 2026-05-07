import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import io

# ==========================================
# 0. 网页配置
# ==========================================
st.set_page_config(page_title="FeliEase 双荧光高通量筛选", layout="wide")

st.title("🧪 FeliEase 筛选分析平台 - 双荧光重复实验版")
st.markdown("""
**更新说明：**
- **双荧光归一化**：自动处理 Firefly(F) 和 Renilla(R) 的比值。
- **3次重复合并**：自动匹配不同 Sheet 或文件的孔位。
- **Z-score 统计**：计算 $Z = \frac{x - \mu}{\sigma}$，评估化合物抑制强度。
""")

# 侧边栏：参数配置
with st.sidebar:
    st.header("⚙️ 实验参数")
    PLATE_HEIGHT = st.number_input("每板行数 (Data Rows)", value=8)
    
    st.divider()
    THRESHOLD_Z = st.number_input("Z-score Hit 阈值 (通常 < -2 或 -3)", value=-2.0)
    
    st.info("默认设置：\n1. 自动根据 Plate_ID 和 Physical_Well 匹配3次重复。\n2. 计算公式：Normalized = F/R; Z = (Ratio - Mean)/Std")

# ==========================================
# 1. 核心计算工具
# ==========================================
def process_single_file(df, plate_height):
    """处理单个原始文件，提取 F 和 R 并归一化"""
    results = []
    # 这里的索引逻辑需根据您的 CSV 布局微调
    # 假设：A列ID, B列内参F, C-L列药物F, M列毒性F... 
    # 此处为简化，假设上传的是清洗后的干净长表，若仍是原始排布，建议使用该函数进行解析
    return df

# ==========================================
# 2. 数据处理逻辑
# ==========================================
st.subheader("1️⃣ 上传数据")
col_u1, col_u2 = st.columns(2)

with col_u1:
    files = st.file_uploader("上传3个重复的实验结果 (CSV/XLSX)", type=["csv", "xlsx"], accept_multiple_files=True)
with col_u2:
    lib_file = st.file_uploader("上传化合物库信息表", type=["csv", "xlsx"])

if len(files) >= 1:
    dfs = []
    for f in files:
        if f.name.endswith('.csv'):
            tmp_df = pd.read_csv(f)
        else:
            tmp_df = pd.read_excel(f)
        dfs.append(tmp_df)
    
    # 假设每个文件已有 Physical_Plate, Physical_Well, Firefly, Renilla 四列
    # 如果没有，此处需要运行您之前的“自动找板”逻辑进行转换
    
    if len(dfs) > 0:
        # 合并重复
        try:
            # 1. 基础合并
            main_df = dfs[0].copy()
            # 关键：计算单次实验的归一化比值
            for i, d in enumerate(dfs):
                # 计算 F/R
                d[f'Ratio_Rep{i+1}'] = d['Firefly'] / d['Renilla']
            
            # 2. 纵向合并并计算均值
            combined = pd.concat(dfs)
            # 按孔位聚合
            final_res = combined.groupby(['Physical_Plate', 'Physical_Well']).agg({
                'Firefly': 'mean',
                'Renilla': 'mean'
            }).reset_index()
            
            # 计算平均 Ratio
            final_res['Avg_Ratio'] = final_res['Firefly'] / final_res['Renilla']
            
            # 3. 计算 Z-score (基于所有药物孔的 Avg_Ratio)
            mu = final_res['Avg_Ratio'].mean()
            sigma = final_res['Avg_Ratio'].std()
            final_res['Z_score'] = (final_res['Avg_Ratio'] - mu) / sigma
            
            # 4. 计算 3 次重复的 CV (变异系数) - 衡量实验质量
            # 这里需要把 Rep1, Rep2, Rep3 横向拼起来算
            
            # 5. 判定 Hit
            final_res['Is_Hit'] = final_res['Z_score'].apply(lambda x: "Yes" if x < THRESHOLD_Z else "No")
            
            # ==========================
            # 关联库信息
            # ==========================
            if lib_file:
                lib_df = pd.read_excel(lib_file) if lib_file.name.endswith('.xlsx') else pd.read_csv(lib_file)
                final_res = pd.merge(final_res, lib_df, on=['Physical_Plate', 'Physical_Well'], how='left')

            # ==========================
            # 可视化
            # ==========================
            st.divider()
            c1, c2, c3 = st.columns(3)
            c1.metric("平均比值 (Mean Ratio)", f"{mu:.4f}")
            c2.metric("标准差 (Std)", f"{sigma:.4f}")
            c3.metric("发现 Hits", len(final_res[final_res['Is_Hit'] == "Yes"]))

            st.subheader("📊 筛选散点图 (Z-score 模式)")
            fig, ax = plt.subplots(figsize=(10, 5))
            colors = ['red' if x == "Yes" else 'lightgray' for x in final_res['Is_Hit']]
            ax.scatter(range(len(final_res)), final_res['Z_score'], c=colors, s=10, alpha=0.5)
            ax.axhline(THRESHOLD_Z, color='blue', linestyle='--', label=f'Threshold (Z={THRESHOLD_Z})')
            ax.set_ylabel("Z-score")
            ax.set_xlabel("Compound Index")
            ax.legend()
            st.pyplot(fig)

            # 展示结果
            st.subheader("📋 筛选结果详表")
            st.dataframe(final_res.sort_values('Z_score').head(50))
            
            # 下载
            output = io.BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                final_res.to_excel(writer, index=False)
            st.download_button("📥 下载全量分析结果", output.getvalue(), "Screening_Zscore_Final.xlsx")

        except Exception as e:
            st.error(f"数据处理失败，请确保文件包含 'Firefly' 和 'Renilla' 列。错误详情: {e}")