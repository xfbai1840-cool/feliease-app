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
# 1. 核心解析函数 (适配“左右双拼”双荧光 CSV)
# ==========================================
def parse_raw_csv(file):
    """
    根据 Discovery 实验室特定版式解析 CSV:
    左侧 F：B列(1)为溶剂, M列(12)为毒性, C-L列(2:11)为药物
    右侧 R：P列(15)为溶剂, AA列(26)为毒性, Q-Z列(16:25)为药物
    计算：Ratio = (Drug_F / Avg_Ctrl_F) / (Drug_R / Avg_Ctrl_R)
    """
    df = pd.read_csv(file, header=None)
    results = []
    
    i = 0
    while i < len(df):
        val_a = str(df.iloc[i, 0]).strip()
        
        # 如果 A 列有内容（例如 ZY-3K-1#-F），说明找到了一板的开头
        if val_a and val_a.lower() != 'nan':
            if i + int(PLATE_ROWS) <= len(df):
                # 清洗 Plate ID（去掉后缀 -F，统一叫 ZY-3K-1#）
                plate_id = val_a.replace("-F", "").strip()
                
                # --- 1. 计算这一板的溶剂对照平均值 (B列和P列) ---
                ctrl_f_vals = pd.to_numeric(df.iloc[i : i+int(PLATE_ROWS), 1], errors='coerce').dropna()
                avg_ctrl_f = ctrl_f_vals.mean() if not ctrl_f_vals.empty else 1.0
                
                ctrl_r_vals = pd.to_numeric(df.iloc[i : i+int(PLATE_ROWS), 15], errors='coerce').dropna()
                avg_ctrl_r = ctrl_r_vals.mean() if not ctrl_r_vals.empty else 1.0
                
                # --- 2. 计算这一板的毒性底线平均值 (M1~M5 和 AA1~AA5) ---
                tox_f_vals = pd.to_numeric(df.iloc[i : i+5, 12], errors='coerce').dropna()
                min_tox_f = tox_f_vals.mean() if not tox_f_vals.empty else 0
                
                tox_r_vals = pd.to_numeric(df.iloc[i : i+5, 26], errors='coerce').dropna()
                min_tox_r = tox_r_vals.mean() if not tox_r_vals.empty else 0
                
                # --- 3. 提取 C~L 列 和 Q~Z 列的药物数据 ---
                for r in range(int(PLATE_ROWS)):
                    for c_idx in range(10): # 每行 10 个药
                        col_f = 2 + c_idx  # 对应 C(2) 到 L(11)
                        col_r = 16 + c_idx # 对应 Q(16) 到 Z(25)
                        
                        raw_f = pd.to_numeric(df.iloc[i+r, col_f], errors='coerce')
                        raw_r = pd.to_numeric(df.iloc[i+r, col_r], errors='coerce')
                        
                        if pd.notna(raw_f) and pd.notna(raw_r):
                            # 计算相对值
                            rel_f = raw_f / avg_ctrl_f
                            rel_r = raw_r / avg_ctrl_r
                            
                            # 计算最终核糖体移码效率 (Ratio)
                            ratio = rel_f / rel_r if rel_r > 0 else np.nan
                            
                            # 毒性判定：只要 F 或 R 低于阳性对照平均值，即判定为有毒
                            is_toxic = "Yes" if (raw_f < min_tox_f or raw_r < min_tox_r) else "No"
                            
                            # 生成孔位 ID (例如 A02, A03)
                            well_id = f"{chr(65+r)}{c_idx+2:02d}"
                            
                            results.append({
                                "Plate_ID": plate_id,
                                "Well": well_id,
                                "Raw_F": raw_f,
                                "Raw_R": raw_r,
                                "Rel_F": rel_f,
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
