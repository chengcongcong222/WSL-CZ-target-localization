# P2 RC2 报告：连续方位 / 平台运动下的运动学候选边界

项目：**B3D-3** 深海第一会聚区单平台拖曳阵目标定位  
工作包：**P2 RC2 kinematic boundary**  
UTC：2026-09-18T10:45:14.453348+00:00

## 0. 停止条件自检

| 条件 | 状态 |
| --- | --- |
| 观测模型正确（完整二维 + atan2） | 是 |
| 真值不泄露到估计器 | 是 |
| 网格 / off-grid 检查 | 是 |
| 参数扫描完成 | 是 |
| 理论边界形成 | 是 |
| 因不可辨识追加传感器/声学特征 | 否 |
| 自动进入 RC3 | 否 |

## 1. 状态与观测

```
x = [r0, θ0, v, ψ]^T   (z=200 m 固定)
平台: (0,0)→航向0°, U=2 m/s; 单次转角在 t=T/2
目标: xt=r0 cosθ0 + v t cosψ, yt=r0 sinθ0 + v t sinψ
观测: θ_obs = atan2(yt−yp, xt−xp) + ν,  ν~N(0,σθ²)
```

候选盒：r0∈[45,60] km, θ0∈[−5°,5°], v∈[1,3] m/s, ψ∈[−15°,15°]。
粗网格：r 1 km × θ 0.5° × v 0.2 m/s × ψ 1°（含 0°），节点数 **114576**。
真值不保证在网格上；估计=粗搜索+局部细化；off-grid 三组对照。

### 分析场景

| 场景 | 真值 | 含义 |
| --- | --- | --- |
| **LIN** | r0=50 km, θ0=0.0°, v=2.0 m/s, ψ=0.0° | 方位恒为0，r–v本征不可辨识；RC2理论下界 |
| **CR5** | r0=50 km, θ0=0.0°, v=2.0 m/s, ψ=5.0° | 存在方位变化率，RC2可部分约束 |
| **CR10** | r0=50 km, θ0=0.0°, v=2.0 m/s, ψ=10.0° | 方位变化率更强 |

主性能曲线用 **CR5**；LIN 作为 RC2 极限边界与甲方可解释的“不可辨识”证据。

## 2. 方法要点

1. 真实生成含噪方位序列，θ0 未知。
2. `cost=Σ wrap(θ_pred−θ_obs)²`；接受集 `cost≤c_min+13.3σθ²`。
3. **主指标**：接受集在 (r0,θ0,v,ψ) 上的宽度、收缩率、真值是否落入接受带。
4. 点估计：top 网格局部细化；若接受集呈宽脊，则用接受集中位数作为诚实点估计。
5. FIM/`HᵀR⁻¹H`/CRLB 在真值处计算，仅用于可观测性，不喂给估计器。
6. 非零 σθ 每格 N_TRIALS 次随机试验；σθ=0 为理想上限。

## 3. 四个核心问题

### Q1 距离–速度–航向何时不可区分？

**定理性事实（本包定量验证）**：当目标与平台相对速度沿视线、且方位变化率为 0 时，
`θ(t)` 与 `(r0,v)` 无关，FIM 秩亏，候选沿 r–v 脊线分布。
冻结基准 LIN 正落入该情形。

- LIN σθ=0.1° T=600 s：接受集 r 宽度 **15.000 km**，后验σr **4.592 km**，v 宽度 **2.000 m/s**，**CRLB(r)：未定义/无界（FIM秩亏）**，cond **∞**，span **0.000°**。
- 注：秩亏时伪逆数值可能给出 CRLB≈0，**不代表距离误差下界为 0**；对外一律记为 UNDEFINED/UNBOUNDED。
- 直航 + 近共线时 ψ 与 r0/v 强耦合；|ψ| 增大或平台转向后，方位 span>0，信息量上升，但在冻结的 σθ/T/转角范围内距离往往仍只有数 km 级不确定度。

### Q2 方位精度需要达到多少？

见 `bearing_accuracy_time_boundary.csv` 与图1（CR5/LIN）。
- CR5 σθ=0.1° T=600 s：后验σr **4.589 km**，r宽度 **15.000 km**，|Δr| **2.861 km**，span **0.120°**，r可辨识比例 **0.00**。
- CR5 σθ=0.5° T=300 s：后验σr **4.604 km**。
- 读图：给定允许距离误差，在 σθ×T 图上找满足 posterior σ_r 或 r 可辨识门限的组合。
- **关键结论**：在 E-STD 冻结量级（r≈50 km, U=v≈2 m/s, T≤1200 s, 转角≤15°）下，即使 σθ=0.01°–0.05°，弱交叉 CR5 的方位总变化仅约 0.01°–0.24°，距离信息量有限；共线 LIN 则对任意 σθ 都不能唯一恢复 r–v。
- 这说明：**纯 RC2 + 本包冻结的平台机动能力，不足以在第一会聚区距离上形成高精度定位**；更优 σθ 只能把脊线“照清楚”，不能消除几何不可观测。更高性能需要更强机动/更长观测/传播增量（RC3）——但 P3 只允许在困难候选上验证传播增量，不回退审计。

### Q3 小幅转向是否真正有效？

见 `maneuver_boundary.csv` 与图3。
- CR5 σθ=0.1° T=600 s：转角 0°→15°，后验σr **4.588→4.418 km**，r宽度 **15→15 km**，r可辨识比例 **0.00→0.00**。
- **口径修正（对外统一）**：小幅转向**改变了局部观测几何**，但在 r≈50 km、T≤1200 s、转角≤15° 条件下，**不足以形成实用的距离约束**。不写“机动已恢复距离可观测性”。sv_min 在直航/近共线时为机器极小量或 0（打印 0.0000），无需纠结“严格满秩 vs 极弱非零”。
- 更长 T + 更好 σθ + 转角时改善更明显（如 T=1200 s、σθ=0.05°、转角15°：后验σr≈2.93 km），仍属 km 级不确定。
- 本轮不评估拖曳阵孔径/波束（P4）。

### Q4 哪些困难候选交给 RC3？

交付 `hard_candidate_pairs.csv` 与 `hard_candidate_pairs_rc3_handoff.csv`。机制：

1. **r–v 补偿**：更远更快 vs 更近更慢，方位接近（LIN 下整条脊线皆是）；
2. **航向近镜像**：±ψ 差在有限 span 下不足；
3. **θ0 与 ψ 联合补偿**。

RC3 **只允许**对这些困难候选对比较 “RC2 only” vs “RC2+传播”，回答传播是否提供独立增量。

## 4. 关键结果摘录

### 4.1 σθ×T（场景对比）

| 场景 | T | σθ° | 后验σr km | r宽度 km | \|Δr\| km | contr_r | 可辨识r | CRLB r | span° |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CR5 | 300 | 0.02 | 4.586 | 15.000 | 2.910 | 0.000 | 0.00 | 0.000 | 0.060 |
| CR5 | 300 | 0.05 | 4.589 | 15.000 | 2.859 | 0.000 | 0.00 | 0.000 | 0.060 |
| CR5 | 300 | 0.1 | 4.593 | 15.000 | 2.790 | 0.000 | 0.00 | 0.000 | 0.060 |
| CR5 | 300 | 0.2 | 4.596 | 15.000 | 2.696 | 0.000 | 0.00 | 0.000 | 0.060 |
| CR5 | 300 | 0.5 | 4.604 | 15.000 | 2.606 | 0.000 | 0.00 | 0.000 | 0.060 |
| CR5 | 600 | 0.02 | 4.638 | 15.000 | 2.962 | 0.000 | 0.00 | 0.000 | 0.120 |
| CR5 | 600 | 0.05 | 4.599 | 15.000 | 2.907 | 0.000 | 0.00 | 0.000 | 0.120 |
| CR5 | 600 | 0.1 | 4.589 | 15.000 | 2.861 | 0.000 | 0.00 | 0.000 | 0.120 |
| CR5 | 600 | 0.2 | 4.592 | 15.000 | 2.813 | 0.000 | 0.00 | 0.000 | 0.120 |
| CR5 | 600 | 0.5 | 4.595 | 15.000 | 2.738 | 0.000 | 0.00 | 0.000 | 0.120 |
| CR5 | 1200 | 0.02 | 4.528 | 15.000 | 3.120 | 0.000 | 0.00 | 0.000 | 0.240 |
| CR5 | 1200 | 0.05 | 4.591 | 15.000 | 3.027 | 0.000 | 0.00 | 0.000 | 0.240 |
| CR5 | 1200 | 0.1 | 4.611 | 15.000 | 2.933 | 0.000 | 0.00 | 0.000 | 0.240 |
| CR5 | 1200 | 0.2 | 4.587 | 15.000 | 2.896 | 0.000 | 0.00 | 0.000 | 0.240 |
| CR5 | 1200 | 0.5 | 4.591 | 15.000 | 2.768 | 0.000 | 0.00 | 0.000 | 0.240 |
| LIN | 300 | 0.02 | 4.590 | 15.000 | 2.922 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 300 | 0.05 | 4.592 | 15.000 | 2.906 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 300 | 0.1 | 4.592 | 15.000 | 2.891 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 300 | 0.2 | 4.596 | 15.000 | 2.745 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 300 | 0.5 | 4.605 | 15.000 | 2.576 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 600 | 0.02 | 4.607 | 15.000 | 2.700 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 600 | 0.05 | 4.593 | 15.000 | 2.901 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 600 | 0.1 | 4.592 | 15.000 | 2.906 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 600 | 0.2 | 4.592 | 15.000 | 2.905 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 600 | 0.5 | 4.595 | 15.000 | 2.801 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 1200 | 0.02 | 4.610 | 15.000 | 2.500 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 1200 | 0.05 | 4.607 | 15.000 | 2.659 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 1200 | 0.1 | 4.598 | 15.000 | 2.854 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 1200 | 0.2 | 4.592 | 15.000 | 2.904 | 0.000 | 0.00 | 0.000 | 0.000 |
| LIN | 1200 | 0.5 | 4.592 | 15.000 | 2.906 | 0.000 | 0.00 | 0.000 | 0.000 |

### 4.2 航向可辨识（σθ=0.1°）

| T | ψ° | \|Δψ\|° | 可辨识率 | CRLB ψ° | span° | r宽度 km |
| --- | --- | --- | --- | --- | --- | --- |
| 300 | -15 | 3.476 | 1.00 | 1551.768 | 0.178 | 15.000 |
| 300 | -5 | 2.566 | 0.75 | 1552.328 | 0.060 | 15.000 |
| 300 | +0 | 2.988 | 0.38 | 5.020 | 0.000 | 15.000 |
| 300 | +5 | 2.521 | 0.62 | 1552.328 | 0.060 | 15.000 |
| 300 | +10 | 0.860 | 1.00 | 1552.117 | 0.119 | 15.000 |
| 300 | +15 | 2.772 | 1.00 | 1551.768 | 0.178 | 15.000 |
| 600 | -15 | 2.263 | 1.00 | 281.574 | 0.356 | 15.000 |
| 600 | -5 | 1.831 | 0.88 | 281.776 | 0.120 | 15.000 |
| 600 | +0 | 1.083 | 0.88 | 1.818 | 0.000 | 15.000 |
| 600 | +5 | 1.652 | 1.00 | 281.776 | 0.120 | 15.000 |
| 600 | +10 | 0.760 | 1.00 | 281.700 | 0.239 | 15.000 |
| 600 | +15 | 2.318 | 1.00 | 281.574 | 0.356 | 15.000 |
| 1200 | -15 | 2.271 | 1.00 | 49.501 | 0.713 | 15.000 |
| 1200 | -5 | 1.501 | 1.00 | 49.571 | 0.240 | 15.000 |
| 1200 | +0 | 0.291 | 1.00 | 0.651 | 0.000 | 15.000 |
| 1200 | +5 | 1.610 | 1.00 | 49.571 | 0.240 | 15.000 |
| 1200 | +10 | 0.507 | 1.00 | 49.544 | 0.478 | 15.000 |
| 1200 | +15 | 2.344 | 1.00 | 49.501 | 0.713 | 15.000 |

### 4.3 机动边界（σθ=0.1°）

| 场景 | 转角° | T | 后验σr km | r宽度 km | \|Δr\| km | sv_min | span° | 可辨识r |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CR5 | 0 | 300 | 4.592 | 15.000 | 2.861 | 0.0000 | 0.060 | 0.00 |
| CR5 | 0 | 600 | 4.588 | 15.000 | 2.882 | 0.0000 | 0.120 | 0.00 |
| CR5 | 0 | 1200 | 4.599 | 15.000 | 2.949 | 0.0000 | 0.240 | 0.00 |
| CR5 | 5 | 300 | 4.592 | 15.000 | 2.760 | 0.0000 | 0.030 | 0.00 |
| CR5 | 5 | 600 | 4.586 | 15.000 | 2.834 | 0.0000 | 0.060 | 0.00 |
| CR5 | 5 | 1200 | 4.485 | 15.000 | 2.866 | 0.0000 | 0.120 | 0.00 |
| CR5 | 15 | 300 | 4.601 | 15.000 | 2.519 | 0.0000 | 0.059 | 0.00 |
| CR5 | 15 | 600 | 4.418 | 15.000 | 3.603 | 0.0000 | 0.118 | 0.00 |
| CR5 | 15 | 1200 | 4.035 | 15.000 | 2.694 | 0.0000 | 0.236 | 0.00 |
| LIN | 0 | 300 | 4.592 | 15.000 | 2.882 | 0.0000 | 0.000 | 0.00 |
| LIN | 0 | 600 | 4.592 | 15.000 | 2.906 | 0.0000 | 0.000 | 0.00 |
| LIN | 0 | 1200 | 4.599 | 15.000 | 2.845 | 0.0000 | 0.000 | 0.00 |
| LIN | 5 | 300 | 4.591 | 15.000 | 2.847 | 0.0000 | 0.030 | 0.00 |
| LIN | 5 | 600 | 4.567 | 15.000 | 3.046 | 0.0000 | 0.060 | 0.00 |
| LIN | 5 | 1200 | 4.568 | 15.000 | 2.261 | 0.0000 | 0.120 | 0.00 |
| LIN | 15 | 300 | 4.548 | 15.000 | 3.198 | 0.0000 | 0.089 | 0.00 |
| LIN | 15 | 600 | 4.521 | 15.000 | 2.495 | 0.0000 | 0.178 | 0.00 |
| LIN | 15 | 1200 | 3.721 | 15.000 | 3.293 | 0.0000 | 0.356 | 0.00 |

### 4.4 off-grid 摘要

| OG | T | σθ | \|Δr\| km | r宽度 km | 真值在r带 | 真值在网格? |
| --- | --- | --- | --- | --- | --- | --- |
| OG1 | 300 | 0.05 | 0.137 | 15.000 | 1.00 | False |
| OG1 | 300 | 0.1 | 0.207 | 15.000 | 1.00 | False |
| OG1 | 300 | 0.2 | 0.262 | 15.000 | 1.00 | False |
| OG1 | 600 | 0.05 | 0.610 | 15.000 | 1.00 | False |
| OG1 | 600 | 0.1 | 0.609 | 15.000 | 1.00 | False |
| OG1 | 600 | 0.2 | 0.569 | 15.000 | 1.00 | False |
| OG1 | 1200 | 0.05 | 0.371 | 15.000 | 1.00 | False |
| OG1 | 1200 | 0.1 | 0.581 | 15.000 | 1.00 | False |
| OG1 | 1200 | 0.2 | 0.604 | 15.000 | 1.00 | False |
| OG2 | 300 | 0.05 | 0.819 | 2.125 | 0.12 | False |
| OG2 | 300 | 0.1 | 1.399 | 10.625 | 1.00 | False |
| OG2 | 300 | 0.2 | 2.243 | 15.000 | 1.00 | False |
| OG2 | 600 | 0.05 | 1.851 | 14.125 | 1.00 | False |
| OG2 | 600 | 0.1 | 2.162 | 14.875 | 1.00 | False |
| OG2 | 600 | 0.2 | 2.366 | 15.000 | 1.00 | False |
| OG2 | 1200 | 0.05 | 3.783 | 15.000 | 1.00 | False |
| OG2 | 1200 | 0.1 | 4.653 | 15.000 | 1.00 | False |
| OG2 | 1200 | 0.2 | 4.723 | 15.000 | 1.00 | False |
| OG3 | 300 | 0.05 | 8.144 | 15.000 | 1.00 | False |
| OG3 | 300 | 0.1 | 7.771 | 15.000 | 1.00 | False |
| OG3 | 300 | 0.2 | 6.768 | 15.000 | 1.00 | False |
| OG3 | 600 | 0.05 | 6.543 | 15.000 | 1.00 | False |
| OG3 | 600 | 0.1 | 6.422 | 15.000 | 1.00 | False |
| OG3 | 600 | 0.2 | 6.574 | 15.000 | 1.00 | False |
| OG3 | 1200 | 0.05 | 5.869 | 15.000 | 1.00 | False |
| OG3 | 1200 | 0.1 | 6.049 | 15.000 | 1.00 | False |
| OG3 | 1200 | 0.2 | 6.139 | 15.000 | 1.00 | False |

三组 off-grid 真值均不在粗网格节点上；覆盖指标与 on-grid 场景同量级，排除“真值写入网格”伪影。

## 5. 图目录

| 文件 | 内容 |
| --- | --- |
| figures/fig1a_sigma_T_range_error.svg |  |
| figures/fig1a_sigma_T_range_poststd_cr5.svg |  |
| figures/fig1a_sigma_T_range_poststd_lin.svg |  |
| figures/fig1a_sigma_T_range_width_cr5.svg |  |
| figures/fig1a_sigma_T_range_width_lin.svg |  |
| figures/fig1b_sigma_T_range_err_cr5.svg |  |
| figures/fig1b_sigma_T_range_err_lin.svg |  |
| figures/fig1b_sigma_T_speed_error.svg |  |
| figures/fig1c_sigma_T_contraction_heatmap.svg |  |
| figures/fig1c_sigma_T_speed_err_cr5.svg |  |
| figures/fig1c_sigma_T_speed_err_lin.svg |  |
| figures/fig1d_sigma_T_contraction_cr5.svg |  |
| figures/fig1d_sigma_T_contraction_lin.svg |  |
| figures/fig1d_sigma_T_psi_error.svg |  |
| figures/fig1e_sigma_T_truth_in_band_cr5.svg |  |
| figures/fig1e_sigma_T_truth_in_band_lin.svg |  |
| figures/fig2a_heading_err_vs_T.svg |  |
| figures/fig2b_heading_identifiable_rate_heatmap.svg |  |
| figures/fig2b_psi_crlb_vs_T.svg |  |
| figures/fig2c_heading_identifiable_rate_heatmap.svg |  |
| figures/fig2c_psi_crlb_heatmap.svg |  |
| figures/fig2d_psi_crlb_heatmap.svg |  |
| figures/fig3a_turn_T_range_error.svg |  |
| figures/fig3a_turn_T_range_poststd_cr5.svg |  |
| figures/fig3a_turn_T_range_poststd_lin.svg |  |
| figures/fig3a_turn_T_range_width_cr5.svg |  |
| figures/fig3a_turn_T_range_width_lin.svg |  |
| figures/fig3b_turn_T_svmin.svg |  |
| figures/fig3b_turn_T_svmin_cr5.svg |  |
| figures/fig3b_turn_T_svmin_lin.svg |  |
| figures/fig3c_turn_T_contraction.svg |  |
| figures/fig3c_turn_T_contraction_cr5.svg |  |
| figures/fig3c_turn_T_contraction_lin.svg |  |
| figures/fig3d_turn_T_range_heatmap.svg |  |
| figures/fig3d_turn_T_range_width_heatmap_cr5.svg |  |
| figures/fig3d_turn_T_range_width_heatmap_lin.svg |  |
| figures/fig3e_hard_pair_bearing_rmse.svg |  |
| figures/fig3e_turn_T_bearing_span_cr5.svg |  |
| figures/fig3e_turn_T_bearing_span_lin.svg |  |
| figures/fig3f_hard_pair_bearing_rmse.svg |  |

## 6. 局限

- 不研究声传播、谱线、阵列方向图。
- z 固定，不讨论深度可观测。
- “不可辨识”不触发新传感器，也不自动进入 RC3。
- E-STD 是理论算例，不是实测海区复现。

## 7. 交付清单

- `results/P2_RC2_kinematic_boundary/bearing_accuracy_time_boundary.csv`
- `results/P2_RC2_kinematic_boundary/figures\fig1a_sigma_T_range_error.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1a_sigma_T_range_poststd_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1a_sigma_T_range_poststd_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1a_sigma_T_range_width_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1a_sigma_T_range_width_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1b_sigma_T_range_err_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1b_sigma_T_range_err_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1b_sigma_T_speed_error.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1c_sigma_T_contraction_heatmap.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1c_sigma_T_speed_err_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1c_sigma_T_speed_err_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1d_sigma_T_contraction_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1d_sigma_T_contraction_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1d_sigma_T_psi_error.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1e_sigma_T_truth_in_band_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig1e_sigma_T_truth_in_band_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig2a_heading_err_vs_T.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig2b_heading_identifiable_rate_heatmap.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig2b_psi_crlb_vs_T.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig2c_heading_identifiable_rate_heatmap.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig2c_psi_crlb_heatmap.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig2d_psi_crlb_heatmap.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3a_turn_T_range_error.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3a_turn_T_range_poststd_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3a_turn_T_range_poststd_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3a_turn_T_range_width_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3a_turn_T_range_width_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3b_turn_T_svmin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3b_turn_T_svmin_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3b_turn_T_svmin_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3c_turn_T_contraction.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3c_turn_T_contraction_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3c_turn_T_contraction_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3d_turn_T_range_heatmap.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3d_turn_T_range_width_heatmap_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3d_turn_T_range_width_heatmap_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3e_hard_pair_bearing_rmse.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3e_turn_T_bearing_span_cr5.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3e_turn_T_bearing_span_lin.svg`
- `results/P2_RC2_kinematic_boundary/figures\fig3f_hard_pair_bearing_rmse.svg`
- `results/P2_RC2_kinematic_boundary/hard_candidate_pairs.csv`
- `results/P2_RC2_kinematic_boundary/hard_candidate_pairs_all.csv`
- `results/P2_RC2_kinematic_boundary/hard_candidate_pairs_rc3_handoff.csv`
- `results/P2_RC2_kinematic_boundary/heading_identifiability_boundary.csv`
- `results/P2_RC2_kinematic_boundary/maneuver_boundary.csv`
- `results/P2_RC2_kinematic_boundary/observability_scan.csv`
- `results/P2_RC2_kinematic_boundary/offgrid_trials.csv`
- `results/P2_RC2_kinematic_boundary/P2_RC2_CONFIG.json`
- `results/P2_RC2_kinematic_boundary/P2_RC2_GPT_SYNC.md`
- `results/P2_RC2_kinematic_boundary/P2_RC2_REPORT.md`
- `results/P2_RC2_kinematic_boundary/p2_results.json`

---

**P2 完成后停止。** 下一轮仅基于困难候选对判断 RC3 增量。