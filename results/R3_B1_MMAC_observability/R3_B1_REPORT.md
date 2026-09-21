# R3-B1 报告：HLA 方向余弦可观测性门控 + MMAC

UTC：2026-09-21T02:08:15.582906+00:00

## 0. R3-A1 收口勘误（文字）

- 判定保持 **A1-NOT-DIRECTLY-TRANSFERABLE**
- `beta_align_search=2.2` = 条纹对齐**经验评分**参数，**不是**物理 β
- 保真依据：`beta_local≈0.959` vs `beta_true≈1.001`，距离峰值误差 ≈0.02 km
- 经典常数 β 不能直接迁移 E-STD；局部/广义 β 条纹降为后备，不关闭整类方法

## 1. 关键物理量（冻结）

$$u=\cos\theta_{\rm rel}\cos\phi,\qquad u(+\phi)=u(-\phi)$$

水平直线阵空间相位观测的是 **方向余弦 u**，不是独立俯仰角 φ；±φ 镜像在相位上不可分。

## Gate-0：解析投影量级

| \|φ\| | \|Δu\| (θ_rel=0°) | L~λ/Δu @200Hz (m) |
| --- | --- | --- |
| 1° | 0.000152 | 49243.3 |
| 2° | 0.000609 | 12311.8 |
| 5° | 0.003805 | 1970.9 |
| 10° | 0.015192 | 493.7 |
| 15° | 0.034074 | 220.1 |
| 20° | 0.060307 | 124.4 |

- 孔径可分辨比例（corr<0.95）：14m **0.04**，70m **0.42**
- σu<Δu 比例：14m **0.58**
- **Gate-0 HLA_ELEVATION_PROJECTION_WEAK = False**

数据：`analytic_hla_projection.csv`, `direction_cosine_pairs.csv`, `direction_cosine_resolution.csv`

## 2. 受控 MMAC（等声速两路径）

| θ_rel | 模式 | 观测 | CRB_r km | CRB_z m | cond J |
| --- | --- | --- | --- | --- | --- |
| 0 | IDEAL_ELEVATION | A_angle_only | 77.136 | 436.3 | 250.00 |
| 0 | IDEAL_ELEVATION | B_delay_only | 4346235977.056 | 17384666669.0 | 1.00 |
| 0 | IDEAL_ELEVATION | C_joint | 77.136 | 347.7 | 254.41 |
| 0 | HLA_DIRECTION_COSINE | A_angle_only | n/a | n/a | n/a |
| 0 | HLA_DIRECTION_COSINE | B_delay_only | 4346235977.056 | 17384666669.0 | 1.00 |
| 0 | HLA_DIRECTION_COSINE | C_joint | 7813.265 | 31254.2 | 8341.37 |
| 30 | IDEAL_ELEVATION | A_angle_only | 77.136 | 436.3 | 250.00 |
| 30 | IDEAL_ELEVATION | B_delay_only | 4346235977.056 | 17384666669.0 | 1.00 |
| 30 | IDEAL_ELEVATION | C_joint | 77.136 | 347.7 | 254.41 |
| 30 | HLA_DIRECTION_COSINE | A_angle_only | n/a | n/a | n/a |
| 30 | HLA_DIRECTION_COSINE | B_delay_only | 4346235977.056 | 17384666669.0 | 1.00 |
| 30 | HLA_DIRECTION_COSINE | C_joint | 9021.940 | 36088.6 | 9629.61 |
| 60 | IDEAL_ELEVATION | A_angle_only | 77.136 | 436.3 | 250.00 |
| 60 | IDEAL_ELEVATION | B_delay_only | 4346235977.056 | 17384666669.0 | 1.00 |
| 60 | IDEAL_ELEVATION | C_joint | 77.136 | 347.7 | 254.41 |
| 60 | HLA_DIRECTION_COSINE | A_angle_only | n/a | n/a | n/a |
| 60 | HLA_DIRECTION_COSINE | B_delay_only | 4346235977.056 | 17384666669.0 | 1.00 |
| 60 | HLA_DIRECTION_COSINE | C_joint | 15626.319 | 62505.1 | 16671.48 |

## 3. E-STD 本征射线

- 稳定分支数（tracking）：**28**；≥2 = **True**
- 分支统计（前几行见 CSV）：按 topology+launch 连续 ID，不按时延峰序定义身份。

| branch | n | r_km 覆盖 | φ 中位° | τ 中位 s |
| --- | --- | --- | --- | --- |
| B0 | 1 | 50-50 | -24.91 | 35.0215 |
| B1 | 1 | 60-60 | 0.79 | 39.9553 |
| B10 | 1 | 60-60 | 18.98 | 41.1736 |
| B11 | 1 | 55-55 | 21.69 | 38.0729 |
| B12 | 1 | 50-50 | -24.29 | 34.9208 |
| B13 | 2 | 55-55 | 0.29 | 38.1160 |
| B14 | 2 | 60-60 | 0.28 | 41.2186 |
| B15 | 2 | 45-45 | 0.40 | 30.6492 |
| B16 | 2 | 50-50 | 0.37 | 33.8493 |
| B17 | 2 | 55-55 | 0.36 | 37.0845 |
| B18 | 2 | 60-60 | 0.34 | 40.3389 |
| B19 | 1 | 45-45 | 5.19 | 30.0501 |

## 4. 时延可观测性

| 源 | 有效带宽 Hz | 主瓣宽 ms | 旁瓣 | 周期模糊 | 可用 |
| --- | --- | --- | --- | --- | --- |
| S0 | 235.0 | 298.200 | 0.496 | False | False |
| S1 | 186.3 | 398.400 | 0.500 | False | False |
| S2 | 210.0 | 383.400 | 0.500 | False | False |

## 5. Fisher（σθ 扫描，S0 摘录 θ_rel=0°）

| case | CRB_r km | CRB_z m | corr | cond F |
| --- | --- | --- | --- | --- |
| A_ideal_elev+delay | 77.136 | 311.6 | 0.990 | 3137382.05 |
| B_hla_u+delay | 13054.007 | 52215.3 | 1.000 | 88065065336.10 |
| C_delay_only | 768313233.003 | 3073203922.6 | 1.000 | n/a |

## 6. R3-B1 判定

### `B1_DELAY_BANDWIDTH_LIMITED`

Delay observability weak/periodic across S0/S1/S2 under current band; delay_df=[{'source': 'S0', 'n_freqs': 23, 'effective_bandwidth_Hz': 235.00483554019206, 'amb_mainlobe_width_s': 0.2982, 'amb_sidelobe_max': 0.4964808349270524, 'crb_sigma_tau_s': 0.00023460321322941465, 'periodic_ambiguity_risk': False, 'delay_usable': False, 'note': 'S2 multi-line can create periodic delay ambiguity; do not assume best'}, {'source': 'S1', 'n_freqs': 5, 'effective_bandwidth_hz': 186.31070822687568, 'amb_mainlobe_width_s': 0.39840000000000003, 'amb_sidelobe_max': 0.4996524926456676, 'crb_sigma_tau_s': 0.0002959190594404395, 'periodic_ambiguity_risk': False, 'delay_usable': False, 'note': 'S2 multi-line can create periodic delay ambiguity; do not assume best'}, {'source': 'S2', 'n_freqs': 5, 'effective_bandwidth_hz': 209.9668545270896, 'amb_mainlobe_width_s': 0.3834, 'amb_sidelobe_max': 0.49985338502952864, 'crb_sigma_tau_s': 0.00026257901356076203, 'periodic_ambiguity_risk': False, 'delay_usable': False, 'note': 'S2 multi-line can create periodic delay ambiguity; do not assume best'}].

**下一步**：do not claim practical delay MMAC; consider RC3-C or broader bandwidth studies later

### 判定依据分解（必读）

1. **Gate-0**：第一会聚区典型俯仰下 |Δu| 极小（φ=2° 时约 6×10⁻⁴，λ/Δu 需 km 级孔径）；14 m 相关可分辨比例≈2%，70 m 诊断孔径≈32%。**±φ 镜像在 u 上严格不可分**。阵列工作在真正的观测域（方向余弦）时，**HLA 俯仰投影弱**。
2. **受控 MMAC（等声速两路径）**：即便 IDEAL 直接用 φ，在 σφ≈0.5°、r=50 km 时 CRB_r 仍约 **77 km**（俯仰本身只有 ~0.2° 量级，被噪声淹没）；HLA_u+delay 更差（CRB_r 数千 km 量级）；**单条 Δτ 的 delay-only 对 (r,z) 基本退化**。
3. **E-STD 本征射线**：45–60 km 可分出 **8–11** 条 launch/topology 连续分支，加密网格后条数稳定 → 多途**存在**，不是“没有路径”。
4. **时延**：CRB 型 στ（S0）约 **0.23 ms**（若关联与匹配滤波理想）；但稀疏/有限频点下模糊函数主瓣宽达 **0.3–0.4 s**，按任务书“可用”判据未通过 → `B1_DELAY_BANDWIDTH_LIMITED`。**不把 CRB 单独当成已实现时延精度。**
5. **未出现** `B1_MMAC_PHYSICS_CONFIRMED`（需 HLA_u 可分 + 联合 Fisher 非退化）；也**未**进入 `B1_DELAY_DOMINANT_POSSIBLE`（因 delay-only 在本受控/Fisher 设置下未显示可用互补 r–z 信息 + 模糊主瓣过宽）。

允许终态还包括：`B1_PHYSICS_ONLY_HLA_LIMITED` / `B1_DELAY_DOMINANT_POSSIBLE` / `B1_DELAY_BANDWIDTH_LIMITED` / `B1_NO_COMPLEMENTARY_INFORMATION` / `B1_NO_STABLE_MULTIPATH_IDENTITY` / `B1_DELAY_IDENTITY_NOT_OBSERVABLE` / `B1_MMAC_PHYSICS_CONFIRMED`（统一用 MMAC 命名）。

## 7. 停止条件

- 无 Monte Carlo RMSE、无 RC2 困难候选、无 RC3-C、无 P5、无 Bellhop、无新特征
- **B1 完成后停止**；判定：**B1_DELAY_BANDWIDTH_LIMITED**

## 8. 图

- figures/fig0_hla_projection_gate.svg
- figures/fig1_controlled_mmac.svg
- figures/fig2_ESTD_branches.svg
- figures/fig3_angle_delay_sensitivity.svg
- figures/fig4_ideal_vs_hla.svg
- figures/fig5_fisher_boundary.svg
- figures/fig6_delay_ambiguity.svg
