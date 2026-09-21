# R3-A1 报告：条纹方法保真 + HLA-MFP 参考修正

UTC：2026-09-20T07:26:06.389540+00:00  ·  输出：`results/R3_A1_fidelity/`

## 0. R3-0 口径修正（冻结）

- R3-0 MFP = frozen-observation matched-field reference ambiguity (not full HLA physical ceiling)
- R3-0 RC3-A = algorithm skeleton only, method fidelity NOT yet verified
- score width != localization RMSE

## 1. RC3-A1 判定：`A1-NOT-DIRECTLY-TRANSFERABLE`

Paper repro passed but E-STD beta_eff NON_CONSTANT (median=0.8235755194946945, IQR=0.8988305024903624, median_by_f_std=0.6855537631860799). Classical constant-beta fringe ranging does not directly transfer to deep first-CZ. | HLA-MFP: single: r宽=7.75km z宽=92.5m; hla_8x2m: r宽=15.00km z宽=100.0m; hla_diag_8x10m: r宽=15.00km z宽=100.0m

**下一步**：STOP RC3-A parameter tuning; proceed to RC3-B (MMAC) next round

## 2. 论文条件复现（PE-WI）

> **R3-A1 收口勘误（文字，不重跑）**：`beta_align_search=2.2` 是条纹对齐评分中的**经验最优参数**，不作为物理 β 估计值。方法保真主要依据：`beta_local≈0.959` 对 `beta_true≈1.001`，以及距离峰值复现误差 ≈0.02 km。判定保持 `A1-NOT-DIRECTLY-TRANSFERABLE`：经典常数 β 模型不能直接迁移到 E-STD；广义/局部 β 条纹整类方法降为后备，不关闭。

| β_true | β_local | β_align | r误差 peak/count/pl | 评分峰宽 | pass |
| --- | --- | --- | --- | --- | --- |
| 1.0013 | 0.9594 | 2.2000 | 0.020 / 1.257 / -5.160 | 0.493 | β:True r:True |

## 3. E-STD β_eff(r,f)

- 判定：**NON_CONSTANT**
- median=0.8236, IQR=0.8988, P5–P95=[-1.8678, 5.0369]
|β|<1.5 占比=0.73，|β|>3 占比=0.14
|f| 方向 median 漂移 std=0.6856

**禁止**：仅扩大 β 搜索范围而不检验 β_eff(r,f) 是否常数。

> **结论**：深海第一会聚区 E-STD 场中 β_eff 随 (r,f) 明显变化 → **经典常数波导不变量模型不能直接迁移**（迁移条件不成立，不是“方法无效”）。

## 4. HLA-MFP 参考（空间孔径）

对比：单通道 / 8×2m / 大孔径诊断（8×10m，非装备要求）。S0 宽带，无噪，环境匹配，未知源按频消去（归一化空间相关）。

| 配置 | 孔径 m | r̂ km | ẑ m | r主瓣宽 km | z主瓣宽 m | corr(r,z) | 真值最优 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| single | 0.0 | 50.00 | 200.0 | 7.75 | 92.5 | -0.5692 | True |
| hla_8x2m | 14.0 | 45.00 | 167.5 | 15.00 | 100.0 | 0.0000 | False |
| hla_diag_8x10m | 70.0 | 50.25 | 250.0 | 15.00 | 100.0 | 0.0000 | False |

解读：
- R3-0 的“MFP”实为**冻结观测（单通道/多频）匹配场参考**，不是完整水平阵物理上限。
- 若大孔径 HLA 的 r/z 主瓣仍宽 → 更接近“当前场景传播信息有限”；若明显变窄 → 空间孔径有潜在增益（理论诊断，非型号指标）。
- **主瓣宽 ≠ 最终 RMSE 下限**（离网点估计误差可以远小于主瓣宽）。

## 5. 图与数据

- figures/beta_eff_ESTD.svg
- figures/HLA_MFP_single.svg / HLA_MFP_hla_8x2m.svg / HLA_MFP_hla_diag_8x10m.svg
- PAPER_REPRO.csv, beta_eff_map.csv, beta_eff_summary.csv, HLA_MFP_COMPARE.csv

## 6. 停止条件

- 未跑 RC2 困难候选 / Monte Carlo（除非 A1-PASS）
- 未进入 RC3-B/C
- 未进入 P5
- **本轮判定：A1-NOT-DIRECTLY-TRANSFERABLE**
