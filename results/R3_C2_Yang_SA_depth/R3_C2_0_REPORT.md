# R3-C2.0 报告：Yang 合成孔径模态可观测性预关口

UTC：2026-09-22T05:17:48.914782+00:00

## 0. C1 冻结

- Zhu-QF：CLOSED；Zhu-QM：`C1_QM_MULTIFREQ_ONLY_CONDITIONAL`（20 m 候选信号 ≠ 定深精度）
- HLA 无额外深度自由度；不再 C1.x

## 1. 文献恢复

- Yang 2015 DOI 10.1121/1.4929748；Erratum 2018 DOI 10.1121/1.5081712
- **网络未取回 Erratum 内容 → `C2_0_PAPER_RECOVERY_BLOCKED`**，不臆造 Yang 公式
- f0 条件：KNOWN_F0 / KNOWN_TRUE_RANGE_INCREMENT

## 2. 孔径定义

**L = 径向相对距离变化 |Δr|**，不是航程。

| L km | v_r(T=600) m/s | v_r(T=1200) m/s |
| --- | --- | --- |
| 0.15 | 0.2500 | 0.1250 |
| 0.30 | 0.5000 | 0.2500 |
| 0.60 | 1.0000 | 0.5000 |
| 1.20 | 2.0000 | 1.0000 |
| 2.40 | 4.0000 | 2.0000 |

## 3. 模态与 η（z_s=200, thr=0.1, L=2.4 km）

| f Hz | n_eff | median η | n 可分对 | D(200/220) | D(200/210) |
| --- | --- | --- | --- | --- | --- |
| 201 | 12 | 0.010 | 0 | 0.0384 | 0.0098 |
| 235 | 12 | 0.009 | 0 | 0.0384 | 0.0098 |
| 283 | 12 | 0.007 | 0 | 0.0384 | 0.0098 |
| 338 | 12 | 0.006 | 0 | 0.0384 | 0.0098 |

## 4. C1 的 338 Hz 现象

 f_hz  C1_corr_200_220  C1_D_abs_implied  n_eff_modes_A0p1  n_eff_modes_A0p2  n_resolvable_eta_L2p4km  median_eta_L2p4km  mode_D_200_220  mode_D_200_210
201.0            0.579             0.421                12                 9                        0           0.010104        0.038384        0.009785
235.0            0.994             0.006                12                 9                        0           0.008641        0.038384        0.009785
283.0            1.000             0.000                12                 9                        0           0.007175        0.038384        0.009785
338.0            0.067             0.933                12                 9                        0           0.006007        0.038384        0.009785

not fully explained — keep as numerical cue only

## 5. 判定

### `C2_0_PAPER_RECOVERY_BLOCKED`

Yang 2015 full text + 2018 Erratum (DOI 10.1121/1.5081712) NOT recovered (PubMed/AIP/DOI network blocked). Per protocol, no Yang beamforming formulas invented. Task-specified modal/aperture physics pre-gate computed below as supporting numbers only.

**下一步**：recover Yang2015+Erratum text before any SA depth estimator; then re-open C2 gate

## 6. 停止

- 不做完整 Yang 深度估计 / Doppler / f0 误差 / MC / P5


## 7. 重要物理限制（预计算）

- 在当前 E-STD ModeModel 中 **φ_m(z) 与频率近似无关**，故四频的模态深度向量 D 相同（≈0.038）——**不能**用它解释 C1 的 338 Hz 选择性。
- 相邻有效模态 η@L=2.4 km ≈ **0.006–0.01 ≪ 1**：**径向孔径 2.4 km 仍远不足以分开模态波数**（Δk_SA=2π/L）。要 η~1 需要 L~|Δk|^{-1} 量级的大孔径。
- 深度签名 D(200/220)≈0.04、D(200/210)≈0.01 → 即便模态可分，对 200 m 近邻深度也弱。

**C2_0_PAPER_RECOVERY_BLOCKED** 为本轮正式判定；上述为任务书公式的支撑量级，不是 Yang 原文复现。
