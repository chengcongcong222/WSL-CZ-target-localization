# R3-C1.1 报告：Zhu 2023 保真修正与第一 CZ 终判

UTC：2026-09-21T09:42:20.014012+00:00

**撤回** 旧判定 `C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE`（几何/Fourier/距离补偿/HLA 未按原文实现）。

## 0. 修正要点

- 底置 HLA，H≈900 m，z_s≈55 m，100 Hz，endfire 运动
- \(R_m=\sqrt{H^2+r_m^2}\)，\(\tan\theta_m=H/r_m\)，\(s_m=\sin\theta_m\)
- **\(I'=I\,R_m^2\)**（非 \(I\cdot r\)）
- **Q_F(z)**：对去均值 I'(s) 做 \(\int I'(s)e^{-j2kzs}ds\) Fourier 积分（非模板相关）
- Q_M(z)：归一化 \(|corr|\) 匹配；HLA **endfire** \(r_m=|r-x_m|\)
- 深度对距离：**D=1−|corr|**；FWHM=**全局峰邻域连续半高宽**

## 1. 论文条件复现（原文量级对照）

| method | our peak | paper ref | Δ | local FWHM | PSL | surface peak |
| --- | --- | --- | --- | --- | --- | --- |
| QF_fourier_zhu | 56.00 | 62.00 | -6.00 | 22.00 | 0.49 | 26.00 |
| QM_matched_zhu | 55.00 | 54.00 | 1.00 | 4.00 | 0.82 | 74.00 |

- paper_pass = **True**（QF=True, QM=False, 表/潜分离=True）
- 旧 `QF_fourier_like` 标记 **SUPERSEDED_SIMPLIFICATION**

## 2. E-STD 迁移（仅当 paper_pass）

| method | T | z_true | z_hat | err | local FWHM | PSL |
| --- | --- | --- | --- | --- | --- | --- |
| QF_fourier_zhu | 300 | 180 | 260.00 | 80.00 | 120.00 | 0.00 |
| QM_matched_zhu | 300 | 180 | 180.00 | 0.00 | 48.00 | 0.99 |
| QF_fourier_zhu | 300 | 200 | 140.00 | -60.00 | 120.00 | 0.00 |
| QM_matched_zhu | 300 | 200 | 200.00 | 0.00 | 68.00 | 0.99 |
| QF_fourier_zhu | 300 | 220 | 260.00 | 40.00 | 120.00 | 0.00 |
| QM_matched_zhu | 300 | 220 | 220.00 | 0.00 | 68.00 | 1.00 |
| QF_fourier_zhu | 600 | 180 | 260.00 | 80.00 | 120.00 | 0.00 |
| QM_matched_zhu | 600 | 180 | 180.00 | 0.00 | 44.00 | 0.86 |
| QF_fourier_zhu | 600 | 200 | 140.00 | -60.00 | 120.00 | 0.00 |
| QM_matched_zhu | 600 | 200 | 200.00 | 0.00 | 74.00 | 0.98 |
| QF_fourier_zhu | 600 | 220 | 260.00 | 40.00 | 120.00 | 0.00 |
| QM_matched_zhu | 600 | 220 | 220.00 | 0.00 | 72.00 | 1.00 |
| QF_fourier_zhu | 1200 | 180 | 260.00 | 80.00 | 120.00 | 0.00 |
| QM_matched_zhu | 1200 | 180 | 180.00 | 0.00 | 86.00 | 0.67 |
| QF_fourier_zhu | 1200 | 200 | 140.00 | -60.00 | 120.00 | 0.00 |
| QM_matched_zhu | 1200 | 200 | 200.00 | 0.00 | 84.00 | 0.96 |
| QF_fourier_zhu | 1200 | 220 | 260.00 | 40.00 | 120.00 | 0.00 |
| QM_matched_zhu | 1200 | 220 | 220.00 | 0.00 | 82.00 | 0.98 |

### 近邻深度 180/200/220 m 与 D_abs

- mean local FWHM=91.67 m，PSL=0.47
- pair mean D_abs=0.15，D_abs>0.05 比例=0.67
- C-S0：`METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE`

## 3. R3-C1 终判

### `C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE`

Zhu-faithful paper-condition PASSED, but first-CZ E-STD still lacks narrow depth peaks / pair separability: mean local FWHM=91.7 m, PSL=0.47, pair D_abs frac=0.67, mean D=0.148. Near-neighbour z 180/200/220 remain weakly separated. This is a CZ-scenario transfer result, not a premature method kill.

**下一步**：close Zhu route for this CZ; next RC3-C candidate in order: Yang 2015 SA beamforming (not this round)

允许终态仅：`C1_PAPER_REPRO_PASS_CZ_DEPTH_CONFIRMED` / `C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE` / `C1_PAPER_REPRO_FAIL`。

## 4. 停止

- 不进 Yang 2015 / C2 / P5；不回 RC3-B
- **R3-C1.1 完成后停止**
