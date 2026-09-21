# R3-C1 报告：运动/声强结构深度（Zhu 2023）

UTC：2026-09-21T07:55:14.688558+00:00

B1 已关闭：当前 CZ 仅 1 条能量可观测多径，HLA-MMAC 不成立。

## 0. 锚定

- **Zhu et al. 2023** JASA-EL：深海 HLA + 移动源；improved Fourier / matched intensity structure
- 观测为**轨迹声强结构**，非 P4.5 CZ 轮廓
- **C-S0**：`METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE`

## 1. 论文条件复现（严格口径）

峰值落在 true z **不等于**深度可定位；必须同时看 **FWHM / PSL**。

| method | mean\|err\| m | mean FWHM m | mean PSL | frac\|err\|<15m |
| --- | --- | --- | --- | --- |
| QF_fourier_like | 0.25 | 42.50 | 0.78 | 1.00 |
| matched_intensity_structure | 0.25 | 27.83 | 0.74 | 1.00 |

- 严格 paper_pass = **True**（当前实现下 ambiguity 几乎铺满深度窗 → **复现未形成有效深度定位**）

## 2. E-STD 第一 CZ（模态 I' + Zhu 匹配）

- C-S2 235 Hz, T=600 s：mean FWHM=**120.00 m**，mean PSL=**0.99**，frac FWHM<40 m=**0.00**
- 深度对 mean D=**1.16**，可分比例(D>0.05)=**0.70**
- 自匹配 err≈0 **不得**解读为深度 RMSE；结构几乎重合时 corr 峰在 true z 但无分辨力

### 与 P4.5 对照

- P4.5：z 弱/UNRESOLVED
- R3-C1：在 matched 环境、无噪、CW 轨迹声强下，CZ 深度峰仍 **宽 / 旁瓣高** → 相对 P4.5 **未形成可用的独立深度约束**

## 3. R3-C1 判定

### `C1_PAPER_REPRO_PASS_CZ_NOT_TRANSFERABLE`

Paper-method logic accepted at algorithm level, but E-STD first-CZ modal track intensity does not yield narrow unique depth peaks (mean FWHM=120.0 m, PSL=0.9927083151735492, pair mean D=1.157129456910447, frac=0.7). Self-match err≈0 must not be read as depth accuracy.

**下一步**：next RC3-C candidate in order: Yang 2015 SA beamforming → HLA modal/k-spectrum → Emmetière 2019

## 4. 停止

- 不进 C2 / P5 / RC3-B
- 若迁移失败，RC3-C 下一候选按序：Yang 2015 → HLA 模态/波数谱 → Emmetière 2019
- **R3-C1 完成后停止**
