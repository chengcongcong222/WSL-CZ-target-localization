# R3-C2.4A MMAR primary-source recovery + E-STD 准入前审计

UTC: 2026-09-24T06:44:02.063185+00:00

## 判定

### `C2_4A_PRIMARY_PDF_REQUIRED`

Liang2018 MMAR 主文 PDF 未取得（Hindawi/Wiley/DOI 均 403）。按规则停止公式恢复与编码。C2.3A 口径已修正；E-STD 模式族/Fourier 稳健性与 AR 分辨率需求已审计（不依赖 MMAR 公式）。

**请用户提供 Liang et al. 2018 PDF（DOI 10.1155/2018/7824671）。**

## C2.3A 口径冻结（不重跑）

- baseline strict：**0/20 有效**（19 EMPTY_UNIQUE_MODE；1 例 283 Hz/z=220 仅 1 unique 且 ẑ=250、margin<1 → NONINFORMATIVE）
- UNIQUE=2 / GROUP=1777：**AGGREGATED_OVER_ALL_SCANS**
- `estd_wavenumber_peaks.csv` 空：`OUTPUT_ARTIFACT_MISSING`
- EQ5 20/20：**ORACLE_EQ5_DEPTH_SIGNATURE_SURVIVES**（非可恢复性能）

冻结句：E-STD 理想深度签名存在，但 L≤2.4 km 普通 Fourier 无法可靠完成谱峰→唯一模态编号。

## 模式族稳健性（无 true zs）

见 `estd_mode_family_classification.csv`、`estd_mode_set_fourier_robustness.csv`。

**FOURIER_MODE_IDENTITY_LIMIT_ROBUST_TO_MODE_SET**（max n_rayleigh_unique=0）

## AR 需要分开的 Δk

见 `ESTD_REQUIRED_WAVENUMBER_RESOLUTION.csv`（由相邻模态 Δk 分位数换算 L_required=2π/Δk）。

## HLA 迁移（仅几何）

见 `MMAR_HLA_MIGRATION_TABLE.csv`。**HLA 不是免费**；后续须区分 MMAR-PAPER vs SINGLE_SENSOR_ABLATION。

## 停止

不编码 AR、不跑 MMAR、不进 P5。等待用户 PDF + GPT 审计。
