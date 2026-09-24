# OBSERVABLE_MODE_TRUTH_AUDIT

UTC: 2026-09-24T16:02:27.087144+00:00

## 区分

- **理论模态存在**：KRAKEN k'_m 在表中
- **可观测真实谱峰** `TRUE/OBSERVABLE MODE PEAK`：生成模型中该 mode 对 beam 的贡献足够大

## 预冻结阈值（看结果前）

`OBS_REL_THRESH = 0.03`（占该 case 最大 modal beam 贡献 |contrib(r0)| 的比例）

依据：低于主模 3% 的贡献在 5% prominence 峰检测下通常不可独立成峰，避免把弱理论模态自动当真峰。

## 检查 primary 指述：1.30–1.36 m⁻¹ 无 true spectral peak

对照 `PAPER_THEORETICAL_MODES` / 本 case `observability` 表。若模型在该窗内仍有强可观测模态，标
`PAPER_TRUE_PEAK_STRUCTURE_NOT_REPRODUCED`。
