# R3-RC3-REANCHOR-3 observation-error boundary

UTC: 2026-09-25T10:12:48.225195+00:00

## 判定

### `PROFILED_RC3_MARGIN_ROBUST_TO_TESTED_ERROR`

ΔJ>0 比例: ε=0.5→1.00, 1.0→1.00, 2.0→1.00; median ΔJ: 2.26/1.90/1.38 dB (S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND; 非实测噪声)

## 轴与样本

ε∈[0.0, 0.25, 0.5, 1.0, 1.5, 2.0] dB · seeds [0, 1, 2, 3, 4] · 30 对 × z_true 180/200/220 · z profile 150:5:250  
标签：`TL_OBSERVATION_ERROR_SENSITIVITY_AXIS`（非海试标定）  
信号：`S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND`

## ΔJ = J_alt* − J_true*（均 profile）

| ε (dB) | median ΔJ | P(ΔJ>0) | P(ΔJ>0.5) |
| ---: | ---: | ---: | ---: |
| 0.5 | 2.261 | 1.000 | 1.000 |
| 1.0 | 1.897 | 1.000 | 0.966 |
| 2.0 | 1.378 | 1.000 | 0.828 |

A05：`RELATIVE_TL_STATIC_RANGE_BLIND_CASE`（已从机制统计剔除）

## 边界

角色不变；REANCHOR-2 的 96.7% 仅 `DECISION_THRESHOLD_SENSITIVITY` / 理想模型，非性能。
