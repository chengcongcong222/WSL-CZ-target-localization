# MMAR_AR_ORDER_AUDIT

UTC: 2026-09-24T09:03:08.658004+00:00

## AR_SAMPLE_COUNT_NOTATION_CONFLICT

| 证据 | 内容 |
| --- | --- |
| 正文 | HLA 有 **2L+1** receivers |
| Eq.(17) | y[i], i=1…**2L+1** |
| Eq.(18) | p often **(2/3)(2L+1)** |
| 仿真 | **L=5 → 11 elements** |

字面解释：N_data=11, p≈7。与 4990/1990 m 运动孔径谱估计物理不一致。

| 解释 | 状态 |
| --- | --- |
| PRINTED_LITERAL: N_data=11, p≈7 | 记录 |
| PHYSICALLY_REQUIRED_MOVING_SEQUENCE: 由 range samples 组成，**论文无独立样本数符号** | 记录 |
| 本轮选定 | **禁止决定** |

## AR estimator

- **modified covariance**（避免 spectral line splitting）
- 峰位置 → k̂_m；**峰高不得作幅度**（large variance）

## AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT

P_AR(l) 的 l→物理 k_r 映射**原文未明确**。不得假定 k=ω/c 或 k=l/Δr。
