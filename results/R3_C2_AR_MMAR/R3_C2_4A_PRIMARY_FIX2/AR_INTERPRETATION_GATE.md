# AR_INTERPRETATION_GATE

UTC: 2026-09-24T09:51:45.113589+00:00

## 记号

| 量 | 符号 | 规则 |
| --- | --- | --- |
| HLA 阵元数 | **2L+1** | 仅阵元 |
| 移动距离样本数 | **N_r** | 独立符号；禁止再用 2L+1 |
| AR 阶数 | p | 打印 (2/3)(2L+1) → `PRIMARY_NOTATION_CONFLICT`；若用 N_r 规则 → **我们的实现解释**，不得反写成论文 |

## 采样

- Δt=1 s **不是** Liang 明文
- 借 Yang/Ref.7 → **`REF7_CONSISTENT_SAMPLING_ASSUMPTION`**
- 空间无折叠 control：Δr < π/k_max（见 SPATIAL_SAMPLING_ALIAS_AUDIT）

## Eq.(19)

- `PRINTED_EQ19_MODULUS_POWER_1`
- `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`
- **先** 可控复指数验证峰**位置**；**禁止**按 D(z) 深度结果择优

## l → k_r

`AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT`

待验证假设：ω_s ≡ ± k_m Δr (mod 2π)；符号/分支由 unit test + paper scale 决定，**本轮不冻结为实现公式**。

## Eq.(24)

- `PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`
- DP → **`OUR_EXACT_SOLVER_FOR_EQ24`**（本轮只写契约）
- 核心指标：**MODE_ORDER_RECOVERY_RATE**（非 Rayleigh unique）
