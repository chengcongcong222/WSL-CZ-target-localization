# MMAR_AR_IMPLEMENTATION_CONTRACT

UTC: 2026-09-24T09:16:56.736527+00:00

## AR 输入

y[i] = B(r_i) S(r_i)  （HLA beam × Eq.8 shading）

不是 raw hydrophone / AR peak / Hankel spectrum。

## 样本与阶数

- `AR_SAMPLE_COUNT_NOTATION_CONFLICT`：印刷 2L+1 与运动序列冲突
- p often (2/3)(2L+1)
- `PAPER_RANGE_SAMPLE_INTERVAL_NOT_EXPLICIT`：Δt 未在 Liang 本节写明
  - 若用 1 s → **`REF7_CONSISTENT_SAMPLING_ASSUMPTION`**（不得写 LIANG_EXPLICIT_DT）

## 谱

- PRINTED_EQ19_MODULUS_POWER_1
- STANDARD_AR_PSD_CONTROL_MODULUS_POWER_2（对照，非原文）
- 先做 controlled complex-exp：两式 **峰位置** 是否一致（不得按深度择优）

## 峰 → k̂ → 幅度

- AR 峰位 → k̂_m
- 幅度 **只** 从 Hankel g(k̂_m)
- 禁止 AR peak height → amplitude

## l → k_r

`AR_NORMALIZED_SPECTRAL_AXIS_MAPPING_NOT_EXPLICIT`

允许待验证假设：l = ± k_r Δr (mod 2π)，须 unit test + paper scale 决定符号/尺度。

## Eq.24 求解

`PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`；DP → `OUR_EXACT_SOLVER_FOR_EQ24`。

## b_m

`FROM_EQ9_DERIVED_BM` + BF_m（Eq.4）。

## δ

`ORACLE_OFFSET_ALIGNMENT`（δ=0）。
