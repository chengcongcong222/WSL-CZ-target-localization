# BAND_EDGE_SCORING_AUDIT

UTC: 2026-09-24T11:51:27.216561+00:00
基线 commit: 1663d16adc2758fc8a56718ed14488d782613acd

## 问题

K_BAND=[1.32,1.50] 用 ±1e-9 硬截断；离散谱网格 Δω=2π/16384 映射到 k 的半格误差

$$
\epsilon_k=\frac{\Delta\omega}{2\Delta r}+10^{-12}
$$

会把 1.319956（k_true=1.32）等 **网格量化** 结果误判为带外 → 0 recall。

## 修正（仅评分）

合法分支：

$$
k_{min}-\epsilon_k \le k_{cand} \le k_{max}+\epsilon_k
$$

- ε_k **只由** Δω 与 Δr 决定（与真值无关）
- **不**把 k_cand 截回 1.32/1.50
- AR 系数 / p / 谱网格 / 峰检测 / prominence **全部不变**

`BAND_EDGE_GRID_QUANTIZATION_CORRECTED`

## 示例

见 `band_edge_examples.csv`（k=1.32/1.50 边界点）。
