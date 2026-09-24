# RANK_VS_CONDITIONING_AUDIT

UTC: 2026-09-24T11:37:52.613113+00:00

## DATA_MODEL_RANK_LIMIT

无噪声 q 个复指数是**低秩模型**：单音 rank≪p 不等于 p=133 在真实多模上不稳。
单音 rank deficiency 记为 `DATA_MODEL_RANK_LIMIT`。

## NUMERICAL_CONDITIONING

独立 diagnostic：固定 seed、ε=1e-6/1e-4 复高斯扰动，看系数/峰位敏感度
（`numerical_conditioning_diagnostic.csv`）。**不计入** MMAR 成功率，**不**用来调 p。

## 性能分支（未改）

p=7 `PRINTED_LITERAL_CONTROL`；p=133 `OUR_MOVING_SAMPLE_INTERPRETATION`。
