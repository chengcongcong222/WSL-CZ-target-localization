# MMAR_OFFSET_RANGE_AUDIT

UTC: 2026-09-24T09:03:08.658004+00:00

论文主张 without knowing absolute source range；仿真 r0=5010 m 当作未知。

## MMAR_INHERITS_OFFSET_RANGE_CAVEAT

Yang 2018 Erratum：合成距离积分与数据距离存在 offset-range δ，
g = b_m φ_m(z_r) exp(i k_m δ)；真实数据需搜索 δ。

因此：

- 不得写 “MMAR 不需要距离”
- paper reproduction：**δ=0 = ORACLE_OFFSET_ALIGNMENT**
- 以后：RC2-conditioned δ search（本轮不做）
