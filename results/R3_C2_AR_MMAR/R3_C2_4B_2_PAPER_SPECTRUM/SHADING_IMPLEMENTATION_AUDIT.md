# SHADING_IMPLEMENTATION_AUDIT

UTC: 2026-09-24T12:40:54.120845+00:00

Liang Eq.(8): S(r)=⟨|B(r)|²⟩^{-1/2}

实现：`OUR_RANGE_AVERAGE_GLOBAL_MEAN` — 对 |B|² 在整个 range 序列上取均值后开方取负一次幂
（等价于尺度补偿的全局归一；主文未指定滑动窗）。

**不按 Fig.3/4 调 window。** 若后续需要局部平滑，单独标 interpretation。
来源：Liang Eq.(8) + Yang generalized-Hankel 实践。
