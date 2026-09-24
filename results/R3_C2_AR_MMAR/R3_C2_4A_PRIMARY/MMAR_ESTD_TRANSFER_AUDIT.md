# MMAR_ESTD_TRANSFER_AUDIT

UTC: 2026-09-24T09:03:08.658004+00:00

| paper assumption | E-STD | 判定 |
| --- | --- | --- |
| range-independent | 近 CZ 近似 | APPROX |
| source depth fixed in aperture | UUV 短窗理论可 | YES_THEORETICAL；**FIXED_DEPTH_DURING_APERTURE_REQUIRED** |
| known radial increments | 理论上限 | YES_UPPER_BOUND |
| known f0 | S2 上限 | YES_S2 |
| HLA | 项目有拖曳 HLA | YES_HARDWARE_CLASS |
| θ̂=θ | 需 RC2 bearing | **NOT_FREE**；本轮 ORACLE_LOOK_DIRECTION |
| source along beam | 一般不成立 | NO_GENERAL |
| far-field r_i+ld sinθ≈r_i | 14 m 阵 @50 km | YES_NEAR |
| shallow modal density | CZ 更密 | **DIFFERENT** |
| Δ≈½ max\|φ\| | Liang 正则 | LIANG_PAPER_REGULARIZER |
| δ | Erratum | ORACLE_OFFSET_ALIGNMENT |
