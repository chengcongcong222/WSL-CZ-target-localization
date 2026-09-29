# RANGE-UB-1C-FIX2

**SINGLE_DEPTH_AUTOCORR_OBSERVABLE_NOT_NUMERICALLY_STABLE**

早时延自相关形状 (0-5s, 5ms bins, 归一化) 1601 vs 3201 收敛: 4/12 cosine>=0.99。
L2 全部通过 (max 0.008)，但 cosine 0.97-0.995，8点略低于 0.99 门槛。

按预注册 gate: 关闭单深度可观测 TDOA 支线。
保留 FIX2 all-eigenray upper bound。
下一步: VERTICAL_APERTURE_TDOA_REFERENCE_UPPER_BOUND。
