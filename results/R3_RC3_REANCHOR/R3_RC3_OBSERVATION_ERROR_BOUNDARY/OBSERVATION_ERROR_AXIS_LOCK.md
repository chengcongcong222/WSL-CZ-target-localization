# OBSERVATION_ERROR_AXIS_LOCK

UTC: 2026-09-25T10:12:48.225195+00:00

`TL_OBSERVATION_ERROR_SENSITIVITY_AXIS`：

ε = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0] dB（每频零均值 TL 形状扰动，扰动后再去均值）

seeds = [0, 1, 2, 3, 4]（有限敏感性，**非** Monte Carlo 性能概率）

信号条件：`S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND`（201/235/283/338 Hz 持续可跟踪；每频独立去均值；不要求源级；不声称 S0/S1）

指标：ΔJ = J_alt* − J_true*（两者均对 z∈150:5:250 profile）

角色不变：`DEPTH_PROFILED_NUISANCE_VARIABLE`；标签 `IDEAL_PROFILED_PROPAGATION_MARGIN_OBSERVED`（非 96.7% 性能）。
