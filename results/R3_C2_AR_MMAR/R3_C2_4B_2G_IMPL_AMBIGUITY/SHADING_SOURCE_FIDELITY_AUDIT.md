# SHADING_SOURCE_FIDELITY_AUDIT

UTC: 2026-09-24T16:28:07.573303+00:00

## 当前

`OUR_RANGE_AVERAGE_GLOBAL_MEAN`：S(r)=const=[mean_aperture |B|²]^{-1/2}

**不能**补偿 1/√r 扩展（Liang Eq.8 / Yang 广义 Hankel：⟨·⟩=range averaging，使 S≈√r）。

## 控制（非 Liang 实现）

`ORACLE_SQRT_R_SPREADING_CONTROL`：S(r)=√(r/r0)

只用于判断“20 dB 多峰缺失是否主因在 shading”；**不得**冒充数据驱动 Eq.(8)。

`PAPER_RANGE_SMOOTHING_WINDOW_NOT_EXPLICIT`（primary/Ref7 未给窗口；禁止找最佳窗）。
