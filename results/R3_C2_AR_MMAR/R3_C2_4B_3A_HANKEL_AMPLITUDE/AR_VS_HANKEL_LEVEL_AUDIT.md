# AR_VS_HANKEL_LEVEL_AUDIT

UTC: 2026-09-25T03:59:05.755902+00:00

`AR_PEAK_LEVEL_NOT_MODAL_AMPLITUDE`

- AR 峰高只作诊断列 `ar_peak_level_NOT_amplitude`
- 模态幅度 **只** 用 Liang Eq.(7) generalized Hankel 在 k̂_m 处：
  g = e^{jπ/4}/√(2π k̂) ∫ B(r) e^{j k̂ r} S(r) dr

两条 shading：
- CURRENT_DATA_SHADING（全局 mean 常数）
- ORACLE_SQRT_R_SPREADING_CONTROL（机制上限）

本轮 **R0 noiseless only**。
