# RANGE_GRID_INTERPOLATION_AUDIT

UTC: 2026-09-28T07:34:51.584019+00:00

代表样本：['A03', 'A10', 'B01', 'B08', 'C01', 'C08']
频率：201/235/283/338 Hz · z=180/200/220 · ref+alt

比较 DIRECT p(f,r(t),z) 与 GRID raw L(f,z,r) → interp → L(r(t))，
分别报 raw TL RMS、trajectory-demeaned TL RMS、demeaned correlation。

## 固定候选 grid

 grid_m  max_demeaned_rms_db  min_corr  pass
   25.0             4.223925  -0.99998 False
   10.0             2.234815  -0.99998 False
    5.0             1.074672   0.00000 False

选定：**RANGE_GRID_PATH_REJECTED_INSUFFICIENT_ACCURACY_USE_DIRECT**

主计算路径不使用 range-grid 插值，而是 F_m(r(t)) 直接求 p。
本审计记录 grid 逼近误差量级，供后续若引入 grid 加速时复用同一门槛
（demeaned RMS ≤ 0.05 dB 且 corr ≥ 0.999）。
