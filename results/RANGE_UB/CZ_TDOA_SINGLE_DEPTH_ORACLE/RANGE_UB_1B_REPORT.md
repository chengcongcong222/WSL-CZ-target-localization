# RANGE-UB-1B — E-STD 单深度 Oracle 多途时延距离上界

UTC: 2026-09-29T03:38:00.102431+00:00

基线：`dda967d77c1870217be8973d279c6e0d4533dad1`

**不是复现徐嘉璘2024论文。** 不使用论文未恢复的公式。

## 问题

> 在 E-STD、单接收深度 z_r=200 m，若理想获得多途到达时延，时延差本身能否打破 r-z 简并？

## 判定

### `SINGLE_DEPTH_ORACLE_TDOA_RANGE_INFORMATION_PARTIAL`

rank1_all=False, unique_min_at_50=True, alias_persists=False; ranks=[np.int64(1), np.int64(1), np.int64(2)]; frac_2plus_arrivals=1.00

## Oracle 假设

- PERFECT_MULTIPATH_DETECTION
- PERFECT_ARRIVAL_ASSOCIATION
- ZERO_DELAY_MEASUREMENT_ERROR

## Arrival 完整性

检查点数 48，≥2 arrivals: 48 (100.0%)

## 排序诊断

 z_true_m  true_rank  true_J  best_false_J  best_false_r_km  best_false_z_s_m  n_scored
    180.0          1     0.0     20.124612             45.0             180.0        93
    200.0          1     0.0      0.000000             50.0             220.0        93
    220.0          2     0.0      0.000000             50.0             200.0        93

## 距离别名审计

 z_true_m  r_track_km  present  J_tau_star    delta_J  best_z_star  n_at_r
    180.0        45.0     True   20.124612  20.124612        180.0       3
    180.0        50.0     True    0.000000   0.000000        180.0       3
    180.0        56.0     True  158.745068 158.745068        180.0       3
    180.0        58.0     True  128.860382 128.860382        180.0       3
    180.0        60.0     True  184.445099 184.445099        180.0       3
    200.0        45.0     True   28.460499  28.460499        180.0       3
    200.0        50.0     True    0.000000   0.000000        200.0       3
    200.0        56.0     True  165.680404 165.680404        180.0       3
    200.0        58.0     True  133.491566 133.491566        180.0       3
    200.0        60.0     True  189.855193 189.855193        180.0       3
    220.0        45.0     True   29.199856  29.199856        180.0       3
    220.0        50.0     True    0.000000   0.000000        200.0       3
    220.0        56.0     True  165.680404 165.680404        180.0       3
    220.0        58.0     True  135.394156 135.394156        180.0       3
    220.0        60.0     True  184.676092 184.676092        180.0       3

## 未做

噪声扫描、RC2/RC3 融合、垂直阵、论文公式复现、平台转向、P5。
