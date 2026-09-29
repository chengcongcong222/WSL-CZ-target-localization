# RANGE-UB-2A — Xu2024 VLA Formula Numerical Closure

UTC: 2026-09-29T07:06:26.331855+00:00

基线：`c6de5bee23bcc4f267b79c7264d1caf895847acf`

## 论文条件

H = 4314.5 m, R = 50 km, c = 1525.0 m/s, VLA 20–1620 m, LFM 100–300 Hz

## 判定

### `XU2024_VLA_RANGE_FORMULA_NUMERICALLY_CLOSED`

theory_slopes=PASS, table1=PASS, table2=PASS, units=PASS

## 理论斜率闭合

slope  computed_ms_per_m  paper_target_ms_per_m  abs_error  gate_5e-4
  k21           0.226334                 0.2263   0.000034       True
  k31          -0.113167                -0.1132   0.000033       True
  k41           0.339502                 0.3395   0.000002       True
  k32          -0.339502                -0.3395   0.000002       True
  k42           0.113167                 0.1132   0.000033       True

## Table 1 闭合

         expr  k_est_ms_per_m  R_calc_km  R_paper_km  err_R_km  err_pct_calc  err_pct_paper  pass
 R≈4H/(c·k21)          0.2508  45.122493       45.12  0.002493      9.755013            9.8  True
R≈-2H/(c·k31)         -0.1056  53.582961       53.58  0.002961      7.165922            7.2  True
 R≈6H/(c·k41)          0.3594  47.231725       47.23  0.001725      5.536550            5.5  True
R≈-6H/(c·k32)         -0.3476  48.835103       48.83  0.005103      2.329793            2.3  True
 R≈2H/(c·k42)          0.1083  52.247097       52.25  0.002903      4.494195            4.5  True

## Table 2 二次相关闭合

         expr  k_est_ms_per_m  R_calc_km  R_paper_km  err_R_km  err_pct_calc  err_pct_paper  formula_matches_paper_error                                                       note  pass
R≈-6H/(c·k32)         -0.3453  49.160388       49.16  0.000388      1.679224            1.7                         True                                                             True
 R≈2H/(c·k42)          0.1088  52.006991       53.84  1.833009      4.013983            4.0                         True paper R value likely typo; formula matches paper error_pct  True
 R≈8H/(c·k43)          0.4620  48.990136       48.99  0.000136      2.019729            2.0                         True                                                             True

## 机制边界

- `XU2024_RANGE_OBSERVABLE_IS_DELAY_DIFFERENCE_VS_RECEIVER_DEPTH_SLOPE`
- `VERTICAL_DEPTH_SAMPLING_IS_ESSENTIAL_TO_PAPER_METHOD`

## 未做

E-STD VLA、BELLHOP、RC2/RC3、噪声、P5。
