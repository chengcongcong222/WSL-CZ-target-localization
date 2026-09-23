# R3-C2.1-PARSER-FINAL

UTC：2026-09-23T08:12:42.907979+00:00

撤销：**不得写 Parser OK**；路线保持未判定。
正式状态见 DECISION；**不评价 Yang**。

## 原始 header 审计（235 Hz 已知硬事实）
               file  f_hz    zs  recl_words  recl_bytes  Nmedia  Ntot  NMat         depths  M_mod  M_prt  file_size   k_off MODE_COUNT_DISCREPANCY  depths_n depth_covers_150_250_step2  z_min  z_max
estd_f235_zs180.mod 235.0 180.0          32       128.0     1.0     2     2 [180.0, 200.0]    749  749.0   102784.0 96768.0                  False       NaN                        NaN    NaN    NaN
estd_f235_zs200.mod 235.0 200.0          32       128.0     1.0     1     1        [200.0]    749  749.0   102784.0 96768.0                  False       NaN                        NaN    NaN    NaN
estd_f235_zs220.mod 235.0 220.0          32       128.0     1.0     2     2 [200.0, 220.0]    749  749.0   102784.0 96768.0                  False       NaN                        NaN    NaN    NaN
     zgrid_f235.mod 235.0 200.0         102         NaN     NaN    51    51            NaN    749    NaN        NaN     NaN                    NaN      51.0                       True  150.0  250.0

## 12 case 矩阵
 f_hz    zs  header_ok  depth_ok  M_ok  k_checksum_ok  phi_ok  overall_ok  M_mod  M_prt MODE_COUNT_DISCREPANCY  case
201.0 180.0       True      True  True           True    True        True    641  641.0                  False   NaN
201.0 200.0       True      True  True           True    True        True    641  641.0                  False   NaN
201.0 220.0       True      True  True           True    True        True    641  641.0                  False   NaN
235.0 180.0       True      True  True           True    True        True    749  749.0                  False   NaN
235.0 200.0       True      True  True           True    True        True    749  749.0                  False   NaN
235.0 220.0       True      True  True           True    True        True    749  749.0                  False   NaN
283.0 180.0       True      True  True           True    True        True    902  903.0                   True   NaN
283.0 200.0       True      True  True           True    True        True    902  903.0                   True   NaN
283.0 220.0       True      True  True           True    True        True    902  903.0                   True   NaN
338.0 180.0       True      True  True           True    True        True   1078 1078.0                  False   NaN
338.0 200.0       True      True  True           True    True        True   1078 1078.0                  False   NaN
338.0 220.0       True      True  True           True    True        True   1078 1078.0                  False   NaN
201.0 200.0       True      True  True           True    True        True    641    NaN                    NaN zgrid
235.0 200.0       True      True  True           True    True        True    749    NaN                    NaN zgrid
283.0 200.0       True      True  True           True    True        True    902    NaN                    NaN zgrid
338.0 200.0       True      True  True           True    True        True   1078    NaN                    NaN zgrid

## 283 MODE_COUNT_DISCREPANCY
见 `kraken_283_mode_count_audit.md`（.mod 902 vs .prt 903）

## 跨 zs 一致性
 f_hz  zs_a  zs_b  phi200_abs_corr  max_abs_dk  k_consistent  phi_consistent
201.0 180.0 200.0              1.0         0.0          True            True
235.0 180.0 200.0              1.0         0.0          True            True
283.0 180.0 200.0              1.0         0.0          True            True
338.0 180.0 200.0              1.0         0.0          True            True

## 判定 `C2_1_PARSER_VALIDATED`

12/12 case ok=True, missing=[], zgrid_cover=True, cross-zs consistent=True, MODE_COUNT_DISCREPANCY_283=True. phi now read only NMat complex from real depth record (no linspace fallback).

本轮无 FWHM / z_hat / Yang pass-fail。
