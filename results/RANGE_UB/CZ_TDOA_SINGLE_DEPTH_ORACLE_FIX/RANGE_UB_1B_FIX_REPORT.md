# RANGE-UB-1B-FIX — Correct BELLHOP Arrival Parser + Full r-z Oracle

UTC: 2026-09-29T04:27:27.153162+00:00

基线：`4c6e8f55f0e3d5849f2f39ca1c81fef1efbc9268`

## 旧 1B 失效原因

1. `.arr` parser 把 phase (col1) 误读为 delay — 正确为 col2
2. 候选 z 只算 180/200/220 而非 150:5:250 全 21 深度
3. sorted-prefix 不是 path-metadata association

## 修复

| 项 | 状态 |
|---|---|
| Parser gate | PASS（delay 正、秒量级） |
| Env integrity | PASS（rigid bottom 与 E-STD 一致） |
| 候选网格 | **651** = 31 r × 21 z_s |
| Path association | PATH_METADATA_ORACLE_ASSOCIATION |
| TDOA | pairwise，单位 ms |

## 判定

### `SINGLE_DEPTH_ORACLE_TDOA_RANGE_INFORMATION_CONFIRMED`

range_identified=True, depth_identified=True, aux=[]; ranks=[np.int64(1), np.int64(1), np.int64(1)]; n_valid=[np.int64(651), np.int64(651), np.int64(651)]

- **RANGE_IDENTIFIABILITY**: PASS
- **DEPTH_IDENTIFIABILITY_AT_TRUE_RANGE**: PASS

## 排序诊断

 z_true_m  true_rank  true_J_ms  n_valid_scored
    180.0          1        0.0             651
    200.0          1        0.0             651
    220.0          1        0.0             651

## 距离别名审计

 z_true_m  r_track_km  present  J_tau_star_ms  delta_J_ms  best_z_star  n_matched  machine_degenerate
    180.0        45.0     True    1144.348594 1144.348594        235.0         64               False
    180.0        50.0     True       0.000000    0.000000        180.0         99                True
    180.0        56.0     True    1274.678527 1274.678527        205.0         63               False
    180.0        58.0     True    1773.283007 1773.283007        150.0         68               False
    180.0        60.0     True    2217.035323 2217.035323        190.0         63               False
    200.0        45.0     True    1146.690106 1146.690106        235.0         65               False
    200.0        50.0     True       0.000000    0.000000        200.0         99                True
    200.0        56.0     True    1323.217140 1323.217140        205.0         64               False
    200.0        58.0     True    1764.053694 1764.053694        150.0         67               False
    200.0        60.0     True    2235.981632 2235.981632        190.0         64               False
    220.0        45.0     True    1122.302105 1122.302105        235.0         65               False
    220.0        50.0     True       0.000000    0.000000        220.0         98                True
    220.0        56.0     True    1292.494834 1292.494834        215.0         63               False
    220.0        58.0     True    1762.705879 1762.705879        150.0         67               False
    220.0        60.0     True    2234.351981 2234.351981        190.0         64               False

## 未做

VLA Oracle、论文 Radon、测量误差扫描、RC2/RC3 融合、平台转向、S1、P5。
