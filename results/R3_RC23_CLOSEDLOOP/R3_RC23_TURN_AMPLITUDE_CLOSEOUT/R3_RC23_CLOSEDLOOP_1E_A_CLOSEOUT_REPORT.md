# R3-CLOSEDLOOP-1E-A-CLOSEOUT

UTC: 2026-09-30T04:37:34.734546+00:00

## Identity

PASS (independent comparison with 1D-FIX2)

## Performance Layers

| Layer | Finding |
|---|---|
| Truth retained | PASS at 0.25/0.5/1.0 dB |
| True rank = 1 | PASS at 0.25/0.5/1.0 dB |
| Strict single-bin | NOT uniform (original FAIL label retained) |
| Set width | 1 km @ 0.25 dB → 2 km @ 0.5 dB → up to 15 km @ 1 dB |

## Working Range Summary

 A_rms_db config  n_combos  strict_single_bin  truth_retained  best_range_50  min_wrong_range_margin  median_r_width  max_r_width  max_dev
     0.00 TRIPLE        12                 12              12             12                0.601456             0.0          0.0      0.0
     0.00   FOUR        12                 12              12             12                0.552272             0.0          0.0      0.0
     0.25 TRIPLE        12                 12              12             12                0.487401             0.0          1.0      1.0
     0.25   FOUR        12                 12              12             12                0.421123             0.0          1.0      1.0
     0.50 TRIPLE        12                 12              12             12                0.384283             0.5          1.0      1.0
     0.50   FOUR        12                 12              12             12                0.318109             0.5          2.0      2.0
     1.00 TRIPLE        12                 12              12             12                0.242041             7.5         15.0     10.0
     1.00   FOUR        12                 12              12             12                0.182269             3.0         15.0     10.0

## Injected vs Residual RMS

 shape  A_raw_rms_db  A_after_window_demean_rms_db  W1_residual_rms_db  W2_residual_rms_db
A-RAMP          0.25                      0.125000            0.126020            0.123954
A-RAMP          0.50                      0.250000            0.252041            0.247908
A-RAMP          1.00                      0.500000            0.504082            0.495816
A-RAMP          2.00                      1.000000            1.008163            0.991632
A-SINE          0.25                      0.110804            0.112199            0.109368
A-SINE          0.50                      0.221609            0.224398            0.218736
A-SINE          1.00                      0.443218            0.448796            0.437473
A-SINE          2.00                      0.886435            0.897593            0.874946

## 判定

**AMPLITUDE_VARIATION_BREAKS_TURN_RANGE_ANCHOR_IN_TESTED_RANGE**

含义：STRICT_SINGLE_RANGE_BIN_CRITERION_NOT_UNIFORMLY_RETAINED

工程口径：在冻结门限下，严格单距离格保持对慢变幅漂敏感。0.25–0.5 dB 条件下最优候选排序仍保持，距离集合仍有明显收缩。尚未得到真实系统的幅度稳定性要求。
