# 1E-A-CLOSEOUT-FIX

UTC: 2026-09-30T06:10:51.412522+00:00

## 修复

strict_single_bin 改用 R_keep（tau=0.5 survivor bins），不再用 R_best。

## Regression

96/96 与原始 1E-A range_anchored 一致。

## 最终统计

 A_rms_db config  n_combos  strict_single_bin  truth_retained  true_rank1  best_range_unique_50  min_wrong_range_margin  median_r_width  max_r_width  max_dev
     0.00 TRIPLE        12                 12              12          12                    12                0.601456             0.0          0.0      0.0
     0.00   FOUR        12                 12              12          12                    12                0.552272             0.0          0.0      0.0
     0.25 TRIPLE        12                 11              12          12                    12                0.487401             0.0          1.0      1.0
     0.25   FOUR        12                  8              12          12                    12                0.421123             0.0          1.0      1.0
     0.50 TRIPLE        12                  6              12          12                    12                0.384283             0.5          1.0      1.0
     0.50   FOUR        12                  6              12          12                    12                0.318109             0.5          2.0      2.0
     1.00 TRIPLE        12                  1              12          12                    12                0.242041             7.5         15.0     10.0
     1.00   FOUR        12                  1              12          12                    12                0.182269             3.0         15.0     10.0

## 三层

optimal rank → survivor-set contraction → strict single-bin anchor

## 判定

**AMPLITUDE_VARIATION_BREAKS_TURN_RANGE_ANCHOR_IN_TESTED_RANGE** = STRICT_SINGLE_RANGE_BIN_CRITERION_NOT_UNIFORMLY_RETAINED

辅助: OPTIMAL_TRUTH_RANK_ROBUST_TO_TESTED_1DB_AMPLITUDE_VARIATION
      RANGE_SET_CONTRACTION_RETAINS_LOCAL_VALUE_THROUGH_TESTED_0P5DB
