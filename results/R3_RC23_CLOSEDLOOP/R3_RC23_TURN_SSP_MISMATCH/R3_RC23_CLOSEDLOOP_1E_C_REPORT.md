# R3-CLOSEDLOOP-1E-C — SSP Mismatch

UTC: 2026-09-30T09:50:34.361217+00:00

## SSP Stress

E1: shift 50m; E2: shift 150m + 2m/s (4G frozen)

## 判定

### `SSP_MISMATCH_DEGRADES_RANGE_RANKING_IN_TESTED_STRESS`

all_strict=False, triple=False, four=False, rank1=False, best50=False

## Working Range

config env  n_cases  strict  truth_ok  rank1  best50  max_width  min_margin
TRIPLE  E0        3       3         3      3       3        0.0    0.601456
TRIPLE  E1        3       0         2      0       0       15.0   -0.899328
TRIPLE  E2        3       0         2      0       0       14.0   -1.322937
  FOUR  E0        3       3         3      3       3        0.0    0.552272
  FOUR  E1        3       0         2      0       0       14.0   -0.943773
  FOUR  E2        3       0         2      0       0       14.0   -0.902284

## Depth Profiling Compensation

config env  z_true  z_star_truth  z_shift   true_J
TRIPLE  E0   180.0         180.0      0.0 0.000000
TRIPLE  E0   200.0         200.0      0.0 0.000000
TRIPLE  E0   220.0         220.0      0.0 0.000000
TRIPLE  E1   180.0         190.0     10.0 1.253686
TRIPLE  E1   200.0         195.0     -5.0 0.776750
TRIPLE  E1   220.0         190.0    -30.0 2.040296
TRIPLE  E2   180.0         190.0     10.0 2.599999
TRIPLE  E2   200.0         170.0    -30.0 0.666473
TRIPLE  E2   220.0         215.0     -5.0 0.763016
  FOUR  E0   180.0         180.0      0.0 0.000000
  FOUR  E0   200.0         200.0      0.0 0.000000
  FOUR  E0   220.0         220.0      0.0 0.000000
  FOUR  E1   180.0         155.0    -25.0 1.759084
  FOUR  E1   200.0         195.0     -5.0 0.832646
  FOUR  E1   220.0         195.0    -25.0 2.290666
  FOUR  E2   180.0         155.0    -25.0 2.738334
  FOUR  E2   200.0         170.0    -30.0 0.701756
  FOUR  E2   220.0         230.0     10.0 1.202992

## 未做

幅漂、频漂、跟踪误差、联合corner、P5。
