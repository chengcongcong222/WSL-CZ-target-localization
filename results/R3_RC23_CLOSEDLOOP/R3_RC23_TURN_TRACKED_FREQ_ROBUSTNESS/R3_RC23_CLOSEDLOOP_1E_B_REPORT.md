# R3-CLOSEDLOOP-1E-B — Tracked Frequency Drift

UTC: 2026-09-30T07:28:04.491825+00:00

## 判定

### `FOUR_REDUNDANCY_IMPROVES_TRACKED_FREQUENCY_DRIFT_ROBUSTNESS`

triple_strict=False, four_strict=True, rank1=True, best50=True

## Working Range

config      D  D_pct  strict  truth_ok  rank1  best_50  max_width  min_margin
TRIPLE 0.0000   0.00       3         3      3        3        0.0    0.601456
TRIPLE 0.0025   0.25       1         3      3        3        6.0    0.228291
TRIPLE 0.0050   0.50       1         3      3        3        6.0    0.228291
TRIPLE 0.0100   1.00       1         3      3        3        6.0    0.228291
TRIPLE 0.0200   2.00       1         3      3        3        6.0    0.228291
  FOUR 0.0000   0.00       3         3      3        3        0.0    0.552272
  FOUR 0.0025   0.25       3         3      3        3        0.0    0.541104
  FOUR 0.0050   0.50       3         3      3        3        0.0    0.541104
  FOUR 0.0100   1.00       3         3      3        3        0.0    0.541104
  FOUR 0.0200   2.00       3         3      3        3        0.0    0.541104

## D=0 Identity

config  z_true_m  ids_match    max_J_err  pass
TRIPLE     180.0       True 8.881784e-16  True
TRIPLE     200.0       True 8.881784e-16  True
TRIPLE     220.0       True 8.881784e-16  True
  FOUR     180.0       True 8.881784e-16  True
  FOUR     200.0       True 8.881784e-16  True
  FOUR     220.0       True 8.881784e-16  True

## 未做

F_NOMINAL、跟踪误差、幅漂、SSP、P5。
