# 1E-B-FIX

UTC: 2026-09-30T07:49:13.619437+00:00

## 频率模型

44 个唯一频率全部构建成功（5B env 方法）。

## D-effect Gate

非零 D 与 D=0 RMS 差异 > 1e-6 dB：PASS

## 判定

### `TRIPLE_AND_FOUR_STRICT_RANGE_ANCHOR_SURVIVES_TRACKED_2PCT_FREQ_DRIFT`

triple_strict=True, four_strict=True, rank1=True, best50=True

## Working Range

config      D  D_pct  strict  truth_ok  rank1  best_50  max_width  min_margin
TRIPLE 0.0000   0.00       3         3      3        3        0.0    0.601456
TRIPLE 0.0025   0.25       3         3      3        3        0.0    1.920934
TRIPLE 0.0050   0.50       3         3      3        3        0.0    1.517021
TRIPLE 0.0100   1.00       3         3      3        3        0.0    1.534708
TRIPLE 0.0200   2.00       3         3      3        3        0.0    1.130036
  FOUR 0.0000   0.00       3         3      3        3        0.0    0.552272
  FOUR 0.0025   0.25       3         3      3        3        0.0    2.382146
  FOUR 0.0050   0.50       3         3      3        3        0.0    2.500552
  FOUR 0.0100   1.00       3         3      3        3        0.0    1.470907
  FOUR 0.0200   2.00       3         3      3        3        0.0    1.580886
