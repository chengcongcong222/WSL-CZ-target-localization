# R3-CLOSEDLOOP-1D — Frequency-Subset Sufficiency of Turn-Assisted Range Anchor

UTC: 2026-09-30T02:57:39.145995+00:00

基线：`e3030060e3bc36fc35e3232e90a2dbb9ead92da9`

## RC2 Identity

n_RC2 = 385 (expected 385), truth_in = True

## Subset Anchor Summary

         subset  n_lines  three_depth_range_anchor  n_anchored
            201        1                     False           0
            235        1                     False           0
            283        1                     False           0
            338        1                     False           0
        201+235        2                     False           0
        201+283        2                     False           1
        201+338        2                     False           0
        235+283        2                     False           1
        235+338        2                     False           1
        283+338        2                     False           0
    201+235+283        3                      True           3
    201+235+338        3                     False           1
    201+283+338        3                     False           1
    235+283+338        3                     False           2
201+235+283+338        4                      True           3

## 判定

### `SOME_THREE_LINE_TURN_RANGE_ANCHOR`

any_single=False, all_single=False, any_double=False, all_double=False, any_triple=True, four=True

## Frequency Complementarity

 freq_hz  n_subsets_containing  n_anchored_subsets
   201.0                     8                   2
   235.0                     8                   2
   283.0                     8                   2
   338.0                     8                   1

## 未做

SSP、幅漂、频漂、噪声、P5。
