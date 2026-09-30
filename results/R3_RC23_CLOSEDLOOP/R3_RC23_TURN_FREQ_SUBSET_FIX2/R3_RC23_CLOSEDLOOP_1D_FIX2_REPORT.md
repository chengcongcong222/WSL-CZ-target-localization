# R3-CLOSEDLOOP-1D-FIX2 — Independent Identity & Complementarity

UTC: 2026-09-30T03:35:27.020040+00:00

## Identity (independent 1C recomputation)

- RC2: PASS (n=385, symmetric_diff=0)
- 235: PASS (max_J_err < 1e-10)
- FOUR: PASS (max_J_err < 1e-10)

## Anchor Summary

         subset  n_lines  three_depth_range_anchor
            201        1                     False
            235        1                     False
            283        1                     False
            338        1                     False
        201+235        2                     False
        201+283        2                     False
        201+338        2                     False
        235+283        2                     False
        235+338        2                     False
        283+338        2                     False
    201+235+283        3                      True
    201+235+338        3                     False
    201+283+338        3                     False
    235+283+338        3                     False
201+235+283+338        4                      True

## 判定

### `SOME_THREE_LINE_TURN_RANGE_ANCHOR`

identity=PASS, any_single=False, any_double=False, any_triple=True, four=True

辅助：`MULTIFREQUENCY_ALIAS_COMPLEMENTARITY_CONFIRMED`

正式口径：`RMS_FUSION_EXPLOITS_COMPLEMENTARY_RANGE_ALIAS_STRUCTURE`

## Cross-Frequency Alias Incompatibility

共 6 个错误距离被三频联合排除但单频均可解释。

 z_true_m  r_km  Jstar_201  Jstar_235  Jstar_283  Jstar_triple  CROSS_FREQUENCY_ALIAS_INCOMPATIBILITY
    200.0  47.0   0.351877   0.341078   0.438410      1.412065                                   True
    220.0  46.0   0.289855   0.498899   0.302961      1.287387                                   True
    220.0  47.0   0.296479   0.146314   0.489166      0.847470                                   True
    220.0  49.0   0.447528   0.138133   0.284108      0.909078                                   True
    220.0  53.0   0.208507   0.335541   0.289120      0.879489                                   True
    220.0  58.0   0.377568   0.222020   0.446018      1.630920                                   True
