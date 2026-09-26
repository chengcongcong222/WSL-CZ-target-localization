# R3-RC3-REANCHOR-5A trackable-line availability

UTC: 2026-09-26T03:58:23.823471+00:00

## 判定

### `ANY_SINGLE_OF_FOUR_TRACKABLE_LINES_SUFFICIENT_IDEAL_CONTROL`

4 个单频子集全部 PASS（理想稳定线谱条件）

## 轴

`S1_BRIDGE_LINE_AVAILABILITY_ONLY` · E0 · 15 非空子集 · 双侧 depth profile · A05 剔除

子集 PASS：fraction(ΔJ>0)≥0.8 且 median ΔJ>0

## 按线数

 n_lines  n_subsets    worst_subset  worst_frac_gt0  worst_median_dJ     best_subset  best_frac_gt0  median_of_subset_medians  n_pass  all_pass
       1          4             201             1.0         1.785184             201            1.0                  1.812488       4      True
       2          6         201;235             1.0         2.035687         201;235            1.0                  2.287620       6      True
       3          4     201;235;283             1.0         2.542427     201;235;283            1.0                  2.570106       4      True
       4          1 201;235;283;338             1.0         2.712532 201;235;283;338            1.0                  2.712532       1      True

## 含义边界

理想稳定线谱下的**结构性**下限；**不是** S1 已验证。
下一阶段：5B 线谱幅度/频率不稳定性。

## 文档修正

- index：4F→**4G**，目录 `..._FIX2`
- 4G 文案：B/C 签名非重复；A01 静态重复属预期（不改 CSV）
