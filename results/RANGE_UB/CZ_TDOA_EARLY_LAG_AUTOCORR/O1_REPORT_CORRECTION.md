# O1_REPORT_CORRECTION

UTC: 2026-09-29T06:20:38.114210+00:00

## 纠正

旧 1C-FIX 摘要写 "O1 convergence PASS 12/12" 有误。

实际 O1 = 10/12：
- zs=180, r=56: matched_fraction=0.754 < 0.8
- zs=200, r=50: matched_fraction=0.455 < 0.8

正确标签：`O1_STRONG_PATH_MEMBERSHIP_NOT_FULLY_CONVERGED`

## 保留现象

匹配上的强路径时延本身非常稳定（RMS 0.02–0.11 ms）。
不稳定来自个别 eigenray 跨越 -20 dB 选择边界。

不改变旧 1C-FIX 主判。
