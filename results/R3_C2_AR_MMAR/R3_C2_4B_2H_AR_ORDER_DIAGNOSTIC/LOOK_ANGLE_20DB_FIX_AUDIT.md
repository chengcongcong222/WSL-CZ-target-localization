# LOOK_ANGLE_20DB_FIX_AUDIT

UTC: 2026-09-24T16:47:18.071304+00:00

## Bug 修复

错误：`Pn = Ps / 10**20` ⇒ SNR=200 dB  
正确：`Pn = Ps / 10**(20/10)` = Ps/100 ⇒ SNR=20 dB

反算 `SNR_AT_R0_DB`：见 `look_angle_20db_fixed.csv`（均应 ≈20）。

## 冻结（仅在修正后成立）

**PAPER_LOOK_ANGLE_AMBIGUITY_NONMATERIAL**（θ=0/30/60/90°，4990 m，p=7，20 dB）

R0 下非材料性结论保留；**20 dB 版本以本表为准**。
