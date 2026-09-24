# R3-C2.4B-1F2 complex modified covariance

UTC: 2026-09-24T11:37:52.613113+00:00

## 判定

### `AR_K_RECOVERY_NUMERICALLY_UNSTABLE`

P1=P1_COMPLEX_EXP_IDENTITY_VALIDATED; F1_single=0.444; multi_complete=0.75; multi_F1=0.942

## 修复

1. **`COMPLEX_BACKWARD_CONJUGATION_REQUIRED`**：Ab=conj(y[...])，bb=−conj(y[n−p])
2. **p=1 硬门** `P1_COMPLEX_EXP_IDENTITY_VALIDATED`：a1=−e^{−jkΔr}；ω=wrap(−kΔr)
3. sign mapping 由解析恒等式冻结（无 majority vote）
4. rank：`DATA_MODEL_RANK_LIMIT` ≠ 数值不稳；另做 ε 扰动 conditioning diagnostic

## 指标（power=1）

见 `SINGLE_TONE_RECOVERY_FINAL.csv` / `MULTITONE_RECOVERY_FINAL.csv`。

## 停止

不进 4B-2 / Hankel / D(z) / FIELD / E-STD / MC / P5。
