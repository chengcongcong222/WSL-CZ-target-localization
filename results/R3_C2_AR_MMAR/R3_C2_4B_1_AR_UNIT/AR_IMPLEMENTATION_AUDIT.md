# AR_IMPLEMENTATION_AUDIT

UTC: 2026-09-24T10:06:36.216185+00:00
基线 commit: c6a49c1b6e033941d6e0df0f4fff1c28a78995a8
来源仅：MMAR_EQUATION_LOCK_FINAL / AR_INTERPRETATION_GATE / AR_ORDER_PREREGISTRATION / SPATIAL_SAMPLING_ALIAS_AUDIT / MMAR_PAPER_CONFIG_FINAL

## AR

- **modified covariance**（forward-backward LS），对应论文选型（避免 line splitting）
- 输入：y[n]=Σ A_q e^{-j k_q n Δr}（真值 k 冻结后生成）
- 只恢复**峰位置** k̂；不用 AR peak height 作幅度

## Eq.(19)

- `PRINTED_EQ19_MODULUS_POWER_1`
- `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`
- 只比峰位置，不按深度择优

## AR order（预注册，未改）

- `PRINTED_LITERAL_CONTROL`: p=7
- `OUR_MOVING_SAMPLE_INTERPRETATION`: p=⌊2N_r/3⌋=133（N_r=200）
- **无额外 p 扫描**

## Eq.(24)

- `OUR_EXACT_SOLVER_FOR_EQ24`（DP）
- 与 brute-force 小规模交叉验证：cost_match=True
- 人工 mode set 上 MODE_ORDER_RECOVERY_RATE=0.950

## 本轮禁止（已遵守）

Hankel g、D(z)、Liang Fig5/6、FIELD、E-STD、MC、P5
