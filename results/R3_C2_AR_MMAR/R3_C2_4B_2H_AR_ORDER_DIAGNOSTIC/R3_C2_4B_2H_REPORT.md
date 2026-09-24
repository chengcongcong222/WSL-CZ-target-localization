# R3-C2.4B-2H AR order / 20dB diagnostic

UTC: 2026-09-24T16:47:18.071304+00:00

## 判定

### `C2_4B2_BLOCKED_BY_PEAK_PICKING_AMBIGUITY`

存在多个 local maxima 但 5% 后仅 1–2 峰 → PAPER_PEAK_PICKING_RULE_NOT_EXPLICIT

## 1. SNR bug 已修

`Pn=Ps/10^(20/10)`；`SNR_AT_R0_DB` 见 look_angle_20db_fixed.csv。  
**PAPER_LOOK_ANGLE_AMBIGUITY_NONMATERIAL**（修正后）。

## 2–3. 极点 vs 选峰（p=7 未改）

见 `P7_POLE_STRUCTURE_R0_VS_20DB.csv`、`P7_LOCAL_MAXIMA_VS_PROMINENCE.csv`。

20 dB：N_POLES / N_LOCAL_MAXIMA / N_5PCT_PEAKS 分列。

## 4–5. AR order

- p=7 = `PRINTED_LITERAL_CONTROL`
- 2N/3 = **`MAX_ORDER_BOUND_SENSITIVITY`**（Ref.[18] 为 p≤2N/3 上限）
- **`PAPER_AR_ORDER_NOT_RECOVERABLE`**

## 停止

无 p sweep；不进 D(z)/Fig.5–7/E-STD/MC/P5。
