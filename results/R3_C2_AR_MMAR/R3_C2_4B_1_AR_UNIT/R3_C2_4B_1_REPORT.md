# R3-C2.4B-1 AR 波数恢复单元测试

UTC: 2026-09-24T10:06:36.216185+00:00

## 判定

### `AR_K_RECOVERY_PARTIAL_WITH_ALIAS_LIMIT`

single_ok=0.78, fold_amb=True, fold_rec=True, nofold=False

## 指标

- 单音 |k̂−k|<1e-3 比例：**0.778**
- 多音全恢复比例：**0.833**
- Eq24 solver cost_match：**True**
- 合成 mode set MODE_ORDER_RECOVERY_RATE：**0.950**

## 两分支 p（未调参）

p=7（PRINTED_LITERAL_CONTROL）与 p=⌊2N_r/3⌋=133（OUR_MOVING_SAMPLE_INTERPRETATION），N_r=200。

## 别名

['FOLDED_BRANCH_AMBIGUOUS', 'FOLDED_BUT_BRANCH_RECOVERABLE']

详见 `ALIAS_BRANCH_RECOVERY.md`。

## 停止

不进 4B-2 / FIELD / E-STD / MC / P5。
