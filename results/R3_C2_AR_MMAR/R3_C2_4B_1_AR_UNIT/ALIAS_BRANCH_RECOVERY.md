# ALIAS_BRANCH_RECOVERY

UTC: 2026-09-24T10:06:36.216185+00:00

模型：y[n]=Σ A_q exp(−j k_q n Δr) ⇒ AR 归一化角频率 ω ≡ ± k Δr (mod 2π)

## 分类

| 类 | 条件 |
| --- | --- |
| `NO_FOLD_DIRECT_RECOVERY` | \|ω\|≤π 且 k=ω/Δr 落在物理带 [1.32,1.50] |
| `FOLDED_BUT_BRANCH_RECOVERABLE` | 折叠，但在已知物理带内唯一可恢复 |
| `FOLDED_BRANCH_AMBIGUOUS` | 带内多候选 / 无法唯一 |

## 单音结果摘要（按 Δr）

见 `SINGLE_TONE_RECOVERY.csv`（k=1.32/1.41/1.50；p=7 与 p=⌊2N_r/3⌋；Eq19 幂次 1 与 2）。

本审计**不**评价 MMAR 深度成败。
