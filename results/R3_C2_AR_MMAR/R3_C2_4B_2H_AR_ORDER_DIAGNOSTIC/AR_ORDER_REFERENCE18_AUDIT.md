# AR_ORDER_REFERENCE18_AUDIT

UTC: 2026-09-24T16:47:18.071304+00:00

| 来源 | 含义 |
| --- | --- |
| Liang Eq.(18) 印刷 | p often (2/3)(2L+1)；L=5 ⇒ **p=7** `PRINTED_LITERAL_CONTROL` |
| Ref.[18] / modified covariance 体系 | **p ≤ (2/3)N** 为相对数据长度 N 的非奇异性/阶数**上限**，非“固定取 2N/3” |

因此 `OUR_MOVING_SAMPLE_INTERPRETATION` 改称：

**`MAX_ORDER_BOUND_SENSITIVITY`**（非“论文可能主实现”；禁止因结果好而采用）。

## 结论

Liang Fig.3/4 实际 simulation p **未由 primary + Ref.[18] 恢复**。

→ **`PAPER_AR_ORDER_NOT_RECOVERABLE`**（合法科学结论）

禁止 p sweep / 用 Fig.3 反推最佳 p。
