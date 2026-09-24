# EQ24_ORDERING_AUDIT

UTC: 2026-09-24T16:02:27.087144+00:00

Liang Eq.(24): k0(1)<…<k0(M0)。

修正：进入 DP 前 **k̂ 与 k' 均按物理 k 升序**；模态 ID 经原始 KRAKEN index 映射回。

SAB 与 MMAR 使用同一规则。KRAKEN 文件内 mode 存储序为 k **降序**，**不得**直接当 Eq24 顺序。

重算 R0/R1 全部 metrics → `MODE_ORDER_METRICS_CORRECTED.csv`。
