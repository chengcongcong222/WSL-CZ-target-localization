# MMAR_MODE_ORDER_MATCHING

UTC: 2026-09-24T09:03:08.658004+00:00

Eq.(24)–(27)：

- 估计向量 k=[k_1…k_{M0}]^T
- 模型 KRAKEN k'=[k'_1…k'_M]^T
- 从 k' 中选 **有序子集** k0，满足 k0(1)<…<k0(M0)
- 最小化 (k−k0)^H (k−k0)

类型：**ORDERED_SUBSET_GLOBAL_MATCH**

不是：independent nearest neighbor / Hungarian / C2.3A Rayleigh unique-only。

## 对 C2.3A 的含义

MMAR **不要求** 每个 AR 峰在 Rayleigh cell 内唯一 mode。
它用高分辨 k̂ + 全局有序子集匹配定 mode order。

后续指标应是 **MODE_ORDER_RECOVERY_RATE**，不是 N_UNIQUE_RAYLEIGH_MODES。
