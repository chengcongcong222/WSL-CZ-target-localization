# C2_1 freeze + audit-code notes (no KRAKEN rerun)

UTC: 2026-09-23T08:46:33.049990+00:00

## 冻结

- `C2_1_PARSER_VALIDATED`：`.mod` 读取链可靠（用户独立二进制复核一致）
- zgrid 150:2:250 为 KRAKEN 真实深度节点，可合法计算 φ_m(z)
- self-built FEM 保持 **NOT_ADMISSIBLE**
- 283 Hz：`M_mod=902`, `M_prt=903`，保留 `MODE_COUNT_DISCREPANCY`；后续以 `.mod` 结构为准

## 审计代码漏洞（本轮不改数据、不重跑 KRAKEN；仅记录）

1. 非 235 Hz `depth_ok` 含 `(... or True)`，使该项退化
2. 跨 z_s 一致性只比了前两个运行，未覆盖 z_s=220（用户三组两两独立复核 corr=1、Δk=0）
3. zgrid case 的 `k_checksum_ok` / `phi_ok` 有硬编码 True 成分

原则：**DECISION 中的 True 本身不是证据，必须追溯 True 的产生方式。**

本轮按 R3-C2.2A 范围**不**改 parser gate 代码；留给后续独立审计补丁。
