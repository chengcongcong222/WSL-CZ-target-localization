# YANG2015_REF7_DEPENDENCY

UTC: 2026-09-23T11:30:59.914058+00:00

Yang2015 明确：*simulated data given in Ref.7 for a downward refractive SSP. (The readers are referred to Ref.7 for details.)*

## 完全复现论文模拟仍缺（Yang2015 未给，不得猜）

1. 具体 SSP 数值剖面（Fig.1 of Ref.7）
2. 水深 H
3. 底质声学参数（密度、声速、衰减）
4. 模态求解/衰减模型细节
5. 产生 p(r,z) 的传播配置（采样几何、源级归一等）
6. Fig.1 所用 true k_m 表（ORACLE 标注）

## 允许状态

- `PAPER_REPRO_NEEDS_REF7` <- **当前**
- `PAPER_REPRO_ENV_COMPLETE`（取得 Ref.7 并锁定后）

## 边界

- Eq.(1)-(10)、A1-A4 **不依赖** Ref.7 即可锁定（已完成）
- 论文条件复现数值结果 **依赖** Ref.7
- 即使 NEEDS_REF7，**YANG_ROUTE 仍为 UNDECIDED**，不得 pass
