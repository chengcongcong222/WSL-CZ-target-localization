# R3-C2.2D 论文单因素终审

UTC: 2026-09-24T05:07:24.472480+00:00

## 冻结

基础工具全部冻结（parser / flp / shd / alpha / Eq1/3/5/6 / env / r / dr / delta=0 / M=9999）。

**incremental 措辞**：仅 **m=2** 做过增量检查且 PASS；聚合 M=1/2/4/8/16/23 均 corr≈1、residual≈5e-4；M=1 的 c_ref 在 M=23 成立。已足够支持 `FIELD_EQ3_MULTIMODE_VALIDATED`。**不**写“m=2–8 全过”。

## 核心事实（zs=50, zr=18）

| Delta | z_hat | D50/Dmax | D50/D43 |
| --- | --- | --- | --- |
| 0.05 | 43 | 0.59 | — |
| 0.10 | 43 | **0.86** | ≈1 |
| 0.20 | **50** | 1.00 | 43 为第二峰 |

真值峰未消失，而是与 43 m **竞争排序**。

## 判定

### `C2_2D_REGULARIZER_SENSITIVE`

zs=50,zr=18 的 43↔50 排序随 Delta 翻转：0.05–0.10 主峰 43（50 为强竞争峰），0.20 主峰 50；zr=70 control 对 Delta 稳定。属于 Eq.(6) 经验正则化下的竞争峰排序，不是 Yang 失败。

机制标签：['C2_2D_FINITE_APERTURE_LEAKAGE', 'C2_2D_REGULARIZER_SENSITIVE', 'C2_2D_SHADING_SENSITIVE', 'PEAK_SELECTION_ID_NOT_PRIMARY_CAUSE', 'RECEIVER_DEPTH_REGULARIZATION_SENSITIVITY']

论文阶段冻结：**`C2_PAPER_METHOD_MECHANISM_REPRODUCED`**

## 范围

不再用 2 m gate 定方法真假。浅海 paper 阶段 **到此结束**。
下一阶段必须是 **R3-C2.3 E-STD 迁移**（本轮不执行）。
YANG_ROUTE_UNDECIDED。
