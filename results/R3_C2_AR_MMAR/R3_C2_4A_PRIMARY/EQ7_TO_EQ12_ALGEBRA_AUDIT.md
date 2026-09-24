# EQ7_TO_EQ12_ALGEBRA_AUDIT

UTC: 2026-09-24T09:03:08.658004+00:00

## Eq.(9) 积分核（印刷）

exp[ j(k_r−k_m) r − α_m r ]

解析积分：

∫ exp[j(k_r−k_m)r − α_m r] dr = Δexp / [ j(k_r−k_m) − α_m ]
= Δexp / [ j( k_r−k_m + j α_m ) ]

故闭式分母应与 **(k_r − k_m + j α_m)** 一致 — Eq.(9) 第二形式分母印刷为 **(k_r−k_m+jα_m)** ✓

## Eq.(10) 印刷形式

a_m = [ e^{ j[(k_r−k_m)−α_m](r0+R) } − e^{ j[(k_r−k_m)−α_m] r0 } ] / ( j√(k_r k_m) sin b(X_m) )

印刷把 α_m 放进 **j[(k_r−k_m)−α_m] = j(k_r−k_m) − j α_m**，
与积分核指数 **j(k_r−k_m)r − α_m r** 的系数 **j(k_r−k_m) − α_m** 不一致。

### 状态：`EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY`

实现只允许 **Eq.(9) 积分形式** 或由 Eq.(9) 重推的闭式。

## Eq.(12) α 下标

印刷：**e^{−α_0 r'}**（α_0 而非 α_m）。
正文/Eq9/Eq11 要求 mode-specific α_m。

### 状态：`EQ12_ALPHA_SUBSCRIPT_CONFLICT`

实现：**FROM_EQ9_DERIVED_BM**，禁止 PRINTED_EQ12_BLIND_COPY。
（对比 Yang2015：b_m 用 e^{−α_m r0}，分母 α_m k_m，sinh(α_m R/2)，另乘 sin b(X_m)。）
