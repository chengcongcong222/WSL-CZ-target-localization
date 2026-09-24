# EQ7_TO_EQ12_ALGEBRA_AUDIT_V2

UTC: 2026-09-24T09:16:56.736527+00:00

## Eq.(9) 积分核

exp[j(k_r−k_m)r − α_m r] → 分母 j(k_r−k_m)−α_m = j(k_r−k_m+jα_m)

与印刷第二形式 (k_r−k_m+jα_m) **一致**。

## Eq.(10) 转录修正

**正确印刷**：sin b(X_m) **乘在** /[j√(k_r k_m)] **之后**，**不在分母**。

```latex
a_m=\frac{e^{j[(k_r-k_m)-\alpha_m](r_0+R)}-e^{j[(k_r-k_m)-\alpha_m]r_0}}{j\sqrt{k_r k_m}}\;\sin b(X_m)
```

## 仍保留内部不一致

`EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY`：指数里 α 进 j[(k_r−k_m)−α_m]=j(k_r−k_m)−jα_m，
与积分核系数 j(k_r−k_m)−α_m 不符。

→ 实现：**FROM_EQ9_DERIVED_BM** / Eq.(9) 积分。

## Eq.(12)

e^{−α_0 r'} → `EQ12_ALPHA_SUBSCRIPT_CONFLICT`；实现用 e^{−α_m r'}。

## sin b(X_m)

`SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR`：

```latex
BF_m=\frac{1}{2L+1}\frac{\sin[(L+1/2)dX_m]}{\sin[(d/2)X_m]}
```
