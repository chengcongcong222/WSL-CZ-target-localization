# YANG2015_BM_DEFINITION

UTC: 2026-09-23T11:30:59.914058+00:00

## 主文 Eq.(5) 印刷定义（视觉锁定）

$$
b_m = \frac{2 e^{-\alpha_m r_0}}{\alpha_m k_m} \sinh\left[\frac{\alpha_m \Delta R}{2}\right] \phi_m(z_s)
$$

- $r_0 \equiv (r_2+r_1)/2$
- $\Delta R = r_2-r_1 = vT$（合成孔径径向跨度）
- $b_m$ 为 **实数**（论文原文：is a real number）

## 因子分解核对

| 因子 | 是否在 b_m 中 |
| --- | --- |
| $\phi_m(z_s)$ | **YES** |
| $\alpha_m$ | **YES**（$e^{-\alpha_m r_0}$ 与 $\sinh(\alpha_m \Delta R/2)$） |
| $k_m$ | **YES**（分母 $\alpha_m k_m$） |
| $r_0$ | **YES** |
| $\Delta R$ | **YES** |
| source-level / 源幅度 | **NO** |
| $1/\sqrt{k_m}$ 独立因子 | **NO**（并入 $1/k_m$，来自 $\sqrt{k_r k_m}|_{k_r=k_m}$） |
| $\phi_m(z_r)$ | **NO**（在 $g=b_m\phi_m(z_r)$ 外侧） |
| range spreading $\sqrt{r}$ | 已在 Eq.(1) 的 S(r)~sqrt(r) 中处理，不在 b_m |

## 与 delta_m（Eq.9）关系

$$
\delta_m = \frac{2 e^{-\alpha_m r_0}}{\alpha_m k_m}\sinh[(\alpha_m \Delta R)/2] = b_m / \phi_m(z_s)
$$

即 **mode shading coefficient**（Ref.13）。

## 禁止

不得再用 $A_m=\phi_m(z_s)\phi_m(z_r)$ 代替 Yang 的 $b_m$。
$A_m$ 仅是去掉 $k_m,\alpha_m,r_0,\Delta R$ 因子后的简化深度乘积，**不是** Eq.(5) 系数。

## Erratum（carry-forward）

$\delta$ 初始距离偏移时：$g(k_m,z_r)=b_m\phi_m(z_r)e^{i k_m \delta}$。

**PAPER_REPRODUCTION_BASELINE: delta = 0**（本轮不重核 Erratum 全文）。
