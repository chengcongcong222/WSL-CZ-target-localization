# MMAR_EQUATION_LOCK_V2

UTC: 2026-09-24T09:16:56.736527+00:00
PRIMARY SHA256: `9364914d34f5408f80065a31cd7fbcb5195c21047c97b2eb118348191b861834`

## 修正记录（相对 V1）

| 项 | V1 错误 | V2 正确 |
| --- | --- | --- |
| Eq.(10) | sin b(X_m) 误放在**分母** | **乘在分式之后** |
| sin b(X_m) | 似为未定义函数 | `SINB_NOTATION_UNDEFINED_USE_EQ4_BEAM_FACTOR` |
| Eq.(29)(30) | 文字占位 | 按 PDF 补全 |
| HLA d | d_m=0.5 m | d=λ/2，数值未由正文钉死 |

## Eq.(1)–(4) HLA

```latex
(1) p_l(r_i,z_r)=\sqrt{2\pi}e^{-j\pi/4}\sum_m \phi_m(z_s)\phi_m(z_r)\frac{\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\}}{\sqrt{k_m[r_i+ld\sin\theta_i]}}
(2) p_l=\sum_m A_m\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\},\ A_m=\sqrt{2\pi}e^{-j\pi/4}\phi_m(z_s)\phi_m(z_r)/\sqrt{k_m r_i}
(3) B(\hat\theta_i)=\frac{1}{2L+1}\sum_{l=-L}^{L}e^{j k l d\sin\hat\theta_i}p_l(r_i,z_r)
(4) B=\frac{1}{2L+1}\sum_m A_m e^{-jk_m r_i-\alpha_m r_i}\underbrace{\frac{\sin[(L+1/2)dX_m]}{\sin((d/2)X_m)}}_{\mathrm{array\ factor\ (Eq.4)}},\
X_m=-(k_m-j\alpha_m)\sin\theta_i+k\sin\hat\theta_i
```

**BF_m（实现用，源自 Eq.4）**：

```latex
BF_m=\frac{1}{2L+1}\cdot\frac{\sin[(L+1/2)dX_m]}{\sin((d/2)X_m]}
```

论文后文的 `sin b(X_m)` **无独立定义** → 使用 BF_m，禁止自造 `sinb()`。

## Eq.(5)–(8)

```latex
(5) Hankel pair p \leftrightarrow g(k_r)
(6) g(k_r)\sim\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int p e^{ik_r r}\sqrt{r}\,dr
(7) g(k_r,z_r)=\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{r_0}^{r_0+R}B(r)e^{ik_r r}S(r)\,dr
(8) S(r)=\langle|B(r)|^2\rangle^{-1/2}
```

## Eq.(9)（积分形式 — 实现主依据）

```latex
g\sim\sum_m \frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}\,BF_m
\int_{r_0}^{r_0+R} e^{j(k_r-k_m)r-\alpha_m r}\,dr
=\sum_m a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m+j\alpha_m}
```

积分核：**exp[j(k_r−k_m)r − α_m r]**（α 在实部衰减项，不在 j[…] 内）。

## Eq.(10) — **印刷原样（乘法因子，非分母）**

```latex
a_m=\frac{ e^{j[(k_r-k_m)-\alpha_m](r_0+R)}-e^{j[(k_r-k_m)-\alpha_m]r_0} }{j\sqrt{k_r k_m}}\;\sin b(X_m)
```

即：

```latex
a_m=\Big[\cdots\Big]/\big(j\sqrt{k_rk_m}\big)\ \times\ \sin b(X_m)
```

**不是** `/ [j√(k_r k_m) sin b(X_m)]`。

仍保留 `EQ10_PRINTED_FORM_INTERNAL_INCONSISTENCY`：α_m 被印进 j[(k_r−k_m)−α_m]。

## Eq.(11)–(16)

```latex
(11) g(k_m,z_r)\sim b_m\phi_m(z_r)
(12) PRINTED: b_m=\frac{2e^{-\alpha_0 r'}}{\alpha_m k_m}\sinh(\alpha_m R/2)\phi_m(z_s)\sin b(X_m),\ r'=r_0+R/2
(13)–(16) g=\Phi b,\ \Phi=\mathrm{diag}(\phi_m(z_r)),\ b=[b_m]
```

Eq.(12) `e^{-α_0 r'}` → `EQ12_ALPHA_SUBSCRIPT_CONFLICT`。

## FROM_EQ9_DERIVED_BM（实现闭式）

令 Δ_k=k_r−k_m，

```latex
I_m=\int_{r_0}^{r_0+R}e^{(j\Delta_k-\alpha_m)r}dr
=\frac{e^{(j\Delta_k-\alpha_m)(r_0+R)}-e^{(j\Delta_k-\alpha_m)r_0}}{j\Delta_k-\alpha_m}
```

k_r=k_m 时：

```latex
I_m=\frac{2e^{-\alpha_m r'}}{\alpha_m}\sinh\left(\frac{\alpha_m R}{2}\right),\quad r'=r_0+R/2
```

故：

```latex
b_m^{\mathrm{impl}}=\frac{2e^{-\alpha_m r'}}{\alpha_m k_m}\sinh\left(\frac{\alpha_m R}{2}\right)\phi_m(z_s)\,BF_m
```

标 **`FROM_EQ9_DERIVED_BM`**；**不用 printed α_0**。实现 α→0 稳定极限：sinh(x)/x 形式。

## Eq.(17)–(19) AR

```latex
(17) y[i]=B(r_i)S(r_i),\ i=1,\ldots,2L+1
(18) y[i]=-\sum_{k=1}^{p}a[k]y[i-k]+u[i],\ p\ \mathrm{often}\ (2/3)(2L+1)
(19) PRINTED: P_{AR}(l)=\frac{\sigma^2}{\left|1+\sum_{k=1}^{p}a[k]e^{-ilk}\right|}
```

`EQ19_DENOMINATOR_POWER_NOT_RESOLVED`：印刷 **|A(l)|^1**；DSP 常规 **|A(l)|^{-2}** 仅作 STANDARD_AR_PSD_CONTROL，非原文。

## Eq.(20)–(23) 深度函数

```latex
(20) D(z)=\varphi(z)\,b b^H\,\varphi^H(z)
(21) \varphi(z)=[\phi_1(z),\ldots,\phi_M(z)]
(22) b=(\Phi+U)^{-1}g
(23) U=\mathrm{diag}([\Delta/\phi_1(z_r),\ldots,\Delta/\phi_M(z_r)]),\
\Delta\ \text{on the order of one-half of max mode function}
```

**LIANG_PAPER_REGULARIZER Δ≈½ max|φ|**（≠ Yang 0.1）。

## Eq.(24)–(27) ordered matching

```latex
(24) \min (k-k_0)^H(k-k_0)\ \mathrm{s.t.}\ k_0(1)<\cdots<k_0(M_0)
(25) k=[k_1,\ldots,k_{M_0}]^T
(26) k'=[k'_1,\ldots,k'_M]^T
(27) k_0\subseteq k'
```

`PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`；若自写 DP：**`OUR_EXACT_SOLVER_FOR_EQ24`**。

## Eq.(28)–(30)

```latex
(28) SNR=10\lg(P_s/P_n)|_{r=r_0}
(29) P=C/C_0,\quad C_0=500,\quad \mathrm{correct}:|\hat z-z_\mathrm{true}|\le 5\ \mathrm{m}
(30) 90\%\ \mathrm{CI}: \ P \pm 1.645\sqrt{P(1-P)/C_0}
```
