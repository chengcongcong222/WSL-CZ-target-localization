# MMAR_EQUATION_LOCK_FINAL

UTC: 2026-09-24T09:51:45.113589+00:00
PRIMARY SHA256: `9364914d34f5408f80065a31cd7fbcb5195c21047c97b2eb118348191b861834`
SUPERSEDES: R3_C2_4A_PRIMARY (V1), R3_C2_4A_PRIMARY_FIX (V2)

状态标记：`PRINTED_PRIMARY` / `PRIMARY_CONFLICT` / `DERIVED_IMPLEMENTATION`

---

## Eq.(1) `PRINTED_PRIMARY`

```latex
p_l(r_i,z_r)=\sqrt{2\pi}\,e^{-j(\pi/4)}\sum_{m=1}^{M}
\phi_m(z_s)\phi_m(z_r)
\frac{\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\}}{\sqrt{k_m[r_i+ld\sin\theta_i]}}
```

## Eq.(2) `PRINTED_PRIMARY`

```latex
p_l(r_i,z_r)=\sum_{m=1}^{M} A_m\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\},
\quad A_m=\sqrt{2\pi}e^{-j\pi/4}\frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_m r_i}}
```

（far-field：r_i+ld\sin\theta_i\approx r_i）

## Eq.(3) `PRINTED_PRIMARY`

```latex
B(\hat\theta_i)=\frac{1}{2L+1}\sum_{l=-L}^{L} e^{j k l d \sin\hat\theta_i}\,p_l(r_i,z_r)
```

## Eq.(4) `PRINTED_PRIMARY`

```latex
B(\hat\theta_i)=\frac{1}{2L+1}\sum_{m=1}^{M} A_m e^{-jk_m r_i-\alpha_m r_i}
\frac{\sin[(L+1/2)d X_m]}{\sin[(d/2) X_m]},
\quad X_m=-(k_m-j\alpha_m)\sin\theta_i+k\sin\hat\theta_i
```

**BF_m（阵因子，Eq.4 原式）**：

```latex
BF_m=\frac{1}{2L+1}\cdot\frac{\sin[(L+1/2)d X_m]}{\sin[(d/2) X_m]}
```

（分母书写必须为 `\sin[(d/2)X_m]`，不得 `sin((d/2)X_m]`。）

## Eq.(5) `PRINTED_PRIMARY`

```latex
p(r;z_s,z_r)=\int_0^\infty g(k_r;z_s,z_r)J_0(k_r r)k_r\,dk_r,
\quad g(k_r;z_s,z_r)=\int_0^\infty p(r;z_s,z_r)J_0(k_r r)r\,dr
```

## Eq.(6) `PRINTED_PRIMARY`

```latex
g(k_r;z_s,z_r)\sim\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{-\infty}^{+\infty} p(r;z_s,z_r)e^{i k_r r}\sqrt{r}\,dr,
\quad k_r r\gg 1
```

## Eq.(7) `PRINTED_PRIMARY`

```latex
g(k_r,z_r)=\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{r_0}^{r_0+R} B(r)\,e^{i k_r r}\,S(r)\,dr,
\quad k_r r_0\gg 1
```

## Eq.(8) `PRINTED_PRIMARY`

```latex
S(r)=\langle |B(r)|^2\rangle^{-1/2}
```

## Eq.(9) `PRINTED_PRIMARY`（严格原文印刷）

```latex
g(k_r,z_r)\sim\sum_{m=1}^{M}
\frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}\,\sin b(X_m)
\int_{r_0}^{r_0+R} e^{j(k_r-k_m)r-\alpha_m r}\,dr
=\sum_{m=1}^{M} a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m+j\alpha_m}
```

积分核：**exp[j(k_r−k_m)r − α_m r]**。此处 **sin b(X_m) 为论文印刷原样**，未替换为 BF_m / mathcal B_m。

---

## Eq.(9) `DERIVED_IMPLEMENTATION`（单独解释层）

```latex
\sin b(X_m)\ \Rightarrow\ BF_m=\frac{1}{2L+1}\cdot\frac{\sin[(L+1/2)d X_m]}{\sin[(d/2) X_m]}
```

标记：**`DERIVED_EQ4_BEAM_FACTOR_INTERPRETATION`**（由 Eq.(4)→Eq.(7) 代数，非论文明文定义 sin b）。

实现公式（`FROM_EQ9_DERIVED_BM`，与 PRINTED_PRIMARY 分离）：

```latex
b_m^{\mathrm{impl}}=\frac{2e^{-\alpha_m r'}}{\alpha_m k_m}\sinh\left(\frac{\alpha_m R}{2}\right)\phi_m(z_s)\,BF_m,
\quad r'=r_0+R/2
```

## Eq.(10) `PRIMARY_CONFLICT`（印刷原样）

```latex
a_m=\frac{ e^{j[(k_r-k_m)-\alpha_m](r_0+R)}-e^{j[(k_r-k_m)-\alpha_m]r_0} }{j\sqrt{k_r k_m}}\;\sin b(X_m)
```

- **sin b(X_m) 为乘法因子，不在分母**（相对 V1 转录已修正）
- 仍冲突：α_m 印入 j[(k_r−k_m)−α_m]，与 Eq.(9) 积分核不符

## Eq.(11) `PRINTED_PRIMARY`

```latex
g(k_m,z_r)\sim b_m\phi_m(z_r)
```

## Eq.(12) `PRIMARY_CONFLICT`（印刷原样）

```latex
b_m=\frac{2 e^{-\alpha_0 r'}}{\alpha_m k_m}\sinh\left(\frac{\alpha_m R}{2}\right)\phi_m(z_s)\sin b(X_m),
\quad r'=r_0+R/2
```

冲突：`e^{-α_0 r'}` 下标 α_0（应为 mode 的 α_m）。

## Eq.(13) `PRINTED_PRIMARY`

```latex
\mathbf{g}=\Phi\cdot\mathbf{b}
```

## Eq.(14) `PRINTED_PRIMARY`

```latex
\mathbf{g}=[g(k_1,z_r),g(k_2,z_r),\ldots,g(k_M,z_r)]^T
```

## Eq.(15) `PRINTED_PRIMARY`

```latex
\Phi=\mathrm{diag}([\phi_1(z_r),\phi_2(z_r),\ldots,\phi_M(z_r)])
```

## Eq.(16) `PRINTED_PRIMARY`

```latex
\mathbf{b}=[b_1,b_2,\ldots,b_M]^T
```

## Eq.(17) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（样本记号）

```latex
y[i]=B(r_i)S(r_i),\quad i=1,2,\ldots,2L+1
```

冲突：`2L+1` 同时用作 HLA 阵元数与 AR 数据长度（见 AR_INTERPRETATION_GATE）。

## Eq.(18) `PRINTED_PRIMARY` + `PRIMARY_CONFLICT`（阶数）

```latex
y[i]=-\sum_{k=1}^{p} a[k]y[i-k]+u[i],
\quad p\ \mathrm{often\ set\ to}\ (2/3)(2L+1)
```

## Eq.(19) `PRINTED_PRIMARY` + `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`（分母幂次对照，非原文冲突）

```latex
P_{AR}(l)=\frac{\sigma^2}{\left|1+\sum_{k=1}^{p}a[k]\exp[-ilk]\right|}
```

对照：`EQ19_DENOMINATOR_POWER_NOT_RESOLVED` → 记为需控制实验的实现歧义；**非 PRIMARY_CONFLICT**。
- `PRINTED_EQ19_MODULUS_POWER_1`
- `STANDARD_AR_PSD_MODULUS_POWER_2_CONTROL`（DSP 常规，非论文）

## Eq.(20) `PRINTED_PRIMARY`

```latex
D(z)=\varphi(z)\,\mathbf{b}\mathbf{b}^H\,\varphi^H(z)
```

## Eq.(21) `PRINTED_PRIMARY`

```latex
\varphi(z)=[\phi_1(z),\phi_2(z),\ldots,\phi_M(z)]
```

## Eq.(22) `PRINTED_PRIMARY`

```latex
\mathbf{b}=(\Phi+U)^{-1}\mathbf{g}
```

## Eq.(23) `PRINTED_PRIMARY`

```latex
U=\mathrm{diag}\left(\left[\frac{\Delta}{\phi_1(z_r)},\frac{\Delta}{\phi_2(z_r)},\ldots,\frac{\Delta}{\phi_M(z_r)}\right]\right)
```

Δ：**on the order of one-half of the maximum value of the mode function**（`LIANG_PAPER_REGULARIZER`）。

## Eq.(24) `PRINTED_PRIMARY` + `PAPER_DOES_NOT_SPECIFY_NUMERICAL_SOLVER`（实现细节未指定，非公式冲突）

```latex
\min_{k_0}(\mathbf{k}-\mathbf{k}_0)^H(\mathbf{k}-\mathbf{k}_0)
\quad \mathrm{s.t.}\quad k_0(1)<k_0(2)<\cdots<k_0(M_0)
```

## Eq.(25) `PRINTED_PRIMARY`

```latex
\mathbf{k}=[k_1,k_2,\ldots,k_{M_0}]^T
```

## Eq.(26) `PRINTED_PRIMARY`

```latex
\mathbf{k}'=[k'_1,k'_2,\ldots,k'_M]^T
```

## Eq.(27) `PRINTED_PRIMARY`

```latex
\mathbf{k}_0\subseteq\mathbf{k}'
```

## Eq.(28) `PRINTED_PRIMARY`

```latex
\mathrm{SNR}=10\lg\frac{P_s}{P_n}\Big|_{r=r_0}
```

## Eq.(29) `PRINTED_PRIMARY`

```latex
P=\frac{C}{C_0},\quad C_0=500,
\quad \mathrm{correct}:\ |\hat{z}-z_{\mathrm{true}}|\le 5\ \mathrm{m}
```

## Eq.(30) `PRINTED_PRIMARY`

```latex
P\pm 1.645\sqrt{\frac{P(1-P)}{C_0}}\quad (90\%\ \mathrm{CI})
```

---

## DERIVED_IMPLEMENTATION 汇总

| 代号 | 内容 | 来源 |
| --- | --- | --- |
| `DERIVED_EQ4_BEAM_FACTOR_INTERPRETATION` | 后文 `sin b(X_m)` 解释为 Eq.(4) 阵因子 BF_m；**论文未独立定义 sin b** | Eq.(4)→Eq.(7) 代数 |
| `FROM_EQ9_DERIVED_BM` | `b_m=\frac{2e^{-\alpha_m r'}}{\alpha_m k_m}\sinh(\alpha_m R/2)\phi_m(z_s)\,BF_m` | Eq.(9) 积分 |
| `OUR_EXACT_SOLVER_FOR_EQ24` | （契约，本轮不跑）DP 有序子序列 | 我方实现 |
