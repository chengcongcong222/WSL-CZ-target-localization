# MMAR_EQUATION_LOCK

UTC: 2026-09-24T09:03:08.658004+00:00

视觉+文本双通道；冲突以 PDF 视觉为准。

## Eq.(1) (P2)

```latex
p_l(r_i,z_r)=\sqrt{2\pi}e^{-j(\pi/4)}\sum_{m=1}^M \phi_m(z_s)\phi_m(z_r)\frac{\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\}}{\sqrt{k_m[r_i+ld\sin\theta_i]}}
```

HLA element normal-mode field with path r_i+ld sinθ_i

## Eq.(2) (P3)

```latex
p_l(r_i,z_r)=\sum_m A_m\exp\{-j(k_m-j\alpha_m)[r_i+ld\sin\theta_i]\},\quad A_m=\sqrt{2\pi}e^{-j\pi/4}\frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_m r_i}}
```

far-field r_i+ld sinθ≈r_i; modal amp A_m

## Eq.(3) (P3)

```latex
B(\hat\theta_i)=\frac{1}{2L+1}\sum_{l=-L}^{L} e^{j k l d \sin\hat\theta_i}\, p_l(r_i,z_r)
```

conventional HLA beamforming weights

## Eq.(4) (P3)

```latex
B(\hat\theta_i)=\frac{1}{2L+1}\sum_m A_m e^{-jk_m r_i-\alpha_m r_i}\frac{\sin[(L+1/2)d X_m]}{\sin((d/2)X_m)},\ X_m=-(k_m-j\alpha_m)\sin\theta_i+k\sin\hat\theta_i
```

array factor; BEAM_FACTOR_FROM_EQ4 = sin[(L+1/2)dX]/sin((d/2)X)

## Eq.(5) (P3)

```latex
p(r;z_s,z_r)=\int_0^\infty g(k_r)J_0(k_r r)k_r dk_r,\quad g(k_r)=\int_0^\infty p(r)J_0(k_r r) r dr
```

Hankel pair

## Eq.(6) (P3)

```latex
g(k_r;z_s,z_r)\sim\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{-\infty}^{+\infty}p(r;z_s,z_r)e^{i k_r r}\sqrt{r}\,dr,\ k_r r\gg1
```

FT approx of Hankel

## Eq.(7) (P3)

```latex
g(k_r,z_r)=\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{r_0}^{r_0+R} B(r)\,e^{i k_r r}\,S(r)\,dr,\ k_r r_0\gg1
```

generalized Hankel on beam output

## Eq.(8) (P3)

```latex
S(r)=\langle |B(r)|^2\rangle^{-1/2}
```

data spreading compensation

## Eq.(9) (P3)

```latex
g\sim\sum_m \frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}\sin b(X_m)\int_{r_0}^{r_0+R}e^{j(k_r-k_m)r-\alpha_m r}dr=\sum_m a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m+j\alpha_m}
```

closed form after Eq4+Eq7; integral kernel exp[j(kr-km)r - αm r]

## Eq.(10) (P3)

```latex
a_m=\frac{e^{j[(k_r-k_m)-\alpha_m](r_0+R)}-e^{j[(k_r-k_m)-\alpha_m]r_0}}{j\sqrt{k_r k_m}\,\sin b(X_m)}
```

PRINTED a_m — algebra audit flagged

## Eq.(11) (P3)

```latex
g(k_m,z_r)\sim b_m\phi_m(z_r)
```

peak at k_r=k_m

## Eq.(12) (P3)

```latex
b_m=\frac{2e^{-\alpha_0 r'}}{\alpha_m k_m}\sinh\left(\frac{\alpha_m R}{2}\right)\phi_m(z_s)\sin b(X_m),\ r'=r_0+R/2
```

PRINTED b_m with e^{-α0 r'} — α0 subscript flagged

## Eq.(13) (P3)

```latex
\mathbf{g}=\Phi\cdot\mathbf{b}
```

matrix form

## Eq.(14) (P3)

```latex
\mathbf{g}=[g(k_1,z_r),\ldots,g(k_M,z_r)]^T
```

spectrum vector

## Eq.(15) (P3)

```latex
\Phi=\mathrm{diag}([\phi_1(z_r),\ldots,\phi_M(z_r)])
```

diag mode functions at zr

## Eq.(16) (P3)

```latex
\mathbf{b}=[b_1,\ldots,b_M]^T
```

modal coefficients

## Eq.(17) (P4)

```latex
y[i]=B(r_i)S(r_i),\quad i=1,2,\ldots,2L+1
```

AR input = beamformed × S; SAMPLE COUNT CONFLICT

## Eq.(18) (P4)

```latex
y[i]=-\sum_{k=1}^{p}a[k]y[i-k]+u[i],\quad p\ \mathrm{often}=(2/3)(2L+1)
```

AR model; order notation CONFLICT

## Eq.(19) (P4)

```latex
P_{AR}(l)=\frac{\sigma^2}{\left|1+\sum_{k=1}^{p}a[k]\exp[-ilk]\right|}
```

AR spectrum; l→k_r mapping NOT_EXPLICIT

## Eq.(20) (P4)

```latex
D(z)=\varphi(z)\,\mathbf{b}\mathbf{b}^H\varphi^H(z)
```

matched-mode depth function

## Eq.(21) (P4)

```latex
\varphi(z)=[\phi_1(z),\phi_2(z),\ldots,\phi_M(z)]
```

replica row

## Eq.(22) (P4)

```latex
\mathbf{b}=(\Phi+U)^{-1}\mathbf{g}
```

regularized solve

## Eq.(23) (P4)

```latex
U=\mathrm{diag}\left(\left[\frac{\Delta}{\phi_1(z_r)},\frac{\Delta}{\phi_2(z_r)},\ldots,\frac{\Delta}{\phi_M(z_r)}\right]\right),\ \Delta\sim\tfrac12\max|\phi|
```

Liang regularization — NOT Yang 0.1

## Eq.(24) (P4)

```latex
\min_{k_0}(\mathbf{k}-\mathbf{k}_0)^H(\mathbf{k}-\mathbf{k}_0)\ \mathrm{s.t.}\ k_0(1)<k_0(2)<\cdots<k_0(M_0)
```

ORDERED_SUBSET_GLOBAL_MATCH

## Eq.(25) (P4)

```latex
\mathbf{k}=[k_1,\ldots,k_{M_0}]^T
```

estimated wavenumbers

## Eq.(26) (P4)

```latex
\mathbf{k}'=[k'_1,\ldots,k'_M]^T
```

KRAKEN/model wavenumbers

## Eq.(27) (P4)

```latex
\mathbf{k}_0\subseteq\mathbf{k}'
```

ordered subset of model k'

## Eq.(28) (P5)

```latex
\mathrm{SNR}=10\lg\frac{P_s}{P_n}\Big|_{r=r_0}
```

SNR at r0

## Eq.(29) (P5+)

```latex
(\text{correct-depth probability — see paper Sec.4 MC})\ \ | \hat z-z_\mathrm{true}|\le 5\ \mathrm{m}
```

MC success criterion |err|<=5m; C0=500

## Eq.(30) (P5+)

```latex
(\text{CI for probability — binomial CI as printed in paper})
```

confidence interval of correct-depth rate
