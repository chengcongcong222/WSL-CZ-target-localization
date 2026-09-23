# YANG2015_EQUATION_LOCK

UTC: 2026-09-23T11:30:59.914058+00:00
PRIMARY SHA256: `b34952ead13fdaabdf9f1e7d87728fc2cb72fd09497240aac20ebdfc47f4b5b4`

视觉核对页：PDF p.3 (JASA 1679) Eq.(1)-(6)；p.5 Eq.(7)；p.7 Eq.(8)-(9)；p.8 Eq.(10)。

## Eq.(1)  ·  1679 (PDF p.3)

```latex
g(k_r,z_r)=\frac{e^{i\pi/4}}{\sqrt{2\pi k_r}}\int_{r_1}^{r_2} p(r,z_r)\,e^{i k_r r}\,S(r)\,dr,\quad k_r r_1\gg 1
```

- 符号：g(k_r,z_r): SA modal beam output / wavenumber spectrum; k_r: steering wavenumber; z_r: receiver depth; p(r,z_r): block pressure at range r; S(r): spreading-loss shading; r1,r2: aperture ends; e^{i\pi/4}: cylindrical-wave prefactor
- 输入：p(r_n,z_r) on radial samples, S(r), k_r grid
- 输出：complex spectrum g(k_r,z_r)
- 角色：synthetic-aperture modal beamforming (generalized Hankel if S=sqrt(r))
- 单水听器：YES (core) · VLA：YES (per hydrophone, then summed)
- 实现：YES

## Eq.(2)  ·  1679 (PDF p.3)

```latex
S(r)=\langle |p(r,z_r)|^2\rangle^{-1/2}
```

- 符号：\langle\cdot\rangle: range averaging/smoothing; |p|^2: intensity
- 输入：pressure intensity vs range
- 输出：shading S(r) compensating spreading loss (~sqrt(r))
- 角色：weighting in Eq.(1) integrand (NUMERATOR of integrand, multiplies p*exp(i k_r r))
- 单水听器：YES · VLA：replaced by Eq.(8)
- 实现：YES

## Eq.(3)  ·  1679 (PDF p.3)

```latex
p(r,z_r)=\sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r}}\,\phi_m(z_s)\phi_m(z_r)\,e^{-i k_m r-\alpha_m r-i\pi/4}
```

- 符号：phi_m(z): mode depth function; k_m: mode wavenumber; alpha_m: mode attenuation (>0); z_s: source depth; M: number of modes; spatial phase e^{-i k_m r} (as printed)
- 输入：modal parameters, r, z_s, z_r
- 输出：normal-mode pressure (range-independent)
- 角色：propagator used to derive Eq.(4)/(5)
- 单水听器：derivation / oracle peaks · VLA：same
- 实现：for algebra check and ORACLE mode peaks

## Eq.(4)  ·  1679 (PDF p.3)

```latex
g(k_r,z_r)\approx\sum_{m=1}^{M}\frac{\phi_m(z_s)\phi_m(z_r)}{\sqrt{k_r k_m}}\int_{r_1}^{r_2} e^{i(k_r-k_m)r-\alpha_m r}\,dr=\sum_{m=1}^{M} a_m\frac{\phi_m(z_s)\phi_m(z_r)}{k_r-k_m-i\alpha_m}
```

- 符号：a_m=[e^{i(k_r-k_m)r_2-\alpha_m r_2}-e^{i(k_r-k_m)r_1-\alpha_m r_1}]/(i k_m) (AS PRINTED)
- 输入：Eq.(3) into Eq.(1) with S~sqrt(r)
- 输出：spectrum as sum of modal Lorentzians
- 角色：closed form after substitution; peaks at k_r=k_m
- 单水听器：YES (equivalent form) · VLA：YES
- 实现：YES - implement FIRST FORM (integral) or Eq.(5) peak form; see algebra check

## Eq.(5)  ·  1679 (PDF p.3)

```latex
g(k_m,z_r)=b_m\,\phi_m(z_r)
```

- 符号：b_m=(2 e^{-\alpha_m r_0})/(\alpha_m k_m)\sinh[(\alpha_m \Delta R)/2]\,\phi_m(z_s); r_0=(r_2+r_1)/2; \Delta R=r_2-r_1=vT
- 输入：spectrum value at modal peak
- 输出：modal peak amplitude b_m phi_m(z_r)
- 角色：peak relation; b_m carries z_s through phi_m(z_s) and range factors
- 单水听器：YES (feeds Eq.6) · VLA：YES (Eq.9 with delta_m)
- 实现：YES

## Eq.(6)  ·  1679 (PDF p.3)

```latex
D(z)=\left|\sum_{m=1}^{M}\phi_m(z)\,\frac{g(k_m,z_r)}{\bar{\phi}_m(z_r)}\right|^2,\quad \bar{\phi}_m^{-1}(z_r)=\frac{\phi_m(z_r)}{\phi_m^2(z_r)+\Delta^2}
```

- 符号：D(z): source-depth ambiguity; \bar{\phi}_m^{-1}: regularized inverse mode depth function at z_r; \Delta: empirical regularizer (order 0.1*max|phi_m|)
- 输入：g(k_m,z_r) at identified modal peaks; nominal phi_m(.)
- 输出：depth ambiguity D(z) on search grid
- 角色：single-hydrophone source-depth estimator (Bartlett class)
- 单水听器：YES (core) · VLA：NO (uses Eq.9 instead)
- 实现：YES

## Eq.(7)  ·  1681 (PDF p.5)

```latex
\tilde{\theta}(n+1)=\tilde{\theta}(n)+K_1 u(n)+K_2\sum_{m=0}^{n} u(m)
```

- 符号：theta~(n): PLL output phase; u(m): phase error data-oscillator; K1=0.001; K2=0.0001 (2nd-order PLL)
- 输入：narrowband complex time series
- 输出：tracked instantaneous phase theta~(t)
- 角色：remove fast random phase (jitter + medium)
- 单水听器：YES for real data · VLA：YES for real data
- 实现：YES for real-data path; optional for pure simulation

## Eq.(8)  ·  1683 (PDF p.7)

```latex
S(r)=\left\langle \sum_{j=1}^{N}|p(r,z_j)|^2\,\Delta z\right\rangle^{-1/2}
```

- 符号：z_j: VLA hydrophone depths; \Delta z: phone spacing; N: #phones
- 输入：VLA pressure vs range
- 输出：depth-integrated intensity shading
- 角色：VLA replacement of Eq.(2); S~sqrt(r) ignoring alpha_m differences
- 单水听器：NO · VLA：YES
- 实现：for VLA path only

## Eq.(9)  ·  1683 (PDF p.7)

```latex
D(z_j)=\left|\sum_{m=1}^{M} g(k_m,z_j)\right|^2=\left|\sum_{m=1}^{M}\phi_m(z_j)\phi_m(z_s)\delta_m\right|^2,\quad j=1,\ldots,N
```

- 符号：delta_m=(2 e^{-\alpha_m r_0})/(\alpha_m k_m)\sinh[(\alpha_m \Delta R)/2] (mode shading coefficient); note delta_m = b_m/phi_m(z_s)
- 输入：per-phone SA spectra at modal peaks g(k_m,z_j)
- 输出：beam power vs hydrophone depth = depth ambiguity
- 角色：VLA model-independent depth estimation (no phi_m needed experimentally)
- 单水听器：NO · VLA：YES
- 实现：for VLA path only

## Eq.(10)  ·  1684 (PDF p.8)

```latex
D_{\mathrm{MMP}}(r=r_s,z)=\left|\sum_{m=1}^{M}\phi_m(z)\phi_m(z_s)\delta_m'\right|^2,\quad \delta_m'\equiv e^{-2\alpha_m r_s}/k_m
```

- 符号：r_s: true source range; delta_m': MMP shading; related D_MMP(r,z)=(r/2pi)|sum[A_m^{rplc}]^* A_m^{data}|^2
- 输入：matched-mode amplitudes
- 输出：MMP range-depth ambiguity at r=r_s
- 角色：relation to matched-mode processing (comparison only)
- 单水听器：NO · VLA：comparison
- 实现：NO for Yang SA depth (documentation only)
