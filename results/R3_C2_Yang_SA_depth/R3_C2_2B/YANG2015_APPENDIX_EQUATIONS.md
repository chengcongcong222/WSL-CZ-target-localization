# YANG2015_APPENDIX_EQUATIONS

UTC: 2026-09-23T11:30:59.914058+00:00

视觉核对：PDF p.9 / JASA 1685。

## A1-A4

### Eq.(A1)

```latex
p(r,z_0,t)=b(t)\sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r}}\phi_m(z_s)\phi_m(z_0)\times e^{-i(2\pi f_c t - k_m r)-\alpha_m r-i\pi/4}+n(t)
```

- 角色：moving-source narrowband time series (as printed; spatial phase +i k_m r after distributing -i)
- 注记：SPATIAL PHASE SIGN OPPOSITE to Eq.(3) e^{-i k_m r}

### Eq.(A2)

```latex
p(r,z_0,t_{n,l})\simeq e^{i\theta(t_{n,l})}\sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r_n}}\phi_m(z_s)\phi_m(z_0)\times e^{-i(2\pi f_c T_n - k_m r_n)-\alpha_m r_n-i\pi/4}\,e^{-i 2\pi f_D (l\Delta t)}+n(t_{n,l})
```

- 角色：Doppler-shifted block model; f_D=f_c(1+v/c); r_n=v T_n; dr=v dt
- 注记：assumes constant Doppler for dominant modes; |b|~1; e^{-alpha_m (l dr)}~1

### Eq.(A3)

```latex
\tilde{p}(r_n,z_0,f)=\sum_{l=1}^{L} e^{i 2\pi f t_{n,l}}\left(e^{-i\hat{\theta}(t_{n,l})}p(r_{n,l},t_{n,l})\right)= e^{i 2\pi (f-f_c)T_n}\sum_{l=1}^{L} e^{i 2\pi f (l\Delta t)} e^{-i 2\pi f_D (l\Delta t)} e^{i\Delta\theta(t_{n,l})}\times \sum_{m=1}^{M}\sqrt{\frac{2\pi}{k_m r}}\phi_m(z_s)\phi_m(z_0)e^{i k_m r_n-\alpha_m r_n-i\pi/4}
```

- 角色：PLL/FFT block pressure estimator (printed modal factor has e^{+i k_m r_n})
- 注记：SPATIAL PHASE +i k_m r_n as printed; conflicts with Eq.(3)

### Eq.(A4)

```latex
\tilde{p}(r_n,z_0,f_D)= L\, e^{i 2\pi (\Delta f) T_n}\, e^{i\overline{\Delta\theta}(T_n)}\, p(r_n,z_0)
```

- 角色：block-by-block pressure recovery; df=f_D-f_c; mean dtheta; p(r_n,z_0) stated to be Eq.(3)
- 注记：if dtheta=0, correct Doppler phase by e^{-i 2 pi (df) T_n} then use p in Eq.(1)

## 完整观测链

```
原始时域 p(r,z0,t) [A1]
  -> Doppler f_D / dr 分块 [A2]
  -> PLL 相位 theta_hat + 低通 [Sec.III / Eq.7]
  -> 块 FFT 压力 p~(r_n,z0,f) [A3]
  -> Doppler 相位校正 -> p(r_n,z0) [A4]
  -> Eq.(1) 波数谱 g(k_r,z_r)
  -> 模态峰识别 (ORACLE / PRACTICAL)
  -> Eq.(6) D(z) -> 归一化 -> argmax
```

## 相位约定冲突（必读）

- Eq.(3) 印刷：空间相位 **e^{-i k_m r}**
- A1/A2/A3 印刷：空间相位 **e^{+i k_m r}**（由 e^{-i(2 pi f_c t - k_m r)} 与 e^{+i k_m r_n} 可见）
- A4 声称 p(r_n,z0) 即 Eq.(3)

**主方法实现必须采用 Eq.(3) 约定**（与 Eq.(1)(4)(5) 代数闭环）。
若从 A3 直接得到 +i k_m r 的压力，须与 Eq.(1) steering 约定一致（或取共轭）后再进谱估计。
此为 **PAPER_APPENDIX_VS_EQ3_PHASE_CONVENTION**，不阻塞主链锁定。
