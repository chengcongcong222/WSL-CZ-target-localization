# YANG2015_IMPLEMENTATION_SPEC

UTC: 2026-09-23T11:30:59.914058+00:00

**本轮禁止运行。** 仅伪代码规范。

## INPUT

- moving CW pressure blocks p(r_n,z_r)（或原始时域走 A1-A4）
- 可选：nominal 环境 -> phi_m(z)
- 可选：true k_m（ORACLE_MODE_ID）

## 单水听器主路径

```
# 0) 压力块（若从时域）
raw time series
  -> [A2] Doppler / range increments r_n = r1 + n * dr
  -> [Eq.7] PLL + low-pass  (real data)
  -> [A3/A4] block pressure p(r_n, z_r)

# 1) shading
S(r_n) = <|p|^2>^{-1/2}          # Eq.(2)

# 2) SA spectrum
g(k_r, z_r) = e^{i pi/4}/sqrt(2 pi k_r)
              * sum_n p(r_n,z_r) e^{i k_r r_n} S(r_n) dr   # Eq.(1)
# 实现等价：Eq.(4) 第一形式积分，禁止未校正的 Eq.(4) 第二形式

# 3) modal peaks
detect peaks in |g(k_r,.)|
mode identification:
  ORACLE_MODE_ID     = match true k_m          # paper simulation
  PRACTICAL_MODE_ID  = spacing + sim guide     # real-data style
g_m <- complex peak value g(k_m, z_r)

# 4) regularized inverse at receiver
phi_m(z_r);  phi_bar_inv = phi_m(z_r) / (phi_m(z_r)^2 + Delta^2)
# Delta = EMPIRICAL_NOT_UNIQUE ~ 0.1 * max|phi_m|  (order-of-magnitude only)

# 5) depth ambiguity Eq.(6)
for z in search_grid:
    D(z) = | sum_m phi_m(z) * g_m * phi_bar_inv_m |^2

# 6) normalized distribution (paper display convention)
Dn(z) = D(z) / sum_z D(z)

# 7) estimate
z_hat = argmax_z D(z)    # or Dn; same argmax
```

## VLA 路径（Sec.IV B）

```
S(r) = < sum_j |p(r,z_j)|^2 dz >^{-1/2}     # Eq.(8)
per phone j: g_j = Eq.(1)
D(z_j) = | sum_m g(k_m, z_j) |^2            # Eq.(9)
P(z_j|z_s) = K D(z_j), K = [sum_j D(z_j)]^{-1}
z_hat = argmax_j D(z_j)
# 不需要 phi_m / SSP
```

## 明确禁止（本轮）

- 运行任何 Yang / E-STD 仿真
- 输出 z_hat / FWHM / PSL / 性能 pass-fail
- delta 扫描、Doppler 误差、MC、AR、P5
- 关闭任何候选算法
- 使用 A_m=phi_s phi_r 或 sqrt(member_energy) 作为 Yang score
