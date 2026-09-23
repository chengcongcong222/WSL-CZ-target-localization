# PAPER_ENV_GEOMETRY_AUDIT

UTC: 2026-09-23T12:44:30.356064+00:00

## Yang2014 Sec.III-C 正文（DIRECTLY_SPECIFIED）

- 热跃层约 10 m：c≈1533 m/s
- 10 m 向下近似线性降至 40 m：c≈1478 m/s
- 40 m 至约 **88 m**：保持 1478 m/s
- bottom: cp≈1650 m/s, rho=1.76 g/cm3, attenuation **0.8 dB / wavelength**
- f=350 Hz, 5 kn, pressure every 1 s, **range increment ≈2.5 m**
- r1=5010 m（生成用）；span 原文 490 m → **勘误 4990 m**
- KRAKEN 生成压力场

## water-bottom interface

| 项 | 结论 |
| --- | --- |
| 正文 | “remains at that value until a depth of **88 m**. The bottom has…” |
| Fig.1 SSP 插图 | 深度轴 0–100 m；剖面在约 88 m 处终止并接底质描述；与正文一致 |
| 判定 | **H = 88 m**，状态 `TEXT_SPECIFIED` + `FIGURE_CONSISTENT` |
| 是否为反调 16-modes | **否**。不因 mode 数调 H。 |

若需极小不确定集：正文/图均指向 88 m，**不**预注册 89/90。图轴 100 m 只是绘图框，不是水深。

## KRAKEN bottom attenuation 换算

论文：0.8 dB / wavelength。
主 baseline 取底质波长：lambda_b = cp/f = 1650/350 m
ap = 0.8 / lambda_b = 0.8 * f / cp ≈ 0.169697 dB/m
（若用水层波长则 ≈0.8*350/1478≈0.189 dB/m；单因素以后再扫，不改主 baseline。）

## paper range coordinate

r = 5010 : 2.5 : 10000 m  （span 4990 m）
**不用** 5 kn 精确 2.5722 m/s 替换论文 2.5 m。
