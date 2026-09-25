# R3-C2.4B-3A Hankel amplitude bridge

UTC: 2026-09-25T03:59:05.755902+00:00

## 判定

### `HANKEL_AMPLITUDE_BRIDGE_VALIDATED`

rank_corr=0.643, dominant_overlap=0.50, mean M0=6.0

## 范围

- **机制验证**，非 Liang Fig.5/6 严格复现
- R0 无噪；p=7；span 4990/1990；z_s=4/50
- AR 只给 k̂；幅度 **Eq.(7) Hankel**
- 两种 shading 控制；不调 window

## 产物

`HANKEL_AT_AR_PEAKS.csv`（复数 g）· `AMPLITUDE_SANITY_CHECK.csv` · `MODE_COUNT_FOR_DZ.csv` · `DZ_INPUT_CONTRACT.md`

## 停止

不计算 b/D(z)/Fig.5–6；不进 E-STD/MC/P5。
