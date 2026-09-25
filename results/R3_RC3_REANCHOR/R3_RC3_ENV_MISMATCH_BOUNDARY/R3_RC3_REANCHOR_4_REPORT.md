# R3-RC3-REANCHOR-4 env-mismatch boundary

UTC: 2026-09-25T12:59:04.576212+00:00

## 判定

### `PROFILED_RC3_MARGIN_SURVIVES_TESTED_SSP_STRESS`

E1 fraction(ΔJ>0)=1.000 median=2.578 dB; E2 fraction=1.000 median=2.551 dB; truth=E_k, template=E0; S2 four-tone; no random TL error

## 轴

E0 nominal · E1 SSP 垂向平移 50 m · E2 平移 150 m + 2 m/s  
`PROJECT_PREDEFINED_SSP_MISMATCH_STRESS_AXIS`（工程敏感性，非海况概率）

**真值 E_k，模板恒为 E0**（ε=0，无随机 TL 误差；S2 四稳定线谱）

## 结果

| env | fraction_tested(ΔJ>0) | median ΔJ |
| --- | ---: | ---: |
| E1 | 1.000 | 2.578 dB |
| E2 | 1.000 | 2.551 dB |

A05 继续 `RELATIVE_TL_STATIC_RANGE_BLIND_CASE`（不进主统计）。

## 口径（相对 REANCHOR-3）

- `fraction_tested` 非概率（5 seeds 非独立 MC）
- 仅 IID 零均值 TL 扰动已测；慢变/环境误差另论（本阶段即环境）
- 2 dB 最差 ΔJ≈0.11 dB（A10）仍为正，但不宽裕

## 禁止

S0/S1、随机 TL、Liang、P5、重挑样本、改 E1/E2。
