# R3-C2.2C-FIX2 报告

UTC: 2026-09-24T04:30:40.403724+00:00

## 根因（已证实）

旧 `.flp` 将 `1 350.0` 写在官方 **M（模态数）** 位置 → FIELD 只叠加 **M=1**。
**不是** Eq.(3) 多模失败，**不是** Yang 失败。

官方格式（field.hlp）：`TITLE / 'RA' / M / NPROF RPROF / NR R / NSD SD / NRD RD / NRR RR`。
频率来自 `.mod`，不写入 `.flp`。NRR 必须 = NRD。

## M 差分

- M=1 vs M=23 相对差 = 1.409e+00（须 > 1e-3）
- M=9999 == M=23：True
- flp_parse_ok = True

## FIELD vs Eq3

FIELD_EQ3_MULTIMODE_VALIDATED

见 `field_vs_eq3_by_modecount.csv`、`field_incremental_mode_check.csv`（c_ref 取自 M=1 固定）。

## 判定

### `C2_2C_FIX2_PAPER_REPRO_PARTIAL`

FLP/FIELD OK (FIELD_EQ3_MULTIMODE_VALIDATED); 3/4 fidelity; shallow/deep may still separate

YANG_ROUTE_UNDECIDED。不进 E-STD / 50–60 km / 2.4 km CZ / δ 搜索 / MC / AR / P5。

## 停止

## FIELD vs Eq3
| M | corr | residual | fixed c_ref residual |
| --- | --- | --- | --- |
| 1 | 1.0000 | 4.9e-4 | 4.9e-4 |
| 23 | 1.0000 | 4.6e-4 | 4.6e-4 |
| 9999 | =M23 | same | same |

M=1 vs M=23 rel diff 1.41; 9999==23 hash. Incremental m=2 corr~1.0. field_status=FIELD_EQ3_MULTIMODE_VALIDATED.

## FIELD paper cases
| zs | zr | z_hat | err | ratio | ok |
| --- | --- | --- | --- | --- | --- |
| 4 | 18 | 5 | 1 | 0.88 | Y |
| 4 | 70 | 5 | 1 | 0.87 | Y |
| 50 | 18 | 43 | 7 | 0.86 | N |
| 50 | 70 | 50 | 0 | 1.00 | Y |

Shallow/deep separable; prominence 1/3/5/10% stable. 3/4 engineering gate => PAPER_REPRO_PARTIAL.
Not a Yang accuracy claim. YANG_ROUTE_UNDECIDED.

## 措辞修正 (R3-C2.2D)
incremental check 仅验证 **m=2** 且 PASS。不得写 modes 2-8 全过。聚合 M=1/2/4/8/16/23 与 fixed c_ref 已足够支持 FIELD_EQ3_MULTIMODE_VALIDATED。
