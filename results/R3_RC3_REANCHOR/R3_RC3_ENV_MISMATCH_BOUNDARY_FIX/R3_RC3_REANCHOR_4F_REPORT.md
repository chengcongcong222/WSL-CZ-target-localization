# R3-RC3-REANCHOR-4F SSP integrity fix

UTC: 2026-09-25T15:15:18.658656+00:00

## 判定

### `PROFILED_RC3_MARGIN_SURVIVES_TESTED_SSP_STRESS`

SSP parser fixed (251 pts); E1 frac(ΔJ>0)=1.000 med=2.569; E2 frac=1.000 med=2.527; truth=E_k template=E0

## 旧 REANCHOR-4

`OLD_REANCHOR4_INVALIDATED_BY_SSP_PARSE_BUG`；意外控制=`UNINTENDED_SSP_PERTURBATION_CONTROL`。

## 修复

- SSP parser：`'R'` 前六列声学行；断言 251 / 0–5000 / 严格递增
- E1/E2 定义不变；anchor 校验 max_err=0.00e+00
- write_env 只改 SSP block；roundtrip max_err=5.00e-05
- KRAKEN 全重跑；forward gate **24** 个 perturbed case

## 环境失配（正确 E1/E2）

| env | fraction_tested(ΔJ>0) | median ΔJ |
| --- | ---: | ---: |
| E1 | 1.000 | 2.569 dB |
| E2 | 1.000 | 2.527 dB |

## 证据链

深度消元通过 → IID TL 通过 → **环境失配：本报告** → S1 仍暂缓。

禁止：S1/S0、随机 TL、Liang、P5、改 E1/E2 定义。
