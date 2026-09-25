# R3-RC3-REANCHOR-4G mode-depth tabulation fix

UTC: 2026-09-25T16:02:28.471191+00:00

## 判定

### `PROFILED_RC3_MARGIN_SURVIVES_TESTED_SSP_STRESS`

tail+tabulation OK; z_true signatures distinct; E1 frac(ΔJ>0)=1.000 med=2.569; E2 frac=1.000 med=2.575

## 修复

- `'R'` 后 tail **逐行复制**（51 receivers 150:2:250）
- `.mod` 含 180/200/220 且 φ 向量不全同
- FIELD 无 `Modes not tabulated` 警告
- 源深签名 RMS(L180−L200)、RMS(L220−L200) 非机器零
- B01/C01/A01 的 y(180/200/220) 非重复

旧 4F：`REANCHOR4F_BLOCKED_BY_MODE_DEPTH_TABULATION`；z=200 单点=`Z200_ENV_STRESS_POSITIVE_CONTROL`。

## 结果

| env | fraction_tested(ΔJ>0) | median ΔJ |
| --- | ---: | ---: |
| E1 | 1.000 | 2.569 |
| E2 | 1.000 | 2.575 |

S1 仍暂缓。
