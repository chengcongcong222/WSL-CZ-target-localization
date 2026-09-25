# ENV_MISMATCH_AXIS_LOCK

UTC: 2026-09-25T12:59:04.576212+00:00

`PROJECT_PREDEFINED_SSP_MISMATCH_STRESS_AXIS`（工程敏感性轴，非海况概率）

| 档 | 定义 |
| --- | --- |
| E0 | nominal E-STD（冻结 zgrid SSP） |
| E1 | `SSP_VERTICAL_SHIFT_50M`：c(z)→c(z−50 m)，端点外推 |
| E2 | `SSP_VERTICAL_SHIFT_150M_PLUS_2MPS`：c(z)→c(z−150 m)+2 m/s |

信号固定 `S2_LIKE_TRACKABLE_FOUR_TONE_UPPER_BOUND`；ε=0（无随机 TL 误差）。

**真值用 E_k；估计器模板只用 E0**（否则不是失配）。
