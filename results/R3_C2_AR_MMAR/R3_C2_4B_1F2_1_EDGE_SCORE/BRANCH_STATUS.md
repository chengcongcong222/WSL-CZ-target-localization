# BRANCH_STATUS

UTC: 2026-09-24T11:51:27.216561+00:00

| 分支 | 角色 | single F1 | multi F1 | multi complete | multi recall |
| --- | --- | ---: | ---: | ---: | ---: |
| p=7 `PRINTED_LITERAL_CONTROL` | **PAPER_PRINTED_PRIMARY_BRANCH** | 1.000 | 0.972 | 0.92 | 0.958 |
| p=133 `OUR_MOVING_SAMPLE_INTERPRETATION` | **OUR_INTERPRETATION_SENSITIVITY_BRANCH** | 1.000 | 0.972 | 0.92 | 0.958 |

## Provenance（进入 4B-2 前按来源冻结，非按结果择优）

- p=7 = 论文打印主分支 → **4B-2 论文复现以 p=7 为主**
- p=133 = 我方解释敏感性分支 → 继续报告，**不**决定 Liang paper reproduction 是否通过
- p=133 小扰动下系数/峰位更敏感 → `P133_CONDITIONING_LIMIT`（见 1F2 diagnostic）

## 冻结（未改）

`P1_COMPLEX_EXP_IDENTITY_VALIDATED`
`OUR_EXACT_SOLVER_FOR_EQ24_VALIDATED`
`AR_CORE_IMPLEMENTATION_VALIDATED`（P1 + 复数 MCOV 共轭）
`P7_PRINTED_LITERAL_BRANCH_ACCEPTED`
`P133_MOVING_SAMPLE_BRANCH_CONDITIONING_LIMITED`
