# Gate2B independent audit acceptance

AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED_ACCEPTED

AUX1_HALF_T5_JOINT_BUDGET_ESTABLISHED accepted.

Parent execution SHA: 333ead6461a7574a469495033d48fd89682aa496
Design SHA: 8c8c0e52df09f4cf1111201396d6a2dec3c5f8b0

Accepted historical execution:2904cells;30000trials/cell;4tests PASS;140702independent checks PASS,0FAIL. These are prior-stage results; this closeout performs no new experiment.

| Family / anchor | Maximum tested passing lambda | Worst P95 (%) | Next registered fail lambda | Fail P95 (%) |
|---|---:|---:|---:|---:|
| T5 / A | 0.5 | 4.0481 | 0.75 | 5.1309 |
| T5 / B | 0.75 | 4.5794 | 1 | 5.4940 |
| T10 / A | 0.75 | 8.3565 | 1 | 11.1964 |
| T10 / B | 0.75 | 8.9366 | 1 | 12.0425 |

T5 A maximum tested joint lambda=.50; T5 B=.75. Independent T10 family A/B=.75. Do not hide stronger tested B evidence behind the half point; keep B lambda1FAIL5.4940%. T10 lambda1both FAIL. Gate2A single-factor maxima are not interchangeable with Gate2B simultaneous vectors.

Gate0: conditional auxiliary geometry; no actual measurement chain. Gate1: ideal baseline/bearing tradeoff. Gate2A: signed single-factor requirements only. Gate2B: two frozen joint static families, discrete ranges and signed corners. No continuous error-box, off-grid or dynamic certificate; physical identifiability and engineered cold start are distinct.

Single-array cold start NOT_ESTABLISHED; depth engineering CLOSED; exact horizontal depth mechanism SUPPORTED_ORACLE_ONLY. Actual auxiliary hardware UNKNOWN; sync/association NOT_NUMERICALLY_VALIDATED; R4-A1-NEW NOT_OPENED; R4=0%.

Client package READY; no customer message sent. Numerical architecture work STOPPED. No automatic geometry/MC/optimizer/acoustic/depth or timing/association/tracking run. Next step is research-lead client confirmation, then reviewed architecture decision.

## Saved evidence references

- [Gate0 conditional geometry](../R4_ROUTE_RESET_AUXILIARY_NODE_GATE0/AUX_GATE0_REPORT.md)
- [Gate1 baseline/bearing tradeoff](../R4_AUX_GATE1_CROSSTRACK_REQUIREMENT/AUX_GATE1_REPORT.md)
- [Gate2A single-factor requirements](../R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET/AUX_GATE2A_REPORT.md)
- [Gate2B joint static requirements](../R4_AUX_GATE2B_JOINT_ERROR_BUDGET/AUX_GATE2B_REPORT.md)

Every joint requirement binds both the execution frontier row and the design vector definition in AUX_REQUIREMENT_PROVENANCE.csv. No new trial draw or numerical reconstruction was generated in this package.
