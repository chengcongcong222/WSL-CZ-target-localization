# R4 master plan

Baseline: `166575043d9c884f6a9ecf62333b79fbcc7e1a46`. R3 is frozen and excluded from R4 management progress.

The research lead owns the route, Gate acceptance, evidence wording and final independent audit. Codex executes the authorized stage autonomously. Writing a script is not scientific completion. No progress is credited before the research lead audits the resulting commit.

| Module | Objective | Management weight |
|---|---|---:|
| A1 | Off-grid generalization and bearing-accuracy baseline | 15% |
| A2 | Platform, sensor and line-spectrum single-factor boundaries | 15% |
| A3 | Non-extreme scenes and joint-perturbation Monte Carlo | 15% |
| A4 | Engineering working-envelope freeze | 5% |
| B1 | Conditional depth identifiability | 10% |
| B2 | Off-grid continuous-depth estimation | 10% |
| B3 | SSP-depth coupling and environmental adaptation | 10% |
| B4 | Depth working-envelope freeze | 5% |
| C | Application metric freeze | 15% |

Only A1 is currently authorized. A2/B/P5 are not opened. R4 starts at 0%.

## A1 scientific questions

Validate genuinely off-grid horizontal observation generation without truth access by the estimator; separate coarse-grid top-1, full survivor envelopes and the quantization floor; measure bearing-sigma dependence with multiple realizations; report either an empirically stable tested region or NO_STABLE_REGION_ESTABLISHED.

The inherited forward/scoring structure and the pre-run design are recorded in `../R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_DESIGN.md`. R3 scripts and numerical evidence must remain unchanged. Results are synthetic matched E0 controls with perfect relative-TL observations, not end-to-end observed-data or ocean performance guarantees.
