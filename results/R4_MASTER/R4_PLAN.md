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

## Accepted structural stop and repair Gate

The research lead independently accepted 4633ff0: A1-1 off-grid pipeline passed; A1 remains blocked and receives no progress credit. R3 remains frozen. A1-FIX develops observation-derived continuous bearing hypotheses and multifrequency search. Even a repair PASS would leave R4 at 0% until the full frozen A1 statistical experiment is completed and audited.

V1 had a degenerate theta-offset global-search initialization. Its partial evidence is retained under INITIAL_IMPLEMENTATION and is excluded from V2 confirmation. V2 corrects that programming defect without changing search budgets/thresholds, then freezes a fresh holdout seed 2026100203 before observations. V2 noisy regression does not recover every matched basin, so A1-FIX remains blocked pending search coverage/convergence work. Do not open A2, B1 or P5.
