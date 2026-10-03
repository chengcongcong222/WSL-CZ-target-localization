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

Only the A1 search-tractability PRE-RUN FREEZE checkpoint is currently authorized. A2/B/SSP/P5 are not opened. R4 remains at 0%.

## A1 scientific questions

Validate genuinely off-grid horizontal observation generation without truth access by the estimator; separate coarse-grid top-1, full survivor envelopes and the quantization floor; measure bearing-sigma dependence with multiple realizations; report either an empirically stable tested region or NO_STABLE_REGION_ESTABLISHED.

The inherited forward/scoring structure and the pre-run design are recorded in `../R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_DESIGN.md`. R3 scripts and numerical evidence must remain unchanged. Results are synthetic matched E0 controls with perfect relative-TL observations, not end-to-end observed-data or ocean performance guarantees.

## Accepted structural stop and repair Gate

The research lead independently accepted 4633ff0: A1-1 off-grid pipeline passed; A1 remains blocked and receives no progress credit. R3 remains frozen. A1-FIX develops observation-derived continuous bearing hypotheses and multifrequency search. Even a repair PASS would leave R4 at 0% until the full frozen A1 statistical experiment is completed and audited.

V1 had a degenerate theta-offset global-search initialization. Its partial evidence is retained under INITIAL_IMPLEMENTATION and is excluded from V2 confirmation. V2 corrects that programming defect without changing search budgets/thresholds, then freezes a fresh holdout seed 2026100203 before observations. V2 noisy regression does not recover every matched basin, so A1-FIX remains blocked pending search coverage/convergence work. Do not open A2, B1 or P5.

## Accepted FIX1 and FIX2 acoustic coverage Gate

The research lead accepted the independent audit of `956f2dfba8e719561641fd135f246a650fe76dca`: continuous forward self-match, truth-independent continuous RC2 and legacy R3 preservation pass. Acoustic global coverage remains blocked. The original nine regression truths and six Q holdout truths are all development material for every subsequent method change.

FIX2 measures likelihood-feasible oracle basin sections and coupled sensitivity without supplying oracle coordinates to search. Three preregistered deterministic range/radial meshes, separated depth-branch beams, continuous local refinement and direct-modal final verification test coverage and convergence. Both raw component results and cumulative candidate retention are exported; retention alone does not certify convergence. An independent Sobol/Nelder-Mead family checks five preselected difficult cases and two further development failures explicitly selected before fresh confirmation.

After method/config/code/budget/threshold freeze, a new seed-defined eight-truth interior off-grid panel is frozen before observation generation, with noiseless and two nominal realizations each. Evidence is isolated under `../R4_A1_FIX2_ACOUSTIC_COVERAGE/`. FIX2 has no management credit; full A1 remains blocked unless all six repair gates pass and the research lead accepts the resulting commit. R3 and prior A1/FIX1 numerical artifacts remain frozen. No A2, depth development, SSP or full bearing-sigma sweep is opened.

## Accepted FIX2 stop and tractability pre-run checkpoint

The research lead accepted `7ab24845e6e1551b75287fefb1ab662e92b395b8` as A1_FIX2_BLOCKED_BY_UNCLOSED_ACOUSTIC_COVERAGE: 2/6 Gates pass, raw coverage remains unclosed. R3/A1/FIX1/FIX2 numerical evidence is frozen. A likelihood-section width is not an optimizer capture basin; finite catalogs without an alias do not establish global uniqueness; search failure is not physical non-identifiability.

The new controlling staged instruction authorizes only external request-budget enforcement, SHGO/DIRECT API/exact-objective/exclusion tests, cost calibration and an independent pre-run freeze commit. Evidence is isolated under `../R4_A1_SEARCH_TRACTABILITY_AUDIT/`. Both solver families have frozen per-case/per-depth admitted-request caps 16/64/256. No noisy development case or fresh panel has run. The checkpoint awaits the research lead's audit and separate release; all scientific tractability Gates remain NOT_EVALUATED. No new scientific progress credit is awarded. A2, depth/B, SSP and P5 remain unopened.
