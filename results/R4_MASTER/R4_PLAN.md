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

The frozen 21-case noisy development at `9cfa99c14bb80d08e6f16210ebf89fd9c541b850` is accepted by the research lead's third independent audit. The current document-only stage is `R4_A1_CLOSEOUT_AND_ROUTE_DECISION`; the A1 search-repair chain is CLOSED and R4_ROUTE_REDESIGN_REQUIRED. No further experiment is authorized. A2/B/SSP/P5 are not opened. R4 remains at 0%.

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

## Accepted infrastructure and unreleased execution harness

The lead independently accepted b4839b295777127ec0c8ade56b76c776db99148b as PRE_RUN_INFRASTRUCTURE_FREEZE_ACCEPTED. Its core objective/controller and frozen 16/64/256 per-depth budgets remain unchanged. The new authorization is only DEVELOPMENT_EXECUTION_HARNESS_FREEZE: non-oracle 21-case manifest, independent raw-branch driver, threshold-hit/witness-cluster semantics and Gate logic, tests and immutable release manifest. It does not authorize real noisy development.

At the harness-freeze checkpoint, development_released=false and noisy_development_runs=0; the second audit and separate development release were pending. Raw top witnesses and discrete subthreshold clusters are not local minima or basin certificates. No management credit is awarded; R4=0%; A2/depth/SSP/P5 unopened.

## Accepted harness and completed frozen development

The lead accepted `c5c9ba2e1258e73486e62df26499e03cce204f27` and authorized only the 21-case noisy development, recorded before search at `d731cf51d72b41dea7fd4c7a48fcf21a623f7bc4`. The unchanged runtime-released harness executed all 126 raw runs / 2646 depth branches. Its decision is `A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED`; raw Gate = False; dual Gate = NOT_REACHED. All failures and raw request/witness evidence are preserved under ../R4_A1_SEARCH_TRACTABILITY_AUDIT/. See DEVELOPMENT_RUN_REPORT.md and DEVELOPMENT_FINAL_DECISION.json. At that development-submission checkpoint the third independent audit was pending (now accepted; see the closeout below); no T4, budget/tolerance change, solver substitution, FIX3/FIX4 or fresh confirmation is authorized. R4 remains 0%; A2/depth/SSP/P5 remain unopened.

## Accepted third audit: A1 closeout and route decision

Accepted development SHA: `9cfa99c14bb80d08e6f16210ebf89fd9c541b850`. The research lead accepted `A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED`. Both SHGO and DIRECT recovered 0/21 at every frozen raw budget; both raw T2/T3 convergence counts were 0/21, so dual agreement was NOT_REACHED. All 126 raw runs were execution-valid, 2646 depth branches were accounted for, and hard-budget enforcement failures were 0. The existing cold reconstruction passed 14479 checks / 0 failures across 232973 unique case/depth/state combinations (maximum J discrepancy 1.1619603057511085e-9 dB). The accepted negative result does not establish physical non-identifiability, absence of a correct match, global uniqueness or universal algorithm impossibility.

FIX2 remains frozen BLOCKED at `7ab24845e6e1551b75287fefb1ab662e92b395b8` (2/6 Gates). Its 18/21 cumulative recovery is not raw convergence: B3 raw is 17/21 and both B2/B3 raw recovered is 16/21. Independent recovery/agreement is 5/7 and 3/7; FIX2 fresh was executed and recovered 8/8 noiseless and 13/16 nominal. Tractability fresh was NOT_REACHED and never executed. A likelihood-section width is neither an optimizer capture basin nor ranging precision.

The current engineering decision is `A1_SEARCH_REPAIR_CHAIN_CLOSED` and `R4_ROUTE_REDESIGN_REQUIRED`, scoped to the tested frozen direct-TL global-search proposal. No T4, increased budget, FIX3/FIX4, new optimizer or fresh confirmation is authorized. R3 and all A1/FIX/FIX2/tractability numerical artifacts remain frozen. R4=0%; A2/depth/B/SSP/P5 remain unopened.

Route priority: **R1 > R2 > R3**. Recommend **RC2 temporal candidate narrowing + RC3 local propagation anchor**; second, redesign a searchable RC3 observable; third, retain the current exact-TL matcher only as oracle/local diagnostic/conditional rescoring. R1 is a recommendation, not a validated estimator. Frozen R3-1B straight-line time accumulation and R3-1C tested turns did not narrow the RC2-only 15 km range domain. A future design must identify an observation-derived narrowing mechanism, an honest cold start/prior source, retained candidate support, and conditional capture/convergence; recursion alone is not presumed to solve the range ridge.

The next stage is a separately defined and authorized **route redesign Gate**, not A2 or another A1 search repair. Its quantitative criteria, observation/initialization scope, budgets and confirmation rules remain NOT_EVALUATED; this closeout does not authorize implementation, experiments or a new panel. Stop after committing and pushing the document-only closeout.

Current entry points: [evidence chain](../R4_A1_CLOSEOUT_AND_ROUTE_DECISION/A1_EVIDENCE_CHAIN.md), [closeout report](../R4_A1_CLOSEOUT_AND_ROUTE_DECISION/A1_CLOSEOUT_REPORT.md), [route decision](../R4_A1_CLOSEOUT_AND_ROUTE_DECISION/R4_ROUTE_DECISION.md), [machine decision](../R4_A1_CLOSEOUT_AND_ROUTE_DECISION/R4_A1_CLOSEOUT_DECISION.json), [GPT sync](../R4_A1_CLOSEOUT_AND_ROUTE_DECISION/GPT_SYNC.md). The ledger is append-only: the new acceptance event supersedes the historical pending-third-audit row without deleting it. R4_PROGRESS.json stays byte-identical as the development snapshot; read the new closeout decision for current route status. Previous freeze manifests remain untouched; their old plan/ledger text identities are checked against the accepted parent Git blobs.
