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

The frozen 21-case noisy development at `9cfa99c14bb80d08e6f16210ebf89fd9c541b850` is accepted by the research lead's third independent audit. The closeout at `206f6c2f46bc3c52890a32004cf54255ed30fcef` is independently accepted. The current document-only stage is `R4_ROUTE_REDESIGN_GATE_0_COLD_START_INFORMATION_AUDIT`; R1 cold start is NOT_ESTABLISHED and tracking is only conditionally plausible. The A1 search-repair chain stays CLOSED and R4_ROUTE_REDESIGN_REQUIRED. No further experiment is authorized. A2/B/SSP/P5 are not opened. R4 remains at 0%.

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

## Accepted closeout and Gate-0 cold-start information audit

The lead accepted `206f6c2f46bc3c52890a32004cf54255ed30fcef` as `R4_A1_CLOSEOUT_AND_ROUTE_DECISION_ACCEPTED`. Its R1 > R2 > R3 ranking was a pre-Gate0 architecture recommendation, not evidence of acquisition feasibility. Gate-0 now distinguishes acquisition from tracking and selects user-defined outcome B: **R1_COLD_START_NOT_ESTABLISHED** and **R1_TRACKING_ONLY_CONDITIONALLY_PLAUSIBLE**. No qualifying cold-start range contraction source was identified in the inherited noisy nominal scope. The audited bootstrap using a previous posterior whose origin is the failed blind direct-TL matcher has `COLD_START_CIRCULAR_DEPENDENCY`; independent legitimate acquisition followed by tracking is a different conditional and acyclic interface, not an available current capability.

RC1's frozen role is a 45--60 km conditioned support; its contraction ratio is NOT_IDENTIFIABLE_FROM_CURRENT_EXPERIMENT, and no narrower first-CZ prior is admitted. R3-1B straight accumulation and all five R3-1C tested turns retain 15 km RC2-only range width. FIX1's 21 noisy continuous-RC2 clouds retain truth in the likelihood but span 14.9008764373--14.9987589195 km (median 14.9734269641 km). These are finite exported samples, not confidence certificates. The 15 noiseless rank-four singleton controls and conditional turn+RC3 coarse positives remain valid exceptions; Gate-0 does not claim bearing has no range information or universal cold-start impossibility.

G0-1/2/3 are not closed in the audited acquisition scope. G0-4 (support retention plus domain contraction metrics) and G0-5 (conditional/local RC3 role) are satisfied only as structural definitions, with no quantitative thresholds or performance PASS. Prediction preserves/propagates prior information and can widen uncertainty; it does not create the first credible range posterior. R1 acquisition is not admitted for quantitative design or implementation. Tracking remains an unvalidated conditional proposal with engineering scope **POST_ACQUISITION_TRACK_REFINEMENT**, requiring independently obtained credible local joint support.

Recommend the research lead's next instruction address **R2 acquisition observable/information-source document design**; retain R1 as a conditional post-acquisition module and Route R3 as oracle/local diagnostic. This refines the old ranking by separating acquisition and tracking; it does not authorize R2 features, experiments, new observation assumptions or any R1 implementation. Cold-start information-source identification is false; the circular-dependency flag is scoped to the unsupported bootstrap, not an existing code defect.

See [source matrix](../R4_ROUTE_REDESIGN_GATE_0/COLD_START_INFORMATION_SOURCE_MATRIX.csv), [acquisition versus tracking](../R4_ROUTE_REDESIGN_GATE_0/R1_ACQUISITION_VS_TRACKING.md), [dependency audit](../R4_ROUTE_REDESIGN_GATE_0/R1_CIRCULAR_DEPENDENCY_AUDIT.md), [Gate-0 report](../R4_ROUTE_REDESIGN_GATE_0/ROUTE_REDESIGN_GATE0_REPORT.md), [decision](../R4_ROUTE_REDESIGN_GATE_0/R1_GATE0_DECISION.json), [GPT sync](../R4_ROUTE_REDESIGN_GATE_0/GPT_SYNC.md). Prior closeout and all historical numerical/code/budget/Gate/manifest evidence remain unchanged; only this plan and append-only ledger change. R4_PROGRESS.json remains a historical development snapshot. Management weights and R4=0% stay unchanged; A2/depth/SSP/P5 remain unopened. This audit creates no numerical experiment, no new threshold or fresh panel, and stops after commit/push pending lead review and a separate next-stage instruction.


## Accepted Gate-0 and Gate-1 R2 observable audit

The research lead accepted Gate-0 `193abb7cd13a38fac91ad1715e6873dd3d6440e9` as **R4_ROUTE_REDESIGN_GATE0_ACCEPTED**. The current documentation stage is **R4_ROUTE_REDESIGN_GATE_1_R2_ACQUISITION_OBSERVABLE_DESIGN_AUDIT**. R1 cold start remains NOT_ESTABLISHED; tracking remains CONDITIONALLY_PLAUSIBLE_POST_ACQUISITION_ONLY. RC3 remains LOCAL_CONDITIONAL_ANCHOR_ONLY. Earlier route rankings are historical recommendations; the latest decisions distinguish acquisition from conditional tracking.

Gate-1 audits all six required families and chooses **CZ_ENVELOPE** first and **INTERFERENCE_FRINGE** second for measurement/information-chain clarification. Both are CONDITIONAL, not performance winners. Parent-candidate classes: A=0, B=5, C=1; candidate_admitted_for_quantitative_design=false. Ordered/differential frequency, multipath TDOA and Doppler are conditional; the current additional HLA spatial range-anchor variant is REJECT_CURRENT_SCOPE. Specific rejected variants do not reject every possible member of a candidate family. No candidate simultaneously closes G1-1 through G1-6. Structural admission does not require completed numerical performance, but cannot rely on unspecified source/signal/receive conditions or oracle information.

Current A1 observations are bearing plus per-frequency/per-window demeaned relative levels at 201/235/283 Hz, not archived raw coherent HLA signals, dense spectra, delay tracks or independent Doppler tracks. Source-level differencing alone does not remove unknown cross-frequency source spectra. The valid Xu2024 source PDF was read from the D-drive original project papers and matches the prior primary-lock SHA; its broadband/VLA depth-slope mechanism does not directly transfer to current HLA. The repository's untracked non-PDF same-name file remains untouched. Mobile-HLA review, bottom-bounce, Yang and Liang papers are audited by required observable and scenario.

Labeled MMAC closure, single-depth aggregate-autocorrelation closure, all-eigenray physical upper bounds and constant-beta CZ transfer limits all remain frozen. Aggregate delay mechanisms do not automatically require long-term path labels, but still require a signal/extractability/stability chain; generalized intensity fringe remains conditional rather than declared useless or available. No absolute source or propagation phase is taken from the forward model.

Recommend a separately instructed **R2_MEASUREMENT_CHAIN_AND_SIGNAL_CONDITION_SPECIFICATION** document stage before any quantitative Gate design. This recommendation authorizes no new signal assumptions, simulation, estimator, optimizer, panel, aperture, VLA, Z maneuver or implementation. No numerical thresholds are set. Gate-1 review is pending; no progress credit is awarded. See [matrix](../R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/R2_OBSERVABLE_ADMISSIBILITY_MATRIX.csv), [measurement chain](../R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/R2_MEASUREMENT_CHAIN_AUDIT.md), [literature](../R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/R2_LITERATURE_TRANSFERABILITY.md), [ranking](../R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/R2_CANDIDATE_RANKING.md), [report](../R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/ROUTE_REDESIGN_GATE1_REPORT.md), [decision](../R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/R2_GATE1_DECISION.json), [GPT sync](../R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/GPT_SYNC.md).

Only this plan and the ledger receive append-only updates. All historical numerical/code/budget/Gate/manifest evidence, closeout and Gate-0 remain unchanged. R4_PROGRESS.json stays the historical development snapshot; management weights and R4=0% remain unchanged. New numerical experiments=0; A2/depth/SSP/P5 remain UNOPENED. Commit, push, verify HEAD==remote/main, then stop pending lead review and a new explicit stage instruction.


## Accepted Gate-1 and Gate-2 measurement / signal specifications

The research lead accepted Gate-1 `9e50b0d949523c85fb35c1168ff77c37ddf37aeb` as R4_ROUTE_REDESIGN_GATE1_ACCEPTED. Current documentation stage: R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CHAIN_AND_SIGNAL_CONDITION_SPECIFICATION. M0 is the frozen synthetic bearing plus 201/235/283 Hz window-centered relative-level interface. M1 actual receive products, ADC sampling, synchronization, calibration, equipment geometry and real target spectrum remain UNKNOWN / NOT_FROZEN; neither model configurations nor proposed project routes establish device capability. M2 source-oracle/cooperative shortcuts are not admitted.

Two minimum passive proposals are specified: CZ_ENVELOPE uses associated multi-line time-intensity/relative-level data and auditable normalization/gain/noise metadata; INTERFERENCE_FRINGE needs a dense adequately sampled target-bearing intensity spectrogram with genuine spectral occupancy and source-trend separation. Neither requires known waveform, absolute source phase, known source spectrum, known emitted f0 or transmitter cooperation. Raw element waveforms are optional for R2 consumers if an auditable beam spectral product is supplied; upstream relative-channel synchronization/beam calibration must still be supported. Dense bandwidth is not mandatory for the selected CZ line-time variant. M0 can express temporal proxies without a new sensor/output; the extension concerns real measured products and their provenance, not an assertion that all representations need new acquisition hardware.

Gate-2 outcome B: R2_REQUIRES_MEASUREMENT_CONTRACT_EXTENSION. Both leading contracts are PROPOSED_CONTRACT; G2-1 is not closed because current product existence and research-lead acceptance of the proposed extension are absent. Other true gate flags express conditional structural nuisance/range/falsifiability reasoning only. candidate_admitted_for_quantitative_design=false. No claim of real UUV spectrum, numerical range contraction or engineering capability is made.

Recommend the separately instructed R2_MEASUREMENT_CONTRACT_EXTENSION_DECISION: research lead decides whether to accept outputs/conditions or seek concrete system evidence. This is not automatic quantitative design, data collection, receiver implementation, procurement or experiment authorization. R1 cold start NOT_ESTABLISHED, post-acquisition tracking conditional, RC3 LOCAL_CONDITIONAL_ANCHOR_ONLY. Historical B1/TDOA/WI and A1 closure are preserved; A2/depth/SSP/P5 remain UNOPENED.

See [receiver audit](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/CURRENT_RECEIVER_CHAIN_AUDIT.md), [measurement matrix](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/R2_MEASUREMENT_CONTRACT_MATRIX.csv), [signal matrix](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/R2_SIGNAL_CONDITION_MATRIX.csv), [CZ contract](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/CZ_ENVELOPE_MEASUREMENT_CONTRACT.md), [fringe contract](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/INTERFERENCE_FRINGE_MEASUREMENT_CONTRACT.md), [decision](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/R2_GATE2_DECISION.json), [report](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/ROUTE_REDESIGN_GATE2_REPORT.md), [GPT sync](../R4_ROUTE_REDESIGN_GATE_2_MEASUREMENT_CONTRACT/GPT_SYNC.md). Only append this plan/ledger. Gate-1 and all historical code/data/thresholds/freezes are unchanged. R4_PROGRESS.json and nine management weights remain unchanged; R4=0%, new numerical experiments=0. Commit, push, verify local HEAD==remote/main, then stop pending independent audit and a new explicit stage instruction.


## Accepted Gate-2 and Gate-3 CZ envelope quantitative design draft

The research lead accepted Gate-2 `dcad14337491ac5bf90acacf0be34a57eafac3f1` as R4_ROUTE_REDESIGN_GATE2_ACCEPTED, and accepted CZ_ENVELOPE Contract A for QUANTITATIVE GATE DESIGN ONLY. INTERFERENCE_FRINGE remains CONDITIONAL_NOT_ADMITTED. This new lead decision supersedes the historical unaccepted-extension management recommendation without changing Gate-2 files. Actual M1 receiver capability remains UNKNOWN / NOT_FROZEN; real UUV multi-line signal is not verified.

Current stage: R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_QUANTITATIVE_DESIGN_DRAFT. Three M0-compatible proposed families are WINDOW_SHAPE_Q2, LINE_BLOCK_MEANS_B3, NORMALIZED_INTENSITY_BLOCKS_B3. Recommend WINDOW_SHAPE_Q2 for static-offset invariance, low-order shape, Contract A compatibility and explicit drift projection, not measured performance. Physical smoothing/block scale selection is not yet justified; q=2/three blocks are reviewable finite drafts. One primary and at most one preregistered sensitivity must be chosen before execution. No additional fringe/TDOA/Doppler/spatial observation enters the main scheme.

N0 static line/window offsets, N1 bounded linear source/gain drift with bound not frozen, and N2 unconstrained nonidentifiable symbolic limit are distinguished. Proposed search retains conservative full 4D support and all range intervals; no truth-motion slice/top1/oracle warm start. Whole-cell bounds and grid/hard budgets remain pre-run requirements, not a promised solved coverage algorithm. Unknown regions are retained and prevent a coverage/convergence claim. Truth is evaluation-only.

Thresholds are PROPOSED_NOT_FROZEN: lenient W_hull<=7.5km/C<=0.5; nominal <=5km/C<=1/3 (recommended); strict <=3km/C<=0.2. All primary cases require joint truth/reference retention, actual range increment and raw grid/coverage closure; components and endpoint stability are proposed as reviewable limits. These are not formal <10% acceptance or downstream RC3 capture guarantees. tau_F/noise, N1 bound and implementation budgets are NOT_FROZEN. Old scientific thresholds are unchanged.

Primary development is the already-used 21 noisy P/Q manifest; 15 noiseless controls are regression only. Old FIX2's 24 holdout cases are also already-used development/regression. Fresh confirmation is NOT_GENERATED and can only follow full representation/nuisance/threshold/budget/search/panel-rule freeze plus new explicit authorization. First-round data are SYNTHETIC_MATCHED_MEASUREMENT_CONTRACT_PROXY, not real receiver validation. Stop on retention loss, no actual range contraction, dense aliases/coverage failure, oracle dependence or nuisance absorption; no unlimited retuning.

See [representations](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/CZ_ENVELOPE_REPRESENTATION_CANDIDATES.md), [nuisance](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/CZ_ENVELOPE_NUISANCE_MODEL.md), [support design](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/CZ_ENVELOPE_SUPPORT_ESTIMATOR_DESIGN.md), [metrics](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/CZ_ENVELOPE_METRICS_AND_BASELINES.md), [thresholds](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/CZ_ENVELOPE_THRESHOLD_PROPOSALS.md), [split](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/CZ_ENVELOPE_DEVELOPMENT_CONFIRMATION_SPLIT.md), [stop](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/CZ_ENVELOPE_STOP_RULES.md), [decision](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/R2_GATE3_DESIGN_DECISION.json), [report](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/ROUTE_REDESIGN_GATE3_DESIGN_REPORT.md), [GPT sync](../R4_ROUTE_REDESIGN_GATE_3_CZ_ENVELOPE_DESIGN/GPT_SYNC.md). Recommended next stage: RESEARCH_LEAD_GATE3_DESIGN_REVIEW_AND_PRE_RUN_FREEZE. Draft review pending, quantitative run NOT_AUTHORIZED, implementation NOT_AUTHORIZED, fresh NOT_GENERATED, new experiments=0. Only append plan/ledger, keep all prior code/data/Gates/freezes and R4_PROGRESS.json unchanged. R4=0%, weights unchanged, A2/depth/SSP/P5 UNOPENED. Commit/push/verify HEAD==remote/main then stop.


## Gate-3A CZ envelope PRE-RUN freeze

Research lead accepted Gate-3 design at e56ba746ebe7370d1dee3d0a181bc348b5fecc5f. Choices: WINDOW_SHAPE_Q2; sensitivity NONE; N0 static per-line/window offsets; NOMINAL; 21 noisy P/Q development only. Policy pushed first at a393c7b45c92c42fd148d9c162faeb5ea35160c7.

Certified-ball whole-cell enclosure, precharged counters, independent raw coarse/fine logic and evaluation sentinel implemented. Fixed non-development calibration only: tau_F=1e-7 dB; coarse widths=(0.25 km,0.625 deg,0.125 m/s,1.875 deg), fine=half; final cell caps=14/29. Forty-seven unit tests pass. Conservative validity does not establish useful bounds or development closure.

[PRE-RUN decision](../R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_PRE_RUN_DECISION.json) and [freeze](../R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_PRE_RUN_DESIGN_FREEZE.json) pending independent audit. Development feature/search=0; entry disabled; fresh NOT_GENERATED. Stop after commit/push. N1, A2/depth/SSP/P5 unopened; R4=0%; prior evidence and weights unchanged.


## Gate-3A PRE-RUN conservative closure correction

Research lead rejected e9316c031d40f11fc07772959d439b965d876cd5: terminal possible cells were incorrectly forced to obtain upper compatibility proof and 14/29 requests were below single-path completion bounds. Old sources/artifacts remain historical and unchanged.

Correction policy pushed first at 98f4e3184752025be73da2ba2bad568bdc0e9992. New engine separates conservative partition completion, uniform compatibility and natural budget completion. Original lower bounds, Q2/N0/tau=1e-7, NOMINAL Gate, domain/grid/profile and primary manifest unchanged. Single-path minima=37/45; corner/midpoint B0/CZ ladders=37/148/439 and 45/180/879.

Outcome: CZ_ENVELOPE_PRE_RUN_SEARCH_CLOSURE_NOT_ESTABLISHED; selected caps={"COARSE": null, "FINE": null}; all-run reference retained=True. Supplementary midpoint-specific geometric minimum=451/579: coarse max439 is below its necessary minimum. Preserve this policy limitation without another tier; no CZ information/physics-impossibility conclusion.

[Corrected decision](../R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_CORRECTED_PRE_RUN_DECISION.json), [fixture summary](../R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_CLOSURE_FIXTURE_SUMMARY.json), [budget limitation](../R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN/CZ_CLOSURE_BUDGET_LIMITATION.md). Sixty-nine tests pass; independent full partition/counter/metadata verification accompanies checkpoint. Development feature/search=0; entry disabled; fresh NOT_GENERATED. Stop after push. R4=0%; N1/A2/depth/SSP/P5 unopened.


## Accepted CZ route closure and authorized R4 B1 design freeze

Research lead accepts parent 20846a9e1ebb10518762587635856c71763da790 as CZ_ENVELOPE_PRE_RUN_SEARCH_CLOSURE_NOT_ESTABLISHED_ACCEPTED. CZ_ENVELOPE_CERTIFIED_GLOBAL_SEARCH_ROUTE_CLOSED: WINDOW_SHAPE_Q2 + N0 + whole-cell certified conservative outer partition + 45-60 km cold-start 4D global support search within preregistered budgets. No C4/F4, increased budgets, grids/tau/enclosure FIX or 21-case CZ development. Mathematical lower-bound mechanism remains historical evidence; no physical impossibility or no-information conclusion. R4-A1 BLOCKED / NOT_COMPLETED, 0/15%, cold-start acquisition implementation not established. R4 remains0%.

New explicit stage R4_B1_CONDITIONAL_DEPTH_PROFILE_IDENTIFIABILITY authorized by lead. Use the complete six historical S7/S8 configurations: triple x3 depths plus four-line x3 depths, not six triple cases; no selection by depth success. Matched E0, frozen historical RMS score, 150:5:250 m, H0 exact horizontal condition plus eight mechanically inherited one-axis grid steps. [Design freeze](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_DESIGN_FREEZE.json), [case provenance](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_CASE_PROVENANCE.csv), [CZ route closeout](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_CZ_ROUTE_CLOSEOUT.md). No new B1 numerical depth evaluation before design commit/push and remote verification. Then execute once unchanged, commit/push and stop for independent audit. B2 and A2/A3/A4/B3/B4/C/P5 NOT_OPENED. No progress credited; R4_PROGRESS.json unchanged.


## R4 B1 conditional-depth execution checkpoint

B1_DEPTH_IDENTIFIABLE_ONLY_UNDER_TIGHT_HORIZONTAL_CONDITIONING; PENDING_RESEARCH_LEAD_AUDIT. Design push 0f93b162b15cd89340a6a2fba9ac26411c0c5fb2 preceded all B1 depth evaluations. Complete six historical configurations (S7 three-frequency and S8 four-frequency, three depths) are preserved. H0 unique truth minima 6/6; H1 within10 32/48, rank<=3 26/48, worst error 50.0 m. Proposed progress 0%; confirmed R4 remains 0%. No full 10% credit is awarded before lead audit.

CZ certified global-search route CLOSED; A1 BLOCKED / NOT_COMPLETED (0/15%), acquisition implementation not established. Conditional depth mechanism is separate from cold start and end-to-end performance. [B1 report](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_REPORT.md) and [decision](../R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY/B1_DECISION.json). B2 and A2/A3/A4/B3/B4/C/P5 NOT_OPENED. Stop after commit/push.


## Accepted conditional B1 and authorized B1A range-speed boundary

Research lead accepts3fc3518161722890b7aa8da7fab524fea7b8f917 as B1_DEPTH_IDENTIFIABLE_ONLY_UNDER_TIGHT_HORIZONTAL_CONDITIONING_ACCEPTED: H0 exact mechanism established, horizontal tolerance envelope not established. B1 remains0/10%, overall R4=0%; no progress credit. Original B1 remains frozen.

New explicit R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY stage: retain all six S7/S8 depth configurations and unchanged per-window RMS score/observations/E0/150:5:250 depth labels. Theta0/psi5 fixed; mechanically inherited dyadic offsets, all81 r-v signed states per configuration,486 profiles/10206 rows. Fixed target +/-0.25 km,+/-0.05 m/s: all150 profiles must have depth error<=10 m AND true rank<=3. Smaller robust rectangles/Pareto frontier are descriptive and never rescue failed target. No continuous interior or sensor-accuracy claim. [Design](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_DESIGN_FREEZE.json), [Gate](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_GATE.md), [lattice](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_LATTICE.csv). Commit/push design and verify remote before any B1A profile; execute once unchanged, push then stop for audit. B2/A2/A3/A4/B3/B4/C/P5 NOT_OPENED. CZ acquisition route remains CLOSED; A1 NOT_COMPLETED; R4_PROGRESS.json unchanged.


## B1A range-speed conditioning execution checkpoint

B1_HORIZONTAL_CONDITIONING_ENVELOPE_BELOW_TARGET; PENDING_RESEARCH_LEAD_AUDIT. Design c3a9744f29ddc657c147539f135ccbecac5a3c1a pushed before all486 new profiles;10206 score rows retained. Frozen target +/-0.25 km, +/-0.05 m/s: 13/150 pass, 137 failures. Proposed B1/overall progress=0%; confirmed0%. Largest tested robust rectangle(s)=[]; this cannot substitute for the primary target. Nonmonotonic witnesses retained. [B1A report](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_REPORT.md), [decision](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_DECISION.json). Finite lattice, fixed theta/psi, matched environment only; no continuous safety/engineering sensor/end-to-end claim.

A1 NOT_COMPLETED, CZ certified-global-search route CLOSED. B2 and A2/A3/A4/B3/B4/C/P5 NOT_OPENED. No post-result lattice additions, score/tau/grid changes, SSP mismatch or numerical rerun. Stop after commit/push for independent audit. R4_PROGRESS.json remains historical and unchanged.


## Accepted negative B1A and authorized B1B saved-data diagnostic

Research lead accepts8a60aeb99297c57e574daf9b796a08114b9c1f3d as B1_HORIZONTAL_CONDITIONING_ENVELOPE_BELOW_TARGET_ACCEPTED. Exact-horizontal mechanism supported; plug-in r/v not robust in tested lattice. B1 NOT_COMPLETED, B1=0/10%, overall0%. No finer steps to seek a passing envelope; B2 not opened. Original B1/B1A unchanged.

Authorized R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT only: saved B1A scores, all six configurations and16 original rectangles. MIN, fixed empirical-nearest-rank Q25, MEDIAN plus preregistered secondary MIN_NO_CENTER. Target route diagnostic D1/D2 MIN6/6 rank<=3/within10m; D3 Q25>=5/6 within10m and>=2/3 S7. No confirmation Gate or B1 completion credit regardless outcome. Every aggregate traced to an old score row; no new propagation/score/node/depth/case. [Policy](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/B1B_DIAGNOSTIC_POLICY.json), [definitions](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/B1B_DIAGNOSTIC_DEFINITIONS.md). Commit/push policy before analysis, then closeout/results commit/push and stop for independent review. Negative diagnostic closes current depth engineering route; positive only warrants separately instructed fresh B1R design. B1R/B2/A2/A3/A4/B3/B4/C/P5 NOT_OPENED; A1 NOT_COMPLETED, CZ acquisition CLOSED; R4_PROGRESS.json unchanged.


## B1B saved joint r-v-z diagnostic and B1 closeout

JOINT_RVZ_PROFILE_ROUTE_NOT_SUPPORTED_BY_EXISTING_EVIDENCE; CURRENT_CONDITIONAL_DEPTH_ENGINEERING_ROUTE_CLOSED; PENDING_RESEARCH_LEAD_AUDIT. Policycb9c0cfa93bf829ca9d8f6312ad14c73283e038e pushed before saved-score analysis. Target MIN/Q25/MEDIAN/no-center summaries={"MIN": {"n_cases": 6, "exact_depth": 6, "within_10m": 6, "rank_le3": 6, "practical_pass": 6, "worst_error_m": 0.0, "worst_true_rank": 1, "minimum_second_margin_db": 0.6524505796925171}, "Q25": {"n_cases": 6, "exact_depth": 0, "within_10m": 3, "rank_le3": 0, "practical_pass": 0, "worst_error_m": 30.0, "worst_true_rank": 13, "minimum_second_margin_db": -0.30153222811184976}, "MEDIAN": {"n_cases": 6, "exact_depth": 0, "within_10m": 1, "rank_le3": 0, "practical_pass": 0, "worst_error_m": 25.0, "worst_true_rank": 16, "minimum_second_margin_db": -0.659365582983999}, "MIN_NO_CENTER": {"n_cases": 6, "exact_depth": 0, "within_10m": 4, "rank_le3": 0, "practical_pass": 0, "worst_error_m": 45.0, "worst_true_rank": 17, "minimum_second_margin_db": -0.7028222103223405}}. New propagation/depth scores/nodes/cases=0; all16 rectangles and six configurations retained. No B1 completion credit, B1=0/10%, R4=0%.

[B1 closeout](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/R4_B1_CLOSEOUT.md), [diagnostic report](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/B1B_REPORT.md). Exact-horizontal mechanism retained; plug-in range/speed not robust, angles stable only at tested scales. B1R/B2/A2/A3/A4/B3/B4/C/P5 NOT_OPENED. A1 NOT_COMPLETED, CZ certified global acquisition CLOSED. Stop after commit/push. R4_PROGRESS.json unchanged.


## Accepted B1B closure and authorized auxiliary-node geometry Gate0

Parent6b8d57ef6082f24720c9ece654b9cd710a687e13 accepted negative B1B; current conditional-depth engineering route CLOSED. Single-array cold-start acquisition NOT ESTABLISHED, CZ certified global route CLOSED; exact-horizontal depth mechanism RETAINED_AS_ORACLE_EVIDENCE. Five-parameter target supported=false. R4 original A/B sequence PAUSED_PENDING_ARCHITECTURE_RESET; R4=0%, B1=0/10%. Mature continuous bearing information retained.

R4_ROUTE_RESET_AUXILIARY_NODE_GATE0 authorized: saved evidence closeout plus one mobile passive node known position/directed bearing geometry only. [Design freeze](../R4_ROUTE_RESET_AUXILIARY_NODE_GATE0/AUX_GATE0_DESIGN_FREEZE.json). All original range/baseline/crossing/noise cells retained; physical triangle feasibility explicit; both roots/mirrors; fixed20000 trials seed2026100501. Design commit/push/remote verification before any geometry run. No propagation/depth/optimizer/tracker/TDOA; no automatic new R4-A1/depth release. A2/A3/A4/B2/B3/B4/C paused/not opened, P5not opened. Execute once, commit/push/verify then stop. R4_PROGRESS.json unchanged.


## Auxiliary-node geometry Gate0 execution

AUXILIARY_NODE_GEOMETRY_CONDITIONALLY_PRACTICAL; PENDING_RESEARCH_LEAD_AUDIT. Known-position paired-directed-bearing geometry only. Physically impossible grid cells retained as infeasible. [Report](../R4_ROUTE_RESET_AUXILIARY_NODE_GATE0/AUX_GATE0_REPORT.md). R4=0%; no automatic architecture, motion or depth restart. Commit/push/verify then stop.


## Accepted AUX Gate0 and authorized Gate1 cross-track requirement

Parent79bc6b5c04e08d6074719e3c045947c1325253c8 accepted conditional geometry, not engineering capability. Single-array cold start NOT ESTABLISHED; depth engineering CLOSED, oracle mechanism retained; AUX1 CONDITIONAL_ARCHITECTURE_CANDIDATE. R4=0%, original A/B sequence paused.

R4_AUX_GATE1_CROSSTRACK_BASELINE_BEARING_REQUIREMENT authorized: MAIN(0,0), AUX(0,+/-B), B2:0.5:10km, R50:1:60km, effective directed-bearing RMS0.05/0.075/0.1/0.15/0.2deg;50000 fixed common Gaussian trials per cell. [Design](../R4_AUX_GATE1_CROSSTRACK_REQUIREMENT/AUX_GATE1_DESIGN_FREEZE.json). Freeze/push/remote verify before MC/tests, execute once unchanged, independent audit/results commit/push then stop. Full-range means11 tested ranges, not continuous guarantee. Actual hardware accuracy UNKNOWN; no numerical CRLB capability promise. No propagation/depth/SSP/TDOA/tracker/5D/off-gridA1/position-bias-time scenarios. R4-A1-NEW/Gate2 not opened; A2/A3/A4/B2/B3/B4/C paused/not opened; P5not opened; R4_PROGRESS.json unchanged.


## AUX Gate1 requirement execution

AUX1_BASELINE_BEARING_TRADEOFF_ESTABLISHED; PENDING_RESEARCH_LEAD_AUDIT. Fixed cross-track deployment, discrete11-range requirement frontier; no hardware or continuous-range guarantee. [Report](../R4_AUX_GATE1_CROSSTRACK_REQUIREMENT/AUX_GATE1_REPORT.md). Single-array cold start not established, depth closed, AUX1 conditional, R4=0%; no R4-A1-NEW/Gate2 restart. Commit/push/verify then stop.


## Accepted AUX Gate1 and authorized Gate2A single-factor nonideal budget

Parent1fde9dac9d848e5c8b9987b010182411fdbc79c6 accepted tradeoff established, ideal geometry only. Single-array cold-start not established, depth engineering closed; AUX1conditionalarchitecturecandidate, hardwareunknown; R4=0%, originalA/Bsequencepaused.

R4_AUX_GATE2_NONIDEAL_MEASUREMENT_ERROR_BUDGET authorized Gate2A only. Anchors5km/.05deg and7km/.075deg, B6/.1deg historicalreferenceonly. Three SINGLE-FACTOR families: complete81signedbias, per-axisnavigation7levels, knownactualdeployment9signedangles;11rangesandmirrors.30000fixedcommontrials seed2026100503. [Freeze](../R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET/AUX_GATE2A_DESIGN_FREEZE.json). Freeze/push beforeMC/tests, executeonce, independentaudit, commit/push/verifythenstop. No jointstress/time-skew/associationMC/propagation/depth/SSP/TDOA/tracker/5D. Client requirements do not claim simultaneous axis tolerances. Actualhardwareunknown. R4-A1-NEW/Gate2Bnotopened;A2/A3/A4/B2/B3/B4/Cpaused/notopened;P5notopened. R4_PROGRESS.jsonunchanged.


## AUX Gate2A execution

AUX1_NONIDEAL_REQUIREMENTS_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED; PENDING_RESEARCH_LEAD_AUDIT. Single-factor signed-bias/navigation/deployment budgets only, not a jointly validated point. [Report](../R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET/AUX_GATE2A_REPORT.md). Actual hardware unknown; R4=0%; R4-A1-NEW/Gate2B not opened, depth closed. Commit/push/verify then stop.


## AUX Gate2A acceptance and Gate2B design

Parent5d8b9c9f2d8efdb74a95da4ea546aae11457a24f accepted SINGLE_FACTOR_REQUIREMENTS. Gate2B explicitly authorized as final geometry-only joint static qualification; two families and fixed half point. [Criteria](../R4_AUX_GATE2B_JOINT_ERROR_BUDGET/AUX_GATE2B_CRITERIA.md). Actual hardware UNKNOWN; R4=0%; A1-NEW closed. Push freeze before MC; stop after execution push.


## AUX Gate2B execution

AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED; PENDING_RESEARCH_LEAD_AUDIT. Joint static ladders complete. [Report](../R4_AUX_GATE2B_JOINT_ERROR_BUDGET/AUX_GATE2B_REPORT.md). Hardware UNKNOWN; time/association not validated; R4=0%; A1-NEW not opened. STOP numerical architecture work; client confirmation next.


## AUX Gate2B acceptance and client confirmation freeze

Parent333ead6461a7574a469495033d48fd89682aa496 independently accepted JOINT_STATIC_REQUIREMENT and HALF_T5_JOINT_BUDGET. Strongest tested T5 A lambda=.50; B=.75; T10 A/B=.75. [Client package](../R4_AUX_CLIENT_CONFIRMATION/GPT_SYNC.md). NUMERICAL_ARCHITECTURE_WORK_STOPPED; CLIENT_CONFIRMATION_REQUIRED. Actual hardware UNKNOWN; time/association not numerically validated; R4-A1-NEW NOT_OPENED; depth CLOSED; R4=0%. No new experiment and no automatic Gate2C; client evidence and research-lead authorization required before any next stage.


## Application pre-research authorization and A1-new PRE-RUN

PROJECT_APPLICATION_PRE_RESEARCH, not hardware acceptance. Newlead instruction supersedes prior numericalclient-confirmationSTOP for this authorizedA1-newstage. HardwareUNKNOWN is nonblockingdesignassumption. Scientificbaseline333ead6 retained; chronologicalparent0eef9ff retained. PrimaryAnchorA,secondaryB;12off-gridtruths x500each;two-stage observationonlyfourparameterdynamicestimator. [FrozenGate](../R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_GATE.md). Credit15onlyifprimaryPROJECTallcasescloses; otherwise0. Alllegacyfrozenartifactsandprogresspreserved; newactiveapplicationprogresswillbesavedafterexecution. Depthnotopened. PushfreezebeforeMC,completebothanchors,pushthenSTOPforlead audit; noautomaticA2.


## R4-A1-new execution

R4_A1_NEW_BASELINE_NOT_ESTABLISHED; PENDING_RESEARCH_LEAD_AUDIT. APPLICATION_PRE_RESEARCH design assumptions; hardware UNKNOWN nonblocking. [Active application progress](R4_APPLICATION_PROGRESS.json). A1-new/overall scientific credit=0%. [Report](../R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_REPORT.md). B secondary only; depth not opened. Commit/push/verify then STOP for audit; no automatic A2.


## A1-new accepted speed bottleneck and diagnostic freeze

Parent e39182ac82b988e71b79d73a06a97076ca76f6bf accepted R4_A1_NEW_BASELINE_NOT_ESTABLISHED_ACCEPTED. Under frozen application panel:range/bearing STRONG12/12;heading PROJECT12/12,STRONG11/12;speed PROJECT0/12. Internal3/4horizontalmetricsestablished,speedbottleneck;notwholeA1credit. NewexplicitDEVELOPMENTdiagnostic D0-D4 reuses12000savedrealizations; no newdraws/panel/baseline/time/turn/RMS/Gate. [Policy](../R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC/SPEED_DIAGNOSTIC_POLICY.json). Pushfreeze before re-solving. R4-A1/R4=0%;depthclosed. Fullanalyses,pushthenSTOP;noautomaticfresh/FIX2.


## Speed bottleneck diagnostic execution

SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO; DEVELOPMENT_ONLY, PENDING_RESEARCH_LEAD_AUDIT. [Report](../R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC/SPEED_DIAGNOSTIC_REPORT.md). [Active speed status](R4_SPEED_DIAGNOSTIC_PROGRESS.json). AcceptedA1range/bearingSTRONG,headingPROJECT; speedunresolved. Reused12000scenes; newMC0; R4-A1/R4=0%. Depthnotopened. PushthenSTOP; freshonlyifseparatelyauthorized.


## Speed information efficiency design freeze

Accepted prior 1200s bottleneck remains scoped to frozen formation/noise/tested estimators. Reuse saved D3 scenes for E0 identity, E1 Gaussian L2, E2 oracle initialization, four-state known-sigma local FIM. [Policy](../R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY/SPEED_INFO_AUDIT_POLICY.json). Push before all resolves/FIM tests. No new scenarios; R4-A1/R4=0%; depth NOT_OPENED; stop after execution commit.


## Speed information efficiency execution

CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB; DEVELOPMENT ONLY, pending independent audit. [Report](../R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY/SPEED_INFO_REPORT.md). [Active status](R4_SPEED_INFORMATION_PROGRESS.json). New MC=0; R4-A1/R4=0%; depth NOT_OPENED. Next is a recommendation only; push then STOP.


## Observability frontier research reset: literature and designs

Accepted parent ee656d49bb7459268a991db6c44ba50d661456db, scoped local/model finding. Same-condition bearing loss/initialization/filter repair CLOSED. PROJECT_APPLICATION_PRE_RESEARCH; hardware UNKNOWN nonblocking. [Report](../R4_OBSERVABILITY_FRONTIER_RESEARCH_RESET/GPT_SYNC.md), [dual status](R4_OBSERVABILITY_FRONTIER_PROGRESS.json). Original joint A1/R4=0%, speed NOT_PASSED. Completed pre-research evidence separately listed, no invented percentage. Literature coverage, TOP3 and TOP2 two-layer designs delivered; scientific calculation, new MC, acoustic synthesis and propagation NONE. FDSL/POMAP full-text and strict comparator formulas pending. Independent audit pending. Commit/push/verify then STOP; no automatic TOP2, extended time, maneuvers, depth, A2, SSP or P5.


## E1 primary acoustic source and formula lock

Parent 0db30f60e781676f9168d380678693a32f1325a1 frontier design ACCEPTED_WITH_E1_PRIMARY_SOURCE_GAP. [E1 lock](../R4_OBSERVABILITY_FRONTIER_RESEARCH_RESET/GPT_SYNC_E1_LOCK.md). Four full PDFs and methods recovered; Sun measurement Eq4 locked, complete strict MFB-AUKF baseline conflicts unresolved. Xu pressure-to-radial-velocity conditional chain locked; CZ phase/sign/covariance not established. A-F structural checks and revised B0/LIT-F/LIT-RV/UUV/NEW contracts frozen. E1_PRIMARY_SOURCE_OR_MODEL_LOCK_INCOMPLETE; NOT_READY. 11 fixed algebra checks only; new scientific/MC/signal/propagation/project-FIM runs=0. E2 NOT_OPENED; FDSL/POMAP comparator pending. R4=0%. Commit/push/verify then STOP; no automatic scientific execution.


## E1-G0 frequency mechanism screen

Source measurement lock ACCEPTED_WITH_ORIGINAL_ESTIMATOR_EXCEPTION at f9b994a; G0 admitted, G1 not admitted. Design46242899 frozen/pushed before execution. [Report](../R4_E1_G0_FREQUENCY_INFORMATION/E1_G0_REPORT.md). 24 nominal mirror scenes,1896 deterministic information records; no MC/audio/propagation. E1_G0_FREQUENCY_INCREMENT_NOT_ESTABLISHED: best registered 3-line1mHz N2 C0 median variance reduction0.125681%, below20%; N0/N5 also below gate. All corners/packages retained. Numerical controls,finite differences,dense covariance reconstruction pass. Mechanism screen completed, no P95 claim; R4-A1/R4=0%. Recommended E2_REVIEW only; G1/E2/depth/A2/SSP/P5 unopened. Push execution commit then STOP for independent audit.


## E2-G0 physical design freeze

E1-G0 accepted at43fa1f57 and closed, no G1. New explicit custom HLA autoproduct observability screen. [Freeze](../R4_E2_G0_HLA_DIFFERENCE_INFORMATION/E2_G0_DESIGN_FREEZE.json):24 scenes, same1200s trajectory, three physical snapshots, two new8×2m HLAs, 150–250Hz exact-frequency modal provider, source/gain/depth nuisance profiling and shared-pair covariance. Physics admission and bounded computation precede information claims. No new MC/audio; R4=0%; stop after execution push.


## E2-G0 physical screen execution

E2_G0_PHYSICS_OR_COVARIANCE_INCOMPLETE. [Report](../R4_E2_G0_HLA_DIFFERENCE_INFORMATION/E2_G0_REPORT.md). One bounded deterministic attempt; null and physical admission records retained. No MC/audio; R4=0%. STOP pending lead independent audit; no automatic extraction/depth/A2/SSP/P5.


## R4 E2 forward fidelity diagnostic
E2-G0 independent audit accepted at4605cf2; original admission FAIL. Authorizes bounded error attribution only; four third-grid calls; no information/MC/audio. R4=0%. Freeze before execution; commit/push then stop after result.

E2 forward diagnostic completed: E2_CA_AND_RAW_FORWARD_UNSTABLE. Original FAIL_UNCHANGED; information NOT_EVALUATED. Four bounded calls; independent audit; STOP. R4=0%.


## E2 paired-frequency bounded recovery
Diagnostic independent audit accepted at e6f1d68f. Authorizes9 registered pairs/13 frequencies; new80001x9 and160001x13, max22 calls. No information/MC/audio. Original admission FAIL_UNCHANGED, R4=0%. Freeze then push before once-only execution; stop after result.

E2 paired recovery completed: E2_PAIRWISE_CA_NUMERICAL_RECOVERABILITY_SUPPORTED. Original admission FAIL_UNCHANGED; information NOT_EVALUATED; R4=0%. Suggested FULL_BAND_CERTIFICATION_REVIEW only; new stage NOT_AUTHORIZED. STOP.
