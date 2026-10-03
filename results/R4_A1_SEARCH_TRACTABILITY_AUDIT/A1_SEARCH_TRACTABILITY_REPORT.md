# A1 search tractability: pre-run checkpoint

Decision: **PRE_RUN_FREEZE_READY_FOR_INDEPENDENT_AUDIT**. The scientific tractability decision is **NOT_EVALUATED**. Stop here for the research lead's independent audit and separate development release.

Baseline: `7ab24845e6e1551b75287fefb1ab662e92b395b8`. The isolated worktree began clean with HEAD==origin/main==baseline. The main workspace's pre-existing unrelated untracked files were preserved. R3/A1/FIX1/FIX2 numerical artifacts are unchanged; FIX2 remains accepted and frozen as A1_FIX2_BLOCKED_BY_UNCLOSED_ACOUSTIC_COVERAGE.

## Frozen implementation and budgets

New files: r4_a1_search_tractability.py, r4_a1_search_tractability_audit.py, tests/test_r4_a1_search_tractability.py. Solver choices remain SciPy SHGO and DIRECT; only observation-conditioned (r0,u_r) and inherited finite nuisance-depth branches are allowed. Exact scores use the direct modal model. The code has no enabled noisy-development driver.

| Solver | Raw T1 | Raw T2 | Raw T3 | Scope |
|---|---:|---:|---:|---|
| SHGO | 16 | 64 | 256 | admitted objective requests per case and depth branch |
| DIRECT | 16 | 64 | 256 | admitted objective requests per case and depth branch |

The 21 inherited depth labels are uniformly enumerated. Per-case aggregate upper bounds are 336/1344/5376 requests per solver; each branch is independent. The cost-only policy measures p95=0.0213035749912 seconds per exact request and selects base=16 from the declared six-hour/factor-two formula. Across both 21-case development families and all raw budgets, the maximum admitted-request count is 296352; the estimated exact-plus-validation time is 3.507 hours, excluding I/O/catalog overhead. This is an estimate, not a runtime guarantee. Budgets were not selected from recovery results.

TRACTABILITY_BUDGET_FREEZE.json records every option, threshold, tolerance, depth label, domain and calibration hash. METHOD_FREEZE.json binds implementation, tests, complete controlling instructions/design, method, budget, input hashes and raw API ledgers. All code/design hashes were checked after freeze. Recovery remains strictly J_exact<0.001 dB.

## Hard cap and evidence

Each native hard-cap probe admits exactly N requests and refuses N+1 before cache lookup or modal evaluation. Cache hits remain counted requests. The four required counts are retained and reconstructed independently. There is no admitted-request or exact-forward overshoot.

Important native API detail: SHGO propagates ObjectiveBudgetExceeded directly; SciPy 1.15.3 DIRECT wraps it in SystemError with the original exception retained as __cause__. Both leave the native call at the first refused callback. Native classes/causal chains are logged, no library patch/workaround is used, unrelated SystemError is re-raised, and cap exhaustion never becomes solver success. The lead should review this explicit causal-propagation interpretation at the checkpoint.

34 new unit tests passed after freeze; 228 independent pre-run integrity/reconstruction checks passed. Twelve API-only cap probes, four physical compatibility probes, and 84 exact cost-only requests are preserved. Compatibility catalogs are evaluated witnesses explicitly not certified local minima; no basin-coverage inference is made from them.

## Mandatory stopping state

Noisy development cases executed: **0**. No development panel/results, fresh panel/observations, dual-solver scientific agreement, alias finding or six-Gate scientific decision has been generated. Scientific Gates remain NOT_EVALUATED at this checkpoint. No future stage is authorized by this report.

R4 remains **0%**. R3/A1/FIX1/FIX2 remain frozen. A2, depth/B, SSP and P5 remain unopened. The pre-run commit must be pushed and execution stops; development requires a separate research-lead instruction after audit.

Pre-run portability correction: the initial freeze commit is `5fa63092007ee10301b3d24b7378aef1333232ff`. Main checkout had 124 historical text files with CRLF/LF-only differences. The correction preserves all historical files and original raw hashes, adds frozen baseline Git-blob/text normalization evidence, and tests conversion acceptance plus content/binary rejection. Frozen budgets, exact objective and zero-development stopping scope are unchanged. The subsequent correction commit is the current checkpoint HEAD.
