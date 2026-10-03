# Research-lead sync: PRE-RUN FREEZE only

Baseline: 7ab24845e6e1551b75287fefb1ab662e92b395b8.

Checkpoint decision: PRE_RUN_FREEZE_READY_FOR_INDEPENDENT_AUDIT; scientific tractability NOT_EVALUATED. The commit containing this file is the independent pre-run freeze commit, before any noisy development result. Use Git HEAD/remote main for its full SHA; a commit does not embed its own SHA.

Please audit TRACTABILITY_BUDGET_FREEZE.json, METHOD_FREEZE.json, HARD_BUDGET_TEST_RESULTS.csv, raw API_REQUEST_LOGS and the three new code/test files. Both solvers use raw per-case/per-depth caps 16/64/256 with fresh state/cache. Boundary tests admit N, refuse N+1, count cache hits and record all four counters. DIRECT's native SystemError preserves ObjectiveBudgetExceeded as cause and exits on the first refusal; this is explicitly disclosed, and exhaustion remains solver_success=false.

Validation after freeze: 34 new tests passed; 228 independent integrity/reconstruction checks passed; 2,989 protected baseline files unchanged. Cost/API calibration only: 84 exact cost requests and 16 native compatibility/cap probes. No noisy development cases (0); no fresh confirmation; no scientific Gates claimed PASS.

R4=0%; A2/depth/SSP/P5 unopened. Execution stops at the checkpoint, awaiting independent audit and separate release. There is no automatic continuation.

Pre-run portability correction: the initial freeze commit is `5fa63092007ee10301b3d24b7378aef1333232ff`. Main checkout had 124 historical text files with CRLF/LF-only differences. The correction preserves all historical files and original raw hashes, adds frozen baseline Git-blob/text normalization evidence, and tests conversion acceptance plus content/binary rejection. Frozen budgets, exact objective and zero-development stopping scope are unchanged. The subsequent correction commit is the current checkpoint HEAD.
