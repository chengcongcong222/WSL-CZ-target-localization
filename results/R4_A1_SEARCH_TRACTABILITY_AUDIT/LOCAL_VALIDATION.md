# Pre-run local validation

Scope: PRE-RUN FREEZE only; no noisy development cases or fresh confirmation.

After METHOD_FREEZE.json creation, `python -m pytest tests/test_r4_a1_search_tractability.py -q` passed **31 tests**. Exact command/output and timestamp are saved in PYTEST_PRE_RUN_FREEZE.json and .log. Tests cover admitted N / refused N+1 boundaries, cache-hit request accounting, native SHGO/DIRECT exception propagation, SHGO local-refinement requests, repeated refused callbacks as enforcement failure, unrelated SystemError rejection, direct-modal self-match and repeated exact identity, saved geometry, bearing/physical feasibility, effective-depth mapping, catalog separation and static/runtime privileged-information exclusion. Runtime sentinels cover both solvers, extraction and exact verification. Self-match is a synthetic isolated unit test, not scientific recovery evidence.

API probes: both solvers at caps 1,3,7,20,64,128 (12 records), plus two repeat physical compatibility probes per solver at cap 32 (4 records). Each admits exactly N requests, logs exactly one refused N+1 attempt, and performs no forward on refusal. SHGO exposes ObjectiveBudgetExceeded directly; DIRECT exposes SystemError with ObjectiveBudgetExceeded as its cause. Budget exits are solver_success=false. This causal behavior is disclosed for independent audit.

Cost calibration: 84 fixed synthetic cost-only requests, all 84 exact forward calls, across 21 inherited depth labels; no noisy observation or historical recovery result is loaded. Four physical compatibility catalogs are reproducible and contain evaluated witnesses, not certified local minima.

`python r4_a1_search_tractability_audit.py audit` independently passed **227 integrity/reconstruction checks**: protected-input/frozen-byte hashes, cost-derived budget policy, all four counters from raw ledgers, exception semantics, saved-state geometry/bearing/direct scores and exact ranking. FROZEN_INPUT_HASHES.csv rechecks 2,989 baseline script/artifact files. New Git attributes preserve frozen bytes on Windows.

These results validate pre-run implementation and accounting. No scientific tractability Gate, recovery rate, global uniqueness or environment robustness beyond the recorded installed environment is claimed. R4=0%; A2/depth/SSP/P5 unopened.
