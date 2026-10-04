# Frozen 21-case development result

Decision: `A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED`. R4 remains 0%. Fresh confirmation was not authorized and did not run. A2, depth development, SSP and P5 remain unopened. This stage stops after submission.

Accepted parent harness: `c5c9ba2e1258e73486e62df26499e03cce204f27`. Runtime release authorization: `d731cf51d72b41dea7fd4c7a48fcf21a623f7bc4`. The source release flags remain false. Algorithm, objective, budgets, solver choices, manifest, Gate logic and tolerances remain byte-identical to the accepted harness. No failed case was retried or removed.

The single frozen execution completed 21 noisy cases x two solvers x three independent raw budgets = 126 case-budget runs, each with 21 inherited depth branches = 2646 branches. All SHGO runs preceded all DIRECT runs. Every branch used a fresh engine, controller and cache with admitted-request caps 16/64/256. Recovery means only that the best feasible exact evaluated witness has J < 0.001 dB. These are inherited synthetic matched E0 controls with perfect relative-TL observations and noisy bearings; the findings have that model and panel scope.

| Solver | T1 threshold hits | T2 threshold hits | T3 threshold hits | Raw T2/T3 convergence |
|---|---:|---:|---:|---:|
| SHGO | 0/21 | 0/21 | 0/21 | 0/21 |
| DIRECT | 0/21 | 0/21 | 0/21 | 0/21 |

Raw budget convergence Gate: `False`. Dual solver agreement: `NOT_REACHED`. Dual Gate: `NOT_REACHED`. The dual comparison is conditional on the raw Gate and cannot receive a scientific PASS from informal partial comparisons.

Exact alias: `NO_IN_TESTED_CATALOGS`. This finding applies only to exported tested catalogs and does not establish global uniqueness. All candidates retain `EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM`; neither optimizer capture basins nor local/global minima are certified. Failure of these capped searches does not establish physical non-identifiability or universal computational intractability.

Optimizer admitted objective requests: 296352; exact forward evaluations: 279114; unique physical states: 292946; cache hits: 3406; refused requests: 2646. Mandatory driver cold validation made 279114 further exact forward evaluations outside the optimizer cap. The frozen alias prefilter selected 0 witnesses for additional driver exact forward reconstruction. Hard-budget enforcement failures: 0. Cap exhaustion is an expected accounting event and does not itself declare search or execution failure.

Mechanical reconstruction passed 14479 checks and failed 0. It reads all 2646 request ledgers, independently reconstructs counters, physical-state hashes, every retained feasible witness, representatives, branch/case aggregates, raw convergence, conditional dual evaluation and strict joint-alias status. Independent geometry and cold direct-modal RMS reconstruction covered 232973 unique case/depth/state combinations in 14754 batched forward calls; maximum J difference = 1.16196030575e-09 dB and maximum bearing cost difference = 4.33680868994e-18 rad^2. This audit-only memo was created after all searches finished and was never supplied to a solver.

Execution-invalid raw runs: 0. Failed recovery raw runs: 126. Failed convergence comparisons: 42. Their full case lists are preserved in DEVELOPMENT_EXECUTION_INVALID_CASES.csv, DEVELOPMENT_FAILED_RECOVERY_CASES.csv and DEVELOPMENT_FAILED_CONVERGENCE_CASES.csv. Raw results, request ledgers, full evaluated witnesses, representative catalogs, exception chains and numerical failures remain preserved.

Post-run unchanged harness tests: ........................................................................ [ 97%]
..                                                                       [100%]
74 passed in 4.45s

Actual execution workspace: `D:\CodexExperiments\WSL-CZ-tractability-20261004`. The original C drive had no free space, so execution used a verified independent D-drive checkout with all frozen identities preserved. The original C checkout requires a separate disk-space check before synchronization. Historical accepted reports and manifests remain intact; this report and DEVELOPMENT_FINAL_DECISION.json are new evidence for the third independent audit.
