# Development execution harness freeze

Checkpoint: **DEVELOPMENT_EXECUTION_HARNESS_FREEZE_READY_FOR_INDEPENDENT_AUDIT**. The research lead accepted the parent as PRE_RUN_INFRASTRUCTURE_FREEZE_ACCEPTED. Scientific tractability remains NOT_EVALUATED. This checkpoint completes only the driver/semantics/Gate freeze and stops for the second independent audit.

Parent PRE-RUN SHA: `b4839b295777127ec0c8ade56b76c776db99148b`. The commit containing this report is the new harness-freeze SHA; it cannot embed its own SHA. Local HEAD and GitHub main are verified equal after push.

| Binding | SHA256 |
|---|---|
| frozen budget | `40cd9c0b2cbc7a6fced83faf97009c25d652223bf7cbb77e7ebc6a981089672d` |
| core objective/controller | `0890fd3c5df680a2d56e13d04709dfe727d09d3b44aabfe4b539cab21de92e2b` |
| development driver | `9ceb0b0b4e5aae9462e3d2999b3566921a58d0839550b0aa7676f3ee70632a34` |
| Gate logic | `b5d205ce6b852fe9941e86b648a4e962a13308fa1fc36937ac28525afce4c66b` |

Core and budget bytes are unchanged. SHGO and DIRECT retain 16/64/256 admitted requests per case/depth, the same options and all accepted tolerances. The new driver is r4_a1_search_tractability_development.py; comparison logic is r4_a1_search_tractability_gates.py; metadata projection/independent checks are r4_a1_search_tractability_harness_audit.py. Every real scientific entry is blocked by DEVELOPMENT_RELEASED=False before loading runtime observations/models. Fresh generation has its own false release flag.

DEVELOPMENT_CASE_MANIFEST.csv: **21/21** unique nominal cases; privileged columns absent; explicit observation case_ids mapping true. The driver runtime loads only the permitted manifest and immutable CASE_OBSERVATIONS.npz under the accepted FIX1 directory. No truth panel, full old result/error dataframe, old recovery/J/state input or future fresh panel is supplied to the optimizer. Manifest selection and archive hashes/provenance are independently reconstructed.

Raw witnesses retain **EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM**. Each raw branch independently creates engine, controller, cache and optimizer call; each case uses all 21 branches and selects the direct-modal raw top witness. Recovery is solely its strict J_exact<0.001 threshold hit. Budget exhaustion and solver_success=false are separate from recovery; execution invalidity separately blocks Gate PASS. All feasible distinct exact states will be retained in *_EVALUATED_WITNESSES.csv before score-first, same-depth frozen witness dedup.

Formal comparison is threshold-hit / discrete witness-cluster convergence, not local stationarity. Raw T2/T3 comparison requires valid recovered top states/depth and T3->T2 cluster containment; dual T3 comparison requires both 21/21 raw-convergence families first, then top state/depth and bidirectional cluster matching. No certified minima files are manufactured. Strict alias reconstruction retains its narrower joint-match criteria and finite-catalog scope. Raw/dual failure cannot reach fresh generation, increase budgets, change tolerances or add a repair algorithm.

Validation after final freeze: **74 tests passed; 312 integrity/reconstruction checks passed**. See LOCAL_VALIDATION.md and raw logs. Protected baseline identities, budget SHA, core SHA and all new freeze hashes match. The parent accepted pre-run method/manifests/report remain archived verbatim; no scientific result has been overwritten.

Noisy development cases: **0**. Development released: **false**. Scientific Gates: **NOT_EVALUATED**. R4: **0%**. A2/depth/SSP/P5: **UNOPENED**. Stop after commit/push; require a second independent audit and separate 21-case release.
