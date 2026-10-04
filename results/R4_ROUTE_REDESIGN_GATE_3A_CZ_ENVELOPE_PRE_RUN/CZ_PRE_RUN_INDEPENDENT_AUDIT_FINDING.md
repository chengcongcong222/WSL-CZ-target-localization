# Independent audit of e9316c031d40f11fc07772959d439b965d876cd5

Research-lead decision: **CZ_ENVELOPE_PRE_RUN_FREEZE_NOT_ACCEPTED**. No 21-case execution is released. R4 remains 0%.

Finding A: the historical engine used `closed = valid and not any(r.status == 'RETAINED_UNRESOLVED' for r in records)`. A terminal possible outer cell therefore could not close coverage. This confuses a completed conservative outer partition with a proof that every retained cell is uniformly compatible.

Finding B: the frozen TINY fixture, about 1e-12 wide on each state axis, has feature upper 3.19847126251623e-6 dB, above unchanged tau_F=1e-7 dB. The upper certificate is too loose to certify uniform compatibility; this does not invalidate the lower-bound rejection mechanism.

Finding C: inherited dyadic split counts are coarse (6,4,4,4), D=18; fine (7,5,5,5), D=22. The lead identified approximately 36/44 requests to reach a terminal feasible path. Counting the root and classification of the last sibling gives exact necessary requests **37/45** for complete queue exhaustion in a single unresolved-path tree: D path internals + D rejected siblings + one terminal retained leaf. The old caps 14/29 are below both versions of the lower bound. This is a pre-development infrastructure contradiction, not a measured difficulty of the 21 cases.

The mathematical lower-bound enclosure, truth isolation, Q2/N0, tau_F, NOMINAL 5 km/component/0.5 km thresholds and primary manifest remain accepted and unchanged. The rejected checkpoint remains historical; its source files and artifacts are preserved byte-for-byte. Correction uses a new engine/source and appended artifacts, never rewrites the old READY decision as accepted.

The correction distinguishes domain_partition_closed, compatibility_certified and search_budget_closed. Terminal possible boxes belong to conservative outer support. Invalid/capped regions remain retained with failed closure. Rejection still requires the original certified whole-cell lower bound; no sample or heuristic replaces it.

The policy registers a domain corner and the domain midpoint, their fixed depths, matched observation generation, a fixture-only inherited noiseless cutoff, the same grids and a finite structure/runtime/storage-derived ladder. If the maximum cap fails on these fixtures, record SEARCH_CLOSURE_NOT_ESTABLISHED and stop; no C4/F4 or development run.
