# Corrected conservative outer-partition completion semantics

The research lead rejected the e9316c PRE-RUN execution design. Its mathematical lower-bound mechanism and scientific configuration remain unchanged. The correction lives in `r4_cz_envelope_support_corrected.py`; the rejected engine and every old artifact remain historical, not edited into an accepted result.

## Three independent concepts

**domain_partition_closed** means every part of the complete initial 4D domain has a leaf classification: certified excluded, or retained at the frozen terminal resolution. Possible terminal boxes are legal conservative outer support; a uniform compatibility proof is not necessary. The corrected implementation also subdivides uniformly compatible nonterminal boxes to the same frozen resolution, so every retained science leaf has the declared maximum widths.

**compatibility_certified** is an optional diagnostic. A cell's bearing upper bound must be below its unchanged cutoff, and its feature upper bound below tau_F for CZ. B0 uses bearing only. It says every point in the box is compatible with a permitted nuisance profile. It is stronger than outer-support membership and is never a condition for partition closure. The result-level diagnostic is true only when nonempty retained support consists entirely of such boxes.

**search_budget_closed** means natural exhaustion of the traversal queue without encountering a denied operation or storage transaction. Natural completion on the final admitted request is allowed. Clearing a queue after a cap or invalid bound is not exhaustion. A failed bound also invalidates execution.

The future Gate phrase `raw global coverage closure required` means exactly **FULL_DOMAIN_CONSERVATIVE_OUTER_PARTITION_CLOSED_AT_FROZEN_RESOLUTION**. It does not mean every retained cell is proven fully compatible.

## Status rules

- REJECTED_BY_CERTIFIED_BOUND: unchanged certain bearing_LB > cutoff, or completed all-depth feature_LB > tau_F.
- RETAINED_CERTIFIED_COMPATIBLE: terminal retained box with the uniform upper certificate.
- RETAINED_POSSIBLE_AT_TERMINAL: terminal retained box not certified excluded, with valid full-cell bounds; the upper certificate can be loose.
- RETAINED_UNRESOLVED_BUDGET: current and all pending boxes preserved after a cap; partition and search closure false.
- RETAINED_INVALID_BOUND: current and pending domain preserved on invalid/incomplete bounds; execution and partition closure false.

BOUNDARY_CENSORED remains an independent flag. Nonterminal non-excluded cells split with the unchanged breadth-first largest-normalized-width rule. Sampling, optimizer results and heuristic probabilities never supply a whole-cell exclusion. The original model/depth/bearing/feature bound sources are used byte-for-byte.

## Support, metrics and evaluation

Scientific W_hull, W_union, components, raw stability and post-seal retention use the union of **certified-compatible and terminal-possible** leaf boxes. A closed empty outer partition is still a scientific failure; hull width is None. Cap/invalid regions are preserved in conservative diagnostics, but cannot enter a scientific Gate: `scientific_cells` and `scientific_metrics` raise when partition closure is false. Thus there is no false contraction by omitting pending regions.

The future same-case B0 has the same outer semantics, domain, resolution and frozen observation-side bearing cutoff. Its partition must also close. The NOMINAL 5 km, 1/3 contraction, component caps, 0.5 km endpoint/width and topology criteria remain unchanged. Stability compares two independently closed outer partitions. A terminal possible cell means only that impossibility was not established at the chosen resolution; it does not certify every point or continuous uniqueness.

Every fixture/resolution/B0/CZ/tier uses a fresh engine, queue, cache and counters. State/depth used to generate the fixed infrastructure observations are not passed to the estimator. Results are immutable objects returned before the guarded separate evaluator receives the generating reference. These fixtures test completion and reference preservation, not scientific CZ performance.

## Structural budget and finite feasibility check

The structural analysis derives each axis' dyadic split count from domain extents and unchanged widths. For depth D, the ideal single-path partition contains D path internals, D excluded siblings and one retained terminal leaf: **2D+1** requests. This gives 37/45. The approximate 36/44 reaches a terminal child but can leave its sibling unclassified. The analytic tests distinguish these cases.

The pushed correction policy derives the three-tier ladders from that necessary minimum, the old measured worst cell cost, a declared 30/60 minute target and 32 MiB record envelope. There is no P/Q input or success-dependent additional tier. Timing targets are advisory estimates, not hard wall-clock guarantees. Storage accounting uses conservative 8192-byte record units, not a measured process-RSS guarantee; request/profile/bound/cache/validation counters and aggregate limits are hard and precharged. Every admitted run remains within the declared conservative record envelope.

Budget selection requires BOTH fixed non-development physical fixtures and BOTH B0/CZ partitions to complete and retain their generating references at a tier. Select the first qualifying tier separately per resolution. Maximum failure yields CZ_ENVELOPE_PRE_RUN_SEARCH_CLOSURE_NOT_ESTABLISHED, selected budget NOT_SELECTED, no development release. A closed partition excluding the reference or an invalid bound yields IMPLEMENTATION_INVALID. An admissible single-path budget is necessary, never a promise of physical fixture closure.

## Diagnostic depth convention

In the fixture CSV, `maximum_depth` is the maximum depth of a dequeued cell, including the active cell whose next request is denied by a hard cap. It is not the depth of every pending frontier leaf. A child may have been created at one level deeper before the cap. Independent partition reconstruction therefore also checks the maximum saved leaf depth and verifies it is no more than one level above the reported dequeued depth. Natural completion processes every leaf; a capped frontier is never described as terminal completion.
