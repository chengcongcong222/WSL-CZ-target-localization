# R4-A1-new application pre-research dynamic baseline

R4_A1_NEW_BASELINE_NOT_ESTABLISHED

Project stage: APPLICATION_PRE_RESEARCH. Hardware UNKNOWN is a nonblocking design-assumption boundary, not an actual capability claim. This new authorization supersedes the previous client-confirmation numerical STOP for this stage only. Scientific baseline333ead6 retained; chronological parent0eef9ff retained.

Primary A is the sole verdict source; B is SECONDARY_CONFIRMATION_ONLY.12off-grid truths x500realizations per anchor;121epochs;0..1200s;dt10s;MAIN2m/s,15deg turn after600s. AUX has identical displacement/velocity to MAIN and a fixed global lateral offset relative to nominal initial LoB0deg plus sampled beta, mirrored by side. No target-guided formation control.

Per realization common/HALF-diff/beta sampled independently Uniform over frozen bounds. This APPLICATION_PRE_RESEARCH_BOUNDED_SYSTEMATIC_ENSEMBLE is an assumed statistical design model, not an observed physical distribution or a claim that Gate2B certified the continuous interior. Navigation: independent epoch-wise Gaussian per node per axis; not a specific INS/GNSS temporal model. Bearing random errors independent epoch/node Gaussian. All variates/seeds frozen before MC; no resampling failures.

Observation-only estimator: positive directed-ray intersection initializer, mathematically invalid intersections excluded, fivefixed IRLS line fits with all valid points retained, then four Cartesian initial-position/velocity parameters optimized using wrapped dual-node bearing residual, soft_l1loss1.5sigma, frozen analytic Jacobian and convergence settings. No truth warm start, parameter-grid starts, acoustic model, depth or nuisance truth. r0/theta0 reference estimated MAIN position at t0; v/psi from fitted global velocity. Truth only in generator/evaluator.

Nearest-rank unconditional errors; optimizer/initializer failures count infinite error in all metrics. Per-case and pooled median/P90/P95/P99 retained. Primary Gate requires all12cases and uses worst case-level P95; no pooled substitution. Local covariance/condition/residual diagnostics are approximate and do not account for systematic/navigation model mismatch; Monte Carlo truth errors determine performance.

| Anchor | Worst range P95 % | Bearing P95 deg | Speed P95 % | Heading P95 deg | STRONG | PROJECT |
|---|---:|---:|---:|---:|---|---|
| A | 1.300624 | 0.085321 | 28.938635 | 2.563957 | False | False |
| B | 1.352744 | 0.132820 | 31.169332 | 2.850499 | False | False |

STRONG:range/speed5%,bearing0.5deg,heading2deg. PROJECT:range/speed10%,bearing1deg,heading5deg. All12cases required. Credit15only if primary PROJECT passes, pending research-lead independent audit; otherwise0. B cannot rescue primary. No application claim of five-parameter/multi-parameter<10% from these four horizontal submetrics.

Independent audit regenerates all draws, reconstructs all scenes with independent geometry, slope intersections and IRLS, cold-recomputes all12000errors/quantiles and decision.48frozen diagnostic runs independently optimized using3-point finite-difference Jacobian; this is saved-scene verification, not extra performance MC.

Single-array conditional depth engineering stays CLOSED; exact-horizontal depth mechanism SUPPORTED_ORACLE_ONLY. Depth may be reconsidered only as a separately authorized augmented-posterior stage after audit. No KRAKEN/BELLHOP/TL/CZ/SSP or depth score. Stop after execution commit/push; A2/A3/A4 are next roadmap stages only, no automatic experiment.

```json
{
  "scientific_parent_SHA": "333ead6461a7574a469495033d48fd89682aa496",
  "commit_parent_SHA": "0eef9ff7d0f40cf05400707b791db8a80e074fac",
  "design_SHA": "5632be4b2b9d9422be67ad0eca203d5ccb72eaaf",
  "project_stage": "APPLICATION_PRE_RESEARCH",
  "primary_application_scenario": "Anchor A",
  "primary_truth_cases": 12,
  "primary_runs": 6000,
  "secondary_runs": 6000,
  "primary_summary": {
    "anchor": "A",
    "n_runs": 6000,
    "failure_rate": 0.0,
    "STRONG_pass": false,
    "PROJECT_pass": false,
    "range_error_worst_case_P95": 0.013006243169259484,
    "range_error_pooled_median": 0.004856252576527137,
    "range_error_pooled_P90": 0.010093771613685507,
    "range_error_pooled_P95": 0.011525126513425875,
    "range_error_pooled_P99": 0.01399859610278917,
    "bearing_error_deg_worst_case_P95": 0.0853208192427844,
    "bearing_error_deg_pooled_median": 0.029110847867880943,
    "bearing_error_deg_pooled_P90": 0.06833214817875076,
    "bearing_error_deg_pooled_P95": 0.08069372653187937,
    "bearing_error_deg_pooled_P99": 0.10410332842861447,
    "speed_error_worst_case_P95": 0.2893863491269,
    "speed_error_pooled_median": 0.07751915235829465,
    "speed_error_pooled_P90": 0.18871373456441162,
    "speed_error_pooled_P95": 0.22775459679345614,
    "speed_error_pooled_P99": 0.30751760068116096,
    "heading_error_deg_worst_case_P95": 2.5639571210910868,
    "heading_error_deg_pooled_median": 0.3471343388450945,
    "heading_error_deg_pooled_P90": 1.1376870340353644,
    "heading_error_deg_pooled_P95": 1.496341107028101,
    "heading_error_deg_pooled_P99": 2.365832223770355
  },
  "secondary_summary": {
    "anchor": "B",
    "n_runs": 6000,
    "failure_rate": 0.0,
    "STRONG_pass": false,
    "PROJECT_pass": false,
    "range_error_worst_case_P95": 0.01352743648528172,
    "range_error_pooled_median": 0.005252796365770336,
    "range_error_pooled_P90": 0.010918744451395386,
    "range_error_pooled_P95": 0.012623136616615621,
    "range_error_pooled_P99": 0.015093923436361555,
    "bearing_error_deg_worst_case_P95": 0.13281969650172462,
    "bearing_error_deg_pooled_median": 0.044027963489263834,
    "bearing_error_deg_pooled_P90": 0.10241882952302381,
    "bearing_error_deg_pooled_P95": 0.1212926737540345,
    "bearing_error_deg_pooled_P99": 0.15768245618248986,
    "speed_error_worst_case_P95": 0.31169331631999375,
    "speed_error_pooled_median": 0.08253330576647275,
    "speed_error_pooled_P90": 0.20725436665791258,
    "speed_error_pooled_P95": 0.24808332299086266,
    "speed_error_pooled_P99": 0.34257692011391033,
    "heading_error_deg_worst_case_P95": 2.850498598912172,
    "heading_error_deg_pooled_median": 0.4630646061885072,
    "heading_error_deg_pooled_P90": 1.4039834892236218,
    "heading_error_deg_pooled_P95": 1.8275930625471724,
    "heading_error_deg_pooled_P99": 3.0900213335381213
  },
  "scientific_decision": "R4_A1_NEW_BASELINE_NOT_ESTABLISHED",
  "R4_A1_progress_percent": 0,
  "R4_overall_percent": 0,
  "progress_scope": "PRE_REGISTERED_SCIENTIFIC_COMPLETION; PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT",
  "hardware_status": "UNKNOWN; NONBLOCKING_DESIGN_ASSUMPTION",
  "depth": "NOT_OPENED",
  "next": "STOP",
  "automatic_next_experiment": false,
  "stop_after_commit_B": true,
  "new_acoustic_propagation": 0,
  "new_depth_score": 0
}
```

## Observed Gate closure

| Anchor | Metric | STRONG cases passed | PROJECT cases passed | Worst case | Worst P95 |
|---|---|---:|---:|---|---:|
| A | range_error | 12/12 | 12/12 | H12 | 1.300624 |
| A | bearing_error_deg | 12/12 | 12/12 | H02 | 0.085321 |
| A | speed_error | 0/12 | 0/12 | H09 | 28.938635 |
| A | heading_error_deg | 11/12 | 12/12 | H04 | 2.563957 |
| B | range_error | 12/12 | 12/12 | H11 | 1.352744 |
| B | bearing_error_deg | 12/12 | 12/12 | H02 | 0.132820 |
| B | speed_error | 0/12 | 0/12 | H09 | 31.169332 |
| B | heading_error_deg | 10/12 | 12/12 | H04 | 2.850499 |

Both anchors have0initializer/solver failures in6000realizations, but speed P95 fails10%inall12cases. Primary A range/bearing meet STRONGinall12; heading meets PROJECTinall12 andSTRONGin11. This is a performance Gate failure under the frozen observation-only estimator and application error model, not evidence of a runtime failure, target nonexistence, physical nonidentifiability or impossibility of every algorithm. No post-result experiment, parameter change or threshold relaxation was undertaken.

The different worst-case P95 values may come from different truth cases; the frozen Gate is four submetric requirements, not a single joint95%event probability guarantee. The proposed static input specification alone is insufficient to claim the frozen dynamic four-parameter PROJECTbaseline.
