# Auxiliary-node geometry Gate0

AUXILIARY_NODE_GEOMETRY_CONDITIONALLY_PRACTICAL

Geometry-only architecture diagnostic; no five-parameter or depth claim.

Requested cells: 1920; feasible: 240; impossible: 1680. Each feasible cell has 20000 Gaussian trials. Every physical branch and mirror retained; summaries take worst branch/mirror.

## Physical feasibility

With MAIN=(0,0), target=(R,0), auxiliary-target distance d obeys d^2-2R cos(alpha)d+R^2-B^2=0. Feasible iff B>=R sin(alpha), d>0. Both roots and reflected placements are frozen. Impossible baseline/crossing combinations are not simulated. All 20-90 degree cells are impossible at the registered far ranges/baselines. No geometry added after freeze.

## Sigma 0.1 degree

| Baseline km | Best feasible-subset angle | Feasible ranges km | Worst P95 % | Best full 50/55/60 km P95 % |
|---|---|---|---|---|
| 0.5 | NONE | NONE | N/A | N/A |
| 1 | NONE | NONE | N/A | N/A |
| 2 | 2 | 50;55 | 14.133082 | N/A |
| 5 | 5 | 50;55 | 5.639378 | 14.646014 |
| 10 | 10 | 50;55 | 2.856294 | 6.009670 |

Full per-range regions (worst both roots/mirrors):
```json
[
  {
    "range_km": 50,
    "baseline_km": 2,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.14133082214185946,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 55,
    "baseline_km": 2,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.1405761807895145,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 50,
    "baseline_km": 5,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.1464601356491147,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 55,
    "baseline_km": 5,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.145849331301152,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 60,
    "baseline_km": 5,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.1452274482981257,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 50,
    "baseline_km": 5,
    "alpha_deg": 5,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.05639377645111139,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": true
  },
  {
    "range_km": 55,
    "baseline_km": 5,
    "alpha_deg": 5,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.055771977414679375,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": true
  },
  {
    "range_km": 50,
    "baseline_km": 10,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.15407941748240758,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 55,
    "baseline_km": 10,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.15268917837019907,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 60,
    "baseline_km": 10,
    "alpha_deg": 2,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.15166525450310667,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": false
  },
  {
    "range_km": 50,
    "baseline_km": 10,
    "alpha_deg": 5,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.060096699319185376,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": true
  },
  {
    "range_km": 55,
    "baseline_km": 10,
    "alpha_deg": 5,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.05945022046439228,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": true
  },
  {
    "range_km": 60,
    "baseline_km": 10,
    "alpha_deg": 5,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.058973586430809986,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": false,
    "le10pct": true
  },
  {
    "range_km": 50,
    "baseline_km": 10,
    "alpha_deg": 10,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.02856294121782497,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": true,
    "le10pct": true
  },
  {
    "range_km": 55,
    "baseline_km": 10,
    "alpha_deg": 10,
    "sigma_deg": 0.1,
    "n_feasible": 4,
    "worst_p95_relative_range": 0.027966740222428982,
    "worst_failure_rate": 0.0,
    "status": "FEASIBLE",
    "le5pct": true,
    "le10pct": true
  }
]
```

Full-range <=5% regions: []
Full-range <=10% regions: [{"baseline_km": 10, "alpha_deg": 5, "sigma_deg": 0.1, "n_feasible_cells": 12, "feasible_ranges_km": "50;55;60", "all_three_ranges_feasible": true, "worst_feasible_p95_relative_range": 0.060096699319185376, "worst_failure_rate": 0.0, "all_feasible_cells_strong": false, "all_feasible_cells_usable": false, "all_feasible_cells_project_level": true, "full_range_usable": false, "full_range_project_level": true}]

## Statistical and engineering scope

Independent zero-mean Gaussian errors between sensors; shared frozen standard-normal draws across cells, mirror sign reversal for exact reflection checks. Cells are not independent trials. Directed target bearings, correct target association and common time are assumed; sensor positions are known exactly. Median/P90/P95/P99 are unconditional nearest-rank order statistics. Parallel, nonfinite and behind-sensor estimates receive infinite errors; finite ill-conditioned cases are retained. Ill-conditioning proxy: 1/abs(sin(observed crossing))>1000. Rates are separate and may overlap.
Range is norm(position-MAIN); cross-range is absolute perpendicular position error to the true MAIN LoB; 2D position error is Euclidean target error. P95 empirical rank uncertainty uses a frozen approximate binomial-normal 95% order-statistic interval; Wilson intervals report fraction <=5%. Primary Gate is empirical P95, as frozen; these intervals do not give selected-grid simultaneous guarantees.
H/FIM inverse is a local linearized lower-bound/approximation, not actual estimator guarantee. At near-parallel geometry nonlinear ratios can have long tails; Monte Carlo quantiles, including failures, govern the Gate. Known positions and unbiased bearings exclude positioning, attitude, shared bias, clock offset, latency, detection and signal-association error.
Baseline 0.5/1 km has no feasible registered angle >=2 deg; N/A is no tested triangle, not proof those baselines are universally impossible. The frozen grid omits their smaller physically possible angles and does not cover every placement. No post-result angle additions. Similarly an infeasible 60 km cell cannot be replaced by a nearer target when claiming full-range performance.
The 5 km practical ceiling is a planning assumption, not verified deployment capability. Hardware, independent directed-bearing accuracy, front/back resolution, self-localization/attitude, time/target association, communication, maneuver limits and acoustic detection remain future gates. An omnidirectional hydrophone alone cannot provide the assumed bearing; no HLA/VLA type is preselected.

## Route boundaries

Two nonparallel directed bearings remove instantaneous horizontal scale ambiguity with known separated nodes. A temporal position sequence could support velocity/heading, but a single epoch cannot observe v/psi, and no tracker or motion accuracy is evaluated. Depth is absent from this model.
Independent horizontal acquisition may justify future conditional-depth validation under a new architecture. Few-percent range is not a validated B1 depth-safe tolerance: B1A/B1B supplied no universally robust envelope. Current single-array depth stays closed, exact-horizontal mechanism retained as oracle evidence.
Outcome A requires full-range <=5% at <=5 km, worst all roots/mirrors. Outcome B admits only conditional consideration (partial-range 5% or full-range 10%); architecture acceptance and hardware scale are research-lead decisions. No automatic R4-A1-NEW or depth restart.

No acoustic propagation/depth scores/global optimizer/SSP/5D estimator/tracker. R4=0%; A2/A3/A4/B2/B3/B4/C paused/not opened; P5 not opened. Commit/push/verify and stop.

## Near-parallel boundary and implementation audit

At sigma0.1 deg and alpha2 deg, every feasible registered region fails even the10% range band: worst-root P95 spans14.06%-15.41%, despite zero finite-forward-intersection failure rate. This is accuracy/conditioning degradation, not execution failure. The maximum sigma0.5-deg behind/nonfinite/parallel combined failure rate is0.215%; these failures remain in unconditional quantiles. No untested continuous crossing-angle transition is inferred.

The independent audit required a documented infrastructure correction to two SHA256 identifiers in the frozen audit script. Original source remains unchanged; corrected source changes only the hash function name and policy-hash field. Both failed audit attempts and hashes are recorded in AUX_AUDIT_IMPLEMENTATION_CORRECTION.json/.md. No primary execution was rerun; no numerical policy, geometry, draws, metric or Gate changed. Independent slope/covariance verification passed33611 checks; all10 unit tests passed. Use r4_aux_gate0_audit_corrected.py --verify.
