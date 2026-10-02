# R4-A1 pre-run design

Baseline: `166575043d9c884f6a9ecf62333b79fbcc7e1a46`. Design frozen on 2026-10-02 before synthetic observations or outcome inspection.

## Reuse and isolation

Reuse the formulas and read-only modal inputs of CLOSEDLOOP-1C/1D/1E: fixed platform turn, target kinematics, direct modal F matrix, per-window/per-frequency demeaning, shared-depth profiling, RC2 cmin+13.3 sigma² and RC3 tau=0.5. Implement a separate R4 harness without importing the old scripts' side-effectful module globals or running their main functions. R3 files are never overwritten.

The old modal files have 2m source-depth samples; the 5m nuisance labels are nearest-sample lookups (e.g. 155m maps to154m). Export this mapping. Source depths 180/200/220m are exact physical samples. Hold them here to isolate horizontal off-grid generalization, not develop a continuous-depth estimator or silently interpolate source modes. All four horizontal truth coordinates are jointly off-grid; range propagation is evaluated at continuous trajectory ranges, never rounded to integer km or a candidate node.

## Panel and noise design

Nine states cover the central nominal neighborhood and non-extreme lower/upper range and speed offsets, both signs of heading and bearing, and both sides of a bearing-grid cell. Two small but nonzero bearing offsets help distinguish near-grid behavior from interior-cell behavior without defining success only on favorable truths. Range47.41–55.36km, speed1.57–2.43m/s and bearing -1.18–1.17deg stay inside the inherited support. OFFGRID_TRUTH_PANEL.csv and R4_A1_CONFIG.json are the frozen inputs.

Noise levels0.02/0.05/0.1/0.2/0.4/0.8deg bracket nominal0.1deg. Thirty independent noise realizations per truth at each positive level give1620 main cases. Nine noiseless controls characterize quantization and selection without stochastic noise. Reuse each seed's standardized noise within a panel across sigma levels for paired comparisons. A three-seed nominal pilot diagnoses cost and candidate loss; it is part of the fixed main design, not a favorable subset. If structural results invalidate local refinement, document that before consuming unnecessary solver time.

## Estimator boundary and integrity

Observation generation receives the truth. Candidate estimation receives only bearing/TL observations, known sigma, immutable candidate grid, platform geometry and E0 forward models. No panel coordinates, nearest-grid IDs, true depths or truth-centered initializations cross that API. Post-hoc evaluation receives the truth and output IDs. Nearest-grid oracle and bracketing-cell vertices are diagnostic computations only. If extra oracle forward scores are evaluated, they are in a separate evaluation-only table and never injected into the accepted cloud.

Replay the on-grid R3 nominal control against stored node IDs and candidate scores at all three depths, using seed20260912 only for that identity check. Demonstrate that perturbing evaluation-only truth cannot change a fixed-observation estimate. Include angular wrapping, topology, missing-cell and nearest-oracle tie tests.

## Metrics and interpretations

Report top1 r/v relative errors and theta/psi shortest circular angular errors in degrees, plus every survivor's worst errors, parameter widths, occupied bins, 4D Manhattan components and range components. For continuous truth, the main results do not report truth-node retention/rank because that node does not exist.

Quantization floor is the coordinatewise minimum over the frozen complete grid. Excess error is actual top1 error minus each coordinate floor. The bracketing coarse cell has up to16 vertices; neighborhood retention means any vertex remains. Report nearest-oracle retention separately. Predeclare coarse-cell localization as the resolution-relative success diagnostic; do not invent the customer's joint-error norm or angle percentages. A stable tested region requires every fixed-panel/seed case at a positive sigma to meet that diagnostic. Empirical P95 from this finite design is not a95% performance guarantee for arbitrary scenes.

Distinguish grid mismatch at the RC2 cutoff, wrong coarse acoustic minima, broad/multicomponent survivor sets and range aliases. GRID_QUANTIZATION_DOMINATED / IDENTIFIABILITY_DOMINATED / MIXED are evidence interpretations, not predetermined labels. Local continuous refinement is optional and must start from estimator-produced coarse cells. If the upstream cloud deletes the relevant cells, local refinement is not presumed to repair a global search failure.

## Computation and evidence

Compute the full frozen-grid bearing predictions once; select clouds independently for every case. Acoustic scores are deterministic per truth observation and candidate, so compute direct modal responses once for the union of accepted IDs, then apply each case's own RC2 cloud and Jmin cutoff. This is caching of exact formulas, not interpolation of a coarse range/TL library. Benchmark first; chunk modal multiplication and cap BLAS threads to avoid oversubscription. Store observation provenance, accepted IDs/costs, per-panel candidate score catalogs and per-case survivor outputs so another audit can reconstruct every result.

No positive progress credit is entered merely for implementation. Execution-complete gates are submitted for independent audit; the lead's final instruction reserves management credit until that audit.
