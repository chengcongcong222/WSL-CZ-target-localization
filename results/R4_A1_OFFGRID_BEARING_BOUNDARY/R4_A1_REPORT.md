# R4-A1 structural stop report

Decision: **R4_A1_BLOCKED_BY_COARSE_RC2_GRID_SELECTION_AND_OFFGRID_MODAL_MISMATCH**. R4 progress remains **0%**. R3 stays frozen.

## Main finding

This frozen coarse pipeline does not generalize successfully to the9 jointly off-grid truths in the executed nominal pilot. Top1 coarse-cell localization succeeds in0/27 cases. RC2 already deletes every bracketing-cell vertex in25/27 pilot cases; only1/27 final sets retains any such vertex. These are resolution-relative diagnostics, not the customer's as-yet-undefined joint error norm.

The cheaper full RC2-only axis confirms a structural selection defect: at sigma0.02 and0.05deg, every one of270 cases loses the truth's bracketing coarse cell; the9 noiseless cases also lose it. Larger noise relaxes the inherited selection cloud and raises neighborhood retention. This inversion does not mean worse bearing measurements improve physical observability. The rule represents bearing noise but not coarse bearing-grid mismatch.

On-grid replay reproduces all385 R3 RC2 IDs and TRIPLE scores/survivors at180/200/220m. Thus this is an R4 off-grid counterexample to extending the frozen tested control, not a reinterpretation or overwrite of R3 numerical evidence.

## Frozen design and executed scope

Baseline `166575043d9c884f6a9ecf62333b79fbcc7e1a46`. Nine jointly off-grid horizontal states:

```text
panel_id  r_km  theta_deg  v_mps  psi_deg  z_true_m                                     role
     P01 50.37       0.17   2.07     5.35       200                     central_joint_offset
     P02 49.63      -0.18   1.93     4.65       180                  central_opposite_offset
     P03 50.71       0.32   2.13     6.42       220          central_other_bearing_cell_side
     P04 47.41      -0.73   1.57    -6.35       180         lower_range_and_speed_nonextreme
     P05 48.68       0.68   2.37     8.62       220                 lower_range_higher_speed
     P06 55.36       1.17   2.43     9.38       200         upper_range_and_speed_nonextreme
     P07 53.72      -1.18   1.73    -8.42       220 upper_range_lower_speed_opposite_heading
     P08 51.24       0.06   1.87     2.28       180             small_bearing_offset_control
     P09 54.61      -0.07   2.23    -3.64       200    small_negative_bearing_offset_control
```

TRIPLE201+235+283Hz, E0 matched modal inputs,2m/s platform,15deg turn at600s,1200s observation,121 bearing/TL samples. Source depths180/200/220m are profiled nuisances, held on physical samples to isolate horizontal gridding. The old5m nuisance labels map to2m modal samples; see NUISANCE_DEPTH_MAPPING.csv. No depth interpolation, depth-estimator development, nav/SSP sweep, joint corner, real HLA extraction or P5.

Pre-run configuration/panel/design hashes are retained in DESIGN_FREEZE_MANIFEST.json. Planned positive levels are0.02/0.05/0.1/0.2/0.4/0.8deg with30 seeds per panel; sigma0 has one noiseless case per panel. Executed:1629 RC2-only cases,27 nominal end-to-end pilot cases with seeds410001–410003, and3 on-grid identity controls. The main1620-case end-to-end matrix was **not run**. PILOT_ACTION_DECISION.md records the structural-stop reason. Acoustic evaluation used only the union of2165 pilot-accepted nodes plus separate evaluation-only oracle scores, not a truth-centered estimator search.

## Leakage and forward integrity

Continuous truth coordinates generate bearing and direct modal propagation along exact continuous ranges. The estimator sees only observations, known sigma and the fixed candidate grid/catalog. Neither truth coordinates, nearest-oracle IDs nor true source depth are estimator inputs. Evaluation-only oracles are stored separately and never injected into acceptance or initialization. OFFGRID_INTEGRITY_AUDIT.csv records truth-feature nonidentity, estimator signatures, fixed-observation evaluation invariance and R3 replay.

The inherited per-window/per-frequency relative-TL demeaning, shared z profile, RC2 cmin+13.3sigma² and RC3 J<=Jmin+0.5 remain fixed. Source selection does not choose new frequencies, turn or support from the truth.

## Quantization versus estimation

Nearest-grid oracle: coordinatewise nearest complete-grid node, with circular angular differences and all ties. Quantization floor is its per-axis minimum error. Bracketing cell is the Cartesian product of floor/ceiling grid coordinates, up to16 vertices; retention of any vertex and retention of nearest oracle are separate. Top1 excess error is actual error minus coordinate floor. Angle errors are degrees, never percentages.

```text
panel_id  floor_rel_r  floor_abs_theta_deg  floor_rel_v  floor_abs_psi_deg
     P01     0.007346                 0.17     0.033816               0.35
     P02     0.007455                 0.18     0.036269               0.35
     P03     0.005719                 0.18     0.032864               0.42
     P04     0.008648                 0.23     0.019108               0.35
     P05     0.006574                 0.18     0.012658               0.38
     P06     0.006503                 0.17     0.012346               0.38
     P07     0.005212                 0.18     0.040462               0.42
     P08     0.004684                 0.06     0.037433               0.28
     P09     0.007142                 0.07     0.013453               0.36
```

The nearest-grid modal score is higher than the selected score in all27 pilot cases, with positive score differences stored in FAILURE_MODE_DIAGNOSTIC.csv. Consequently, merely reinserting that oracle node would not establish correct top1 localization. Coarse source propagation mismatch/alias selection and upstream RC2 loss both appear. Interpretation: **MIXED**. This does not prove intrinsic continuous physical non-identifiability; a reliable continuous-horizontal search has not been validated. Local refinement at a selected coarse minimum is deferred rather than assumed to recover the deleted support.

## Nominal pilot output errors

The following are empirical linear-interpolation P50/P95 for27 cases with3 seeds, **only at0.1deg**. Relative r/v entries are fractions; angular entries are degrees. They are not a multi-sigma engineering accuracy curve or a95% ocean performance guarantee.

```text
       metric  top1_P50  top1_P95  survivor_P50  survivor_P95
        rel_r  0.075596  0.160093      0.175975      0.208946
abs_theta_deg  0.180000  0.230000      0.180000      0.230000
        rel_v  0.069519  0.273885      0.234568      0.734104
  abs_psi_deg  5.650000  8.650000      6.580000      8.650000
```

All individual top1 errors, coordinate floors/excesses, survivor worst errors and spans are in CASE_LEVEL_RESULTS.csv / SURVIVOR_SET_QUALITY.csv. 4D Manhattan components and range-projection components are reported per case. z_star is a nuisance diagnostic, not a depth estimate. truth-node retained/rank are not off-grid performance metrics.

## RC2 bearing-noise axis

```text
                scope  sigma_deg  n_cases  nearest_grid_oracle_retention_fraction  truth_cell_vertex_retention_fraction  n_RC2_P50  n_RC2_P90  n_RC2_P95  n_RC2_worst_observed
NOISELESS_RC2_CONTROL       0.00        9                                0.000000                              0.000000        1.0        1.0       1.00                   1.0
   RC2_ONLY_FULL_AXIS       0.02      270                                0.000000                              0.000000        3.0        9.0      14.00                  17.0
   RC2_ONLY_FULL_AXIS       0.05      270                                0.000000                              0.000000       15.0       66.0      76.55                 147.0
   RC2_ONLY_FULL_AXIS       0.10      270                                0.003704                              0.088889      103.0      373.7     418.30                 490.0
   RC2_ONLY_FULL_AXIS       0.20      270                                0.148148                              0.211111      632.0      943.3     970.00                1228.0
   RC2_ONLY_FULL_AXIS       0.40      270                                0.355556                              0.533333     1886.0     2491.6    2672.20                2936.0
   RC2_ONLY_FULL_AXIS       0.80      270                                0.859259                              0.985185     4798.0     5521.0    5618.00                5826.0
```

This table covers only RC2 candidate selection. Do not read it as final output accuracy or a stability region. Better bearing precision narrows the tolerated residual cloud enough to exclude coarse representations of the continuous truth. First observed selection failure is already at the lowest positive tested sigma0.02deg and exists in noiseless controls; no sensor-degradation crossover is established.

## Gate status and next route

- A1-1: off-grid observation/estimator separation validated in execution; independent audit pending.
- A1-2: quantization floors and pilot excess/envelope diagnosed; full-stage discrimination incomplete, interpretation MIXED.
- A1-3: blocked; the full end-to-end bearing-to-output statistical boundary was not executed.
- A1-4: NO_STABLE_REGION_ESTABLISHED in the executed evidence; unscored sigma levels are not inferred failures or successes.

Recommendation: open an **A1-FIX Gate**, not A2, to validate truth-independent treatment of coarse bearing-grid mismatch before hard pruning and a continuous-horizontal scoring/search strategy. Keep the existing truth panel/seeds as an unchanged regression challenge; do not tune thresholds or initialize around truth. Only after that Gate should the full planned end-to-end matrix resume. Management credit remains0%; code and a pipeline replay do not equal A1 completion.

## Reproduction and audit artifacts

`python r4_a1_offgrid_bearing.py pilot` reproduces the structural pilot selection.
`python r4_a1_offgrid_bearing.py run` executes the documented reduced structural-stop path.
`python r4_a1_offgrid_bearing.py report` regenerates summaries without propagation.

PILOT_OBSERVATIONS.npz stores exact noisy bearing and perfect relative-TL observations; RC2_ACCEPTED_CLOUDS.npz stores accepted IDs/costs in CSR order; CANDIDATE_SCORE_CATALOG.csv.gz stores deterministic per-panel scores on the pilot union; SURVIVOR_NODES.csv.gz stores each case's final nodes. These artifacts reconstruct case-level minima, IDs and metrics. GRID_STATE_TABLE.csv.gz freezes ID-to-coordinate correspondence. INPUT_SOURCE_MANIFEST.csv identifies read-only modal/control inputs.

`python r4_a1_result_audit.py` independently reconstructs the pilot bearing geometry and circular residuals, accepted node IDs, final scored outputs, quantization floors and error envelopes from saved artifacts. It also verifies source/design hashes and recomputes a small modal-score sample. See [LOCAL_VALIDATION.md](LOCAL_VALIDATION.md) for scope and commands. This local consistency validation does not award scientific Gate credit.

The R4 evidence directory preserves exact bytes through its local .gitattributes so design hashes survive Windows checkout. The source manifest records both execution-byte SHA256 and LF-normalized SHA256 for the legacy CSV source; modal hashes are exact binary hashes. Input-source hashes were captured after execution; only the design manifest is a pre-run freeze.
