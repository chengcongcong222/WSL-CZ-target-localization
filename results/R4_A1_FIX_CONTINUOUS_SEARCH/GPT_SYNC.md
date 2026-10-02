# R4-A1-FIX continuous search

**A1_FIX_BLOCKED_BY_FINITE_CONTINUOUS_SEARCH_FAILURE; R4 progress 0%.** This is a repair Gate, not the full sensor-accuracy experiment. R3 and structural-stop evidence remain frozen. A1-1 was independently accepted at baseline 4633ff0.

## Method and frozen scope

See METHOD_FREEZE.md and R4_A1_FIX_CONFIG.json. Starts are generated from bearing observations, the full bounded state space and independent Sobol sequences. Two global acoustic DE islands search the continuous bearing likelihood region; up to sixteen diverse spline-refined basins and up to four exact-polished basins are preserved. Collapsed populations can produce fewer separated starts; actual counts are exported. Exact modal scores determine final outputs. Neither truth coordinates nor truth-cell vertices enter either estimator API. Continuous source depth is not estimated: the shared inherited nuisance profile is unchanged.

The six holdout states follow seed 2026100203, fixed interior ranges, rejection of on-grid coordinates and the frozen generation rule. Their panel and method hash were frozen before any holdout observation. Development: 9 truths x (noiseless + one nominal seed). Holdout: 6 truths x (noiseless + two nominal seeds). A V1 theta-offset initialization defect was corrected and a fresh holdout frozen; V1 evidence remains in INITIAL_IMPLEMENTATION. Search budgets and scientific thresholds were not changed to improve panel success. No full multi-sigma MC is run.

## Continuous self-match and legacy controls

| panel_id | bearing_cost_noiseless | J_TRIPLE | z_star_label_m | bearing_tolerance | acoustic_tolerance_db | pass_self_match |
| --- | --- | --- | --- | --- | --- | --- |
| P01 | 0 | 9.810473e-15 | 200 | 1e-20 | 1e-10 | True |
| P02 | 0 | 3.08075e-14 | 180 | 1e-20 | 1e-10 | True |
| P03 | 0 | 1.340504e-14 | 220 | 1e-20 | 1e-10 | True |
| P04 | 0 | 1.111331e-14 | 180 | 1e-20 | 1e-10 | True |
| P05 | 0 | 1.136085e-14 | 220 | 1e-20 | 1e-10 | True |
| P06 | 0 | 1.067933e-14 | 200 | 1e-20 | 1e-10 | True |
| P07 | 0 | 4.785249e-14 | 220 | 1e-20 | 1e-10 | True |
| P08 | 0 | 1.484275e-14 | 180 | 1e-20 | 1e-10 | True |
| P09 | 0 | 2.06434e-14 | 200 | 1e-20 | 1e-10 | True |

Tolerance: bearing cost <1e-20 rad squared; exact profiled acoustic J <1e-10 dB. Truth evaluation is an oracle forward-integrity control only. R3's exact 385-node IDs and all TRIPLE scores/survivors at three depths are checked against frozen evidence. All 27 legacy A1 failure cases are replayed from frozen observations/catalogs; legacy scoring remains BASELINE_LEGACY.

## Regression results

| case_id | rc2_cmin | Jmin | top1_rel_r | top1_abs_theta_deg | top1_rel_v | top1_abs_psi_deg | recovered_matched_acoustic_basin |
| --- | --- | --- | --- | --- | --- | --- | --- |
| P01_s0_seed0 | 0 | 1.644818e-11 | 1.410647e-16 | 0 | 6.693519e-14 | 3.694822e-13 | True |
| P01_s0.1_seed410001 | 0.000287994 | 7.830926e-12 | 8.46388e-16 | 4.428102e-11 | 2.076063e-12 | 1.323144e-09 | True |
| P02_s0_seed0 | 0 | 1.502024e-10 | 2.577024e-15 | 0 | 6.43124e-14 | 3.694822e-13 | True |
| P02_s0.1_seed410001 | 0.0003465852 | 0.9565184 | 0.03424413 | 0.04614817 | 0.003745017 | 1.16643 | False |
| P03_s0_seed0 | 0 | 1.296635e-11 | 2.802377e-16 | 0 | 6.796858e-14 | 3.979039e-13 | True |
| P03_s0.1_seed410001 | 0.0002587273 | 6.584841e-12 | 0 | 4.263256e-13 | 8.256306e-14 | 4.16378e-11 | True |
| P04_s0_seed0 | 0 | 2.577323e-11 | 2.997438e-16 | 0 | 2.079016e-14 | 7.105427e-13 | True |
| P04_s0.1_seed410001 | 0.0003127648 | 2.863739 | 0.1057827 | 0.07178723 | 0.06077171 | 6.414827 | False |
| P05_s0_seed0 | 0 | 1.043275e-10 | 1.605581e-15 | 0 | 1.80259e-13 | 1.392664e-12 | True |
| P05_s0.1_seed410001 | 0.0003345002 | 2.468742 | 0.03290067 | 0.04811987 | 0.005642059 | 1.140028 | False |
| P06_s0_seed0 | 0 | 1.228342e-10 | 5.133979e-15 | 0 | 2.288065e-13 | 1.875833e-12 | True |
| P06_s0.1_seed410001 | 0.0003507874 | 3.164496 | 0.1702863 | 0.04369743 | 0.166636 | 1.700857 | False |
| P07_s0_seed0 | 0 | 3.252412e-10 | 4.364838e-15 | 0 | 5.236659e-14 | 6.82121e-13 | True |
| P07_s0.1_seed410001 | 0.0003411763 | 2.82687e-11 | 1.058143e-15 | 2.759748e-11 | 1.880063e-12 | 8.62002e-10 | True |
| P08_s0_seed0 | 0 | 1.43693e-11 | 0 | 0 | 9.131139e-14 | 1.98952e-13 | True |
| P08_s0.1_seed410001 | 0.0002362221 | 5.176917e-12 | 1.386695e-16 | 4.831691e-12 | 6.839449e-14 | 1.019771e-10 | True |
| P09_s0_seed0 | 0 | 1.162015e-10 | 2.211907e-15 | 0 | 3.405348e-14 | 1.136868e-13 | True |
| P09_s0.1_seed410001 | 0.0004001031 | 3.420892 | 0.08443922 | 0.02858491 | 0.1006882 | 1.971503 | False |

## Holdout results

| case_id | rc2_cmin | Jmin | top1_rel_r | top1_abs_theta_deg | top1_rel_v | top1_abs_psi_deg | recovered_matched_acoustic_basin |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Q01_s0_seed0 | 0 | 2.98582e-10 | 1.952721e-15 | 0 | 5.253065e-14 | 2.557954e-13 | True |
| Q01_s0.1_seed420001 | 0.0003565861 | 1.003801 | 0.1798214 | 0.005395629 | 0.001295376 | 1.054368 | False |
| Q01_s0.1_seed420002 | 0.0004529755 | 0.8897175 | 0.03624825 | 0.03761986 | 0.002575963 | 1.925037 | False |
| Q02_s0_seed0 | 0 | 6.909689e-11 | 3.847411e-16 | 0 | 9.580162e-14 | 7.673862e-13 | True |
| Q02_s0.1_seed420001 | 0.0003461412 | 2.471263e-11 | 2.564941e-16 | 1.020339e-11 | 1.169477e-12 | 4.685319e-10 | True |
| Q02_s0.1_seed420002 | 0.000450388 | 2.780646e-11 | 1.28247e-16 | 1.460876e-11 | 2.401543e-12 | 9.550263e-10 | True |
| Q03_s0_seed0 | 0 | 1.357924e-11 | 2.609819e-15 | 0 | 3.651219e-14 | 8.526513e-14 | True |
| Q03_s0.1_seed420001 | 0.0003983493 | 7.01898e-12 | 2.747178e-16 | 1.432454e-11 | 2.637404e-13 | 2.399645e-10 | True |
| Q03_s0.1_seed420002 | 0.0003700999 | 7.993255e-12 | 4.120767e-16 | 2.950173e-11 | 3.645659e-13 | 2.914931e-10 | True |
| Q04_s0_seed0 | 0 | 2.106926e-11 | 1.304194e-15 | 0 | 1.544664e-14 | 1.705303e-13 | True |
| Q04_s0.1_seed420001 | 0.0003207714 | 0.001682408 | 4.220414e-07 | 0.01952176 | 0.0008443693 | 0.4903585 | False |
| Q04_s0.1_seed420002 | 0.0003524499 | 1.649136e-11 | 0 | 4.14957e-12 | 1.596984e-13 | 8.927259e-11 | True |
| Q05_s0_seed0 | 0 | 1.210281e-11 | 1.414566e-16 | 0 | 1.070954e-14 | 5.684342e-14 | True |
| Q05_s0.1_seed420001 | 0.000376076 | 1.021326 | 0.03446841 | 0.03962856 | 0.01243396 | 3.224984 | False |
| Q05_s0.1_seed420002 | 0.0003017158 | 7.590827e-12 | 0 | 1.455192e-11 | 2.812156e-13 | 3.952039e-10 | True |
| Q06_s0_seed0 | 0 | 9.497016e-11 | 1.970375e-15 | 0 | 5.106133e-15 | 2.842171e-14 | True |
| Q06_s0.1_seed420001 | 0.0003279162 | 1.285725e-11 | 1.313583e-16 | 1.443823e-11 | 5.354145e-13 | 4.332321e-10 | True |
| Q06_s0.1_seed420002 | 0.0003918848 | 1.190822e-11 | 1.313583e-16 | 1.952571e-11 | 6.675904e-13 | 5.302638e-10 | True |

## Interpretation and limits

Gate flags: `{"CONTINUOUS_FORWARD_SELF_MATCH_VALIDATED": true, "TRUTH_INDEPENDENT_RC2_CONTINUATION_VALIDATED": true, "CONTINUOUS_MULTIFREQUENCY_SEARCH_VALIDATED": false, "LEGACY_R3_CONTROL_PRESERVED": true, "FROZEN_OFFGRID_REGRESSION_IMPROVED": false, "HOLDOUT_OFFGRID_STRUCTURAL_GENERALIZATION_CONFIRMED": false}`.

Recovered exact acoustic basins: regression noiseless 9/9, nominal 4/9; holdout noiseless 6/6, nominal 8/12. Continuous-truth bearing-likelihood coverage is 36/36. These fixed small panels are structural tests, not Monte Carlo accuracy boundaries.

Noiseless bearing solutions use independently full-rank Cartesian geometry and bounded nonlinear multistart refinement; distance/region extent and all convergence starts are exported. For noisy cases the cutoff is unchanged, but acts on exact continuous states. RC2 sampling and acoustic global search remain finite and do not certify exhaustive modes or a confidence envelope. Survivor envelopes describe exported candidates only. DE reaching its fixed generation budget is explicitly recorded, not counted as optimizer convergence. Local convergence fractions are in raw logs.

`n_sampled_acoustic_modes` in raw case data is a finite separated-candidate proxy at the frozen deduplication resolution, not a certified count of continuous local or physical modes. MODE_INVENTORY.csv distinguishes the certified unique noiseless bearing solution from non-exhaustive noisy hypothesis clusters. DE `nfev` denotes vectorized objective batches, not individual state evaluations.

An exact synthetic truth has near-zero acoustic score, so a positive optimized minimum, while truth remains inside the bearing region, is evidence of finite-search failure. It does not establish physical aliasing. No two distinct continuous states with indistinguishable observations have been demonstrated. A recovered acoustic basin requires exact J <0.001 dB, not just an improvement in top1 error.

CASE_LEVEL_RESULTS.csv reports angular circular errors, relative r/v errors, coordinatewise best-candidate errors, candidate envelopes, bearing/acoustic costs, sampled basin counts and convergence. Coarse quantization floors are reference only. HOLDOUT_RESULTS.csv is structurally confirmatory for this fixed interior design, not an engineering bearing tolerance or general ocean guarantee.

## Reproduce

`python r4_a1_fix_continuous.py run` (checkpoints allow deterministic resume); `python r4_a1_fix_continuous.py report` regenerates summaries. Forward interpolation is independently checked at fixed random ranges and final candidates use exact scoring. API signatures, frozen hashes and legacy isolation are locally checked; independent scientific audit remains pending.

`python r4_a1_fix_audit.py` reconstructs all final acoustic scores, bearing costs, errors and candidate envelopes from separately regenerated observations, verifies legacy challenge bytes, exports optimizer starts/convergence and checks exact-alias diagnostics. `python r4_a1_fix_figures.py` creates scientific figures from saved data only. See LOCAL_VALIDATION.md for final check counts and limitations.

Recommended next Gate if blocked: strengthen truth-independent acoustic basin coverage and verify global-search convergence, for example with physically resolved hierarchical range/profile searches or independent global solvers. Preserve V2 regression/holdout failures as development evidence and preregister a new confirmation panel after any method change. Do not inflate the bearing threshold or begin full A1 statistics before repair validation.
