# B1B saved joint r-v-z diagnostic

**JOINT_RVZ_PROFILE_ROUTE_NOT_SUPPORTED_BY_EXISTING_EVIDENCE**; ROUTE_DIAGNOSTIC only, PENDING_RESEARCH_LEAD_AUDIT. R4-B1 and overall remain0%. Parent 8a60aeb99297c57e574daf9b796a08114b9c1f3d was independently accepted negative; policy cb9c0cfa93bf829ca9d8f6312ad14c73283e038e was pushed before this analysis. New propagation evaluations=0, new depth scores=0, new cases/nodes/depth labels=0. Every aggregate value is exactly a selected existing B1A score with original CSV line/value and support-node provenance.

All original six configurations and16 rectangles are retained. S7 triple-frequency is primary interpretation; S8 four-line is the historical companion, not independent statistical confirmation. MIN is minimum across the saved rectangle; Q25 is the empirical inverse-CDF nearest-rank order statistic ceil(n/4), not interpolated percentile or a continuous optimizer. MEDIAN uses ceil(n/2), with odd support-node counts. Ties select ascending score/dr/dv. MIN_NO_CENTER removes only (0,0) and remains secondary. Q25/MEDIAN horizontal provenance is the selected order-statistic carrier, not an argmin; their argmin columns stay blank. This rule was frozen before execution; no quantile shopping or new rectangle.

Target rectangle R0.25 km/V0.05 m/s: {"MIN": {"n_cases": 6, "exact_depth": 6, "within_10m": 6, "rank_le3": 6, "practical_pass": 6, "worst_error_m": 0.0, "worst_true_rank": 1, "minimum_second_margin_db": 0.6524505796925171}, "Q25": {"n_cases": 6, "exact_depth": 0, "within_10m": 3, "rank_le3": 0, "practical_pass": 0, "worst_error_m": 30.0, "worst_true_rank": 13, "minimum_second_margin_db": -0.30153222811184976}, "MEDIAN": {"n_cases": 6, "exact_depth": 0, "within_10m": 1, "rank_le3": 0, "practical_pass": 0, "worst_error_m": 25.0, "worst_true_rank": 16, "minimum_second_margin_db": -0.659365582983999}, "MIN_NO_CENTER": {"n_cases": 6, "exact_depth": 0, "within_10m": 4, "rank_le3": 0, "practical_pass": 0, "worst_error_m": 45.0, "worst_true_rank": 17, "minimum_second_margin_db": -0.7028222103223405}}. D1=True, D2=True, D3=False. D3 requires Q25 within10m in at least5/6 and at least2/3 S7 cases; MEDIAN/leave-center cannot replace it. D4 all16 rectangles reported without demanding monotonicity. These route-admission diagnostics are not a B1 scientific completion or confirmation Gate.

| Case | Aggregation | H0 z-hat | H0 margin (dB) | Profiled z-hat | True rank | Profiled margin (dB) |
|---|---|---:|---:|---:|---:|---:|
| S7_Z180 | MIN | 180 | 1.11621 | 180 | 1 | 0.916792 |
| S7_Z180 | Q25 | 180 | 1.11621 | 210 | 10 | -0.212804 |
| S7_Z180 | MEDIAN | 180 | 1.11621 | 205 | 13 | -0.489271 |
| S7_Z180 | MIN_NO_CENTER | 180 | 1.11621 | 190 | 9 | -0.270465 |
| S7_Z200 | MIN | 200 | 0.982204 | 200 | 1 | 0.949073 |
| S7_Z200 | Q25 | 200 | 0.982204 | 205 | 5 | -0.0940771 |
| S7_Z200 | MEDIAN | 200 | 0.982204 | 205 | 4 | -0.169926 |
| S7_Z200 | MIN_NO_CENTER | 200 | 0.982204 | 245 | 14 | -0.659233 |
| S7_Z220 | MIN | 220 | 0.853971 | 220 | 1 | 0.738455 |
| S7_Z220 | Q25 | 220 | 0.853971 | 210 | 5 | -0.178731 |
| S7_Z220 | MEDIAN | 220 | 0.853971 | 205 | 9 | -0.31801 |
| S7_Z220 | MIN_NO_CENTER | 220 | 0.853971 | 230 | 6 | -0.137802 |
| S8_Z180 | MIN | 180 | 1.20415 | 180 | 1 | 0.898012 |
| S8_Z180 | Q25 | 180 | 1.20415 | 210 | 6 | -0.171302 |
| S8_Z180 | MEDIAN | 180 | 1.20415 | 195 | 16 | -0.659366 |
| S8_Z180 | MIN_NO_CENTER | 180 | 1.20415 | 170 | 8 | -0.301536 |
| S8_Z200 | MIN | 200 | 0.9953 | 200 | 1 | 0.973308 |
| S8_Z200 | Q25 | 200 | 0.9953 | 215 | 13 | -0.220876 |
| S8_Z200 | MEDIAN | 200 | 0.9953 | 215 | 7 | -0.120499 |
| S8_Z200 | MIN_NO_CENTER | 200 | 0.9953 | 155 | 17 | -0.702822 |
| S8_Z220 | MIN | 220 | 0.794481 | 220 | 1 | 0.652451 |
| S8_Z220 | Q25 | 220 | 0.794481 | 210 | 13 | -0.301532 |
| S8_Z220 | MEDIAN | 220 | 0.794481 | 205 | 13 | -0.445253 |
| S8_Z220 | MIN_NO_CENTER | 220 | 0.794481 | 230 | 16 | -0.391372 |


Target wrong-depth best joint margins span 0.652451-0.973308 dB. MIN true-depth argmin uses center in 6/6 configurations. At wrong depths, noncenter-node improvement over same-depth H0 occurs in 108/120 case/depth entries. Adjacent-depth argmin jumps per case: [19, 19, 19, 20, 19, 20]; range jumps [13, 14, 13, 16, 15, 15], speed jumps [15, 15, 12, 18, 14, 13]. Depth labels numerically equivalent to global minimum per case: [1, 1, 1, 1, 1, 1].

The maps demonstrate finite-node range/speed compensation where noncenter nodes reduce the same wrong-depth cost, and discrete depth-dependent branch changes. They do not establish a smooth continuous depth-range/speed tradeoff or causal decomposition; no curve is fitted. Numerical-equivalence uses the inherited 4.30873114965e-10 dB tolerance, not an invented statistical interval. Positive wrong-depth margins report sampled separation only, never continuous uniqueness. Center is in every primary rectangle and is the exact generation state: preserving its near-zero true score is expected and cannot establish an observation-derived engineering route. Q25, MEDIAN and leave-center diagnose dependence on this exact node; none represents a calibrated uncertainty likelihood or confidence interval.

The previous plug-in target had13/150 practical passes, and all16 rectangles failed universally; the smallest had7/54. B1A preserved369 nonmonotonic witnesses. H0 unique true minima6/6 and positive margins/curvature remain accepted. B1A range/speed conditioning is not robust; prior theta+/-0.5 and psi+/-1 single-factor tests are comparatively stable at those scales, without a joint4D claim. New fine conditioning fractions, B2 depth interpolation, SSP/receiver-depth errors, optimizer or propagation are not run.

25619 analysis integrity checks pass; 12 unit tests pass; separate saved-row audit reconstructs every aggregate selection, metric, compensation map,16-rectangle summary and route decision. Historical files/inputs and immutable source hashes are protected. Three main aggregations have6048 curve rows; secondary no-center has2016; all384 profiles and8064 source mappings are complete. Metadata and file hashing are the only other input reads; analysis imports no forward/depth-score module.

Current engineering status: CURRENT_CONDITIONAL_DEPTH_ENGINEERING_ROUTE_CLOSED. A positive route-admission result would only recommend a separately instructed B1R_LOCAL_JOINT_RVZ_SUPPORT_VALIDATION fresh design. A negative result closes the current conditional-depth engineering route and preserves EXACT_HORIZONTAL_DEPTH_MECHANISM as oracle/mechanism evidence. B1R/B2/A2/A3/A4/B3/B4/C/P5 are NOT_OPENED, and CZ certified acquisition remains CLOSED. Stop after analysis commit/push for independent review.

See [curves](B1B_PROFILED_DEPTH_CURVES.csv), [no-center](B1B_MIN_NO_CENTER_DEPTH_CURVES.csv), [source provenance](B1B_SCORE_SOURCE_PROVENANCE.csv), [compensation map](B1B_RV_COMPENSATION_MAP.csv), [all rectangles](B1B_RECTANGLE_DEGRADATION_SUMMARY.csv), [B1 closeout](R4_B1_CLOSEOUT.md).
