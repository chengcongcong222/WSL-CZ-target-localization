# H3-G0 finite vertical mechanism result

Final scoped decision: H3_G0_VERTICAL_DEPTH_INFORMATION_PRESENT_CONDITIONALLY.
High-level class: H3_G0_VERTICAL_DEPTH_INFORMATION_PRESENT_CONDITIONALLY.
Design SHA: cfab8f500e12af45e576e46fb9312a2530cae2ab. Execution SHA is commit B; verified remote SHA is reported after commit.

Independent saved-evidence reconstruction: 1272 PASS / 0 FAIL. Maximum cold pressure relative difference 0; Jacobian 1.50230303e-16; effective depth QR/SVD 2.7598781e-12.

Registered input/gradient arrays: 36. Information rows: 720. Finite F label diagnostics: 4500.

## Calibration and horizontal profiling

|Noise|Resource|Calibration|Profiled Iz min/max (m^-2)|Oracle-horizontal Iz min|Ranks|Worst local scale (m)|
|---|---|---|---|---|---|---|
|RELATIVE|B0|C0|2755.147 / 307473.3|2800.274|5|0.01905143|
|RELATIVE|B0|C1a|2687.906 / 298999.1|2787.372|5|0.01928826|
|RELATIVE|B0|C1b|1893.569 / 201117.6|1999.746|5|0.0229805|
|RELATIVE|B0|C2|0 / 0|0|0|inf|
|RELATIVE|B0|B3_KNOWN_SOURCE_AND_CALIBRATION|25199.58 / 754831.3|25920.95|5|0.00629946|
|RELATIVE|B1|C0|11813.8 / 892757|11893.18|5|0.009200368|
|RELATIVE|B1|C1a|11699.69 / 867756.3|11827.93|5|0.009245125|
|RELATIVE|B1|C1b|7722.412 / 505461.6|7869.754|5|0.01137951|
|RELATIVE|B1|C2|0 / 0|0|0|inf|
|RELATIVE|B1|B3_KNOWN_SOURCE_AND_CALIBRATION|43926.76 / 1537506|44807.64|5|0.004771286|
|RELATIVE|B2|C0|2861.327 / 330247.8|2910.905|5|0.0186946|
|RELATIVE|B2|C1a|2794.448 / 321365.9|2897.445|5|0.01891699|
|RELATIVE|B2|C1b|1975.716 / 215940.1|2087.543|5|0.02249768|
|RELATIVE|B2|C2|0 / 0|0|0|inf|
|RELATIVE|B2|B3_KNOWN_SOURCE_AND_CALIBRATION|37778.46 / 1123718|38950.28|5|0.005144911|
|ABSOLUTE_FLOOR|B0|C0|5448.498 / 15982.29|5532.084|5|0.01354758|
|ABSOLUTE_FLOOR|B0|C1a|5244.821 / 15640.01|5476.736|5|0.01380812|
|ABSOLUTE_FLOOR|B0|C1b|3643.011 / 10310.67|3928.681|5|0.01656799|
|ABSOLUTE_FLOOR|B0|C2|0 / 0|0|0|inf|
|ABSOLUTE_FLOOR|B0|B3_KNOWN_SOURCE_AND_CALIBRATION|55031.68 / 95853.88|57907.5|5|0.004262787|
|ABSOLUTE_FLOOR|B1|C0|27390.49 / 77195.36|27605.22|5|0.00604227|
|ABSOLUTE_FLOOR|B1|C1a|26488.43 / 70603.72|26870.13|5|0.006144293|
|ABSOLUTE_FLOOR|B1|C1b|16951.16 / 40815.7|17431.79|5|0.007680692|
|ABSOLUTE_FLOOR|B1|C2|0 / 0|0|0|inf|
|ABSOLUTE_FLOOR|B1|B3_KNOWN_SOURCE_AND_CALIBRATION|98957.63 / 183459.1|102928.7|5|0.003178889|
|ABSOLUTE_FLOOR|B2|C0|5655.805 / 17310.83|5743.941|5|0.01329697|
|ABSOLUTE_FLOOR|B2|C1a|5447.726 / 16889.13|5682.286|5|0.01354854|
|ABSOLUTE_FLOOR|B2|C1b|3788.557 / 11316.89|4080.772|5|0.01624662|
|ABSOLUTE_FLOOR|B2|C2|0 / 0|0|0|inf|
|ABSOLUTE_FLOOR|B2|B3_KNOWN_SOURCE_AND_CALIBRATION|82432.11 / 143560.4|87035.37|5|0.003482987|

Values above are deterministic nominal local information diagnostics; if the numerical Gate fails, they must not be interpreted as established depth information or estimator precision. Any INF depth scale denotes a null/unidentifiable direction; it is never reported as zero variance. Full singular values, weakest direction and nullspace loading are retained per row.

C1b unknown per-element/per-frequency gains are fixed across all three snapshots, and source is arbitrary at each frequency/snapshot. C1a fixes per-element gain across every frequency. C2 spans the full response and gives zero information. Arbitrary free modal gains analytically absorb depth: all modal source-200 coefficients are nonzero in the registered caches. Constrained group gain is NOT_ADMITTED: no physical constraint or group count has been invented. C1e environment/pose uncertainty is NOT_EVALUATED.

## Vertical gain and numerical closure

|Scene|Noise|B1/B2 profiled depth information|B1 10m signal diagnostic|
|---|---|---|---|
|H01_M|RELATIVE|3.908666|878.7725|
|H01_M|ABSOLUTE_FLOOR|4.474304|1301.966|
|H01_P|RELATIVE|3.908666|878.7725|
|H01_P|ABSOLUTE_FLOOR|4.474304|1301.966|
|H06_M|RELATIVE|3.718802|3718.296|
|H06_M|ABSOLUTE_FLOOR|3.868981|1707.356|
|H06_P|RELATIVE|3.718802|3718.296|
|H06_P|ABSOLUTE_FLOOR|3.868981|1707.356|
|H12_M|RELATIVE|2.340749|7109.582|
|H12_M|ABSOLUTE_FLOOR|3.606617|2020.29|
|H12_P|RELATIVE|2.340749|7109.582|
|H12_P|ABSOLUTE_FLOOR|3.606617|2020.29|

Primary numerical rows STABLE: 12/12.
|Diagnostic|Worst primary value|Frozen limit|
|---|---|---|
|response_to_noise1_ratio|0.0522481399|0.5|
|Jacobian_mesh_relative|0.00711952724|0.02|
|depth_step_relative|0.000415465661|0.02|
|depth_information_relative|0.0515245138|0.1|

Other calibration/resource rows numerically stable: 180/180. Raw two-grid field differences are retained separately; new H3 relative-response tests do not change old E2's failed 0.2% raw-field Gate. All finite labels use exact cached source/receiver entries; local depth derivatives at 200m are step comparisons, not continuous-depth validation.

## Evidence and limits

Shared-reference covariance and direct source-nuisance information are cross-checked with dense covariance; data-processing inequality is checked against original known-source raw fields. All controls are in CONTROL_CHECKS.csv and read-only reconstruction in INDEPENDENT_AUDIT_CHECKS.csv.

F retains every 25 horizontal error nodes × five depth labels at 160001. These are nominal tangent nuisance projections of finite log/phase contrasts with source and C1b gain jointly removed across time. They are not nonlinear nuisance optimization, estimator recovery, full RC2 support, branch certification or a continuous-depth guarantee. U/U0: H3_FULL_HORIZONTAL_SUPPORT_NOT_ESTABLISHED. Six labels represent three identical-acoustics mirror pairs; no evidence multiplier is applied.

Noise is a proper-complex ideal-field delta-method design at 1%/5%, with separate frozen absolute floor. LOW_ENERGY_RECEIVING_DIAGNOSTIC.csv records weak field regions and checks whether local high-SNR approximation may be stressed. Unknown source level forbids interpreting these assumptions as actual SNR. No robustness conclusion to SSP/pose or signal extraction is drawn.

Proposed next REVIEW only: H3_EXTRACTION_REVIEW. No next execution authorized. No new KRAKEN/FIELD/MC/audio; no full-band E2. H3_EXTRACTED_OBSERVABLE NOT_OPENED. R4=0%. Commit B/push/verify and STOP.


Finite-horizontal conditional rankings: B1 source200 first at 62/150 relative-noise nodes and 64/150 fixed-floor nodes; duplicates reduce to 31/75 and 32/75 distinct-main-geometry nodes. These are descriptive local-chart rankings, not recovery rates or support coverage. H12 fixed-floor depth information retains approximately8.07% of the relative model; 5% reference floor can exceed the weakest field amplitude. Read INTERPRETATION_LIMITS.md before interpreting local information or considering any extraction review.
