# B1A range-speed conditioning boundary

**B1_HORIZONTAL_CONDITIONING_ENVELOPE_BELOW_TARGET**. Execution valid, PENDING_RESEARCH_LEAD_AUDIT. Proposed B1/overall progress=0%; confirmed R4 remains0%. B2 stays NOT_OPENED.

Accepted parent B1 3fc3518161722890b7aa8da7fab524fea7b8f917: exact-horizontal depth mechanism ESTABLISHED, conditioning-tolerance envelope NOT_ESTABLISHED. The B1A design was pushed first at c3a9744f29ddc657c147539f135ccbecac5a3c1a; full fixed dyadic Cartesian lattice was then executed once, without changing cases, score, frequency/window weighting, depth labels, numerical rules or Gate. Existing B1 observations were read byte-for-byte, not regenerated. Historical B1/R3 and the closed CZ acquisition branch remain unchanged.

Six configurations are S7 triple (201/235/283 Hz) and S8 four-line (adds338 Hz), each at180/200/220 m: three source-depth scenarios, not six independent statistical samples. Theta=0 degree, psi=5 degree fixed. Range offsets0,+/-0.125,+/-0.25,+/-0.5,+/-1 km and speed0,+/-0.025,+/-0.05,+/-0.1,+/-0.2 m/s derive from old1 km/0.2 m/s steps by1/8,1/4,1/2,1. All486 profiles and10206 depth rows are retained; no failed point or quadrant removed.

Primary criterion: depth-label error<=10 m AND exact saved-value true-depth rank<=3. Secondary strict criterion: exact true label, rank1 and positive second-best margin greater than numerical tolerance. Strict performance is diagnostic and never substitutes for the primary criterion.

| Axis | Absolute offset | Exact | <=5 m | <=10 m | Rank<=3 | Practical pass | Worst error (m) | Min second margin (dB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| range | 0.125 | 1/12 | 3/12 | 4/12 | 3/12 | 1/12 | 45.0 | -1.63121 |
| range | 0.25 | 0/12 | 1/12 | 7/12 | 0/12 | 0/12 | 45.0 | -2.35295 |
| range | 0.5 | 0/12 | 1/12 | 3/12 | 0/12 | 0/12 | 70.0 | -1.29542 |
| range | 1.0 | 0/12 | 4/12 | 4/12 | 2/12 | 1/12 | 50.0 | -3.21498 |
| speed | 0.025 | 0/12 | 4/12 | 5/12 | 0/12 | 0/12 | 60.0 | -1.30993 |
| speed | 0.05 | 0/12 | 1/12 | 2/12 | 0/12 | 0/12 | 60.0 | -2.41248 |
| speed | 0.1 | 0/12 | 4/12 | 6/12 | 0/12 | 0/12 | 45.0 | -2.70914 |
| speed | 0.2 | 0/12 | 1/12 | 4/12 | 0/12 | 0/12 | 35.0 | -1.80843 |


Frozen primary target |dr|<=0.25 km, |dv|<=0.05 m/s includes150 profile points (25 signed lattice nodes x6 configurations). PASS requires all150, not an average. Actual practical pass=13/150, failures=137, worst error=60.0 m, worst rank=21. All16 candidate rectangles, including axes and interior signed nodes, are explicitly tabulated. Largest nondominated TESTED_ROBUST_RECTANGLE(s): []. A smaller passing rectangle cannot change the target decision or award10%.

| Configuration | All-lattice practical pass | Target practical pass | Target failures | Target worst depth error (m) | Target worst true rank |
|---|---:|---:|---:|---:|---:|
| S7 | 21/243 | 7/75 | 68 | 60.0 | 21 |
| S8 | 22/243 | 6/75 | 69 | 55.0 | 21 |


NON_MONOTONIC_DEPTH_CONDITIONING_RESPONSE observed=True; exhaustive same-case/same-closed-orthant coordinatewise smaller-failure/larger-success witnesses=369, including170 with the other axis fixed. Every witness is saved; sign-specific recovery islands are not smoothed into a monotonic boundary. This is a tested lattice boundary only: no spline, fitted boundary, between-node continuous safety guarantee or precise continuous tolerance claim. Four scientific figures are driven only by saved complete data. Heatmap blocks represent individual nodes on lattice indices; nonuniform physical dyadic spacing is explicitly labeled, not interpreted as uniform continuous cells.

B1 angle/heading evidence remains only the previously tested single-factor +/-0.5-degree theta and +/-1-degree psi controls (12/12 exact each). No angle sweep is repeated or expanded; these controls are not multiplied into a joint four-parameter robustness claim.

Numerical validation: repeat every full profile after fresh model reload, independently recompute all21 scores via the accepted independent struct/scalar-geometry/compensated-mode-sum/math.fsum path. Maximum cold repeat difference=0 dB; independent difference=4.07820444082e-11 dB. Frozen tol=max(128*eps*max(1,max_abs_J),10*max(new reconstruction residual,accepted B1 observation/score residual)) gives4.30873114965e-10 dB, not fitted to margins. 19177 runtime checks pass, including7273 historical file identities;630 score rows at old B1 origin and range/speed endpoints match accepted B1. 25 unit tests pass. Separate stdlib-only saved-score/metric/rectangle/frontier/nonmonotonic/tolerance reconstruction checks provide a second audit without physical reruns.

This identifies an upstream horizontal-state conditioning requirement in a synthetic matched E0 mechanism experiment, not platform sensor accuracy, an attainable horizontal estimate, cold-start acquisition, end-to-end observed depth accuracy, environmental robustness or final project depth specification. A passing tested rectangle would remain finite node evidence, not certification of all continuous points. S7 and S8 are both required by Gate; no four-frequency-only rescue. Failure below the frozen target is preserved without a finer lattice or downgraded target. Source-depth labels150:5:250 map to nearest stored modal samples; effective source/receiver depths are saved, with no off-grid z or interpolation.

CZ acquisition route CLOSED_FOR_CERTIFIED_GLOBAL_SEARCH; A1 remains BLOCKED/NOT_COMPLETED. No C4/F4, changed enclosure/threshold/grid/budget, 21-case CZ development, SSP/bottom/receiver-depth mismatch, bearing/platform/noise sweep or Monte Carlo. B2 and A2/A3/A4/B3/B4/C/P5 NOT_OPENED. Stop after execution commit/push and wait research-lead independent audit.

See [profiles](B1A_DEPTH_PROFILES.csv), [points](B1A_POINT_METRICS.csv), [rectangles](B1A_ROBUST_RECTANGLES.csv), [decision](B1A_DECISION.json), [validation](B1A_VALIDATION.json), [range plot](figures/FIG_A_RANGE_ONLY_DEPTH_ERROR.svg), [speed plot](figures/FIG_B_SPEED_ONLY_DEPTH_ERROR.svg), [depth-error map](figures/FIG_C_JOINT_WORST_DEPTH_ERROR.svg), [rank map](figures/FIG_D_JOINT_WORST_TRUE_DEPTH_RANK.svg).
