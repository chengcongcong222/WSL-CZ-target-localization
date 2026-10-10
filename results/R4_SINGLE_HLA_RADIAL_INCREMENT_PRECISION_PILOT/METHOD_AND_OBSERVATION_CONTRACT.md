# HLA-H2 frozen conditional feature contract

Parent: dd13a0a32a4bd279c03c272b9ee8b86ceea14730. This is an adapted feature-level experiment, not a replication of a paper or actual HLA feature extraction. Source amplitude/phase, emitted frequency and propagation have not been solved. R4=0%.

## Geometry and conditions

h=[r_km,theta_deg,v_mps,psi_deg], domain [45,60] km × [-5,5] deg × [1,3] m/s × [-15,15] deg. Angles are counterclockwise from +x, not the north-axis convention in the historical literature lock. Target p0=1000r(cos theta,sin theta), CV velocity v(cos psi,sin psi). Platform speed 2 m/s along +x through 600 s, then 15 deg; trajectory continuous. Navigation is exact, signed bearing error is independent Gaussian 0.1 deg, system bearing bias zero. Actual single-HLA left/right ambiguity resolution is not certified.

121 bearings at 0:10:1200 s. Six nonoverlapping 200 s windows at endpoints 0:200:1200. qbar=(rho(end)-rho(start))/200 in m/s uses exact endpoint distances. It is signed NET displacement rate, not instantaneous radial speed or window mean absolute speed. Within-window sign crossings and near-zero net rates are reported for all 48 geometry/windows and never excluded.

All eight H1 off-grid states are the DEVELOPMENT_GEOMETRY_PANEL. Archived z labels are metadata only. No extra array, auxiliary node, depth observation, H1 acoustic fusion or altered physical domain.

## Observation data and unknowns

B: y_beta=wrap(beta+sigma_beta e_beta). No q.
S: y_q=qbar+sigma_q e_q, signed calibrated stronger control.
U: u=abs(qbar+sigma_q e_q), only nonnegative observations, calibrated zero. Folded-normal observation family, not abs(qbar)+Gaussian.
Biased data: u=abs(qbar+b_true+sigma_q e_q). Producer sets b_true=+0.10 m/s at sigma=.05. U-mismatch ignores b; UB-matched knows only b in [-.20,.20], shared across all six windows. The bias is a conditional surrogate, not a complete physical source-frequency model.

Sigma_q .02/.05/.10/.20 m/s; 32 replicates × 8 geometries. B256; S1024; U1024; biased U-mismatch256; UB256: 2816 configurations. Only 256 independent base innovation pairs; all scales/methods share each pair. SeedSequence([2026101002,zero_based_geometry,zero_based_replicate,stream]) and NumPy PCG64; stream0 bearings, stream1 radial. Noise between bearing/radial and radial windows is independent by ASSUMPTION, not by validated receiver extraction. S data live in a separate file; unsigned estimator API has beta, radial, sigma, kind only. Truth and innovations stay in generator/evaluator files, and never enter initialization or branch selection.

Primary: U, sigma=.05, zero bias. Stronger S cannot rescue the primary result. No main eight-scene performance runs before Commit A.

## Scores and conditional coverage

T_beta=chi2.ppf(.975,121); T_q=chi2.ppf(.975,6), numerical values in DESIGN_FREEZE.json. No subtraction of fitted parameters. Accept BOTH J_beta<=T_beta and J_radial<=T_q. S uses signed residual. U uses abs(q)-u residual. UB completely minimizes abs(q+b)-u over the registered b domain at fixed terminal h.

At the fixed true h, the Gaussian squared-noise events each have probability .975; the inequality ||a|-|b||<=|a-b| and a bias domain containing true b can only reduce the corresponding true-state radial residual. The union bound gives at least .95 ideal-continuous-set conditional inclusion for correctly modeled S/U/UB, and .975 for B alone. This is NOT coverage of the exported finite terminal set. U score is not an exact chi-square distribution, folded-normal MLE, LR interval or global support certificate. U-mismatch has the wrong zero-point model and does not inherit the stated biased-data coverage.

## Observation-only estimator

Enumerate all 64 lexicographic {-1,+1}^6 expansions of NOISY unsigned observations; no true sign, monotonic-q or branch pruning. Each branch fits the single CV trajectory to all original bearings and six independent window features. Smooth branch residual q+b-s*u; signed S uses its signed observations.

Endpoint observed bearing vectors e_j and d_j=200 sum_{i<=j} unfolded_q_i yield
p0+v*t_j-rho0*e_j+b*t_j*e_j=platform(t_j)+d_j*e_j.
Scaled SVD rcond=1e-12 gives an initializer only. Physical scales [52500 m,52500 m,2 m/s,2 m/s,52500 m,(.2 m/s)]. Rank, condition, nonpositive rho0 and clipping logged; no configuration discarded. Nonlinear distance consistency, angular measurement error and correlated cumulative increments prevent treating this initializer as unbiased estimation.

U/UB/mismatch: two starts per branch (linear geometric and domain center [52.5,0,2,0], UB b0). S: these two starts, observed signed branch. B: bearing-only Cartesian initializer plus 16 fixed range strata 45+15*(i+.5)/16; theta0,v2,psi0. Total expected optimization records 203008.

SciPy least_squares TRF, linear loss, normalized physical box, analytic Jacobian; max_nfev100; ftol=xtol=gtol=1e-9. Bounds enforced. Distance uses km/m factor1000; angle columns degree/radian factor pi/180. No post-result solver repair. Physical residuals reconstructed at all terminals. UB stores original b/J_branch/J_folded and the independent fixed-h complete piecewise quadratic profile b/J_profile separately. Profiling changes b only, never h.

Choose minimum J_beta+J_profile among accepted records, tie by terminal_id. None accepted => ABSTAIN, rejected best diagnostic kept separately. A finite nonconverged optimizer point may be accepted if it satisfies physical compatibility; status, bound touches and exceptions remain visible. No search-found result is a continuous global optimum certificate. Dedup scales [.005 km,.001 deg,.005 m/s,.05 deg,.005 m/s bias] only label terminal clusters; raw records remain.

## Metrics and predeclared descriptive decisions

Errors: relative range and speed, absolute smallest circular theta/psi differences. Every scene/method/sigma: accepted count, abstain/exception/nonconvergence, true compatibility, search failure, conditional and INF-inclusive unconditional median/P90/P95/max, denominators, Wilson95 output-rate interval. Nearest rank; n32 P95=31st order statistic, exploratory. NOT_EVALUATED is missing, never a zero/failure. Across scenes use equal-scene means of individual medians and worst individual-scene P95, not pooled replacement.

Primary positive signal >=5/8 scenes, each full32 observations, output>=.90, and >=.30 unconditional-median improvement vs paired B in at least r/v/psi. B median0 cannot use a ratio. B medianINF only qualifies a direction with finite U P90 strictly better than B P90; do not use infinity ratios.
Sign-information-required: U primary <5 and S .05 >=5 qualifying scenes.
Reference-critical: >=5 scenes with UB output rate at least10pp below clean U or any r/v/psi unconditional median at least1.30 times the finite positive U median.
Velocity-gain-with-range-limit: >=5 primary scenes with >=.30 median velocity improvement and any primary scene unconditional range P95 >.10.
Continuous-search-limited: any true-compatible configuration with no accepted searched terminal. Otherwise limited-at-tested-precisions if primary <5. Implementation invalid/partial supersede scientific classifications. Criteria are exploratory route signals, not client acceptance or proven real success probability.

Reference lines rangeP95 .10; speedP95 .10; headingP95 5deg; bearingP95 1deg, each tested bin separately. No interpolation of a critical precision, no joint-four-parameter Gate. Report all finite terminal envelopes and farthest wrong terminals; these are NOT confidence sets or full continuous outer support.

## Verification, budgets and stop

Preflight seed2026101012, fixed mathematical controls only; independent Cartesian geometry and finite differences, analytic Jacobian tolerance1e-5, scores tolerance1e-7 relative, all-sign bias profile, net displacement/integrated radial rate, sign flip/wrap, self-match, abstain/INF/count/API/summary controls.

One execution after A pushed and local/remote verified. Single process, BLAS4, memory8GiB, total main+cold wall14400s; main stops at12600s reserving1800s cold; delivery<=512MiB. Configuration order geometry→replicate→B→each sigma S/U→mismatch→UB. Flush all raw terminals and completed estimates, preserve partial and remaining NOT_EVALUATED. Independent audit every terminal (geometry, all scores, b minima, acceptance), configuration/branch counts, selected output and all data innovations; eight registered U/.05/rep0 configurations get independent terminal FD-Jacobian audit with no new optimization. Report metrics/Gates independently rebuilt.

A freezes code/inputs/thresholds/starts/budget. B publishes once, verifies remote and STOP. After A no estimator, inputs, threshold or start changes; implementation/budget issues retain original data and stop. No KRAKEN/FIELD/BELLHOP, audio or raw complex spectrum, H1 rewrite, H3, next mechanism, enhanced architecture or real extraction. New work is feature MC plus pure geometry. Old R4 remains0%; COLD_START_ACOUSTIC_CAPABILITY_NOT_ESTABLISHED.
