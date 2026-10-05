# Frozen AUX Gate1 cross-track requirement boundary

Parent79bc6b5c04e08d6074719e3c045947c1325253c8 Gate0 accepted conditional. Current single-array five-parameter engineering target NOT SUPPORTED; AUX1 conditionally practical geometry, not known engineering capability. New stage only cross-track requirement map.

MAIN=(0,0), target=(R,0), AUX=(0,+/-B). B=2:0.5:10km (17 values), R=50:1:60km (11 values), effective directed-bearing RMS=0.05/0.075/0.10/0.15/0.20deg (5 values); both mirrors,1870 cells. Alpha follows atan(B/R); no per-range placement/angle adjustment. Exact positions, simultaneous time, correct association/directed bearing, ideal reference main-LoB alignment. No position/bias/time nuisance sweep.

50000 trials per mirror cell, one PCG64seed2026100502 frozen before design push. One shared50000x2 standard-normal table; Gaussian zero-mean independence between nodes. Cell and mirror comparisons are not independent repeats. Mirror sign reverses both perturbations to isolate symmetry. Estimator uses known positions plus two noisy directed LoBs only. Truth sets scene/evaluation, never clips or initializes the estimator.

All failures retained as infinite error: parallel determinant abs<=1e-12/nonfinite/behind either sensor. Finite ill-conditioned cases retained; report1/abs(sin(observed crossing))>1000 separately. Nearest-rank unconditional median/P90/P95/P99 for relative range, absolute range,2D error and cross-range. Range=distance(MAIN,estimated target). No deletion/outlier trimming.

For each fixed B,sigma: max empirical P95 over all11 ranges and mirrors; full-range requires all cells execution-valid and all P95 finite. T5<=5%, T10<=10%. Not continuous50-60km guarantee. Frontier=minimum TESTED B reaching each threshold, or NOT_REACHED_WITHIN_B<=10KM. FixedB5/10 requirement=max TESTED sigma passing; preserve next tested failure bracket, no interpolation/extrapolation. Seed/Gate/grid/threshold unchanged after results.

Outcome A: at sigma0.10deg someB<=5km T5. Outcome B: no A but some sigma>=0.05deg B<=10km T5; tradeoff geometry only. Outcome C: no T5 combination. No hardware or five-parameter progress even on A/B. R4-A1-NEW, Gate2 not opened; depth closed pending independent horizontal acquisition; R4=0%.

Independent slope intersection cold check and independently derived local covariance. FIM inverse is local linearized lower-bound/approximation, not actual-estimator guarantee. Local P95=1.959963984540054*sigma_r/R. Check sigma/sin(alpha) scaling descriptively, not Gate. Finite-MC P95 rank intervals approximate binomial-normal95%; Wilson95% fraction<=5/10% intervals; no selected-grid simultaneous guarantee. No unconfirmed SNR/snapshots/aperture used to promise sensor accuracy.

Freeze Commit A, push/verify remote before any new MC or numerical tests. Execute once unchanged, independent audit, Commit B/push/verify then STOP. No propagation/CZ/depth/SSP/TDOA/5D/tracker/off-gridA1. No new grid points after results.
