# Failure scope and interpretation

All60 providers succeeded; independent reconstruction passed. No observed lower-window discrepancy: the1490/1480 full complex fields and centered SSP secants agree exactly within the stored numerical representation. All lower boundary margins exceed10m/s. Full mode counts and assignment results are retained; no modes were deleted to enforce equality.

The B label H3_P1_MODAL_WINDOW_NOT_CERTIFIED is caused by the preregistered joint lower/upper1m/s boundary-margin guard: all60 records fail its upper(CHIGH=1800) margin. The minimum upper margin is0.15054m/s. It does NOT mean the two CLOW windows differed, that negative SSP deleted modes, or that an omitted contribution has been observed. No wider-CHIGH calculation was performed, so the omitted upper-window contribution remains unbounded. This stage does not relax the pre-execution guard after observing the result.

The C label is independently supported by144/144 two-step and144/144 positive/negative secant checks failing their frozen5%/20% thresholds. The complete3-scenario,3-resource,2-noise,2-grid,2-window results are retained. Both source-only and C1b tangent diagnostics are present in SSP_SECANT_LOCALITY.csv.

Finite raw-field SSP secants followed by a nominal linear nuisance projection are the preregistered representation here. A finite common phase change and curvature can affect this representation before projection. This does not prove that every possible nonlinear source-invariant coordinate has nonlocal response. No alternative coordinate or smaller step was tested, and no such method is certified by these data. These failures are limits of this frozen numerical response contract, not negative evidence of physical depth identifiability or absence of SSP information.

The registered source-invariant principal-phase diagnostic itself passed. The low-energy and absolute-floor weighting controls retain weak received samples and their projected direction changes; amplitudes are model units, not measured source levels or SNR. No local C1e-P0 depth information was computed.

STOP: do not lower CLOW again, widen CHIGH, reduce steps or add frequencies/grids automatically. No current propagation-derivative-route continuation is opened. The next recommendation is a separately authorized RC2-horizontal-support/depth-set design review.

## Additional recorded aggregates

{
  "solver_calls": 60,
  "solver_total_elapsed_s": 471.7960000000894,
  "solver_call_max_s": 13.766000000061467,
  "same_count_across_tested_offsets_and_windows_by_frequency": {
    "150.0": [
      478
    ],
    "200.0": [
      638
    ],
    "250.0": [
      797
    ]
  },
  "all_sampled_assignments_matched": true,
  "min_sampled_shape_correlation": 0.27496303070088773,
  "lower_margin_checks_PASS": true,
  "upper_margin_checks_PASS": false,
  "grid_secant_relative_max": 0.008835052938637724,
  "window_secant_relative_max": 0.0,
  "step_checks_failed": 144,
  "forward_backward_checks_failed": 144,
  "principal_phase_checks_failed": 0
}
