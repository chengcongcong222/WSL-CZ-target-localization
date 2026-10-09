# H3-P1 unified modal window and local SSP response

Parent: 7aaab2b99150c11b6d9a18092d4a97076d1c64af
Design SHA: f827bd38fc86e0b12ff679956b3b0fab0f559d2e
Execution SHA: commit containing this report; final verified local/remote SHA stored after push.

## Result and scope
- H3_P1_MODAL_WINDOW_NOT_CERTIFIED
- H3_P1_SSP_SECANT_LOCALITY_NOT_ESTABLISHED

This is a frozen three-frequency numerical screen, not a source-depth/SSP identifiability or depth-accuracy claim. OriginalP0 remains H3_P0_NUMERICAL_OR_FORMULA_INCOMPLETE. True source-depth derivative accuracy remains ACCURACY_UNCERTIFIED. C1e-P0 depth survival NOT_EVALUATED; complete RC2 support NOT_ESTABLISHED; R4=0%.

## Provider
New KRAKEN calls 60/60; generated files 60/60. Every window includes its own new nominal baseline. No old1500m/s field is used as a new-window nominal. Frequencies150/200/250Hz; grids80001/160001; CLOW1490/1480m/s; CHIGH1800m/s; water-speed shifts0,±0.01,±0.02m/s. All other frozen environment fields and13 exact source/receive sample depths are unchanged.

CHIGH has not received an independent truncation certification. Agreement between the two lower windows is evidence only within this tested contract, not a universal propagation completeness theorem. The SSP step sizes are locality probes, not a measured ocean error range.

## Frozen diagnostics
{
  "failed_checks_by_type": {
    "LOCAL_STEP": 144,
    "LOCAL_FB": 144,
    "WINDOW_MARGIN": 60
  },
  "two_step_C1b_relative_max": 0.2775716771812911,
  "forward_backward_C1b_relative_max": 1.4784434172486591,
  "source_removed_phase_max_rad": 0.4515005815171616,
  "CLOW_margin_min": 10.07543564252569,
  "CHIGH_margin_min": 0.1505407165896031,
  "GRID_raw_relative_max": 0.0005873553477335053,
  "GRID_source_response_relative_max": 0.00033754408642788927,
  "GRID_response_noise1_ratio_max": 0.05781259976070951,
  "WINDOW_raw_relative_max": 0.0,
  "WINDOW_source_response_relative_max": 3.0256464851385373e-16,
  "WINDOW_response_noise1_ratio_max": 0.0
}

All three H01/H06/H12 scenarios and B0/B1/B2 resources, two meshes, two windows and both inherited noise assumptions are retained. Full raw-field, source-removed spatial-projector and fixed-gain-profiled response comparisons appear in RELATIVE_RESPONSE_NUMERICAL_STABILITY.csv. Secants use only the new nominal baseline and are compared at both steps; no source-depth D1 enters them.

C1b removes unknown complex source per frequency/time and complex element/frequency gain fixed across all snapshots. Principal phase is measured on element/reference source-invariant ratios; it is not a certified unwrap or infinitesimal derivative. Low-pressure locations and absolute-noise-floor sensitivity remain in the per-case tables.

## Complete modal sets
MODE_BOUNDARY_AND_CONTINUITY.csv records counts, phase speeds and boundary margins. MODE_MATCHING.csv records sampled-shape Hungarian assignment, correlations and apparent unmatched indices. Matching based on13 depth samples is not an exact identity theorem. Count differences alone do not fail a Gate. No modes are discarded.

WINDOW_FIELD_CONTRIBUTION.csv rebuilds unmatched and low-correlation subsets while retaining the complete full field for all scientific comparisons. Assignment is only diagnostic; SSP secants are differences of full fields. Full-window comparisons, both grids, both secant steps and nuisance-profiled reproducibility jointly determine admission.

## Independent review
3555 checks; 3555 PASS; 0 FAIL. Cold pressure relative maximum 1.3774630051071827e-15. Independent binary parser, Cartesian geometry and termwise modal summation reconstruct every saved field, with a separate nuisance chart and pivoted QR for response and secant checks. Scientific Gate failures are distinct from implementation reconstruction failures.

## Stop
No FIELD calls, MC, receive recordings or continuous depth localization. No automatic lower window, smaller SSP step, new frequency/grid or original±1m/s stress trial. If numerical support is accepted, only a separately designed full-band research-lead review is suggested. Otherwise stop the current propagation derivative route and prioritize a separate RC2-horizontal-support/depth-set design decision. No next stage is opened here.

## Essential failure-scope qualification

The B label is an upper(CHIGH) boundary-guard failure, not an observed difference between the CLOW windows. Their full complex fields and SSP secants coincide. The C label independently records failure of the frozen raw-field-secants/nominal-projection locality tests. These results do not disprove locality in every nonlinear source-invariant coordinate. See FAILURE_SCOPE.md and LOCALITY_RESOURCE_AND_NOISE_SUMMARY.csv.

## Supplemental mode-diagnostic cold review

A further1638 checks independently reconstructed normalized sampled-shape matching and complete low-correlation subset contributions:1638 PASS,0 FAIL. Combined with3555 main checks:5193 checks,0 FAIL. The supplementary script is MODE_DIAGNOSTIC_RECHECK.py; run from the repository root via python -c "import runpy;runpy.run_path('results/R4_H3_P1_UNIFIED_WINDOW_SSP/MODE_DIAGNOSTIC_RECHECK.py',run_name='__main__')". This saved-data review used0 additional solvers and did not change any Gate or classification.
