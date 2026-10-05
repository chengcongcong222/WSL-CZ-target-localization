# AUX Gate1 cross-track baseline and bearing requirement

AUX1_BASELINE_BEARING_TRADEOFF_ESTABLISHED

Parent79bc6b5c04e08d6074719e3c045947c1325253c8 accepted conditional Gate0. Gate1 uses a new frozen deployment, not a retuned Gate0 triangle. MAIN=(0,0), target=(R,0), AUX=(0,+/-B); alpha=atan(B/R). AUX angle is never adjusted per target range.

## Requirement frontier

| Effective bearing RMS deg | Minimum tested B for5% km | Minimum tested B for10% km |
|---|---|---|
| 0.05 | 3.0 | 2.0 |
| 0.075 | 4.5 | 2.5 |
| 0.1 | 6.0 | 3.0 |
| 0.15 | 9.0 | 4.5 |
| 0.2 | NOT_REACHED_WITHIN_B<=10KM | 6.0 |

All11 tested ranges50:1:60km and both mirrors must be valid/finite; worst nonlinear MC P95 governs T5/T10. Finite sample requirement on a discrete grid, no continuous range certificate, no interpolated B/sigma threshold.

Fixed-baseline bearing requirements:
```json
[
  {
    "baseline_km": 5,
    "largest_tested_sigma_for_5pct_deg": 0.075,
    "next_tested_failing_sigma_for_5pct_deg": 0.1,
    "worst_p95_at_5pct_boundary": 0.04413373847551559,
    "largest_tested_sigma_for_10pct_deg": 0.15,
    "next_tested_failing_sigma_for_10pct_deg": 0.2,
    "worst_p95_at_10pct_boundary": 0.08894963902966377
  },
  {
    "baseline_km": 10,
    "largest_tested_sigma_for_5pct_deg": 0.15,
    "next_tested_failing_sigma_for_5pct_deg": 0.2,
    "worst_p95_at_5pct_boundary": 0.044613563497423785,
    "largest_tested_sigma_for_10pct_deg": 0.2,
    "next_tested_failing_sigma_for_10pct_deg": "NOT_BRACKETED_ABOVE_TESTED_MAX",
    "worst_p95_at_10pct_boundary": 0.059488856021590306
  }
]
```

Primary sigma0.10-deg selection:
```json
{
  "baseline_km": 6.0,
  "alpha_min_deg": 5.710593137499642,
  "alpha_max_deg": 6.84277341263094,
  "worst_range_km": "60",
  "worst_p95_relative": 0.049110591639527475
}
```

## Sampling and propagation of error

Each of1870 registered mirror cells evaluates50000 trials. The shared single50000x2 PCG64 standard-normal table compares cells; cell-trials are not independent trials across geometry. Mirrors reverse both Gaussian perturbations and are symmetry controls. No repeat seed, new scenario or post-result grid addition. Independent nodes have zero-mean Gaussian effective directed-bearing error of equal RMS; truth enters scene/evaluation only, never the intersection estimator.
Parallel/nonfinite/behind-sensor intersections receive infinite errors and remain in unconditional median/P90/P95/P99. Finite poorly conditioned intersections remain. Failure-rate metrics may overlap. Empirical nearest-rank quantiles, no outlier trimming. P95 rank intervals are approximate binomial-normal diagnostics; Wilson intervals describe proportions<=5/10%. No simultaneous selected-grid or deployment guarantee.
The local linearized covariance check uses an independently derived sensitivity. sigma_r= sigma_rad*sqrt(R^4+(R^2+B^2)^2)/B. At B/R small this gives P95 approximately1.96*sqrt(2)*sigma_rad/sin(alpha). Ratio/scaling diagnostics use actual MC, not surrogate Gate results. Local inverse FIM is a linearized lower-bound/approximation, not estimator guarantee.
Nonlinear/linear P95 ratio range: 1.007706--1.121796. MC P95 divided by sigma_rad/sin(alpha): 2.793313--3.109404. Baseline monotonicity violations: 0. See AUX_SCALING_DIAGNOSTIC.json.

## Ideal geometry and measurement requirement

MAIN/AUX positions exact, timestamps simultaneous, association correct, directed bearing resolved. Perfect cross-track alignment to the reference main LoB is a Gate1 ideal deployment assumption; uncertainty of initial bearing-guided placement/relative alignment is not simulated here. Node positioning, attitude, time offset, systematic bearing bias/correlation are Gate2 topics and NOT_OPENED.
Effective directed-bearing RMS is a total input requirement, not beamformer-only CRLB or known receiver capability. Future measurement budget must include DOA random error, calibration residuals, heading/attitude and relative alignment. Gaussian zero-mean independence excludes systematic/common biases in this Gate. No numerical allocation to components is inferred, no SNR/aperture/snapshot guess.
Existing project materials leave actual element count/aperture/navigation/heading/bearing accuracy UNKNOWN / CLIENT_CONFIRMATION_REQUIRED. Historical model arrays and0.1-degree synthetic noise are not device acceptance. See AUX_BEARING_MEASUREMENT_REQUIREMENT.md and source-bound original excerpts.

## Architecture boundary

Single-array cold start NOT ESTABLISHED; single-array depth engineering CLOSED; exact-horizontal mechanism oracle evidence retained. AUX1 remains CONDITIONAL_ARCHITECTURE_CANDIDATE, even if the requirement tradeoff is established geometrically. No actual auxiliary bearing chain, node pose/time/association/deployment evidence supplied. No velocity/heading/depth precision claim.
R4-A1-NEW NOT_OPENED; depth CLOSED_PENDING_INDEPENDENT_HORIZONTAL_ACQUISITION; AUX Gate2 NOT_OPENED. R4=0%, original A/B sequence paused; A2/A3/A4/B2/B3/B4/C paused/not opened; P5 not opened. No acoustic TL/CZ/depth/SSP/TDOA/5D estimator/tracker. Commit/push/remote verify and STOP.

## Boundary interpretation and verification

At primary effective directed-bearing RMS0.10deg, minimum tested cross-track baseline is6.0km for T5 and3.0km for T10; worst range60km. Selected T5 crossing spans5.710593--6.842773deg. These are minimum tested values on0.5km spacing, not continuous minima; no5.75/6.25km additions or interpolated thresholds.

At fixed B5km, largest tested RMS passing T5 is0.075deg (next tested0.10deg fails); largest passing T10 is0.15deg (next0.20deg fails). At B10km, T5 largest tested0.15deg (0.20deg fails); T10 passes at highest tested0.20deg, with no tested upper failure bracket. Sigma0.05deg T10 already passes at the smallest registered baseline2km, so the lower boundary is left-censored by this grid; no smaller-baseline requirement is asserted.

All1870 cells are execution-valid and P95 finite; all85 fixed(B,sigma) summaries include all11 ranges and both mirrors. Eight tests and98240 independent cold geometry checks pass, zeroFAIL. Supplementary checks bind the decision JSON to CSV frontiers and verify original hardware evidence hashes. No primary MC rerun or policy/source correction was needed. No R4 credit or automatic next-stage release.
