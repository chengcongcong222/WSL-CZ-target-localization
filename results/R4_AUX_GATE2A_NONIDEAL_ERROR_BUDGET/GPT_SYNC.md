# AUX Gate2A nonideal measurement requirements

AUX1_NONIDEAL_REQUIREMENTS_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED

All reported budgets are SINGLE FACTOR, not a simultaneous deployable error budget. No joint bias/navigation/deployment stress was run.

```json
{
  "bias_requirements": [
    {
      "anchor": "A",
      "baseline_km": 5.0,
      "sigma_random_deg": 0.05,
      "max_tested_abs_common_5pct_deg": 0.1,
      "max_tested_abs_diff_5pct_deg": 0.05,
      "max_tested_abs_common_10pct_deg": 0.1,
      "max_tested_abs_diff_10pct_deg": 0.1
    },
    {
      "anchor": "B",
      "baseline_km": 7.0,
      "sigma_random_deg": 0.075,
      "max_tested_abs_common_5pct_deg": 0.1,
      "max_tested_abs_diff_5pct_deg": 0.05,
      "max_tested_abs_common_10pct_deg": 0.1,
      "max_tested_abs_diff_10pct_deg": 0.1
    }
  ],
  "position_requirements": [
    {
      "anchor": "A",
      "baseline_km": 5.0,
      "sigma_random_deg": 0.05,
      "max_tested_sigma_pos_5pct_m": 50.0,
      "max_tested_sigma_pos_10pct_m": 100.0
    },
    {
      "anchor": "B",
      "baseline_km": 7.0,
      "sigma_random_deg": 0.075,
      "max_tested_sigma_pos_5pct_m": 50.0,
      "max_tested_sigma_pos_10pct_m": 200.0
    }
  ],
  "deployment_requirements": [
    {
      "anchor": "A",
      "baseline_km": 5.0,
      "sigma_random_deg": 0.05,
      "max_tested_abs_beta_5pct_deg": 20.0,
      "max_tested_abs_beta_10pct_deg": 20.0
    },
    {
      "anchor": "B",
      "baseline_km": 7.0,
      "sigma_random_deg": 0.075,
      "max_tested_abs_beta_5pct_deg": 20.0,
      "max_tested_abs_beta_10pct_deg": 20.0
    }
  ]
}
```

The two anchors are frozen by the research lead: A B5km/random0.05deg, B B7km/random0.075deg. B6km/random0.10deg is preserved only as a historical threshold reference; no new nonideal sweep for it. New zero-control quantiles use a new frozen30k sample table and may differ from Gate1 finite-MC results; they do not replace accepted Gate1 evidence.

## Model and interpretation

Random bearing terms independent zero-mean Gaussian at each anchor sigma. Fixed b1=b_common-b_diff, b2=b_common+b_diff; b_diff is half the difference, so a0.05deg bound corresponds to0.10deg b2-b1. Every signed81 combination kept. Common bias is a shared reference-like parametrization; differential bias is relative mismatch-like, not a verified physical decomposition.
Position error per-axis1-sigma in metres: both estimated node positions receive independent2D Gaussian offsets; observations use true geometry. Range output is norm(p_hat-MAIN_est), compared with true MAIN-to-target R; position error is norm(p_hat-target_true). Navigation-error common translations can cancel in relative range while remaining in absolute2D position; both metrics reported.
Deployment beta modifies true AUX=B[cos(90+beta),sin(90+beta)] and reflected placement. Estimator knows actual position exactly, so this is geometry tolerance, not navigation error. Beta is not converted into a bearing bias. Other nonideal factors zero in each single-factor family.
Exactly30000 trials per registered cell; one frozen30000x6 Gaussian table shared across cells. First2 coordinates bearing, next2 MAIN position, last2 AUX position; position arrays only used in POSITION family. Mirrors reflect noise and navigation y components. Bias mirror compares(c,d,side+) with(-c,-d,side-). Mirrors are controls, not independent trials.
All11 discrete50:1:60km ranges and mirrors included. Full-range P95=max over22 cells, never pooled or averaged across signs. Nearest-rank unconditional quantiles; parallel/nonfinite/behind errors assigned infinity, no outlier removal. Approximate P95 rank/Wilson intervals diagnose finite-MC uncertainty, not simultaneous assurance.
Axis frontiers require every registered inner signed node through the reported absolute cap to pass; common axis diff=0, differential axis common=0. Full bias map must be consulted for combinations: separate maximum axis limits cannot be applied simultaneously. Position and beta caps likewise preserve all lower tested magnitudes/signs. No interpolation or continuous tolerance certificate. A cap at the largest tested magnitude is right-censored; larger values are untested.

## Hardware and next decision

Actual auxiliary bearing array/capability, navigation/attitude, clock and calibration specifications remain UNKNOWN / CLIENT_CONFIRMATION_REQUIRED. Simulated arrays, synthetic0.1deg noise and ideal navigation assumptions are not hardware evidence. See source-bound AUX_HARDWARE_EVIDENCE_AUDIT.md and AUX_CLIENT_REQUIREMENTS.md.
Independent zero-mean random components may enter variance addition only if their independence is substantiated. Systematic biases stay separate; correlated random terms require covariance. No bias is silently folded into Gaussian RMS.
Numerical nonzero single-factor budgets do not admit an actual measurement chain. No Gate2B joint budget, time-skew/association MC, R4-A1-NEW, acoustic/depth/SSP/TDOA/tracker/5D. Single-array depth stays closed; R4=0%. Commit/push/verify then stop for research-lead audit and client confirmation.
