# Frozen gauge and nuisance contract

Inference sees observed complex coefficient samples and independent background records only. No modal fields,
nominal q/g, source draws or generating noise covariance enter the PRIMARY statistical interface.
Oracle covariance is supplied only to a separately labeled call.

Unknown q varies over frequency/template. Fixed unknown complex g varies over element/frequency, not template.
Multiplying g_f by c_f and dividing q_tf by |c_f|^2 leaves covariance unchanged; phase is unobservable.
For arbitrary fixed per-element h_f, (g_f h_f)*(p_tf/h_f) also leaves all three template responses unchanged.
The extractor therefore estimates a gain-distorted response projector, not an uncalibrated physical field.
No source/gain truth, calibration prior or pseudoinverse-based zero variance is supplied.

The full raw and CSD likelihood witnesses use identical synthetic non-truth nuisance candidates fixed in code.
They test the Gaussian covariance family and sufficiency, not a fit of the physical propagation model.
All full channel-frequency blocks, including shared Hann frequency noise, enter the model.
The primary likelihood is K logdet C + K trace(C^-1 Rhat), with separate background likelihood.
Sample CSD may be singular. Only MODEL/background covariance is factored; no sample CSD inverse or ridge.

The conditional proper-Gaussian source model differs from H3-G0's deterministic unknown complex means.
No Fisher efficiency ratio between these distinct experiments is admissible.
The normalized eigenvector statistic has no established finite-sample likelihood and receives no Fisher efficiency.

For a single separable component, source-depth amplitude is absorbed by unknown q.
For free per-mode complex gains, changing nonzero source modal coefficients is compensated mode by mode.
For C2 free element/frequency/template responses, any source-depth change is likewise absorbed.
These are exact zero-information controls, not evidence against constrained propagation models.
Constrained propagation C1g, environment/pose C1e and complete RC2 support remain NOT_ESTABLISHED.

The 256 background samples support an unconstrained proper complex full covariance because B exceeds
the largest selected joint dimension 156 (union archive dimension 208). Positive definiteness is still checked.
Background samples and source/observation noise use distinct frozen RNG substreams.
No background regularization or observed-data-dependent covariance restriction is allowed.

Per-frequency compressed direction is released only above a registered energy threshold.
With independent complex Gaussian background B, E[Rn_hat^-1]=B/(B-N) Rn^-1;
Markov yields the conservative per-block alpha=0.01 threshold. This is not an optimized detector.
Undetected cells retain failures and no inferred direction. Quality thresholds are descriptive extractor
criteria chosen before data, not a changed project accuracy Gate. Detected-only bias and variance cannot
be interpreted as unconditional estimator risk; all failure denominators remain.

K and package comparisons, noise families/scales, resources and two meshes are paired.
Only 48 independent innovation pools exist. The original eight noise channels are exactly shared;
four vertical and four co-depth supplemental channels have independent artificial noises.
Three template labels do not establish time-domain windows or 64-second stationarity.
