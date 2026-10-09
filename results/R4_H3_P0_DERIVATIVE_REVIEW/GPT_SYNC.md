# H3-P0 saved depth derivative and modal-window review

Parent: ea021e88344167fc1fafcefe4f367631ff67efbe
Design SHA: f4e497af57f5fe7ab54de25c838b8374b64861bb
Execution SHA: commit containing this report; final local/remote SHA recorded separately after push.

## Frozen outcome
- H3_P0_DERIVATIVE_STENCIL_CONSISTENCY_SUPPORTED
- H3_P0_DEPTH_DERIVATIVE_ACCURACY_UNCERTIFIED
- H3_P0_MODAL_WINDOW_REDESIGN_REQUIRED

Original P0 remains NUMERICAL_OR_FORMULA_INCOMPLETE. Its source-column 2% Gate still has20/36 failures;126/126 original G0 replay is not revoked. SSP depth information NOT_EVALUATED; R4=0%.

## Three separate diagnostics
All five derivatives and all ten pairwise comparisons were retained across H01/H06/H12, B0/B1/B2, both meshes and both noise models. Effective information includes both registered sigma levels.

Stencil: projected h2 fitted ratio range 3.7634183–3.8556472, median 3.8375352; expected4 under a smooth asymptotic central-difference expansion. Maximum h2 vector residual 0.065095039. This is a truncation-pattern diagnostic, not a true-error bound.

R1/R2 maximum projected-direction difference 0.0023898496; maximum effective-information difference 0.0035506688. R1 effective information / oldD1 ranges 1.0281877–1.0723646.

Grid: high-order same-method maximum nuisance-profiled direction difference 0.0055528283; maximum effective-information difference 0.0040786972. Complete same-method grid comparisons appear in GRID_NUMERICAL_DIAGNOSTIC.csv.

Effective information is a conditional local tangent diagnostic. It is not achieved depth accuracy, an SSP-robust bound, a finite-sample estimator result or a global branch certificate. Null directions receive INF scale, never pseudoinverse zero variance.

## Sampling, cancellation and weighting
All source and receiver depths are exact saved samples. Source200m lies at a C-linear SSP interpolation knot; the left/right sound-speed slopes differ. This does not by itself prove derivative inaccuracy, but the smoothness needed for an unqualified fourth-order remainder bound has not been certified. Mode shapes are serialized complex float32, and R1/R2 share samples. No independent analytic derivative or rigorous remainder bound is available.

CANCELLATION_AND_LOW_AMPLITUDE.csv retains per-case minimum pressure, modal cancellation and stencil cancellation for each method. Modal cancellation measures sum of absolute contributions / magnitude of their complex sum; stencil cancellation measures coefficient-weighted pressure magnitude / derivative magnitude. A large ratio identifies sensitivity, not independently measured roundoff error. ABSOLUTE_FLOOR weights are sqrt(2)*|pressure|/REF, so weak-pressure observations contribute less than under RELATIVE whitening; these are frozen design assumptions, not measured SNR.

## Modal-window risk
All52 prepared positive/negative SSP environments and26 nominal environments were read; no perturbed modes were generated. Negative profiles have minimum1499m/s while CLOW=1500m/s; positive minimum1501m/s. CHIGH=1800m/s exceeds the sampled water speeds, but this alone is not a mode-convergence certificate. The narrowest nominal lower-window margin is 0.095416221m/s; 106 nominal entries lie within1m/s across both grids/frequencies.

The [official KRAKEN manual](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html) documents C-linear interpolation and exclusion of slower modes by nonzero CLOW. Therefore the present contract lacks a certificate that the ±1m/s SSP parameter difference uses a comparable physically relevant mode family. This is a design risk, not observed loss or predicted perturbed counts.

Any future contract must preregister a common window covering every SSP offset with boundary margin, or independently bound the contribution of excluded modes. It must check boundary/count/contribution convergence and rebuild the nominal provider using that same window. New and old nominal/perturbed providers must not be mixed. No new window is executed or authorized here.

## Independent verification and stop
Independent binary parsing, termwise modal summation, a separate nuisance chart and pivoted QR rebuilt 19257 checks: 19257 PASS, 0 FAIL. Max pressure discrepancy 1.0202342e-15; derivative 1.0296132e-15; effective information 2.5642669e-12. The original20 depth-column failures remain reconstructed.

New KRAKEN/FIELD/MC/audio=0. Internal consistency does not certify true derivative accuracy. The current local Fisher environment route stops at ACCURACY_UNCERTIFIED; a redesigned propagation contract is a separate future decision. No automatic SSP solve or source-depth expansion follows. This review contributes no negative physical evidence about SSP-eroded depth information, which has not been evaluated.

## Cancellation and noise detail

See CANCELLATION_AND_NOISE_INTERPRETATION.md for per-resource noise-weighting and stencil information changes, sample indices and cancellation ratios. The52 full environment contracts were additionally reconstructed:251 water-speed rows each shifted exactly?1m/s; all other fields unchanged. Decision classification was independently recomputed.
