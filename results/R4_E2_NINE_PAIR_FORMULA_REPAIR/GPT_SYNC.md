# E2 analytic formula repair: pre-development stop

Final: E2_FORMULA_REPAIR_FAILED. Formula Gate FAIL.
Parent 2574fd59cc04ae49d2fb9e7669a79b2d80376caa; design 45c747d0aebbf38b59d2f2fdfc84906d6f001af3.
Execution SHA is the commit containing this document. Verify remote main after push.

The one authorized repair execution reached all 48 scene/mesh analytic controls,
then stopped before development regression because 20 of 48 nominal pressure identity
checks exceeded the frozen 1e-12 relative bound. Maximum: 1.46661321453e-11.
The implementation mixed original Cartesian element ranges and the polar-coordinate
roundtrip used by the old pilot. Maximum range change was 1.45519152284e-11 m.
The new and old pressures were independently reconstructed using their respective
geometries. Their residuals were 0 and 3.52458087255e-15.
Thus this stop is a frozen numerical-consistency guard failure in this implementation,
not a demonstration that the physical mechanism lacks information. The frozen
threshold, source files, and execution marker were not changed and no rerun was made.

P1 primary physical rank: 48/48 PASS, ranks
[4]; maximum cofactor-null projected response
1.42755132269e-17.
Independent geometric SVD and QR nuisance elimination verified the same rank bound,
including both depth stencils; maximum null response 1.40714816194e-17.
P2 measured pre-control ranks: [5]. This is not a full P2 scientific validation.
Maximum analytic-vs-old-FD identifiable Jacobian discrepancy:
0.00222312540383.
Independent per-mode scalar analytic horizontal derivative discrepancy: 1.42436022305e-16.
Pre-control old/new effective-information changes: {"P1": {"range_relative_change": 0.03319022700490536, "depth_relative_change": 0.3802758549823611}, "P2": {"range_relative_change": 0.004092456403175905, "depth_relative_change": 0.0002065840178119294}}.
These are formula diagnostics only; they do not replace the cancelled regression.

P0/P1/P2 main C0/C1/C2 development comparison, scene-level data-processing regression,
unknown-depth range information stability, depth/velocity weak-direction stability,
relative/fixed-floor scientific sensitivity and candidate regression were NOT EXECUTED.
The information and stability CSVs explicitly contain stop statuses without estimates.
The registered synthetic plane-wave, gain-absorption and shared-frequency controls
were run; full covariance and scientific information closure remains unevaluated.
No application claim can be made from P2's pre-control rank or its diagnostic information.

Validation: 6327 PASS / 20 FAIL. Independent stop audit:
5042 PASS / 0 FAIL.
Every frozen input and original pilot artifact was hash-checked unchanged.
Original E2-G0: FAIL_UNCHANGED. Original pilot: IMPLEMENTATION_INVALID.
Original E2 scientific information: NOT_EVALUATED. R4=0%.
All new diagnostics: FORMULA_REPAIRED_DEVELOPMENT_ONLY.
New KRAKEN/Monte Carlo/audio: 0/0/0; development regressions: 0.
The single allowed repair opportunity is closed. No budget, depth step, SVD threshold
or pressure identity threshold was changed after observing results.
Next: H3_REVIEW (review only), requiring separate authorization for execution.
Do not launch full-band certification, H3, A2, SSP or P5.
