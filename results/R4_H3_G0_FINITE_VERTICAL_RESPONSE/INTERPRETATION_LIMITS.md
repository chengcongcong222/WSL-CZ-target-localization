# Interpretation of saved evidence

The frozen local classification is H3_G0_VERTICAL_DEPTH_INFORMATION_PRESENT_CONDITIONALLY. It applies to the exact environment, exact source-200 center, ideal raw complex-field noise and C1b channel gains fixed over three snapshots. It does not establish globally reliable depth estimation.

## Finite horizontal conditioning is a material limit

At each of the 25 registered horizontal offset nodes, sort the five finite depth labels by the already saved C1b nominal tangent log/phase contrast. The source200 label ranks first in:

|Resource|Relative noise|Fixed floor|
|---|---|---|
|B0 original HLA|50/150|34/150|
|B1 vertical additions|62/150|64/150|
|B2 co-depth resource control|44/150|38/150|

These denominators include six identifiers, each with 25 nodes. M/P have duplicate main-acoustic geometry: B1 unique-geometry counts are 31/75 and 32/75. These are descriptive conditional rankings, not recovery rates: each node imposes a different horizontal state and the finite contrast uses local tangent calibration profiling, not full nonlinear nuisance optimization. The data demonstrate that the strong center-local depth signature does not guarantee preservation of the correct depth label under finite horizontal offsets. U/U0 remains H3_FULL_HORIZONTAL_SUPPORT_NOT_ESTABLISHED.

## Local information is an upper bound

At 1% relative ideal complex noise, B1/C1b depth information profiled over all horizontal tangents is 7722.41–505461.60 m^-2. C1b retains 56.62–67.24% of C0 depth information. Relative-noise B1/B2 information ratios are 2.34075–3.90867; fixed-floor ratios are 3.60662–4.47430. Thus the registered vertical contrast exceeds pure channel duplication within the frozen local model.

The very small formal inverse-sqrt-information depth scales (approximately 1.4–11.4mm for the relative design) are nominal Fisher diagnostics only. They are not achieved errors, continuous-depth guarantees, receiver signal-extraction results, or application commitments. Source derivatives are finite cached depth steps; random draws and estimators were not run.

At B1/C1b the scaled Jacobian has rank5 in both meshes without forcing rank. Its weakest state direction is predominantly heading psi (absolute normalized loading 0.99295–0.99828), coupled with bearing/speed. Identifiable condition numbers reach 5.2971e6. This experiment retains original HLA spatial bearing content; it does not use the old E2 P1 nuisance model that separately removed instantaneous bearing directions and imposed rank<=4.

## Noise and propagation limits remain open

The fixed absolute floor can suppress relative-noise depth information: H12 B1/C1b falls from 505461.60 to 40815.70 m^-2, retaining about 8.07%. The weakest nominal receiver amplitude is 3.63438e-6. Fixed-floor noise/nominal-amplitude reaches 0.29784 at the 1% reference setting, and 1.48919 at the 5% setting. Nonlinear phase/response extraction in such weak channels is not validated by local Fisher or delta-method calculations. These are declared source-level design assumptions, not real SNR.

Completely free complex modal gains absorb the depth tangent analytically. All 36 regular-source field configurations reproduce this absorption, and C2 gives rank0/INF variance. No independently justified restricted CZ group-gain model is admitted. C1e environment/deployment uncertainty remains NOT_EVALUATED. Any subsequent extraction review must address these restrictions and horizontal support; no extraction execution is authorized here.

There were no new KRAKEN/FIELD/MC/audio calls. The next H3_EXTRACTION_REVIEW is a review proposal only. E2 remains paused with its exact formula-guard scope, original E2 admission FAIL_UNCHANGED, extracted observable NOT_OPENED, R4=0%. Stop after execution commit and remote verification.
