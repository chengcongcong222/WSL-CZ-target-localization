# CAP06 resource and research gaps

The physical gates and the receiver algorithm must be distinguished. Increasing node spacing cannot guarantee a usable extra observation. CAP06 assumes the same per-channel SNR for both nodes; it does not establish a CZ link budget. Calibration transmissions at known directions, raw complex samples, noise-background samples and stable channel responses are additional resources.

| Condition | Mathematical necessity | Tested hypothetical resource | Still unverified |
|---|---|---|---|
| Independent depth/range power information | Differential gain cannot span the full propagation contrast | G0-G4 toy, actual finite-sample power extraction | Same-CZ power average, real depth/range response, environment covariance |
| HLA directional information | Finite spatial phase slope; no truth-provided side | 8 elements,14 m,200 Hz,0/10/20 dB; both mirror branches | Real CZ multimodal receiver field, hardware calibration |
| Independent orientation reference | Unknown channel phase and attitude cannot be silently fixed | Four known-angle2048-sample40 dB transmissions/node | How these calibrations could be obtained and maintained at sea |
| Horizontal support | Observable trajectory variations must survive nuisance freedom and branch unions |31 snapshots over1200 s; B0 and B3; bounded Cartesian outer search | Continuous optimum, statistical simultaneous angular coverage, actual association |
| Target depth | An independent depth-sensitive observable is required | None admitted to P3 | NOT_EVALUATED; cannot infer from extra bearing node |

Useful necessary scales are symbolic rather than invented accuracy specifications. For an ideal white-noise plane wave with complex amplitude profiled, angular information is proportional to (2 pi f/c)^2 SNR sum_m (x_m-xbar)^2 sin(alpha)^2. Endfire is locally weak; a small axial-cosine error need not be a small angle error. At exact endfire first-order angular sensitivity vanishes. The frozen1001-point cosine grid and finite snapshot span are also processing limits, not a physical impossibility theorem.

For toy relative intensity d=Dp with common source eliminated, effective state sensitivity is the part of D J that survives projection against differential gain, environment and other states. Arbitrary differential gain makes that component zero. More snapshots alone do not remove this structural degeneracy. Shared reference noise is correlated; the unweighted toy matrix is only algebraic mechanism evidence.

For two independently directional nodes, larger baseline may improve geometrical intersection, but bearing calibration, source association, navigation and range-dependent reception remain necessary. No numerical hardware minimum follows from this campaign. The finite Pareto CSV compares the same matched model and registered resources only; it must not be mixed with historical oracle-bearing T3 precision or H3 local complex information.

No partner input is required to finish this hypothetical study. Hardware, ocean reception, calibration availability and real signal stability remain UNKNOWN, rather than being silently assumed established. No external messages are sent. Future investigation requires a new authorized contract; the present campaign stops after final delivery.
