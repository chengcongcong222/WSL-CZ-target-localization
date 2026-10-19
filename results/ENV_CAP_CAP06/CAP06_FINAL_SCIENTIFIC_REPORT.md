# CAP06 final scientific report

The authorized receiver campaign is complete. New receiver evidence is separated from unresolved first-CZ transfer. No historical gate changes.

## O2 power ratios

O2_FIRST_CZ_FORWARD_NOT_ADMITTED. 960 toy sample records passed through power estimation, background subtraction and log ratios. Independent audit: 1947 checks, 0 failures. This is no depth-dependent first-CZ forward certification.

| SNR | Valid contrasts / total | Nonpositive channel power |
|---|---:|---:|
| -20 dB | 3603/11520 | 6742 |
| 0 dB | 11520/11520 | 0 |
| 10 dB | 11520/11520 | 0 |

Free differential gain gives algebraic rank zero. Linear drift can erode independent toy directions. Source cancellation does not cancel independent channel gain. Toy unweighted information is not a noisy-ratio Fisher estimate. Shared reference covariance, temporal correlation and missingness remain explicit.

## O6 extracted receiver products

| Holdout mode | SNR | Present epochs | Extracted | Covered | Axis median / P95 deg |
|---|---:|---:|---:|---:|---|
| CAL | 0 | 22648 | 22644 | 22643 | 1.11415 / 5.35751 |
| UNCAL | 0 | 22648 | 2968 | 1904 | inf / inf |
| CAL | 10 | 22648 | 22336 | 22316 | 0.393532 / 3.38229 |
| UNCAL | 10 | 22648 | 0 | 0 | inf / inf |
| CAL | 20 | 22648 | 21421 | 21378 | 0.184649 / inf |
| UNCAL | 20 | 22648 | 0 | 0 | inf / inf |

These epoch counts pool three nodes from256 independent base tracks, paired across qualities and resources. They are not independent acoustic trials. Epoch Wilson intervals are descriptive only. The main waveform bytes are shared between resources. All mirror angles and up-to3 spectral candidates remain unassigned, including leakage peaks. Truth frequency is evaluation-only. CAL uses extra known-direction signals; UNCAL receives no true bias or hard prior.

All fixed stress cases and failures are in P2_EXTRACTION_SUMMARY.csv. A specific frontend failure does not prove physical absence of information. No stress data retune the main gate.

## Information and conditional support

Local design matrices profile source coefficient and compare calibrated versus freely unknown node orientation. They are plane-wave upper bounds, not achieved accuracy. Pseudoinverse is only for nuisance projection; no null direction is assigned zero variance. Depth remains unobserved.

P3 was NOT_EVALUATED because the frozen receiver gate failed. No injected directed angles or truth-centered candidates replace unavailable receiver evidence.

## Resource requirements and decision

S0 has8 channels; S2 has16. Baselines are5 km or15.297 km. Equal20 dB input SNR at both nodes is assumed, not demonstrated CZ detection. Four known-angle high-SNR calibration transmissions per node, shape knowledge, raw baseband, noise-only background, gain stability, navigation/time bounds and association handling are substantive resources. Target depth is unobserved.

No P3 configuration was admitted, so no empirical accuracy-resource Pareto frontier was established. The Pareto deliverable lists resource costs with NOT_EVALUATED dominance and performance. It is not a global hardware minimum or real CZ accuracy frontier. Unknown multimodal coherence, seabed and SSP remain outside this plane-wave model. NUM03 remains nominal ray geometry. No first-CZ receiver transfer is established.

Retain the conditional receiver benchmark, calibration costs and ambiguity requirements. Do not expand baselines, SNR, angles or closed numerical chains automatically. O2 requires a certified same-CZ power-average contract. O6 requires a credible multipath receiver and calibration/association evidence before transferring old0.05-degree assumptions. No new5D result is available.

## Audit and execution

{'P1': {'PASS': True, 'checks': 1947, 'FAIL': 0}, 'P2': {'PASS': True, 'checks': 35975, 'FAIL': 0}, 'P3': {'PASS': True, 'checks': 1, 'FAIL': 0}}. Reconstruction does not replace physical admission.

Scientific wall=1406.0s; peak process RSS=0.367GiB; raw=3.467GiB; external solvers=0. All started units execute once.

Original main, historical gates and R4=0% remain. FINAL STOP.

## Sealed-output failure stratification (post hoc descriptive)
The preregistered CAL20 main gate failed in both development (5323/5632 covered,94.5135%) and holdout (21378/22648,94.3924%), below95%. CAL20 target extraction is94.5823%; the median axis error0.18465degrees coexists with unconditional P95=INF. No switch toCAL10 was allowed.

The complete saved-data stratification shows20326/20326 target-present holdout epochs without a distractor covered at every quality level. With a distractor, coverage is2317/2322 at0dB,1990/2322 at10dB,1052/2322 at20dB. These clustered epochs are not independent trajectories. This supports a failure concentrated in this processor's interference/model-fit handling as signal level rises; it is not an independent causal experiment and does not justify changing the residual threshold. Low-SNR acceptance sets are broader, so higher coverage alone does not prove better accuracy. The favorable subset cannot replace the frozen main panel.

At20dB the coherent-echo stress retains41/320 target angular sets; phase drift64/320; bad channel0/320; interchannel clock stress0/320. These are limits of this registered processor and synthetic model, not impossibility results for calibration or array methods. Uncalibrated target extraction is0/22648 at10/20dB. Noise-only false alarm control is0/320 per quality; non-target tones during target absence are a distinct association issue.

The positive evidence is actual calibrated single-arrival sample-to-axial-candidate extraction and the explicit gain/weak-energy boundary. The principal unresolved problems are a legitimate first-CZ receiver forward model forO2 and interference/calibration/association conditions forO6. No horizontal localization or depth performance was executed. See P2_SNAPSHOT_PHYSICAL_LIMITS.md for absent cross-window/cross-node coherent propagation.
