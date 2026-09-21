# B3D-3 Gate Bottleneck Audit — 2026-09-09

## Purpose

This is a derived audit of the existing `B3D3_ANCHOR_ADMISSION_TABLE.csv`. It does **not** replace the official admission table and does **not** upgrade any candidate to an admitted physical anchor. Its purpose is to identify which gate first collapses the mature-model candidate pool, so that the next B3D-3 run spends effort on the binding uncertainty rather than adding more features.

## Evidence level

- Propagation source: existing mature `bellhopcxx 2-D` run in the B3D-3 package.
- Environment in that run: reconstructed South-China-Sea-type SSP, **not observed SSP**.
- This audit: deterministic post-audit of the 26 rows in the existing admission table.
- No new Bellhop/Kraken/RAM propagation solve was executed in this audit.

## Gate waterfall

| Stage | Gate | Survivors | Interpretation |
|---|---|---:|---|
| 0 | All mature-model relative-delay candidates | 26 | Starting pool in existing admission table |
| 1 | Adequate common range window | 25 | Only one candidate is eliminated solely by a short common window |
| 2 | Numerical delay stability | 23 | Two additional candidates fail the existing 0.075 ms numerical-change diagnostic |
| 3 | Arrival resolvability at the existing B=900 Hz criterion | 9 | Unresolved neighboring arrivals are a major but not universal blocker |
| 4 | Not weak relative to the strongest single arrival | 3 | Only three candidates remain attractive before environment robustness |
| 5 | Stable association in all 8 diagnostic environment perturbations | **0** | First gate that annihilates the pre-environment shortlist |
| 6 | Environment delay shift <= 0.15 ms | 0 | Also failed by every candidate, but not needed to obtain zero after the 8/8 association gate |
| 7 | Observed/traceable SSP available | 0 | Administrative/evidence gate; absent in the existing run |
| 8 | Full waveform extraction passed | 0 | Signal-processing gate; absent for every candidate in the existing admission table |

### Key conclusion

The current `0 admitted` result is **not** explained by a lack of ideal local `(r,z)` sensitivity alone. The first decisive physical-robustness bottleneck, after continuity/numerics/resolvability/amplitude screening, is **cross-environment path-family association**. Therefore the next B3D-3 run should prioritize an observed/traceable SSP ensemble and family identity tracking, not invent additional pairwise delays.

Importantly, this conclusion does not depend on the exact 0.15 ms environment-drift threshold: after the first four gates, the three remaining candidates already fail the stricter requirement of being associated in all 8 diagnostic environments.

## Pre-environment shortlist

Three candidates survive the first four gates:

| key | reference | common range (km) | delay over window (ms) | nearest competing arrival, median (ms) | max numerical change (ms) | environment association | max environment delay shift (ms) | median rel. level (dB) | Jacobian cells |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| S0B1_s1r0 | S0B1_s0r1 | 52.0–54.25 | -3.3150 to -2.9641 | 1.5145 | 0.0114 | 4/8 | 0.8087 | -19.99 | 6 |
| S1B2_s1r0 | S0B1_s0r1 | 52.0–54.5 | 979.4655 to 1015.9187 | 47.9431 | 0.0038 | 7/8 | 5.3253 | -18.15 | 11 |
| S3B2_s0r1 | S0B1_s0r1 | 52.0–54.5 | 1074.7109 to 1119.5831 | 51.4488 | 0.0077 | 7/8 | 5.9128 | -18.91 | 11 |

None is admitted. The first is a millisecond-scale short relative delay but has only 4/8 environment associations. The latter two are approximately 1 s long-delay boundary-reflection combinations and have 7/8 association, but their environmental delay movement is several milliseconds and they still require an observed-SSP re-run and waveform extraction validation.

## Negative audit decisions

1. Do not count the 26 algebraic candidates as 26 independent anchors.
2. Do not count the two ~1 s long-delay candidates as two independent anchors until their sensitivity directions and path identities are shown to be non-redundant under the observed-SSP ensemble.
3. Do not weaken failed path labels by silently dropping the failing eighth environment.
4. Do not reinterpret the three-candidate shortlist as operational admission.
5. Do not use the old shallow 18N/110E WOA23 reference asset for this project; its ~102–110 m water depth is incompatible with ~200 m source/receiver depth.

## Next run target

The next mature propagation pass should start with the three shortlisted path combinations and their reference family, using a traceable deep-water WOA23/SCSPOD14/Argo SSP ensemble. Full-field arrivals must still be generated so that competing-path ambiguity can be detected, but admission analysis can focus on whether these three identities survive and whether any genuinely new path family adds a distinct singular direction.
