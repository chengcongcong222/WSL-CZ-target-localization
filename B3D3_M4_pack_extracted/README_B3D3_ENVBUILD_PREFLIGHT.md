# B3D-3 deep-environment construction and preflight — 2026-09-09

## Status

**Purpose:** preflight the sensitivity of first-CZ-like deep-turning propagation geometry before the two P1 path families are rerun with a mature propagation solver.

**Evidence level:** `DIAGNOSTIC_NOT_ADMISSION`.

This directory does **not** contain a new Bellhop/Kraken/RAM run and does **not** close the real-anchor admission gate.

## Evidence anchors used

Upper-ocean checkpoints are retained from the project's previous WOA23 extraction at 18°N, 114°E, April climatology:

- c(0 m) = 1538.81 m/s
- c(100 m) = 1523.82 m/s
- c(700 m) = 1486.84 m/s
- c(1000 m) = 1483.92 m/s

An independent published South-China-Sea summer observed mean SSP reports approximately:
- minimum c ≈ 1484 m/s around 1180 m,
- deep/bottom c ≈ 1529.4 m/s.

Independent South-China-Sea acoustic experiments report deep-water settings around ~3910 m and mean ~4305 m. These are used only to define a **regional diagnostic bathymetry bracket**, not as the exact bathymetry at 18°N,114°E.

## Diagnostic construction

The 0–1000 m WOA-derived checkpoints are joined with shape-preserving PCHIP interpolation. The ~1180 m and deep-speed observations are cross-study constraints. The following stress cases are created:

- E_D3900_C1529
- E_D4100_C1529
- E_D4300_C1529
- E_D4100_C1525
- E_D4100_C1534

Intermediate points are explicitly marked `DERIVED_DIAGNOSTIC_INTERPOLATION`.

## Custom ray-integral preflight

A range-independent Snell-law integral is used for a source and receiver at 200 m, sweeping downward grazing angles 0.5–15°. It computes the first deep-turning return to 200 m.

This is deliberately labeled `CUSTOM_DIAGNOSTIC_RAY_INTEGRAL_NOT_BELLHOP`.

Key return-range envelopes:

| Environment | Return-range envelope (km) | Angles landing in 50–60 km |
|---|---:|---:|
| E_D3900_C1529 | 44.08–48.09 | 0 |
| E_D4100_C1529 | 46.66–50.66 | 11 |
| E_D4300_C1529 | 49.24–53.23 | 58 |
| E_D4100_C1525 | 49.28–52.72 | 48 |
| E_D4100_C1534 | 44.29–48.81 | 0 |

### Falsification-relevant result

Holding the deep endpoint speed at 1529.4 m/s while changing diagnostic water depth from 3900 m to 4300 m moves the upper edge of the first deep-turning return envelope from
**48.09 km to 53.23 km**, a shift of
**5.14 km**.

At 4100 m water depth, changing only the diagnostic deep endpoint speed from 1525 to 1534 m/s moves the envelope maximum from
**52.72 km to 48.81 km**, a shift of
**3.90 km**.

Therefore, the existing P1 failure mode (7/8 path association) is physically credible: deep structure and bathymetry can move/terminate deep-turning families on kilometre scales. The preflight does **not** identify `S1B2_s1r0` or `S3B2_s0r1`; only the mature solver with topology/eigenray association can do that.

## Next admissible computation

Only these two P1 families are to be rerun:

1. `S1B2_s1r0`
2. `S3B2_s0r1`

The run must use:
- full-depth traceable SCS SSP(s),
- site/track bathymetry with provenance,
- mature Bellhop/BellhopCXX (and ideally a second solver cross-check),
- topology/eigenray continuity rather than symbol-name matching,
- then full-channel broadband extraction,
- then admitted-only Jacobian/SVD and Core/Guard revalidation.

## Public sources

- ODB WOA23 API source and usage: https://github.com/cywhale/woa23
- NOAA WOA23: https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/
- Gao, Xu & Li (2022), Chinese Physics B, DOI 10.1088/1674-1056/ac6014
- SCS deep-water acoustic experiment / CTD context: https://pmc.ncbi.nlm.nih.gov/articles/PMC5856022/
- IOA/CAS SCS deep-water experiment context: https://english.ioa.cas.cn/as/201512/t20151230_158285.html

## Conclusion

This preflight **does not upgrade B3D-3 from HOLD**.

It does justify a stricter next step: full-depth SSP and bathymetry are not optional nuisance inputs for the P1 families. They are likely first-order determinants of whether those families exist and remain associable in the 50–60 km window.
