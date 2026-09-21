# B3D3 M4 Local Regression Pack — 2026-09-12

Purpose: verify that three historical online validation results can be reproduced or independently recomputed in the new permanent `Ubuntu-B3D3` environment.

This pack is NOT physical anchor evidence and does NOT freeze a project site.

## Regression targets

### R1 — Snell diagnostic preflight
Historical diagnostic-only custom Snell-integral result:
- E_D3900_C1529 upper return envelope ≈ 48.09 km
- E_D4300_C1529 upper return envelope ≈ 53.23 km
- depth 3900→4300 m shifts upper envelope by ≈ +5.14 km
- E_D4100_C1525 upper ≈ 52.72 km
- E_D4100_C1534 upper ≈ 48.81 km
- cdeep 1525→1534 m/s shifts upper envelope by ≈ -3.90 km

This is DIAGNOSTIC_NOT_ADMISSION. It was a custom range-independent Snell integral, not Bellhop/Kraken.

### R2 — P1 nominal local Jacobian
At nominal reconstructed environment, z=200 m, r=53 km:
- S1B2_s1r0:
  dτ/dr ≈ -14.7476 ms/km
  dτ/dz ≈ -0.19627 ms/m
  scaled [1 km,20 m] gradient ≈ [-14.7476,-3.9254] ms
- S3B2_s0r1:
  dτ/dr ≈ -18.1428 ms/km
  dτ/dz ≈ +0.11177 ms/m
  scaled ≈ [-18.1428,+2.2354] ms
- acute angle ≈ 21.929°
- joint unwhitened singular values ≈ [23.3928,4.45370] ms
- condition ≈ 5.25244

This supports nominal local nonredundancy only. It is not real-ocean anchor admission.

### R3 — candidate gate waterfall
Historical mature-model candidate gate sequence:
26 → 25 → 23 → 9 → 3 → 0

Three pre-environment survivors:
- S0B1_s1r0 (P2/control)
- S1B2_s1r0 (P1)
- S3B2_s0r1 (P1)

Final admission count remained 0 because environment/association/extractability gates were not fully passed.

## Evidence language
Use:
- REGRESSION_PASS
- REGRESSION_PARTIAL
- REGRESSION_FAIL

Do NOT use:
- PHYSICAL_ANCHOR_PASS
- REAL_OCEAN_VALIDATED
- PROJECT_SITE_VALIDATED

Project site remains TBD.
