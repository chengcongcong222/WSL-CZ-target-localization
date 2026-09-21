# MiMo Stage M4 — historical validation regression

M1/M2/M3 have passed. Execute M4 only.

Goal: independently reproduce/check three historical validation results in the permanent `Ubuntu-B3D3` environment. This is a regression/migration audit, NOT new anchor admission.

Do not:
- freeze a project site
- run formal P1 physical re-identification
- use 18N114E as project site
- change old distros
- replace historical source files with newly invented values
- call the Snell diagnostic “Bellhop”

## Inputs

Use the provided `B3D3_M4_regression_pack_20260912.zip`.
Copy/unzip it under the new project, for example:
`~/projects/cz-target-localization/reference/M4_legacy_regression/`

Preserve original files read-only if convenient.

Create new code/results separately:
- `scripts/regression/`
- `results/M4_legacy_regression/`

## R1 — Snell diagnostic regression

Read:
- `B3D3_SSP_DIAGNOSTIC_ENSEMBLE.csv`
- `B3D3_RAY_RETURN_PREFLIGHT.csv`
- `B3D3_RAY_RETURN_SUMMARY.csv`
- `README_B3D3_ENVBUILD_PREFLIGHT.md`

First audit the files and understand what was actually computed.

Then independently implement or reconstruct the range-independent Snell-invariant quadrature used for the diagnostic deep-return preflight. Do not just copy the reference output into a new CSV.

Recompute, as closely as the retained inputs permit:
- E_D3900_C1529
- E_D4100_C1529
- E_D4300_C1529
- E_D4100_C1525
- E_D4100_C1534

Historical target summaries:
- D3900/C1529 upper ≈ 48.09 km
- D4100/C1529 upper ≈ 50.66 km
- D4300/C1529 upper ≈ 53.23 km
- D4100/C1525 upper ≈ 52.72 km
- D4100/C1534 upper ≈ 48.81 km

Important:
This is `DIAGNOSTIC_NOT_ADMISSION`.
Never call it Bellhop.

Compare:
- return-range envelopes
- number of valid samples
- depth sensitivity trend
- deep-c sensitivity trend

If exact algorithm details are missing, distinguish:
- exact reproduction
- independent approximate reconstruction
Do not force equality.

## R2 — P1 nominal Jacobian regression

Read:
- `B3D3_P1_LOCAL_JACOBIAN_20260911.csv`
- `B3D3_P1_LOCAL_JACOBIAN_AUDIT_20260911.md`

Audit whether the CSV contains sufficient raw neighboring delay values to recompute the derivatives independently.

If sufficient, recompute from raw values rather than copying listed metrics.

Target:
S1B2_s1r0:
- dτ/dr ≈ -14.7476 ms/km
- dτ/dz ≈ -0.19627 ms/m
- scaled gradient [1 km,20m] ≈ [-14.7476,-3.9254]

S3B2_s0r1:
- dτ/dr ≈ -18.1428 ms/km
- dτ/dz ≈ +0.11177 ms/m
- scaled ≈ [-18.1428,+2.2354]

Joint:
- acute gradient angle ≈ 21.929°
- determinant ≈ -104.1845 ms²
- singular values ≈ [23.3928,4.45370] ms
- condition ≈ 5.25244

Use explicit formulas and save them in code/report.

Do NOT reproduce the old conditional whitening as a required regression unless exact rho/noise assumptions are present in the retained files. If not fully specified, keep it out of PASS criteria.

## R3 — gate waterfall regression

Read:
- `B3D3_GATE_WATERFALL.csv`
- `B3D3_PRE_ENVIRONMENT_SHORTLIST.csv`
- `B3D3_PRE_ENVIRONMENT_SHORTLIST_ENRICHED.csv`
- `B3D3_GATE_BOTTLENECK_AUDIT.md`
- `B3D3_P1_ANCHOR_DECISION.csv`
- `B3D3_STAGE_STATUS.json`

Independently audit row counts and gate logic from the retained tables.

Confirm or reject:
26 → 25 → 23 → 9 → 3 → 0

Confirm identities of the three pre-environment survivors:
- S0B1_s1r0
- S1B2_s1r0
- S3B2_s0r1

Check whether the final zero is genuinely caused by later environment/association/admission gates rather than an accidental filtering/counting bug.

Do not invent candidate→specific failed-environment mapping if the old exports do not contain it.

## Regression policy

For each R1/R2/R3 output one of:
- REGRESSION_PASS
- REGRESSION_PARTIAL
- REGRESSION_FAIL

R1 can PASS if the retained data/algorithm are independently reconstructed and key envelopes/trends match within a clearly justified numerical tolerance.
If algorithm details are incomplete but trends match, use PARTIAL.

R2 should normally be strict because it is local arithmetic from retained data.

R3 should normally be strict for counts/IDs if the retained tables are complete.

## Required outputs

Create:
- `M4_REPORT.md`
- `M4_STATUS.json`
- `R1_snell_regression.csv`
- `R1_snell_regression_summary.json`
- `R2_p1_jacobian_recomputed.csv`
- `R2_p1_jacobian_metrics.json`
- `R3_gate_waterfall_recomputed.csv`
- `R3_gate_audit.json`
- `SHA256SUMS_M4.txt`

Also preserve code/scripts used to recompute these results.

## Final report

1. Files audited
2. R1 method and result
3. R2 formulas/raw source/result
4. R3 filtering/count result
5. Discrepancies
6. Any legacy evidence that should be downgraded
7. Final regression verdict
8. Remaining blockers before M5

Machine-readable JSON:
{
  "stage": "M4",
  "project_site": "TBD",
  "R1_snell": "REGRESSION_PASS|REGRESSION_PARTIAL|REGRESSION_FAIL",
  "R2_p1_jacobian": "REGRESSION_PASS|REGRESSION_PARTIAL|REGRESSION_FAIL",
  "R3_gate_waterfall": "REGRESSION_PASS|REGRESSION_PARTIAL|REGRESSION_FAIL",
  "p1_nominal_nonredundancy_supported": true/false,
  "waterfall_26_25_23_9_3_0_supported": true/false,
  "legacy_candidate_specific_env_failure_mapping_available": false,
  "old_distros_modified": false,
  "M4_status": "PASS|CONDITIONAL_PASS|HOLD",
  "blockers": [],
  "next_recommended_action": "..."
}

If a historical result cannot be regenerated from retained evidence, report that honestly; do not infer missing raw data.
