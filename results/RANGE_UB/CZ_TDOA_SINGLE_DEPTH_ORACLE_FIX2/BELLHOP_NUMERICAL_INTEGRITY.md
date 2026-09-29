# BELLHOP_NUMERICAL_INTEGRITY

UTC: 2026-09-29T04:52:02.586751+00:00

## BOX correction

Before: `0.0 r_max_m H_DEPTH` → Box%z=70000m, Box%r=5000km (WRONG)
After: `0.0 5001.0 70.0` → STEP=auto, ZBOX=5001.0m, RBOX=70.0km

## .prt audit (21 runs)

- All return_code=0: True
- Zero insufficient_storage: True
- ZBOX≈5001m: True

Status: PASS
