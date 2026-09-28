# INTEGRITY_GATES

UTC: 2026-09-28T08:40:16.269226+00:00

## 1. RC2 cloud matches P2 definition — PASS
- coarse grid N=114576 (expect 114576)
- truth on grid: True
- truth in RC2 accepted: True
- cost = sum wrap², accept = cost ≤ cmin + 4.051412e-05

## 2. True grid node not artificially injected — PASS
- true node located by nearest-grid match to frozen truth (50 km, 0°, 2 m/s, 5°)
- estimator never receives truth

## 3. Acoustic observation only from truth — PASS
- L_obs generated from truth trajectory at z_true
- RC2 generation uses bearing only; acoustic enters only in RC3 scoring

## 4. RC3 zero-information control — PASS
- if J(h)≡const, then {"J≤Jmin+τ"} keeps all → RC2+RC3 = RC2
- survivors would be 1326 = n_RC2

## 5. MAIN / REFERENCE share same RC2 cloud — PASS
- single RC2 run; both acoustic configs scored on identical accepted indices

## Summary

All integrity gates PASS.
