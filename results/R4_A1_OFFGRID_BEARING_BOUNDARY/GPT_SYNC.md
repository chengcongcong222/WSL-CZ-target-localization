# R4-A1 structural stop

**R4_A1_BLOCKED_BY_COARSE_RC2_GRID_SELECTION_AND_OFFGRID_MODAL_MISMATCH**

9 jointly off-grid truths; full RC2-only axis:1629 cases,30 seeds per positive sigma/panel.
End-to-end:27 nominal0.1deg pilot cases only; full MC paused.
RC2 bracketing-cell retention:0/270 at0.02deg,0/270 at0.05deg;0/9 noiseless.
Pilot top1 coarse-cell success:0/27. R3 on-grid replay passes.

       metric  top1_P50  top1_P95  survivor_P50  survivor_P95
        rel_r  0.075596  0.160093      0.175975      0.208946
abs_theta_deg  0.180000  0.230000      0.180000      0.230000
        rel_v  0.069519  0.273885      0.234568      0.734104
  abs_psi_deg  5.650000  8.650000      6.580000      8.650000

Interpretation MIXED; no intrinsic continuous-identifiability conclusion.
A1-3 incomplete; NO_STABLE_REGION_ESTABLISHED. R4 progress0%.
Recommend A1-FIX quantization/continuation Gate before further end-to-end MC. R3 frozen; P5 not opened.
