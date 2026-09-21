# R3-B1-FINAL-CHECK — GPT 同步

- **判定：B1_NO_STABLE_MULTIPATH_IDENTITY**
- B1: PERMANENTLY_CLOSED
- Only 1 continuation-continuous AND modal-energy-supported branches at center.
- 下一步：B1 permanently closed; RC3-B MMAC not established

M_observable=1 delays=0 delay2d_allowed=False
asserts_pass=False nesting_ok=True
S0_usable=True S1/S2_amb=True

## Fisher (σθ=0.1°)
                 case  sigma_theta_deg  theta_rel_deg  M_observable  n_obs  rank  forbidden  crb_r_m  crb_z_m  corr_rz  cond_F  sv_min
IDEAL_ELEVATION+delay              0.1            0.0             1      1     1       True      NaN      NaN      NaN     NaN     NaN
          HLA_u+delay              0.1            0.0             1      1     1       True      NaN      NaN      NaN     NaN     NaN
           delay-only              0.1            0.0             1      0     0       True      NaN      NaN      NaN     NaN     NaN
IDEAL_ELEVATION+delay              0.1           30.0             1      1     1       True      NaN      NaN      NaN     NaN     NaN
          HLA_u+delay              0.1           30.0             1      1     1       True      NaN      NaN      NaN     NaN     NaN
           delay-only              0.1           30.0             1      0     0       True      NaN      NaN      NaN     NaN     NaN
IDEAL_ELEVATION+delay              0.1           60.0             1      1     1       True      NaN      NaN      NaN     NaN     NaN
          HLA_u+delay              0.1           60.0             1      1     1       True      NaN      NaN      NaN     NaN     NaN
           delay-only              0.1           60.0             1      0     0       True      NaN      NaN      NaN     NaN     NaN

## Branches
branch_id       topology   phi_deg     tau_s  continuation_5pt  modal_energy_supported        status_OBS
       C0      refracted -0.303597 33.330537              True                    True OBSERVABLE_BRANCH
       C1 surface_bounce 13.147062 33.813456             False                   False    NOT_OBSERVABLE

关闭 B1；无 B2/RC3-C/P5。
