# R3 FINAL PARAMETER-CONTRACTION AUDIT

UTC: 2026-09-30T10:28:22.783014+00:00

## Stage Contraction

                     stage  n_candidates  fraction_of_initial  r_width  theta_width  v_width  psi_width
           S0_INITIAL_GRID        114576             1.000000     15.0         10.0      2.0       30.0
            S1_RC2_W1_600S          1326             0.011573     15.0          0.0      2.0       13.0
     S2_RC2_RC3_W1_MAIN235           121             0.001056     15.0          0.0      0.2        8.0
S3_RC2_W1W2_STRAIGHT_1200S           495             0.004320     15.0          0.0      2.0       11.0
   S4_RC2_RC3_W1W2_MAIN235            43             0.000375     15.0          0.0      0.2        3.0
             S5_RC2_TURN15           385             0.003360     15.0          0.0      2.0        8.0
        S6_TURN15_PLUS_235            29             0.000253     15.0          0.0      0.2        2.0
     S7_TURN15_PLUS_TRIPLE             2             0.000017      0.0          0.0      0.0        1.0
       S8_TURN15_PLUS_FOUR             2             0.000017      0.0          0.0      0.0        1.0

## Final Survivors (S7/S8)

n=2, r=50, theta=0, v=2, psi=4-5°

## Truth-Relative Error

                stage  z_true_m  n_survivors  max_rel_err_r  max_abs_err_theta_deg  max_rel_err_v  max_rel_err_psi  worst_survivor_psi_deg
   S6_TURN15_PLUS_235     180.0            2           0.08                    0.0            0.0              0.0                     5.0
   S6_TURN15_PLUS_235     200.0           29           0.20                    0.0            0.1              0.2                     6.0
   S6_TURN15_PLUS_235     220.0           15           0.20                    0.0            0.0              0.2                     4.0
S7_TURN15_PLUS_TRIPLE     180.0            2           0.00                    0.0            0.0              0.2                     4.0
S7_TURN15_PLUS_TRIPLE     200.0            2           0.00                    0.0            0.0              0.2                     4.0
S7_TURN15_PLUS_TRIPLE     220.0            2           0.00                    0.0            0.0              0.2                     4.0
  S8_TURN15_PLUS_FOUR     180.0            2           0.00                    0.0            0.0              0.2                     4.0
  S8_TURN15_PLUS_FOUR     200.0            2           0.00                    0.0            0.0              0.2                     4.0
  S8_TURN15_PLUS_FOUR     220.0            2           0.00                    0.0            0.0              0.2                     4.0

## 结论

- 距离锚定：r=50 单格
- 四维候选：psi 仍有 1° bin 宽度
- universal_lt10pct_set_bound = NOT_ESTABLISHED
- DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION
