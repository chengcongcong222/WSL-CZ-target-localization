# R3-B1 — GPT 同步稿

- UTC: 2026-09-21T02:08:15.582906+00:00
- **判定：B1_DELAY_BANDWIDTH_LIMITED**
- Delay observability weak/periodic across S0/S1/S2 under current band; delay_df=[{'source': 'S0', 'n_freqs': 23, 'effective_bandwidth_Hz': 235.00483554019206, 'amb_mainlobe_width_s': 0.2982, 'amb_sidelobe_max': 0.4964808349270524, 'crb_sigma_tau_s': 0.00023460321322941465, 'periodic_ambiguity_risk': False, 'delay_usable': False, 'note': 'S2 multi-line can create periodic delay ambiguity; do not assume best'}, {'source': 'S1', 'n_freqs': 5, 'effective_bandwidth_Hz': 186.31070822687568, 'amb_mainlobe_width_s': 0.39840000000000003, 'amb_sidelobe_max': 0.4996524926456676, 'crb_sigma_tau_s': 0.0002959190594404395, 'periodic_ambiguity_risk': False, 'delay_usable': False, 'note': 'S2 multi-line can create periodic delay ambiguity; do not assume best'}, {'source': 'S2', 'n_freqs': 5, 'effective_bandwidth_Hz': 209.9668545270896, 'amb_mainlobe_width_s': 0.3834, 'amb_sidelobe_max': 0.49985338502952864, 'crb_sigma_tau_s': 0.00026257901356076203, 'periodic_ambiguity_risk': False, 'delay_usable': False, 'note': 'S2 multi-line can create periodic delay ambiguity; do not assume best'}].
- 下一步：do not claim practical delay MMAC; consider RC3-C or broader bandwidth studies later

## Gate-0
 abs_phi_deg  delta_u_theta0  L_req_200Hz_m  L_req_300Hz_m
         1.0        0.000152   49243.345269   32828.896846
         2.0        0.000609   12311.773889    8207.849259
         5.0        0.003805    1970.934286    1313.956191
        10.0        0.015192     493.672859     329.115239
        15.0        0.034074     220.108052     146.738701
        20.0        0.060307     124.362891      82.908594

HLA_PROJECTION_WEAK=False

## 控制受控 CRB_r (km)
theta_rel_deg  obs_mode              obs_set     
0.0            HLA_DIRECTION_COSINE  A_angle_only             inf
                                     B_delay_only    4.346236e+12
                                     C_joint         7.813265e+06
               IDEAL_ELEVATION       A_angle_only    7.713585e+04
                                     B_delay_only    4.346236e+12
                                     C_joint         7.713585e+04
30.0           HLA_DIRECTION_COSINE  A_angle_only             inf
                                     B_delay_only    4.346236e+12
                                     C_joint         9.021940e+06
               IDEAL_ELEVATION       A_angle_only    7.713585e+04
                                     B_delay_only    4.346236e+12
                                     C_joint         7.713585e+04
60.0           HLA_DIRECTION_COSINE  A_angle_only             inf
                                     B_delay_only    4.346236e+12
                                     C_joint         1.562632e+07
               IDEAL_ELEVATION       A_angle_only    7.713585e+04
                                     B_delay_only    4.346236e+12
                                     C_joint         7.713585e+04

## 延迟源
source  n_freqs  effective_bandwidth_Hz  amb_mainlobe_width_s  amb_sidelobe_max  crb_sigma_tau_s  periodic_ambiguity_risk  delay_usable                                                                  note
    S0       23              235.004836                0.2982          0.496481         0.000235                    False         False S2 multi-line can create periodic delay ambiguity; do not assume best
    S1        5              186.310708                0.3984          0.499652         0.000296                    False         False S2 multi-line can create periodic delay ambiguity; do not assume best
    S2        5              209.966855                0.3834          0.499853         0.000263                    False         False S2 multi-line can create periodic delay ambiguity; do not assume best

停止：B1 后停止；不自动 B2/RC3-C/P5。
