# R3-B1.1 — GPT 同步稿

- UTC: 2026-09-21T03:26:14.655763+00:00
- 旧判定撤销：B1_DELAY_BANDWIDTH_LIMITED
- **新判定：B1_DELAY_AMBIGUITY_LIMITED**
- Local delay mainlobe ms-scale OK (S0_DENSE 0.00520000000000001 s) but global ambiguity remains (S0_23 comb side=0.9999551387640602, S0_DENSE side=0.4907037600350636). CRB vs ambiguity must stay separated.
- 下一步：do not kill MMAC; resolve global delay ambiguity or restrict to local peak tracking

## 时延修正
    source  n_freqs  freq_spacing_Hz  T_amb_period_s  local_mainlobe_width_s  first_null_width_s  global_max_sidelobe  crb_sigma_tau_s  effective_bw_Hz                                 represents
  S0_DENSE      226         1.000000        1.000000                  0.0052              0.0082             0.490704         0.000244       225.997788                       continuous broadband
S0_23POINT       23        10.227273        0.097778                  0.0050              0.0080             0.999955         0.000235       235.004836 discrete comb control (NOT true broadband)
        S1        5        36.000000        0.027778                  0.0064              0.0982             0.792748         0.000296       186.310708                               sparse lines
        S2        5        35.000000        0.028571                  0.0056              0.0178             0.981494         0.000263       209.966855                               sparse lines

S0_DENSE local mainlobe=0.00520000000000001 s; S0_23 comb side=0.9999551387640602; S0_DENSE side=0.4907037600350636

## 分支
n_branch=12 stable=4 max_resid=0.09861705982609692 m

## HLA / Fisher
{"delay_local_ok": true, "delay_global_ok_S0_DENSE": false, "delay_comb_ambiguity_S0_23": true, "hla_resolved_corr_14m": 0.0, "hla_resolved_corr_70m": 0.43333333333333335, "HLA_ELEVATION_PROJECTION_WEAK": true, "HLA_projection_note": "Corrected gate uses actual E-STD branch Delta-u + steering corr; prior False flag from composite CRB approximation is NOT frozen.", "estd_ideal_crb_r_m": NaN, "estd_hla_crb_r_m": NaN, "estd_delay_crb_r_m": NaN, "estd_ideal_ok": false, "estd_hla_ok": false, "estd_delay_fish_ok": false}

停止：无 B2/RC3-C/P5。
