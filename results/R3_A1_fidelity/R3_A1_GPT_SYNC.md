# R3-A1 — GPT 同步稿

- UTC: 2026-09-20T07:26:06.389540+00:00
- **判定：A1-NOT-DIRECTLY-TRANSFERABLE**
- Paper repro passed but E-STD beta_eff NON_CONSTANT (median=0.8235755194946945, IQR=0.8988305024903624, median_by_f_std=0.6855537631860799). Classical constant-beta fringe ranging does not directly transfer to deep first-CZ. | HLA-MFP: single: r宽=7.75km z宽=92.5m; hla_8x2m: r宽=15.00km z宽=100.0m; hla_diag_8x10m: r宽=15.00km z宽=100.0m
- 下一步：STOP RC3-A parameter tuning; proceed to RC3-B (MMAC) next round

## PE-WI 复现
             scenario  beta_true_mean  beta_true_std  beta_local_slope  beta_local_slope_std  beta_align_search  r_true_km  r_hat_peakmatch_km  err_r_peakmatch_km  r_hat_fringe_count_km  err_r_fringe_count_km  n_fringe_crossings  r_hat_powerlaw_km  err_r_powerlaw_km  range_error_km  fringe_score_width_km  fringe_score_width_powerlaw_km  pass_beta  pass_range                                                                         note
PE_WI_PAPER_CONDITION        1.001281        0.00072          0.959415              0.843322                2.2       25.0               25.02                0.02              26.256507               1.256507                  10              19.84              -5.16            0.02               0.493333                        0.493333       True        True beta_true uses fringe-slope convention; score_width is NOT localization RMSE

## β_eff (E-STD)
{"status": "NON_CONSTANT", "median": 0.8235755194946945, "iqr": 0.8988305024903624, "p5": -1.8677970185272932, "p95": 5.036863866285385, "std_all": 2.7215931135472693, "median_by_f_std": 0.6855537631860799, "local_std_mid_f": 0.33860232832294834, "frac_abs_lt_1p5": 0.7311715481171548, "frac_abs_gt_3": 0.13598326359832635}

## HLA-MFP
        config  n_elements  spacing_m  aperture_m                        note  r_hat_km  z_hat_m  err_r_km  err_z_m  r_mainlobe_km  z_mainlobe_m   corr_rz  truth_is_best        J_min
        single           1        0.0         0.0                单通道/当前R3-0参考     50.00    200.0      0.00      0.0           7.75          92.5 -0.569204           True 1.110223e-16
      hla_8x2m           8        2.0        14.0          基准水平阵 8×2m（孔径约14m）     45.00    167.5     -5.00    -32.5          15.00         100.0  0.000000          False 1.555147e-08
hla_diag_8x10m           8       10.0        70.0 大孔径诊断参考 8×10m（孔径约70m，非装备要求）     50.25    250.0      0.25     50.0          15.00         100.0  0.000000          False 4.334881e-08

## R3-0 口径
- MFP=冻结观测参考，非完整HLA物理上限
- RC3-A骨架≠论文复现完成
- score width ≠ RMSE

停止：无 hard-pair MC（除非PASS），无 R3-B/C，无 P5。
