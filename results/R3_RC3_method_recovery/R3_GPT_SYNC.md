# R3-0 — GPT 同步稿

- UTC: 2026-09-20T07:03:13.529234+00:00
- 重锚定：RC3 文献方法正式验证；P4.5≠RC3上限；不进入P5

## 冻结主线

- **上限基准**：宽带匹配场MFP → r,z
- **RC3-A距离主线**：会聚区宽带干涉条纹/波导不变量 → r,Δr
- **RC3-B距离–深度主线**：多途角+时延联合MMAC → r,z
- **RC3-C深度主线**：HLA模态/波数谱+合成孔径深度 → z

## MFP 上限（摘录）

source truth_tag  r_true_km  z_true_m  r_hat_km  z_hat_m  err_r_km  err_z_m         J_min  r_mainlobe_km  z_mainlobe_m  n_local_minima_r  truth_is_global_opt  corr_rz_basin  n_freqs
    S0   ON_GRID       50.0     200.0      50.0    200.0       0.0      0.0 -2.220446e-16           7.75          92.5                 1                 True      -0.563123       24
    S0  OFF_GRID       52.3     213.0      52.5    210.0       0.2     -3.0  2.361480e-05           7.50          92.5                 1                False      -0.645085       24
    S0 OFF_GRID2       47.6     187.0      47.5    187.5      -0.1      0.5  2.291616e-05           6.75          87.5                 1                 True      -0.428992       24
    S1   ON_GRID       50.0     200.0      50.0    200.0       0.0      0.0  0.000000e+00           6.25          82.5                 1                 True      -0.644342        5
    S1  OFF_GRID       52.3     213.0      52.5    210.0       0.2     -3.0  2.841852e-05           6.00          90.0                 1                False      -0.714341        5
    S1 OFF_GRID2       47.6     187.0      47.5    187.5      -0.1      0.5  2.198347e-05           6.00          85.0                 1                 True      -0.436177        5
    S2   ON_GRID       50.0     200.0      50.0    200.0       0.0      0.0  0.000000e+00           6.25          85.0                 1                 True      -0.626083        5
    S2  OFF_GRID       52.3     213.0      52.5    210.0       0.2     -3.0  2.796713e-05           7.75          92.5                 1                False      -0.630346        5
    S2 OFF_GRID2       47.6     187.0      47.5    187.5      -0.1      0.5  2.062429e-05           6.00          87.5                 1                 True      -0.409382        5

## RC3-A 条纹复现（摘录）

            case                                                           note  r_true_km  z_true_m  beta_hat  beta_score  r_hat_km  err_r_km  fringe_align_width_km  n_freqs                                                                                        method
PAPER_LIKE_IDEAL E-STD modal field, matched env, noiseless — algorithm fidelity       50.0     200.0       2.5    1.132150     54.25      4.25                   0.00       46 WI-like fringe alignment after range-frequency de-trend (paper-like beta search + range warp)
         OFFGRID                                                 off-grid truth       52.3     213.0       2.5    1.514182     51.50     -0.80                   0.25       46 WI-like fringe alignment after range-frequency de-trend (paper-like beta search + range warp)
        S2_LINES                          fringes sampled on S2 line freqs only       50.0     200.0       2.5    0.706468     45.00     -5.00                   2.50        5 WI-like fringe alignment after range-frequency de-trend (paper-like beta search + range warp)

## 停止

本轮仅 R3-0 + RC3-A 首轮复现启动；**不**自动 R3-B/C，**不** P5。
