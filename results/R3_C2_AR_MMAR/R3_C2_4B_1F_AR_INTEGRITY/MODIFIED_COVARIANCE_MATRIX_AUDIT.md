# MODIFIED_COVARIANCE_MATRIX_AUDIT

UTC: 2026-09-24T10:18:32.049426+00:00

Forward (n=p..N-1): A=[y[n-1]…y[n-p]] → −y[n]
Backward (n=p..N-1): A=[y[n-p+1]…y[n]] → −y[n-p]  （与 forward **独立**）

2(N_r−p) at N_r=200,p=133: **134** 方程 / **133** 参数 → 仍偏紧（rank 见 csv）。

  p  N_r  n_forward  n_backward  n_equations  expected_2_Nr_minus_p  n_params  rank         cond  backward_dup_forward  sufficient_independent
  1  200        199         199          398                    398         1     1 1.000000e+00                 False                    True
  7  200        193         193          386                    386         7     2 1.321878e+14                 False                   False
133  200         67          67          134                    134       133     2 3.030928e+16                 False                   False

backward_dup_forward 必须为 False（修复后）。