# R3-C2.1 报告

UTC：2026-09-23T07:32:46.685640+00:00

- 求解器：**KRAKEN (AT Win10 预编译)**，退出自研 FEM
- Erratum δ 相位已知；本轮 **δ=0** paper-like
- Yang 公式：**YANG_FORMULA_RECOVERY_PARTIAL**

## 频率边界 (L=2.4 km)
 f_hz  mean_fwhm_2p4  mean_abs_err  mean_psl  mean_n_peaks
201.0      19.333333     25.333333  0.435237      0.666667
235.0       2.000000     33.333333  0.119551      0.333333
283.0       1.333333     33.333333  0.159012      0.333333
338.0       1.333333     33.333333  0.079130      0.333333

## 判定 `C2_1_YANG_APERTURE_LIMITED`

Mature solver still cannot form enough modal/depth peaks at L≤2.4km. [{'f_hz': 201.0, 'mean_fwhm_2p4': 19.333333333333332, 'mean_abs_err': 25.333333333333332, 'mean_psl': 0.4352366913681409, 'mean_n_peaks': 0.6666666666666666}, {'f_hz': 235.0, 'mean_fwhm_2p4': 2.0, 'mean_abs_err': 33.333333333333336, 'mean_psl': 0.11955145154280102, 'mean_n_peaks': 0.3333333333333333}, {'f_hz': 283.0, 'mean_fwhm_2p4': 1.3333333333333333, 'mean_abs_err': 33.333333333333336, 'mean_psl': 0.15901230843058112, 'mean_n_peaks': 0.3333333333333333}, {'f_hz': 338.0, 'mean_fwhm_2p4': 1.3333333333333333, 'mean_abs_err': 33.333333333333336, 'mean_psl': 0.0791304347825986, 'mean_n_peaks': 0.3333333333333333}] | YANG_FORMULA_RECOVERY_PARTIAL: SA mode-sum + modal peak depth score with δ=0; not a full paper estimator.

## snapshot vs SA
 f_hz    zs    L_m  z_hat  fwhm      psl  n_groups     kind           note
201.0 180.0 1200.0  180.0  48.0 0.999635       1.0      NaN            NaN
201.0 180.0 2400.0  180.0  48.0 0.999635       1.0      NaN            NaN
201.0 180.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
201.0 200.0 1200.0  206.0  10.0 0.306075       1.0      NaN            NaN
201.0 200.0 2400.0  206.0  10.0 0.306075       1.0      NaN            NaN
201.0 200.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
201.0 220.0 1200.0  150.0   0.0 0.000000       0.0      NaN            NaN
201.0 220.0 2400.0  150.0   0.0 0.000000       0.0      NaN            NaN
201.0 220.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
235.0 180.0 1200.0  150.0   0.0 0.000000       0.0      NaN            NaN
235.0 180.0 2400.0  150.0   0.0 0.000000       0.0      NaN            NaN
235.0 180.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
235.0 200.0 1200.0  200.0   6.0 0.358654       1.0      NaN            NaN
235.0 200.0 2400.0  200.0   6.0 0.358654       1.0      NaN            NaN
235.0 200.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
235.0 220.0 1200.0  150.0   0.0 0.000000       0.0      NaN            NaN
235.0 220.0 2400.0  150.0   0.0 0.000000       0.0      NaN            NaN
235.0 220.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
283.0 180.0 1200.0  150.0   0.0 0.000000       0.0      NaN            NaN
283.0 180.0 2400.0  150.0   0.0 0.000000       0.0      NaN            NaN
283.0 180.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
283.0 200.0 1200.0  200.0   4.0 0.477037       1.0      NaN            NaN
283.0 200.0 2400.0  200.0   4.0 0.477037       1.0      NaN            NaN
283.0 200.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
283.0 220.0 1200.0  150.0   0.0 0.000000       0.0      NaN            NaN
283.0 220.0 2400.0  150.0   0.0 0.000000       0.0      NaN            NaN
283.0 220.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
338.0 180.0 1200.0  150.0   0.0 0.000000       0.0      NaN            NaN
338.0 180.0 2400.0  150.0   0.0 0.000000       0.0      NaN            NaN
338.0 180.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
338.0 200.0 1200.0  200.0   6.0 0.434765       1.0      NaN            NaN
338.0 200.0 2400.0  200.0   4.0 0.237391       1.0      NaN            NaN
338.0 200.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture
338.0 220.0 1200.0  150.0   0.0 0.000000       0.0      NaN            NaN
338.0 220.0 2400.0  150.0   0.0 0.000000       0.0      NaN            NaN
338.0 220.0    0.0    NaN   NaN      NaN       NaN snapshot no SA aperture

停止：无 δ 误差 / Doppler / f0 / MC / P5。