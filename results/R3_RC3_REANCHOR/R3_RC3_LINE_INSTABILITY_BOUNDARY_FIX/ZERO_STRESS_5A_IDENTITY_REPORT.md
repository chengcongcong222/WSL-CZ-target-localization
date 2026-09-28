# ZERO_STRESS_5A_IDENTITY_REPORT

UTC: 2026-09-28T07:34:51.584019+00:00

## Gate

`ZERO_STRESS_5A_IDENTITY_GATE` = **PASS**

A_rms=0, D_f=0 时 5B-FIX 与 5A raw `SUBSET_PROFILED_MARGIN_RESULTS.csv` 逐行比较
(pair_id × z_true × config)，量：J_true* / J_alt* / ΔJ。

浮点容差：1e-9 dB。全局最大绝对误差 = 7.105e-15 dB。

## 按 config

config  n_rows  max_err_Jt   max_err_Ja   max_err_dJ  all_pass
   201      90         0.0 2.442491e-15 2.442491e-15      True
   235      90         0.0 3.330669e-15 3.330669e-15      True
   283      90         0.0 5.995204e-15 5.995204e-15      True
   338      90         0.0 3.108624e-15 3.108624e-15      True
  FOUR      90         0.0 7.105427e-15 7.105427e-15      True

## 结论

零扰动严格退化到 5A。observable = 轨迹窗去均值 `SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE`。
