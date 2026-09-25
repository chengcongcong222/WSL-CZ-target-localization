# R3-RC3-REANCHOR-2 depth-profiled discrimination

UTC: 2026-09-25T08:56:26.990583+00:00

## 冻结角色

科学 `DEPTH_ROLE_MIXED_OR_UNRESOLVED` · 系统 `DEPTH_PROFILED_NUISANCE_VARIABLE` · 主状态 [r,θ,v,ψ]

## 方法

30 对困难样本；y = SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE；
J_alt_fixed = J(h_alt, z_true)；J_alt_profiled = min_z J(h_alt, z)。

阈值（预注册）：accept 若 J≤{0.5,1.0,1.5} dB RMS，主阈值 **1.0** dB。

## 结果（thr=1.0）

| | fixed z | **profiled z** | Δ |
| --- | ---: | ---: | ---: |
| 错误候选排除率 | 0.967 | **0.967** | +0.000 |

标签：`PROFILED_PROPAGATION_REJECTS_MAJORITY_OF_HARD_ALTS`

## 机制

mechanism  n_cases  median_J_alt_fixed  median_J_alt_profiled  rej_fixed_z_thr0.5  rej_profiled_thr0.5  incremental_rej_thr0.5  rej_fixed_z_thr1.0  rej_profiled_thr1.0  incremental_rej_thr1.0  rej_fixed_z_thr1.5  rej_profiled_thr1.5  incremental_rej_thr1.5
        A       30            4.448839               3.542300                 0.9                  0.9                     0.0                 0.9                  0.9                     0.0                 0.9                  0.8                    -0.1
        B       30            3.554176               2.629880                 1.0                  1.0                     0.0                 1.0                  1.0                     0.0                 1.0                  0.8                    -0.2
        C       30            3.382595               2.597404                 1.0                  1.0                     0.0                 1.0                  1.0                     0.0                 1.0                  0.8                    -0.2

## 边界

- 无深度 RMSE；z* 仅补偿诊断
- 不用旧 P3 97%（`HISTORICAL_P3_SIMPLIFIED_RESULT`）
- 非 Liang D(z)

## 结论方向

在允许 z∈[150,250] profile 后，传播仍对错误水平轨迹有**可观排除力**（见 rej_profiled）。
若增量相对 fixed-z 有限，说明深度消元是正确建模而非“白送信息”。
