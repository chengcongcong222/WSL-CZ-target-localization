# R3-CLOSEDLOOP-1A — 单窗口 RC2→RC3 完整候选集收缩

UTC: 2026-09-28T08:40:16.269226+00:00

基线：`3dd1f350b48110744ebe49447267a0e5bc26194f`

## 判定

### `SET_LEVEL_RC3_INCREMENT_CONFIRMED`

MAIN tau=0.5: true_kept=True count>=50%=True axes>=2/3 @>=20%=True; counts=[np.int64(122), np.int64(121), np.int64(137)]/1326; r_contr=[np.float64(0.0), np.float64(0.0), np.float64(0.0)]; v_contr=[np.float64(0.8), np.float64(0.9), np.float64(0.6)]; psi_contr=[np.float64(0.385), np.float64(0.385), np.float64(0.385)]

## 设定

| 项 | 值 |
|---|---|
| 场景 | CR5：r0=50 km, θ0=0°, v=2 m/s, ψ=5° |
| 窗口 | T=600 s, σθ=0.1°, turn=0°, dt=10 s (n_obs=61) |
| 种子 | 20260912 · `SINGLE_WINDOW_DEMO_REALIZATION` |
| 搜索盒 | 16×21×11×31 = **114576** 节点 |
| RC2 接受 | cost ≤ cmin + 13.3 σ² = 4.0514e-05 |
| RC3 | DIRECT_F_MATRIX_ON_TRAJECTORY，z profile 150:5:250 m |
| 声学 | MAIN=235 Hz；REFERENCE=FOUR；E0 无扰动 |

## RC2 候选云

| 指标 | 值 |
|---|---|
| n_RC2 | **1326** / 114576 |
| r 宽度 | 15.0 km |
| θ0 宽度 | 0.00° |
| v 宽度 | 2.00 m/s |
| ψ 宽度 | 13.00° |
| 真值在网格 | True |
| 真值在 RC2 | True |

## MAIN=235 Hz, τ=0.5 dB 集合收缩

config  z_true_m  tau_db  n_RC2  n_RC3  count_contraction  r_width_km  r_contraction  theta_width_deg  theta_contraction  v_width_mps  v_contraction  psi_width_deg  psi_contraction  true_retained  true_delta_J
  MAIN     180.0     0.5   1326    122           0.907994        15.0            0.0              0.0                0.0          0.4            0.8            8.0         0.384615           True           0.0
  MAIN     200.0     0.5   1326    121           0.908748        15.0            0.0              0.0                0.0          0.2            0.9            8.0         0.384615           True           0.0
  MAIN     220.0     0.5   1326    137           0.896682        15.0            0.0              0.0                0.0          0.8            0.6            8.0         0.384615           True           0.0

## REFERENCE=FOUR, τ=0.5 dB

   config  z_true_m  tau_db  n_RC2  n_RC3  count_contraction  r_width_km  r_contraction  theta_width_deg  theta_contraction  v_width_mps  v_contraction  psi_width_deg  psi_contraction  true_retained  true_delta_J
REFERENCE     180.0     0.5   1326     99           0.925339        15.0            0.0              0.0                0.0          0.0            1.0            7.0         0.461538           True           0.0
REFERENCE     200.0     0.5   1326     94           0.929110        15.0            0.0              0.0                0.0          0.0            1.0            7.0         0.461538           True           0.0
REFERENCE     220.0     0.5   1326    102           0.923077        15.0            0.0              0.0                0.0          0.0            1.0            7.0         0.461538           True           0.0

## 排序诊断（MAIN, z_true=200）

config  z_true_m  true_rank  true_J  true_z_star  best_false_J    J_p10  J_median    J_p90  J_min  n_acc local_label  local_r0_km  local_theta0_deg  local_v  local_psi_deg  local_J local_is_true
  MAIN     200.0          1     0.0        200.0      0.008116 0.575312  1.863003 3.797663    0.0   1326         NaN          NaN               NaN      NaN            NaN      NaN           NaN
  MAIN     200.0         -1     NaN          NaN      0.000000      NaN       NaN      NaN    0.0   1326     LOCAL_1         50.0               0.0      2.0            5.0 0.000000          True
  MAIN     200.0         -1     NaN          NaN      0.008116      NaN       NaN      NaN    0.0   1326     LOCAL_2         58.0               0.0      2.0            4.0 0.008116         False
  MAIN     200.0         -1     NaN          NaN      0.009051      NaN       NaN      NaN    0.0   1326     LOCAL_3         56.0               0.0      2.0            6.0 0.009051         False

## 完整性

见 `INTEGRITY_GATES.md`。零信息控制、真值不注入、RC2 与 P2 一致。

## 关键结构发现

- **真值 rank=1**（MAIN 与 REFERENCE，三个 z_true 均是），J_true*=0 为自匹配。
- **候选数收缩 ~90%**：1326 → 121（MAIN τ=0.5）。集合级增量成立。
- **v 轴收缩最强**：τ=0.5 时宽度缩 60–90%；τ=0.25 时塌缩到单值 v=2.0。
- **ψ 轴中等收缩**：约 38%（13° → 8°）。
- **r 轴不收缩**：幸存者仍横跨 45–60 km。第二局部极小在 r=58 km（J≈0.008），第三在 r=56 km（J≈0.009），与 J_min=0 极为接近。单窗口会聚区 TL 形状在不同距离上存在近简并，这是**结构性发现**，不是实现问题。
- 幸存机制：12 个 near-truth cluster，其余 109 个为跨距离近简并（归入 other）。

## 标签

- `SINGLE_WINDOW_DEMO_REALIZATION`（非概率性能）
- `RC3_SET_THRESHOLD_SENSITIVITY`（预冻结 τ，非设备标定阈值）
- `WEAKEST_TESTED_SINGLE_LINE_MAIN` / `FOUR_LINE_REFERENCE`
- `DIRECT_F_MATRIX_ON_TRAJECTORY`

## 未做

多窗口递推、S1 迁移、line dropout、combined corner、SSP、幅频漂移、Liang、P5。
