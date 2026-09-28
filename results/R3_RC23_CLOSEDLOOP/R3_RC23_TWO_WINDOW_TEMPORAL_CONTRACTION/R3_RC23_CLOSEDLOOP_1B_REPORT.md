# R3-CLOSEDLOOP-1B — 双窗口时间累积与距离脊线审计

UTC: 2026-09-28T09:30:36.300105+00:00

基线：`8d5b9b7e74a97639daf5351050b01d4d5b868bf0`

## 判定

### `MULTIWINDOW_RC3_SET_INCREMENT_RANGE_UNRESOLVED`

MAIN tau=0.5: true_kept=True count>=50%=True r>=20%=False; cumRC2 r_width=15.0km (W1=15.0); counts=[np.int64(41), np.int64(43), np.int64(41)]/495; r_contr=[np.float64(0.0), np.float64(0.0), np.float64(0.0)]; RC2_time_r_contr=0.000

## 设定

| 项 | 值 |
|---|---|
| 场景 | CR5 匀速直线 t=0–1200 s |
| W1 | 0:10:600 s（严格复现 1A） |
| W2 | 610:10:1200 s |
| RC2 累计 | cost_W1+cost_W2，全 114576 节点重算 |
| RC3 | 分窗去均值 + 共享 z profile，J₁₂ 联合残差 |
| 声学 | MAIN=235 Hz，REFERENCE=FOUR，E0 |

## W1 Identity

W1 RC2 count=1326（1A=1326），node IDs 一致，RC3 最大误差 8.88e-16 dB。

## 累计 RC2（1200 s）

| 指标 | W1 | W1+W2 |
|---|---:|---:|
| count | 1326 | 495 |
| r 宽度 | 15.0 km | 15.0 km |
| v 宽度 | 2.00 | 2.00 |
| ψ 宽度 | 13.00° | 11.00° |

RC2 时间累积本身对 r 的收缩 = **0.0%**

## MAIN=235 Hz 双窗口 RC3，τ=0.5

config  z_true_m  tau_db  n_RC2_cum  n_RC2RC3_cum  count_contraction_vs_cumRC2  r_width_km  r_contraction_vs_cumRC2  v_width_mps  v_contraction_vs_cumRC2  psi_width_deg  psi_contraction_vs_cumRC2  true_retained  true_delta_J
  MAIN     180.0     0.5        495            41                     0.917172        15.0                      0.0          0.0                      1.0            3.0                   0.727273           True           0.0
  MAIN     200.0     0.5        495            43                     0.913131        15.0                      0.0          0.2                      0.9            3.0                   0.727273           True           0.0
  MAIN     220.0     0.5        495            41                     0.917172        15.0                      0.0          0.0                      1.0            3.0                   0.727273           True           0.0

## 距离别名持久性（MAIN, z_true=200）

 r_track_km  present_in_cumRC2  n_at_r     best_J12  delta_J12  best_v  best_psi_deg  best_z_star  survives_tau05  is_true_r
       45.0               True      27 4.527005e-02   0.045270     2.0           4.0        240.0            True      False
       50.0               True      30 3.395205e-13   0.000000     2.0           5.0        200.0            True       True
       56.0               True      32 3.339010e-02   0.033390     2.0           4.0        230.0            True      False
       58.0               True      34 2.450419e-02   0.024504     2.0           4.0        220.0            True      False
       60.0               True      37 4.323237e-02   0.043232     2.0           4.0        180.0            True      False

## SURVIVORS_BY_RANGE（τ=0.5）

config  z_true_m  tau_db  r0_km  n_cumRC2_at_r  n_survivors_at_r  any_survivor  best_J_at_r  delta_J_at_r
  MAIN     200.0     0.5   45.0             27                 3          True 4.527005e-02      0.045270
  MAIN     200.0     0.5   46.0             28                 3          True 2.701007e-02      0.027010
  MAIN     200.0     0.5   47.0             27                 2          True 4.133768e-02      0.041338
  MAIN     200.0     0.5   48.0             25                 2          True 4.254565e-02      0.042546
  MAIN     200.0     0.5   49.0             28                 2          True 1.568004e-02      0.015680
  MAIN     200.0     0.5   50.0             30                 2          True 3.395205e-13      0.000000
  MAIN     200.0     0.5   51.0             31                 2          True 7.575130e-02      0.075751
  MAIN     200.0     0.5   52.0             30                 4          True 7.794578e-02      0.077946
  MAIN     200.0     0.5   53.0             31                 2          True 4.296739e-02      0.042967
  MAIN     200.0     0.5   54.0             33                 3          True 2.620993e-02      0.026210
  MAIN     200.0     0.5   55.0             34                 3          True 3.434823e-02      0.034348
  MAIN     200.0     0.5   56.0             32                 3          True 3.339010e-02      0.033390
  MAIN     200.0     0.5   57.0             32                 3          True 5.733344e-02      0.057333
  MAIN     200.0     0.5   58.0             34                 3          True 2.450419e-02      0.024504
  MAIN     200.0     0.5   59.0             36                 3          True 1.782149e-02      0.017821
  MAIN     200.0     0.5   60.0             37                 3          True 4.323237e-02      0.043232

## 排序诊断（MAIN, z_true=200）

config  z_true_m  true_rank       true_J  true_z_star  best_false_J  best_false_r_km  best_false_v  best_false_psi_deg    J_p10  J_median    J_p90        J_min  n_acc
  MAIN     200.0          1 3.395205e-13        200.0       0.01568             49.0           2.0                 5.0 0.787819   2.42449 4.134565 3.395205e-13    495

## 完整性

见 `INTEGRITY_GATES.md`。

## 关键结构发现

**双窗时间推进未能打破距离脊线。**

1. **RC2 累计（1200 s bearing）本身不收缩距离**：r 宽度仍 15 km。方位观测的时间积累对绝对距离无贡献。
2. **RC3 双窗也不收缩距离**：τ=0.5 时 41–43 个幸存者仍横跨 45–60 km。
3. **每一个 r 网格点都有幸存者**（每点 2–4 个），45/56/58/60 km 的假极小 δJ₁₂ 均 < 0.05 dB。
4. **深度 profile 吸收了距离差**：各距离假解靠不同 z* 补偿（45 km→z=240, 56→230, 58→220, 60→180）。
5. **v/ψ 仍然强力收缩**：v 塌缩到单值 2.0 m/s，ψ 从 13° 缩到 3°。

### 科学结论

$$
\text{source-level-free relative TL shape} + \text{straight-line continuous bearing}
\not\Rightarrow \text{absolute range resolution}
$$

与 A05 静态距离盲区一致。仅靠当前 observable + 直航，绝对距离缺少锚点。

下一步应转向**有针对性的距离锚定因素**（RC1 会聚区绝对距离先验，或小幅平台转向），而非继续堆声学鲁棒性。

## 标签

- `SINGLE_WINDOW_DEMO_REALIZATION_EXTENDED`
- `PER_WINDOW_DEMEAN_SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE`
- `DIRECT_F_MATRIX_ON_TRAJECTORY`

## 未做

平台转角、absolute TL、已知源级、RC1 先验、新声学特征、S1、P5。
