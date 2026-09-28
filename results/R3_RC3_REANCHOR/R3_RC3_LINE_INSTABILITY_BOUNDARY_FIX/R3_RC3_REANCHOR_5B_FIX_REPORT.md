# R3-RC3-REANCHOR-5B-FIX

UTC: 2026-09-28T07:34:51.584019+00:00

## 判定

### `TRACKABLE_LINE_INSTABILITY_ROBUST_IN_TESTED_RANGE`

ZERO_STRESS=PASS; amp1=True amp05=True F_TR@2%=True F_NOM@1%=False aux=['TRACKABLE_LINE_REQUIRES_FREQUENCY_TRACKING']

## Gate 链

1. **ZERO_STRESS_5A_IDENTITY_GATE = PASS**
   - 逐行对齐 5A raw（450 行 × J_true*/J_alt*/ΔJ）
   - 全局最大误差 = 7.105e-15 dB（容差 1e-9）
2. **AMPLITUDE_BRANCH_INTEGRITY = PASS**
   - COMMON `[+1,+1,+1,+1]` vs LINE_SPECIFIC `[+1,-1,+1,-1]`
   - A=0 完全一致；A>0 注入数组不同
3. **RANGE_GRID_INTERPOLATION_AUDIT**
   - 状态：`RANGE_GRID_PATH_REJECTED_INSUFFICIENT_ACCURACY_USE_DIRECT`
   - 25/10/5 m 均未达 demeaned RMS≤0.05 dB 且 corr≥0.999
   - **主路径不使用 range-grid**，而是 F_m(r(t)) 直接 |p|；grid 加速正式废弃

## observable

$$\tilde L_f(t)=L_f(t)-\overline{L_f(t)}$$

- 幅漂：L_prop + a_f(t) 后再轨迹窗去均值
- 频漂：先拼完整 L_f(t)(t)，再整条线时间窗去均值
- F_TRACKED 模板用观测 f(t)；F_NOMINAL 模板恒 f_0

## 相对旧 5B

| 项 | 旧 5B | 5B-FIX |
|---|---|---|
| TL 去均值 | 全 44–61 km R_GRID | 轨迹时间窗 |
| 零扰动 vs 5A | 不一致 | 逐行一致 |
| LINE_SPECIFIC | 字符串未匹配，未执行 | multiplier 显式 [±1] |
| 判定 | REANCHOR5B_SUPERSEDED_PENDING_OBSERVABLE_FIX | `TRACKABLE_LINE_INSTABILITY_ROBUST_IN_TESTED_RANGE` |

## 辅助标签

- `TRACKABLE_LINE_REQUIRES_FREQUENCY_TRACKING`（若 F_NOMINAL@1% 未过）
- `FREQUENCY_DRIFT_FORWARD_MODEL_SPOTCHECK_PASSED`（自旧 5B 保留）
- `WEAKEST_IDEAL_SINGLE_LINE_235HZ`（单列）

## 未做

combined corner、旧 S1 频率迁移、S0、SSP、IID TL、intermittent、Liang、P5。
