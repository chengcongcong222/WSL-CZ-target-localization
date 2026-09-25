# R3-RC3-REANCHOR-1 深度角色信息价值审计

UTC: 2026-09-25T04:19:11.212460+00:00

## 判定

### `DEPTH_ROLE_MIXED_OR_UNRESOLVED`

分布跨阈值；frac_strong=0.00 frac_non=0.25 med50=0.7088021089831597

## 测量

`SOURCE_LEVEL_FREE_RELATIVE_TL_SHAPE`（4 频各自去均值 TL）· E-STD · 30 对 P2/P3 困难样本 · z_true∈{180,200,220} · z∈150:5:250

核心：ρ_z = J_profile/J_fixed（错误轨迹可否靠改深度冒充真轨迹）。

## 机制汇总

mechanism  n_cases  n_informative  median_rho_z  frac_nonmaterial  frac_moderate  frac_strong  median_rho_50m  median_rho_25m  median_rho_10m
        A       30             27      0.750697          0.407407       0.592593          0.0        0.794492        0.849170        0.921654
        B       30             30      0.700610          0.200000       0.800000          0.0        0.700610        0.724240        0.752851
        C       30             30      0.634919          0.166667       0.833333          0.0        0.634919        0.671927        0.702224

## 边界

- z_alt* = `WRONG_TRAJECTORY_DEPTH_COMPENSATION_DIAGNOSTIC`（非深度 RMSE）
- 先验宽 = `ORACLE_DEPTH_PRIOR_WIDTH_CONTROL`
- 不做 Liang D(z) / 五维精度 / P5

## 含义（待 GPT 确认）

若 NUISANCE：RC3 目标改为「不知精深时用传播排除错误 r,v,ψ」，z 作隐变量。
