# R3_FINAL_WORKING_ENVELOPE

UTC: 2026-09-30

Canonical sync: 2026-10-02。R3 已完成 GPT 独立参数审计并冻结；下列口径限定于已有测试证据。

## 搜索空间 vs 已验证真值

**搜索空间** = 候选网格；**已验证真值** = 测试用例。两者不同。

## Geometry

```text
first CZ (deep ocean)
receiver/tow depth ≈ 200 m
tested source truth depths = 180, 200, 220 m
candidate depth profiling = 150:5:250 m (21 nodes)
candidate range grid = 45:1:60 km (16 nodes)
validated anchor truth range = 50 km
```

**禁止**把 45–60 km candidate grid 写成"已验证所有真实距离 45–60 km"。

## Motion

```text
platform speed ≈ 2 m/s
target truth speed = 2 m/s
truth heading = 5°
turn = 15° at t=600 s
total observation = 1200 s (two windows)
bearing noise sigma = 0.1°
one frozen seed = 20260912
```

不能承诺：任意目标运动、所有转向角、所有观测时长。

## Spectral Condition

成功最小子集：**201+235+283 Hz**

项目选定 ideal trackable tones，不是已证明的真实 UUV 通用谱线。

## Source-level Requirement

每频、每窗口去均值：

```text
no absolute source level required
no cross-frequency absolute amplitude calibration required
```

但仍依赖可跟踪的线谱轨迹。

## 环境

```text
E-STD Munk-like SSP
H = 5000 m
rigid bottom
E0 = nominal (matched environment)
```

## 已验证工作包络汇总

| 轴 | 测试范围 | 结果 |
|---|---|---|
| 频率子集 | 15 subsets | 201+235+283 最小成功 |
| 幅漂 | A=0–2 dB 原阶段测试；candidate-level质量审计覆盖0–1 dB | truth rank在tested 1 dB保持；strict single-bin从0.25 dB起不再统一保持，0.5 dB仍明显收缩 |
| 频漂 | D=0–2% tracked；perfect tracking control | TRIPLE与FOUR到2%均保持strict range anchor |
| SSP | E1/E2 truth vs E0 template | E1总体rank 2–24，E2总体rank 6–29；距离排序降级 |
| 转向 | δ=15° | 已验证 |

`SSP_MODEL_MATCH_IS_A_CRITICAL_DEPENDENCY_IN_TESTED_STRESS` 描述当前测试依赖，不构成实际海洋误差源的普遍排序。

## 最终 nominal 四维候选质量

TRIPLE、FOUR与三个真实深度的六个matched case均保留 `(50 km,0°,2 m/s,4°)`、`(50 km,0°,2 m/s,5°)`，真值为唯一全局最优。距离、方位与速度均收缩至单一离散网格值，psi仍有1°离散歧义。z作为profile nuisance，不作为测深输出。

within frozen 114576-node RC1-conditioned discrete grid: final S7/S8 survivor count = 2; excluded fraction = 99.998254%。RC1_contraction_ratio = NOT_IDENTIFIABLE。

全survivor-set最坏航向相对诊断为20%；`final_nominal_universal_lt10pct_set_bound=NOT_ESTABLISHED`。top-1精确落在真值与集合保证分开报告。甲方联合误差范数未定义，theta真值为0°、z为nuisance；`DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION`，不构成正式验收判定。
