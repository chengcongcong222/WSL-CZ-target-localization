# R3_FINAL_WORKING_ENVELOPE

UTC: 2026-09-30

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
| 幅漂 | A=0–2 dB | 0.25–0.5 dB 集合有效 |
| 频漂 | D=0–2% tracked | 两配置到 2% 全保持 |
| SSP | E1/E2 | **降级排序** |
| 转向 | δ=15° | 已验证 |
