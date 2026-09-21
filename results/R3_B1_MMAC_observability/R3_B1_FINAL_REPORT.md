# R3-B1-FINAL-CHECK 报告

UTC：2026-09-21T06:06:52.832098+00:00

**B1 自本报告起永久关闭**，不再派生 B1.x；不进入 B2 / RC3-C / P5。

## 0. 撤销

- 撤销 R3-B1.2 的 `B1_MMAC_PHYSICS_CONFIRMED`（delay-only 观测维数错误 + 疑似分支切换）
- 保留正结果：S0 连续宽带局部时延 ms 级；S1/S2 全局模糊；HLA 俯仰投影弱

## 1. 可观测分支（continuation + 模态能量）

- 五点 stencil 无 shoot-fill；连续性：topology + |Δφ|<3° + |Δτ|<50ms + |Δlaunch|<5°
- 模态能量容差：**15 ms**（与宽带主瓣/第一零点同量级）
- **M（OBSERVABLE_BRANCH）= 1**：['C0']

| id | topo | φ° | τ s | cont5 | energy | OBS |
| --- | --- | --- | --- | --- | --- | --- |
| C0 | refracted | -0.30 | 33.3305 | True | True | OBSERVABLE_BRANCH |
| C1 | surface_bounce | 13.15 | 33.8135 | False | False | NOT_OBSERVABLE |

诊断摘录见 `observable_branch_table.csv` 的 notes（含 13°→3.6° 类跳变）。

## 2. 时延观测构造（硬断言）

- delay-only 独立维数 = M−1 = **0**
- delay-only 估计 (r,z) **要求 M≥3** → 本轮 **允许=False**
- 搜索窗覆盖实际 Δτ：**±0.600 s**（max|Δτ|=0.500 s）
- S0 使用 0.25 Hz 密频积分，**不用** 1 Hz 采样制造 1 s 周期

| 源 | τ0 | 局部主瓣 ms | PSL | status |
| --- | --- | --- | --- | --- |
| S0_DENSE | 0.000 | 5.30 | 0.217 | **GLOBAL_DELAY_USABLE_IN_TEST_WINDOW** |
| S0_23POINT | 0.000 | 5.10 | 1.000 | **GLOBAL_AMBIGUOUS_COMB_CONTROL** |
| S1 | 0.000 | 6.30 | 0.910 | **GLOBAL_AMBIGUOUS** |
| S2 | 0.000 | 5.70 | 0.982 | **GLOBAL_AMBIGUOUS** |

## 3. Fisher + 不变量

- nesting / assert 总通过：**False**（nesting_ok=True）

| case | σθ | θ° | M | n_obs | rank | forbidden | CRB_r m | CRB_z m |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| IDEAL_ELEVATION+delay | 0.1 | 0.0 | 1 | 1 | 1 | True | n/a | n/a |
| HLA_u+delay | 0.1 | 0.0 | 1 | 1 | 1 | True | n/a | n/a |
| delay-only | 0.1 | 0.0 | 1 | 0 | 0 | True | n/a | n/a |
| IDEAL_ELEVATION+delay | 0.1 | 30.0 | 1 | 1 | 1 | True | n/a | n/a |
| HLA_u+delay | 0.1 | 30.0 | 1 | 1 | 1 | True | n/a | n/a |
| delay-only | 0.1 | 30.0 | 1 | 0 | 0 | True | n/a | n/a |
| IDEAL_ELEVATION+delay | 0.1 | 60.0 | 1 | 1 | 1 | True | n/a | n/a |
| HLA_u+delay | 0.1 | 60.0 | 1 | 1 | 1 | True | n/a | n/a |
| delay-only | 0.1 | 60.0 | 1 | 0 | 0 | True | n/a | n/a |

断言明细：`fisher_invariant_checks.csv`（含 n_delay=M−1、rank 限制、Fisher 半正定嵌套、CRB 单调）。

## 4. B1 最终判定（永久）

### `B1_NO_STABLE_MULTIPATH_IDENTITY`

Only 1 continuation-continuous AND modal-energy-supported branches at center.

**下一步**：B1 permanently closed; RC3-B MMAC not established

## 5. 关闭声明

- B1 状态：**PERMANENTLY_CLOSED**
- 不进入 B2 / RC3-C / P5；不派生 B1.x
- 不增加新物理特征
