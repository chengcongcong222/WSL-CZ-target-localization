# P4-G4 修正回合报告

UTC：2026-09-20T04:43:28.628249+00:00

**范围**：仅修正 RC2→RC3 候选关系、SNR 统计尺度与展示口径；**不**重跑深度扫描、P2/P3/Bellhop；**不**启动 P5。

## 修正后 G4 状态：`G4_PARTIAL_BOUNDARY_VALIDATED`

Nested filtering holds for all revised cells (n_RC2+RC3<=n_RC2). Standardized residual J uses fixed threshold across SNR. Display is retain-first. Integrated RC2+RC3 statuses: OK r=1, INVALID=1, UNRESOLVED=1.

## 1. 结构修正：强制候选嵌套

$$\mathcal{C}_{RC2+RC3}=\mathcal{C}_{RC2}\cap\mathcal{C}_{RC3}$$

实现：先在标准化方位残差上得到 $\mathcal{C}_{RC2}$，再**仅在该集合内**用标准化声学残差继续筛选。

- source×SNR 全部格点 nested_ok：**True**
- C-L/M/U nested_ok：**True**

### 层间候选数（修正后）

| 场景 | RC1先验 | RC2 | RC2∩RC3 | 排除率 | nested |
| --- | --- | --- | --- | --- | --- |
| C-L | 48400 | 4180 | 819 | 0.804 | True |
| C-M | 48400 | 1075 | 1 | 0.999 | True |
| C-U | 48400 | 105 | 1 | 0.990 | True |

旧版 C-U 中 445→5122 的反向增加已消除：修正后 RC3 只能减少或保持候选数。

## 2. SNR 统一标准化残差

固定同一形式，**不**按 SNR 重调阈值：

$$J=\sum\frac{|y-\hat{y}|^2}{\sigma_n^2}$$

- 方位：$\sigma_n=\sigma_\theta$，阈值 $\Delta\chi^2\approx13.28$（固定）
- 声学特征：$\sigma_n^2 \propto 10^{-\mathrm{SNR}/10}/F$，阈值固定；只改变测量噪声尺度，不改接受规则

### source×SNR（RC2+RC3，retain 优先）

| source | SNR | n RC2→RC3 | r状态 | r宽度 | r收缩 | retain_r |
| --- | --- | --- | --- | --- | --- | --- |
| S0 | 20 | 975→1 | OK | 0.000 | 1.000 | 1.00 |
| S0 | 10 | 985→1 | INVALID_NOT_COVERED | n/a | n/a | 0.00 |
| S0 | 0 | 1040→2 | INVALID_NOT_COVERED | n/a | n/a | 0.00 |
| S0 | -10 | 1075→1070 | UNRESOLVED | 15.000 | 0.000 | 1.00 |
| S1 | 20 | 1070→1 | INVALID_NOT_COVERED | n/a | n/a | 0.00 |
| S1 | 10 | 945→2 | OK | 0.000 | 1.000 | 1.00 |
| S1 | 0 | 1070→6 | INVALID_NOT_COVERED | n/a | n/a | 0.00 |
| S1 | -10 | 1085→999 | UNRESOLVED | 15.000 | 0.000 | 1.00 |
| S2 | 20 | 1070→1 | INVALID_NOT_COVERED | n/a | n/a | 0.00 |
| S2 | 10 | 930→1 | INVALID_NOT_COVERED | n/a | n/a | 0.00 |
| S2 | 0 | 1070→53 | UNRESOLVED | 14.000 | 0.067 | 1.00 |
| S2 | -10 | 1075→1073 | UNRESOLVED | 15.000 | 0.000 | 1.00 |

在 OK 格点上，SNR 反常现象已通过统一 $J$ 尺度得到抑制或标注为 UNRESOLVED/INVALID。

## 3. retain 优先于 width

若 `truth retention=0`：

- 状态记 **`INVALID_NOT_COVERED`**
- **不报告**“宽 0”
- **不**把 contraction 画成 1
- 核心图 5 使用深红 INVALID 块，而不是高性能柱

### 三档综合（RC2+RC3，修正展示）

| 场景 | r | θ | z | v | ψ | n RC2→RC3 |
| --- | --- | --- | --- | --- | --- | --- |
| C-L | UNRESOLVED (w=15.000) | w=0.000, c=1.000 | UNRESOLVED (w=100.000) | UNRESOLVED (w=2.000) | UNRESOLVED (w=30.000) | 4180→819 |
| C-M | **INVALID_NOT_COVERED** | w=0.000, c=1.000 | UNRESOLVED (w=n/a) | **INVALID_NOT_COVERED** | **INVALID_NOT_COVERED** | 1075→1 |
| C-U | w=0.000, c=1.000 | w=0.000, c=1.000 | **INVALID_NOT_COVERED** | **INVALID_NOT_COVERED** | w=0.000, c=1.000 | 105→1 |

## 4. 仅改口径、不重跑的结论

### 小机动（冻结表述）

> 0–15° 小转向改变运动学候选结构，但在本研究的远距离、有限观测时间条件下，**未观察到稳定单调的距离或航向性能提升**；其作用更多体现为候选几何变化，而非直接高精度测距。

### 环境 E0/E1/E2（非核心性能）

> 环境扰动可显著改变候选结构；当前三点**不足以**形成单调鲁棒性边界。E1 出现的 r/z 零宽尖峰**不得**解释为“轻度失配反而提升性能”。

## 5. 不受本次修正影响的项目结论（冻结）

- RC2: long-range pure bearing cannot provide practical ranging; collinear r-v ridge is intrinsic.
- RC3: under frozen theory propagation + pre-locked hard pairs, CZ provides independent increment on range/range-trajectory candidates.
- Depth: z=200/220 m essentially unconstrained in current model; depth remains UNRESOLVED.
- Maneuver wording (no re-run): 0–15° small turns change kinematic candidate structure, but no stable monotonic range/course performance gain in this far-range limited-T study.
- Environment (no re-run): mismatch can change candidate structure; three points E0/E1/E2 insufficient for a monotone robustness boundary; E1 spikes are not evidence that mild mismatch improves performance.

## 6. 本回合输出

- `results/P4_G4_performance_boundary/P4_G4_REVISION.md`
- `results/P4_G4_performance_boundary/boundary_source_snr_revised.csv`
- `results/P4_G4_performance_boundary/integrated_CL_CM_CU_revised.csv`
- `results/P4_G4_performance_boundary/rc1_rc2_rc3_contraction_revised.csv`
- `results/P4_G4_performance_boundary/figures/core_fig2_source_snr.svg`
- `results/P4_G4_performance_boundary/figures/core_fig5_five_state.svg`
- `results/P4_G4_performance_boundary/figures/core_fig6_stage_contraction.svg`

## 7. 停止条件

- 未启动 P5
- 未重跑深度 / P2 / P3 / Bellhop / T×σθ / 机动扫描
- 允许终态：`G4_PARTIAL_BOUNDARY_VALIDATED` 或 `G4_BOUNDARY_LOGIC_NOT_RESOLVED` → 本轮为 **G4_PARTIAL_BOUNDARY_VALIDATED**
