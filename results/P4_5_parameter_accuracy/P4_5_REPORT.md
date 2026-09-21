# P4.5 报告：五维参数统一估计精度基准

UTC：2026-09-20T06:16:37.298488+00:00  ·  目录：`results/P4_5_parameter_accuracy/`

本阶段**不**开发新方法；把成果从“候选排除率”转换为 **先验 → RC2 → RC2+RC3** 的参数 RMSE/MAE/覆盖率。

## 0. 实验设计

- 每场景 **N_truth=30** × **N_noise=20** = 600 次估计
- 真值连续 **off-grid** 采样；M1/M2 共用真值、方位噪声、候选空间与优化框架
- M0=先验中心参考；M1=仅方位；M2=方位+CZ（时间轮廓特征，未知 S(f) 消去）
- 点估计：粗网格 argmin → 局部细化；CI：J≤J_min+Δχ²
- 覆盖率不足 → INVALID/OVERCONFIDENT，不把小 RMSE 称高精度

## 1. 主表 `parameter_accuracy_baseline_vs_method.csv`

| 场景 | 参数 | 先验std | RC2 RMSE | RC2覆盖 | RC2+RC3 RMSE | RC2+RC3覆盖 | 改善 | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B-L | r (km) | 4.330 | 5.107 | 1.00 | 5.982 | 0.40 | -0.171 | INVALID_COVERAGE |
| B-L | theta (deg) | 2.887 | 0.110 | 0.05 | 0.129 | 0.04 | -0.175 | INVALID_COVERAGE |
| B-L | z (m) | 28.868 | n/a | 1.00 | 21.730 | 0.70 | n/a | M2=OK |
| B-L | v (m/s) | 0.577 | 0.779 | 0.96 | 0.611 | 0.33 | 0.216 | INVALID_COVERAGE |
| B-L | psi (deg) | 8.660 | 6.516 | 0.15 | 8.591 | 0.13 | -0.318 | INVALID_COVERAGE |
| B-M | r (km) | 4.330 | 5.934 | 1.00 | 5.849 | 0.40 | 0.014 | INVALID_COVERAGE |
| B-M | theta (deg) | 2.887 | 0.079 | 0.01 | 0.103 | 0.03 | -0.315 | INVALID_COVERAGE |
| B-M | z (m) | 28.868 | n/a | 1.00 | 24.422 | 0.66 | n/a | M2=UNRESOLVED_OR_WEAK |
| B-M | v (m/s) | 0.577 | 0.790 | 0.89 | 0.598 | 0.36 | 0.243 | INVALID_COVERAGE |
| B-M | psi (deg) | 8.660 | 4.938 | 0.46 | 7.392 | 0.27 | -0.497 | INVALID_COVERAGE |
| B-U | r (km) | 4.330 | 7.541 | 0.14 | 4.467 | 0.04 | 0.408 | INVALID_COVERAGE |
| B-U | theta (deg) | 2.887 | 0.063 | 0.00 | 0.076 | 0.00 | -0.203 | INVALID_COVERAGE |
| B-U | z (m) | 28.868 | n/a | 1.00 | 23.726 | 0.07 | n/a | INVALID_COVERAGE |
| B-U | v (m/s) | 0.577 | 0.555 | 0.69 | 0.466 | 0.02 | 0.161 | INVALID_COVERAGE |
| B-U | psi (deg) | 8.660 | 3.237 | 0.56 | 2.926 | 0.01 | 0.096 | INVALID_COVERAGE |

## 2. 论文式表述（B-M 主结果）

- **r**：RMSE=5.849 但覆盖=0.40 → **INVALID_COVERAGE**，不称高精度。
- **v**：RMSE=0.598 但覆盖=0.36 → **INVALID_COVERAGE**，不称高精度。
- **theta**：RMSE=0.103 但覆盖=0.03 → **INVALID_COVERAGE**，不称高精度。
- **psi**：RMSE=7.392 但覆盖=0.27 → **INVALID_COVERAGE**，不称高精度。
- **z**：RC2 结构性弱/不可辨识；CZ 后 RMSE=24.422 m，95%宽=26.400，覆盖=0.66。

## 3. CZ 直接敏感度（B-M 切片）

| 参数 | Δ@+9.21 | 解读 |
| --- | --- | --- |
| r | 11.0000 | 距离敏感度有限 |
| v | 0.7000 | 单独弱，多来自 r–v 耦合 |
| z | n/a | 深度有一定敏感度 |
| psi | n/a | 航向有 CZ 信息 |
| theta | n/a | θ 有 CZ 信息 |

**结论口径**：CZ 直接观测贡献主要来自 **距离 r**；速度改善更多来自距离轨迹耦合；深度/航向改善有限。**不得**用“候选排除率”反推参数精度。

## 4. 图与文件

- figures/figA_rmse_baseline_vs_method.svg
- figures/figB_rmse_improvement.svg
- figures/figC_rc3_parameter_sensitivity.svg
- results/P4_5_parameter_accuracy/P4_5_CONFIG.json
- results/P4_5_parameter_accuracy/parameter_accuracy_baseline_vs_method.csv
- results/P4_5_parameter_accuracy/parameter_accuracy_trials.csv
- results/P4_5_parameter_accuracy/parameter_coverage.csv
- results/P4_5_parameter_accuracy/rc3_1d_sensitivity.csv
- results/P4_5_parameter_accuracy/rc3_rv_sensitivity.csv
- results/P4_5_parameter_accuracy/rc3_rz_sensitivity.csv

**P4.5 完成后停止；未进入 P5。**

---

## 6. 结果解释（必须与 P3 排除率区分）

P4.5 用**同一真值 / 同一噪声 / 同一候选空间**做了 3×600 次估计。主结论如下。

### 6.1 RC2 基线（冻结结论再次量化）

| 场景 | r RMSE (RC2) | r 95%宽 (RC2) | v RMSE (RC2) | 解读 |
| --- | --- | --- | --- | --- |
| B-L 近共线 | ≈5.1 km | ≈15 km（≈先验全宽） | ≈0.78 m/s | 结构性弱 / UNRESOLVED |
| B-M 开放几何 | ≈5.9 km | ≈14.7 km | ≈0.79 m/s | 距离–速度仍宽 |
| B-U 有利条件 | ≈7.5 km | ≈2.5 km 但覆盖仅≈0.14 | ≈0.56 m/s | 尖峰集过窄 → INVALID_COVERAGE |

**可写**：RC2 下距离 RMSE 约 5–8 km 量级，95% 候选宽常与先验同量级；**不能**写“RC2 已能测距”。

### 6.2 加入 CZ 后（M2）

- B-M：r RMSE 5.93→5.85 km（改善仅 ~1%），覆盖 0.40 → **INVALID_COVERAGE**
- B-M：v RMSE 0.79→0.60 m/s（~24%），覆盖 0.36 → 误差有降但**覆盖不足**，不称高精度
- B-U：r RMSE 7.54→4.47 km（~41%），覆盖 0.035 → 同样 **INVALID_COVERAGE**
- θ/ψ：M2 未稳定优于 M2 之外的 RC2；ψ 在 B-M 反而变差
- z：RMSE ~20–24 m，多为 UNRESOLVED / 弱约束

### 6.3 CZ 单参数敏感度（回答“RC3 看到了谁”）

B-M 切片，Δ@J_min+9.21 半高宽：

| 参数 | 半高宽 | 含义 |
| --- | --- | --- |
| r | **~11 km** | J_CZ(r) 峰**不尖锐** |
| v | ~0.7 m/s | 仍宽 |
| z / ψ / θ | 无稳定 +9.21 截断（曲线平） | 单独几乎不敏感 |

**严谨结论**：

> 1. P3 的“困难候选排除 ~97%”是**多维联合假设判别**，**不能**换算成“距离精度提高 97%”。  
> 2. P4.5 显示：在统一 RMSE/覆盖基准下，**当前理论传播模型 + 时间轮廓特征，尚未把距离压到高覆盖的小误差**；覆盖不足时标 INVALID，不得包装成高精度。  
> 3. CZ 的信息方向更接近“约束与距离轨迹相关的联合几何”，而不是给出 sharp 的单参数测距。  
> 4. 与冻结结论一致：RC2 不实用测距；深度仍弱；θ 主要来自方位。

### 6.4 甲方/论文允许的表述

- “RC2 下距离 RMSE≈X km，95% 候选宽≈Y km；加入 CZ 后 RMSE≈Z km，但 95% 覆盖仅≈C，故仅能说明候选结构变化，不能宣称达到高精度定位。”
- “RC3 在预锁定困难候选上可显著排除错误联合假设（P3），但全网格参数 RMSE 改善有限且覆盖率不足（P4.5）——两者层次不同。”
- 禁止：“排除 97% ⇒ 测距提高 97%”。

**P4.5 完成后停止。下一步是否进入 P5，由任务书决定。**
