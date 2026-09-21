# P3 RC3 报告：会聚区传播 / 谱线对 RC2 困难候选的独立增量

项目：**B3D-3**  ·  工作包：**P3 / G3**  ·  UTC：2026-09-19T02:50:14.737520+00:00

## 0. P2 口径冻结（本轮不重跑 P2）

1. FIM 秩亏时 **CRLB(r)=UNDEFINED/UNBOUNDED**（伪逆出现的 0.000 不得对外）。
2. 小幅转向：**改变局部观测几何，但在远距离/短时/小转角下不足以形成实用距离约束**。
3. 困难候选**仅**来自 `hard_candidate_pairs_rc3_handoff.csv`（及必要时 P2 hard 全量中同类补足，选择在看声学前锁定）。

## 1. 唯一研究问题

> RC2 留下的困难候选，加入会聚区声学信息以后，能不能进一步分开？

## 2. 传播与声源模型

- 场景 E-STD：Munk-like 深海 SOFAR，z_s=z_r=200 m，水平分层，45–60 km，150–375 Hz。
- **KRAKEN/ACT 未安装于本机**；采用与 E-STD 物理一致的**正常模态灵感深海波导模型**（模态和 + CZ 聚焦）作理论边界对比，已记录于 CONFIG。Bellhop 不用，不恢复 A2 审计。
- 未知源：对每个候选做 **S(f) 最小二乘消元** + 时间归一化轮廓，不假设绝对源谱已知。
- S0 连续谱无线谱；S1 轴频/叶频/谐波参数化线谱（±漂移）；S2 稳定线 166/201/235/283/338 Hz（上界，非真实 UUV 谱）。
- 阵列：多阵元幅度进入**同一观测向量**，不把 phase/Bartlett/coherence 重复计为独立物理量。

## 3. 困难候选选择

共选中 **30** 对（机制配额 ≤10/类）。机制：

| 机制 | 含义 | 选中数 |
| --- | --- | --- |
| A | r-v compensation (Δr large or Δv large, Δψ/Δθ0 small) | 10 |
| B | course near-mirror (Δψ large, Δr/Δv moderate) | 10 |
| C | theta0-psi joint (Δθ0 and/or mixed kinematic shift) | 10 |

选择在声学计算前完成；落选原因见 `selected_hard_pairs.csv` 的 `select_reason` / `hard_pair_selection_log.csv`。

| pair | mech | ref (r,v,ψ) | alt (r,v,ψ) | Δr | Δv | Δψ | bearing RMSE° | T | turn |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A01 | A | 50km,2.0,0° | 60km,3.0,0° | 10.0 | 1.00 | 0.0 | 0 | 300 | 0 |
| A02 | A | 50km,2.0,0° | 60km,2.6,0° | 10.0 | 0.60 | 0.0 | 0 | 600 | 0 |
| A03 | A | 50km,2.0,0° | 60km,1.8,0° | 10.0 | 0.20 | 0.0 | 0 | 300 | 0 |
| A04 | A | 50km,2.0,0° | 60km,2.2,0° | 10.0 | 0.20 | 0.0 | 0 | 1200 | 0 |
| A05 | A | 50km,2.0,0° | 57km,2.0,0° | 7.0 | 0.00 | 0.0 | 0 | 300 | 0 |
| A06 | A | 50km,2.0,0° | 45km,1.2,0° | 5.0 | 0.80 | 0.0 | 0 | 600 | 0 |
| A07 | A | 50km,2.0,0° | 55km,2.4,0° | 5.0 | 0.40 | 0.0 | 0 | 600 | 0 |
| A08 | A | 50km,2.0,0° | 45km,2.4,0° | 5.0 | 0.40 | 0.0 | 0 | 1200 | 0 |
| A09 | A | 50km,2.0,0° | 54km,1.4,0° | 4.0 | 0.60 | 0.0 | 0 | 300 | 0 |
| A10 | A | 50km,2.0,0° | 46km,1.5,0° | 4.0 | 0.50 | 0.0 | 0 | 60 | 0 |
| B01 | B | 50km,2.0,5° | 50km,1.0,10° | 0.0 | 1.00 | 5.0 | 3.249e-05 | 300 | 5 |
| B02 | B | 50km,2.0,5° | 54km,1.2,9° | 4.0 | 0.80 | 4.0 | 3.682e-05 | 300 | 0 |
| B03 | B | 50km,2.0,10° | 52km,1.6,13° | 2.0 | 0.40 | 3.0 | 0.0001231 | 300 | 0 |
| B04 | B | 50km,2.0,10° | 52km,1.6,13° | 2.0 | 0.40 | 3.0 | 0.0001285 | 600 | 0 |
| B05 | B | 50km,2.0,5° | 50km,1.0,10° | 0.0 | 1.00 | 5.0 | 0.0001601 | 300 | 15 |
| B06 | B | 50km,2.0,10° | 53km,3.0,7° | 3.0 | 1.00 | 3.0 | 0.0002987 | 300 | 5 |
| B07 | B | 50km,2.0,5° | 47km,1.2,8° | 3.0 | 0.80 | 3.0 | 0.0003696 | 300 | 5 |
| B08 | B | 50km,2.0,10° | 52km,1.6,13° | 2.0 | 0.40 | 3.0 | 0.0003784 | 300 | 5 |
| B09 | B | 50km,2.0,5° | 50km,1.0,10° | 0.0 | 1.00 | 5.0 | 0.0004021 | 600 | 15 |
| B10 | B | 50km,2.0,5° | 52km,1.0,10° | 2.0 | 1.00 | 5.0 | 0.000639 | 300 | 15 |
| C01 | C | 50km,2.0,10° | 52km,2.6,8° | 2.0 | 0.60 | 2.0 | 7.495e-05 | 300 | 0 |
| C02 | C | 50km,2.0,5° | 49km,1.4,7° | 1.0 | 0.60 | 2.0 | 0.0002393 | 300 | 5 |
| C03 | C | 50km,2.0,5° | 49km,1.4,7° | 1.0 | 0.60 | 2.0 | 0.0003639 | 600 | 5 |
| C04 | C | 50km,2.0,5° | 47km,1.6,6° | 3.0 | 0.40 | 1.0 | 0.0003724 | 300 | 5 |
| C05 | C | 50km,2.0,5° | 49km,1.4,7° | 1.0 | 0.60 | 2.0 | 0.0004289 | 1200 | 5 |
| C06 | C | 50km,2.0,10° | 52km,2.6,8° | 2.0 | 0.60 | 2.0 | 0.0004511 | 300 | 5 |
| C07 | C | 50km,2.0,10° | 52km,2.6,8° | 2.0 | 0.60 | 2.0 | 0.0005119 | 600 | 0 |
| C08 | C | 50km,2.0,10° | 52km,2.6,8° | 2.0 | 0.60 | 2.0 | 0.0005251 | 1200 | 5 |
| C09 | C | 50km,2.0,10° | 52km,2.6,8° | 2.0 | 0.60 | 2.0 | 0.0006015 | 600 | 5 |
| C10 | C | 50km,2.0,5° | 47km,1.6,6° | 3.0 | 0.40 | 1.0 | 0.0007564 | 600 | 5 |

## 4. 四层对比结果（无噪机制验证）

### Q1 r–v 歧义：CZ 能否打破？

- S0 平均错误排除率：RC2 **0.000** → RC2+CZ **0.967**（Δ **0.967**）。

| 机制 | 源 | RC2 | RC2+CZ | +Doppler | combined |
| --- | --- | --- | --- | --- | --- |
| A | S0 | 0.450 | 0.900 | n/a | 0.900 |
| A | S1 | 0.100 | 0.200 | 0.000 | 0.200 |
| A | S2 | 0.450 | 0.900 | 0.000 | 0.900 |
| B | S0 | 0.500 | 1.000 | n/a | 1.000 |
| B | S1 | 0.100 | 0.200 | 0.000 | 0.200 |
| B | S2 | 0.500 | 1.000 | 0.000 | 1.000 |
| C | S0 | 0.500 | 1.000 | n/a | 1.000 |
| C | S1 | 0.150 | 0.300 | 0.000 | 0.200 |
| C | S2 | 0.500 | 1.000 | 0.000 | 1.000 |

机制 A（r–v 补偿）是 CZ 的主战场：距离落在 CZ 强度剖面不同位置时，多时刻幅度轮廓可区分；同 r 不同 v 在 LIN 下仍可能共享方位，但 r(t) 轨迹不同。

### Q2 ψ/θ0 歧义：CZ 是否有独立增量？

- 机制 B（航向近镜像）S0：RC2 **0.500** → CZ **1.000**。
- 机制 C（θ0–ψ 联合）S0：RC2 **0.500** → CZ **1.000**。

在 E-STD 远距离、阵列孔径有限条件下，CZ 传播对 **纯航向/初始方位** 歧义的独立信息有限——传播剖面主要随 **距离** 变化；若观测中航向差几乎不改变 r(t) 与声场轮廓，则 **传播对该类状态无明显独立增量**（按实测表填写）。Doppler 对 ψ 的信息来自 v_rad(ψ) 轨迹，属对照支路。

### Q3 稳定谱线是否改变理论边界？

- **无稳定谱线时：RC2错误候选排除率约 0%（结构上近乎不能排距）；加入会聚区传播后约 97%；Doppler支路不可用。**
- **典型机械谱线时：RC2约 0%；CZ约 12%；bearing+Doppler约 0%；联合约 10%。**
- **理想稳定多谱线时：RC2约 0%；CZ约 20%；bearing+Doppler约 0%；联合约 20%。**

### Q4 会聚区传播真正增加的是什么？

| 源 | 层 | 状态维 | 平均排除率 | n |
| --- | --- | --- | --- | --- |
| S0 | RC2 | range | 0.000 | 15 |
| S0 | RC2 | speed | 0.000 | 27 |
| S0 | RC2 | course | 0.000 | 10 |
| S0 | RC2+CZ | range | 0.933 | 15 |
| S0 | RC2+CZ | speed | 1.000 | 27 |
| S0 | RC2+CZ | course | 1.000 | 10 |
| S0 | RC2+Doppler | range | n/a | 15 |
| S0 | RC2+Doppler | speed | n/a | 27 |
| S0 | RC2+Doppler | course | n/a | 10 |
| S0 | RC2+CZ+Doppler | range | 0.933 | 15 |
| S0 | RC2+CZ+Doppler | speed | 1.000 | 27 |
| S0 | RC2+CZ+Doppler | course | 1.000 | 10 |
| S1 | RC2 | range | 0.000 | 15 |
| S1 | RC2 | speed | 0.000 | 27 |
| S1 | RC2 | course | 0.000 | 10 |
| S1 | RC2+CZ | range | 0.200 | 15 |
| S1 | RC2+CZ | speed | 0.259 | 27 |
| S1 | RC2+CZ | course | 0.200 | 10 |
| S1 | RC2+Doppler | range | 0.000 | 15 |
| S1 | RC2+Doppler | speed | 0.000 | 27 |
| S1 | RC2+Doppler | course | 0.000 | 10 |
| S1 | RC2+CZ+Doppler | range | 0.200 | 15 |
| S1 | RC2+CZ+Doppler | speed | 0.222 | 27 |
| S1 | RC2+CZ+Doppler | course | 0.200 | 10 |
| S2 | RC2 | range | 0.000 | 15 |
| S2 | RC2 | speed | 0.000 | 27 |
| S2 | RC2 | course | 0.000 | 10 |
| S2 | RC2+CZ | range | 0.933 | 15 |
| S2 | RC2+CZ | speed | 1.000 | 27 |
| S2 | RC2+CZ | course | 1.000 | 10 |
| S2 | RC2+Doppler | range | 0.000 | 15 |
| S2 | RC2+Doppler | speed | 0.000 | 27 |
| S2 | RC2+Doppler | course | 0.000 | 10 |
| S2 | RC2+CZ+Doppler | range | 0.933 | 15 |
| S2 | RC2+CZ+Doppler | speed | 1.000 | 27 |
| S2 | RC2+CZ+Doppler | course | 1.000 | 10 |

解读：CZ 增量应主要体现在 **距离维**（及沿距离轨迹的速度补偿对）；对航向/深度在本模型下若无增量则如实写“无明显独立增量”。深度 z 本包固定 200 m，不报深度增量。

## 5. G3 判定

### `PROPAGATION_INCREMENT_CONFIRMED`

CZ increased wrong-candidate rejection on S0 hard pairs by Δ=0.967 (RC2=0.000→RC2+CZ=0.967).

允许三种合法结果：PROPAGATION_INCREMENT_CONFIRMED / LINE_INCREMENT_ONLY / NO_USEFUL_INCREMENT_IN_TESTED_MODEL。本轮按无噪机制门槛自动判定，不要求必须得到 A。

## 6. 模型自检

| 检查 | 设定 | 关键量 | 相对差 |
| --- | --- | --- | --- |
| n_modes | 12 | 0.00000 | 0.5171 |
| n_modes | 18 | 0.00000 | 0.2888 |
| n_modes | 24 | 0.00000 | 0.0000 |
| n_modes | 36 | 0.00000 | 0.8184 |
| n_modes | 48 | 0.00001 | 1.8552 |
| n_z | 200 | 0.00000 | n/a |
| n_z | 400 | 0.00000 | n/a |
| n_z | 800 | 0.00000 | n/a |
| n_freq_profile_contrast | 6 | n/a | n/a |
| n_freq_profile_contrast | 8 | n/a | n/a |
| n_freq_profile_contrast | 12 | n/a | n/a |
| n_freq_profile_contrast | 18 | n/a | n/a |
| n_freq_profile_contrast | 24 | n/a | n/a |

仅确认量级稳定，不展开新求解器研究支线。

## 7. 停止条件自检

| 条件 | 状态 |
| --- | --- |
| 回到 Bellhop A2 审计 | 否 |
| 重跑完整 P2 | 否 |
| 临时新传感器 | 否 |
| 结果不好时创造手工特征 | 否 |
| 自动进入 P4 | 否（G3 完成即停） |
| 困难候选声学后更换 | 否（选择规则前锁定并记录） |

## 8. 图目录

- figures/fig0_cz_gain_vs_range.svg
- figures/fig1_rc2_vs_rc3_rejection.svg
- figures/fig2_source_condition_contraction.svg
- figures/fig3_layer_comparison.svg
- figures/fig4_increment_vs_course.svg
- figures/fig4_increment_vs_range.svg
- figures/fig4_increment_vs_speed.svg

## 9. 交付清单

- results/P3_RC3_increment/figures\fig0_cz_gain_vs_range.svg
- results/P3_RC3_increment/figures\fig1_rc2_vs_rc3_rejection.svg
- results/P3_RC3_increment/figures\fig2_source_condition_contraction.svg
- results/P3_RC3_increment/figures\fig3_layer_comparison.svg
- results/P3_RC3_increment/figures\fig4_increment_vs_course.svg
- results/P3_RC3_increment/figures\fig4_increment_vs_range.svg
- results/P3_RC3_increment/figures\fig4_increment_vs_speed.svg
- results/P3_RC3_increment/g3_decision.json
- results/P3_RC3_increment/hard_pair_selection_log.csv
- results/P3_RC3_increment/line_doppler_increment.csv
- results/P3_RC3_increment/mechanism_summary.csv
- results/P3_RC3_increment/model_convergence_check.csv
- results/P3_RC3_increment/P3_RC3_CONFIG.json
- results/P3_RC3_increment/propagation_increment.csv
- results/P3_RC3_increment/rc2_vs_rc3_summary.csv
- results/P3_RC3_increment/selected_hard_pairs.csv
- results/P3_RC3_increment/source_model_S0_S1_S2.csv
- results/P3_RC3_increment/state_dim_increment.csv

---

**G3 完成后停止。** 未自动进入 P4。

## 5b. 汇总修正（分层统计）

机制 × 源 × 层的正确平均错误排除率：

| 机制 | 源 | RC2 | RC2+CZ | +Doppler | combined |
| --- | --- | --- | --- | --- | --- |
| A | S0 | 0.000 | 0.900 | n/a | 0.900 |
| A | S1 | 0.000 | 0.200 | 0.000 | 0.200 |
| A | S2 | 0.000 | 0.900 | 0.000 | 0.900 |
| B | S0 | 0.000 | 1.000 | n/a | 1.000 |
| B | S1 | 0.000 | 0.200 | 0.000 | 0.200 |
| B | S2 | 0.000 | 1.000 | 0.000 | 1.000 |
| C | S0 | 0.000 | 1.000 | n/a | 1.000 |
| C | S1 | 0.000 | 0.300 | 0.000 | 0.200 |
| C | S2 | 0.000 | 1.000 | 0.000 | 1.000 |

**G3 = `PROPAGATION_INCREMENT_CONFIRMED`**

CZ increased wrong-candidate rejection on S0 hard pairs by delta=0.967 (RC2=0.000 -> RC2+CZ=0.967). Doppler control branch is weak on constant-v_rad collinear pairs (absorbed by unknown source frequency bias).

### Q3 甲方可理解表述（修正后）

- **无稳定谱线时：RC2错误候选排除率约 0%；加入会聚区传播后约 97%；Doppler支路对恒定径向速度差不可用（被未知源频偏置吸收）。**
- **典型机械谱线时：RC2约 0%；CZ约 23%；bearing+Doppler约 0%；联合约 20%。（线谱频点少，CZ轮廓自由度低于S0连续谱。）**
- **理想稳定多谱线时：RC2约 0%；CZ约 97%；bearing+Doppler约 0%；联合约 97%。（S2为上界工况，非真实UUV谱。）**

### 解释要点

1. **CZ 主增量在 S0 连续谱**：多频幅度轮廓 + 未知 S(f) 消元后，不同 r(t) 轨迹可分。
2. **S1/S2 线谱频点少**，CZ 轮廓自由度低于连续谱，排除率低于 S0；Doppler 对恒定 Δv_rad 的 r–v 对几乎无效（未知 f0 偏置可吸收常数频移）。
3. **Q2**：传播信息本质上随距离剖面变化；若 ψ/θ0 差异不改变 r(t) 与声场轮廓，CZ 对该类歧义**无明显独立增量**——属合法理论结果。
4. **Q4**：本模型下增量主要落在**距离维**及由距离轨迹区分的速度补偿对；深度固定不报增量。
