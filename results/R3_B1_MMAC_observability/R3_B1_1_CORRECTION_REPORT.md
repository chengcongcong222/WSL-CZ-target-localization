# R3-B1.1 修正报告：MMAC 观测关口

UTC：2026-09-21T03:26:14.655763+00:00

## 0. 撤销与口径

- 旧判定 `B1_DELAY_BANDWIDTH_LIMITED` **撤销**，暂记 `B1_DECISION_PENDING_CORRECTION` 后由本修正给出终态。
- 旧文件保留；`R3_B1_DECISION_SUPERSEDED.json` 记录替代关系。
- 受控 direct+surface 模型保留为 **50 km 弱几何 sanity check**，**不是** E-STD 第一 CZ MMAC 性能。
- `HLA_ELEVATION_PROJECTION_WEAK=False` **不冻结**；修正后以实际分支 Δu + 导向相关为准。

## 1. 修正后时延模糊（核心）

| 源 | 代表 | 间距 Hz | T_amb s | 局部主瓣 s | 第一零点宽 s | 全局旁瓣 | CRB στ s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| S0_DENSE | continuous broadband | 1.000 | 1.0000 | 0.00520 | 0.00820 | 0.491 | 0.000244 |
| S0_23POINT | discrete comb control (NOT true broadband) | 10.227 | 0.0978 | 0.00500 | 0.00800 | 1.000 | 0.000235 |
| S1 | sparse lines | 36.000 | 0.0278 | 0.00640 | 0.09820 | 0.793 | 0.000296 |
| S2 | sparse lines | 35.000 | 0.0286 | 0.00560 | 0.01780 | 0.981 | 0.000263 |

**必须分开**：

- 局部时延分辨（S0_DENSE 主瓣 ≈ **0.00520 s**）与 CRB（~0.2–0.3 ms）一致量级；
- 全局无模糊：S0_23POINT 人为梳状 T_amb≈**0.0978 s**，旁瓣 **1.000** —— **不能**代表连续宽带 S0；
- S0_DENSE 全局旁瓣 **0.491**。

禁止再用“所有 A≥0.5 的最大 τ 跨度”当主瓣宽（旧 300–400 ms 为 **计算错误**）。

## 2. 精化本征射线（根求解）

- 准入：|z_hit−z_r|≤**0.1 m**（二分根求解 + 射线积分）
- branch_id：launch/φ/τ/topology **连续关联**，禁止 `round(launch/2°)`
- 分支 id 数 **12**；stable_for_fisher **4**；最大 residual **0.099 m**

| branch | z_s | n_ranges | r覆盖 km | launch drift° | φ drift° | τ drift s | max resid m | stable |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| E0 | 190.0 | 1 | 50.0-50.0 | 0.000 | 0.000 | 0.00000 | 0.006 | False |
| E1 | 190.0 | 1 | 50.0-50.0 | 0.000 | 0.000 | 0.00000 | 0.071 | False |
| E2 | 200.0 | 6 | 45.0-60.0 | 3.281 | 3.269 | 10.56102 | 0.081 | True |
| E3 | 200.0 | 3 | 45.0-49.5 | 2.360 | 2.394 | 2.87438 | 0.099 | True |
| E4 | 200.0 | 8 | 45.0-60.0 | 3.957 | 3.957 | 9.93371 | 0.059 | True |
| E5 | 200.0 | 4 | 49.0-60.0 | 4.204 | 4.188 | 7.14116 | 0.073 | True |
| E6 | 200.0 | 1 | 60.0-60.0 | 0.000 | 0.000 | 0.00000 | 0.074 | False |
| E7 | 200.0 | 1 | 60.0-60.0 | 0.000 | 0.000 | 0.00000 | 0.023 | False |
| E8 | 200.0 | 1 | 60.0-60.0 | 0.000 | 0.000 | 0.00000 | 0.073 | False |
| E9 | 210.0 | 1 | 50.0-50.0 | 0.000 | 0.000 | 0.00000 | 0.017 | False |
| E10 | 210.0 | 1 | 50.0-50.0 | 0.000 | 0.000 | 0.00000 | 0.017 | False |
| E11 | 210.0 | 1 | 50.0-50.0 | 0.000 | 0.000 | 0.00000 | 0.028 | False |

正确表述：**粗扫描显示第一 CZ 存在多径候选；稳定因果分支以本表 residual/coverage 为准。**

## 3. 实际分支上的 HLA 方向余弦

- corr<0.95 可分比例：14 m **0.000**，70 m **0.433**
- **HLA_ELEVATION_PROJECTION_WEAK = True**
- Corrected gate uses actual E-STD branch Delta-u + steering corr; prior False flag from composite CRB approximation is NOT frozen.

## 4. E-STD 实际分支 Fisher（非 controlled_paths）

| case | rank | CRB_r m | CRB_zs proxy m | cond F |
| --- | --- | --- | --- | --- |

## 5. ray–modal 时延交叉检查（sanity）

| r km | modal peaks | ray admitted taus s |
| --- | --- | --- |
| 49.0 | 0 | [32.66817, 33.16912, 33.16913, 33.23237] |
| 50.0 | 0 | [33.33054, 33.81346] |
| 51.0 | 0 | [33.99292, 34.45896, 34.45898] |

## 6. 修正后判定

### `B1_DELAY_AMBIGUITY_LIMITED`

Local delay mainlobe ms-scale OK (S0_DENSE 0.00520000000000001 s) but global ambiguity remains (S0_23 comb side=0.9999551387640602, S0_DENSE side=0.4907037600350636). CRB vs ambiguity must stay separated.

**下一步**：do not kill MMAC; resolve global delay ambiguity or restrict to local peak tracking

允许终态：`B1_MMAC_PHYSICS_CONFIRMED` / `B1_PHYSICS_ONLY_HLA_LIMITED` / `B1_DELAY_DOMINANT_POSSIBLE` / `B1_DELAY_AMBIGUITY_LIMITED` / `B1_NO_STABLE_MULTIPATH_IDENTITY` / `B1_NO_COMPLEMENTARY_INFORMATION`。

## 7. 停止

- 不进入 B2 / RC3-C / P5 / Monte Carlo RMSE
- 不增加新物理特征
- **B1.1 完成后停止；判定 `B1_DELAY_AMBIGUITY_LIMITED`**
