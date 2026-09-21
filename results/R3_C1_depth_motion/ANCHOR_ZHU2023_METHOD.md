# ANCHOR — Zhu 2023 深海 HLA 移动源深度估计

UTC：2026-09-21T07:52:08.343563+00:00

## 文献

Zhu F., Li F., Zhang Y., et al. *Moving source depth estimation in deep ocean direct arrival zone with a horizontal line array.* JASA Express Letters, 2023, 3(9):096003.

备用锚点：T.C. Yang 2015, synthetic aperture beamforming for moving source depth (PubMed 26428805)。

## 方法要点（本轮实现口径）

- **几何**：深海直达区；水平阵 endfire；源沿阵轴方向连续运动、深度恒定。
- **观测**：沿轨迹的声强/干涉结构 I(r) 或 I'(r)=I·r（距离补偿）。
- **机理**：直达 + 海面镜像路径干涉；相位差 Δφ=k(Ls−Ld) 随轨迹与深度变化 → 声强调制。
- **方法 A（优先）**：improved Fourier integral 型深度权 Q_F(z)：用候选深度预测的干涉核与 I' 的匹配。
- **方法 B（论文内对照）**：matched sound intensity structure：J(z)=corr(I'_obs, I'_pred(·|z))。
- **变量**：可用 sin(arrival angle) 或等价轨迹变量；本轮用 r 轨迹上的 I' 相关（Zhu-form，非 P4.5 CZ 轮廓）。
- **原文条件**：CW/可跟踪稳定频率；移动形成空间采样；深度搜索得到模糊峰。

## paper-condition 场景

- isovelocity c=1500 m/s；direct + surface image；endfire 移动；CW。
- 验证：true z 峰、track 变长峰变尖、表面/水下可分趋势。

## 观测条件分级

| 级 | 定义 | 本轮 |
| --- | --- | --- |
| C-S2 | 稳定线 201/235 Hz + 联合 | 主线 |
| C-S1 | 带内可跟踪机械线 | 对照 |
| C-S0 | 无线谱 | **METHOD_NOT_APPLICABLE_WITHOUT_TRACKABLE_LINE** |

禁止用 P4.5 简化 CZ 轮廓代替本方法观测形式。
