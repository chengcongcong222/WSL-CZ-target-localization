# R3_FINAL_SCIENTIFIC_STORY

UTC: 2026-09-30

Canonical sync: 2026-10-02。基线 `d3394e3cc3c2f03aaa88e9ff24e6b3f4e8bf2e6f` 的参数候选集审计已获 GPT 独立审计通过，R3 正式冻结。本次只同步文档。

## 问题1：为什么单阵直航难以定距？

RC2 压缩方位并维持运动学候选；直航 1200 s 后 r-width 仍为 15 km。RC3 进一步压缩 v/ψ，但留下距离、profiled depth 与运动补偿脊线。以 z_true=200 m 为例，S4 的43个候选形成一个4D连通分量，16个距离格仍全部保留。

## 问题2：小转向解决了什么？

15° 转向形成约 300 m 运动虚拟基线，与 RC3 传播约束协同，在部分深度切断距离脊线。但单频仍有深度依赖（z=200 最难）。

## 问题3：为什么多频有效？

不同频率的错误 r–z 候选不一致。201/235/283 Hz 各自保留不同别名集（Jaccard 0.0–0.6），共享同一 profiled depth 后，6 个跨频不兼容点被联合排除。

## 问题4：什么条件下能出现距离锚？

已测试：50 km 真值、15° 转向、201+235+283 Hz（或 FOUR）、E0 匹配环境、1200 s 双窗口。三深度均只保留 50 km 距离格。

## 问题5：当前测试揭示哪些边界？

幅漂、tracked frequency drift 与 SSP 失配需要分别描述：

- 幅漂：truth rank 在 tested 1 dB 下保持；strict single-bin 从0.25 dB起不再统一保持。0.5 dB仍保留明显的候选集合收缩。
- tracked frequency drift：在 perfect tracking 条件下，TRIPLE 与 FOUR 到已测试的2%频漂均保持 strict range anchor。这不代表未知多普勒或频率跟踪误差达到2%时也鲁棒。
- SSP：E1 在两配置、三个深度中的真值总体 rank 为2–24，E2为6–29；最佳距离排序退化，距离候选空间重新打开。

采用 `SSP_MODEL_MATCH_IS_A_CRITICAL_DEPENDENCY_IN_TESTED_STRESS`。这些结果没有建立对实际海洋误差源的普遍重要性排序。

依据已有[幅漂最终指标](../R3_RC23_CLOSEDLOOP/R3_RC23_TURN_AMPLITUDE_CLOSEOUT_FIX/AMPLITUDE_CASE_METRICS_FINAL.csv)、[tracked频漂FIX2指标](../R3_RC23_CLOSEDLOOP/R3_RC23_TURN_TRACKED_FREQ_ROBUSTNESS_FIX2/TRACKED_FREQ_ANCHOR_CASES_FIXED.csv)与[SSP case表](../R3_RC23_CLOSEDLOOP/R3_RC23_TURN_SSP_MISMATCH/SSP_MISMATCH_ANCHOR_CASES.csv)同步上述口径，本轮未重算这些结果。

## 问题6：最终候选质量是什么？

在当前 matched synthetic control 中，TRIPLE 与 FOUR、三个真实深度均得到同一最终集合：

\[
\Omega_{\rm final}=\{(50\mathrm{km},0^\circ,2\mathrm{m/s},4^\circ),
(50\mathrm{km},0^\circ,2\mathrm{m/s},5^\circ)\}.
\]

| 参数 | 最终 survivor set | 当前结论 |
|---|---|---|
| r | 50 km | 单一离散网格值 |
| theta | 0° | 单一离散网格值 |
| v | 2 m/s | 单一离散网格值 |
| psi | 4°、5° | 剩余1°离散歧义 |
| z | profiled nuisance | 不作为测深结果 |

within frozen 114576-node RC1-conditioned discrete grid:

```text
final S7/S8 survivor count = 2
excluded fraction = 99.998254%
RC1_contraction_ratio = NOT_IDENTIFIABLE
```

这不是对任意先验空间的排除率。实际分支结构见[最终候选收缩 DAG](R3_FINAL_TECHNICAL_ARCHITECTURE.md)。

六个 matched case 中真值 `(50,0,2,5°)` 均为唯一全局最优，J_truth=0、n_strictly_better=0、n_tied=1。top-1 精确落在真值节点，与全 survivor-set 的保证分开报告。

逐节点与评分证据保留在[最终候选表](../R3_FINAL_PARAMETER_CONTRACTION/FINAL_SURVIVOR_STATE_TABLE_FIXED.csv)与[三深度审计表](../R3_FINAL_PARAMETER_CONTRACTION/PARAMETER_CONTRACTION_BY_DEPTH_FIXED.csv)。

集合仍保留4°航向节点，最坏航向相对诊断为20%；因此 `final_nominal_universal_lt10pct_set_bound=NOT_ESTABLISHED`。这不应写成当前 top-1 估计误差20%。甲方联合误差范数尚未定义，theta_truth=0°时通常的相对百分比误差未定义，z又是 nuisance parameter。本项继续标记 `DIAGNOSTIC_ONLY_NOT_REQUIREMENT_VERIFICATION`，不作为正式验收通过或不达标判定。

## 最终落脚

$$
\boxed{
\text{匹配控制场景中距离、方位、速度收缩至单格；剩余核心状态歧义为航向4°/5°}
}
$$

环境失配会重新打开距离候选空间。这构成当前716项目申请中可以使用的技术结论与边界。R3 在这些已测试条件下冻结，P5 未开启。
