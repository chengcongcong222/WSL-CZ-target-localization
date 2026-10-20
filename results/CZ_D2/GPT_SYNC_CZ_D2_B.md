# CZ-D2 Commit B：前置完整性门失败，未进入正式实验

设计SHA：`dc91982e67921eeb7a6aac4e12c482a980d5ad89`；分支 `research/cz-d2`。研究负责人于2026-10-20接受设计并附条件授权执行。当前结论：

`CZ_D2_PREFLIGHT_STOP_HISTORICAL_EXCLUSION_REGISTRY_INCOMPLETE`

## 阻断证据

Commit A的216条排除记录来自12个来源。新执行前核查发现，既有AUX-D1归档 `results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/evaluation_only/calibration_private.npz` 中还有 `h_truth`，shape=(512,4)，单位顺序为km、deg、m/s、deg。

这512个状态是历史独立阈值校准所用的生成真值，不是AUX-D1复用的12个H状态，也不是本轮新生成数据。原 `aux_d1_generate.py`、任务SPEC、1024行双资源校准分数、私有清单与THRESHOLD_FREEZE共同确定其用途与身份。归档在本设计父版本中已经存在，字节与原SHA256封存一致。

只读比较得到：512个状态均无注册表精确匹配，也均未被注册邻域（含镜像）或旧粗网格离格保护覆盖。独立标准库ZIP/NPY解析与NumPy读数一致，第二套纯Python比较得到同样的512个遗漏。证据复核55项通过、0项失败；**这不意味着历史完整性Gate通过**，而是确认遗漏确实存在。

遗漏发生在上一轮由Codex编制排除注册表时。当前注册表和设计合同均原样保留，没有擅自加入512条再继续运行。扫描在首项确定遗漏处停止，不能声称除此之外已经没有其他遗漏；下一次若获准修正，需要完成其余历史来源完整性检查。

## 为什么停止

负责人本次授权明确要求：发现历史登记不完整就STOP报审，不得自行重新选取状态。冻结 `PANEL_AND_HOLDOUT_POLICY.md` 也将额外旧生成真值列为停止条件。因此先核验历史目录这一前置条件后，即取消后续生产实现、面板物化和正式计算。没有为了绕过门槛而将校准真值排除出“已见”范围。

本轮尚未生成新面板，故未证明任何新状态与遗漏旧状态实际碰撞；失败是排除合同不完整，而非已经发生真值泄漏，更不是B3或P1的性能负结果。

## 执行数量与统计

|项目|实际状态|
|---|---|
|新面板状态/新观测|0 / 0|
|正式主实验/正式冷审计|0 / 0|
|优化/Fisher/传播求解/MC|全部0|
|计划方法配置|48000|
|已执行/未执行|0 / 48000|
|四参数P50/P90/P95/P99|NOT_EVALUATED|
|拒判率、触界率、P0/P1敏感性|NOT_EVALUATED|
|生产模型与权限Gate|未实施，不能借旧程序宣称通过|

`CONFIGURATION_EXECUTION_STATUS.csv`逐项登记48000个计划键；其中D01–D12仅为预留ID，不表示面板已经生成。未执行项不记INF、不记为求解失败，也不计成拒判。只有未来实际运行且失败的样本才适用冻结INF规则。

## 交付与下一动作

- `HISTORY_COMPLETENESS_PREFLIGHT.json`：来源字节身份、规则比较及硬门裁定依据。
- `UNREGISTERED_HISTORICAL_CALIBRATION_STATES.csv`：512个旧校准状态的取证副本，**不是修订后的排除表**。
- `preflight_history_completeness.py`：只读核验与独立NPY解析代码，不导入或运行旧生成器。
- `PREFLIGHT_EVIDENCE_REVIEW.json`：取证复核范围，明确非正式科研冷审计。
- `COMMIT_B_GATE_STATUS.csv`、`CZ_D2_B_FINAL_DECISION.json`：已通过、失败及未执行项分别记账。

本次仅提交Commit B前置STOP证据，不宣称完成精度补强。后续需要研究负责人审查这份遗漏，明确授权排除来源修订与重新预检后才能继续；不需要更换研究问题，也不应改seed、阈值或估计器来回避遗漏。未运行的正式科学调用没有被消耗，但本次STOP不能自行解除。

推送核对后STOP；不启动CZ-D3。main、CZ-D1、所有历史Gate、原8.566%条件结果和R4=0%保持不变。深度、真实接收能力与五维性能仍未评价。
