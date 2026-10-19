# 716第一会聚区单拖曳阵项目：三轮阶段交付总控与CZ-D1 Codex任务书

版本：2026-10-19 / DRAFT FOR RESEARCH-LEAD AUDIT

## 0. 任务定位

**编号规范**：本交付包三轮使用`CZ-D1`、`CZ-D2`、`CZ-D3`，以免与仓库历史研究阶段`R3/R4`或路线候选`R1/R2/R3`发生冲突。

**交付目标**：把既有研究结果转换为“在什么观测资源、测量精度、环境/信号前提下，目标五维状态可达到什么输出精度或受限程度”的可追溯证据体系，支持短期向甲方提交条件明确的阶段性技术报告。五维状态为 `x=[r, theta, z, v, psi]`（距离、方位、深度、速度、航向）。不以“114576→2候选”代替实际连续状态精度。

**主研究对象**：深海第一会聚区，目标约50–60 km、目标与接收阵深度约200 m，单平台单水平拖曳阵，典型8元、2 m元间距、14 m物理孔径。机动、补充声学观测、垂向扩展、倾斜阵、辅助节点均作为资源增强条件，分别计账；不得混称原单HLA性能。

**三轮安排**：
1. **CZ-D1（证据归并，不新增物理仿真）**：搜集全项目所有已执行/未执行/关闭实验，编制输入观测精度—五参数输出精度证据主表，覆盖机动、辅助声学观测、垂向观测及辅助节点；记录结果级别、适用范围和失败/证据缺口。
2. **CZ-D2（针对性精度补强，CZ-D1审计后另行授权）**：按CZ-D1预定义选择规则确定最接近交付的主方案，优先检验离网格性能、拒判与误差尾部，预冻结实验、统计及成功失败Gate。
3. **CZ-D3（甲方交付封版）**：出具技术报告、精度—观测条件对照表、方案与必要资源、科研图件和汇报材料；不得把条件性仿真宣传为设备能力或海试精度。

## 1. CZ-D1观测架构分类：必须覆盖

|ID|观测架构|归并对象与比较注意事项|
|---|---|---|
|A0|原单HLA基线|R3/R4-A1/H1/H2/SP1等；原8元×2 m，保留原实验各自真实配置|
|A1|原单HLA+正常小幅机动|0°/2°/5°/10°/15°；区分RC2-only与RC2+RC3；多数增强组自身已经含600 s后15°转向，不能据此推断纯机动增益|
|A2|原单HLA+附加声学特征|多频会聚区包络、强度/干涉、径向信息、特征提取（RX/CAP06）等；严格区分观测值注入和真实接收链提取|
|A3|HLA+短垂向扩展|H3局部信息、模态数值与接收提取；新增垂向通道及跨深度同步/增益校准代价|
|A4|倾斜阵或其他垂向几何|T2-TILT；若前向数值Gate未通过，则主性能必须NOT_EVALUATED，不能写成精度为0或阵型失败|
|A5|HLA+辅助节点|AUX-D1、T3-W1（B0/B1/B2/B3）以及辅助节点误差预算；B0已经是双节点配置，绝非纯单HLA|
|A6|已实际测试的联合增强|只有同一个冻结实验真正联合施加并可分离的配置才能填写；其他组合标NOT_EVALUATED，不拼接独立实验数值|

每个实验归并时必须记录是**同数据配对对照**还是不同实验非可比结果，是否满足单因素因果判别；未完成正反方位假设检验时，保留左右舷模糊。

## 2. 指标定义（所有表格共同遵守）

- 输入：海况/SSP、水深和传感器/源深、工作频带/源谱S0/S1/S2、阵元数/间距/新增通道或节点、机动量与观测时长、采样历元、测向精度（随机/系统偏差分开）、导航误差、校准误差、相干/同步/关联前提、声学特征偏差与提取成功率、接收数据来源（理想注入/合成波形/真实录音）、oracle权限。
- 输出：距离`|rhat-r|/r`（%）；方位`wrap(theta_hat-theta)`最小角距（°）；深度`|zhat-z|`（m）；速度同时报绝对误差（m/s）和相对误差（%）；航向最小圆周角距（°）。另列**五维联合成功率**，前提是甲方明确联合标准或预冻结逐参数阈值；没有定义则`NOT_DEFINED`。
- 统计：样本数、独立基础创新数、复用/配对结构、分场景及最坏场景P50/P90/P95/P99、各参数条件/无条件分位、拒判率、估计输出率、正确候选覆盖与全域捕获状态。失败或拒判的无条件精度按冻结处理规则（通常INF），不能仅报告成功条件P95。离网格/网格内应当分开。
- 证据类别：`CONTINUOUS_ESTIMATOR_ERROR`、`GRID_CANDIDATE_CONTRACTION`、`FISHER_LOCAL_INFORMATION`、`INJECTED_FEATURE_VALUE`、`RECEIVER_EXTRACTION`、`NOT_EVALUATED`。不同类别的“数值好坏”不得混作同一精度比较。
- 核实状态：`LEAD_AUDITED`、`COLD_AUDIT_ONLY`、`PENDING_LEAD_AUDIT`、`FAIL_GATE`、`SOURCE_MISSING`、`NOT_EVALUATED`。先后状态不能静默升级。
- 输出宣称级别：`MECHANISM`→`IDEAL_CONDITIONAL`→`RECEIVER_CHAIN_CONDITIONAL`→`ENVIRONMENT_STRESSED`→`ENGINEERING/SEA_TRIAL`，后一级必须有新增证据，不得从较低级别直接推断。
- 需要在原始记录中区分：**物理不可辨识、数值求解失败、冻结Gate失败、真实观测缺失、设计未执行**；不能把这些都写成“方法不可行”。

## 3. CZ-D1证据总表的必需字段

建议建立 `CZ_D1_EVIDENCE_MASTER`，一实验/质量档/架构对照/统计口径为一行，至少包含：

`experiment_id, phase, architecture_id, version_sha, evidence_path, result_file, scenario, r_true_km, theta_true_deg, z_true_m, v_true_mps, psi_true_deg, N_elements, element_spacing_m, vertical_channels, aux_count, aux_baseline_km, turn_deg, turn_start_s, observation_s, epoch_count, signal_S_class, frequencies_Hz, sigma_bearing_deg, bearing_bias_bound_deg, navigation_sigma_m, channel_gain_contract, feature_sigma_and_unit, extraction_source, environmental_contract, matched_or_mismatched, oracle_flags, estimator_family, offgrid_yes_no, n_total, n_independent, n_rejected, rejection_rule, metrics_scope, r_P50_percent, r_P95_percent, theta_P95_deg, z_P95_m, v_P50_percent, v_P95_percent, v_P95_mps, psi_P95_deg, joint_status, correct_support_retention, cold_audit_status, lead_audit_status, gate, claim_scope, missing_inputs, comments`。

数字单元只放能从原始输出重算的**同一定义**数字；无法证明的数值留空且配套状态字段，禁止填0或编造区间。需另外提交中文简明“甲方版”精度—条件表与技术版全量证据台账。

## 4. 已核对的关键仓库入口（查最新父SHA后引用准确固定版本）

仓库：`chengcongcong222/WSL-CZ-target-localization`。在2026-10-19一次只读核查中，远程最新可见`main`为`a9a80191ee4d0f0cecd9db296198606b63a68080`，执行时必须重新核对，而不是静默认定它仍是HEAD。

- R3工作包络：`results/R3_FINAL_CLOSEOUT/R3_FINAL_WORKING_ENVELOPE.md`。114,576→2仅为理想匹配、冻结网格候选筛选，非连续精度。
- R3机动：`results/R3_RC23_CLOSEDLOOP/R3_RC23_SINGLE_TURN_RANGE_ANCHOR/GPT_SYNC.md` 和 `.../R3_RC23_TURN_FREQ_SUBSET_FIX2/GPT_SYNC.md`。RC3辅助下有条件收缩，RC2-only远距冷启动结论不得被覆盖。
- 单HLA离网格负例：`results/R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_REPORT.md`（0/27粗网格定位成功的受限试验）。
- H2径向特征：`results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/GPT_SYNC.md`。主U档速度等权平均中位0.7538%，最坏几何无条件P95为INF（包含拒判），条件尾部仍差；特征为假设注入而非实际声学提取。
- H3垂向：`results/R4_H3_G0_FINITE_VERTICAL_RESPONSE/GPT_SYNC.md`，后续数值/SSP审计`results/R4_H3_P1_UNIFIED_WINDOW_SSP/GPT_SYNC.md`。2.34–4.47倍只为局部深度信息增益，非深度估计误差倍数。
- 倾斜阵：`results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/GPT_SYNC.md`，冻结前向Gate失败，测试性能NOT_EVALUATED。
- 辅助节点：`results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/GPT_SYNC.md`。U0→U5的条件四参数收益，速度最坏P95仍31.1631%。
- 宽基线：`results/R4_T3_W1_WIDE_BASELINE_DEMO/GPT_SYNC_T3_W1.md`。Q1/B3: 距离P95 0.7026%、方位0.1413°、速度8.5659%、航向1.2215°（均最坏场景），深度NOT_EVALUATED；B3是15.297 km级双节点，不属于原单HLA；48,000个配置重用6000组基础创新，不能当48,000个独立观测样本。
- CAP06等部分最新本地成果：**远程当前未找到同名CAP06完整证据文件**。必须查Codex本地工作目录及数据盘，确认文件、SHA与原始输出后才能升级证据；不存在的或不可访问的文件标SOURCE_MISSING，不得按口头数字替代独立审计。

## 5. CZ-D1需回答的三项决策问题

**Q1：原单HLA可独立承诺的参数精度分别是什么？** 没有连续离网格/可提取接收链支持的，应如实标`NOT_ESTABLISHED`，并列出最接近的条件性结果与主要阻塞点。

**Q2：小幅机动、附加声学特征、垂向扩展、辅助节点分别改善了什么？** 用同合同配对差值说话，分别列出对r/theta/z/v/psi的P95增量；若只有候选集、信息量或不同配置间不可比结果，报告对应本来量纲与需要的对照实验。

**Q3：哪个技术方案最适合CZ-D2短期精度补强？** 按甲方原单HLA主线一致性、已有完整接收—估计链、证据可靠性、离网格可靠性、错误尾部、可实现成本、补强计算量综合排序。可以为首交付保留一个增强架构的“条件工程参考”，但不能用其掩盖原单阵证据缺口。排序需写正反证据，结果不足可不选赢家。

## 6. CZ-D1正式Codex执行指令

1. 只读取既有研究仓库和授权本地原始档案，先核对`HEAD`、`remote/main`及数据目录，形成完整Stage/Commit/Status清单，**不重写历史结果**。
2. 逐份阅读重要`DESIGN_FREEZE`、`GPT_SYNC`、`DECISION`、冻结统计CSV及原始数据哈希，把不同资源/状态/统计定义如实写入CZ-D1主表；系统梳理机动、附加观测、垂向、倾斜及辅助节点。不能从Git文档“通过”字样推断科学已审计。
3. 检查实验输入精度与五维输出精度是否真的同实验闭合。遇到缺少原始数据、标注含混或统计不可比，记录问题和精确路径，禁止静默补值。
4. 尽量使用已存CSV执行**只读复算**（合法访问的情况下），重点复核R3、R4-A1、H2、H3、AUX-D1、T3-W1的样本分母、P95、INF、配对关系和适用边界。只读复算不等同新模型仿真。
5. 交付：`CZ_D1_EVIDENCE_MASTER.csv`、`CZ_D1_ARCHITECTURE_METRICS.md`、`CZ_D1_MANEUVERS_AND_AUGMENTATIONS.md`、`CZ_D1_FIVE_STATE_CLAIM_MATRIX.md`、`CZ_D1_NEGATIVE_GATES_AND_GAPS.md`、`CZ_D2_SELECTION_CANDIDATES.md`、`CZ_D1_SOURCE_PROVENANCE.json`、`GPT_SYNC_CZ_D1.md`。如生成甲方简表，同步保留完整溯源路径。
6. 仅在新CZ_D1目录内生成必要脚本与报告，不修改历史冻结目录。新数值仿真、KRAKEN/FIELD、估计器再优化、采样新噪声、重新选择旧优胜样本、扩大阵列参数扫描一律禁止。未做的测试写`NOT_EVALUATED`。
7. 输出CZ-D1 Gate：`EVIDENCE_INVENTORY_COMPLETE`、`METRIC_DEFINITIONS_RECONCILED`、`SOURCE_GAPS_DISCLOSED`、`CZ_D2_CANDIDATE_RANKING_READY`逐项说明是否通过。允许有合理的`PENDING`或`NO_CANDIDATE_ADMITTED`，不得为了关Gate而“补出”精度。
8. 本轮结束提交Git并推送，报告真实parent/design/result SHA和remote/main核对结果，之后STOP。**不自动开启CZ-D2**，等待GPT独立研究审计与负责人授权。

## 7. CZ-D2、CZ-D3依赖条件

CZ-D2只有在上述主表溯源和独立审计完成后才能另行编写冻结实验；优先离网格、拒判、无条件P95/尾部与左右舷，不在未知接收可提取性下借助虚构特征。若单阵完全没有可短期闭合的精度链，可将其作为“研究边界”，同时交付双节点条件性四参数示例，不能宣称五维成功。

CZ-D3对外报告只发布达到对应证据级别的指标；逐参数约束和统计定义明确。报告中须并列展示原单HLA与增强架构，标明设备资源、运行条件、理论上限、统计范围及需要甲方后续确认的输入信息。未经确认不发布“五维联合精度<10%”或“实际海试达到该指标”。

**CZ-D1到此为止；其输出需独立审计后决定第二轮方案。**
