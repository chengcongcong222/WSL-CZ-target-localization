# CZ-D2-HISTORY-REPAIR-01 同步报告

**结论：受控修订已封存；来源完整性UNRESOLVED，正式仿真继续STOP。**

父STOP提交：`3c7291bd6d4683796deccb352d1e8c105db63d87`。原设计：`dc91982e67921eeb7a6aac4e12c482a980d5ad89`。本报告的结果提交SHA以包含本文件的Git提交为准，最终回传核对远端分支，不改main。

- V2为744条：原216条完整保留，新增AUX-D1校准512条，以及本次新发现的T2-G0静态校准16条。原表、原设计与原STOP不变。
- 冻结快照14,582份文件完成元数据及来源发现；4,818份相关外部私有/校准档案SHA核验通过。不能把发现扫描称为全量语义完整性通过。
- 独立只读预检786项通过、0失败，覆盖已登记来源身份、V2继承、新增真值、镜像/邻域/粗网格和原合同未变。它不是正式科学冷审计。
- 来源完整性未闭合：早期MFP确实使用了3组径深真值，但没有定义方位和运动状态。当前四维/静态二维排除合同没有规定仅距离真值的处理方式。未补造theta，未新增一维距离排除，也未将该来源擅自忽略。其后待裁定来源仍明确留为未完成。
- 保留H2严格浮点相等检查的1项FAIL；独立补充诊断显示最多2 ULP，已被原邻域覆盖，历史SHA一致。补充诊断曾有列优先读取错误，已保留原输出并单列勘误，未影响V2或正式实验。

**需要负责人裁定的事项**：仅r,z的历史生成数据，是否属于本轮四维“已见状态”的适用来源；若纳入，需要明确独立的一维投影排除政策。任何新规则都需另行授权。该裁定完成后仍需补齐剩余来源语义审查，不能直接恢复48,000配置。

新目标面板、随机创新、公开观测、主优化、Fisher、Monte Carlo、正式冷审计均为0。48,000配置全部NOT_EVALUATED。CZ-D3未启动；历史Gate、main、R4=0%保持。

入口：`SOURCE_DISCOVERY_AND_COVERAGE.md`、`HISTORICAL_STATE_EXCLUSIONS_V2.csv`、`SOURCE_COMPLETENESS_GAPS.json`、`INDEPENDENT_REPREFLIGHT_RESULT.json`、`REPAIR_CHANGE_AND_EXCEPTION_RECORD.md`、`REPAIR_FINAL_DECISION.json`。提交并核对科研分支后STOP。
