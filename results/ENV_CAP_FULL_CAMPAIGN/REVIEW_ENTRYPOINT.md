# ENV-CAP整体审阅入口

建议先读[最终同步](GPT_SYNC_ENV_CAP_FINAL.md)与[最终科研判定](ENV_CAP_FINAL_DECISION.json)，再读[完整科研报告](ENV_CAP_FINAL_SCIENTIFIC_REPORT.md)。本轮是负结果及解析资源分支交付，不是新CZ精度PASS。

关键核查：

1. [E1最终裁定](E1_FINAL_ADJUDICATION.json)与[冷审计](E1_INDEPENDENT_COLD_AUDIT.json)：274项中12FAIL，原9个临时标签保留但正式准入为0；[差异诊断](E1_COLD_DISCREPANCY_DIAGNOSTIC.csv)不修复FAIL。
2. [E2最终裁定](E2_FINAL_ADJUDICATION.json)：整数溢出控制失败；[源与公式身份](E2_PAPER_METHOD_IDENTITY.md)和[勘误](E2_YANG_ERRATUM_FORMULA_AUDIT.md)不借oracle条件。
3. [性能缺失表](E3_PERFORMANCE_BY_CASE.csv)：30单元、零尝试，不填零误差；[资源表](E4_RESOURCE_INVERSION.csv)不将历史分位档当ENV-CAP可行性。
4. [申请可用文字](APPLICATION_USABLE_CLAIMS.md)及[禁用表述](APPLICATION_FORBIDDEN_CLAIMS.md)，检查有没有将研究问题写成已完成能力。
5. [审计范围汇总](INDEPENDENT_AUDIT_SUMMARY.json)、[来源冷核对](E0_SOURCE_INDEPENDENT_COLD_AUDIT.json)、[E2/E4只读冷核对](E2_E4_INDEPENDENT_COLD_AUDIT.json)、[预算](BUDGET_USAGE.json)和[历史台账未改](HISTORICAL_LEDGER_PRESERVATION.json)。只读PASS不改变科学FAIL。

图表入口：[七组科研图](figures/)，PNG供浏览、PDF供排版；[图数据血缘](FIGURE_DATA_PROVENANCE.csv)与[大文件原始索引](raw_local_archive_index.csv)逐一绑定D盘输入。新科学与冷审计代码在scripts/，任务包原字节在task_package/。

科研分支research/env-cap-20261016，最终SHA由外层回报和D盘ENV_CAP_FINAL_GIT_VERIFICATION.json核对。main与R4保持原样，STOP=true。最终审核应判断失败作用域、条件证据和下一研究必要性，不要求当前失败被包装为性能通过。
