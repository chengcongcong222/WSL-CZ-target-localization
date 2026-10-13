# T2-N3 同步报告

设计SHA：`56eb1f6acdf2f469ee08618f0cc72f3c09b5392f`。执行SHA由最终回传列示，避免文件自引用。
任务SHA256：`5082da133db2ab54de93c71cd2179f8b2df531cc5dbefbcc86335994c714e127`。

本轮判定：`T2_N3_PARTIAL_NUMERIC_NOT_ADMITTED`。

|阶段|完成范围与结论|
|---|---|
|N3-1 原物理合同|COMPLETE；8684条存档模态分类，1728条复相干带贡献记录|
|N3-2 主解析|1/2已执行配置通过；预登记12项|
|N3-2 粗网格|0/0已执行；预登记4项，零执行不称通过|
|N3-3 E-STD底界|0/144向量；NOT_EVALUATED|

原SSP逐字节核验：水深5000 m，自由海面V、刚性海底R、CVWT和实际PRT的Thorp；c_min=1500.0 m/s，c_max=1635.1811 m/s，声道轴1300 m。水柱密度1 g/cm³。所有新增1800–2400 m/s模态的相速度高于水柱最大声速，其kz²在全部原SSP节点为正。此为名义水柱无内部转向的判定，不是实测海底回波证书。

模态带在复数域完整相加并独立重建；相消和弱场均保留。各带范数之比可能大于1，不能称为独立声线能量比例。

独立解析参考已复现并冷核对：16条计数、240条窗口差异、80条Hankel/远场差异。它们描述两个等声速理想波导，不是原Munk型SSP或真实第一会聚区。

已执行KRAKEN解析场最大空间形状差异：0.22149397%（门槛0.2%），最大水平波数相对偏差5.2303846e-08（门槛2e-5）。最坏：n3_PR_R_f201_h2400_n160001、P04、B1_V4。逐配置模态数量、原始复场差、全五探针四资源共用一个复比例的诊断全部保存。

冻结停止原因：`[{"call_order": 2, "stem": "n3_PR_R_f201_h2400_n160001", "reason": "ANALYTIC_HARD_GATE_FAIL", "future_calls": "CANCELLED_NO_RETRY"}]`。没有重试、修改环境、增加预算或改门槛；后续配置明确标记NOT_EVALUATED。

冷复核：22649项，0 FAIL，PASS=True。使用独立MOD解析、本征根brentq、J0-iY0径向项及QR投影。冷复核通过只认证结果和Gate记录可重建，不能改变科学Gate失败。

新增KRAKEN 2/19；SCOOTER/FIELDS/FIELD、波形、DFT、MC、FIM、优化、N1模板重评分均0。求解期壁钟29.61 s，进程树峰RSS 68.02 MiB；科学全链预算另见PIPELINE_BUDGET.json。

原G0/N1/N2全部冻结；N1阵型证据仅限旧1800窗口，不移植到本轮。CHIGH全域完整性未认证。T2-A/B接收和定位性能NOT_EVALUATED，真实CZ1传播NOT_CERTIFIED，实际底质UNKNOWN，R4=0%。

后续建议：若解析Gate失败，先研究本次已保存的模型/工具数值合同差异；若全部通过且底界敏感，则须获得有来源的底质与传播合同，再讨论阵型优选。任何后续工作均需独立任务，本轮不自动执行。

负责人审计PENDING。STOP。

工具边界定义参考：[官方KRAKEN手册](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html)、[边界条件](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node10.html)。

交付补充：40条已执行解析场比较中36条通过、4条超限，均为CHIGH=2400的P04，分别对应A0_0、A_P4M、B1_V4、B2_DUP；最大0.22149397%发生于B1_V4。两次程序调用均成功且无Thorp、Im(kr)=0；硬门失败仅指本次解析场精度未达到冻结阈值，不能扩大为KRAKEN实现错误或物理模型不存在。

辅助诊断作用域：ANALYTIC_MODAL_ERRORS.csv中的单模态函数绝对差采用共同复比例对齐，包含幅度因子，不能作为仅相位校准后的归一化函数精度证书。该因子从未进入主复场、波数或空间形状Gate；原始MOD和全部因子保留。本轮保留这一辅助口径限制，不改冻结代码或主科学结果。

完整运行与冷复核壁钟55.203秒，进程树峰RSS约125.58 MiB，均在冻结预算内。LOCAL_ARCHIVE_INDEX.csv是报告生成时快照，其中报告进程日志随后完成写入；最终字节身份以FINAL_LOCAL_ARCHIVE_INDEX.csv为准，该输出索引更新不涉及科学结果、求解器重试或Gate修改。