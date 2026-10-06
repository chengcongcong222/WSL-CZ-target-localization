# R4_OBSERVABILITY_FRONTIER_RESEARCH_RESET 交接

基线：ee656d49bb7459268a991db6c44ba50d661456db。
接受 CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB，限已保存双节点方位/1200s/当前几何噪声/正确四状态/局部Gaussian及已测估计器。26.2082%为局部1.96σ等效，不是实际P95下界。

本轮完成：历史实现覆盖清点；全文/章节/摘要阅读范围分级；未知量扣除公式；关闭路线范围复核；TOP3假设；TOP2两层最小试验设计；增量指标与成本口径。
未运行科学数值计算，未生成新MC、接收声学数据、传播分支或确认集。文档一致性校验不计科学实验。

推荐顺序：
1. E1未知源频率/漂移的双节点速度信息及接收频率提取链。
2. E2源谱剖面化的多通道差频传播表示，以真分支保留/错误分支削减为目标；原FDSL/POMAP完整对照待正文公式补齐，不把自定义处理器标作论文复现。
3. H3只保留深度区间/组统计研究假设，不打开第三项或depth执行。

主要核查结论：
- Yang2015必须连同2018勘误：offset δ影响相干深度评分；已知CW不等于未知源频。
- literature/2018.pdf实际上是Yang2014的2015勘误，490→4990m；不是新2018方法。
- Liang2018已知运动速度；AR波数/复Hankel幅度分工不能混淆。历史到Hankel桥，未完成严格Fig.5/6深度链。
- Cockrell2010原文浅海β=1、宽带、径向运动；历史解析二模态控制不等于整篇复现。
- Xu2024范围方法依赖跨接收深度斜率；Xu2023短VLA深度方法是底反射区，不可直接搬到第一CZ/HLA。
- RF2018可借鉴CFO/频率提取机制，载频/信标及主动轨迹条件不同。
- FDSL2019、POMAP2021/2023与2026 few-element FDSL目前只有摘要/出版身份核实；全文获取受限。禁止宣称已恢复其完整公式或few-element就是当前HLA。

[覆盖表](LITERATURE_IMPLEMENTATION_COVERAGE.csv)；[互补分析](OBSERVABLE_INFORMATION_COMPLEMENTARITY.md)；[路线范围](CLOSED_ROUTE_SCOPE_REVIEW.md)；[TOP3](TOP3_RESEARCH_HYPOTHESES.md)；[TOP2](TOP2_MINIMAL_EXPERIMENT_DESIGNS.md)；[指标](BREAKTHROUGH_METRICS_PROPOSAL.md)；[来源阅读记录](PRIMARY_SOURCE_READ_LOG.csv)。

双状态：
- 原联合指标：A1/R4=0%，速度未通过，历史Gate与结论不改。
- 预研交付：双节点水平基线、速度归因/效率审计、本轮文献覆盖/假设/设计完成；无虚增百分比、无新方法性能PASS。

hardware UNKNOWN不阻碍预研；time synchronization和target association仍NOT_NUMERICALLY_VALIDATED。纯方位同条件修补关闭；A2、depth执行、SSP、P5、时间延长均未开放。
本轮交付待研究负责人独立审计。提交、推送并核对HEAD=remote/main后STOP，不自动运行TOP2。
