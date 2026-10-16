# ENV-CAP-NUM03：最终同步

**NUM03_REFERENCE_AND_NOMINAL_FIRST_RETURN_GEOMETRY_SUPPORTED**。本轮完成并恢复STOP。

- 父SHA：`392710ecaa69e02b6247068cb83d062a0e610995`。
- 设计SHA：`a641df8e72adf7acec4a71e295c06c78b7458079`。
- 最终SHA：本报告所在Commit B；完整身份见远端`research/env-cap-num03`。

## 分层结论

|标签|结果与作用域|
|---|---|
|NUMERICAL_REFERENCE_VALIDATED|True；全部72,018登记条目独立复核，最大距离差9.55e-09 m|
|FIRST_RETURN_GEOMETRY_SUPPORTED|9/9名义单元；返回源深平面的折射回返与折叠|
|REPRESENTATION_SENSITIVITY_CHARACTERIZED|True；C-linear/PCHIP首候选移动3.94–32.21 m，次级候选也保留|
|ENERGY_ATTRIBUTION_UNRESOLVED|保持；路径共存不自动否定几何，也未量化贡献|
|YANG_MMAR_PERFORMANCE_NOT_EVALUATED|保持；无模态提取、测深或接收性能新结论|

冷审计**73,820项，0失败，PASS=True**。另有30组60/80位参考，独立双精度与MP80最大距离差5.35e-05 m；84次登记DOP853对照，最大与Snell距离差0.00548 m。合成开发16,113项有效数值控制通过；开发测试脚本的一次语法失败和修改原因完整保留，未用9个SCS结果调试。

这次对旧失败角使用了新身份、新方法和新冻结合同，未调用或修改旧科学程序；旧6项失败保持原样，新对照见NUM03_HISTORICAL_FAILURE_RECHECK.csv。数值参考、首候选精度、拓扑和表示差异分别交付，不用一个PASS替代。

全部9单元仍只属于一个WOA气候资料来源族。新几何证据不能扩展为真实海域、所有SSP、全波传播、目标测深或工程环境鲁棒性。建议下一步只设计一次有界的传播可观测量研究；当前数值修补专项收口，不自动建设相干提供器。

主分析与冷审计各一次；累计科学壁钟468.85 s，采样RSS峰值300.0 MiB。新增外部声学求解／MC／接收SHD均0。main、ENV-CAP、REC01、CZROB-02、T2/FUNC2及R4=0%不变。

详细见[科研去留报告](NUM03_RESEARCH_DECISION.md)、[9单元敏感性表](NUM03_NINE_CASE_SENSITIVITY.csv)、[数学推导](ANALYTIC_REFERENCE_DERIVATION.md)及[独立冷审计](NUM03_INDEPENDENT_COLD_AUDIT.json)。
