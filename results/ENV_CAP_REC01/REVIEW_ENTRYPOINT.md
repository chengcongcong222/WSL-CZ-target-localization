# REC01审阅入口

先读GPT_SYNC_ENV_CAP_REC01.md，再读REC01_FINAL_DECISION.json。代码、规则和输入身份锁定在Commit A，正式主分析与独立冷审计均各一次；执行记录和局部原始记录路径见RUN_STATE.json与RAW_DATA_INDEX.csv。

关键核查：31,716个全域能量窗是否完整；旧12项失败是否全部映射；9个条件单元是否仍只算一个资料族；边界路径未分离的限制是否明确；R2新通过是否被错误升级为旧E2通过；相干提供器未准入时是否停止Yang/MMAR性能。FIGURE_DATA_PROVENANCE.csv给出科研图的数据来源。

main、旧ENV-CAP、T2/FUNC2及R4台账不改。本分支STOP，不自动合并main。
