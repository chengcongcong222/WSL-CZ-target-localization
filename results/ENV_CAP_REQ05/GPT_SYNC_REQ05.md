# REQ05同步报告

结论：**OBSERVATION_RESOURCE_CONDITIONS_ESTABLISHED_NO_NEW_PERFORMANCE_ADMISSION**。本轮完成观测资源论证，没有新增定位/测深精度。父提交d658f660d6bd6c9b3a585f0e5b8812c49ee8b0c4；设计SHA：80516e049cafdf10dc65896f1cd78276286c390e；最终SHA为包含本报告的Commit B，推送后另核对。分支research/env-cap-req05。

六类观测都已给出测量模型、源/通道/环境未知量、零空间及五维状态作用范围。18组无量纲代数控制和124项独立实现冷审计通过，0失败；这是公式与来源身份的核验，不是新CZ精度证书。冷审计由同一代理的独立代码实施，外部科研审核待负责人完成。

主要结论：相对强度可消去共同源级，但逐时自由差分增益仍使其零信息；互谱消去源相位也不等于通道已校准；相对时延、频率和相干孔径各有不可忽略的参考/路径自由度。辅助水平方位没有深度列。处理不能制造垂向孔径、第二节点视线或丢失的复数据。

保留两个**不同目标下并列的条件候选**：

|候选|为什么保留|必须先解决|当前边界|
|---|---|---|---|
|同步多接收深度相对强度|检验公共源级消除后是否仍有r/z信息|同步垂向数据、差分增益/AGC、噪声底及可解释功率前向|没有功率比正结果；不继承H3复响应局部倍数|
|辅助节点方位测量链|已有AUXD1/T3W1条件水平Demo|独立偏置依据、报告导航、侧向多支、同步和目标关联|不再扫基线；真实接收能力UNKNOWN，z未评价|

不设总分选优，不承诺这两个候选均能准入。其余O1单频扫描暂停，O3旧SP1/T2链不重开，O4宽带时延和O5源频/相干孔径等待具体输入条件。NUM03名义几何、RC3两航向候选、H3局部信息和全部负结果各自保留原范围。

交付重点：OBSERVATION_MODELS_AND_NULLSPACES.md；四层OBSERVATION_EVIDENCE_LAYERS.csv；30项FIVE_STATE_INFORMATION_MATRIX.csv；13项RESOURCE_NECESSARY_CONDITIONS.csv；19项PARTNER_DATA_REQUEST.csv及716交流稿；NEXT_EXPERIMENT_RECOMMENDATIONS.md。没有把必要条件变成最低设备数量或实际精度，没有向合作方发送消息。

新增声学求解0，MC0，接收数据生成0。正式代数分析和冷审计合计约0.406s，采样进程树峰值约22.8MiB；完整来源/代码/矩阵可复核。所有项目文件在D盘；历史main和R4=0%不变。两个候选均PROPOSED_NOT_FROZEN、NUMERICAL_EXECUTION_NOT_AUTHORIZED；提交推送核对后STOP。
