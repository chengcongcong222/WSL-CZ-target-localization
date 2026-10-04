# GPT 同步：R2 Gate-1

研究负责人接受 parent Gate-0 `193abb7cd13a38fac91ad1715e6873dd3d6440e9` 为 `R4_ROUTE_REDESIGN_GATE0_ACCEPTED`。本文件所在提交是本轮 Gate-1，提交后核对 local HEAD==remote/main。

本轮仅 observable admissibility + information-chain 审计：六类齐全，Class A=0、B=5、C=1；`candidate_admitted_for_quantitative_design=false`。Priority 1=CZ_ENVELOPE，Priority 2=INTERFERENCE_FRINGE，二者只获 CONDITIONAL，排序是先澄清什么，不是已发现可用距离锚。

核心缺口是从实际接收链到 feature 的来源：现有 A1 只有 bearing + 三频、每频/每窗口去均值的 relative level。原始同步 HLA 信号、dense intensity、可提取到达差、独立 Doppler frequency track 未作为该 API 输入。模型 complex pressure 不提供工程可观测 absolute propagation phase；同频阵元相位与跨频/时间相位分开。

D 盘原始项目文献已读取：Xu2024 有效原文 SHA 与旧 lock 一致，方法依赖 broadband + vertical receiver-depth slope；bottom-bounce 原文、mobile-HLA review、Yang/Liang 原文均按 signal/geometry 逐项迁移。TDOA aggregate 不需要长期 path label 的机制可以讨论，但不能解除旧 single-depth autocorr closure。fringe amplitude-only 可讨论，三个稀疏频点和 constant-β 非直接迁移限制保持；没有“fringe useless/已可用”。源级差分只消共同标量 gain，不消 per-line 源谱；Doppler 只可为运动辅约束；当前 HLA spatial range variant 拒绝。

下一阶段建议 `R2_MEASUREMENT_CHAIN_AND_SIGNAL_CONDITION_SPECIFICATION`，仍为新的负责人指令下的文档任务，未授权 quantitative design、实现或实验。R1 cold start NOT_ESTABLISHED，tracking post-acquisition 条件合理而未验证；RC3 LOCAL_CONDITIONAL_ANCHOR_ONLY。

审阅：[矩阵](R2_OBSERVABLE_ADMISSIBILITY_MATRIX.csv)、[接收链](R2_MEASUREMENT_CHAIN_AUDIT.md)、[文献迁移](R2_LITERATURE_TRANSFERABILITY.md)、[排序](R2_CANDIDATE_RANKING.md)、[报告](ROUTE_REDESIGN_GATE1_REPORT.md)、[机器判定](R2_GATE1_DECISION.json)、[完整性](VALIDATION.json)。master plan 只追加，ledger 只追加；historical R4_PROGRESS.json 保留。

新数值实验 0；R4=0%；A2/depth/SSP/P5 UNOPENED。推送后停止，等待独立审计及新的明确阶段指令。
