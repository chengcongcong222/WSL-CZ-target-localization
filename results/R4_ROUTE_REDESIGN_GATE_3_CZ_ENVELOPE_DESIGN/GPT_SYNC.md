# GPT sync：CZ envelope Gate-3设计

parent Gate-2 `dcad14337491ac5bf90acacf0be34a57eafac3f1` 已被负责人接受。Contract A accepted for quantitative Gate design ONLY；fringe未准入，M1 receiver/UUV真实条件未验证。本轮交付是draft，不是pre-run freeze或run许可。

候选E1 WINDOW_SHAPE_Q2、E2 LINE_BLOCK_MEANS_B3、E3 NORMALIZED_INTENSITY_BLOCKS_B3均M0-compatible，定义从三频/121样本/旧两window relative level出发。E1推荐只因静态offset invariance、最低trend+curvature与可审计drift投影；没有计算feature或选择性能最佳者。block尺度没有物理充分性证明，显式SCALE_SELECTION_NOT_YET_JUSTIFIED。

N0/N1/N2分开；N1 bound未冻结，N2任意源级可吸收所有传播信息的limit只做符号推论。tau_F需noise/numeric元数据，禁止用truth或旧J门限替代。输出joint support与全部range intervals；保守unknown cells计入width，不能只用point cloud/top1。full-cell enclosure、grid、hard预算均待pre-run选择。

建议nominal：所有主casejoint truth/reference保留，W_hull≤5km/C≤1/3，rangecomp≤3/joint broad comp≤6，双grid端点/width drift≤0.5km并保持预注册coarse拓扑；lenient7.5km，strict3km，全部PROPOSED_NOT_FROZEN。没有RC3capture保证，非正式<10%指标。

旧21noisy P/Q主development，15noiseless仅regression；旧FIX2 24case全部归已见development。fresh NOT_GENERATED；不调整历史dataset/threshold/code。停止规则覆盖丢truth、无range缩域、alias/search未闭合、truthmotion依赖与nuisance吸收。

建议下一阶段仅 `RESEARCH_LEAD_GATE3_DESIGN_REVIEW_AND_PRE_RUN_FREEZE`。quantitative_run_authorized=false，新数值实验0，R4=0%，A2/depth/SSP/P5 UNOPENED。只append master，推送后停止。

阅读 [design report](ROUTE_REDESIGN_GATE3_DESIGN_REPORT.md)、[representation definitions](CZ_ENVELOPE_REPRESENTATION_CANDIDATES.md)、[threshold proposals](CZ_ENVELOPE_THRESHOLD_PROPOSALS.md)、[machine decision](R2_GATE3_DESIGN_DECISION.json)、[integrity validation](VALIDATION.json)。
