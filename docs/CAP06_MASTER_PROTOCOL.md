# CAP06 高级Codex全阶段自主研究启动协议（Git预置版）

研究背景：深海第一会聚区必须存在；目标状态 x=[r,theta,z,v,psi]。研究以单平台单拖曳水平阵为核心，不预设单阵一次解五维。海深约5 km、第一CZ约50–60 km、源收深约200 m是代表性条件，不是普适物理或甲方验收线。

## 0. 固定版本与限制

- 从REQ05最后提交 dddad1f66e6d61668ecd5af38cad3e3af9fc9245 承接；本分支 research/env-cap-cap06 是设计种子分支，不是科学执行已经开始。
- 正式新结果只写 results/ENV_CAP_CAP06/ 和独立D盘新目录。旧ENV-CAP/REC01/CZROB02/NUM03/CZO04/REQ05及T2/FUNC2、main和R4=0%全部保持，不重启已关闭的旧科学实验。
- **用户没有任何可对接合作方人员，也没有可取得的新专有资料**：不得索取/等待合作方数据，不外部联系。自主定义有依据、可复现的合成接收信号和校准误差；不能把设定数值冒称实测设备能力。
- 先读REQ05的OBSERVATION_MODELS_AND_NULLSPACES、RESOURCE_AND_RESEARCH_DECISION、NEXT_EXPERIMENT_RECOMMENDATIONS及原始冷审计，再读NUM03、CZO04、REC01、R4 H3、T3-W1关键结论。冻结输入Git blob、数据SHA、物理单位及错误使用边界。
- 本任务包目标是**真正研究非oracle测量链和有条件的信息增量**；不能把NUM03射线几何、CZO04非相干I类强度、旧T3注入的有向方位直接升级为海洋接收成功。

## 1. 证据等级和四组资源

L0解析零空间；L1公开环境来源；L2名义首CZ射线；L3同一CZ环境下数值可认证的新声学前向；L4理想信息；L5从合成接收信号非oracle提取的观测；L6冻结留出条件估计；L7真实硬件/海试本轮不预期。

S0：一个平台8元14m HLA，主阵中心深约200m。
S1：S0加四个同时不同深度接收通道，例如196/198/202/204m，需独立增益/AGC校准；不默认已有硬件。
S2：S0加一个独立8元HLA测向节点；仅旧B0 5km和B3 15.297km两种冻结几何，不再扫描更多基线。
S3：S1+S2；仅两条增量观测独立通过Gate，才允许五维互补性测试。资源成本分别列通道、阵元、节点、布放、导航、同步、通信、时频校准和处理开销。

## 2. P0 单阵基准+冻结Commit A

建立生成端truth_private、接收侧receiver_public、封存评价evaluation_private三个不互通的数据域，冻结阵形、被动源、信号带宽、频率、噪声、增益、AGC、导航、平台小转向、第一CZ条件环境、误差和关联模型。历史G/H/T3是开发材料，必须登记新的状态/噪声留出，不把匹配环境称独立海域。

零控制：平面波同轴HLA无附加r/z信息；y镜像不自动消歧；完全自由单频源级/差分通道gain无独立信息；弱场/纯噪声不能给可靠目标角；任意未知有向方位/真实偏置/真实标签绝不喂估计器。

合成程序与代码可在Commit A前利用独立解析控制开发迭代、保留修改日志。冻结所有物理/统计Gate与预算后，建立 Commit A（代码SHA、输入SHA、样本seed、状态域、拒判定义、误差阈值），此后每个SCIENCE_STARTED单元只有一次正式执行和独立冷审计。

## 3. P1 O2 同时不同接收深度相对功率比（深度/距离支线）

物理模型是相对接收功率 d_i=10log10[(I_i-N_i)/(I_ref-N_ref)]；共源发射功率在合适条件下消去，但差分gain、AGC、噪声底、源方向性、环境、r/z/运动状态的耦合仍可能吞噬信息。低能量或扣噪后非正值为缺测，不能强取log；共参考产生相关噪声，不能作独立Fisher叠加。

冻结已校准增益、未知恒定差分增益、未知线性漂移、逐历元自由gain、只有共同AGC五类控制。自由gain时额外信息必须为零；其余只有在与所有待估状态及环境干扰联合剖面后才可能有信息增量。

最关键Gate：旧I类BELLHOP非相干强度不能当作窄带稳定CW真实接收功率，H3局部复响应倍数不能借作功率比增益。若没有经独立数值检查的同一首CZ场景物理一致功率/复场前向，只能开展明确标注的解析/合成玩具处理器和信息上限，标签 O2_FIRST_CZ_FORWARD_NOT_ADMITTED；第一CZ深度精度与L5/L6保持 NOT_EVALUATED，不伪造结果。P1失败不阻断P2。

## 4. P2 O6 辅助测向接收链（水平支线，独立于P1）

从带噪阵列时序或复基带自己实现方位前端；不能把T3旧Q1/Q2有向角直接注入。先做8元14m HLA单到达平面波正控制、纯噪声零控制，再登记固定多途干扰、通道幅相误差、阵形变化、SNR、漏检、假警、AGC与导航偏置。

同节点共x-z平面的左右镜像及前后不确定性均作为多支输出，不许根据真位置/真方位符号/真目标关联预选。接收产品只有候选轴角、时标、质量标记和带噪报告导航。独立校准样本是额外资源，有无偏置先验分开；跨节点相关导航/时偏、错关联和拒判都需要统计。

B0/B3资源对照配对共用主阵的随机创新，不扩基线。信号前端若只用远场单平面波，严格标 ARRAY_FRONTEND_CONDITIONAL_NOT_CZ_MULTIPATH；该产品可以支持条件水平信息研究，但**不是**第一CZ多途声场取得的实际方位质量，不能报告首CZ波形级端到端P95。该水平模型z列为零；不能把历史T3的水平四参Demo改写成五维定位。

## 5. CP1 后P3条件融合与留出

P1/P2各交一次正式主实验、一次独立冷审计，CP1推科研分支。一个失败另一继续。仅通过接收及未知量Gate的资源配置进入P3；P1失败则只许S0/S2水平4D条件实验，深度未评；P2失败则禁止任何oracle有向角输入；两者都失败则P3不得发明性能统计。

估计器保留完整方位镜像、多途身份和离网格未决域：区别 retained/excluded/undecided，禁止看真值补格点、假定正确目标关联、用真z/e或有限Uniform偏置support冒充独立校准。生成/估计代码隔离；只用合法输入，使用联合噪声协方差。冻结独立留出，统计拒判、错支、未决、正确分支保留、无条件P95、95%事件区间和最坏状态（若MC/声场Gate不合法则统一NOT_EVALUATED）。CP2推执行结果或负结果。

## 6. P4 资源条件反推/科研结题

S0/S1/S2/S3必须公平计增量成本：新增垂向/第二节点、部署与校准、源谱稳定性、噪声相关、时钟/导航/目标关联、持续观测、阵列处理预算。精度等级（例如距离1/5/10%、深度5/10/20m、速度0.1/0.25/0.5m/s）仅为科研示例，不是甲方验收值。

逐条区别 MATHEMATICAL_NECESSARY、SIMULATION_CONDITIONAL_INFORMATION、TESTED_CONDITIONAL_FEASIBLE、TESTED_INFEASIBLE_UNDER_REGISTERED_MODEL、UNKNOWN_OR_UNTESTED。只有同级证据/有限集合可以比较Pareto，不写最低真实装备，也不把增加节点或功率归一化自身称创新。

形成完整科研报告、负结果、图及数据血缘、源码/日志/独立冷审计、设备假设/成本向量、申请可用与禁止措辞、FINAL_DECISION、RUN_STATE、GPT_SYNC_CAP06_FINAL、STOP。

## 7. Git/预算/滚动约束

科学累计壁钟≤72h，进程树RSS≤8GiB，新D盘数据≤30GiB，新Git≤120MiB，新增声学求解器≤48次，新增MC配置≤20000；这些是上限不是必须使用额度。合理更小预算必须在Commit A登记。主分析与冷审计相互独立；审计程序不可直接import主科学判定函数。

跨会话更新RUN_STATE.json+SESSION_HANDOFF.md，记录冻结SHA、每个SCIENCE_STARTED单元、输入/代码SHA、已推检查点及下一条确定性命令。已经正式启动的科学单元失败，不靠调高SNR、增加预算、换模型或新种子事后救活。缺合作方数据不能成为自动终止所有合法条件研究的理由；真正缺少第一CZ前向时只能暂停相应物理性能声称。

常规CP1和CP2推送后自动继续，只有违反硬背景、需要未登记新硬件/算法、超预算、科研造假/原历史覆盖等重大问题才请求负责人。FINAL提交、push并核对远端HEAD、STOP；禁止自动merge main，R4=0%不改变。当前文件只是设计种子，尚未开展CAP06科学实验。

## 8. 当前必读的历史入口

- results/ENV_CAP_REQ05/OBSERVATION_MODELS_AND_NULLSPACES.md
- results/ENV_CAP_REQ05/RESOURCE_AND_RESEARCH_DECISION.md
- results/ENV_CAP_REQ05/NEXT_EXPERIMENT_RECOMMENDATIONS.md
- results/ENV_CAP_CZO04/CZO04_RESEARCH_DECISION.md
- results/ENV_CAP_NUM03/NUM03_RESEARCH_DECISION.md
- results/ENV_CAP_REC01/REC01_COHERENT_FORWARD_FEASIBILITY.md
- results/R4_T3_W1_WIDE_BASELINE_DEMO/GPT_SYNC_T3_W1.md
- results/R4_H3_G0_FINITE_VERTICAL_RESPONSE/INTERPRETATION_LIMITS.md
- results/R4_ROUTE_REDESIGN_GATE_1_R2_OBSERVABLE/R2_MEASUREMENT_CHAIN_AUDIT.md

如果你有完整CAP06_FULL_CAMPAIGN_20261019_v1.zip，则把其中04—10详细合同及12_STAGE_SPEC.json一起登记进入Commit A，运行包内manifest/plan/state工具。此文件可作为在未提供ZIP时的自足启动合同；不得因ZIP缺席而任意放宽本文件的科学或审计要求。
