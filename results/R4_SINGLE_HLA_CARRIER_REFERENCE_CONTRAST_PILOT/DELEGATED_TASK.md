# HLA-RX1Z：已存接收记录的频率参考与跨窗相对运动量小试

版本：2026-10-10 / v1.0
性质：APPLICATION_PRE_RESEARCH / SAVED_DATA_NEW_REPRESENTATION_PILOT
任务标识：HLA_RX1Z_CARRIER_LEDGER_AND_INTERWINDOW_CONTRAST

## 0. 这是一项新任务，不是重发、续跑或重启RX1

仓库：chengcongcong222/WSL-CZ-target-localization
唯一parent：41463d9aa71da1ae6ebed9b025de22cb3d679e68
RX1设计：2873786350c23cd5bc90bfb5808c3dbbd2ed728c
新结果目录：results/R4_SINGLE_HLA_CARRIER_REFERENCE_CONTRAST_PILOT/
新代码前缀：hla_rx1z_

旧任务HLA_RX1_CODEX_TASK_20261010.md的SHA256：
0d4a7841d5cbaee946dec6131e5efef6342addc06c145f37d78aa5552edd8ea8

上述旧任务已经在parent提交结果。不得再运行hla_rx1_execute.py，不得补跑DRIFT，不得把本任务当作旧任务恢复许可。用户转交本文件后，只授权本文件的有限存档处理范围。

旧结论保持：FEATURE_MAPPING_NOT_ESTABLISHED、PARTIAL_EXECUTION、IMPLEMENTATION_OR_PROVIDER_LIMIT；STABLE 128条已生成，DRIFT 128条未生成。STABLE使用1500 m/s群延迟后备假设，动态模态物理映射未认证。新分析不能替旧提供器补发认证。

原单HLA资源不变；不增加垂向通道或辅助节点，不研究深度收紧。H1/H2/E1及历史H3结果、原R4=0%全部保留。

## 1. 新问题与交付目标

RX1在每个200 s窗独立估计接收谱线中心offset，然后乘exp(-i*2*pi*offset*t)，最后从该残余序列的相关谱峰映射q_eff。源频偏与平均运动频移可能同时被移除。旧结果已经提供了每窗f_rx_hz和rx_center_offset_hz。

本轮只问：
1. 逐窗重定心究竟保留/消除了什么？缺失的是绝对速度零点，还是连跨窗变化也被后续表示忽略？
2. 已保存的、从接收侧得到的六窗载频序列，其跨窗差值是否与某种相对径向变化有可重复的对应？
3. 对应即使成立，是否主要只反映已知平台转向，而非对目标速度有额外区分？

本轮最大声明是“在RX1存档信号合同下，有/没有可利用的跨窗接收频率变化，以及其几何映射诊断”。绝不声称恢复绝对径向幅值、未知发射频率已经解决、第一CZ接收提取已经通过，或H2速度精度已经实现。

预期交付必须有实际频率轨迹、差分与误差/信号尺度表，而不只是一次文档审计和大量PASS计数。

## 2. 与RX1的实质区别

| 项目 | 旧RX1 | 本轮RX1Z |
|---|---|---|
| 输入 | 新生成接收记录 | 只用128条已存STABLE记录及其接收侧频率输出 |
| 主表示 | 逐窗去接收中心后取ACF残余峰 | 保留实际接收中心频率，取五个相邻窗差分 |
| 目标 | 径向幅值提取与H2映射 | 频率参考信息归属、相对变化映射与收益边界 |
| 新传播/信号 | 曾有30次KRAKEN和128条信号 | 一律0 |
| DRIFT | 未生成 | 继续NOT_EVALUATED，不补齐 |
| H2回灌 | 未准入 | 本轮仍不做 |

本轮不能只是删除DC峰、换主峰或改相关窗来“救”旧q_eff，也不能把旧REF结果借给接收侧新表示。

## 3. 只读材料与输入可用性

必读并登记SHA：
- 旧目录的GPT_SYNC.md、DECISION.json、DESIGN_FREEZE.json、METHOD_INPUT_LOCK_READABLE.md；
- hla_rx1_extract.py、hla_rx1_generate.py；
- 旧目录的EXTRACTED_FEATURES.csv、RECORD_EXECUTION.csv、RECEIVER_DATA_MANIFEST.json、WINDOW_FEATURE_MAPPING.csv、MEASUREMENT_CONDITION_TABLE_READABLE.md；
- 原本机任务SHA与旧目录DELEGATED_TASK.md；
- H2的MEASUREMENT_TO_SIGNAL_GAP.md和FOUR_PARAMETER_RESULTS.md，仅用于边界对照。

读取路径和实际列名必须从真实文件核对，不猜字段语义。原始接收记录路径以manifest为准，当前预期在D:/ProjectStorage/WSL-CZ/HLA_RX1/receiver_records。不得假装当前进程拥有不存在的D盘文件。

若CSV是LFS指针、编码不同或读为空，先检查真实字节/编码和仓库实际对象；不得把“读取为空”当作科学空结果。允许读取本机同SHA的原始文件，不生成替代观测。报告工作树换行与Git blob字节差异时分别登记，不改旧文件。

需要的数据不存在或SHA不符时，先保留已确认的来源和缺口；本轮不得调用旧生成器恢复它们。若只能读取完整EXTRACTED_FEATURES而没有原始记录，可完成存档频率差分，但把原始样本复核标MISSING_RAW_RECORDS，不能声称已完成接收链复核。

## 4. 推导与信息权限（先做，篇幅不超过两页）

在单有效分量简化模型中写出：
x(t)=A exp{i[2*pi*(F0-F_LO)-k*q]*t}。
若完全用接收估计中心去旋转，则恒定运动对应的频移与F0-F_LO一起被移去，残余可为常数。该例只证明此模型及处理的零点歧义，不是全部多模态/多通道方法的No-Go。

必须区分四项：
- 已知接收器LO/数字参考（合法公开元数据）；
- 从接收数据估计的谱线中心（混合源与运动）；
- 未知发射频率/源跨时相位（主处理不可读）；
- 名义单分量中的运动频移（不是直接观测）。

用代数说明公共数字参考改变不应改变“恢复到同一物理参考后的接收频率”；但已知LO不能代替已知F0。不能“把去掉的频率加回来”就宣布知道绝对q。

对稳定未知源的一阶单分量关系，仅作为映射候选：
F_rx,j ~= F0*(1-qbar_j/c).
因此Delta F_rx ~= -(F0/c)*Delta qbar。
这里有未知比例、传播相速度、模态与时间权问题；不得将其当作已建立的第一CZ真公式。

若源频率随时间漂移，上式还会多出源变化项。仅做符号推导说明混淆；本轮没有DRIFT接收数据，不能据此给出DRIFT的数值精度。

## 5. 接收侧主表示：先冻结输出，再读取评价真值

### 5.1 六窗频率账本

接收处理进程只读public receiver metadata和接收频率记录/波形，不读private_source_frequency、true q、运动状态、源深、真实相位或模态身份。record_id仅用于索引。

从EXTRACTED_FEATURES中取method=RX、lag_s=200的记录，每个record/line/window只保留一个原始载频估计：
Fhat[l,j] = public_LO[l] + saved_center_offset[l,j]。
核对它与保存f_rx_hz一致。100 s的重复载频记录是同一估计的重复展示，不作为额外样本；REF行不进入新主表示。

每记录应有3条线x6窗=18个载频值。先构造RECEIVER_FREQUENCY_LEDGER.csv，记录原时间支持、原检测比、原主要峰及质量，不改谱峰选择。不根据真值、REF或预期运动选另一条峰。

载频有效条件只用“数值有限且原检测比>=6”；该阈值沿旧检测规则，不调参。旧残余相关accepted标记单独保留，不将它自动等同于载频准确或不准确。

如果所有输入齐全，本轮有128条STABLE记录（20dB/5dB各64条），2304个线窗载频记录。原计划256条接收记录的另外128条仍未评价。

### 5.2 唯一新主特征

每条频率线取五个相邻窗差分：
d[l,j] = Fhat[l,j+1] - Fhat[l,j], j=0..4。
主单位Hz。两端均有有效载频才输出；其余明确NOT_AVAILABLE，保留完整计划分母。不得用插值或真值补窗。

同时给出辅助“表观径向变化标度”：
Fbar[l] = mean_j Fhat[l,j]（要求六窗都有有效值）；
v_eq[l,j] = -c_ref*d[l,j]/Fbar[l]。
Fbar来自同一接收记录，c_ref沿旧公开元数据。v_eq只叫APPARENT_RADIAL_CHANGE，不叫绝对径向速度，不叫H2径向幅值，不用F0代替Fbar。若六窗不齐，Hz差分可以仍输出，标度转换不补值。

保留原符号的频率差。频率升降本身可观测，但其符号不等于目标绝对接近/远离；不得据它给H2的六个绝对q指定符号。

本轮不新增学习器、拟合偏置、运动模型约束或多峰跟踪优化，不把RF/声学多普勒公式强行套到多模态相关峰。

### 5.3 复算与数字参考检查

若128条原记录可用，只做一次同设置的载频复算。继承8s Hann/1s步长和原固定合并，不改变采样、平滑、频带、峰检、窗长；比较保存的频率中心。可用隔离的新纯函数复现，禁止导入会运行旧生成/求解器或写旧目录的脚本。

在按record_id预先排序的首条20dB记录和首条5dB记录上，做纯数字参考变换检查：各通道乘exp(-i*2*pi*delta_LO*t)，同时更新LO；delta_LO固定0.031Hz。将输入恢复到旧公开物理参考后再使用同一合并和滤波。不得因为LO改了就换物理波束或频带。该检查是同一物理记录的坐标变换，不是新声源数据、不是源频率被识别的证据。

若复算差异影响主结论，记录IMPLEMENTATION_LIMIT并停止相关结论；不当场换峰算法重跑主面板。

以上接收侧输出完成并写SHA后，才开放评价侧读取真值文件。无需Git再多插入一轮用户确认；用输出manifest证明顺序。

## 6. 评价：诊断是否有相对运动信息，尤其避免只测到平台转向

评价器可读旧几何真值，但不能把其结果传回特征处理器。

主要比较对象：H2定义六窗qbar形成的五个有符号差分D*qbar，而不是abs(qbar)，也不是六窗qbar本身。辅助比较对象是旧实际STFT样本时间支持上算出的对应平均径向变化；两种对象分别报告，不能选择较好者作为主值。

必须列出STFT实际time_start/time_end（首尾裁切等），未证明相关/谱峰的时间权时不要称为精确映射。表格标题统一使用MAPPING_DIAGNOSTIC_UNDER_RX1_FALLBACK_MODEL。

对每几何/SNR/频带分别给出：
- 计划数、有效数、缺失数；
- 原频率中心轨迹和五个Hz差分；
- v_eq与D*qbar的偏差、MAE/RMSE、绝对误差中位数/最大值；
- 真实变化的幅度尺度，以及零变化预测器的绝对误差，判断“小误差”是否仅因目标量接近零；
- 有效输出条件统计与缺失计INF的完整分母统计；
- 所有线窗和失败留存，不用只选最好的频率。

同一窗口参与相邻两个差分，误差天然共享；同源三频与同一记录的五个差分都不是独立航迹。每几何/SNR只有8个实现，尾分位仅描述，不承诺P95概率。不得把差分误差直接称为H2的独立Gaussian sigma_q。

转向分组必须预先固定：六窗编号0..5，差分j=2跨过600s转向；j=0,1,3,4在同一平台航段内。分别报告五项合计、单独跨转向项、其余四项。不允许只展示跨转向的大变化而称为目标速度精确可估。

仅作评价侧的代数解释：q=e(t)^T(v_target-v_platform)。观测到跨转向变化可能主要来自已知平台运动，缺失绝对常数后目标速度未必得到强约束。本轮不拟合H2状态，不用oracle方位扣掉平台项后把余量冒称接收观测。

三频融合如需展示，只能附加固定等权中位数，并保留所有逐线结果；它不是主判据，不根据数据选择“好线”。若要计算更复杂融合，留待下阶段，不能本轮加分。

## 7. 路线输出（不要以检查计数替代科学结论）

按实际结果给三个相互独立的状态字段：

A. REFERENCE_INFORMATION_ACCOUNTED：若单窗重定心的信息含义、公开参考与源参考权限及原记录复算一致，则该项成立。这不依赖恢复q成功。

B. INTERWINDOW_FREQUENCY_FEATURE：状态为REPRODUCIBLE / TRACKING_UNSTABLE / DATA_INCOMPLETE。记录接收载频跨窗变化能否复现、是否有跳峰或近零误判。跳峰不靠事后手工换峰消除。

C. RELATIVE_MOTION_MAPPING：最多为SUPPORTED_IN_SAVED_MODEL_DIAGNOSTIC / NOT_SUPPORTED_IN_SAVED_MODEL / UNRESOLVED。主对照STABLE20dB，5dB只作压力对照。结果应比较真实变化与误差、跨转向与同航段，以及各几何一致性，不预设必须达到某个H2误差档。小误差但比零变化预测无增量时，不判为有用的运动映射。

所有结果并列强制保留：
- ABSOLUTE_RADIAL_ZERO_NOT_ESTABLISHED；
- FIRST_CZ_DYNAMIC_PROVIDER_UNCERTIFIED_FROM_RX1；
- DRIFT_NOT_EVALUATED；
- H2_FEEDBACK_NOT_AUTHORIZED；
- R4=0%。

若相对量诊断有价值，只推荐后续单独研究它对运动参数的有效增量，不直接赋予H2收益。若只有跨转向项有信号或载频跳峰主导，停止这条处理版本，保留限制，不在本轮再换多普勒、谱峰或滤波器。

## 8. 最小交付

1. NEW_TASK_AND_PARENT.json：本任务SHA、旧任务SHA、唯一parent、差异清单；
2. DESIGN_FREEZE.json；SOURCE_AND_PERMISSION_LEDGER.md（<=2页）；
3. RECEIVER_FREQUENCY_LEDGER.csv；INTERWINDOW_FEATURES.csv；FEATURE_OUTPUT_MANIFEST.json（在评价真值读取前写）；
4. RELATIVE_MAPPING_BY_GEOMETRY.csv；TURN_VS_SAME_LEG.csv；
5. 最多两张结果图：预先按ID选定记录的载频轨迹/旧残余峰对照；全几何相对变化与误差尺度，不只画成功案例；
6. VALIDATION.json；DECISION.json；OUTPUT_MANIFEST.json；
7. GPT_SYNC.md（约1000中文字符）：旧RX1的事实、本轮新增计算、频率变化是否有用、哪种量测仍未知、能否接H2，以及完整计数和SHA。

任务与报告均用明确UTF-8写入；不要经PowerShell默认ASCII管道创建中文。启动提示词含本任务SHA，必须核对后执行。本轮来源仅以上存档和必要代数推导，不做新一轮文献普查。

## 9. 执行预算与防重复

预算：单计算进程、BLAS<=4线程、峰值内存<=4GiB、总计算/复核wall<=1800s、新增仓库<=64MiB。预算是硬上限，不是预计耗时。

新KRAKEN/FIELD/传播求解=0；新声源/接收记录=0；新噪声抽样=0；新状态优化=0。仅允许既有记录的确定性重表示和标量几何评价。不得调用旧RX1/H1/H2入口，禁止恢复缺失的DRIFT。

流程：
1. HEAD与remote/main必须均为parent。若不同，读取新增提交作用域后报告，不擅自合并、回滚或沿旧parent执行。
2. 核对本任务新SHA与旧SHA不同；检查新目录无已完成同SHA执行记录。若确已执行该新任务，给出真实执行SHA并停止。
3. 新模块预检、固定输入索引/权限/评价定义，Commit A推送。此次转交任务即授权预检通过后执行，不另外等待普通确认。
4. 先写EXECUTION_STARTED.json，再执行一次接收侧处理并冻结输出；随后执行评价与必要冷复核。
5. 旧文件/源码/波形字节不得覆盖。每次局限记录PARTIAL和原因；不补信号、不调整参数救结果。
6. Commit B推送、核对remote/main后STOP。回传设计SHA、执行SHA、本任务SHA及最关键量化结果。禁止自动开展H2回灌、漂移提供器修补或硬件增强。

本轮的价值是把“旧处理丢掉了什么”与“现存记录还保留什么可观测变化”分开，避免反复重跑一个无法建立绝对速度零点的处理流程。

## 来源（截至parent的已读仓库资料）

- results/R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT/GPT_SYNC.md
- results/R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT/DECISION.json
- results/R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT/METHOD_INPUT_LOCK_READABLE.md
- hla_rx1_extract.py；hla_rx1_generate.py
- results/R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT/OUTPUT_MANIFEST.json
- results/R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT/EXTRACTED_FEATURES.csv（已确认载频/offset/窗口等字段存在）

本文件提出的跨窗差分与其评价是新设计，不是上述文件已经完成的新方法结果。
