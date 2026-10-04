# GPT 同步：Gate-2 接收与信号合同

研究负责人接受 Gate-1 `9e50b0d949523c85fb35c1168ff77c37ddf37aeb`；本轮 Outcome B：`R2_REQUIRES_MEASUREMENT_CONTRACT_EXTENSION`。CZ_ENVELOPE 与 INTERFERENCE_FRINGE 都是 PROPOSED_CONTRACT，尚无获准 quantitative design 的候选。

M0 只是已用合成 bearing/三频 relative level；M1 的实际 waveform/PSD、sync/calibration、line detector、Doppler/dense spectrum 没有产品证据。原始项目方案的技术路线和资料可用假设不能当设备能力；第六节和既有问题清单仍要求确认阵列、目标频谱/录音与观测时长。本次绑定这些资料的身份并保持原文不动。

Contract A 选择多条可关联谱线的时间強度变化，不强制 dense PSD，要求可审计 gain、源漂移、归一化及质量记录。稳定 per-line 源偏置可同频时间消去，任意时变源谱不能；要在未知运动/深度/环境下保留 CZ coarse 区段差别，否则 FAIL。M0 可以表达一部分 temporal proxy，但不能替代实际测量证据。

Contract B 提出 target-bearing dense intensity spectrogram，要求足够真实谱占用、time/frequency axes 和 source/gain trend separation；不能插值三频冒充 dense fringe，不能用 truth radial motion 造 range coordinate。slope/spacing 需联合 range、motion、环境/β 与 family/order；没有唯一距离声明。

两份草案均不需 known waveform、absolute source phase、known source spectrum、known f0 或发射端配合。上游 beam 的相对阵元同步/标定不能被省略，也不等于 source phase 已知。若 producer 能给可复核 spectral 产品，R2 不强制读取所有 raw channels。没有授权设备能力扩大。

G2-1 未闭合：负责人尚未接受拟议输出扩展。其他 G2 项是条件结构推理/可证伪设计，不是实测或科学 PASS。建议下一阶段仅 `R2_MEASUREMENT_CONTRACT_EXTENSION_DECISION`，由负责人选择接受扩展或补证。当前不启动 quantitative design、实验、采集或实现。

[TDOA、差频、Doppler、spatial 最低合同](R2_MEASUREMENT_CONTRACT_MATRIX.csv) 保留全部既有边界：MMAC/single-depth 不重开，Doppler/spatial 为辅助，constant-β 不直接迁移。参见 [报告](ROUTE_REDESIGN_GATE2_REPORT.md)、[接收审计](CURRENT_RECEIVER_CHAIN_AUDIT.md)、[信号条件](R2_SIGNAL_CONDITION_MATRIX.csv)、[决策](R2_GATE2_DECISION.json)、[完整性](VALIDATION.json)。

本轮新数值实验 0，R4=0%，A2/depth/SSP/P5 UNOPENED；只追加 master plan/ledger。提交推送后停止。
