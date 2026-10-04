# Route redesign Gate-2：接收与目标信号合同

parent Gate-1：`9e50b0d949523c85fb35c1168ff77c37ddf37aeb`，研究负责人接受 `R4_ROUTE_REDESIGN_GATE1_ACCEPTED`。本轮结论 **`R2_REQUIRES_MEASUREMENT_CONTRACT_EXTENSION`（Outcome B）**。两个草案已具体定义且符合被动单阵结构，但没有证明当前 receiver 已提供所需产品，也没有负责人接受输出扩展；`candidate_admitted_for_quantitative_design=false`。

这不是宣称所有 R2 都需要新增传感器，也不否定 M0 的相对 level 可研究 temporal representations。缺的是实际接收与信号合同。CZ envelope 不强制 dense bandwidth，fringe 才需要足够谱占用与 dense/adequate sampling；二者不强制下游保存阵元 waveform。上游同步/相对校准与 source absolute phase 不混同。

## 三层接收证据

M0 `CURRENTLY_SUPPORTED` 仅为已用模型输入：bearing + 201/235/283 Hz relative levels + 每频每窗口去均值。原 NPZ/API 没有 raw waveform、PSD、源谱或独立 Doppler track；10 s 的特征时间轴不是音频采样率。

M1 真实 multi-channel/beam waveform、dense PSD、line tracking、gain/phase calibration、sampling/sync 能力都未冻结。原始 V0.2 项目稿及问题清单仍要求确认阵列、噪声/频谱/录音与时间尺度；技术路线文字不是设备/数据验收。M2 known waveform/source phase/source spectrum/carrier、active ping、transmitter sync 均不为本轮被动合同所要求或准入。

## 两个最小草案与其他四项

| 候选 | 最低合同 | 当前状态与范围 |
|---|---|---|
| CZ_ENVELOPE | 目标 beam 上多条可关联谱线的时间强度/relative-level 产品，frequency/time、归一化、gain/噪声/失锁元数据 | PROPOSED_CONTRACT；M0 有 temporal proxy，实际测量/源漂移条件未冻结 |
| INTERFERENCE_FRINGE | 目标 beam 的 dense intensity spectrogram，adequate source occupancy、known axes、source/gain trend 与质量记录 | PROPOSED_CONTRACT；三频不能补成 dense observation |
| MULTIPATH_TDOA | broadband/transient/autocorrelation-rich 接收内容、时间分辨、稳定可提取 aggregate lag | CONDITIONAL_SIGNAL_DEPENDENT；不必 known waveform；不必长期 path labeling；当前源信号证据未知，旧 single-depth/MMAC 关闭不变 |
| ORDERED_DIFFERENTIAL_FREQUENCY | 多条稳定 component、可关联相对级与 source-spectrum nuisance model | CONDITIONAL_SIGNAL_DEPENDENT；M0 统计可表达，跨频差不自动 source-invariant |
| DOPPLER | received frequency tracks、时间参考/共享 fractional change 或受约束 source drift，未知 f0 留 nuisance | CONDITIONAL_SIGNAL_DEPENDENT；多条线共享变化不自动给绝对速度/距离，当前独立 stream 未冻结 |
| ARRAY_SPATIAL | 如有 synchronized/calibrated raw complex HLA，可讨论 bearing/spatial coherence 辅助 | PROPOSED_CONTRACT（辅助）；当前额外 range-anchor variant 仍 REJECT_CURRENT_SCOPE，不扩孔径/阵型 |

known waveform 的 matched-filter 变体、known phase、已知源谱与 known carrier 的捷径另列为 NOT_ADMISSIBLE_PASSIVE_SCOPE；不能成为候选暗中依赖。

## nuisance 和 range mechanism

同频时间去均值能消去稳定、频率相关的常数源偏置；不能消任意时间变化的源谱或恢复跨频绝对谱级。单次跨频差只消共同标量，保留未知源谱差。smooth trend、stable ratios、源/传播时尺度分离均为 REQUIRED_SIGNAL_CONDITION，不是事实。

CZ 在条件传播环境中有具物理距离尺度的 edge/broad-maximum 区段；要从 observation-domain 时间形状、未知运动与 nuisance 共同保留多个 coarse range 区间。若只是“第一 CZ 内”、源漂移可任意解释、短窗平坦或环境/深度可补偿，则不能形成新锚。Fringe 的 spacing/curvature 可条件映射到绝对距离尺度，但孤立 slope 只给 range-rate/range/β 的组合；保留 order/family ambiguity，不能用 truth R(t)。

## G2-1–G2-6 结构结果

| Gate | CZ / fringe | 解释 |
|---|---|---|
| G2-1 Raw-data existence | 未闭合 / 未闭合 | 产品需求已定义，现状 UNKNOWN/NOT_FROZEN；提出合同的授权不等于接受该输出扩展 |
| G2-2 Passive scope | 结构满足 / 结构满足 | 不要求 waveform/source phase/active cooperation |
| G2-3 Signal plausibility | 条件合理 / 条件合理 | 可讨论非合作辐射内容；真实 UUV 稳定线谱/dense occupancy 没有证明 |
| G2-4 Nuisance treatment | 模型角色明确 / 模型角色明确 | source level/spectrum、gain、depth、environment 联合或边缘化；不能全部无限自由 |
| G2-5 Range mechanism | 条件机制明确 / 条件机制明确 | envelope 区段尺度 / fringe spacing-curvature；不是已获 absolute range |
| G2-6 Falsifiability | 可定义 / 可定义 | nuisance 后无 coarse range 信息、仅 RC1、源谱伪结构或只剩比例组合时，可 FAIL |

True 仅代表 structural design reasoning；没有当前系统 capability PASS，没有 quantitative validation，不能用表格 true 数量给 R4 加分。既然 G2-1 未接受，不选 Outcome A；两个 truth-free 被动合同可以具体定义，也不按“不现实发射配合”选择 Outcome C。

建议下一阶段 **`R2_MEASUREMENT_CONTRACT_EXTENSION_DECISION`**：由研究负责人决定是否接受所列 receiver 输出及条件，或要求补具体系统证据。此为建议，不由 Codex 自行采购、接入系统、采信号或开始 quantitative design。若扩展获接受，再另行明确 quantitative design 阶段；仍不自动实验。

交付：[receiver audit](CURRENT_RECEIVER_CHAIN_AUDIT.md)、[measurement matrix](R2_MEASUREMENT_CONTRACT_MATRIX.csv)、[signal matrix](R2_SIGNAL_CONDITION_MATRIX.csv)、[Contract A](CZ_ENVELOPE_MEASUREMENT_CONTRACT.md)、[Contract B](INTERFERENCE_FRINGE_MEASUREMENT_CONTRACT.md)、[decision](R2_GATE2_DECISION.json)、[GPT sync](GPT_SYNC.md)、[validation](VALIDATION.json)。Gate-1/全部历史文件保持不变；master 只追加。

新数值实验 0；无 signal generation/STFT/PSD/feature extraction/detector/estimator/传播运行、MC、optimizer、panel 或新性能阈值。R4=0%；A2/depth/SSP/P5 UNOPENED。提交、推送、核对 HEAD==remote/main 后停止。
