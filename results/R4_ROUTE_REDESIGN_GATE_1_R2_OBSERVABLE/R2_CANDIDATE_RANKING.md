# R2 Gate-1 候选排序

当前没有 Class A。Priority 1 为 **CZ_ENVELOPE**，Priority 2 为 **INTERFERENCE_FRINGE**；两者均 `CONDITIONAL` / Class B。排序用于决定先澄清哪条信息链，不是性能评分，不授权 quantitative Gate 设计、实现或实验。

| 优先序 | 候选 | admission / class | 排序理由及限制 |
|---|---|---|---|
| 1 | CZ_ENVELOPE | CONDITIONAL / B | 距当前 relative-level 输入最近，不需 known waveform 或 absolute phase；coarse shape 有希望替代细振荡匹配。但源/增益稳定性、short-window edge 可见性及归一化后的绝对 range relevance 未闭合 |
| 2 | INTERFERENCE_FRINGE | CONDITIONAL / B | intensity-only、单 HLA/单点可讨论，不需 source phase；period/slope/局部 β 可构成 structured ambiguity。三个频点没建立 dense fringe，UUV 有效带宽、未知运动与 β/环境转移是关键缺口 |
| 3 | ORDERED_DIFFERENTIAL_FREQUENCY | CONDITIONAL / B | 少频统计靠近现有输入，可辅助 A/C proposal；跨频源谱偏置未消除，ordering 会丢信息，不能只是 exact-TL 改名 |
| 4 | MULTIPATH_TDOA | CONDITIONAL / B | 独立时延信息的物理上界最强，但当前 tonal/level-only chain 没有可靠 lag；原文 VLA 方法不匹配，旧 aggregate single-depth 支线关闭，缺口比前两项大 |
| 5 | DOPPLER | CONDITIONAL / B | stable-tone frequency 可作为运动辅约束，不需 absolute source phase；未知 f0/drift 和短时 motion excitation 限制，不能单独 anchor absolute range |
| 6 | ARRAY_SPATIAL | REJECT_CURRENT_SCOPE / C | inter-sensor phase 条件可测但主要 bearing；所测小 aperture 与 modal identity 没有新 range 锚。以 MFP 重做 exact global matching 不满足 representation 改造要求 |

## Class A / B / C 的明确边界

Class A=`R2_CANDIDATE_ADMITTED_FOR_QUANTITATIVE_GATE_DESIGN`，须 G1-1–G1-6 同时有结构依据。它不要求本轮证明 numerical accuracy，却至少要有实际 observable/合理信号条件、处理 nuisance 后的 range mechanism，以及可区别于原 exact matching 的 search structure。当前没有谁达到这条线；不能为了产出 top-2 就填一个 Class A。

Class B=`R2_CONDITIONAL_PENDING_MEASUREMENT_CHAIN`。A/C/D 缺 source/gain/representation 的保留证据；B 缺真实 lag 与稳定性；E 缺独立 Doppler 条件。Class B 可以有多个缺口，不只缺硬件 raw signal。结构信息不足也不等于物理不可能。

Class C=`R2_REJECT_CURRENT_SCOPE`。当前 F range-anchor、B 的 Xu2024 VLA 直接移植、C 的 constant-β 直接移植和 D 的未经源谱处理的单次 ordering 都属于此类具体 variant；后三者并不把母候选的条件机制全部关闭。

## 下一阶段建议

推荐 **`R2_MEASUREMENT_CHAIN_AND_SIGNAL_CONDITION_SPECIFICATION`**，仍为负责人另行授权的文档阶段。先说明真实 UUV 有哪些可采集谱/时域条件，接收同步、阵形、gain、line tracking、源谱稳定性有哪些可审计支持；再判断 A 的 envelope 是否保留绝对 range 结构、C 是否有足够频率采样及合法 motion coordinate。不得静默假定 broadband、known f0、known waveform、coherent source、VLA、Z maneuver 或先前已定位。

若这些来源仍不能成立，下一步应判定需要改变 signal/measurement contract 或重新设计 R4 acquisition，而不是选择第三 optimizer。未授权采购、接收系统实现、特征提取或试验。本轮不设 band、SNR、aperture、平滑窗口、proposal width 或 accuracy 数值门限；不开发候选。

完整信息链见 [接收链审计](R2_MEASUREMENT_CHAIN_AUDIT.md)，证据迁移见 [文献审计](R2_LITERATURE_TRANSFERABILITY.md)。R4=0%；A2/depth/SSP/P5 UNOPENED。
