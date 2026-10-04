# Route redesign Gate-1：R2 acquisition observable 审计

parent Gate-0：`193abb7cd13a38fac91ad1715e6873dd3d6440e9`；研究负责人接受 `R4_ROUTE_REDESIGN_GATE0_ACCEPTED`。本轮六类 observable/information-chain 设计审计已经完成，但**没有候选获准 quantitative Gate design**，不标记 scientific PASS。

Priority 1=`CZ_ENVELOPE`，Priority 2=`INTERFERENCE_FRINGE`，均为条件候选；其选择依据是与现有被动幅度链接近、不需绝对相位，并有 coarse/structured representation 的机制期待。不存在已验证的 km/百米 acquisition 锚，也没有估计 feature landscape 或新增 threshold。

## 审计范围与结果

| 候选 | observable reality | Admission | 分类 |
|---|---|---|---|
| CZ envelope | CONDITIONALLY_MEASURABLE | CONDITIONAL | B |
| multipath TDOA | NOT_CURRENTLY_AVAILABLE（所需 lag，不是声学原理不存在） | CONDITIONAL | B；Xu VLA 直接移植拒绝 |
| interference fringe | NOT_CURRENTLY_AVAILABLE（dense striation；少量 level 可用） | CONDITIONAL | B；constant-β 直接移植拒绝 |
| ordered/differential frequency | CONDITIONALLY_MEASURABLE | CONDITIONAL | B；单次未知源谱 ordering 不放行 |
| Doppler | CONDITIONALLY_MEASURABLE（稳定 tone 条件） | CONDITIONAL | B；仅运动辅约束 |
| array spatial | CONDITIONALLY_MEASURABLE（同步/校准原始通道条件） | REJECT_CURRENT_SCOPE | C；当前额外 range-anchor 方案 |

`MEASURABLE` 不能由 forward model 内部 pressure/phase 推出。bearing、三个频点的 demeaned relative level 是现有 API；真实 UUV 的 signal condition 与 raw receive contract 尚未建立。共源同频的阵元相对相位、跨频源相位、时间相位连续性、绝对传播相位逐项区分，不偷模型相位。

TDOA 的 labeled MMAC 与 unlabeled aggregate lag 不是同一问题；后者仍不能绕过旧 single-depth extractability 关闭。Xu2024/海底反射区全文显示 broadband 与 vertical depth sampling 条件，不能直接移植成 single-HLA UUV。有效 Xu2024 PDF 已在 D 盘原始项目资料重新找到，SHA 与后续 primary lock 相同；库内非 PDF 同名文件保留不动。

fringe 有 amplitude-only 测量机制，不需绝对 phase。三个非均匀频点没有证成可提取的 fringe period，constant-β 在 E-STD 的迁移限制仍在；不宣布整类无效，也不宣称现已可用。CZ envelope smoothing 则需证明在 short observation、source drift、未知 motion/depth/environment 下保留 coarse absolute range，而不只是第一 CZ band。D 差分不会自动消未知源谱；E frequency shift 主要测 range-rate；F small-aperture 相位主要测 bearing，depth/modal mechanisms 不补首个距离锚。

## G1-1–G1-6 结构条件

| Gate | 放行所需依据 | 当前结论 |
|---|---|---|
| G1-1 Observable reality | cold-start 可测统计或明确合理 signal condition | 实际接口/源信号条件未闭合；模型 upper bound 不能替代 |
| G1-2 Truth independence | 无 truth、oracle path/phase/mode、已知 source position | 能定义 truth-free 条件方案；不得借 privileged 输入填其他 Gate |
| G1-3 Absolute-range relevance | absolute range 或破除 r-v-psi 的独立信息 | delay/fringe/envelope 有潜在机制，但当前可用输入+nuisance removal 未建立 acquisition relevance |
| G1-4 Nuisance robustness | source level/phase/spectrum、depth、environment 均处理 | constant per-line 去均值有既有支持；时变源谱/增益及转移条件缺口保留 |
| G1-5 Searchability | 有 coarser/smoother/lower-dimensional/ambiguity-structured 理由 | A/B/C/D/E 有结构期待，未验证；F 若仍 exact MFP 就不满足；平滑可能丢信息 |
| G1-6 Compatibility | 单平台/HLA、first CZ、UUV、短观测、小转向、无 Z | 原文 broadband VLA 方案不兼容；generalized A/C/B/E 各需信号条件，不能静默换场景 |

这些是 logical design 条件，不是数值 Gate；没有新增 band/SNR/aperture/window/proposal-width 门限，也没有给出 retention/contraction 或 coverage 的新结果。各候选的逐 Gate 缺口见 [decision](R2_GATE1_DECISION.json)。Class A=0，Class B=5，Class C=1（按母候选，具体 rejected variants 另列）。

## 交付与停止

推荐下一阶段仅为 `R2_MEASUREMENT_CHAIN_AND_SIGNAL_CONDITION_SPECIFICATION` 文档，先冻结真实可观测量与 source/gain/band/motion 条件，再由负责人决定是否允许 quantitative Gate design。本次排序不是实验启动许可。R1=`NOT_ESTABLISHED` / `CONDITIONALLY_PLAUSIBLE_POST_ACQUISITION_ONLY`；RC3=`LOCAL_CONDITIONAL_ANCHOR_ONLY`。

交付：[admissibility matrix](R2_OBSERVABLE_ADMISSIBILITY_MATRIX.csv)、[measurement chain](R2_MEASUREMENT_CHAIN_AUDIT.md)、[literature transfer](R2_LITERATURE_TRANSFERABILITY.md)、[ranking](R2_CANDIDATE_RANKING.md)、[decision](R2_GATE1_DECISION.json)、[GPT sync](GPT_SYNC.md)、[validation](VALIDATION.json)。只追加 master 当前状态和 ledger；不改历史数据、代码、预算、Gate、manifest、Gate-0 或 closeout。

新数值实验 0；R4=0%；A2/depth/SSP/P5 UNOPENED。不得新仿真、BELLHOP/KRAKEN、MC、optimizer/panel、TDOA/fringe estimator、T4/FIX3/FIX4，不重开 B1/range-UB。提交推送后停止。
