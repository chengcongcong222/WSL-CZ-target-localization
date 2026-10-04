# Route redesign Gate-0：冷启动信息来源审计

parent closeout：`206f6c2f46bc3c52890a32004cf54255ed30fcef`，研究负责人已接受 `R4_A1_CLOSEOUT_AND_ROUTE_DECISION_ACCEPTED`。本轮无新数值实验、仿真、Monte Carlo、optimizer、forward recomputation、重跑 R3/A1 或新 panel。

## 判定与作用域

采用用户定义的**结论 B**：

```text
R1_COLD_START_NOT_ESTABLISHED
R1_TRACKING_ONLY_CONDITIONALLY_PLAUSIBLE
```

在当前继承含噪 nominal bearing、单平台既定机动能力、45–60 km 条件距离域及已有 RC3 observable 下，没有识别出足以支持 local-anchor acquisition 的合法 cold-start 缩域来源。这里 `cold_start_information_source_identified=false` 指没有找到符合该作用域且具备所需缩域/保留证据的来源，不是说所有观测毫无距离信息，也不是物理不可辨识证明。

推荐下一路线：**R2 acquisition observable/information-source design**，R1 仅保留为 `POST_ACQUISITION_TRACK_REFINEMENT` 条件子模块；Route R3 保留 oracle/diagnostic。当前 acquisition 研究优先讨论 R2，再保留 R3 的诊断用途；不把条件 tracking 与 acquisition 放在同一性能排名中。任何下一轮设计/实现/实验都须负责人另行授权。

## S1–S6 来源审计

| Source | Available at cold start? | Truth-free? | Adds absolute-range information? | Existing positive evidence | Existing negative evidence | Expected role |
|---|---|---|---|---|---|---|
| RC1 first-CZ prior | 冻结场景 band 有；窄观测 prior 未建立 | 场景支持域可以；窄 prior 无来源 | 有条件 domain restriction；无更窄 acquisition anchor 证据 | RC1-conditioned 45–60 km support | RC1 contraction ratio 未识别；未测试 wider pre-RC1 prior；环境/多 CZ/depth 不确定性 admission 未闭合 | 广域支持/模型条件 |
| continuous bearing | 是 | 是 | 理想 full-rank geometry 可有；nominal 足够缩域未建立 | FIX1 truth coverage 36/36；noiseless 15/15 singleton | 21 nominal 导出 range span 14.9009–14.9988 km；原 A1 coarse selection blocker 不能转成物理否定 | theta、运动耦合与候选可行域 |
| straight-line temporal accumulation | 当前合法采样时间内是 | 是 | 所测 range contraction 0 | R3-1B 候选 1326→495，psi 宽 13→11° | r 宽仍 15 km；更少候选不等于更窄距离支持 | motion consistency，不能作已验证距离锚 |
| small platform turn | 所测 0–15° 是 | 是 | 几何潜力/局部协同有；noisy RC2-only anchor 未建立 | 部分 turn+RC3 配置有收缩，noiseless rank-four 可解 | R3-1C RC2-only 五转角仍 15 km；主判 ridge persists | 几何激励；不等于 Z maneuver/large-baseline |
| target speed bound | 是，继承 1–3 m/s | 来源合法的物理 bound 是 | 限制可行尺度；非独立 range measurement | 排除不合法速度状态 | 与 RC2 联合时仍有宽 r-v compensation；未做独立 bound ablation | candidate admissibility |
| target heading continuity | 常速/常航向模型可用，heading 未知 | 无 truth heading 时是 | 未建立独立 absolute-range anchor | psi 与 temporal consistency 有约束 | heading/velocity 可共同补偿 range；没有已知航向测量或独立缩域证书 | motion coupling/continuation |
| previous posterior | 首次 acquisition 不可用 | 仅当有合法独立来源 | 可承载既有范围信息；不创造首个 anchor | 能定义条件 tracking 接口 | 来源若回到失败 global matcher 则循环；未有实际 tracking admission | post-acquisition 条件输入 |
| RC3 direct-TL local score | 可评分；可靠局部初始化未提供 | scorer API 是；oracle geometry 不可作初始化 | 强局部 range-sensitive matching，不是已建立 global proposal | FIX2 窄 sections、独立近零恢复 | 冻结 SHGO/DIRECT 各 T1/T2/T3 0/21；无 global capture 证书 | LOCAL_CONDITIONAL_ANCHOR_ONLY |

机器表见 [COLD_START_INFORMATION_SOURCE_MATRIX.csv](COLD_START_INFORMATION_SOURCE_MATRIX.csv)，含 source paths 与各项 admission 状态；以上没有未验证项被写为科学 PASS。

## 必须保留的反例、范围与缺口

1. S1：第一 CZ 不能直接变成精确 range；旧 RC1 文档明确只在条件域内算排除率。固定 E0 不是环境不确定性证书；unknown source depth 被 profile，不能替换为已知。multi-CZ/SSP 不确定性尚未审计，本轮不开放。
2. S2：FIX1 15 个 noiseless singleton 是 bearing 正证据，不是 noisy cold-start 通过。21 noisy cases truth retained 与 sampled range span 近 15 km 同时成立；有限 cloud 不是 confidence region。
3. S3：0–15° 单次小转向不是 Z maneuver 或更大 baseline。MAIN tau=0.5、15° 时三个 depth 的 RC2+RC3 range width 为 4/15/15 km；FOUR reference 同场景有单距离格，是有条件 coarse control，不能隐藏。它没有修复后来 off-grid nominal acquisition failure，也不等于 0.001 dB capture Gate。
4. S4：合法 v∈[1,3]、psi∈[-15,15] 已包含在所测 broad domain 中。bounds 与连续运动学筛选状态，不等于 range measurement；没有独立 ablation 就不把其边际贡献写成严格零。
5. S5：prediction 传播已有信息；过程误差还可扩大 support。第一 posterior 未建立时，用它解释冷启动是 `COLD_START_CIRCULAR_DEPENDENCY`，不是已有代码 bug。
6. S6：5.37 mm 是局部 likelihood section，不是 attraction basin 或必须达到的先验宽度。本轮不设“多窄才可 local”的新数字，不把 local discrimination 等同 local convergence/capture 已通过。

来源：[RC1 role](../R3_FINAL_PARAMETER_CONTRACTION/RC1_SUPPORT_ROLE.md)、[R3-1B](../R3_RC23_CLOSEDLOOP/R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION/R3_RC23_CLOSEDLOOP_1B_REPORT.md)、[R3-1C](../R3_RC23_CLOSEDLOOP/R3_RC23_SINGLE_TURN_RANGE_ANCHOR/R3_RC23_CLOSEDLOOP_1C_REPORT.md)、[Original A1](../R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_DECISION.json)、[FIX1](../R4_A1_FIX_CONTINUOUS_SEARCH/R4_A1_FIX_REPORT.md)、[FIX2](../R4_A1_FIX2_ACOUSTIC_COVERAGE/R4_A1_FIX2_REPORT.md)、[tractability](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_RUN_REPORT.md)。全部保持历史原判。

## G0-1–G0-5 逻辑必要条件

| Gate | 本轮审计 | 含义 |
|---|---|---|
| G0-1：合法、truth-free、冷启动可用的缩域来源 | NOT_MET_IN_CURRENT_NOISY_SCOPE | 有允许观测/先验，但没有识别出足够缩域且保持正确 support 的 admitted source |
| G0-2：不循环依赖失败 global search | NOT_MET_FOR_PREVIOUS_POSTERIOR_BOOTSTRAP | previous support 独立来源未建立；若借 failed matcher 来产生所需 prior，循环成立；外部 acquisition 条件 tracking 不循环 |
| G0-3：不直接复用已有完全不缩 range 的方案 | NOT_MET_FOR_TESTED_STRAIGHT_ACCUMULATION_AND_RC2_SMALL_TURN | 这两项旧表已给出 15 km 不变；RC1 更窄 subband 无证据，不借未验证先验绕过 |
| G0-4：可定义 support retention 与 domain contraction 指标 | SATISFIED_AT_DEFINITION_LEVEL_ONLY | 可定义未来的 support inclusion、range projection/joint components/保留与缩域联合报告；没有数值门限或新评估结果 |
| G0-5：RC3 限定为 conditional/local anchor | SATISFIED_AT_ROLE_DEFINITION_LEVEL_ONLY | 已给出角色限制；局部捕获与 tracking 仍未验证 |

所有逻辑必要条件没有同时满足，R1 cold-start 不进入 quantitative Gate 设计/实现。G0-4/G0-5 的定义可行不是 scientific PASS，不增加 R4 进度。

未来若另行授权指标设计，必须分别检查：obs-derived support 是否保留 evaluation-only reference；domain contraction 是否伴随 retained support；projection width、多个连通分量和 joint coupling 如何记录；有限 samples 如何与 continuous support 的证书区分。不得只看 top-1 或把 truth 交给 estimator；本轮不生成数据、不设新阈值。

## 交付及停止

只新增本目录的七项 artifacts、更新 master plan、追加 ledger。旧 closeout 与全部 R3/A1/FIX/FIX2/tractability 数值、源码、预算、Gate、manifest 保持不变。旧 R4_PROGRESS.json 仍是 historical development snapshot；当前读取 [R1_GATE0_DECISION.json](R1_GATE0_DECISION.json)。本轮 hash/表一致性证据见 [VALIDATION.json](VALIDATION.json)。

R4=0%；A2/depth/SSP/P5 UNOPENED。禁止新 optimizer、仿真/MC、重跑 R3/A1、T4、FIX3/FIX4、fresh panel。推送本轮审计后停止，不启动 R2 observable 开发或 R1 tracking 实现。
