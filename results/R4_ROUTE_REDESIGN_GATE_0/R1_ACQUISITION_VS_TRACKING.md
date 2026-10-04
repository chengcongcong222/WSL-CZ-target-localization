# R1：Acquisition 与 Tracking 的边界

阶段 `R4_ROUTE_REDESIGN_GATE_0_COLD_START_INFORMATION_AUDIT`；接受 parent closeout `206f6c2f46bc3c52890a32004cf54255ed30fcef`。本轮只审计既有项目证据和逻辑必要条件，新增数值实验 0。

结论：`R1_COLD_START_NOT_ESTABLISHED`；`R1_TRACKING_ONLY_CONDITIONALLY_PLAUSIBLE`。当前没有识别出能在继承含噪场景中从宽距离域可靠进入 RC3 local-anchor 区域的合法冷启动缩域来源。这个结论不是所有场景 cold start 不可能，也不否定无噪声 bearing 唯一解。

| 项目 | R1-A：Cold-start acquisition | R1-B：Tracking / continuation |
|---|---|---|
| 输入状态 | 无可信 previous posterior；当前 RC1-conditioned range support 45–60 km | 已有可信且足够局部的 previous posterior/candidate support；来源必须独立、可审计 |
| 允许信息 | 场景/第一 CZ 先验、continuous bearing、已知平台轨迹、合法速度范围、运动学、当前 RC3 observable | 上述信息与合法前时刻 support；prediction → candidate continuation → local rescoring → measurement update |
| 当前缺口 | 未建立既保留正确区域、又将宽域缩至可捕获 local-anchor 域的 bootstrap | narrow support 的来源、局部捕获、模型误差、失锁和恢复策略尚未验证 |
| 当前证据 | RC1 无窄域 admission；直航/小转向粗候选域未收缩距离；FIX1 noisy continuous RC2 导出跨度仍近 15 km；global capture 失败 | 已有运动候选维持与局部 acoustic 区分证据，使条件式 refinement 值得讨论；没有完整 tracking 实验 |
| 工程 claim | 不允许 `COLD_START_PASSIVE_LOCALIZATION` | 只能暂用 `POST_ACQUISITION_TRACK_REFINEMENT` 的条件性提案，不能声称已实现或已通过 |

## 可观察信息与范围约束

RC1 的正式 [support role](../R3_FINAL_PARAMETER_CONTRACTION/RC1_SUPPORT_ROLE.md) 写明：所有排除率都在冻结的 114576-node RC1-conditioned grid 内，range support 为 45:1:60 km，没有测试更宽 pre-RC1 prior，`RC1_CONTRACTION_RATIO = NOT_IDENTIFIABLE_FROM_CURRENT_EXPERIMENT`。因此第一会聚区是条件场景/支持域，不是已经观测到的精确距离。没有一份已接受证据给出更窄、且对未知环境/深度有效的 first-CZ cold-start prior。

先验的物理解释依赖传播模型和环境条件。当前 E0/固定 modal inputs 是匹配模型，不是带观测误差的 SSP 后验；depth 在继承 labels 上作 nuisance profiling，不能预先固定为真值。多个 CZ 的分区/序号误判、环境变化下 band 的偏移以及先验保留率未在本链条验证。不能宣称一定有多 CZ alias，也不能假定第一 CZ 归属已由测量无歧义确定。SSP/depth 仍不开放。

## bearing 正证据与含噪限制必须同时保留

[FIX1 CASE_LEVEL_RESULTS.csv](../R4_A1_FIX_CONTINUOUS_SEARCH/CASE_LEVEL_RESULTS.csv) 和 [MODE_INVENTORY.csv](../R4_A1_FIX_CONTINUOUS_SEARCH/MODE_INVENTORY.csv) 显示：15 个 noiseless 控制均有 rank-four bearing geometry、导出 RC2 range width 0，标记唯一 noiseless solution；21 个 nominal noisy case 均保留 truth-in-likelihood，但导出 RC2 range span 为 14.9008764373–14.9987589195 km，中位 14.9734269641 km。此汇总只读取既有表；导出 cloud 是有限采样，不能当作完整 confidence support 或严格连续域边界。

因此：continuous bearing 在合适几何、无噪声条件下可以含绝对距离信息；原 A1 coarse-cell rejection 不是 bearing 完全没有信息的证明。当前含噪、bounds 与 likelihood 条件下没有建立足够缩域能力，truth retention 也不等于 domain contraction。轴向 range span 不能单独描述耦合 support 的维数或体积，但现有结果没有进一步建立沿薄脊线的可审计 local-anchor capture。

直航 bearing 的尺度补偿可从既有运动学作逻辑说明：固定平台速度 u 时，relative trajectory 为 p0+(v-u)t；同时把 p0 和 (v-u) 乘同一正比例，在变换后的速度仍合法时保留 bearing。速度/航向 bounds 能限制允许比例，却不是独立距离测量。此符号说明不是新仿真，且不适用于把转向平台也判为同样精确对称；小转向能改变几何，现有 noiseless full rank 正是不能忽略的例外。

## 当前可接受的 tracking 条件

可信 previous support 必须由独立、合法、observation-derived acquisition 提供，并保留多峰、r-v-psi 耦合、深度 nuisance 和误差预算。FIX2 recovered state、truth、oracle basin 均不能作为工程输入，也不能把外部 acquisition 静默添加到当前单平台场景。

Prediction 传递既有支持，并不凭空创造距离信息；运动/过程误差还可能扩大支持。需要 measurement update 与实际信息增量才可能进一步收缩。previous posterior 可信且局部是必要条件，不是 tracking 自动成功的充分条件。

[FIX2 局部几何与独立恢复](../R4_A1_FIX2_ACOUSTIC_COVERAGE/R4_A1_FIX2_REPORT.md) 支持条件性局部区分的机制角色；5.37 mm likelihood section 不给出 attraction basin、初始化容差或连续捕获证书。tracking/refinement 当前仅 `CONDITIONALLY_PLAUSIBLE`，实现和数值 Gate 都未放行。

下一步建议负责人转 R2 的 acquisition observable/information-source 文档设计；R1 可留作后续取得合法初始 support 后的条件子模块。R1-B 的合理性不能证明 R1-A 成立。
