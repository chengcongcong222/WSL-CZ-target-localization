# R4 路线决策

排序：**R1 > R2 > R3**。推荐 `RC2_TEMPORAL_CANDIDATE_NARROWING_PLUS_RC3_LOCAL_ANCHOR`；状态 `R4_ROUTE_REDESIGN_REQUIRED`。这是待负责人决策的架构提案，没有实现、实验结果或进度信用。

## 三条路线

| 优先级 | 路线与结构 | 现有证据为何支持讨论 | 未知与停止条件 |
|---|---|---|---|
| 1 | R1：continuous bearing、已知平台运动、temporal recursion、state prediction、previous-time posterior/candidate continuation 维持多候选；RC3 在缩小域内局部锚定、排歧、校正 | FIX1 修复 continuous RC2；FIX2 有强局部匹配结构；本轮未建立宽域 blind capture。符合 candidate maintenance ↔ propagation anchoring 的分工 | RC2 缩窄绝对距离域未建立；不得假设已有可靠 previous posterior。若只能删除正确候选、借 truth 起步或持续保留宽距离脊线，R1 前置 Gate 不闭合 |
| 2 | R2：重设计低维 phase/fringe-derived、frequency-dependent、ordered、envelope/differential 或 coarse range-anchor 特征 | exact relative-TL 结构过尖，当前通用搜索未捕获阈值区域，应重审 observable 与 proposal 结构 | 平滑不代表保留判别信息；可能引入别名及 nuisance 耦合。不能假定 relative-TL 幅度数据包含可恢复的相位；新观测须明确数据/传感器来源并另行授权 |
| 3 | R3：现有 direct-TL matcher 仅为 oracle mechanism validation、local diagnostic、idealized upper-bound control、conditional posterior rescoring；退出 global proposal generator 角色 | 保留局部正证据与 tractability 负证据，放弃当前 blind global recovery 工程承诺 | 不单独构成工程 estimator；oracle 是理想机制参照，不是 CRLB、已证明的性能上界或无条件 posterior certificate |

R3 可作为未来 R1/R2 的诊断用途；排第三指它不能独立解决工程定位。R1 排第一基于架构与现有证据，不代表它已通过或实测优于 R2。

## R1 分工及前置约束

```text
RC1（若使用：明确可观测先验来源与不确定性）
  → RC2 candidate maintenance / temporal prediction
  ↔ RC3 local propagation anchoring
  → temporal update / 多候选继续维持
```

RC2 维护 `(r,v,psi,theta...)` 的观测支持域，保留多峰、运动耦合与遗漏风险；RC3 提供条件传播锚定。必须同时检查“域缩小”和“正确区域未被删除”。top-1 方差变小不等于 support coverage，不得制造毫米先验来绕过 capture 问题。

必须带入历史负证据：

- [R3-1B](../R3_RC23_CLOSEDLOOP/R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION/R3_RC23_CLOSEDLOOP_1B_REPORT.md)：直航 bearing 从 600 s 累积至 1200 s，候选 1326→495，但 range width 仍 15 km。时间推进本身未建立绝对距离缩域。
- [R3-1C](../R3_RC23_CLOSEDLOOP/R3_RC23_SINGLE_TURN_RANGE_ANCHOR/R3_RC23_CLOSEDLOOP_1C_REPORT.md)：所测 0–15° 转向，RC2-only range width 仍 15 km；部分 RC2+RC3 配置有协同收缩，但主判 `SMALL_TURN_RANGE_RIDGE_PERSISTS`，不能升级为 general off-grid 工程先验。
- 这些结论限定于各自冻结场景、粗候选和门限；不重开 R3，也不外推 continuous RC2 普遍不可能。它们要求 R1 给出新的、可审计的信息来源与缩域证据。

冷启动必须明确：first-time support 由哪些允许观测、平台轨迹、合法先验产生？previous posterior 从哪里来？若上一时刻仍依赖当前失败的 blind global search，递推不能自动修复初始化。若讨论 RC1 会聚区先验，必须保留宽度、多区歧义及模型误差，不能把它变成隐藏真值距离。

## 待负责人定义的 route redesign Gate

下列仅是下一份设计的决策清单。本轮不设新实验阈值、不实现、不执行。未来数值准则、panel、预算、起步条件、失败处理与独立确认规则须在执行前单独冻结。

| 待决策项 | 设计必须回答 | 当前状态 |
|---|---|---|
| 信息来源与冷启动 | 哪些观测/平台运动/RC1 或运动先验真正新增距离信息？prior/posterior 的误差、来源、可用时刻是什么？ | NOT_EVALUATED |
| RC2 支持与缩域 | 如何同时检查 coverage 和域/成本收缩，保留多峰及 r/v/psi 耦合，防止只留错误 top-1？ | NOT_EVALUATED |
| RC3 条件捕获 | 如何建立 observation-only capture/convergence，而非 oracle section width 或 truth initialization？沿用 exact matching 不得偷改 0.001 dB Gate | NOT_EVALUATED |
| 成本与拒判 | proposal、预测、scoring、local refinement、validation 如何计量？如何显式处理初始化失败与假锁定？ | NOT_EVALUATED |
| 独立确认 | 所有已看 panels 都是历史/development；何时冻结新确认设计，如何完整保留失败？ | NOT_EVALUATED；本轮不生成 panel |

若 R1 无法提出 observation-only 缩域机制，负责人应转向 R2 或重设计 R4，不能默认放行 T4、FIX3/FIX4 或新 optimizer。R2 改 observable 时是一条新路线；原 A1 Gate/阈值保持原判，不能事后改阈值回填 PASS。

## 本轮效力

`A1_SEARCH_REPAIR_CHAIN_CLOSED`；`R4_ROUTE_REDESIGN_REQUIRED`；R4=0%；A2/depth/B/SSP/P5 UNOPENED。下一轮须是新的 route redesign Gate 设计任务；设计接受不自动授权数值执行。提交本轮后停止。
