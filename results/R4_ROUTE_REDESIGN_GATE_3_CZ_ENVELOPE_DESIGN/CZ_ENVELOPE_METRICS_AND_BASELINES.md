# 指标、基线与 future table schema

本轮只定义指标，不填新科学数值。baseline cold-start range 域为 45–60 km，宽 15 km。已知历史 [R3-1B 表](../R3_RC23_CLOSEDLOOP/R3_RC23_TWO_WINDOW_TEMPORAL_CONTRACTION/RC2_CUMULATIVE_SET_METRICS.csv) 中 candidate 1326→495 而 width 15→15 km；candidate 数减少不能替代 range contraction。

| 指标 | 未来定义 | 必须同时报告的限制 |
|---|---|---|
| M1 SUPPORT_RETENTION | truth/reference horizontal state ∈ exported joint support；另报 range projection retention | truth 只给 evaluator；empty/missing/invalid 不删；noiseless bearing singleton 仅 regression |
| M2 RANGE_PROJECTION_WIDTH | W_hull=max(r)-min(r)；另报 interval union length W_union | main Gate 用 hull 防两端小岛；空集不能算零宽成功 |
| M3 CONTRACTION_RATIO | C_r=W_hull/15 km；另报与同 case B0 RC2 width 的条件增量 | 如 B0 已先缩域，不能把其贡献全记 CZ；15 km为固定prior尺度 |
| M4 CONNECTED_COMPONENT_COUNT | closed range-interval union数及 face-adjacent 4D outer-cell component数 | resolution作用域、unresolved连接、inner witness islands另报；不事后merge/delete |
| M5 JOINT_MOTION_COUPLING | retained x 的 r-v-psi projections、每range段内motion support、成分间相关/补偿 | theta保留；z仅profile nuisance；不fix truth motion/depth作主结果 |
| M6 LANDSCAPE_COMPLEXITY | conservative broad-region数、全域bounds/closure、两级grid endpoints/topology稳定性；可用score roughness诊断 | low-dimensional F不保证参数域平滑；hard-cap未闭合不能当coverage收敛 |

## roughness 与配对机制诊断

若下一轮明确授权，可用同一 observation-only 预注册候选 lattice、固定 v/theta/psi 枚举及z-profile，对 s_F 的归一化 total variation/second difference 作纯诊断；尺度标准化方式必须事前写入，零 dynamic range 标记 flat/uninformative，不能因 roughness 低就获 PASS。不在 truth motion 上选择 slice。range-axis诊断不能代替4D closure。

若需要 paired exact pointwise reference，只复用同一候选幅度预测、同一grid和调用计数；不新增 optimizer、recovery polish 或尝试复活 B1。它是额外需明确授权的 mechanism diagnostic，不是本轮执行或必需重算旧失败。不能与历史异域/异预算 J 或 search sample cloud 作定量平滑度排名。

## 三类基线

B0：RC1 prior + 同 case observation-only RC2 support，不加 envelope。报告prior15km与actual B0宽度/retention/成分，不假定每case都恰15km。

B1：旧 exact direct-TL 路线的历史负对照 [tractability decision](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_FINAL_DECISION.json)，SHGO/DIRECT T3均0/21，execution-valid；保留 NOT_ESTABLISHED，不重跑/调optimizer，不推物理不可辨识。

Proposed：选定单一 CZ envelope F；与 B0 比 support retention +实际range缩域 +保守覆盖/复杂度。N0 primary 与 N1 robustness 分开；静态 exact matched existence/control 不给estimator seed。不得把noise-free bearing singleton数当acoustic acquisition证据。

## 输出表结构（定义而非空实验结果）

`CASE_SUPPORT_METRICS.csv`：case_id、split_role、representation_id、nuisance_id、source_provenance、execution_valid、coverage_status、truth_joint_retained_eval_only、truth_range_retained_eval_only、truth_cell_resolved_eval_only、B0_width_hull_km、width_hull_km、width_union_km、C_r_prior、C_r_given_B0、n_range_components、n_joint_components、unresolved_cell_count、motion_support_summary、profile_scope、gate_tier、failure_reason。

`RETAINED_RANGE_INTERVALS.csv`：case_id、representation/nuisance_id、component_id、r_lo/r_hi、support_status、grid_level、boundary_censoring；所有 intervals 输出，禁止仅主峰。

`RETAINED_JOINT_CELLS.csv`：cell_id、4D lows/highs、bearing/feature_bounds、depth/drift scope、retained/rejected/unresolved、bound provenance、grid level。estimator产物不含truth列。

`SEARCH_AUDIT.csv`：objective requests、exact forward state/calls、cache hits、blocked calls、hard-cap limit、closure、raw fine/coarse一致性、effective z mappings、all budget overhead。wrapper严格计数，不用optimizer maxev选项当hard cap。

`LANDSCAPE_DIAGNOSTICS.csv`：same-grid region/component、normalized roughness、flat-feature标志、endpoint/topology drift、reference role/provenance。未做reference时标NOT_RUN，不填0冒充结果。

上述schema必须在执行前和truth/evaluation分离规则一起冻结；本轮未创建这些CSV结果或运行任何计算。
