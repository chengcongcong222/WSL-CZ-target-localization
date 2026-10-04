# Development / confirmation 与数据来源

本轮 **fresh_panel_generated=false**，没有创建新truth、source-drift或measurement panel，也没有读取case feature数值做粗化。只核查旧配置/清单/文件身份。

## 已见材料的永久角色

| 旧材料 | 后续角色 | 禁止解释 |
|---|---|---|
| [P development 与 Q旧留出统一清单](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_CASE_MANIFEST.csv) | 21个nominal noisy cases作为primary development/mechanism draft | Q不再fresh，tractability同批也不算独立确认 |
| [CASE_OBSERVATIONS](../R4_A1_FIX_CONTINUOUS_SEARCH/CASE_OBSERVATIONS.npz) 与 [case table](../R4_A1_FIX_CONTINUOUS_SEARCH/CASE_LEVEL_RESULTS.csv) | 同一P/Q的36个archive cases；其中15个noiseless仅regression/control，21 noisy为上行主清单 | 无噪声RC2 singleton不能证明envelope提供距离 |
| [FIX2 old fresh holdout](../R4_A1_FIX2_ACOUSTIC_COVERAGE/FRESH_HOLDOUT_PANEL.csv)、[observations](../R4_A1_FIX2_ACOUSTIC_COVERAGE/FRESH_HOLDOUT_OBSERVATIONS.npz)、[results](../R4_A1_FIX2_ACOUSTIC_COVERAGE/FRESH_HOLDOUT_RESULTS.csv) | 已见24个case（8 controls+16 noisy）为secondary development/regression；若纳入执行，须全量预注册，不选成功子集 | 文件名FRESH不再代表对新representation的新留出 |
| OFFGRID/P/FIX/FIX2各panel及历史失败/independent solver | 整体属于已接触development背景，不得改写其历史失败 | 不用恢复state、mode list或truth附近seed生成proposal |

primary建议21个完整nominal清单；secondary旧FIX2全24的执行选择由pre-run预算冻结决定。两组分别报告，不跨组反复选feature调至通过。已使用数据可帮助未来development，但任何数据驱动修改都必须作为新版本完整披露，不能再用当前draft声称无事后调参。

## 三条数据通道

Estimator：只读取observation-only whitelist（bearing、relative levels/频率时间、平台/合同/quality元数据）和固定模型先验。manifest的case_id仅用于映射与日志，不能通过P09等ID选择算法分支；不读truth columns、old recovered starts或评价排序。

Generator：如果获明确实验授权，用冻结truth产生`SYNTHETIC_MATCHED_MEASUREMENT_CONTRACT_PROXY`；若未来复用现有档案，记录其旧生成SHA、forward/pressure-floor/depth映射与无raw-signal事实。生成truth可用，estimator不可用；actual receive waveform/PSD若不存在不能填入source标签。

Evaluator：estimator输出封存以后才读truth，计算joint/projected retention与width/components；retention检查不得反馈score cutoff、grid、bounds或nuisance。old true effective z只用于评价相容性，不提供预测深度。

## 未来确认的顺序（不执行）

1. 负责人选择representation family/唯一primary尺度、唯一sensitivity（或无）、nuisance及signal/noise合同。
2. 冻结tau_F、range阈值档、所有support/bound/grid规则、hard预算计数、case完整清单、所有输入SHA、停止规则及pre-run工程检查；当前这些仍是草案/未闭合项。
3. 取得明确的development执行授权；本轮没有。development失败按冻结规则停止，禁止先创建fresh再找补。
4. 如果development所选全部Gate闭合，且另获fresh确认授权，才创建新off-grid panel。representation/nuisance/threshold/budget/search rule和panel生成规则、样本数、seed所有权在生成前冻结。样本数和分层设计由预算/科学coverage目标决定，当前`NOT_FROZEN`，不能看到新feature后定。
5. fresh truth独立于旧P/Q/FIX/FIX2，全量保留边界/内域和失败；不得用truth挑“好看”的CZ时段/可匹配shape。确认失败不回头改同一panel。

本轮最多起草未来matched机制结论范围：如果以后通过，说明选定表示在匹配传播、有限z-profile和冻结source/nuisance合同下有coarse acquisition信息，不是已验证真实UUV/receiver、连续深度、海洋环境鲁棒性或工程系统。
