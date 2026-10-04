# Gate-3 CZ envelope quantitative design draft

负责人已接受Gate-2 `dcad14337491ac5bf90acacf0be34a57eafac3f1` 并将Contract A准入 **`ACCEPTED_FOR_QUANTITATIVE_GATE_DESIGN_ONLY`**。fringe仍 `CONDITIONAL_NOT_ADMITTED`；真实receiver M1 `UNKNOWN / NOT_FROZEN`。本轮设计已完成，**没有执行、没有新科学数值，quantitative_run_authorized=false**。

单因素问题：在不读truth、保留joint motion/nuisance的前提下，coarse temporal representation能否对45–60km support提供实际range缩域，同时不丢正确state？旧direct-TL仅历史负对照或将来明确授权的机制reference，不复活optimizer路线。

| 候选 | 粗化定义 | M0 | 主要风险 |
|---|---|---|---|
| E1 WINDOW_SHAPE_Q2（推荐） | 每line/旧window两个orthonormal trend/curvature coefficients | compatible | source linear drift吸收trend，quadratic自由漂移可吸收整个feature；真实range机制可能丢失 |
| E2 LINE_BLOCK_MEANS_B3 | 每line/旧window三块coarse means及样本权重 | compatible | block尺度缺物理验证，残存振荡或抹去事件 |
| E3 NORMALIZED_INTENSITY_BLOCKS_B3 | 每line normalized intensity，再equal-line average和三块coarse means | compatible | 非线性尖峰、frequency信息丢失、drift不能在feature后简化 |

推荐E1基于measurement合同、N0 offset invariance、低维coarse机制与显式nuisance投影，不基于运行结果。建议q=2仅draft；多尺度选择尚未获得物理充分性，尤其block尺度 `SCALE_SELECTION_NOT_YET_JUSTIFIED`。负责人下一轮选择唯一primary及最多一个sensitivity，不能跑三项挑最好。

N0 static offsets为basic matched control；N1 bounded linear drift的D尚未冻结，未知gain合并nuisance；N2 arbitrary time variation为`NON_IDENTIFIABLE_NUISANCE_LIMIT`纯符号推论。z为有限旧profile nuisance，禁止truth-motion/depth主曲线，SSP/depth development不开放。

Support输出为保守joint cells与range-interval union，M1 joint retention、M2 hull+union width、M3width/15、M4joint/rangecomponents、M5motioncoupling、M6coverage/roughness/grid stability。empty support不得计零宽成功；unknown cell不能被删除；point sampling或局部convergence不是global coverage证书。cell bounds的可实现性、resolution和hard预算仍是pre-run待选项，不因为F低维就宣布已解搜索。

三档阈值均PROPOSED_NOT_FROZEN：lenient W_hull≤7.5km/C≤0.5；nominal≤5km/C≤1/3（推荐）；strict≤3km/C≤0.2。所有主case joint retention100%，另有component/endpoint稳定性建议。它们是coarse角色建议，不是716<10%验收或RC3capture保证；tau_F必须另由noise/numeric metadata确定，不用Jmin或truth调。

Primary development为旧21个noisy P/Q清单；15个noiseless只regression；旧FIX2的24个“fresh”已全部视为secondary development/regression。未来fresh必须在representation/nuisance/threshold/budget/search全部冻结并得到明确授权后生成；当前NOT_GENERATED。

下一阶段建议 `RESEARCH_LEAD_GATE3_DESIGN_REVIEW_AND_PRE_RUN_FREEZE`。负责人还需选择primary/sensitivity、threshold档、tau_F/noise、N1范围、full-cell enclosure与terminal grid/hard预算、development/confirmation完整规则。本轮准入design不等于这些项已冻结，也不自动授权implementation或实验。

交付：[representations](CZ_ENVELOPE_REPRESENTATION_CANDIDATES.md)、[nuisance](CZ_ENVELOPE_NUISANCE_MODEL.md)、[support design](CZ_ENVELOPE_SUPPORT_ESTIMATOR_DESIGN.md)、[metrics](CZ_ENVELOPE_METRICS_AND_BASELINES.md)、[thresholds](CZ_ENVELOPE_THRESHOLD_PROPOSALS.md)、[split](CZ_ENVELOPE_DEVELOPMENT_CONFIRMATION_SPLIT.md)、[stop](CZ_ENVELOPE_STOP_RULES.md)、[decision](R2_GATE3_DESIGN_DECISION.json)、[sync](GPT_SYNC.md)、[validation](VALIDATION.json)。仅append master plan/ledger，所有历史/Gate-2不变。

新数值实验0；无feature calculation、forward/direct-modal、signal generation、KRAKEN/BELLHOP、MC、optimizer、fresh truth panel或estimator代码。R4=0%；A2/depth/SSP/P5 UNOPENED。commit/push并核对HEAD==remote/main后停止。
