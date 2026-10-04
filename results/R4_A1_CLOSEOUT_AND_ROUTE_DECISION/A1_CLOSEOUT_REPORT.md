# A1 最终科学收口

接受 `A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED`；正式关闭 `A1_SEARCH_REPAIR_CHAIN_CLOSED`。R4 路线为 `R4_ROUTE_REDESIGN_REQUIRED`，R4=0%。本轮新增数值实验为 0。

第三次独立审计接受基线：`9cfa99c14bb80d08e6f16210ebf89fd9c541b850`。审计接受的是结果和证据完整性，不是把失败的科学 Gate 改成 PASS。来源见 [证据链](A1_EVIDENCE_CHAIN.md)。

## execution-valid 的负结果

冻结的 21 个含噪 synthetic matched E0 controls 只向 estimator 提供允许的观测和 manifest；relative-TL 为理想匹配模型观测，bearing 含噪。本结论不具备实海或场景总体概率作用域。

每个 solver/case/budget 的 raw 搜索独立执行；每 case 有 21 个继承 nuisance depth labels，未开放连续 depth estimation。T1/T2/T3 admitted-request caps 为每 case/depth 的 16/64/256。没有累计继承低预算成功、truth start、oracle 邻域或失败后重跑。

| 冻结证据 | 结果 |
|---|---:|
| SHGO T1/T2/T3 threshold-region recovery | 0/21、0/21、0/21 |
| DIRECT T1/T2/T3 threshold-region recovery | 0/21、0/21、0/21 |
| SHGO / DIRECT raw T2/T3 convergence | 0/21 / 0/21 |
| Raw convergence Gate / dual agreement Gate | FAIL / NOT_REACHED |
| Execution-valid raw case-budget runs | 126/126 |
| Depth branches / hard-budget enforcement failures | 2646 / 0 |
| 既有独立重建检查 | 14479 PASS，0 FAIL |
| 既有独立冷重建 case/depth/state | 232973 |
| 最大 J 重建差 | 1.1619603057511085e-9 dB |
| Tested catalogs exact alias | NO_IN_TESTED_CATALOGS |
| Tractability fresh confirmation | NOT_REACHED；未执行 |

来源：[final decision](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_FINAL_DECISION.json)、[reconstruction summary](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_RECONSTRUCTION_SUMMARY.json)。本轮仅核对保存的表/摘要，没有重新运行 forward reconstruction、原 harness tests 或优化器。

冻结 raw CSV 的直接汇总：SHGO T3 最好 J 的 case 范围为 0.7287017462–4.1630791420 dB，DIRECT 为 0.8582561358–3.9996842622 dB。40/42 个 T3 case 在 0.8–4 dB；最小值也高于 0.001 dB Gate 约 729 倍。它们是已有表的描述汇总，不是新实验或新 Gate。

原始表：[SHGO](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_SHGO_RESULTS.csv)、[DIRECT](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_DIRECT_RESULTS.csv)。threshold-hit 是可行 exact evaluated witness 满足严格 J<0.001 dB。所有候选保持 `EVALUATED_WITNESS_NOT_CERTIFIED_LOCAL_MINIMUM`；没有认证局部极小、捕获盆地或连续空间 coverage。

## 存在正确匹配，与能够找到它，是两个问题

FIX1 self-match 与 FIX2 tested controls 给出正确匹配存在的证据。FIX2 的 P09/Q06 独立族近零恢复、Q01 主族成功而独立族失败，说明不同有限搜索族会漏掉不同区域。FIX2 的 18/21 与本轮 0/21 属于不同搜索结构、预算、聚合方式，不能视为同一算法的预算退化曲线。

本轮结论是：冻结 SHGO/DIRECT observation-only 路线在所测工程预算内没有稳定进入正确阈值区域。它没有建立物理不可辨识、正确解不存在、所有算法不可能或四维全局唯一；也不能外推为 RC3 传播机制失效。

5.37 mm 是 FIX2 的局部 likelihood section，不是捕获域、传感器精度、毫米全局网格要求或已证明的先验精度要求。finite-catalog exact alias 未发现，只适用于导出的有限候选集合。

## 正式停止与路线

`NOT_ENGINEERING_TRACTABLE_UNDER_TESTED_FROZEN_SEARCHES` 仅标记当前 direct-TL global-search 提案在测试/冻结作用域内没有建立工程可实现性，不是普遍计算不可能的证明。

关闭 A1 search-repair chain。禁止 T4、加预算、FIX3/FIX4、新 optimizer、tractability fresh confirmation、A2、depth/B、SSP、P5。旧失败不得删除、重新分类或当作 fresh confirmation。

建议 `R1 > R2 > R3`：先讨论 RC2 temporal candidate narrowing + RC3 local propagation anchor，再讨论 searchable observable；当前 exact matcher保留为 oracle/local diagnostic/conditional rescoring。R1 是待验证路线，必须面对 R3 直航时间累积与小转向的距离脊线负证据，见 [路线决策](R4_ROUTE_DECISION.md)。下一步 route redesign Gate 须另行定义、冻结、授权；本轮不运行、不开放 A2、不增加进度。

## 历史与管理文件

本轮只新增此目录、更新 `R4_MASTER/R4_PLAN.md`、追加 `R4_MASTER/R4_EVIDENCE_LEDGER.csv`。旧数值证据、算法、预算、Gate、release flags 与 manifest 全部保留。

`R4_MASTER/R4_PROGRESS.json` 保留为 9cfa99c development 历史快照；其 pending-third-audit 不再描述当前状态。当前状态读取 [closeout decision](R4_A1_CLOSEOUT_DECISION.json) 和更新后的 plan，进度仍 0。旧 development freeze 中 plan/ledger 的文本绑定按接受的 parent Git blob 验证；本轮授权修改这两个管理文件，不回写旧 freeze。

输入绑定及完整性核对见 [CLOSEOUT_VALIDATION.json](CLOSEOUT_VALIDATION.json)。此前 74 harness tests 和 14479 reconstruction checks 是已有证据，不能写成本轮新跑的检查。
