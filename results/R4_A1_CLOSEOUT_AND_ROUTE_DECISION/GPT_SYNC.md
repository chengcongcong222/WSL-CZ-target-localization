# GPT 同步：A1 收口与路线决策

parent accepted development SHA：`9cfa99c14bb80d08e6f16210ebf89fd9c541b850`。本文件所在 Git commit 即 closeout SHA；不在正文嵌入自身 SHA。推送后核对 local HEAD == remote/main。

- A1 search tractability：`NOT_ESTABLISHED`；第三次独立审计已接受。
- A1 search repair chain：`CLOSED`。
- physical non-identifiability / global uniqueness：均 `NOT_ESTABLISHED`。
- 当前 direct-TL global search：`NOT_ENGINEERING_TRACTABLE_UNDER_TESTED_FROZEN_SEARCHES`，保留 tested/frozen scope。
- 推荐：RC2 temporal candidate narrowing + RC3 local propagation anchor；`R1 > R2 > R3`。
- R4 route status：`ROUTE_REDESIGN_REQUIRED`；R4=0%；A2/depth/SSP/P5 UNOPENED。
- 本轮新数值实验：0；无 optimizer、simulation、forward reconstruction、fresh generation 或 harness test 执行。

SHGO/DIRECT 各 T1/T2/T3 0/21；各 raw T2/T3 convergence 0/21；dual NOT_REACHED。126 raw runs execution-valid，2646 branches 闭合，hard-cap failures 0；既有 14479 reconstruction checks PASS/0 FAIL，232973 unique state/depth 冷重建，max J delta 1.1619603057511085e-9 dB。finite catalogs 无 exact alias；tractability fresh 未执行。

FIX2 保留 BLOCKED：2/6 Gates PASS；cumulative 9/21→17/21→18/21，B3 raw 17/21，both-raw 16/21；independent 5/7、agreement 3/7；已执行 FIX2 fresh noiseless 8/8、nominal 13/16。5.37 mm 是局部 likelihood section，不是捕获域、精度或网格要求。

负责人下一轮决策重点：另行定义 route redesign Gate。R3-1B 直航累积与 R3-1C 所测小转向 RC2-only 均未缩窄 15 km 距离域，因此 R1 必须回答真实缩域信息来源、truth-free 冷启动及 candidate coverage；不能预设精确 prior/previous posterior，不能循环依赖失败的 global search。R2 只定义 observable 方向，相位不可假定来自现有 relative-TL 幅度数据。R3 仅为 oracle/local diagnostic/conditional rescoring；三条路线均未实现或验证。

审阅入口：[证据链](A1_EVIDENCE_CHAIN.md)、[报告](A1_CLOSEOUT_REPORT.md)、[路线与前置 Gate](R4_ROUTE_DECISION.md)、[机器判定](R4_A1_CLOSEOUT_DECISION.json)、[完整性检查](CLOSEOUT_VALIDATION.json)。master plan 更新，ledger 只追加事件。旧 R4_PROGRESS.json 与 freeze 保留为 9cfa99c 阶段快照，不回写历史。

提交后停止。无 T4、加预算、FIX3/FIX4、新 optimizer、fresh confirmation、A2、depth/B、SSP、P5 授权。新明确指令到达前不启动路线实验。
