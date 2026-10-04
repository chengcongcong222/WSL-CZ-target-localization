# A1 全链条证据收口

阶段：`R4_A1_CLOSEOUT_AND_ROUTE_DECISION`。接受基线：`9cfa99c14bb80d08e6f16210ebf89fd9c541b850`。研究负责人于 2026-10-04 接受第三次独立审计及负结果；本轮仅整理已有证据，新增数值实验为 0。

## 时间与因果链

| 阶段与接受 SHA | 已建立的证据 | 未闭合项与后果 | 原始来源 |
|---|---|---|---|
| Original A1；`4633ff091874bc27d54017a9aff3b87cd4eed97d` | off-grid 生成与 A1-1 管线通过；1629 个 RC2-only case、27 个 nominal pilot 保留 | coarse RC2 会排除真值邻域，coarse acoustic 排序也出现错误；完整 bearing 边界未建立，结构性停止 | [A1 判定](../R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_DECISION.json)、[设计](../R4_A1_OFFGRID_BEARING_BOUNDARY/R4_A1_DESIGN.md) |
| A1-FIX；`956f2dfba8e719561641fd135f246a650fe76dca` | continuous forward self-match、truth-independent continuous RC2、R3 preservation 通过；continuous truth bearing coverage 36/36 | nominal 回归 4/9、Q 留出 8/12；acoustic coverage 失败。V1 初始化缺陷已修复并留档，最终 V2 负结果不能归因于该缺陷 | [FIX 判定](../R4_A1_FIX_CONTINUOUS_SEARCH/R4_A1_FIX_DECISION.json)、[报告](../R4_A1_FIX_CONTINUOUS_SEARCH/R4_A1_FIX_REPORT.md) |
| A1-FIX2；`7ab24845e6e1551b75287fefb1ab662e92b395b8` | matched likelihood sections 与耦合敏感性被表征；独立族在 P09/Q06 找到主族漏掉的近零区域 | 2/6 Gates PASS；raw coverage/convergence、独立一致性、fresh recovery 未闭合；冻结 BLOCKED | [FIX2 判定](../R4_A1_FIX2_ACOUSTIC_COVERAGE/R4_A1_FIX2_DECISION.json)、[报告](../R4_A1_FIX2_ACOUSTIC_COVERAGE/R4_A1_FIX2_REPORT.md) |
| Tractability PRE-RUN；`b4839b295777127ec0c8ade56b76c776db99148b` | 外部 hard-cap、exact objective/API、排除真值输入及计算成本冻结 | 当时没有 noisy development；基础设施通过不等于科学 Gate 通过 | [接受的 PRE-RUN 报告](../R4_A1_SEARCH_TRACTABILITY_AUDIT/ACCEPTED_PRE_RUN_REPORT.md) |
| Harness；`c5c9ba2e1258e73486e62df26499e03cce204f27`；运行授权 `d731cf51d72b41dea7fd4c7a48fcf21a623f7bc4` | 固定 21 noisy cases、独立 raw 分支、witness 语义、Gate 和 16/64/256 caps；仅 development 被单独放行 | fresh 未放行；阈值 witness 不是 certified minimum | [harness 报告](../R4_A1_SEARCH_TRACTABILITY_AUDIT/A1_SEARCH_TRACTABILITY_REPORT.md)、[运行授权](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_RELEASE_AUTHORIZATION.json) |
| Development；`9cfa99c14bb80d08e6f16210ebf89fd9c541b850`；第三次审计已接受 | 两个确定性全局 solver 家族、observation-only、direct-modal exact scoring；126 raw runs execution-valid，2646 branches accounting 闭合 | SHGO/DIRECT 各 T1/T2/T3 都 0/21，各 raw T2/T3 convergence 0/21；dual NOT_REACHED；`A1_SEARCH_TRACTABILITY_NOT_ESTABLISHED` | [运行报告](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_RUN_REPORT.md)、[最终判定](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_FINAL_DECISION.json)、[冷重建摘要](../R4_A1_SEARCH_TRACTABILITY_AUDIT/DEVELOPMENT_RECONSTRUCTION_SUMMARY.json) |
| 本轮 closeout | 接受上述作用域内结论，关闭 search-repair chain，建议 R1 > R2 > R3 | `R4_ROUTE_REDESIGN_REQUIRED`；R4=0%；不启动路线实现或实验 | [收口报告](A1_CLOSEOUT_REPORT.md)、[路线决策](R4_ROUTE_DECISION.md)、[机器判定](R4_A1_CLOSEOUT_DECISION.json) |

因果链：off-grid/coarse selection blocker → continuous RC2 修复 → acoustic coverage 仍未闭合 → FIX2 窄截面与搜索族分歧 → 冻结 SHGO/DIRECT 未捕获阈值区域 → tractability 未建立 → A1 search-repair chain CLOSED。

## FIX2 的数字应怎样读

| 项目 | 冻结结果 | 解释边界 |
|---|---:|---|
| B1/B2/B3 cumulative nominal recovery | 9/21 → 17/21 → 18/21 | B3 raw 只有 17/21；累计保留不是独立收敛证据 |
| B2、B3 各自 raw 均恢复 | 16/21 | 不是 18/21；coverage/convergence Gate 未闭合 |
| 独立族 recovery / 与主族 agreement | 5/7 / 3/7 | P09、Q06 独立成功而主族失败；Q01 相反；P02 两者均失败 |
| FIX2 fresh noiseless / nominal | 8/8 / 13/16 | 24 cases 已执行；科学 Gate 未闭合，不是程序没跑完；三项含噪失败保留 |
| 开发 noiseless | 15/15 | `NOISELESS_RC2_UNIQUE` / `RC2_NUMERICAL_SINGLETON; NOT_ACOUSTIC_BASIN_WIDTH`；主要是 regression/control |
| J<0.001 dB 距离截面宽度 | 约 1.08–21.2 mm；median 5.37 mm | 局部 likelihood section，不是 attraction basin、全局网格要求或测距精度 |
| 最快归一化方向 range loading | median absolute ≈0.9961 | 距离主导并与速度耦合；独立轴截面不能拼成四维捕获域体积 |

FIX2 fresh 曾执行，tractability fresh 没有执行，不能混用。FIX1 的 Q holdout 在后续方法开发中已成为 development，不能再次承担新方法确认作用。

## 三个不同层面

1. **物理/信息：**tested synthetic matched controls 中存在正确 acoustic match，包括独立族捕获的近零区域。物理不可辨识未建立；全局唯一也未建立。
2. **搜索可实现性：**当前 observation-only direct-TL matcher 在测试的冻结 SHGO/DIRECT 预算内未建立稳定 threshold-region capture。限定于模型、panel、搜索结构、预算及 Gate。
3. **工程路线：**研究负责人关闭继续加密网格、增加随机起点、增大 SHGO/DIRECT budget、增加 FIX 的支线。这是基于负证据及冻结停止规则的路线决定，不是所有可能算法都不可能的定理。

旧数值结果、报告与 manifest 保留。历史报告中的 pending/next recommendation 是当时状态，由本轮决定向前覆盖，不回写历史。
