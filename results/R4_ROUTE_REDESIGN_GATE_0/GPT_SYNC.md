# GPT 同步：R1 Gate-0 冷启动信息审计

已接受 parent closeout：`206f6c2f46bc3c52890a32004cf54255ed30fcef`。本文件所在提交为 Gate-0 SHA；提交后核对 local HEAD == remote/main。

结论 B：`R1_COLD_START_NOT_ESTABLISHED`；`R1_TRACKING_ONLY_CONDITIONALLY_PLAUSIBLE`。cold-start information source identified=false。circular dependency=true，限定为 acquisition 依赖 previous posterior，而该 prior 最终又依赖已失败 global matcher 的 bootstrap 分支；不是现有实现 bug，也不代表独立取得 support 后的 tracking 必然循环。

关键既有证据：

- RC1 support role 明确 45:1:60 km 条件域，没有更宽 pre-RC1 对照，contraction ratio 未识别；无更窄 first-CZ prior admission。
- R3-1B：1326→495 candidates，range width 仍 15 km。R3-1C 所测五转角 0/2/5/10/15° 的 RC2-only width 都 15 km。
- FIX1 continuous RC2：15 noiseless controls full rank、sampled width 0；21 noisy cases truth retained，但 sampled range span 14.9008764373–14.9987589195 km，median 14.9734269641 km。有限 samples 不是 confidence certificate；不能否定 bearing 在理想条件下的 range information。
- R3-1C 部分 turn+RC3 coarse reference 有正证据；主判与后来 off-grid noisy failure 均保留，不能作为隐含初始 range prior。
- FIX2 local matched information 强，但 attraction/capture 证书未建立；SHGO/DIRECT 各 budget 0/21 的宽域失败不产生 first posterior。

G0-1/2/3 未闭合；G0-4 可定义 support retention + domain contraction，G0-5 可限定 local role，二者只是定义层面，不是验证通过。R1 acquisition 不放行 quantitative design 或 implementation；R1-B 只有来源合法、可信且足够局部的 previous support 才条件性合理，claim 限为 POST_ACQUISITION_TRACK_REFINEMENT。

推荐负责人下一轮定义 **R2 acquisition observable/information-source 文档设计**。R1 留作条件 tracking/refinement 子模块，Route R3 留 oracle/local diagnostic；R2 没有因此获得实现或实验授权。

审阅：[信息来源表](COLD_START_INFORMATION_SOURCE_MATRIX.csv)、[acquisition/tracking](R1_ACQUISITION_VS_TRACKING.md)、[循环依赖](R1_CIRCULAR_DEPENDENCY_AUDIT.md)、[Gate-0 报告](ROUTE_REDESIGN_GATE0_REPORT.md)、[机器判定](R1_GATE0_DECISION.json)、[完整性](VALIDATION.json)。master plan 更新，ledger 仅追加，不改旧 closeout 或历史证据。

新数值实验 0；R4=0%；A2/depth/SSP/P5 UNOPENED。提交后停止，等待下一轮明确指令。
