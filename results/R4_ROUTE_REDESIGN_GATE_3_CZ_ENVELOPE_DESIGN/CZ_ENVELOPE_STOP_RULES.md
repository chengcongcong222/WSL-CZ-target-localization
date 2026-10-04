# CZ envelope 停止与判定规则草案

当前是设计阶段，所有运行未授权；以下仅供负责人pre-run选择并冻结。

| 停止项 | 未来检查 | 判定与处理 |
|---|---|---|
| S1 retention | joint truth/reference不能在完整冻结主case清单保持；projection-only保留不足 | `CZ_ENVELOPE_COARSE_RANGE_INFORMATION_NOT_ESTABLISHED`；保留失败，不按真值扩大band后重算 |
| S2 no contraction | 所选primary档不满足actual range-domain contraction，或仅candidate数量下降而hull仍15km | 同上；不能以平均/最优case遮盖，也不能失败后降档 |
| S3 aliases/search | coarse representation仍超密oscillation/alias、两级raw support不稳定或global bounds/cap未闭合 | 不声明searchable；有效运行但覆盖未闭合记`SEARCH_COVERAGE_UNCLOSED`，不是物理不可能 |
| S4 oracle dependency | range信息只在truth motion/true depth固定后存在，或使用source/phase/oracle seed | 当前cold-start机制未建立；任何泄漏作为protocol-invalid另标，停止并审计 |
| S5 nuisance absorption | 冻结N1一加入range特征被吸收，或只能用不现实窄drift prior保住 | 对N1合同不成立；N0最多保留ideal-control范围，不包装robust PASS |

主headline依据预注册primary representation+tier+nuisance范围。若选择N0-only first gate，N1缺D则必须写`ROBUSTNESS_NOT_ESTABLISHED`，不能称源漂移稳定；N2人为无限自由度是符号limit，不对真实target作普遍否定。

执行完整性、bound正确性或数据映射坏掉时记`EXECUTION_INVALID / PROTOCOL_INVALID`；不能把软件错误当负科学结论。修复也须在独立审计/新授权下，冻结所有失败证据；不偷偷增加预算、容差或representation。

每case严格hard cap计入所有prediction、depth/drift profile、validation、reference开销，超限阻断并保留unknown cells。两个预注册raw grid levels结束后停止；不T4式增预算、不第三solver、不连续smoothing sweep，不把sensitivity胜者换成primary。不能重开A1/FIX3/FIX4/B1/single-depth/fringe/TDOA/Doppler/spatial实验。

若未来闭合，仅写 `CZ_ENVELOPE_MATCHED_COARSE_ACQUISITION_INFORMATION_ESTABLISHED_IN_FROZEN_SCOPE`，仍需fresh确认、真实接收/信号验证及负责人后续阶段选择。R4管理信用由负责人另行审计，不由draft或单轮机制成功自动增加；当前R4=0%，A2/depth/SSP/P5 UNOPENED。

本轮交付推送即停止，不运行本文的任何检查或新科学试验。
