# GPT sync: Gate-3A PRE-RUN

研究负责人已接受 Gate-3 设计 e56ba746ebe7370d1dee3d0a181bc348b5fecc5f。本轮只实现并冻结 WINDOW_SHAPE_Q2、N0 和 NOMINAL PRE-RUN；未运行 21-case feature/search。

校准政策先以 a393c7b45c92c42fd148d9c162faeb5ea35160c7 独立提交并推送。固定四个域角点、三个 depth 标签及解析单测；tau_F 只按预注册公式产生，最终为 1e-7 dB。它是数值重建 allowance，不是实测 acoustic noise。

Arb/Acb 全单元包络对固定实数模态模型成立，允许 loose bounds。数学保守性与开发集覆盖/收缩是两个不同检查项：本轮没有证明 coarse/fine 网格能在预算内闭合，也没有 CZ 科学 PASS。

预算预占和评价守卫补强后，在最终代码上按原政策重新校准；最终 full-cell cost 约 3.275 s，按事前公式导出 coarse/fine 请求上限 14/29。首次实现校准的 16/32 未进入最终 freeze；没有改变公式、夹具、阈值、表示或源数据。两次校准均未接触开发观测。

最终 47 项单测通过，包含 20 类要求。未知/松弛/耗尽单元保留；coarse/fine raw 各自从全域、新队列、新缓存、新计数开始。B0 不调用声学；空支持失败；评价只能在 estimator 之外运行。开发入口保持 DEVELOPMENT_RELEASED=False。

PRE-RUN 判定见 CZ_PRE_RUN_DECISION.json，输入、代码、测试、tau、grid、budget 和停止规则绑定见 CZ_PRE_RUN_DESIGN_FREEZE.json；完整性检查见 CZ_PRE_RUN_INTEGRITY.csv。

提交推送后立即停止，等待独立审计。R4=0%；A1 搜索修复支线保持关闭；fresh NOT_GENERATED；N1、A2、depth development、SSP、P5 均不开放。真实接收机 M1 与真实 UUV 稳定多线证据仍 UNKNOWN/NOT_VERIFIED；5 km 不表示最终测距 <10% 或 RC3 capture。
