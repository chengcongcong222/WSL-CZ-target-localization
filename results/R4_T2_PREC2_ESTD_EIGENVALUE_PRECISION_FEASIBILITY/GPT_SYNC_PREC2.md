# T2-PREC2 负责人同步

设计SHA：`e7b19f081f587d5c815a6e0b9c5eb3c44a0d5aa5`。执行SHA由最终外层回传。父提交：`ac43a8702282fc6c80663319cb6757a30b084e09`。
任务MD SHA256：`f8eef00f759af543b6cb58199fdf255d0a696ca49ed0bd8401fe0070bc856efe`。更正SPEC SHA256：`71acc2be68c503334f8842472406936dfaf2799440c05c8920cfafb8c02ba16d`。

判定：`PREC2_SPARSE_OUTPUT_LOSS_AND_SUBSET_LOSSLESS_REAL_K_CORROBORATED`。主程序PREC2_DECISION.json原样封存；最终分类见PREC2_FINAL_ADJUDICATION.json。

仅更正获批的源SHA登记及受影响元数据哈希链，原/更正ZIP、原SPEC、首次STOP、从头110项预检记录均随提交留存。原任务MD及科学代码未改变。

|来源|模态数|打印行|逐位float32匹配|60km最大相位差(rad)|
|---|---:|---:|---:|---:|
|n2_K2_f201|980|31|31|0.00174405286|
|n2_K2_f235|1146|31|31|0.00178305103|
|n2_K2_f283|1380|30|30|0.00336368041|
|n2_K3_f201|980|31|31|0.00174405286|
|n2_K3_f235|1146|31|31|0.00178305103|
|n2_K3_f283|1380|30|30|0.00336368041|

P0：184/184条、18/18源文件；表示关系可重建，不是内部全部波数精度证书。两网格合计7012个MOD记录是重复网格视图。
P1：15/15个冻结模态，R4内部收敛标志通过15/15；最大R4相邻外推根差1.2820033923333085e-10，设计预算1.6666666666666667e-09。

|201Hz固定模态|两档Prüfer根差|Prüfer细根–FDM R6|预算内一致|
|---|---:|---:|---|
|1|2.72226686e-13|1.77191595e-13|True|
|481|1.59045777e-10|5.90520965e-11|True|
|961|2.01444861e-10|1.12072573e-11|True|

正式冷审计：766项，0FAIL；3个预登记IVP根实际完成3。零完成不称通过。最大IVP/FDM差5.90520965459973e-11。
主程序1次、冷审计1次，重试0；科学链壁钟182.281s，峰RSS82.14MiB，新增结果381581字节。禁止的声学求解、MC、波形、DFT、FIM与候选反演均0。

异常/停止原因：无科学执行异常；登记流程完成后STOP。

本轮仅提供稀疏输出精度损失证据和有限无损实根数值一致性。Thorp虚部精度、全部有损复波数、完整声场、Windows exe高精度全模态输出能力、CHIGH及真实海底均未认证。

N2/N3/PREC1作用域及原G0/N1历史结论保持原样；T2-A/B实际接收性能NOT_EVALUATED，R4=0%。负责人审计PENDING。STOP。

有限实根若获得独立一致性，只建议研究负责人另行设计有来源的高精度全模态输出/侧车与完整复波数验证，不能直接替换虚部或旧声场。

公开源码与旧exe构建未溯源；模态窗口和实际底质合同独立保留。全模态侧车、E-STD完整复场、海底替换、接收信号和定位性能均未授权。R4=0%。STOP。


Delivery note: 766 actual cold checks = 18 hashes + 6 mode counts + 4x184 scalar checks + 6 IVP checks. Prior package synthetic tests are not the actual archive audit. Runtime result bytes are sampled before report and manifest generation; final delivery remains below64MiB. Execution/remote SHA are supplied after Commit B, avoiding a self-referential commit identifier.
