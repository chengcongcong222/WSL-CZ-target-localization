# PREC2研究来源与状态

- 私有项目原始记录：`chengcongcong222/WSL-CZ-target-localization`，固定父提交`ac43a8702282fc6c80663319cb6757a30b084e09`。
- 历史 N2：`results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/SCIENCE_OUTPUT_FREEZE.json`、`SOLVER_PLAN.json`、`modal_archive/n2_K{2,3}_f{201,235,283}.{mod,prt}`、`environments/n2_K{2,3}_f{201,235,283}.env`。本任务机器JSON列出所有18个原始文件SHA256。
- 历史 PREC1：`results/R4_T2_N3_PREC1_MOD_PRECISION_ATTRIBUTION/GPT_SYNC_PREC1.md`（仅理想等声速验证，有别于本轮）。
- 外部公开源码：https://raw.githubusercontent.com/oalib-acoustics/Acoustics-Toolbox/main/Kraken/kraken.f90 。源码文本显示双精度本征值与`CMPLX(k)`写MOD，但没有证明此文件与封存Windows exe为同一构建，不能宣称已取得可执行文件的完整内部双精度侧车。
- 外部官方说明：https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html 。KRAKEN对吸收扰动处理，KRAKENC复数运算；本轮选择损耗剥离后的自伴实特征值诊断，刻意**不**认证完整复特征值的虚部及实际传播场。
- 主动防止混淆：GitHub连接器只读PRT/MOD抽样核查得到184/184是**预研**；用户本地正式主分析、独立冷审计尚未运行。附的CSV不得作为正式程序的输入；程序须直接读取历史二进制。
- 不声称真实海域传声损失、CZ边界、左右舷分辨、五参数精度、硬件优劣。R4=0%。
