# T2-PREC2 研究简报：非均匀SSP的水平波数有效精度与来源验证

## 研究衔接

基线 remote/main：`ac43a8702282fc6c80663319cb6757a30b084e09`。旧T2-G0的`FORWARD_NUMERIC_NOT_ADMITTED`、N2的`T2_N2_NUMERIC_BRACKETING_NOT_ESTABLISHED`、N3的`T2_N3_PARTIAL_NUMERIC_NOT_ADMITTED`不改判。PREC1仅支持等声速理想波导下的MOD波数单精度损失机理。

原E-STD是5000 m水柱、251个C-linear SSP节点、压力释放海面V、理想刚性海底R、`CVWT`/Thorp、恒密度。MOD储存的复波数同时包含实部与扰动吸收虚部。不能将无损高精度实波数认作完整有损复波数，也不能拿旧1800窗口与2400窗口混合作为一组真值。

## 本会话已新增的直接证据（非Codex正式冷审计）

已从GitHub父提交直接读取N2的K2(160001)和K3(240001)两网格、201/235/283Hz的六套PRT与MOD：模态总数分别980、1146、1380；两套各打印31、31、30个PRT抽样模态。**184/184**个PRT实波数在舍入为`float32`后与MOD实部逐位一致；全部抽样PRT实波数本身与MOD实部不等。按60 km换算，三频最大相位差量级分别0.001744、0.001783、0.003364 rad。相应原始Git输入、SHA256在机器冻结配置中明确登记。`PRT`并未打印全部模态，不得对缺失模式作插值补全。

公开KRAKEN源码曾使用双精度本征值及`CMPLX(k)`单精度写MOD的步骤；这是来源机制解释，但没有源码到当前Windows exe的构建证明。官方使用说明指出KRAKEN的吸收以扰动处理，KRAKENC在复域计算；只能在严格标明“无损自伴实特征值”的作用域内检验独立数值方法。

## 新一轮冻结研究范围

1. **PREC2-P0：六套既有PRT/MOD/ENV的全量SHA与逐位舍入身份**。输出全部抽样模态对照，不将PRT当作全模态侧车文件。
2. **PREC2-P1：独立无损实波数小样本诊断**。对每个频率，从K2 PRT抽样模态中按冻结列表位置取五个模态；采用恒密度Dirichlet/Neumann Sturm–Liouville算子，四个垂向离散密度`N=20k/40k/80k/160k`，在lambda=kr²层做Richardson外推。记录R4两层的差，不能把外推收敛指标称为严格误差上界。
3. **PREC2-P2：独立冷审计**。重新读取历史MOD/PRT，用另一种PRT/MOD解析实现校验全部抽样；在201Hz选取三个冻结模态，以Prüfer相位积分/DOP853和根求解检验FDM估计，分别记录算法步长/容差依赖和不一致样本。

设计工程指标：距离60000m时单模态相位绝对扰动不超过`1e-4 rad`，对应`|delta kr| <= 1.6666666667e-9 m^-1`。这是人工冻结的**波数精度设计预算**，并非阵列0.2%空间形状误差的充分条件。独立小样本达标，也不能保证全部980/1146/1380个模态精确。

## 结果与治理

- **禁止**：新的KRAKEN/KRAKENC/FIELD/FIELDS/SCOOTER、重跑旧G0、声压/波形/DFT/MC/FIM、低场样本剔除、任意改变海底条件、把已有complex64强制转换为complex128后称作高精度。
- **保留负结果**：任何指定模态未收敛、PRT/MOD舍入不一致、常速或源文件校验异常立即记录并按冻结Gate STOP，不事后增网格、更换排名位置或改变容差。
- **下一步**：只有这次明确建立“有限实模态高精度取得路径”，才能另行考虑经源码溯源的全模态双精度侧车，或独立求解完整复波数。CHIGH与真实底质仍需独立治理。
- **论文主张**：只能报告“存档输出表示损失的实证”与“非均匀名义SSP部分实波数的条件数值一致性”；T2-A/B接收性能继续`NOT_EVALUATED`，R4=0%。

参考：KRAKEN公开手册 <https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html>；公开源码 <https://raw.githubusercontent.com/oalib-acoustics/Acoustics-Toolbox/main/Kraken/kraken.f90>；历史N2/N3/PREC1按父提交锁定。
