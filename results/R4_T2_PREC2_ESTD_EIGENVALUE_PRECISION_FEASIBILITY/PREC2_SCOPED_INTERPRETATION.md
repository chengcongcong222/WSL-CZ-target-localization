# T2-PREC2 科学作用域

本轮只核验原始E-STD档案的输出表示关系，并计算原251节点C-linear SSP、恒密度、水深5000m、V海面/R刚性底下的无损自伴实特征值小样本。原ENV CVWT/Thorp保持原字节，无损算子明确剥离吸收；未将结果写回旧MOD或合成有损场。

最终限定判定：PREC2_SPARSE_OUTPUT_LOSS_AND_SUBSET_LOSSLESS_REAL_K_CORROBORATED。18份源输入、184条稀疏PRT打印数值、15个预登记FDM根和3个固定Prüfer根分别记账，不混成全模态认证。实际完成数为18/18源、184/184打印记录、15/15FDM根、3/3IVP根（完成/预登记）。

单模态设计预算为1.66666666667e-09 m^-1，即60km时1e-4rad；该数值是冻结设计量，不是空间形状0.2%门槛的充分条件。Richardson网格一致性为内部诊断，DOP853两档及另一路连续算子检验仅覆盖3个根。没有严格余项界、全部有损实/虚部精度证书或真实海洋传播认证。

PRT仍是每套30/31条稀疏打印，不能插值填满980/1146/1380模态，也不能将complex64转complex128解释成恢复丢失位数。原MOD虚部只保留并记录打印/存储数值和符号对应；完整虚部精度未评价。

源码上下文：[公开KRAKEN源码](https://raw.githubusercontent.com/oalib-acoustics/Acoustics-Toolbox/main/Kraken/kraken.f90)有CMPLX(k)写盘；[官方手册](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html)区分扰动吸收和复数求解。这些资料不建立封存Windows二进制与公开源码的构建等价关系。

历史保留：N2扩窗空间形状差异最大99.9902%及98.7434%，网格0.3712%；N3仅2次解析配置，2400窗P04最大0.22149%触发停止；PREC1的40组理想波导归因与1787波数舍入证据保持。旧N1优势证书只限旧1800窗口。旧G0/N1/N2/N3/PREC1表格和代码未改变。

CHIGH全域完整性、真实海底、全部复模态、T2-A/B接收与定位表现继续未认证或NOT_EVALUATED。R4=0%。STOP。
