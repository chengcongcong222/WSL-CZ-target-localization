# H3-G1：接收可提取性与水平支持审查

## 判定

主结论：H3_DEPTH_SUPPORT_AND_PROPAGATION_UNRESOLVED。

并列缺口：
- H3_EXTRACTION_SOURCE_CONDITION_INCOMPLETE；
- H3_EXTRACTION_CALIBRATION_CONDITIONAL；
- PHYSICAL_PROPAGATION_NUISANCE_CONTRACT_NOT_ESTABLISHED；
- FULL_RC2_HORIZONTAL_SUPPORT_NOT_ESTABLISHED。

未授予H3_EXTRACTION_MINIMAL_PILOT_READY的无条件准入。可以提出条件化有限频谱样本试验的冻结草案，不将源/硬件资料未知变为实装验收失败。当前是申请阶段预研设计。

## 保留的正向证据

H3-G0局部深度信息接受，1257公式/1272独立重建/180数值稳定检查及B1/B2增益2.34–4.47倍保持原作用域，不覆盖原数据。主B1的有限表联合MIN即使移除真实水平中心，仍在两种噪声各3/3独立几何排200 m最低；B0/B2均0/3。分数并非近零、不是有限样本似然，也不是完整RC2支持或估计成功。

固定错误水平节点的条件排序仍为B1 62/150、64/150，镜像只对应3个几何。中心包含时的MIN只能视为oracle构造，不重复历史B1B的工程成功误读。详见FINITE_LABEL_PROFILE_SCOPE.md及逐行来源CSV。

## 主要未闭合条件

公共源相位可从同窗空间交叉谱消除，源功率/噪声偏置、共享参考/窗口相关、漂移及低SNR仍需统计模型。无需目标跨三个快照保持绝对相位；需要阵内共钟与跨快照相对通道响应稳定。

13传播频点不证明实际UUV有13条可用频率。当前没有源谱占有、窗内运动/场稳定、噪声空间频率协方差及校准漂移证据。自由模态gain仍吸收深度；没有独立物理依据选择受限组模型。RC2历史粗网格损失、有限连续云、静态预算和动态点估计均不能提供完整H3水平支持。

H12固定底的强信息损失保留；5%噪声可超过最弱接收幅度，不用线性化log/phase把它包装为可靠提取。

## 下一次最小试验及资源

先建议条件化FINITE_SAMPLE_FREQUENCY_DOMAIN_STATISTIC接口试验：三个几何；B0/B1/B2；固定单频200、三频150/200/250和13频资源包；K=1/8/32；1%/5%相对及固定底；独立256个噪声背景样本；16次有限样本实现。完整联合CSD及原始复谱保留共享频率噪声。Gaussian源协方差族与确定性线谱模型分别处理，不能跨模型假称信息损失。

主接口用完整CSD似然，归一化响应为附加带不确定度特征；无法认证的压缩信息不冒充Fisher效率。设计上限5184诊断单元、单次1小时、≤1 GB，需下一阶段独立冻结和授权。本轮只写草案，不执行。

## 管理状态

parent SHA：0cfa6b6c7c6ff231d602522144c5a7b992cf03ad。
review SHA：包含本文件的唯一审查提交；远端核对值由最终回报给出，文件不自引用自己的SHA。
H3-G0：CONDITIONAL_LOCAL_DEPTH_INFORMATION_ACCEPTED。
H3实际提取：NOT_ESTABLISHED。
完整水平支持测深：NOT_ESTABLISHED。
E2：PAUSED_WITH_FORMULA_REPAIR_FAILURE；原E2-G0 FAIL_UNCHANGED。
速度信息瓶颈与原H3-A路径保留，本轮不执行它们。
R4=0%；新KRAKEN/FIELD/MC/audio=0；下一执行NOT_AUTHORIZED。
单次审查提交推送、核对remote/main后停止。
