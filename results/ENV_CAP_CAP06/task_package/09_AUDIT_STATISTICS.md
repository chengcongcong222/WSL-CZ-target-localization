# CAP06 独立审计、统计核对及失败归档

## 为什么单靠一个“PASS计数”不够

以往阶段冷审计对代码复建有价值，但同一模型的两套程序可能共享错误物理假设。新阶段务必分开核查：`PHYSICS_ADMISSION`、`IMPLEMENTATION_RECONSTRUCTION`、`EXTRACTION_IDENTIFIABILITY`、`GLOBAL_SUPPORT`、`GENERALIZATION`五类门。某个数值审计0失败不得替代物理认证。

## 预执行确认

- 源SHA、环境/阵元、信号生成、目标运动、先验边界、增益与噪声模型、留出种子、前端参数、搜索预算、资源配置、Gate和阈值：在Commit A冻结。没有不可逆源身份或真实域先验的试验STOP。
- 合成正控制/零控制可在冻结前开发迭代，保留修补历史；科学代码冻结后每单元运行一次，失败不得看到标签后修改。
- 每个单元`SCIENCE_STARTED.json`存 `unit, codeSHA, inputSHA, configSHA, seed, budget, model_level`。不能因运行超时在旧单元续命加预算。

## 冷审计最低范围

|阶段|必须以独立实现重构的量|关键负面检查|
|---|---|---|
|P0|阵元坐标、源/收/海底深度、波矢符号、单位、零模型Jacobian、数据域读权限|真参数泄漏、远场HLA的虚假r/z信息、镜像符号泄漏|
|P1|功率扣噪、比值、噪声/共参考协方差、投影矩阵、差分增益模型；如通过前向，独立求解数值采样和声功率统计|低能量log偏差、I→CW冒充、差分自由增益非零信息、真环境/z喂估计器|
|P2|阵内复基带仿真与谱估计独立正控制、无符号阵列响应、增益/姿态/导航/时偏的公开—真值隔离、完整错支计数|偷用真实舷侧、真实目标关联、Uniform触界当校准、把旧T3 0.05°当获得能力|
|P3|全部候选/保留/未决/排除、样本分母、拒判、条件和无条件误差、P95及抽样独立性、联合协方差|离网格真值补点、同源派生数据被独立相加、成功者子集替代全部样本|
|P4|从原始P3统计重建资源Pareto、成本口径、各正例对应Gate和数据血缘|解析必要条件偷换经验可行；不同物理证据层级混排|

冷审计代码**不得import主科学判定函数**；可读取注册输入、冻结规范与独立数学公式。记录 `check_id, source, error, tolerance, PASS, interpretation`。对科学不可执行路线，独立核对NOT_EVALUATED而不是生成空的PASS。

## 性能统计

- 量纲：距离绝对/相对（真r非0）、深度绝对m、速度m/s、方位/航向圆周误差deg。每个状态、质量、资源和左右舷给样本数。
- 若确有独立留出，报样本二项事件95%区间（Wilson或Clopper–Pearson）、经验P50/P90/P95/P99并注明对应样本数和P95尾部置信性。真值关联失败和前端拒判必须进入无条件统计，不可仅成功样本P95。
- 同一基础随机数的S0/S1/S2/S3是**配对设计**，不是独立重复；不同质量/状态共用噪声时统计须按配对和种子簇统计。
- 设计与结果区分：无新随机观测仅`MEAN_SEPARATION`; 无非oracle提取仅`ORACLE/TOY`; 只有特定前向可证时可以报告相应层级L5/L6。
- 不因看到漂亮的G02或B3案例而选择P3 holdout；不把历史已见G/H轨迹当fresh。

## 闭环防错

必须生成`EVIDENCE_LEVEL_MATRIX.csv`，每条图和每项结论标注明确的`physical_model_id, input_SHA, processor_SHA, estimator_SHA, statistical_denominator, evidence_level, blockers`。任何无来源数字默认不可进入最终报告。独立审计先验是必要而非充分条件；重大物理矛盾不可“数值重建0FAIL”而忽略。
