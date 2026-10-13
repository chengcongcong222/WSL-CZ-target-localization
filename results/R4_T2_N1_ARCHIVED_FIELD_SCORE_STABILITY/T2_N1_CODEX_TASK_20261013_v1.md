# T2-N1 — 已封存双网格声场的候选评分敏感性与阵型比较数值证据

**任务状态：新一轮已授权的“一次只读数值诊断”实施任务，尚未执行。** 任务基线为 `chengcongcong222/WSL-CZ-target-localization` 的 `6b60c39f52b68a361b2c005e74126b9f68401b77`。本任务不取代、不重新运行、不重新准入 T2-TILT-G0；新目录：`results/R4_T2_N1_ARCHIVED_FIELD_SCORE_STABILITY/`。所有大型临时中间量留在 `D:\ProjectStorage\WSL-CZ\T2_N1`，不纳入 Git。

## 0. 研究问题与可允许结论

仅回答：T2-G0双网格同一静态候选的有限空间形状差异（最坏 `0.0073441923068232945`；冻结门槛 `0.002`），**是否足以改变已登记的匹配分数、远距离错误候选的排序及两项阵型对照的有限数值结果**？

研究对象是“存档数值表示与候选排序的敏感性”，**不是阵型真实性能**。不能把高网格 `160001` 当作物理真值；它只作为一个固定参考数值表示。两个网格可以共同遗漏高相速度模态，原 `CHIGH=1800` 上界完整性 `UNCERTIFIED`；两网格吻合不能证明物理正确。`FORWARD_NUMERIC_NOT_ADMITTED`、A/B `NOT_EVALUATED`、H3-P1 SSP局部性失败及R4=0%保持原样。实际接收、源谱、波形、DFT、噪声、校准、判决覆盖、定位精度和动态五参数**均为NOT_EVALUATED**。不引入任何5km辅助节点或新的垂向结构。

本轮只读取已有 `A0_0, A_P05, A_M05, A_P10, A_M10, A_P4M, A_M4M, B1_V4, B2_DUP` 九资源在 `primary/fine` 两个完整目录和 `80001/160001` 两个模态离散网格下的**已存档复场副本**。T2-A主对照固定为 `A_P4M` 对 `A0_0`（同8通道）；T2-B固定为 `B1_V4` 对 `B2_DUP`（同12通道）。其他倾斜正负角仅作为敏感性对照，不得事后择优替代主方案。

## 1. 不可变来源和身份冻结

先 `git fetch` 并检查 `git rev-parse HEAD` 与 `git ls-remote origin refs/heads/main` 都等于上述parent；若远端或工作区出现相冲突的新提交，**STOP并报告，不改写别人的结果**。现有旧T2-G0文件只读。

必须依次读取/绑定：

1. `results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/GPT_SYNC.md`, `DECISION.json`, `FORWARD_DECISION.json`, `FORWARD_CATALOGUE_NUMERICAL_STABILITY.csv`, `FORWARD_WEAK_FIELD_AND_FAILURE_COUNTS.csv`。
2. `results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/OUTPUT_RECOVERY_INPUT_SEAL.json`, `LOCAL_ARCHIVE_INDEX.csv`, `REPLICA_FREEZE.json`, `PROVIDER_MODE_WINDOW_AND_MESH.csv`, `FIELD_STORAGE_QUANTIZATION_COLD_DIAGNOSTIC.csv`。
3. `results/R4_T2_TILT_VS_SHORT_VERTICAL_G0/CATALOGUE_primary.npy`, `CATALOGUE_fine.npy`, `ARRAY_GEOMETRY_TABLE.csv`, `T2_TILT_G0_EXPERIMENT_SPEC_20261013_v1.json`。
4. `t2_tilt_g0_core.py`中 `shape_error` 与 `score/scores_many`；只读，不允许修改旧函数、来源CSV或旧结果判定。
5. D盘 `T2_TILT_G0/replicas/{RESOURCE}_{CAT}_n80001.npy` 与 `_n160001.npy` 共36个复场文件，以及 `T2_TILT_G0/SHAPE_ERROR_{RESOURCE}_{CAT}.npy` 共18个旧误差数组。精确路径由旧存档索引和源文件核查，不凭文件名推测新副本。

前置只允许进行**文件存在、尺寸、dtype/shape、头信息、原SHA256及目录同一性检查**；不得在Commit A前预先计算新评分、排序、探针优势或新摘要。必要的大文件使用 `np.load(...,mmap_mode='r',allow_pickle=False)`。输入不足或SHA不符则登记 `ARCHIVE_INPUT_INCOMPLETE_OR_MUTATED`，STOP，禁止按真值、旧估计值或者新传播求解重建缺失副本。不要调用用户以往传输目录中名为任务包的旧启动文件。

原登记实际指标：`9*(3355+12705)*3=433620` 个逐频空间向量，其中 **8012** 个 `e>0.002`；54/54资源-目录-频率最大值条件未通过，最大 `0.0073441923068232945`。新的逐点审计必须复算并检查这些历史计数，不能偷换为“54/54个向量都失败”。原FIELD 201/235/283 Hz 原登记0.1526%/0.2252%/0.2170%，保存坐标冷诊断0.1169%/0.1771%/0.1142%，二者只可旁栏列示，不得代替两网格认证或改判FIELD原检查。

## 2. A提交：冻结任务合同、实施代码与预算

任务文件字节SHA256在交接时写入 `TASK_SOURCE_LOCK.json`。在A前实现下述诊断、只读独立审计和确定性测试，但不读取指标矩阵做正式排名。登记 `N1_CONFIG_FREEZE.json` 包括以下全部预注册项目、确切输入SHA、复数归一化公差和输出字段；运行只读前置测试后**一次Commit A并推送**，保存全SHA。任一A后的核心公式、分组、探针、阈值、排序判据修改须STOP并告知研究负责人，不能自行另开科学轮。

本轮预算：一台机器、单个科学进程，线程上限4，科学诊断墙钟 ≤1800s、峰RSS ≤4 GiB、新增D盘 ≤2 GiB、新增Git内容 ≤64 MiB；前置公式自检（不读旧大数组）≤20个小维度确定性矢量；**新增KRAKEN/FIELD/BELLHOP/海洋传播求解=0，新增MC/源噪声抽样/合成录音或复谱=0，优化器/FIM/真实观测生成=0。**运行中预算触发即保存 `PARTIAL_BUDGET` 和已完成项目，不能改变预设分母，STOP。

任务中的固定探针、误差分箱、半径和统计规则不允许在看数据后调整。`N1_CONFIG_FREEZE.json` 必须显式声明：**原 `0.002` 是历史G0前向准入门槛，N1既不改变门槛，也不设置替代认证阈值。**

## 3. 一次数据域诊断B：空间形状误差分布（不删任何候选）

每个阵型、目录、频率，对相同候选索引的两个已存复声压向量计算

\[
 u_c^{(H)}=h_c^{(160001)}/\|h_c^{(160001)}\|,\quad
 u_c^{(L)}=h_c^{(80001)}/\|h_c^{(80001)}\|,
 \qquad e_{c,f}=\sqrt{\max\{0,1-|[u_c^{(H)}]^H u_c^{(L)}|^2\}}.
\]

`e` 是**经过每频公共复比例消除**的射影/空间形状差；原存档的 `SHAPE_ERROR` 必须逐候选复核，独立双精度求法的最大绝对差不大于 `1e-10`（极低差异处应使用稳定的 `||h_H-c h_L||/||h_H||` 求法，避免 `sqrt(1-corr²)` 在e≈0的消减误差；见独立审计）。任何全零、NaN、Inf向量保留其ID并标模型输入异常，不得删除或置0。

完整报告分组：

- 9资源 ×2目录 ×3频，N、max、min、median、P90/P95/P99、超过`0.002`的数量与比例；每资源/目录聚合3频的向量数，不将候选或频点当独立目标几何。
- 按每个资源/目录/频率的高网格向量范数分10等计数档（相等分位以候选ID稳定破同分）记录超限数、误差分位及范数分布；**弱场一律不删除**。
- 固定分层：距离 `[45,50),[50,55),[55,60]` km，`abs(theta)=0` / `0<abs(theta)<5` / `=5` deg，源深标签 `180/190/200/210/220m`，并保留逐频及九资源对照。
- 记录 `e` 与 `log10(||h_H||)` 的Spearman等级相关（仅为非因果描述）、每组最坏5个候选的完整ID与(r,|theta|,z,f,范数,e)。不额外筛样本，也不宣称弱场必然导致误差。

**必须区分**：此前54组“组内最大值失败”与8012个逐频向量超限；若全量计算不能复现这两项，科学部分STOP且错误不写成阵型无效。

## 4. 理论诊断：匹配评分对数值声场差异的通用上界

原T2-G0的源系数剖面化空间匹配，定义单位观测向量 `w_{k,f}`（任意但固定、可以含噪声和未知源）和归一化模板向量 `u_{c,f}^{(m)}`：

\[
 S^{(m)}(c\mid w)=\frac1{3K}\sum_{f=1}^3\sum_{k=1}^K|[u_{c,f}^{(m)}]^H w_{k,f}|^2,
 \quad D^{(m)}=1-S^{(m)}.
\]

在相同**固定观测向量**与两网格候选下，单位秩一投影满足

\[
\|u_Hu_H^H-u_Lu_L^H\|_2=\sqrt{1-|u_H^Hu_L|^2}=e_{c,f},
\quad |S^{(H)}(c\mid w)-S^{(L)}(c\mid w)|\le E(c):=\tfrac13\sum_{f=1}^3 e_{c,f}.
\]

该界对所有单位复`w`成立（包括跨时间段的不同w）；`K=8`并不能让最坏界自动除以√8。对 `D` 同样成立。此处 `E(c)` 是**跨两张已存数值网格的评分敏感性上界**，不是观测噪声标准差、不是海洋传播误差分布，也不是向已通过旧G0的阈值转化的理由。若实际接收重新生成使观测w也变化，此界本身不覆盖w变化；更不界定两个网格共同遗漏的CHIGH模态。

实施时分别从存档 `SHAPE_ERROR` 和独立Gram/SVD方法算 `E(c)`；数值零附近用稳定内积/投影残差路线。严格复核 `E`上界涵盖每个固定探针的全部候选评分差。输出：`SCORE_PERTURBATION_BOUNDS_BY_RESOURCE.csv`与全目录候选E的P50/P90/P95/P99/max、`E>0.002`数量（这是额外描述，**不改旧门槛**）。

## 5. 锁定的无噪声模板探针：错误距离候选排序（不是接收性能）

无需生成声音或另行求传播；每个探针的“观测”只是选取**高网格档案内的一个候选模板**，不是物理测量、不是目标真值、不是新噪声蒙特卡洛。5个固定探针按 `(r km, absθ deg, z m)` 顺序：

1. `P01=(45.5,0.0,180)`；2. `P02=(50.0,0.0,200)`；3. `P03=(55.0,2.5,200)`；4. `P04=(58.0,4.5,210)`；5. `P05=(59.5,5.0,220)`。

全部必须在 `primary` 与 `fine` 两个**独立目录**均有且只有一个精确候选ID；禁止靠最近真值补标签。分别对9资源×2目录×5探针：令固定观察`w_f=u_{p,f}^{(H)}`（八段完全相同，仅为代数模板，无采样含义），用两张数值网格的完整候选计算 `S^{(H)}_p(c)` 和 `S^{(L)}_p(c)`；tie用严格评分、后较小candidateID，记录最高峰、距离排序变化、评分变化和最危险远距离对手。

对每个探针，**事前固定**错误距离目录 `C_wrong={c: |r_c-r_p|>=1.0 km}`，在θ和z全范围剖面。定义每网格的无噪声模板“远距离区分裕量”：

\[
 G^{(m)}_R(p)=S^{(m)}_R(p\mid w_R)-\max_{c\in C_{wrong}} S^{(m)}_R(c\mid w_R).
\]

已存场中`H`网格探针的`S_H(p)=1`属数学必然，不能把高网格的“正确top1”称为定位成功。标记`G`负值/零值、最强错误分支距离、同峰近同分；分支只指预先固定的距离排除1km区间，不代表所有连续物理分支。

令每资源每探针 `B_R(p)=E_R(p)+max_{c∈C_wrong} E_R(c)`，则由逐候选界严格推出

\[
|G^{(H)}_R(p)-G^{(L)}_R(p)|\le B_R(p).
\]

注意：**比较时，两网格共用同一个高网格探针向量w_R**。若把低网格另当“真实信号”重新生成，就不得沿用这个上界。

对两个固定资源对照，各自按对应阵列的高网格探针构造不同向量w（不同通道数不互比样本字节），定义`Adv^{(m)}_A=G_{A_P4M}^{(m)}-G_{A0_0}^{(m)}`，`Adv_B=G_{B1_V4}^{(m)}-G_{B2_DUP}^{(m)}`。数值敏感性界为 `B_t+B_control`。若高网格Adv的绝对值严格大于该和，可标 `SIGN_ROBUST_ACROSS_THE_TWO_ARCHIVED_MESHES_IN_THIS_TEMPLATE_PROBE`并登记方向；否则`SIGN_UNRESOLVED_BY_BOUND`。另记录直接两网格是否实见符号翻转，`same sign`本身不等于被上界证明。

每一探针都独立，不能从五个模板重复测量虚报统计显著性。两张目录的同坐标探针也不是独立接收试验。必须输出`TEMPLATE_PROBE_RANKING.csv`（9×2×5=90行）与`PAIR_ADVANTAGE_NUMERIC_INTERVALS.csv`（2×2×5=20行）。**不允许使用旧门槛进行接受集合、真值保留、定位P95或客户端硬件优劣判定**。

## 6. 必须单独处理的物理／统计范围

- 同维度配对：T2-A `A_P4M` vs `A0_0`均8通道；T2-B `B1_V4` vs `B2_DUP`均12通道。不能将A与B当同资源，也不比较跨资源不同噪声权重的“工程精度”。
- `e/E/B`只是两次离散化的**条件数值误差预算**，完全不涵盖 `CHIGH`遗漏、不同传播求解器/SSP、源频率、姿态位置不确定性或阵型变形。原6 KRAKEN、3 FIELD调度与原0.2%失败完整封存。
- 3频分别赋予每频未知的公共复源谱，归一化后的投影评分不获得跨频相干收益。这里仅考察G0冻结处理器，不复现博士论文全部WSF/海底反射区4/5章。
- **单线阵有左右舷镜像**；本轮目录只使用非负`abs(theta)`代表项，不能证明镜像消歧。实际目标/平台运动、UUV速度航向、深度恢复、源信号提取均不评价。
- FIELD量化复核是原登记失败的补充原因分析，不是一种替代FIELD判据。数据不相等也不能通过更换精度定义来掩盖。

## 7. 独立审计、输出和停止规则

独立审计必须与主评分实现分离（不调用主评分函数），从读出的hi/lo存档复场直接重建至少全部433620个射影误差与e存档、完整E字段、至少全部90组模板探针在两个网格下的参考/最强远候选评分与20项配对优势。用直接投影二次型、`np.vdot`及小规模Hermitian矩阵2-范数验证解析上界，包含零场、安全数值、同相缩放、相位旋转、重复通道控制。

要求比较数值容差：存档shape_error的相对路线最大绝对差 `<=1e-10`，独立模板评分差`<=1e-10`，所有解析不等式可容许数值往返误差 `1e-10`。极小量时不能直接用`1-|inner|²`截断误报为真零；必须有稳定残差形式的替代交叉验证。冷复核不可只读回自己的输出重新计数，而应重新读封存原复场。

输出至少：

- `N1_DESIGN_FREEZE.json`, `INPUT_IDENTITY_AND_SCOPE.json`, `ANALYTIC_PROJECTOR_IDENTITY_CHECKS.csv`；
- `FIELD_SHAPE_DISTRIBUTION.csv`, `FIELD_NORM_DECILE_STRESS.csv`, `STRATIFIED_FAILURE_COUNTS.csv`, `TOP_ERROR_CANDIDATES.csv`；
- `SCORE_PERTURBATION_BOUNDS_BY_RESOURCE.csv`, `TEMPLATE_PROBE_RANKING.csv`, `PAIR_ADVANTAGE_NUMERIC_INTERVALS.csv`；
- `COLD_AUDIT.csv`, `VALIDATION.json`, `DECISION.json`, `GPT_SYNC.md`, `NEXT_RESEARCH_DECISION.md`；
- 可选两幅图（误差强度/弱场分层、远距区分裕量区间），不得画精度P95或声学定位成功率。大于10MB的中间数组保留D盘，Git只保留证据摘要、每组完整分母和哈希。

本轮决策字段必须同时包含 `HISTORICAL_G0_FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED`、`A_NOT_EVALUATED`、`B_NOT_EVALUATED`、`CHIGH_UNCERTIFIED`、`R4_percent=0`。

在这些强制字段外，另行给出三种**分析性而非原G0升级**标签：

- `ARCHIVED_TEMPLATE_NUMERIC_COMPARISON_ROBUST_IN_RECORDED_PROBES`：若某一主资源对在primary和fine目录五个探针全部达到上述优势符号的严格上界判别且优势方向一致；具体到A或B，不能合并。
- `ARCHIVED_TEMPLATE_RANKING_OR_ADVANTAGE_SENSITIVE`：若任一探针/资源实际出现最强错误对手ID变化、G符号翻转或某主对Adv符号翻转；完整逐项报告。此标签不等于物理不可行。
- `ARCHIVED_TEMPLATE_NUMERIC_IMPACT_UNRESOLVED`：未达到前两项证明时仍如实报告未决；可能与前一敏感标签并存。所有标签仅针对存档模板和冻结数值目录。

建议下一轮须根据结果明确分叉：若数值歧义确实影响关键错误分支，另行设计**有独立物理基准**的传播求解认证，不允许自动无限加密网格；若只在模板探针下排序保持稳定，则**仍必须另行认证CHIGH和实际观测可用性**才谈A/B性能；如果存档身份失败则停止并请求原始文件。任何结果**不得**把T2-A/B由NOT_EVALUATED改判通过或失败。

执行过程：A提交并推送（冻结输入、公式、探针、预算）→ **唯一一次**本轮数据诊断 → 完整独立冷复核 → B提交、推送、remote/main全SHA核对 → STOP。B阶段只增加结果/审计文件，不修改已冻结科学代码。发现新问题只保存异常并STOP，不自动做N2、KRAKEN、FIELD、接收数据模拟、优化器、垂向架构变更。同步回复只给本轮新证据，不再重贴整个三档历史。
