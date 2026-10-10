# 716项目 HLA-H2：原单HLA径向运动信息的条件价值与精度需求

版本：2026-10-10 v1.0
性质：APPLICATION_PRE_RESEARCH / CONDITIONAL_MEASUREMENT_VALUE_PILOT
不是径向速度提取成功声明，不是已授权增加实际传感器，不是五维验收。

## 0. 一次任务的目标与边界

仓库：chengcongcong222/WSL-CZ-target-localization
起始parent SHA：dd13a0a32a4bd279c03c272b9ee8b86ceea14730
新目录：results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT/
新代码前缀：hla_h2_

用户转交本文件和启动提示词后，执行授权仅覆盖：实现下述条件观测模型、设计冻结Commit A、一次固定数值试验、独立复核、结果Commit B和STOP。不自动运行声学波形/传播求解器、不改H1、不增加硬件，也不运行后续信号提取。

这一轮必须用离网格连续状态给出数值结论：

> 假如原HLA接收信号经过后续待研究的处理，能够提供少量径向净位移速率特征，那么在1200 s、小幅转向条件下，哪些水平参数有望改善？特征需要多高精度？接近/远离符号不明、速度特征零点未知时，收益是否还在？

“假如有这种特征”和“已经从原HLA提取出这种特征”完全分开。本轮只验证前者，回答是否值得投入后一条提取链。即便全部正向，最大结论也只能是CONDITIONAL_FEATURE_VALUE_SUPPORTED，不能叫原HLA新增精度已经实现。

原单HLA是当前主线。4 m垂向扩展、5 km辅助节点保持独立存档，不使用它们的观测/轨迹/指标。深度不设收紧目标；本轮为二维水平运动几何模型，没有声学传播，因此不把深度设为真值以换取水平收益。二维径向距离与实际传播/斜距的映射尚未认证，必须注明。

继承H1：有限网格公共源幅条件增量保留；M0主比较最优格点同样正确但绝对残差拒绝；M1离网格0/32，事后真值32/32可接受；线特异源变化仍失配。停止H1的旧声学全局搜索修补。本轮不是增加H1起点、换H1优化器或给H1一个真实距离邻域。

## 1. 先读的有限资料与证据分层

只读以下文件及必要直接依赖，不重新清点全部项目：

1. results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/GPT_SYNC.md
2. 同目录 DECISION.json、OFFGRID_SUMMARY.csv、OFFGRID_TRUTH_EVALUATION_ONLY.csv
3. results/R4_OBSERVABILITY_FRONTIER_RESEARCH_RESET/E1_ACOUSTIC_MEASUREMENT_MODEL_LOCK.md
4. 同目录 E1_SOURCE_PROVENANCE.json
5. results/R4_E1_G0_FREQUENCY_INFORMATION/GPT_SYNC.md
6. results/R4_RC2_OBSERVATION_SUPPORT/OBSERVATION_AND_ERROR_CONTRACT.md

文献依据：
- Xu et al., 2018, 方位和径向速度联合的浅海目标运动分析方法，DOI 10.15949/j.cnki.0371-0025.2018.03.007。仓库已锁定有条件的复声压相关→径向速度幅值→方位联合链，尚未在当前CZ实际运行。
- Li et al., 2019, 联合方位-径向速度的粒子滤波目标运动分析，DOI 10.15949/j.cnki.0371-0025.2019.04.013。是同类量测的另一估计器，不是自动增加新物理信息。
- Wang et al., 2025, 利用水声目标辐射线谱的方位-径向速度联合距离估计方法，声学学报50(5):1095–1107，DOI 10.12395/0371-0025.2024010。
  出版社正文：https://www.cpsjournals.cn/article/doi/10.12395/0371-0025.2024010
  PDF：https://www.cpsjournals.cn/data/article/sxxb/preview/pdf/10.12395/0371-0025.2024010.pdf
  本次规划读到了出版社HTML和PDF提取文本；PDF页图接口失败，未独立做完整原文公式图像锁定。因此不直接采用其解析测距公式，也不宣称完整复现该论文。其浅海、固定接收阵、较长观测等条件不能迁移为本项目性能。

下面的有限时间径向增量模型、分支枚举和连续批处理估计器，是本任务显式定义的改编试验，不是三篇论文之一的完整复现。“方位+径向速度”及线性最小二乘本身不作为新颖性主张。

旧E1-G0研究的是双节点等条件下未知源频率/参考的局部信息，不是本轮单HLA条件径向增量的经验性能。但不能利用本轮假设的速度零点绕过E1的未知源频率问题。只有将本轮定位为特征精度需求试验，这种分层才合法。

## 2. 冻结几何与观测时刻

### 2.1 水平状态和平台

h=[r_km, theta_deg, v_mps, psi_deg]。
物理域：r∈[45,60] km，theta∈[-5,5] deg，v∈[1,3] m/s，psi∈[-15,15] deg。
方位从+x轴逆时针计，使用atan2(y,x)。与文献北向坐标的差别必须写明。

目标水平位置：
p_T(t)=1000r[cos(theta),sin(theta)]+v t[cos(psi),sin(psi)]。

平台：t<=600 s以2 m/s沿+x直航；600 s后以2 m/s沿15 deg方向航行。
位置沿用H1的确定性解析式。导航精确、方位系统偏差为0，均为条件假设。

方位：t=0,10,...,1200 s，共121个有符号方位。
继承sigma_beta=0.1 deg，运算以rad进行。本轮不认证真实HLA左右舷解模糊。

### 2.2 只给六个径向特征，不给121个独立速度观测

窗口端点t_j=0,200,400,600,800,1000,1200 s，共6个不重叠200 s窗口；600 s转向位于窗口边界。
rho(t;h)=||p_T(t;h)-p_O(t)||，单位m。
定义
qbar_j(h)=[rho(t_j;h)-rho(t_{j-1};h)]/200，单位m/s。

这是“200 s内净径向位移除以时间”的特征，不是径向速度瞬时值，也不是平均绝对速度。必须通过完整端点几何计算，禁止替换成中点瞬时q后仍声称模型精确。若窗口内q变号，|qbar|与平均|q|不同：保留并单独标记，不能删除该场景。

窗口用于假设后续相关处理的更新尺度，200 s不是已证明的源相干时间。当前不存在从相关谱峰精确导出这个窗口特征的实现。实际提取链若给出另一个时间加权量，应在以后独立定义，不能借用本轮结果冒充。

## 3. 四种估计模型和对照

### B：原条件方位基线

y_beta,k=wrap(beta_k(h_true)+sigma_beta*e_beta,k)。
只使用方位和平台，不使用任何径向量。

### S：有符号、零点已校准的条件上限

y_q,j=qbar_j(h_true)+sigma_q*e_q,j。
估计器看到带噪y_q，不看到真q。
标签：SIGNED_CALIBRATED_FEATURE_UPPER_CONTROL。
本组假定获得接近/远离方向及准确速度零点，属于比主组更强的特征条件。

### U：只有幅值、零点已校准——主研究组

u_j=abs(qbar_j(h_true)+sigma_q*e_q,j)。

估计器只看到非负u，不看到折叠前的符号。不能使用truth、谱条纹未认证符号或假想辅助阵获取符号。
必须强调：它是有符号高斯误差观测取绝对值的folded-normal条件族；并不是abs(qbar)+独立高斯，也不是直接生成无误差abs(qbar)。
标签：UNSIGNED_CALIBRATED_FEATURE_CONDITIONAL。

### UB：只有幅值，且存在未知常量零点偏置——关键边界组

u_j=abs(qbar_j(h_true)+b_true+sigma_q*e_q,j)。
生成端固定b_true=+0.10 m/s；估计器只知道设计域b∈[-0.20,+0.20] m/s，不读取b_true。这个范围是预研假设，不是实测仪器/声源指标。

在同一组偏置数据上运行：
- U-mismatch：错误地假定b=0，用于看清忽略偏置的代价；
- UB-matched：将同一个b跨全部6个窗口联合估计/剖面。

b是一个便于解释的“径向特征零点误差”代理，不能宣称已精确覆盖声源漂移、SSP误差、频率参考及相干性。允许后续真实提取显示这种模型仍不成立。

独立条件：方位噪声与径向特征噪声独立、不同窗口径向误差独立；只是本轮的特征统计合同，不是同一录音派生统计量已经独立。不能将相关多频输出当多份独立q观测。

## 4. 固定样本矩阵：主要是离网格，不重复H1真值在格点的问题

### 4.1 几何

使用H1 OFFGRID_TRUTH_EVALUATION_ONLY.csv中全部8个连续水平状态，不挑成功状态、不改变坐标。本批已见状态只叫DEVELOPMENT_GEOMETRY_PANEL；不叫新的独立确认集。
深度字段仅保留为来源元数据，不进入本轮生成器/估计器几何。

### 4.2 噪声与配置

每几何32个新特征噪声实现，先冻结种子，后生成。
sigma_q∈{0.02,0.05,0.10,0.20} m/s。
主比较：sigma_q=0.05 m/s、b=0、组U。

- B：8×32=256配置，仅生成并估计一次。
- S和U：8×32×4×2=2048配置。
- 偏置边界，仅sigma_q=0.05：8×32×2(U-mismatch与UB-matched)=512配置。
合计2816个方法配置，不把B复制到每个误差档重复计数。

每个(geometry,replicate)使用一组长度121的e_beta和长度6的e_q，在全部方法、sigma_q和偏置条件之间配对复用；共有256组独立基础创新，不是2816个独立样本。
推荐固定seed主值2026101002，用SeedSequence([seed,g,rep,stream])区分beta和q；预检seed2026101012。登记PRNG版本和生成规则，存全部实际观测数组。

本轮允许新增特征级MC；禁止声学时域录音、传播求解器、真实源频率注入、H1声级与q伪独立融合。

## 5. 接受规则：让未知符号进入模型，不用真值挑一支

使用相同方位兼容条件：
J_beta(h)=sum[wrap(beta(h)-y_beta)/sigma_beta]^2 <= T_beta。
T_beta=chi2.ppf(0.975,121)。

S的径向条件：
J_S(h)=sum[(qbar(h)-y_q)/sigma_q]^2 <= T_q。

U的径向条件：
J_U(h)=sum[(abs(qbar(h))-u)/sigma_q]^2 <= T_q。

UB的径向条件：
J_UB(h)=min_{b∈[-.2,.2]} sum[(abs(qbar(h)+b)-u)/sigma_q]^2 <= T_q。

T_q=chi2.ppf(0.975,6)。本轮没有因估计h/b而减自由度。
联合接受为J_beta<=T_beta AND J_radial<=T_q，而不是事后改为Jmin+任意门槛。

概率含义必须从“固定真实状态处的噪声事件”推导：
- signed合同中两个标准化噪声平方和各以概率0.975不超过对应阈值。
- ||a|-|b||<=|a-b|保证folded残差在真值处不大于原始高斯误差；b域包含真实b时剖面不会增加真值代价。
- 用union bound，正确模型的理想兼容集合具有至少0.95的条件真值包含概率；不需要额外减估计参数自由度。

这一概率属于理想的连续兼容集合，不属于有限起点导出的候选集。若真值兼容但求解器没有找到，必须单独标SEARCH_FAILURE。B单独是0.975条件集合，不冒称所有方法有完全相同名义覆盖。

U的残差不是标准chi-square精确分布，而是上述保守事件推导；不得称其为folded-normal边缘最大似然或精确LR置信区间。使用分支最小二乘只是一个已说明的计算策略。

## 6. 一个有结构的连续估计器，不重开声学全局搜索

### 6.1 显式枚举六个folded观测的符号

U/UB完整枚举s∈{-1,+1}^6，共64个观测展开分支。每支使用q_obs=s*u。
枚举的是折叠前有噪声观测的符号，不是强行规定真实q符号。接近零时噪声符号可以与真实q相反，不能以“q应单调”删除这些分支。

固定分支s内的光滑目标：
J_s(h,b)=J_beta(h)+sum[(qbar(h)+b-s*u)/sigma_q]^2。
U取b=0；UB联合估计b。S只用一个已观测有符号分支。

全局CV轨迹与已知平台转向跨六窗一致，不按窗口分别输出六个独立目标状态。

### 6.2 用径向增量积分形成观测侧初始化

对端点t_j，在现有121个方位里取同一端点观测，e_j=[cos(y_beta(t_j)),sin(y_beta(t_j))]。
定义d_0=0，d_j=200*sum_{i<=j} q_obs,i。
由rho_j≈rho_0+d_j-b*t_j与p_T=p_O+rho e，得到初始化用线性关系：

p_0+v*t_j-rho_0*e_j+b*t_j*e_j = p_O(t_j)+d_j*e_j。

用缩放最小二乘求(p_0x,p_0y,vx,vy,rho_0[,b])；它只给起点。
注意角度有噪声，累积d的噪声相关，rho_0与||p_0-p_O0||的非线性关系还没强制，因此该线性解不能直接作为无偏估计或精度证明。最终评分使用原始方位+六个未累积径向观测，避免重复计数。

把线性p0/v转回h，裁剪到物理域只用于初始化；非正rho0、秩不足、裁剪、条件数均记录，不删配置。秩不足时用固定缩放SVD最小范数解并记录。

每个分支只用两个预先规定的起点：线性几何起点、物理域中心[52.5,0,2,0]，UB中心b=0；若两者相同仍保留来源标记，可复用运算但不能重复统计独立成功。

B用原纯方位Cartesian线性初始化加固定16个距离分层中心，各起点不使用q；登记16个范围中心和相同theta/v/psi中心，允许重复结果但不以truth择优。

### 6.3 最终连续优化与导数

使用固定scipy least_squares TRF、linear loss；h归一到[0,1]^4，UB多一个归一b。
max_nfev=100，ftol=xtol=gtol=1e-9，不在结果后提高预算。
实现解析Jacobian：
- 目标位置对r/theta/v/psi的导数；
- beta的atan2链式导数；
- qbar端点距离差导数：
  d qbar_j/d h=[e(t_j)^T d p_T(t_j)/dh-e(t_{j-1})^T d p_T(t_{j-1})/dh]/200。
角度参数的deg/rad因子、r的km/m因子必须显式测试。

每个终点重新用独立几何计算J_beta、原分支J_s、folded J_U；UB另对固定h的b做完整一维分段二次最小化核验：断点为b=-qbar_j，加b域端点；每段符号固定，取该段内二次函数最小点与边界。保留原终点b/原分支代价，和重新剖面的b/代价，不把两种记录混为同一个审计对象。

接受判断按第5节。排序先在已接受候选里选最小J_beta+J_radial；没有接受候选则输出ABSTAIN，并另存所有失败终点及其最小代价作为诊断。不得把拒绝点的误差放进“已接受精度”主表。

64个分支全部保留，不能先按truth、cost或符号只保留若干分支。最终同一h/b的近重复终点可以去重用于绘图，判定和来源表不丢失。去重阈值预注册[r .005km,theta .001deg,v .005m/s,psi .05deg,b .005m/s]，不能据此称物理分支可分辨尺度。

完整符号枚举不等于每个分支的连续全局最优已认证；两起点也不等于排尽局部解。未找到接受解不直接推断集合为空。

## 7. 必须交付的数值指标

### 7.1 四参数误差主表

每个几何、方法、sigma_q/偏置条件分别报告：
- 已接受输出数、拒判数、求解器异常数；
- 真实状态是否属于理想兼容集（只在估计完成后做）；
- 真值兼容但没有任何可接受搜索结果的次数；
- 距离/速度相对误差，方位/航向最小圆周绝对误差；
- 已接受输出的中位/P90/P95/最大误差，并显式分母；
- 将拒判/数值失败记INF后的无条件中位/P90/P95，避免只统计成功；
- 未接受最优点误差另外放SEARCH_DIAGNOSTICS，不冒充有效定位精度。

分位数用nearest-rank，样本32时P95是31st顺序统计量，明确探索性；附输出率Wilson区间。跨8几何汇总用每几何同权，不用漂亮pool代替最坏几何。

### 7.2 多分支结果

对所有接受终点，保存符号展开分支、h、b、代价、是否边界、原始起点等。
报告各参数的终点包络、分离的候选簇及最远错误分支（后者仅评估端）。这些命名为FINITE_TERMINAL_ENVELOPE，不能叫95%置信区间或完整连续外包络。

本轮不做H1式完整四维粗网格排除，也不从单点结果报告“集合宽度0=误差0”。成功信号主要根据离网格估计与失败率，而不是候选数压到1。

### 7.3 精度需求图

固定四档sigma_q逐档展示r/theta/v/psi误差；报告在已测档位内，哪些档位满足参考线：距离P95 10%、速度P95 10%、航向P95 5deg、方位P95 1deg。参考线不是本轮必须通过的Gate，逐参数分别标，不要求四项同时通过。

只在实际测过的档位上说“这档达到/未达到”。不在四点间插值出高精度的临界sigma，不称最优阈值。P95不稳定或INF原样保留。

如果U改善速度但距离仍宽/误差大，登记VELOCITY_GAIN_WITH_RANGE_LIMIT，不判整轮失败。如果只有S改善，登记SIGN_INFORMATION_REQUIRED，不把S性能转给U。如果UB损失明显，登记RADIAL_REFERENCE_CONDITION_CRITICAL。

### 7.4 物理接收条件对照（不执行接收仿真）

写一页MEASUREMENT_TO_SIGNAL_GAP.md：
1. 条件qbar特征、论文相关法的局部速度/幅值、真实单HLA接收数据是不同层次；本轮未建立映射。
2. p(f,t)=S(f,t)H(f,t)的跨时乘积保留S(f,t)S*(f,t+tau)；公共源幅能消除不意味着跨时源相位可消除。不能用裸H代替接收复声压。
3. 未知发射频率不能用真实f0解调或直接换算为q。本轮b控制也不等于完整源相位/频率误差模型。
4. 展示单分量控制：当p(t)=A exp{i(omega-k q)t}，联合改变q和omega可保持p不变。它是简单模型的源/运动混淆示例，不是整个多模态方法族No-Go。
5. 200 s、201/235/283 Hz下可列c_ref/(f*T)的普通谱栅格参考尺度，c_ref=1500 m/s仅示意；它不是实际sigma_q，不是精度下界，零填充不能自动创造分辨力。源相干、相速度、窗内变化、弱场、共享频带相关仍未验证。
6. 不能仅因某档假设sigma_q带来漂亮定位结果，就宣布该档可由真实原HLA获得。

## 8. 前置检查与独立审计

Commit A之前只做有限解析/代码控制，不跑8几何主试验，不看主科学排名后改设计：
- 单位、坐标、端点增量等于积分q的解析/高精度几何对照；
- 直线共速、转向、近零q、窗口内变号、方位圆周跨缝的固定数学控制；
- 无噪声观测下线性初始化恒等与最终正确自匹配；
- 未知符号的min_{s} squared residual等于(abs(q)-u)^2；
- 不强制真实q符号与展开噪声符号一致；
- UB分段二次b剖面与密集独立b控制/解析候选一致；
- 解析Jacobian与中心差分在固定物理内点/边界附近对照，最大归一化差<=1e-5；
- 最终直接物理残差独立复算差<=1e-7相对、单位正常；
- 非空、空集、拒判、INF分位数、重复候选的汇总单元测试；
- 构造器API不读取truth/noise创新/隐藏符号/b_true；评估在全部构造完成后读取。

设计文件记录actual input hashes、派生方位/径向阈值、分支数、起点、版本、时间/内存预算、各文件路径和种子。主执行只读取观测文件与模型合同，真值在单独生成模块和单独评估模块中使用。

Commit B前：
- 冷复算每条导出终点的几何、原分支代价、folded/profiled代价、接受状态；
- 全2816配置的来源/数量/方法/误差档对账，所有失败不丢；
- 复算全部误差、分位数、输出率和Gate；
- 核对所有64符号分支均已运行，异常/预算终止明确；
- 固定8个配置（各几何rep0，U主档）用另一Jacobian实现复核所存终点的梯度/代价，不另找更好解；
- 对每个主配置事后计算true h,b的兼容性，只用于区分SEARCH_FAILURE与统计拒绝。

不把PASS计数当成创新或精度。前置实现失误可在Commit A之前修正，但不可读取主矩阵后调参。Commit A之后若需改变估计器、观测、阈值、起点或合同，停止并交付原结果；不得在本阶段补跑。

## 9. 描述性研究信号和停止分类

主判断只看U、sigma_q=0.05、b=0；S为更强特征控制，不可救主组。

正向信号建议：至少5/8几何同时满足：
- 接受输出率>=90%；
- 无条件中位误差在r、v、psi至少一项相对同观测B改善>=30%。
  B中位为0时该项不能用比值；为INF时只报告从无输出到输出，不靠INF比例自动通过；需要该项P90有限且优于B，或者另一个有限基线项通过。

各参数分开，深度不参与。该门槛是预注册探索性信号，不是0.9真实成功率统计证明，也不是甲方指标。

分类允许并列：
- HLA_H2_UNSIGNED_RADIAL_FEATURE_VALUE_SUPPORTED_CONDITIONALLY；
- HLA_H2_SIGN_INFORMATION_REQUIRED；
- HLA_H2_RADIAL_REFERENCE_CONDITION_CRITICAL；
- HLA_H2_VELOCITY_GAIN_WITH_RANGE_LIMIT；
- HLA_H2_BENEFIT_LIMITED_AT_TESTED_PRECISIONS；
- HLA_H2_CONTINUOUS_SEARCH_LIMITED；
- HLA_H2_IMPLEMENTATION_OR_INPUT_INVALID；
- HLA_H2_PARTIAL_EXECUTION。

真实接收提取始终NOT_EVALUATED，COLD_START_ACOUSTIC_CAPABILITY_NOT_ESTABLISHED。所有正向结论必须含INJECTED_FEATURE_ONLY和签名/零点/噪声条件。

本轮所有等级都在Commit B之后STOP。不根据S/U/UB输赢换论文方法、不追加更精细sigma、不延长观测、不调用H1声学处理做联合性能，不转向增强架构。

## 10. 预算与产物

单进程，BLAS最多4线程；内存上限8 GiB；全部主计算+必要冷复核总wall上限14400 s；仓库交付<=512 MiB；不推大无关缓存。预算是停止上限，不是运行时间承诺。
按照geometry、rep、sigma、method的固定顺序执行，异常写checkpoint。不能因为某几何较难而跳过或增加专用起点。预算到达保留PARTIAL_EXECUTION，未做项为NOT_EVALUATED。
本轮KRAKEN/FIELD/BELLHOP=0、音频/原始复谱生成=0；允许新特征级Monte Carlo和纯几何前向计算，资源表不要把它们写成“没有任何新数值”。

必要产物（可合并同内容，避免文档膨胀）：
- DESIGN_FREEZE.json，METHOD_AND_OBSERVATION_CONTRACT.md，PREFLIGHT.json
- OBSERVATIONS.npz（估计器可见数据，不含真值/隐符号），生成端truth/创新分开存档
- ALL_TERMINALS.csv.gz，ESTIMATES.csv，CONFIGURATION_RESULTS.csv
- FOUR_PARAMETER_RESULTS.md，PRECISION_REQUIREMENT_TABLE.csv
- SIGN_AND_BIAS_BOUNDARY.csv，SEARCH_FAILURE_DIAGNOSTICS.csv
- MEASUREMENT_TO_SIGNAL_GAP.md
- VALIDATION.json，OUTPUT_MANIFEST.json，DECISION.json
- GPT_SYNC.md（建议1000–1600中文字，重点为数字与科学作用域）
- 三张图：四参数误差-径向精度关系；signed/unsigned差异；未知偏置影响。失败、拒判及未做条件必须在图中可见，不裁掉最差场景。

GPT_SYNC必须回答：
1. 原单HLA只是假定增加了什么可提取特征？是否实际提取（必须写没有）？
2. 四参数在U主档下各自中位/P95/失败率怎样？
3. 哪些sigma档有价值、哪些没有？
4. S与U的差别是否说明符号是关键条件？
5. 共同零点未知后收益保留多少？
6. 输出集合/包络到底是否连续全域认证（必须写不是）？
7. 建议实际信号提取优先还是另换机制？只建议，不启动。

## 11. Commit A / Commit B纪律

- 首先HEAD与remote/main必须等于上述parent；不相等时停止并报告差异，不自动rebase/合并。
- 不导入会在顶层运行、建目录或写历史文件的旧脚本。
- Commit A提交实现、来源作用域、预检和冻结，不含主实验结果；推送核对后只执行一次。
- Commit B提交全部结果和审计；历史H1/E1/E2/H3/R3数值保持原字节，只允许在管理台账新增一条有准确作用域的记录。
- 原R4=0%不回填。本轮条件特征试验不获得原硬Gate进度。
- 最终回报完整parent/设计/执行SHA、remote核对、数据入口和STOP。
