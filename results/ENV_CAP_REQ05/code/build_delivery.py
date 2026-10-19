import pathlib,json,csv,hashlib,time,subprocess
R=pathlib.Path(r'D:\XiaomiMiMoProjects\WSL系统汇聚区目标定位');O=R/'results/ENV_CAP_REQ05';W=pathlib.Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_REQ05')
def put(n,s):(O/n).write_text(s.strip()+'\n',encoding='utf-8')
def table(n,rows):
 with (O/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
obs=[
 dict(id='O1',name='单频声级时间包络',mathematical='受约束源/增益下可能有独立距离方向；逐时自由时严格为零',archive='CZO04存档I类200Hz模板及有限目录；仅一几何及镜像保留完整条件证据',extractable='NOT_ESTABLISHED；没有UUV功率/AGC/统计相关性验证',conditional='需要源级/增益时间约束、功率统计与环境导数证据',hardware='原HLA可采功率；处理不能创造源级约束',evidence='CZO04;NUM03',next='STOP_SAME_MODEL_SCAN'),
 dict(id='O2',name='同时不同接收深度的相对强度',mathematical='共同源级可消去；差分增益、r-z与环境耦合仍在',archive='新CZ仅有条件期望能量背景；H3复响应局部增益不能视为功率比增益',extractable='NOT_ESTABLISHED；无真实同步垂向功率比提取证据',conditional='共同激励、差分增益/AGC及噪声底受约束、深度/姿态和环境模型',hardware='名义同深HLA缺垂向采样；需实际已知非共深阵元或额外同步垂向接收',evidence='H3G0;H3G1;H3G2E;REC01;CZO04',next='CONDITIONAL_CANDIDATE_DEPTH_NOT_EXECUTABLE_YET'),
 dict(id='O3',name='同频互谱/相对相位',mathematical='公共相位可消去；未知逐元复响应及平面波控制限制独立r-z信息',archive='H3名义局部正信息；SP1全距离集合未收紧；G2E压缩提取不可靠',extractable='SYNTHETIC_STATIC_CONTROL_ONLY；真实/动态CZ提取NOT_ESTABLISHED',conditional='同步逐元复数据、标定、噪声协方差、认证前向与完整候选',hardware='原8通道可具采集能力；只交方位/功率则无法恢复复信息；处理不能增加实际孔径',evidence='H3G0;H3G1;H3G2E;SP1;FUNC2',next='DEFER_NO_REOPEN_SP1_OR_T2'),
 dict(id='O4',name='多途相对时延',mathematical='公共发射时刻可消去；路径独立导数可能增加r-z信息',archive='NUM03几何及文献条件；无本项目已分离宽带到达与路径能量证书',extractable='NOT_ESTABLISHED；单CW周期模糊及路径混叠保留',conditional='有效带宽/时间标记、相关峰与路径身份、环境/底质约束',hardware='是否需新硬件取决于现有带宽和数据；软件不能生成不存在的信号带宽',evidence='H3_COVERAGE;NUM03;REC01',next='DEFER_NO_BROADBAND_RECEIVER_EVIDENCE'),
 dict(id='O5',name='运动频率与相干孔径',mathematical='相干/频率变化可承载运动和模式信息；未知频率/漂移/相位规范可吞噬',archive='Yang浅海机制与MMAR部分链；E1增量很小；RX1Z相对运动映射未支持',extractable='NOT_ESTABLISHED_IN_CZ；不能将生成k/速度/源频当输入',conditional='受控源/时频参考、径向采样/相干窗、模式身份与可靠复场',hardware='记录本身可能用原HLA；稳定参考、足够径向孔径无法由重处理保证',evidence='PAPER_YANG;PAPER_ERRATUM;PAPER_MMAR;METHOD_ID;E1;RX1Z;FUNC2',next='DEFER_NO_RENAMED_DOPPLER_OR_MMAR_PERFORMANCE'),
 dict(id='O6',name='辅助节点独立方位',mathematical='非共线几何与持续运动可增加水平信息；此模型无z列',archive='AUXD1配对档案和T3W1有界偏置下水平Demo',extractable='HYPOTHETICAL_BEARINGS_ONLY；真实有符号方位/同步/关联UNKNOWN',conditional='双节点观测与报告轨迹、偏置依据、无符号多支、导航协方差与关联',hardware='需第二测向节点及部署/通信/时间资源；不要求跨节点载波相干',evidence='RC2;SUPPORT;AUXD1;T3W1;CLIENT',next='CONDITIONAL_CANDIDATE_HORIZONTAL_CHAIN_ONLY')]
table('OBSERVATION_EVIDENCE_LAYERS.csv',obs)
state=[
('O1',['CZO04有限离散条件增量，非一般锚','仅经运动和传播间接耦合','干扰需剖面；本轮不评深度精度','通过r(t)耦合，漂移可吸收','通过r(t)耦合，镜像保留']),
('O2',['可能与z独立，也可能耦合；NOT_EVALUATED','非单独方位证据，依赖接收几何','假设性差分深度敏感；功率比尚未验证','仅多时运动联合；未验证','仅多时运动联合；镜像未解']),
('O3',['平面波零；复杂波导可能；SP1未收紧','原HLA方向结构但镜像保留','H3复响应局部条件增益；非精度','动态响应条件性；非同窗直接量','动态响应条件性；镜像未解']),
('O4',['路径差导数条件性；路径标签未知','依赖阵列/路径方向，非自动独立','路径差导数条件性；环境耦合','需跨时或其它观测联合','需跨时或其它观测联合']),
('O5',['offset与模式相位耦合，未建立一般锚','依赖方位/姿态或既有HLA','Yang/MMAR模式依赖；第一CZ NOT_EVALUATED','源频/漂移/时钟干扰严重；E1小增量','与LOS运动联合；非独立直接观测']),
('O6',['条件水平Demo；依赖导航和偏置','有向方位假设；无符号分支保留','STRUCTURAL_ZERO_IN_HORIZONTAL_BEARING_MODEL','条件Demo；偏置先验影响弱方向','条件Demo；关联和左右舷未实测'])]
table('FIVE_STATE_INFORMATION_MATRIX.csv',[dict(observation=o,state=s,status=v,layer='SCOPED_EVIDENCE_NOT_COMMON_PRECISION_SCALE') for o,vals in state for s,v in zip(['r','theta','z','v','psi'],vals)])
resources=[
('R01','r','O1','源级+增益不可逐时完全自由','MATHEMATICAL_NECESSARY','常数或受约束漂移是假设；CZO04仅有限目录条件成立','功率/AGC时序、噪声相关、源稳定资料','不能反推出设备最低0.25dB'),
('R02','r,z','O2','差分响应至少一方向不落入增益/状态/环境未知量空间','MATHEMATICAL_NECESSARY','同步不同深度；稳定差分增益是假设','逐元深度、校准及同窗功率统计','不存在已验证最少阵元数/垂向孔径'),
('R03','r,z','O2','不能用共同源消除替代差分增益标定','MATHEMATICAL_NECESSARY','恒定增益多时可能保留变化；静态比例可能全被吸收','跨通道AGC、温漂、幅频响应时序','H3复场2.34-4.47倍不是功率比收益'),
('R04','theta,r,z','O3','公共源剖面后需保留非平面波或几何差异','MATHEMATICAL_NECESSARY','原8元同步复数据已有资源潜力；真实接口未知','逐元复基带、共钟、相位响应、阵形','更换处理器不能增加14m真实孔径'),
('R05','r,z','O3','增益/模式幅度不能全自由；前向和数值需独立可信','MATHEMATICAL_NECESSARY','H3名义局部；T2停止；SP1集合未收紧','环境不确定性与独立复场参考','不再开发旧数值修补链'),
('R06','r,z','O4','可分辨相对延迟和可解释路径身份','MATHEMATICAL_NECESSARY','宽带/多频提供延迟尺度属假设','源自相关、可用带宽、同步原始记录、路径证据','1/B不是通用最小带宽或测距精度'),
('R07','v,psi,z','O5','源频/源相位/时频参考需受约束；径向采样可确定','MATHEMATICAL_NECESSARY','已知源频或其它径向参考均为额外信息','频率稳定资料、时钟、相干时间和运动资料','不能将H2注入径向特征当接收输出'),
('R08','z','O5','模态身份、复幅度、深度函数及offset可辨','MATHEMATICAL_NECESSARY','浅海正控制与部分MMAR链；新CZ未建立','认证复场或独立实测模式参考','2pi/L只是Fourier尺度，不是普适资源下界'),
('R09','r,theta,v,psi','O6','独立非共线观测与受控偏置/导航、关联','MATHEMATICAL_NECESSARY','AUXD1/T3W1条件仿真，不是实装','同目标双节点方位、报告轨迹、侧向分支、偏置校准','节点数量2不是所有方法最低数量'),
('R10','r,v,psi','O6','宽基线收益需接收可用性与相对布放证据','HYPOTHETICAL_CONFIGURATION','T3W1 B3=(3,±15)km条件结果；不继续扫描','部署范围/成本、持续探测及通信记录','几何45-60km不证明第一CZ可检测'),
('R11','z','O6','水平方位不提供z观测，需另一深度敏感机制','STRUCTURAL_LIMIT','该水平观测模型严格z零列','如需深度，独立垂向或传播观测合同','不能填写零深度误差或五维通过'),
('R12','all','all','共享噪声与共同源导致派生量相关','MATHEMATICAL_NECESSARY','先保留完整联合数据/协方差','同步窗口、噪声样本、丢包和处理元数据','不得把同一录音派生特征Fisher直接相加'),
('R13','r,z','all_acoustic','首CZ存在与传播可提取性分别认证','EVIDENCE_BOUNDARY','NUM03同一WOA族名义几何支持','可追溯SSP时空范围、海深/底质及不确定性','不能以9单元宣称9独立海域验证')]
table('RESOURCE_NECESSARY_CONDITIONS.csv',[dict(id=i,target=p,observation=o,condition=c,evidence_kind=k,current_config=v,needed_data=d,forbidden_inference=f) for i,p,o,c,k,v,d,f in resources])
questions=[
('D01','阵列','逐元数量与三维坐标/坐标轴/单位；名义与实测阵形','坐标CSV、机械图、误差说明','识别实际水平/垂向孔径；不能用处理制造几何','O2;O3;O6'),
('D02','姿态','拖阵姿态、倾角、弯曲、接收深度及其时间同步记录','姿态/深度时间序列、校准和协方差','剖面接收位姿混淆；不预设4m垂向已知','O2;O3'),
('D03','通道','幅相频响、通道间差异、温漂和稳定时长','校准扫频与日期、误差及相关性','约束差分gain，而非假定公共源消除全部未知','O2;O3'),
('D04','AGC','AGC作用位置、逐元还是共同、增益时序/限幅记录','原始增益档位、线性标定、饱和/量化标记','任意未知差分AGC会消灭相对功率信息','O1;O2;O3'),
('D05','数据接口','是否保留逐元同步原始时域/复基带，或仅波束/方位','匿名样例、采样率、位宽、数据字典','复相位丢失后无法由方位/功率恢复','O3;O4;O5'),
('D06','互谱','互谱矩阵与功率的计算窗、窗函数、重叠、平均和噪声扣除','处理代码/配置、完整CSD含非对角项','核对共享参考、有限样本偏差和时间相关','O1;O2;O3'),
('D07','信号','有效频带、线谱占有、宽带特征、目标关联与持续性','频谱/时频图及相应原始片段；含无信号负例','传播计算频点不等于可用UUV谱线','O1;O2;O3;O4;O5'),
('D08','源稳定','是否有源频先验、频漂/幅度/相位稳定性资料及来源','测量方法、统计区间；无资料明确UNKNOWN','不向接收器隐式提供真实源频/相位','O1;O5'),
('D09','噪声','背景与目标窗的空间/时间/频率协方差、低能量和缺测','同步噪声样例、状态标签、噪声底定义','SNR和有效样本数须可解释；不只给总RMS','all'),
('D10','阵内时钟','采样共钟、通道时延、时钟偏差与漂移','时钟架构、校准记录及稳定性','阵内互谱不同于跨节点载波相干','O3;O4;O5'),
('D11','导航','主辅报告轨迹、参考系、更新时间及绝对/相对误差','带时间戳轨迹、系统偏差与跨节点协方差','CEP不直接替代每轴1sigma；禁止用真实位置初始化','O6'),
('D12','测向','方位产品是否有符号；左右/前后歧义怎样输出','同记录多分支、置信/拒判、旁瓣和误差证据','不以真实方向或镜像标签完成消歧','O3;O6'),
('D13','测向偏置','随机误差与固定航向/安装偏差能否分离，有何独立校准','完整误差预算、校准前后记录、偏置先验来源','soft_l1或约束边界不代表系统偏置被消除','O6'),
('D14','同步关联','跨节点时间偏差/漂移、同目标关联、并发目标处理','时间基准、关联标签依据、错误/丢失记录','不把正确同步和关联作为无成本事实','O4;O6'),
('D15','环境','SSP原样本/测站时空/测深与插值方法及误差','CTD/气候资料元数据、空间时间变异范围','一个气候族不提供跨海域鲁棒性','O1;O2;O3;O4;O5'),
('D16','边界','海底地形、底质、吸收和边界来源及不确定性','测绘、地声参数来源；无资料明确未知','射线几何支持不认证复场/多途能量','O2;O3;O4;O5'),
('D17','垂向资源','可否同步取得不同接收深度；现有设备/拟增元/布放限制','可用资源清单、姿态测量、功耗/通道/风险/成本','先确认条件，不要求先购买四元或验收米级精度','O2'),
('D18','辅助资源','可用辅助测向节点、相对布放/续航/通信/接收覆盖','能力范围和成本、数据共享方式、缺测记录','不把15km级基线和0.05度视为实际规格','O6'),
('D19','时间与机动','可用观测时长、平台机动约束、时窗平稳性','连续航迹和信号段；机动/断续标记','实际径向孔径不等于平台路程或采样数','O1;O3;O5;O6')]
table('PARTNER_DATA_REQUEST.csv',[dict(id=i,category=c,question=q,requested_format=f,scientific_purpose=p,observations=o,current_status='UNKNOWN_OR_NOT_PROVIDED',numeric_acceptance_threshold='NOT_SET') for i,c,q,f,p,o in questions])
put('PARTNER_TECHNICAL_REQUEST_716.md',r'''
# 716所及合作团队：第一会聚区观测资料确认稿

当前处于项目申请预研阶段。我们希望先确认可使用的观测和资源，再确定值得验证的方法与精度范围；本表不构成设备验收或采购要求。缺少现成设备参数是正常的，可分别标记“已有数据”“设计规格”“拟研究能力”“暂无资料”，并注明来源与不确定性。本文仅形成待交流材料，本轮未发送给任何合作方。

已建立的基础是名义第一CZ折返几何、若干条件水平定位Demo及局部深度信息。实际UUV原始接收提取、一般距离锚和五维联合精度均未建立。希望优先提供以下三类材料：

1. 一份同步多通道接收样例及完整处理元数据：逐元坐标/深度/姿态，时域或复基带，频谱、互谱、AGC、校准和噪声窗。只有方位产品时请明确说明，避免假定可以恢复已经丢弃的幅相信息。
2. 一份方位—报告导航—时间戳的配对样例：若有双节点则提供两侧同目标记录，包含侧向歧义、固定偏差、失配、拒判和缺测。没有双节点数据时可先提供资源能力说明，无需为本轮专门采集或生成。
3. 与样例对应的环境和资源依据：SSP原样本及时间空间来源，海深/底质与误差；同步垂向通道、辅助节点、部署、通信与能耗的可用范围。

详细19项见PARTNER_DATA_REQUEST.csv，均未设数值验收线。历史0.05度、25m导航、4m垂向以及15km级基线均为特定仿真条件，不作为本次要求。请不要将波束宽度、单快照测向误差、航迹平均误差、随机误差和固定安装偏差混成一个数字；时间、节点和通道相关性应单独说明。

收到资料后，优先检验两个独立问题：同步不同深度的相对强度是否保留源级消除后的深度/距离敏感性；辅助节点方位增益能否在实际偏置、左右舷和关联条件下成立。这两项目前是有条件候选，不是本轮启动的新试验。没有资料时可以继续形成明确标注的假设性研究设计，但不能升级为真实设备性能。
''')
recommendations=[
 dict(id='C_DEPTH',observation='O2',name='同步不同接收深度相对强度：差分增益条件下的机制验证',objective='z/r独立敏感性，不预设测深精度',priority='PARALLEL_CONDITIONAL_NOT_RANKED',reason='直接回应CZO04公共源级退化；不要求源绝对相位；H3说明垂向几何值得保留但不证明功率比有效',required=['D01','D02','D03','D04','D06','D07','D09','D15','D16','D17'],go_requirements=['源在同步深度通道间可建共同激励模型','差分增益和噪声背景有可登记约束','有可解释期望功率/观测对应关系且r-z/e同时作为未知','前向可信度或独立数据参考足以分辨候选变化与误差'],no_go=['逐时差分AGC/增益完全未知','仅同深单波束数据且没有新增垂向资源','将I类平均能量直接冒充窄带CW接收功率','仅oracle深度/环境或局部Fisher证明全域可测深'],new_hardware='可能需同步不同深度接收；最少阵元/孔径NOT_ESTABLISHED',status='PROPOSED_NOT_FROZEN_NUMERICAL_EXECUTION_NOT_AUTHORIZED'),
 dict(id='C_HORIZONTAL',observation='O6',name='辅助节点方位链：偏置、同步、目标关联和侧向多分支验证',objective='r/theta/v/psi条件增量能否落到真实或明确非oracle观测链；z不评',priority='PARALLEL_CONDITIONAL_NOT_RANKED',reason='AUXD1和T3W1已有配对条件Demo；下一不确定性在测量链与偏置先验而非继续扩大基线',required=['D01','D07','D09','D11','D12','D13','D14','D18','D19'],go_requirements=['可审查的主辅方位/报告轨迹/时间戳和同目标关联','有向方位有证据或完整无符号候选分支','固定偏置与导航相关性的独立依据或显式未知量范围','实测数据或明确的非oracle提取合同；与既有合成方位结果区分'],no_go=['用真值位置/舷侧/偏置补齐估计输入','只有不含原始观测与误差定义的精度宣传数值','用未知接收覆盖支持更远布放扫描','将既有Uniform边界约束当作实际校准'],new_hardware='第二被动测向节点；实际可用性/部署成本UNKNOWN；无最小基线认证',status='PROPOSED_NOT_FROZEN_NUMERICAL_EXECUTION_NOT_AUTHORIZED')]
put('NEXT_EXPERIMENT_RECOMMENDATIONS.json',json.dumps(recommendations,ensure_ascii=False,indent=2))
put('NEXT_EXPERIMENT_RECOMMENDATIONS.md',r'''
# 最多两个后续候选，按目标并列

当前证据不足以把深度和水平两个不同目标排成同一优劣次序。我们保留O2与O6作为**条件候选**，不指定总冠军，不给它们分配同一精度分值，也不授权数值执行。

## 候选一：同步不同接收深度的相对强度

它直接检验公共源级消除后还剩什么传播信息，且不必先要求跨时源绝对相位。这是资源/观测条件变化，不是把CZO04再换一个平滑器。核心反例已确定：逐时自由差分增益仍会吸收所有信息；即便增益固定，单静态比例也可能被完全吸收。新的贡献应来自受约束差分增益下、不能被r-z和环境耦合解释的变化。

H3G0提供选择垂向资源的动机，但其2.34–4.47倍增量是**复响应局部深度信息**，不能给本功率比候选背书。REQ05尚未证明抛弃相位后仍有足够信息。此候选需先取得D01–04/D06–07/D09/D15–17，登记噪声底、源公共性、时间相关和差分增益；再由新合同判断能否使用独立可信期望功率或接收数据。原HLA名义同深，不具备软件合成同时垂向采样的能力。

No-Go：只有共深数据却宣称垂向新增信息；差分AGC任意未知；以静态I类SHD替代未经验证的CW接收功率；锁定真实源深或环境取得收缩；前向数值误差无法与候选差异分开。出现这些情况应保留数学潜力并停止性能声明，不回到T2求解器修补。

## 候选二：辅助方位的实际测量链和偏置条件

AUXD1/T3W1已有比较充分的水平条件Demo。下一有价值的问题是：正确关联、左右舷候选、导航和偏置先验能否在可审查输入中成立，而不是再搜索更远的辅助位置。只读取合成公开档案能够复核算法，不能把该档案升级为真实海试输入。

需取得D11–14/D18–19等资料。若只能输出轴角，应保留所有离散支路；若偏置没有独立校准，应剖面或给出可追溯范围，不能提供真实偏置。两节点不要求载波相干，但同步、报告位置、同目标关联和通信仍需证据。该观测没有z列，不承担测深目标。

No-Go：依据真值完成舷侧/目标关联；只提供无法追溯误差定义的精度数字；以扩大基线替代偏置或接收链验证；没有接收可用性却将45–60km几何范围当探测覆盖。资料未到位可保留申请阶段假设设计，不宣称设备失败，也不启动又一轮理想方位MC。

## 其余机制为何暂缓

O1已达到同类模型收口点。O3原HLA同步互谱值得索取和保存，但SP1与H3G2E已有明确限制、可靠复场仍缺；没有新增数据或未知量合同，不重跑旧匹配。O4缺宽带接收和路径身份/能量证据。O5源频/径向采样、模式身份和复场链尚未建立，E1/RX1Z负结果不能换名规避。它们不是被判物理不可能，也不是证明O2严格优于O3。

本次选择依据是不同目标下可检验的条件缺口，不是已验证收益最大。收到资源资料后，候选仍可能都不准入；这也是有效的研究选路结果。
''')
# Numerical audit presentation generated strictly from frozen saved results.
a=json.loads((O/'ALGEBRA_RESULTS.json').read_text());cold=json.loads((O/'INDEPENDENT_COLD_AUDIT.json').read_text());assert cold['PASS']
put('NULLSPACE_AUDIT.md','# 小型代数审计\n\n18组矩阵全部是预登记的无量纲反例/条件示例，绝不是六种真实观测的五维Fisher测量值。精确有理数消元独立核对信息矩阵、秩和零空间，'+str(cold['checks'])+'项通过，0失败。独立实现由同一代理编写，未冒称外部科研审核。\n\n|控制|有效秩|不可辨维数|\n|---|---:|---:|\n'+'\n'.join('|'+z['id']+'|'+str(z['rank'])+'|'+str(5-z['rank'])+'|' for z in a['controls'])+'\n\n完整数值矩阵在ALGEBRA_RESULTS.json，精确零向量在EXACT_NULLSPACES.json，输入在ALGEBRA_SPEC.json。AUX_HORIZONTAL_UPPER_BOUND仅用四个抽象独立水平导数说明z仍为零，不表示实际辅助轨迹已取得四维可靠秩。\n\n共享参考D=[[-1,1,0],[-1,0,1]]作用于单位原始协方差时，DD^T=[[2,1],[1,2]]；忽略非对角项会改变信息。主程序正确协方差下的对比信息与原始空间公共源剖面相等，损失矩阵由有理数验证半正定。平面波和运动方位解析导数与独立差分交叉核对。\n\n这里没有由局部秩推断全域唯一性，没有真实环境Jacobian、实际定位CRLB或接收提取精度；镜像、路径/模式置换、离散别名仍按观测合同保留。')
put('FIVE_STATE_INFORMATION_MATRIX.md','# 五维状态与信息证据\n\n符号：条件=需要完整未知量剖面与可靠前向/观测；未评=本项目没有相应能力证明；零=仅在所列观测模型下的结构零。不同层级不能当共同精度榜。\n\n|观测|距离r|方位theta|深度z|速度v/航向psi|\n|---|---|---|---|---|\n|O1单频包络|有限条件离散增量，非一般锚|间接耦合|干扰，未建立精度|运动耦合，漂移可吞噬|\n|O2垂向相对强度|条件，未评|需其它几何|条件，功率比未评|多时联合，未评|\n|O3互谱/相对相位|平面波零；波导条件|原HLA方向能力，镜像|H3局部复响应条件|动态联合条件，未评|\n|O4多途时延|路径差条件|几何相关，未评|路径差条件|跨时联合，未评|\n|O5频率/相干孔径|offset耦合|需方位/姿态|模式链条件，CZ未评|未知频/参考显著退化|\n|O6辅助方位|条件Demo|有向方位假设|水平模型零|条件Demo；偏置限制|\n\n逐状态完整措辞见FIVE_STATE_INFORMATION_MATRIX.csv；数学/存档/提取/标定四层见OBSERVATION_EVIDENCE_LAYERS.csv。')
put('RESOURCE_AND_RESEARCH_DECISION.md',r'''
# REQ05研究去留：把新增资源对应到可解除的退化

正式结论：**OBSERVATION_RESOURCE_CONDITIONS_ESTABLISHED_NO_NEW_PERFORMANCE_ADMISSION**。

本轮完成六类观测的统一模型、五维状态—干扰参数审查、18组最小代数控制、资源必要条件和合作方资料清单。独立有理数程序124项复核通过。没有新的声场、MC、接收样本、定位精度或最小设备配置认证。

最需要保留的结构结论是：**消除公共未知量后，独立状态信息还必须超出剩余差分未知量和其它待估状态的变化空间。**同时垂向功率比可以消除共同源级，却仍受差分增益限制；同频互谱可以消除共同源相位，却仍受通道复响应限制；多途差分可以消除发射时刻，却仍受路径身份、时延偏差及环境限制。重复采样、归一化和更换优化器不能自动解除这些规范自由度。

原8元HLA若能提供逐元同步复基带，具有采集空间互谱、功率和频率的基本资源；实际接口未知。只保存波束或方位产品时，不能逆推出丢掉的复相位。原名义同深阵不具有同时垂向采样，延长处理或增加软件通道不能替代；原单节点也不能凭处理得到第二节点独立视线。另一方面，附加通道和基线只是可能解除退化的资源，不能直接换算成保证精度。

## 历史成果保留到什么程度

NUM03恢复同一WOA来源族九个名义单元的首深折返几何，不认证真实海域、可检测性或复场。CZO04的自由源级零信息与漂移损失成立；两项较完整离散收紧来自一组镜像几何。RC3原名义有限网格的r=50km、theta=0、v=2m/s、psi=4或5度两候选继续保留；后来的负结果不抹去它，但也不把它当连续全域测距。

H3G0垂向相对同通道共深对照的局部复响应深度信息增量保留，G1/G2E的传播、完整水平支持及提取缺口也保留。Yang浅海正控制和MMAR部分谱/幅度链不是第一CZ非oracle测深。T2/FUNC2和Snell修补链继续关闭。

AUXD1补齐了可复用的合成配对方位输入，T3W1表明宽基线在登记偏置条件下可改善水平弱方向；它不证明设备可获得那些观测。无独立偏置先验时信息明显退化，方位误差也不是所有配置都改善。无符号方位/导航/同步/关联的证据不能由真实标签替代。所有旧百分比阈值仍只是历史轮次研究筛选，不是项目统一验收线。

## 最多两个候选与投入顺序

保留两个按目标并列的条件候选：O2同步多深度相对强度，针对深度/距离敏感性；O6辅助节点方位链，针对水平参数的实际条件成立性。前者没有直接功率比正结果，只因能明确检验源级退化后的新自由度而值得保留；后者已有几何Demo，下一投入应解决测量链而非继续布放扫描。不存在可信统一权重证明谁更优。

先取得资料，再决定是否能够冻结新的有限验证：O2优先查差分AGC、同步深度、噪声底和期望功率模型；O6优先查有符号/多支方位、独立偏置资料、报告导航和目标关联。资料未给出时是申请阶段条件未知，不是实际设备不合格。可写明确假设的研究设计，但不发放真实性能结论。

本轮没有给出最少阵元数、最小垂向孔径、最低带宽、最短观测时间或必要辅助距离。已知的只是必要约束类型和若干历史条件配置；资源下界需可靠的剖面导数、误差范围、全域候选及提取证据共同建立。

## 对项目申请的表述

可表述：第一会聚区多参数约束由传播敏感性与源、通道、环境未知量共同决定；已建立条件几何和局部/水平信息基础，提出以独立观测解除退化、再反推资源条件的研究路径。拟分别验证受控差分增益下的垂向相对强度，以及含实际偏置与多分支的协同方位。

不可表述：已有可靠第一CZ测深、原HLA五维定位、两个候选均已有效、9个独立海域验证、4m或15km为最低资源、仿真假设等于合作方设备规格。所有新增性能仍NOT_EVALUATED。

本阶段交付后STOP；下一候选均为PROPOSED_NOT_FROZEN / NUMERICAL_EXECUTION_NOT_AUTHORIZED。main、历史Gate、原STOP与R4=0%保持不变。
''')
# One review record per material inferential hazard, with primary evidence and limits.
review=[
('公共源消除是否被误作通道校准','O2/O3公式保留差分gain，free-gain零控制通过','H3_NUISANCE;CZO04','PASS'),
('局部秩是否被当成全域唯一','所有矩阵标明无量纲；镜像/模式路径置换独立保留','SUPPORT;SP1','PASS'),
('H3复响应是否被误当功率比正证据','候选O2明确NOT_EVALUATED；收益不继承','H3G0;H3G2E','PASS'),
('浅海文献是否迁移成第一CZ精度','方法身份/源频/速度/offset及提供器缺口单列','PAPER_YANG;PAPER_ERRATUM;PAPER_MMAR;METHOD_ID','PASS'),
('理论可辨是否变实际可提取','四层证据表逐项区分；UNKNOWN保留','RX_CONTRACT;RX1Z','PASS'),
('新候选是否绕开旧STOP','不重扫O1/SP1/频率/T2/基线；候选需新输入','FUNC2;CZO04;T3W1','PASS'),
('协同观测是否暗中解决深度/左右舷','z零列，signed前提及无符号分支单列','AUXD1;T3W1;CLIENT','PASS'),
('资料清单是否成设备验收/发信','仅待交流稿，19项无阈值，未发送','CLIENT','PASS'),
('历史RC3条件正例是否被负结果覆盖','两候选精确范围保留，非一般定位','RC3','PASS'),
('是否把同族/配对样本当独立海域','NUM03一WOA族、CZO04镜像不增加独立几何','NUM03;CZO04','PASS'),
('独立冷审计是否冒称外部审核','不同代码路径Fraction；同一代理身份披露','ALGEBRA_SPEC','PASS')]
table('RESEARCH_CLAIM_REVIEW.csv',[dict(issue=i,disposition=d,evidence=e,status=s,review_type='SEPARATE_PASS_BY_SAME_AGENT_NOT_EXTERNAL_REVIEWER') for i,d,e,s in review])
run=json.loads((O/'ANALYSIS_RUN.json').read_text());ar=json.loads((O/'AUDIT_RUN.json').read_text())
final={'stage':'ENV_CAP_REQ05','decision':'OBSERVATION_RESOURCE_CONDITIONS_ESTABLISHED_NO_NEW_PERFORMANCE_ADMISSION','design_SHA':'80516e049cafdf10dc65896f1cd78276286c390e','parent_SHA':'d658f660d6bd6c9b3a585f0e5b8812c49ee8b0c4','result_SHA':'CONTAINING_COMMIT','branch':'research/env-cap-req05','observables':6,'state_rows':30,'resource_conditions':len(resources),'partner_questions':len(questions),'conditional_candidates':[r['id'] for r in recommendations],'candidate_ranking':'PARALLEL_DIFFERENT_OBJECTIVES_NO_TOTAL_ORDER','algebra_controls':18,'cold_audit_PASS':cold['PASS'],'cold_checks':cold['checks'],'external_scientific_review':'PENDING','new_performance':'NOT_EVALUATED','hardware_status':'UNKNOWN','coherent_CZ_provider':'NOT_ESTABLISHED','R4_percent':0,'main_unchanged':'a9a80191ee4d0f0cecd9db296198606b63a68080','historical_gates_unchanged':True,'formal_runs':{'analysis':1,'cold_audit':1},'budget':{'scientific_wall_s':run['wall_s']+ar['wall_s'],'peak_sampled_process_tree_RSS_bytes':max(run['peak_sampled_tree_RSS'],ar['peak_sampled_tree_RSS']),'external_acoustic_calls':0,'MC':0,'received_data_generated':0},'status':'STOP_AFTER_DELIVERY','next_execution':'NOT_AUTHORIZED'}
put('REQ05_FINAL_DECISION.json',json.dumps(final,ensure_ascii=False,indent=2))
put('GPT_SYNC_REQ05.md',r'''
# REQ05同步报告

结论：**OBSERVATION_RESOURCE_CONDITIONS_ESTABLISHED_NO_NEW_PERFORMANCE_ADMISSION**。本轮完成观测资源论证，没有新增定位/测深精度。父提交d658f660d6bd6c9b3a585f0e5b8812c49ee8b0c4；设计SHA：80516e049cafdf10dc65896f1cd78276286c390e；最终SHA为包含本报告的Commit B，推送后另核对。分支research/env-cap-req05。

六类观测都已给出测量模型、源/通道/环境未知量、零空间及五维状态作用范围。18组无量纲代数控制和124项独立实现冷审计通过，0失败；这是公式与来源身份的核验，不是新CZ精度证书。冷审计由同一代理的独立代码实施，外部科研审核待负责人完成。

主要结论：相对强度可消去共同源级，但逐时自由差分增益仍使其零信息；互谱消去源相位也不等于通道已校准；相对时延、频率和相干孔径各有不可忽略的参考/路径自由度。辅助水平方位没有深度列。处理不能制造垂向孔径、第二节点视线或丢失的复数据。

保留两个**不同目标下并列的条件候选**：

|候选|为什么保留|必须先解决|当前边界|
|---|---|---|---|
|同步多接收深度相对强度|检验公共源级消除后是否仍有r/z信息|同步垂向数据、差分增益/AGC、噪声底及可解释功率前向|没有功率比正结果；不继承H3复响应局部倍数|
|辅助节点方位测量链|已有AUXD1/T3W1条件水平Demo|独立偏置依据、报告导航、侧向多支、同步和目标关联|不再扫基线；真实接收能力UNKNOWN，z未评价|

不设总分选优，不承诺这两个候选均能准入。其余O1单频扫描暂停，O3旧SP1/T2链不重开，O4宽带时延和O5源频/相干孔径等待具体输入条件。NUM03名义几何、RC3两航向候选、H3局部信息和全部负结果各自保留原范围。

交付重点：OBSERVATION_MODELS_AND_NULLSPACES.md；四层OBSERVATION_EVIDENCE_LAYERS.csv；30项FIVE_STATE_INFORMATION_MATRIX.csv；13项RESOURCE_NECESSARY_CONDITIONS.csv；19项PARTNER_DATA_REQUEST.csv及716交流稿；NEXT_EXPERIMENT_RECOMMENDATIONS.md。没有把必要条件变成最低设备数量或实际精度，没有向合作方发送消息。

新增声学求解0，MC0，接收数据生成0。正式代数分析和冷审计合计约0.406s，采样进程树峰值约22.8MiB；完整来源/代码/矩阵可复核。所有项目文件在D盘；历史main和R4=0%不变。两个候选均PROPOSED_NOT_FROZEN、NUMERICAL_EXECUTION_NOT_AUTHORIZED；提交推送核对后STOP。
''')
put('STOP.json',json.dumps({'stage':'ENV_CAP_REQ05','status':'STOP','reason':'AUTHORIZED_REVIEW_COMPLETE; NO_NEXT_NUMERICAL_STAGE_AUTHORIZED','historical_STOPs':'PRESERVED','R4_percent':0,'next_allowed_action':'RESEARCH_LEAD_REVIEW'},indent=2))
# Preserve exclusive execution timestamps without regenerating science.
put('EXECUTION_ORDER.json',json.dumps({'design_SHA':final['design_SHA'],'analysis_started_unix':float((W/'analysis.started').read_text()),'audit_started_unix':float((W/'audit.started').read_text()),'development_precheck':'DEVELOPMENT_PRECHECK.json','outputs':'ALGEBRA_RESULTS.json and INDEPENDENT_COLD_AUDIT.json'},indent=2))
print('Wrote report, matrices and conditional recommendations')
