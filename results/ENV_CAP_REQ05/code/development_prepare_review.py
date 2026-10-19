import pathlib,json,csv,hashlib,subprocess,time,sys,importlib.util
R=pathlib.Path(r'D:\XiaomiMiMoProjects\WSL系统汇聚区目标定位');O=R/'results/ENV_CAP_REQ05'
def put(n,s):(O/n).write_text(s.strip()+'\n',encoding='utf-8')
sources=[
('RC2','results/P2_RC2_kinematic_boundary/P2_RC2_GPT_SYNC.md','单阵共线和弱机动条件；不是所有方位轨迹的不可行证明'),
('SUPPORT','results/R4_RC2_OBSERVATION_SUPPORT/GPT_SYNC.md','保守但宽的U0支持；旧U输入缺口不覆盖后来新生成的AUX档案'),
('RC3','results/R3_FINAL_CLOSEOUT/R3_FINAL_DECISION.json','名义有限网格两航向候选；环境匹配与可提取性边界'),
('RC3_DEPTH','results/R3_RC3_REANCHOR/R3_RC3_DEPTH_ROLE_AUDIT/GPT_SYNC.md','深度角色混合或未解决'),
('YANG_DIAG','results/R3_C2_Yang_SA_depth/R3_C2_2D/R3_C2_2D_GPT_SYNC.md','浅海正则化与泄漏敏感性'),
('YANG_ID','results/R3_C2_Yang_SA_depth/R3_C2_3A/R3_C2_3A_GPT_SYNC.md','傅里叶峰模式身份有限'),
('MMAR_PAPER','results/R3_C2_AR_MMAR/R3_C2_4B_2_PAPER_SPECTRUM/GPT_SYNC.md','有限论文谱机制对照'),
('MMAR_PEAK','results/R3_C2_AR_MMAR/R3_C2_4B_2H_AR_ORDER_DIAGNOSTIC/GPT_SYNC.md','峰选规则歧义'),
('MMAR_HANKEL','results/R3_C2_AR_MMAR/R3_C2_4B_3A_HANKEL_AMPLITUDE/GPT_SYNC.md','幅度桥接不等于完整非oracle测深'),
('METHOD_ID','results/ENV_CAP_FULL_CAMPAIGN/E2_PAPER_METHOD_IDENTITY.md','Yang/MMAR身份与完整性'),
('ERRATUM','results/ENV_CAP_FULL_CAMPAIGN/E2_YANG_ERRATUM_FORMULA_AUDIT.md','距离offset相位与共享差分协方差'),
('PAPER_YANG','literature/yang2015.pdf','PDF第2-3页；已知源频与径向孔径'),
('PAPER_ERRATUM','literature/2018-2.pdf','2018勘误，模式相位offset'),
('PAPER_MMAR','literature/liang2018_mmar_7824671.pdf','第2-3页；已知速度、方位简化、AR与Hankel不同职责'),
('H3_COVERAGE','results/R4_H3_DEPTH_OBSERVABLE_REVIEW/H3_PRIMARY_METHOD_COVERAGE.md','含多途与VLA文献条件，不移植论文精度'),
('H3_NUISANCE','results/R4_H3_DEPTH_OBSERVABLE_REVIEW/H3_OBSERVABLE_AND_NUISANCE_MATRIX.md','源、通道、环境、水平支持联合未知量'),
('H3G0','results/R4_H3_G0_FINITE_VERTICAL_RESPONSE/GPT_SYNC.md','名义理想场局部深度信息2.34-4.47倍；非估计精度'),
('H3G1','results/R4_H3_G1_EXTRACTION_SUPPORT_REVIEW/GPT_SYNC.md','传播与完整水平支持未解决'),
('RX_CONTRACT','results/R4_H3_G1_EXTRACTION_SUPPORT_REVIEW/RECEIVED_SIGNAL_REQUIREMENTS.md','有限频谱和实际录音分层；源/校准/相干窗未知'),
('H3G2E','results/R4_H3_G2E_CONDITIONAL_EXTRACTION/GPT_SYNC.md','压缩响应提取不可靠'),
('E1','results/R4_E1_G0_FREQUENCY_INFORMATION/GPT_SYNC.md','未知源频/漂移/时频参考后增量很小，不是所有频率机制为零'),
('RX1Z','results/R4_SINGLE_HLA_CARRIER_REFERENCE_CONTRAST_PILOT/GPT_SYNC.md','跨窗频差不能映射有效相对运动；动态提供器未认证'),
('SP1','results/R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT/GPT_SYNC.md','单HLA三频静态目录未收紧；数值和网格限制'),
('AUXD1','results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT/GPT_SYNC.md','新配对方位档案、条件总增益；有限终点非支持域'),
('T3W1','results/R4_T3_W1_WIDE_BASELINE_DEMO/GPT_SYNC_T3_W1.md','宽基线收益依赖测量链和偏置先验；不得重扫基线'),
('CLIENT','results/R4_AUX_CLIENT_CONFIRMATION/AUX_CLIENT_CONFIRMATION_MATRIX.csv','实际硬件、同步、关联UNKNOWN；旧数值不是设备验收线'),
('NUM03','results/ENV_CAP_NUM03/NUM03_RESEARCH_DECISION.md','名义首回返几何支持，非相干声场和探测认证'),
('REC01','results/ENV_CAP_REC01/GPT_SYNC_ENV_CAP_REC01.md','条件总TL和几何；相干提供器未就绪'),
('CZO04','results/ENV_CAP_CZO04/CZO04_RESEARCH_DECISION.md','单频源/增益自由度退化；有限离散条件增量'),
('CZO04_SYNC','results/ENV_CAP_CZO04/GPT_SYNC_CZO04.md','设计、执行、审计和范围'),
]
# Locate final FUNC2 report by filename, never infer a newer result from a stale task.
func=list(R.glob('results/*FUNC2*/GPT_SYNC*'))
assert len(func)==1,func
sources.append(('FUNC2',str(func[0].relative_to(R)).replace('\\','/'),'数值修补链关闭；不否定方法族'))
rows=[]
for i,p,scope in sources:
 b=(R/p).read_bytes(); old=subprocess.check_output(['git','-c','safe.directory='+str(R).replace('\\','/'),'show','d658f660d6bd6c9b3a585f0e5b8812c49ee8b0c4:'+p],cwd=R)
 assert b==old or b.replace(b'\r\n',b'\n')==old.replace(b'\r\n',b'\n'),p
 rows.append(dict(id=i,path=p,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),parent_blob_sha256=hashlib.sha256(old).hexdigest(),checkout_crlf_only=(b!=old),scope=scope))
with (O/'EVIDENCE_INDEX.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
put('EVIDENCE_README.md','# 固定来源\n\n父提交 d658f660d6bd6c9b3a585f0e5b8812c49ee8b0c4。表中各文件与父提交Git字节逐一核对；未重跑历史实验。PDF保存原始字节，仅复读方法条件和勘误。DOI在线访问本轮失败，未将网页内容当作新证据。\n\n'+'\n'.join('- **'+r['id']+'** ['+r['path']+'](../../'+r['path'].removeprefix('results/')+')：'+r['scope'] if r['path'].startswith('results/') else '- **'+r['id']+'** ['+r['path']+'](../../'+r['path']+')：'+r['scope'] for r in rows))
# Correct relative link root: REQ05 is two levels below repository root.
s=(O/'EVIDENCE_README.md').read_text();
for r in rows:
 if r['path'].startswith('results/'):s=s.replace('(../../'+r['path'].removeprefix('results/')+')','(../../'+r['path']+')')
put('EVIDENCE_README.md',s)
put('DESIGN_CONTRACT.md',r'''
# ENV-CAP-REQ05 冻结合同

父提交：d658f660d6bd6c9b3a585f0e5b8812c49ee8b0c4。分支：research/env-cap-req05。用户本次指令授权文献/证据与数学审查、独立复核、推送；不授权下一项性能试验。

目标：识别新增观测解除的未知量退化及其必要资源，不设新的项目精度PASS，不计算新CZ定位精度。历史main、Gate和STOP不变，R4=0%。原8元HLA为已有研究资源定义，实际设备型号/数据接口仍UNKNOWN。

## 冻结范围和方法

六类观测全部审查，状态顺序[r,theta,z,v,psi]。统一白化Jacobian与未知量投影，零空间保留；离散镜像/模式置换与局部秩分开。时间、通道、频率的共享噪声不能独立相加。无真值初始化、场生成、优化器、MC或传播求解。

ALGEBRA_SPEC登记18个无量纲、整数Jacobian示例与精确预期秩。它们只测试公式、退化反例和必要条件，不代表真实CZ五维信息矩阵。含常数/线性/自由源级，静态/固定/自由差分增益，平面波相位，自由通道，公共/自由路径时间，源频漂移，自由源相位及辅助方位上限。另验共享参考协方差、数据处理、CSD源相位消除、平面波与运动方位解析导数。所有参数和容差见JSON；正式执行前不依赖真实案例结果调节。

独立冷程序不得导入主程序。以Fraction有理数消元及精确投影重建信息/秩/零空间；独立标量有限差分交叉检查方位公式；逐项验证冻结代码、规格及历史来源SHA。独立实现由同一执行代理编写，不能冒称外部研究人员复审。新实验正式分析一次、冷审计一次，无科学重试。

## 选路规则

先记录直接证据的对象和条件，再作偏序选择，不设加权总分。后续机制至多两个：必须有明确剩余状态敏感性、不同于已关闭路线的新增观测条件，以及可列出的输入/No-Go。不因某条路线历史数字较漂亮而将深度局部Fisher和水平MC精度排名。资料不全只允许条件推荐或并列，不声称最优或设备验收失败。

无新增资源的同一单频扫描、旧SP1网格加密、RX频率重命名、旧T2/FUNC2与Snell修补均不选作下一实验。原HLA互谱采集能力可以审查，但没有新增条件时不能据此重开SP1。辅助节点不继续布放扫描；若推荐，限定真实测量链和偏置/侧向歧义证据。垂向强度的建议不能继承H3复响应信息数值。

## 预算与停止

数学分析+冷审计总壁钟上限600s；单进程树4GiB；新增D盘数据128MiB、Git16MiB；外部声学求解0、MC0、原始接收生成0。代码入口run_once.py带排他启动标记和进程预算监控。开发仅可解析合成自测，记录保留；不重新运行任何历史科学检查。

代码/输入身份、代数或预算失败：封存错误，停止机制数值声明，仍完成数据缺口报告。通过也只准入代数与证据审查，不准入下一性能实验。Commit A冻结后一次运行和一次冷审计；Commit B提交/推送/核对后STOP。
''')
put('OBSERVATION_MODELS_AND_NULLSPACES.md',r'''
# 观测、状态与未知量审查

本报告为模型与证据审查。x=[r0,theta0,z,v,psi]；q(t)=(r0 cos theta0,r0 sin theta0)+tv(cos psi,sin psi)，接收位置p_i(t)来自报告导航，不能取生成真值。环境e含SSP/底质/海深，阵位姿也可能未知。下文P代表传播**接收声级响应**（负TL），不是把TL正负号混用。新控制是无量纲代数，不是名义CZ声场。

## 统一判据

固定协方差Sigma下，W Sigma W^T=I，A=W Jx，B=W Jeta。局部有效信息F=A^T(I-BB^+)A；某方向u可辨需Au不在col(B)。判断单个r或z须进一步剖面所有其它状态列。奇异方向为不可辨/无界，不赋零方差。若协方差自身依赖状态或源参数，均值Jacobian不是完整Fisher，须使用完整似然；本轮不假造该项。未知量的独立先验是额外信息，不能算作新接收信息。

变换T的数据处理不等式要求同一原始似然与未知量、正确变换协方差。线性D有Sigma_D=D Sigma D^T；共享参考带来非零交叉项。非线性比值/相位在弱场、噪声偏差、相位绕转时不能直接套高SNR导数。一个局部正秩也不能排除全局别名。

## O1 单频声级时间包络

观测y_k=P(x,t_k,e)+a_k+g_k+n_k。C0中a+g为未知常数，投影去1；C1为常数+线性，投影去[1,t]；C2逐时自由，Jeta=I，F=0。源级和增益的单独值本身也有规范自由度，只能辨其和。若P的距离导数在常数、漂移及其它状态/环境导数张成空间中，距离增量为零。

CZO04是实际已有的条件模板证据：324确定性均值条件、主C1仅G02及镜像具有数值/集合一致的离散收紧，不能当覆盖率。29个可解释配对案例中漂移后新增距离信息中位比8.09e-6，范围不能扩展为全环境结论。全自由源级的退化是解析结构，不靠该比例成立。参考CZO04、NUM03。

## O2 同时刻不同接收深度的相对强度

在同一源激励、可比较的窄带功率平均下，I_i=|S|^2 |g_i|^2 H_i(x,e)+N_i。噪声背景扣除并可解释后，d_i=10log10[(I_i-N_i)/(I_0-N_0)]=10log10(H_i/H_0)+gamma_i-gamma_0。共同源功率消去；源方向性随路径/节点不同不自动为公共系数。有限样本比值有偏、接近噪声底可不可定义，须保留缺测而非强制取log。

静态单快照且每个差分增益自由时，未知量与d维数相同，可吸收全部空间强度形状。多快照且差分增益真正恒定时，时间变化可能保留，但静态部分仍不可辨；逐时自由差分增益使全部信息归零。接收深度/姿态、环境和目标深度改变也可混淆距离，不得固定真z/e。至少要出现不被这些干扰吞噬的独立导数方向；这不推出最少阵元数、最小孔径或已达到米级精度。

原HLA名义同深，不能靠软件制造同时刻垂向采样。沿程改变阵深也不自动等于同时垂向观测（源/增益/运动不同）。H3G0的B1/B2增益来自完整复响应，不是功率比验证；新第一CZ存档I类能量只提供条件期望能量背景，不能直接当CW实测功率。候选深度强度比的信息目前NOT_EVALUATED。

## O3 同频阵元互谱和相对相位

窄带单目标模型Y=S diag(g)h(x,e)+n。充分统计是否为CSD取决于源/噪声似然，不能统一假定。扣除噪声项后C_ij=E|S|^2 g_i g_j* h_i h_j*；公共源相位消去，公共源功率和通道复响应仍在。相位arg C_i0=arg(h_i h_0*)+arg(g_i g_0*)，无需跨快照公共源相位连贯，但阵内共钟及相对响应必须受控。

未知复S消去两实维规范；一个快照中任意非零g_i可取g_i'=g_i h_i(x)/h_i(x')，使任何新候选产生相同Y。固定跨时g不能随意逐时重拟合，因此可能保留动态信息。自由模态复幅度也能吸收源深系数，除非其约束有独立物理根据。共同参考协方差必须保留；归一化不会制造原始复谱没有的信息。

远场单平面波h_i=exp(i k a_i sin theta)且公共幅度未知时，原水平阵控制仅有方向导数，没有r/z导数。真实波导多模态或近场可能另有空间敏感性，但要经前向与提取验证。水平各向同性、共x-z平面阵对源y→-y距离不变；HLA、倾斜及短垂直阵都可能保留该镜像。局部Jacobian不会看见这种离散对称；保留多支，不能借真值选舷侧。

已有SP1在原8元14m三频静态样本中三方法均未建立完整距离包络收紧；不等于所有互谱物理信息为零。H3G0是名义局部信息，G2E压缩提取仍不可靠；旧T2精度链关闭。原8元若能交付同步逐元复基带，已有采集此类量的基本通道资源；仅有波束功率或方位产品无法逆推出丢掉的复相位。参考SP1、H3G0/G1/G2E。

## O4 多途相对时延

路径模型Y(f)=S(f)sum_l A_l exp(-i2pi f(t0+tau_l(x,e)))。若多途已可分辨/关联，则Delta tau_l=tau_l-tau_0消除发射时刻及同一接收时钟的公共偏置；不同接收器TDOA仍保留差分时钟。未知S的自相关主瓣/周期性、路径重叠、色散、路径置换、底反射相位和环境误差限制能否提取这些时延；不允许把射线路径到达时刻直接送入接收器。

J_Delta=D J_tau；每路径自由延迟误差可吸收全部状态效应。单一CW只给模2pi相位，不能单独确定不模糊延迟；常规相关分辨尺度约1/B只是带宽设计量级，不是通用CRLB或最少带宽认证。若不同路径具有非共线的r/z导数才可能解除耦合，同时必须保留模式/路径身份分支。既有CZ时延文献的大深度跨度接收和宽带激励条件不能归给原14m HLA。NUM03只给几何，不给已分离到达峰和能量。

## O5 运动频率与相干径向孔径

单有效分量一阶模型f_rx=f_s(t)(1-q/c_eff)+b_clock(t)，q为路径相关径向率。未知常f_s与常径向率存在乘性规范；源漂移/参考漂移能吸收低阶运动变化。两个节点共享源频可提供差分变化，但各自参考漂移仍是未知量；旧E1没有支持在其完整未知量模型下获得显著增量。多模态的c_eff不能任意换成1500并宣称精确q。

相干孔径Y(t)=S(t)sum_m A_m(z,e)exp[-ik_m r(t)]。任意时间源相位可吸收单通道的相位轨迹；相干处理需受约束的源/时频参考及可确定的径向采样位置。距离offset使模式系数乘exp(i k_m delta)，跨模相位变，故功率谱不变不能推出相干深度评分不变。常规Fourier尺度Delta k~2pi/L是孔径设计参考，不是AR不能超越的普适下界；模式身份/幅度/环境错误仍需检查。

Yang2015原文PDF第2-3页要求已知发射频率用于径向孔径，并依赖名义模态深度函数；2018勘误保留offset。Liang2018第2-3页使用已知速度和简化真实方位，AR峰用于k、Hankel用于复幅度。原文例子及历史部分复现不能自动变为第一CZ非oracle提取。RX1Z频差可重复，但同航段变化不优于零预测；不把它回灌H2。参考PAPER_YANG/PAPER_ERRATUM/PAPER_MMAR、METHOD_ID、E1、RX1Z。

## O6 辅助节点独立方位

beta_nk=atan2(q_y-p_ny,q_x-p_nx)+b_c+s_n b_d+epsilon_nk。节点时偏、姿态、导航偏差和目标关联均需额外参数/离散假设。令d=q-p，J_beta=(-d_y,d_x)/||d||^2 * dq/dx；dq/dr=(cos theta,sin theta)，dq/dtheta=r(-sin theta,cos theta)，dq/dz=0，dq/dv=t(cos psi,sin psi)，dq/dpsi=tv(-sin psi,cos psi)。固定偏置的干扰列为1和s_n，任意逐时偏置则可吸收全部方位。

两个非平行的已标定有向方位在静态二维几何可约束水平位置，多时刻可进一步约束运动；近共线条件很差。没有垂直方向观测，z列严格为0。未知固定偏置、导航相关误差可能侵占弱方向；Uniform边界给出的约束不是数据自己学到的标定。无符号轴角有额外离散分支，双节点本身不保证其消歧。

AUXD1有新配对公开方位档案；T3W1宽基线有条件改善r/v/psi，但方位P95并非总改善、偏置边界依赖强。它是新增观测+几何+登记处理器的条件总收益，不是纯孔径实验，更不是实际设备。辅助节点无需与主节点保持声学载波相干，却需可用时间对准、导航和同目标关联。参考T3W1、AUXD1、CLIENT。

## 资源反推的可用数学形式

在明确的局部固定协方差模型下，剖面后距离方向a_r=(I-P_[B,A_other])A_r，若a_r=0，单纯降低白噪声或增加重复样本不能解除该规范。若非零，局部尺度1/||a_r||仍需前向导数、相关噪声和全域别名验证，不能当已达到精度。对有限候选的剖面间隔d_min，只有共同误差范数小于d_min/2等充分条件才能保证该有限目录最近模板选择；当前没有真实观测误差界，故不提供设备最低规格。

所有模型结论限于写明的源、环境和误差假设；数学必要条件并不充分。第一CZ必须存在是物理背景，5000m海深/50–60km/200m均不是一般硬约束。NUM03的九单元仍为同一WOA来源族；传播路径能量归因、相干场和实际可检测性未认证。
''')
print('Evidence sources',len(rows))
