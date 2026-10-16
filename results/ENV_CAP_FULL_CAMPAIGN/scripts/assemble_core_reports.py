"""Assemble ENV-CAP reports from sealed outputs. No scientific computation."""
from env_provider import *
import shutil,datetime

def md(name,s):(O/name).write_text(s.strip()+'\n',encoding='utf-8',newline='\n')
now=datetime.datetime.now(datetime.timezone.utc).isoformat()
case=list(csv.DictReader((O/'E1_CZ_ENVIRONMENTS.csv').open()));cold=json.loads((O/'E1_INDEPENDENT_COLD_AUDIT.json').read_text());recold=json.loads((O/'E2_E4_INDEPENDENT_COLD_AUDIT.json').read_text())
js(O/'ENV_CAP_FINAL_DECISION.json',{'campaign':'ENV_CAP_20261016','parent_SHA':'a9a80191ee4d0f0cecd9db296198606b63a68080','branch':'research/env-cap-20261016','decision':'NO_GO_CURRENT_CZ_METHOD_PERFORMANCE_DEMO','campaign_route':'PRE_REGISTERED_NEGATIVE_AND_ANALYTIC_RESOURCE_CLOSEOUT','E0':'PASS_SOURCE_AND_CONTRACT_FREEZE','E1':'ADMISSION_WITHHELD_COLD_RECONSTRUCTION_MISMATCH','E1_main_provisional_CZ2_count':9,'E1_provisional_source_families':1,'E1_final_admitted_CZ2_count':0,'E1_independent_checks':274,'E1_independent_FAIL':12,'E2':'IMPLEMENTATION_ZERO_CONTROL_FAILED_INT32_OVERFLOW','E2_stored_shallow_depth_control_PASS':12,'E2_total_control_PASS':25,'E2_total_control_FAIL':1,'E3':'NOT_EVALUATED','E4':'ANALYTIC_SCALES_AND_HISTORICAL_RESOURCE_REFERENCE_ONLY','E4_new_CZ_conditional_feasible_count':0,'E5':'REPORT_AND_NEGATIVE_EVIDENCE_ARCHIVAL','readonly_E2_E4_cold_checks':126,'readonly_E2_E4_cold_FAIL':0,'new_MC_configurations':0,'new_received_pressure_records':0,'coherent_CZ_provider':'NOT_ADMITTED','Yang_MMARDemo':'NOT_EVALUATED_FIRST_CZ','T2_FUNC2_historical_STOP':'UNCHANGED','R4_percent':0,'five_state_precision':'NOT_ESTABLISHED','hardware':'UNKNOWN','time_sync_and_association':'NOT_VALIDATED','automatic_next_experiment':False,'final_STOP_required':True,'not_inferred':['Physical absence of a firstCZ in SCS','Physical impossibility of Yang/MMAR or T2','All high-resolution algorithms limited by ordinary Fourier Rayleigh spacing','Historical T3 signed bearings acquired by current HLA equipment']})
source=[]
for p in sorted((W/'sources').rglob('*')):
 if not p.is_file() or p.suffix.lower() not in ['.nc','.json','.csv','.npy']:continue
 source.append({'source_id':p.relative_to(W/'sources').as_posix(),'path_or_url':str(p),'source_name':p.name,'source_kind':'ARCHIVED_ENVIRONMENT_SOURCE_OR_DIAGNOSTIC','version_or_date':'WOA23/ETOPO2022/archived2026-09; see linked metadata','size_bytes':p.stat().st_size,'sha256':sha(p),'coordinate_system':'depth_positive_down_for_SSP; ETOPO_elevation_positive_up','lat':'17 or18 or25 permetadata','lon':'118 or114 or140 permetadata','depth_max_m':'PER_FILE_METADATA_NOT_ASSUMED5500_VALID','access_status':'LOCAL_COPIED_READ_ONLY','quality_notes':'Original-source-grid, legacy NCSS and reconstructed products remain distinct; filename alone is not scientific admission.'})
source.extend([{'source_id':'NOAA_WOA23','path_or_url':'https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/','source_kind':'OFFICIAL_PRODUCT_METADATA','access_status':'WEB_READ_2026-10-16','quality_notes':'Climatological analyzed T/S, not CTD; physical input bytes are the SHA-indexed extracts.'},{'source_id':'ETOPO2022','path_or_url':'https://www.ncei.noaa.gov/products/etopo-global-relief-model','source_kind':'OFFICIAL_PRODUCT_METADATA','access_status':'WEB_READ_2026-10-16','quality_notes':'Bathymetry is not sediment acoustic property or local survey certification.'},{'source_id':'BELLHOP_OFFICIAL','path_or_url':'https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node61.html','source_kind':'SOLVER_PRIMARY_DOCUMENTATION','access_status':'WEB_READ_2026-10-16','quality_notes':'Ray/beam tracing, inputunits, ray thenTL, beam-doubling convergence; not modal completeness or true ocean precision.'},{'source_id':'INDEPENDENT_25N140E','path_or_url':str(W/'sources/independent_25N140E/ACCESS.json'),'source_kind':'FIXED_NEW_SOURCE_ATTEMPT','access_status':'T_S_REQUEST_HTTP400_BATHY_AVAILABLE','quality_notes':'This request failure does not prove regional WOA data do not exist. No site substitution or invented deep SSP.'}])
table(O/'E1_ENV_SOURCE_INDEX.csv',source)
neg=[]
for c in case:
 if c['CZ_class']!='CZ2_RAY_TL_STABLE':neg.append({'case_id':c['case_id'],'scope':c['CZ_class'],'reason':c['reason'],'physical_impossibility_inferred':False,'new_solver_retry':False})
neg.extend([{'case_id':'SCS17_ANNUAL_z200','scope':'COLD_RECONSTRUCTION_FAIL','reason':'Zero-energy11.375km I0 window vs sliding sum/dBfloor; global28.888dB discrepancy. Focus-band difference2.84e-14dB does not retroactively pass global audit.','physical_impossibility_inferred':False},{'case_id':'ENV_SCS_WOA_18N114E','scope':'SOURCE_INCOMPLETE_FOR_NEW_ADMISSION','reason':'Cross-study deep endpoint diagnostic and legacy NCSS; no admitted full-depth original-source identity.'},{'case_id':'ENV_INDEPENDENT_REGIONAL_25N140E','scope':'T_S_REQUEST_UNAVAILABLE','reason':'Fixed request HTTP400, bathymetry obtained; no acoustic inference from missing data.'}]);table(O/'E1_CZ_NEGATIVE_EXAMPLES.csv',neg)
needs=[('AREA_AND_SEASON','intended basin, coordinates, season, source/receiver depth bounds','FirstCZ location and ray-family topology','Client/site not selected; climatology only'),('SSP','full-depth T/S/CTD, pressure/location/time, uncertainty and temporal covariance','admitted CZ and coherent modal prediction','WOA is climate mean; missing deep data must not be silently completed'),('BATHYMETRY','track bathymetry and uncertainty, sound-channel depth vs bottom','zero-/bottom-bounce families and detection coverage','ETOPO source-point flat extension is conditional'),('SEDIMENT','sourced sound speed/density/attenuation/layering/roughness or justified ranges','coherent-field boundary completeness','Unknown; rigid bottom is mathematical control'),('CW_SOURCE','available tonal lines, f0 knowledge, coherence, amplitude/phase drift','Yang/MMAR radial-coordinate and phase extraction','Source frequency/gain cannot be handed to estimator as truth'),('RECORDING','time-synchronized raw channels, rate, calibration, units, receive bandwidth','actual nonoracle observables and shared covariance','No new CZ pressure records in campaign'),('ARRAY','actual element XYZ/depth, attitude/roll, phase/gain vsfrequency','geometry and depth-sensitive relative field','8/14m base; vertical and transverse aperture separately costed'),('MOTION_NAV','platform/target prior, navigation covariance and maneuver constraints','available relative radial aperture Lr and horizontal branches','platform speed*time is not Lr'),('TIME_REFERENCE','clock offset/drift and any independent source-frequency reference','phase / radial feature gauge','Not equipment-validated'),('SIDE_ASSOCIATION','signed-bearing provenance, ambiguity handling and same-target association','joint horizontal support','T3 assumes signed bearings, cannot certify HLA side resolution'),('AUXILIARY','node count, maintainable5/7/15.3km geometry, telemetry/nav/time/association','historical T3 design transfer','No automatic actual capability from simulation'),('PRIORITIES','which parameters and research grades matter; tolerate refusal/intervals','resource Pareto and experiment scope','No universal5%/10% acceptance; depth can remain interval/unknown')]
table(O/'E4_NEEDS_FROM_CLIENT.csv',[{'item_id':i+1,'category':a,'needed_evidence':b,'reason':c,'current_limit':d,'hardware_or_data_status':'UNKNOWN_OR_NOT_PROVIDED','necessity_not_procurement_requirement':True} for i,(a,b,c,d) in enumerate(needs)])
ledger=[{'evidence_id':'ENV_CAP_E0','stage':'E0','status':'PASS','level':'SOURCE_IDENTITY_AND_CONTRACT','source':'E0_PRECHECK.json','claim':'32 package identities,28 historical sources; independent source selection fixed','limit':'Not science admission'},{'evidence_id':'ENV_CAP_E1_MAIN','stage':'E1','status':'PROVISIONAL_NOT_ADMITTED','level':'RAY_ENERGY_DIAGNOSTIC','source':'E1_CZ_ENVIRONMENTS.csv','claim':'9 provisionalSCS cases in1 source family','limit':'FinalE1 audit hold; flat rigid nominal model'},{'evidence_id':'ENV_CAP_E1_COLD','stage':'E1','status':'FAIL_UNCHANGED','level':'RAW_NUMERICAL_RECONSTRUCTION','source':'E1_INDEPENDENT_COLD_AUDIT.json','claim':'274 checks,12FAIL; cancellation/floor sensitivity retained','limit':'Not evidence of physicalabsence'},{'evidence_id':'ENV_CAP_E2_PAPER','stage':'E2','status':'STORED_CONTROL_RECONFIRMED','level':'HISTORICAL_SHALLOW_ORACLE','source':'E2_PAPER_DEPTH_CONTROL_CURVES.csv','claim':'12 archived depth-score comparisons match','limit':'Not CZ/nonoracle extraction or new eigenfunction certification'},{'evidence_id':'ENV_CAP_E2_CONTROL','stage':'E2','status':'IMPLEMENTATION_INVALID','level':'CONTROL_FAILURE','source':'E2_FINAL_ADJUDICATION.json','claim':'int32 mirror overflow caused1FAIL /26','limit':'No physicalmechanism failure inferred; no retry'},{'evidence_id':'ENV_CAP_E3','stage':'E3','status':'NOT_EVALUATED','level':'NONE','source':'E3_PERFORMANCE_BY_CASE.csv','claim':'30 cells masked with causes,0 attempts','limit':'No precision/coverage estimates'},{'evidence_id':'ENV_CAP_E4','stage':'E4','status':'LIMITED','level':'L0_AND_HISTORICAL_REFERENCE','source':'E4_RESOURCE_INVERSION.csv','claim':'conditional design scales and historical8point resource tradeoffs','limit':'No newCZ feasiblecells or minimal hardware'},{'evidence_id':'ENV_CAP_COLD_E2E4','stage':'E5','status':'READONLY_RECONSTRUCTION_PASS','level':'SCOPE_LIMITED_AUDIT','source':'E2_E4_INDEPENDENT_COLD_AUDIT.json','claim':'126 checks0FAIL','limit':'Never repairs E1/E2 failed gates'}];table(O/'ENV_CAP_EVIDENCE_LEDGER.csv',ledger)
md('E1_CZ_PHYSICS_EVIDENCE.md',r'''
# E1物理证据与最终准入边界

本轮主程序的9个CZ2临时标签没有成为最终环境族准入：独立冷审计274项中12项失败，E1_FINAL_ADJUDICATION将正式已准入数保守记录为0。八个SCS单元逐项一致是可保留的诊断事实，不用于绕过全阶段审计规则。

声学几何采用范围不变C-linear SSP。若向下发射掠射角为alpha，水平慢度p=cos(alpha)/c(zs)。深层折返点满足c(zt)=1/p；一次返回同源深的水平距离为R=2 integral[zs,zt] p*c(z)/sqrt(1-p²*c²(z)) dz。每段线性声速可解析积分，折返点的可积奇性不通过粗深度网格绕开。独立冷审计再用z=zt-u²的数值积分核对。BELLHOP的第一向上穿越点只在此前不触海面、不触海底时被接受。由发射角—回返距离曲线的内部最早折叠候选联合TL能量增强检验，区别于单独挑最高TL峰。

固定频率200Hz，范围10–120km，三组源深150/200/250m，接收带随源深取±25m。I0/I1/I2束数2001/4001/8001、步长50/25/12.5m，深度100–300m每2.5m，距离每125m。先在接收带内对线性能量平均，再做固定9点滑动平均。焦区±5km的两级相邻数值差必须中位≤1dB、P95≤2dB，峰位差≤1km、宽度差≤2km；增强需超过max(0.5dB,2倍数值P95)。这些是本轮条件门槛，不是通用第一CZ定义。

17N118E三套输入均来自WOA23原始四角T/S，使用GSW的p_from_z、SA_from_SP、CT_from_t和sound_speed。共同有效至4000m，海深3909.75m在数据内插范围，January上层与Winter(Jan–Mar)/Annual深层的拼接、全Annual子例分别列出。三者属于同一个气候资料来源，不是独立海域。源点ETOPO水深被保持为平底、底边界为刚性，是本轮名义合同，不是实际海底参数或整段地形。

主计算临时峰位为48.375–51.375km，半显著度宽3.375–5.75km。以JW、源200m为例，解析回返折叠约47.756km，接收带能量峰49.125km、宽4.75km。回返折叠与有限频率、有限深度带能量峰不是同一个量。无损几何能量结果未包含实海洋吸收／底损的不确定性，因此其TL不能用于装备检测预算。字段和所有36条TL曲线保留，不据此发布真实CZ定位性能。

冷审计差异主要来自全零／极弱窗。SCS年平均200m的最大28.888dB差在11.375km；焦区内两种重建差约2.84e-14dB。这是对失败作用域的只读定位，不以缩小审计区间来改判通过。E-STD三例既没有满足本轮内部折叠选择，也未闭合TL数值门，部分记录含类似dB地板敏感性；不得据此认定E-STD所有射线族均无CZ。常速三个控制在登记角域内无折射回返，不排除反射路径。

首CZ已存在是后续主研究硬要求。当前只有条件物理诊断，未形成完整审计准入环境；相干复模态、独立传播精度、真实底质和跨海域推广均未认证。

来源：[WOA23官方元数据](https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/)、[ETOPO2022](https://www.ncei.noaa.gov/products/etopo-global-relief-model)、[BELLHOP官方说明](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node61.html)。具体物理输入与原始输出以本目录CSV/JSON和D盘SHA索引为准。
''')
md('E2_PAPER_METHOD_IDENTITY.md',r'''
# Yang与Liang方法身份

Yang2015，JASA138:1678–1686，DOI10.1121/1.4929748，研究运动源形成的径向合成孔径与模态深度结构。论文要求可相干处理的CW、足够径向位移，已知发射频率用于径向运动换算；它不是单次14m HLA自动联合恢复五维状态。正文、2014参考环境与2018勘误PDF均按title/DOI及SHA识别。原浅海约88m、350Hz、源4/50m与接收18/70m、4990m孔径是正控制，不是第一CZ。

本轮不重跑KRAKEN浅海环境；直接读取存档MOD、真实模式处复谱g与深度评分，另写循环求和重建。四组几何×三种正则比例共12个评分比较一致；独立冷核对从原MOD重建全部曲线。该成果只确认历史oracle计算仍可重建，不扩展为新的传播精度或真实模态提取能力。

Yang2018 Erratum，DOI10.1121/1.5081712，对距离offset的相位作修正；全文PDF已阅读。起点对齐不能从“无需绝对距离”自动推出，见公式审计。原文幅度正则与数据shading也不能靠真实传播幅度代替后宣称已提取。

Liang2018 MMAR，Mathematical Problems in Engineering7824671，DOI10.1155/2018/7824671，AR用于波数位置，广义Hankel谱用于复幅度，再经模式顺序匹配形成深度评分。AR峰高不是模态幅度。原文含已知径向速度和简化方位条件；不能把它们悄悄借给未知目标问题。历史无噪模式次序及Hankel幅度桥接的条件成果保持，但20dB多峰结构、自动峰选和完整D(z)尚未闭合。本轮身份为PAPER_IMPLEMENTATION_PARTIAL，没有声称严格复现图5–7。

新CZ前向仅有射线／非相干能量场，既不输出认证完整模式，也不具新SSP与底边界下的高精度复函数。旧E-STD模态窗口与有损函数限制独立存在。因此Yang/MMAR在新环境的模态密度、模式提取率和深度精度全部NOT_EVALUATED。即使新E2镜像代码没有失败，这个L3缺口也仍存在。

E2的总控制25/26通过，失败的静态HLA镜像控制由int32平方溢出造成，NaN完整保留。它限制当前实现准入，既不撤销存档12项一致性，也不否定镜像的解析恒等关系。没有修补重跑。

主要原文：仓库literature/yang2015.pdf、2018-2.pdf、liang2018_mmar_7824671.pdf。版本与字节见E2_UNIT_FREEZE.json。出版信息可参阅[2018勘误书目信息](https://pubmed.ncbi.nlm.nih.gov/30599696/)；实际公式依据本地原文而非摘要。
''')
md('E2_YANG_ERRATUM_FORMULA_AUDIT.md',r'''
# Yang公式与未知量审计

采用exp(+i*omega*t)约定，单模接收项写作sqrt(2pi/(k*r))*phi_s*phi_r*exp(-i*k*r-alpha*r-i*pi/4)，alpha≥0。理想sqrt(r) shading和在真实k处的积分给g=phi_s*phi_r*exp(-alpha*r1)*(1-exp(-alpha*L))/(alpha*k)；alpha→0的极限为phi_s*phi_r*L/k。新代码以expm1计算，八个固定存档模式的有限积分差小于冻结1e-6限值；这只是单模公式控制，不是多模自动识别或完整传播认证。

若积分坐标相对真实距离存在delta，模式峰变为g_m(delta)=g_m(0)*exp(i*k_m*delta)。其模平方不变，但不同模式间的相位差增加(k_m-k_n)*delta，故相干深度评分一般变化。本轮两模式控制的功率差约2.22e-16，相干和平方变化约0.6436。必须搜索／剖面delta或提供独立参考，不能用真实起始距离对齐后称未知距离已解决。相位预算在资源表中只是事先登记的0.1rad设计坐标。

静态单通道每频率一个自由公共复源系数可以吸收任意标量传播响应；若每个模式另有完全自由复增益，源深对模式系数的影响也可被吸收。多通道、跨时或受约束增益的情形需要重新建立联合未知量合同，不由这些控制一概否定。

单有效分量、恒定径向率模型中，接收频率f0*(1-q/c)具有源频与q的规范不确定性；改变q同时改变未知f0可保持接收频率。该恒定模型的零信息不等于所有非匀速、多模态或多节点观测都零信息，旧E1的很小正信息仍保持原结论。

五个相邻窗差分由六个窗噪声构成，单位原始噪声时协方差对角2、邻对角-1，相关系数-0.5。保留该协方差后变换信息不超过原始信息。不能将同一录音的bearing、radial和modal派生特征独立相加；缺联合协方差时只能做互斥对照。

静态共线／共x-z平面阵在水平各向同性介质中的左右镜像距离解析相等。当前数值控制却因Windows默认int32产生平方溢出而失败，IMPLEMENTATION_INVALID如实保留；此处只陈述解析关系，不用重新计算来修复本轮Gate，也不声称任意主动机动之后仍必然相同。

本轮未构造新的五维Fisher矩阵，未给不可辨方向伪逆零方差。新环境的未知水平状态与深度联合支持、源频、offset及传播未知量剖面均NOT_EVALUATED。
''')
md('APPLICATION_USABLE_CLAIMS.md',r'''
# 可用于项目申请的表述

本项目拟研究第一会聚区环境、可提取传播观测与有限接收资源如何共同决定多参数约束能力，采用环境准入、独立信息、实际接收提取、多分支估计及资源反推逐级检验。

已有Yang浅海模式／深度控制和MMAR部分波数—幅度机制基础；它们揭示了模态身份、径向相干孔径、源频和距离offset的条件依赖。ENV-CAP对浅海存档评分的12项重建一致，但没有取得新第一CZ深度精度。

基于可追溯WOA23资料与ETOPO点深重建的17N118E名义环境给出首折射会聚区的条件诊断线索。应同时披露本轮TL冷重建未闭合，环境族准入暂缓；不得把图中49km量级的诊断焦区当成真实海域已验证结果。

历史T3研究已提供双节点有符号方位条件下的资源—水平误差参照，较宽基线改善距离和速度弱方向，同时存在偏置支持依赖、方位尾部和实际测量链缺口。该结果可支撑对协同几何与系统偏置共同作用的研究问题，不能证明新增基线本身是算法创新。

申请拟突破的问题可包括：可追溯环境族中传播特征的可提取边界；未知源与offset、校准及环境干扰下独立状态方向的保留；面向多分支与拒判的资源—信息设计。以上是拟研究内容，当前未全部完成Demo验证。
''')
md('APPLICATION_FORBIDDEN_CLAIMS.md',r'''
# 当前证据不支持的表述

- 不能写“已认证深海第一CZ环境族并实现Yang/MMAR深度精度”；正式CZ2准入数为0，E2实现控制失败，E3未执行。
- 不能将9个临时正标签称为9个独立环境／海域，或把一套WOA来源的子例称跨海域稳健。
- 不能将气候平均、跨研究拼接、平底刚性底或几何范围内称为实测环境、真实探测覆盖或海底完整性证明。
- 不能将已知源频、真实模式、真实径向位移、delta=0、真实水平状态的控制称为接收侧已取得观测。
- 不能将普通Fourier的2pi/Lr尺度称为全部高分辨算法不可逾越的下界，也不能将14m物理阵长视为1200s的径向孔径。
- 不能将历史T3有符号方位误差表称为原单HLA已解决左右舷、真实系统已校准、同步关联已验证或深度有信息。
- 不能将经验分位档反推成已证明的最小装备、全分支覆盖、真实95%成功率或统一5%/10%甲方验收线。
- 不能将完整性核对PASS改写成E1/E2科学Gate PASS，或将实现失败上升为Yang、MMAR、T2或单HLA物理不可行。
- 不能宣称完成新五维联合定位、米级实际测深、新CZ Monte Carlo性能或为R4计入完成百分比。
''')
print('Core reports, evidence/source registries and final scientific decision written.')
