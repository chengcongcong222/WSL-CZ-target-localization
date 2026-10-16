"""NUM03 delivery from saved data only: no main/reference numerical functions are called."""
from pathlib import Path
import json,csv,hashlib,time,threading,collections
import numpy as np
import psutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_NUM03');start=time.time();done=threading.Event();peak=[0]
def monitor():
 while not done.wait(.1):
  p=psutil.Process();peak[0]=max(peak[0],sum(q.memory_info().rss for q in [p]+p.children(recursive=True) if q.is_running()))
threading.Thread(target=monitor,daemon=True).start()
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,a):p.write_text(json.dumps(a,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def text(p,a):p.write_text(a,encoding='utf-8',newline='\n')
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def table(p,a):
 with p.open('w',encoding='utf-8',newline='') as f:
  d=csv.DictWriter(f,list(dict.fromkeys(k for v in a for k in v)));d.writeheader();d.writerows(a)
f=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'));p=json.loads((O/'NUM03_PRIMARY_RESULTS.json').read_text(encoding='utf-8'));cold=json.loads((O/'NUM03_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'));checks=cold['checks_detail'];rows=read(O/'NUM03_RETURN_CURVES.csv');refs=read(O/'NUM03_FULL_DOMAIN_INDEPENDENT_REFERENCE.csv');mp=read(O/'NUM03_MP_REFERENCE.csv');ode=read(O/'NUM03_ODE_COMPARISON.csv');first=read(O/'NUM03_FIRST_CANDIDATE_INDEPENDENT.csv');old=json.loads((R/'results/ENV_CAP_CZROB02/CZROB02_PRIMARY_RESULTS.json').read_text(encoding='utf-8'));summary=[];domains=[]
numchecks=[a for a in checks if any(s in a['check'] for s in ['_full_angle','_coeff','_MP_','9_cases_72018'])];numeric=bool(len(refs)==72018 and all(a['PASS'] for a in numchecks) and json.loads((O/'DEVELOPMENT_GATE.json').read_text(encoding='utf-8'))['PASS'])
for c in p['cases']:
 name=c['case_id'];oc=next(v for v in old['cases'] if v['case_id']==name);out={'case_id':name,'source_z_m':c['source_z_m'],'source_family':c['source_family'],'representation_range_shift_m':c['representation_range_shift_m'],'representation_turn_shift_m':c['representation_turn_shift_m']};geos=[]
 for kind in ['C_LINEAR','PCHIP']:
  mr=c['methods'][kind];rr=[a for a in refs if a['case_id']==name and a['SSP_method']==kind];pr=[a for a in rows if a['case_id']==name and a['SSP_method']==kind];valid=[a for a in rr if a['status']=='FINITE'];cc=[a for a in checks if a['check'].startswith(name+kind)];cp=bool(all(a['PASS'] for a in cc));fc=mr['first_refined_candidate'];op=oc['methods'][kind]['folds_count'];geo=bool(cp and mr['primary_numeric_PASS'] and mr['angle_resolution_PASS'] and fc and fc.get('strict_interior_minimum') and mr['ODE_PASS'] is not False);geos.append(geo)
  out.update({kind+'_fold_count':mr['finite_grid_fold_count'],kind+'_CZROB02_saved_fold_count':op,kind+'_fold_count_change':mr['finite_grid_fold_count']-op,kind+'_geometry_supported':geo,kind+'_ODE_runs':mr['ODE_count'],kind+'_ODE_PASS':mr['ODE_PASS'],kind+'_refined_range_km':fc['range_m']/1000 if fc else None,kind+'_deep_turn_z_m':fc['turn_z_m'] if fc else None,kind+'_refined_angle_deg':fc['angle_deg'] if fc else None,kind+'_full_domain_max_reference_delta_m':max((float(a['main_reference_delta_m']) for a in valid if a['main_reference_delta_m']),default=None)})
  # Stored independent full-domain values independently reproduce coarse-grid metadata, no new integration.
  candidates=[]
  for j in range(2,len(rr)-2,2):
   a,b,d=rr[j-2],rr[j],rr[j+2]
   if all(v['status']=='FINITE' and v['zero_boundary']=='True' for v in (a,b,d)) and float(b['range_m'])<min(float(a['range_m']),float(d['range_m'])) and 10000<=float(b['range_m'])<=120000:candidates.append(b)
  cf=min(candidates,key=lambda a:float(a['range_m'])) if candidates else None;coarse_match=bool(cf and mr['coarse_grid_candidate'] and abs(float(cf['range_m'])-mr['coarse_grid_candidate']['range_m'])<=f['independent_range_tolerance_m'])
  domains.append({'case_id':name,'SSP_method':kind,'registered_angle_rows':len(rr),'finite_reference_rows':len(valid),'max_main_reference_delta_m':max((float(a['main_reference_delta_m']) for a in valid if a['main_reference_delta_m']),default=None),'max_turn_delta_m':max((float(a['turn_delta_m']) for a in valid if a['turn_delta_m']),default=None),'adaptive_warnings':sum(int(a['warnings']) for a in valid),'main_algorithm_counts':json.dumps(dict(collections.Counter(a.get('algorithm','UNRESOLVED') for a in pr))),'coarse_grid_reference_metadata_match':coarse_match,'coarse_first_delta_m':abs(float(cf['range_m'])-mr['coarse_grid_candidate']['range_m']) if cf and mr['coarse_grid_candidate'] else None,'coarse_fine_100m_diagnostic_PASS':mr['angle_resolution_PASS']})
 out['FIRST_RETURN_GEOMETRY_SUPPORTED']=all(geos);summary.append(out)
table(O/'NUM03_NINE_CASE_SENSITIVITY.csv',summary);table(O/'NUM03_FULL_DOMAIN_SUMMARY.csv',domains)
historical=[]
for a in read(R/'results/ENV_CAP_CZROB02/CZROB02_NUMERICAL_FAILURE_ATTRIBUTION.csv'):
 cid=a['case_id'];kind=a['SSP_method'];idx=int(a['angle_index']);b=next(v for v in refs if v['case_id']==cid and v['SSP_method']==kind and int(v['angle_index'])==idx);m=next(v for v in mp if v['case_id']==cid and v['SSP_method']==kind and v['reason']=='PREREGISTERED_CZROB02_FAILED_ANGLE');historical.append({'case_id':cid,'SSP_method':kind,'angle_deg':a['angle_deg'],'old_CZROB02_reference_delta_m':a['independent_reference_delta_m'],'NUM03_main_reference_delta_m':b['main_reference_delta_m'],'NUM03_double_MP80_delta_m':m['double_MP80_delta_m'],'NUM03_MP60_MP80_decimal_delta_m':m['MP60_MP80_delta_m'],'old_failure_unchanged':True})
table(O/'NUM03_HISTORICAL_FAILURE_RECHECK.csv',historical)
valid=[a for a in refs if a['status']=='FINITE'];maxdelta=max(float(a['main_reference_delta_m']) for a in valid if a['main_reference_delta_m']);maxturn=max(float(a['turn_delta_m']) for a in valid if a['turn_delta_m']);maxmp=max(float(a['double_MP80_delta_m']) for a in mp);maxmpconv=max(float(a['MP60_MP80_delta_m']) for a in mp);maxfirst=max(float(a['range_delta_m']) for a in first);maxfirstz=max(float(a['turn_delta_m']) for a in first);maxode=max(float(a['range_difference_m']) for a in ode);maxodeiv=max(float(a['invariant_error']) for a in ode);geo=sum(a['FIRST_RETURN_GEOMETRY_SUPPORTED'] for a in summary);shift=[a['representation_range_shift_m'] for a in summary];zshift=[a['representation_turn_shift_m'] for a in summary];usage=[json.loads((O/(x+'_PROCESS_USAGE.json')).read_text(encoding='utf-8')) for x in ['ANALYSIS','COLD']];wall=sum(a['elapsed_s'] for a in usage);rss=max(a['tree_peak_RSS_bytes'] for a in usage);positive=bool(numeric and cold['PASS'] and geo==9);decision='NUM03_REFERENCE_AND_NOMINAL_FIRST_RETURN_GEOMETRY_SUPPORTED' if positive else 'NUM03_REFERENCE_OR_GEOMETRY_NOT_CLOSED_NUMERICAL_REPAIR_ROUTE_STOPPED'
# Scientific figures are from saved values; no field or receiver data generated.
fig,axes=plt.subplots(3,3,figsize=(13,10),layout='constrained')
for ax,c in zip(axes.flat,p['cases']):
 for kind,color in [('C_LINEAR','black'),('PCHIP','tab:blue')]:
  q=[a for a in rows if a['case_id']==c['case_id'] and a['SSP_method']==kind and a['status']=='FINITE' and a['zero_boundary_deep_return']=='True'];ax.plot([float(a['angle_deg']) for a in q],[float(a['range_m'])/1000 for a in q],c=color,lw=.8,label=kind);fc=c['methods'][kind]['first_refined_candidate']
  if fc:ax.plot(fc['angle_deg'],fc['range_m']/1000,'o',color=color,ms=3)
 ax.set_xlim(.1,20);ax.set_ylim(10,120);ax.set_title(c['case_id'].replace('SCS17_',''));ax.set_xlabel('Launch angle (deg)');ax.set_ylabel('Source-depth return range (km)');ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=8);fig.suptitle('NUM03 nominal zero-boundary deep-return geometry; registered 10–120 km window')
for suf in ['png','pdf']:fig.savefig(O/f'figures/NUM03_return_geometry.{suf}',dpi=160)
plt.close(fig)
fig,axes=plt.subplots(3,3,figsize=(13,10),layout='constrained')
for ax,c in zip(axes.flat,p['cases']):
 for kind,color in [('C_LINEAR','black'),('PCHIP','tab:blue')]:
  q=[a for a in refs if a['case_id']==c['case_id'] and a['SSP_method']==kind and a['status']=='FINITE'];ax.plot([float(a['angle_deg']) for a in q],[max(1e-12,float(a['main_reference_delta_m'])) for a in q],c=color,lw=.8,label=kind)
 ax.axhline(.01,c='tab:red',ls='--',label='Frozen 0.01 m limit');ax.set_yscale('log');ax.set_xlim(.1,20);ax.set_ylim(1e-12,1e-1);ax.set_title(c['case_id'].replace('SCS17_',''));ax.set_xlabel('Launch angle (deg)');ax.set_ylabel('Independent range difference (m)');ax.grid(alpha=.2)
axes.flat[0].legend(fontsize=7);fig.suptitle('All registered angular entries independently checked; display floor 1e-12 m')
for suf in ['png','pdf']:fig.savefig(O/f'figures/NUM03_numerical_closure.{suf}',dpi=160)
plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained');xx=np.arange(9);labels=[a['case_id'].replace('SCS17_','').replace('_z','/') for a in summary];axes[0].bar(xx,shift,color='tab:blue');axes[0].set_ylabel('C-linear / PCHIP first-fold shift (m)');axes[1].bar(xx-.18,[a['C_LINEAR_fold_count'] for a in summary],.36,label='C-linear');axes[1].bar(xx+.18,[a['PCHIP_fold_count'] for a in summary],.36,label='PCHIP');axes[1].set_ylabel('Registered finite-grid local minima');axes[1].legend()
for ax in axes:ax.set_xticks(xx,labels,rotation=45,ha='right');ax.grid(axis='y',alpha=.2)
fig.suptitle('Two SSP representations; one climatological source family')
for suf in ['png','pdf']:fig.savefig(O/f'figures/NUM03_representation_sensitivity.{suf}',dpi=160)
plt.close(fig)
table(O/'FIGURE_DATA_PROVENANCE.csv',[{'figures':'figures/NUM03_return_geometry.png/.pdf','input':'NUM03_RETURN_CURVES.csv; NUM03_CANDIDATE_CATALOGUE.csv','sha256':' ; '.join(h(O/x) for x in ['NUM03_RETURN_CURVES.csv','NUM03_CANDIDATE_CATALOGUE.csv']),'scope':'Nominal zero-boundary branch; receiver field not evaluated'},{'figures':'figures/NUM03_numerical_closure.png/.pdf','input':'NUM03_FULL_DOMAIN_INDEPENDENT_REFERENCE.csv','sha256':h(O/'NUM03_FULL_DOMAIN_INDEPENDENT_REFERENCE.csv'),'scope':'All72018 registered rows; display floor only'},{'figures':'figures/NUM03_representation_sensitivity.png/.pdf','input':'NUM03_NINE_CASE_SENSITIVITY.csv','sha256':h(O/'NUM03_NINE_CASE_SENSITIVITY.csv'),'scope':'Two representations, one source family; finite-grid topology'}])
labels_status={'NUMERICAL_REFERENCE_VALIDATED':numeric,'FIRST_RETURN_GEOMETRY_SUPPORTED':geo==9,'REPRESENTATION_SENSITIVITY_CHARACTERIZED':bool(numeric and len(summary)==9 and len(first)==18 and all(a['PASS']=='True' for a in first)),'ENERGY_ATTRIBUTION_UNRESOLVED':True,'YANG_MMAR_PERFORMANCE_NOT_EVALUATED':True}
final={'stage':'ENV_CAP_NUM03','parent_SHA':f['parent_SHA'],'design_SHA':'a641df8e72adf7acec4a71e295c06c78b7458079','branch':f['branch'],'status':'STOP','decision':decision,'labels':labels_status,'cold_audit':{'PASS':cold['PASS'],'checks':cold['checks'],'FAIL_count':cold['FAIL_count']},'geometrically_supported_cases':geo,'input_cases':9,'source_families':1,'scope':'Same WOA23_17N118E climatological nodes, C-linear/PCHIP, fixed depth and ideal boundaries; source-depth-plane ray return; not reception or full-wave accuracy','registered_angular_entries':len(rows),'precision_evidence':{'max_primary_independent_range_delta_m':maxdelta,'max_fixed_angle_turn_delta_m':maxturn,'MP_controls':len(mp),'max_independent_double_MP80_delta_m':maxmp,'max_MP60_MP80_decimal_delta_m':maxmpconv,'max_first_candidate_refinement_range_delta_m':maxfirst,'max_first_candidate_refinement_turn_delta_m':maxfirstz,'ODE_runs':len(ode),'max_ODE_Snell_delta_m':maxode,'max_ODE_invariant_error':maxodeiv},'representation_range_shift_m':[min(shift),max(shift)],'representation_turn_shift_m':[min(zshift),max(zshift)],'resource':{'science_wall_s':wall,'sampled_peak_RSS_bytes':rss,'new_acoustic_solver_calls':0,'new_MC':0,'new_receiver_or_SHD':0},'next_recommendation':'DESIGN_ONE_BOUNDED_PROPAGATION_OBSERVABLE_STUDY_ONLY' if positive else 'STOP_THIS_NUMERICAL_REPAIR_DIRECTION_AND_DOCUMENT_APPLICATION_CONDITIONS','automatic_coherent_provider_development':False,'automatic_Yang_MMAR_performance':False,'historical_R4':'0%','all_historical_Gates_STOP_main':'UNCHANGED'}
js(O/'NUM03_FINAL_DECISION.json',final)
lines=['|单元|C-linear首折叠 km|PCHIP首折叠 km|位置差 m|C/P登记极小值数|','|---|---:|---:|---:|---:|']
for a in summary:lines.append(f"|{a['case_id']}|{a['C_LINEAR_refined_range_km']:.6f}|{a['PCHIP_refined_range_km']:.6f}|{a['representation_range_shift_m']:.3f}|{a['C_LINEAR_fold_count']}/{a['PCHIP_fold_count']}|")
text(O/'NUM03_RESEARCH_DECISION.md',f'''# NUM03：最终科研去留

正式结论：**{decision}**。冷审计PASS={cold['PASS']}，{cold['checks']:,}项、{cold['FAIL_count']}失败。本轮完成后STOP；不追加数值修补轮次，不自动建设相干提供器。

## 1. 数值参考闭合到什么程度

合成开发有效数值检查16,113项通过；测试驱动v002的解析失败已保留，没有隐藏或算作科学结果。积分器与独立参考从首版v001起未修改。实际SCS科学在Commit A之后才开始，主分析和冷审计各一次。

完整登记角网格9单元×2表示×4001角，共72,018条，全部由独立斜率、根和局部平方根适应积分重建。主冷回返距离最大差{maxdelta:.6g} m，固定角转向深度最大差{maxturn:.6g} m；预冻结要求分别0.01 m、0.0001 m。30组高精度对照从节点重新建立60/80位参考，独立双精度与MP80最大差{maxmp:.6g} m，60/80位十进制差最大{maxmpconv:.6g} m。独立首候选细化距离最大差{maxfirst:.6g} m，转向位置最大差{maxfirstz:.6g} m。

登记代表案例执行{len(ode)}次DOP853轨迹，最大ODE—Snell距离差{maxode:.6g} m，Snell不变量相对误差最大{maxodeiv:.6g}。非代表单元没有假称也执行了ODE。全部轨迹、步长对照和原始索引保留。

NUMERICAL_REFERENCE_VALIDATED={numeric}，只覆盖这些登记条目、规定的离网格候选与高精度控制，不是所有连续角、任意SSP或严格区间余项的数学证明。Gauss吻合、QUADPACK误差估计及高精度一致性均按各自作用域使用，不能把零求积估计理解成全部误差为零。

## 2. 首深折返与SSP表示

FIRST_RETURN_GEOMETRY_SUPPORTED覆盖{geo}/9个名义单元。这里的距离是**射线返回源深平面的首深折叠距离**，不是定位估计误差，也不是新接收阵获得的测距。表中转向深度同样不是目标深度估计。

{chr(10).join(lines)}

两表示的首细化候选移动{min(shift):.3f}–{max(shift):.3f} m，深转向位置移动{min(zshift):.3f}–{max(zshift):.3f} m。这支持名义几何中的首候选在本次C1对照下保留，不能简单解释成C-linear梯度跳变制造了全部首折叠。它仍不证明真实海洋结构；PCHIP只连续一阶导数，不是实测光滑真值。

REPRESENTATION_SENSITIVITY_CHARACTERIZED={labels_status['REPRESENTATION_SENSITIVITY_CHARACTERIZED']}。主/冷独立全域曲线都保留全部登记极小值与括区细化。部分次级极小值与旧CZROB-02计数不同；新旧目录并列，旧结论不追溯改判。两表示的候选数量也不同，因此不能用首位置相近抹平多分支敏感性。没有声称离散点间所有潜在极小值已被严格排除。

## 3. 能量和可提取性分别记账

本轮没有新BELLHOP、SHD或接收波形。原REC01支持C-linear条件下的总TL增强，旧RAY支持不同路径族几何共存，但没有各族能量分配。**ENERGY_ATTRIBUTION_UNRESOLVED**继续成立，它不再被当作折射几何存在的自动否决条件。

**YANG_MMAR_PERFORMANCE_NOT_EVALUATED**继续成立。没有新相干模态、已提取复幅度、测深误差或未知源/校准后的性能结果。Snell数值精度不是200 Hz全波相位准确性，亦不代表模态完整性；T2/FUNC2保持关闭。

## 4. 去留与项目申请

{'本轮可收束首深回返积分的数值修补专项。只建议由负责人另行设计一次有界传播可观测量研究，检验折射和边界路径共存时哪些量能够被接收处理提取；不直接开放相干提供器开发或大规模Yang/MMAR扫描。' if positive else '本轮未闭合的标签保持未认证。按照最后一次专项约定，停止这条数值修补方向，使用现有条件证据开展申请与资源边界论证；不另起普通NUM/CZROB轮次。'}

申请中可用的是：同一WOA23_17N118E气候来源族下，原节点、平底、理想边界和两种表示中的条件几何/敏感性证据；不是9个独立海域、真实探测覆盖或工程鲁棒性。研究问题可以聚焦环境不确定性下的多参数可辨识性、可提取传播信息和观测资源需求边界，而不要求先承诺四维5%/10%或五维同精度。

科学壁钟{wall:.2f} s，采样RSS峰值{rss/1024**2:.1f} MiB；新增声学求解0、MC0，所有新增材料在D盘。main及所有历史Gate、STOP与R4=0%不变。提交推送后STOP。
''')
text(O/'GPT_SYNC_NUM03.md',f'''# ENV-CAP-NUM03：最终同步

**{decision}**。本轮完成并恢复STOP。

- 父SHA：`{f['parent_SHA']}`。
- 设计SHA：`{final['design_SHA']}`。
- 最终SHA：本报告所在Commit B；完整身份见远端`{f['branch']}`。

## 分层结论

|标签|结果与作用域|
|---|---|
|NUMERICAL_REFERENCE_VALIDATED|{numeric}；全部72,018登记条目独立复核，最大距离差{maxdelta:.3g} m|
|FIRST_RETURN_GEOMETRY_SUPPORTED|{geo}/9名义单元；返回源深平面的折射回返与折叠|
|REPRESENTATION_SENSITIVITY_CHARACTERIZED|{labels_status['REPRESENTATION_SENSITIVITY_CHARACTERIZED']}；C-linear/PCHIP首候选移动{min(shift):.2f}–{max(shift):.2f} m，次级候选也保留|
|ENERGY_ATTRIBUTION_UNRESOLVED|保持；路径共存不自动否定几何，也未量化贡献|
|YANG_MMAR_PERFORMANCE_NOT_EVALUATED|保持；无模态提取、测深或接收性能新结论|

冷审计**{cold['checks']:,}项，{cold['FAIL_count']}失败，PASS={cold['PASS']}**。另有30组60/80位参考，独立双精度与MP80最大距离差{maxmp:.3g} m；{len(ode)}次登记DOP853对照，最大与Snell距离差{maxode:.3g} m。合成开发16,113项有效数值控制通过；开发测试脚本的一次语法失败和修改原因完整保留，未用9个SCS结果调试。

这次对旧失败角使用了新身份、新方法和新冻结合同，未调用或修改旧科学程序；旧6项失败保持原样，新对照见NUM03_HISTORICAL_FAILURE_RECHECK.csv。数值参考、首候选精度、拓扑和表示差异分别交付，不用一个PASS替代。

全部9单元仍只属于一个WOA气候资料来源族。新几何证据不能扩展为真实海域、所有SSP、全波传播、目标测深或工程环境鲁棒性。{'建议下一步只设计一次有界的传播可观测量研究；当前数值修补专项收口，不自动建设相干提供器。' if positive else '本方向停止继续数值修补，转入申请条件与资源边界论证。'}

主分析与冷审计各一次；累计科学壁钟{wall:.2f} s，采样RSS峰值{rss/1024**2:.1f} MiB。新增外部声学求解／MC／接收SHD均0。main、ENV-CAP、REC01、CZROB-02、T2/FUNC2及R4=0%不变。

详细见[科研去留报告](NUM03_RESEARCH_DECISION.md)、[9单元敏感性表](NUM03_NINE_CASE_SENSITIVITY.csv)、[数学推导](ANALYTIC_REFERENCE_DERIVATION.md)及[独立冷审计](NUM03_INDEPENDENT_COLD_AUDIT.json)。
''')
text(O/'PROJECT_APPLICATION_CONDITIONAL_CLAIMS.md','''# 申请论证：三类条件证据分别使用

|证据|可支持的表述|仍不能支持|
|---|---|---|
|REC01射线与总TL、NUM03几何数值参考|同一气候来源族、名义边界下的首深折返几何与条件传播环境；数值实现和敏感性可复核|真实海域检测距离、现场能量占比、接收可提取性|
|CZROB-02及NUM03两种表示对照|首候选位置与多分支结构对环境表示的影响是明确研究问题；旧失败与新证据分开|全部连续SSP变化的鲁棒性、全波相位或模态认证|
|历史T3-W1双节点水平Demo|有符号方位及登记导航/偏置条件下的资源—水平参数信息增量|实际单拖阵左右舷解除、目标深度或五维精度|

拟研究的创新问题可表述为：第一会聚区环境不确定性下，哪些传播结构能转化为可提取观测；这些观测在未知源、校准和宽水平候选支持下如何改变多参数可辨识性，并需要怎样的环境先验、相干孔径与协同资源。

这是研究计划和条件基础，不是宣称已经完成全部创新方法。数值参考闭合与传播信息可提取性必须分开；不把历史四参数5%/10%研究筛选指标写成普遍验收标准。硬件、实际同步与目标关联仍需独立证据，R4不因本轮几何研究增加进度。

来源：本目录NUM03_RESEARCH_DECISION.md；../ENV_CAP_REC01/GPT_SYNC_ENV_CAP_REC01.md；../ENV_CAP_CZROB02/CZROB02_RESEARCH_DECISION.md；../R4_T3_W1_WIDE_BASELINE_DEMO/GPT_SYNC_T3_W1.md。历史结果均只读，本轮不重算T3。
''')
js(O/'STOP.json',{'stage':'ENV_CAP_NUM03','status':'STOP','decision':decision,'labels':labels_status,'no_more_automatic_numerical_repair':True,'new_solver_calls':0,'new_MC':0,'no_automatic_coherent_or_performance_work':True,'historical_R4':'0%','next_stage_requires_new_authorization':True,'time_unix':time.time()})
js(O/'RUN_STATE.json',{'stage':'NUM03','phase':'STOP','formal_runs':1,'cold_runs':1,'development_versions':['v001','v002_parse_failure','v003','v004'],'frozen_code_modified_after_A':False,'science_source_data_reused_not_regenerated':True})
text(O/'REVIEW_ENTRYPOINT.md','''# NUM03审核入口

先读GPT_SYNC_NUM03.md及NUM03_RESEARCH_DECISION.md。分开核查数值参考、首回返几何、两表示敏感性、未解能量和未评价性能。开发失败记录不删除，旧CZROB-02失败不改判。

完整全角域主/冷曲线、30组60/80位控制、全部候选与首候选独立细化、DOP853轨迹和哈希索引均交付。图的数据和范围见FIGURE_DATA_PROVENANCE.csv。所有距离精度首先是积分比较，转向深度不是目标深度估计。完整连续角域无遗漏及真实海洋适用性没有由此证明。

最终身份由设计Commit A和本报告所在Commit B锁定；远端核对后STOP，不合并main、不自动启动下一阶段。
''')
idx=[{'scope':a['source'],'path':a['path'],'bytes':Path(a['path']).stat().st_size,'sha256':a['sha256']} for a in f['execution_inputs']]
for q in sorted(W.iterdir()):
 if q.is_file() and q.suffix in ['.json','.npz','.log','.py']:
  idx.append({'scope':'NUM03_local_execution_or_operation','path':str(q),'bytes':q.stat().st_size,'sha256':h(q)})
  if q.suffix in ['.json','.log']:(O/'inputs'/q.name).write_bytes(q.read_bytes())
  elif q.suffix=='.py':(O/'inputs'/('OPERATIONAL_'+q.name)).write_bytes(q.read_bytes())
table(O/'RAW_DATA_INDEX.csv',idx)
done.set();js(O/'DOCUMENTATION_PROCESS_USAGE.json',{'purpose':'Saved-result summaries/rendering only; no scientific execution','wall_s':time.time()-start,'sampled_peak_RSS_bytes':peak[0]})
print(json.dumps({'decision':decision,'cold_PASS':cold['PASS'],'cold_checks':cold['checks'],'cold_FAIL':cold['FAIL_count'],'max_full_reference_delta_m':maxdelta,'max_MP_delta_m':maxmp,'geometry_cases':geo,'first_candidate_shift_m':[min(shift),max(shift)],'ODE_runs':len(ode),'science_wall_s':wall,'science_peak_MiB':rss/1024**2,'documentation_peak_MiB':peak[0]/1024**2}))
