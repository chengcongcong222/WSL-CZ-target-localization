"""T3-W1 postprocessing only: read sealed outputs, plot distributions, close scope.
No observations, fits, Fisher matrices or scientific tests are rerun.
"""
from pathlib import Path
import csv,json,hashlib,math,shutil,subprocess,datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
ROOT=Path.cwd();OUT=ROOT/'results/R4_T3_W1_WIDE_BASELINE_DEMO';LOCAL=Path('D:/ProjectStorage/WSL-CZ/T3_W1')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(name,data):
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:
  fields=list(dict.fromkeys(k for r in data for k in r));w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(data)
def git(*a):return subprocess.check_output(['git','-c','safe.directory='+ROOT.as_posix(),*a],cwd=ROOT).decode().strip()
def number(x,percent=False):
 if x is None:return 'NA'
 x=float(x)
 if not math.isfinite(x):return 'INF'
 return f'{x*(100 if percent else 1):.4f}'

def main():
 cfg=read(OUT/'T3_W1_DESIGN_FREEZE.json');audit=read(OUT/'INDEPENDENT_COLD_AUDIT.json');budget=read(OUT/'EXECUTION_BUDGET.json');raw_decision=read(OUT/'T3_W1_DECISION.json');seal=read(OUT/'MAIN_OUTPUT_SEAL.json');original=(LOCAL/'LEDGER_BEFORE.bin').read_bytes();ledger=ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'
 for item in seal['files']:assert sha(item['path'])==item['sha256'],item['path']
 for item in cfg['inputs']:assert sha(ROOT/item['path'])==item['sha256'],item['path']
 for item in cfg['frozen_code']:assert sha(OUT/item['path'])==item['sha256'],item['path']
 assert budget['main_calls']==budget['cold_calls']==1 and budget['scientific_retries']==0
 assert ledger.read_bytes()==original
 summary=rows(OUT/'GLOBAL_STATISTICS.csv');cases=rows(OUT/'CASE_STATISTICS.csv');infos=rows(OUT/'LOCAL_INFORMATION.csv');geom=rows(OUT/'GEOMETRY_AND_RESOURCE_AUDIT.csv');runs=rows(OUT/'RUN_ERRORS.csv');estimates=rows(OUT/'ESTIMATOR_OUTPUTS.csv');pairs=rows(OUT/'PAIRED_CHANGES.csv');cold_ok=audit['PASS']
 report_rows=[];info_rows=[];increments=[];directions=[];resources=[]
 for bname,scheme in cfg['geometries'].items():
  gg=[r for r in geom if r['quality']=='Q1' and r['geometry']==bname];resources.append({'geometry':bname,'forward_km':scheme['forward_km'],'lateral_km':scheme['lateral_km'],'node_separation_km':math.hypot(scheme['forward_km'],scheme['lateral_km']),'relative_to5km_resource':math.hypot(scheme['forward_km'],scheme['lateral_km'])/5,'aux_range_min_km':min(float(r['aux_target_min_m']) for r in gg)/1000,'aux_range_max_km':max(float(r['aux_target_max_m']) for r in gg)/1000,'outside45to60km_epoch_locations':sum(int(r['aux_outside45to60km_epochs']) for r in gg),'registered_locations':2904,'nodes':2,'measurements_per_record':242,'hardware_DOA_navigation_sync_association':'UNKNOWN_OR_NOT_VALIDATED','cost_scope':'geometricdeploymentresourceonly;nomonetaryclaim;inrangeisnotdetectability'})
 for qname in cfg['qualities']:
  for bname in cfg['geometries']:
   ss=next(r for r in summary if r['config']==qname+'_'+bname);cc=[r for r in cases if r['config']==qname+'_'+bname and r['side_group']=='BOTH'];ii=[r for r in infos if r['quality']==qname and r['geometry']==bname];resource=next(r for r in resources if r['geometry']==bname)
   report_rows.append({'quality':qname,'geometry':bname,'n':int(ss['n']),'accepted':int(ss['accepted']),'failure_rate':float(ss['failure_rate']),'range_worst_case_P95_percent':100*float(ss['range_rel_worst_case_unconditional_P95']),'bearing_worst_case_P95_deg':float(ss['bearing_deg_worst_case_unconditional_P95']),'speed_abs_worst_case_P95_mps':float(ss['speed_abs_mps_worst_case_unconditional_P95']),'speed_worst_case_P95_percent':100*float(ss['speed_rel_worst_case_unconditional_P95']),'heading_worst_case_P95_deg':float(ss['heading_deg_worst_case_unconditional_P95']),'speed_pooled_median_percent':100*float(ss['speed_rel_unconditional_median']),'speed_pooled_P99_percent':100*float(ss['speed_rel_unconditional_P99']),'max_case_bias_boundary_rate':max(float(r['bias_boundary_rate']) for r in cc),'max_case_state_boundary_rate':max(float(r['state_boundary_rate']) for r in cc),'aux_range_min_km':resource['aux_range_min_km'],'aux_range_max_km':resource['aux_range_max_km']})
   free=[r for r in ii if r['control'].startswith('BIAS_FREE')];known=[r for r in ii if r['control'].startswith('BIAS_KNOWN')];info_rows.append({'quality':qname,'geometry':bname,'free_bias_worst_1p96sigma_speed_percent':100*max(float(r['speed_rel_1p96sigma']) for r in free),'known_bias_worst_1p96sigma_speed_percent':100*max(float(r['speed_rel_1p96sigma']) for r in known),'min_free_rank':min(int(r['rank']) for r in free),'max_scaled_design_condition':max(float(r['condition_design']) for r in free),'scope':'localGaussianequivalentnotempiricalP95;uniformsupportnotFisherprior'})
  for case in sorted({r['case_id'] for r in cases}):
   base=next(r for r in cases if r['config']==qname+'_B0' and r['case_id']==case and r['side_group']=='BOTH');candidate=next(r for r in cases if r['config']==qname+'_B3' and r['case_id']==case and r['side_group']=='BOTH');directions.append({'quality':qname,'case_id':case,'median_speed_improved':float(candidate['speed_rel_unconditional_median'])<float(base['speed_rel_unconditional_median']),'P95_speed_improved':float(candidate['speed_rel_unconditional_P95'])<float(base['speed_rel_unconditional_P95']),'P95_range_improved':float(candidate['range_rel_unconditional_P95'])<float(base['range_rel_unconditional_P95']),'P95_heading_improved':float(candidate['heading_deg_unconditional_P95'])<float(base['heading_deg_unconditional_P95']),'speed_median_ratio':float(candidate['speed_rel_unconditional_median'])/float(base['speed_rel_unconditional_median']),'speed_P95_ratio':float(candidate['speed_rel_unconditional_P95'])/float(base['speed_rel_unconditional_P95'])})
  for bname in ['B1','B2','B3']:
   for case in sorted({r['case_id'] for r in infos}):
    for side in ['-1','1']:
     for control in ['BIAS_FREE_NO_PRIOR_LOCAL_DATA_INFORMATION','BIAS_KNOWN_INFORMATION_UPPER_CONTROL']:
      b=next(r for r in infos if r['quality']==qname and r['geometry']=='B0' and r['case_id']==case and r['side']==side and r['control']==control);c=next(r for r in infos if r['quality']==qname and r['geometry']==bname and r['case_id']==case and r['side']==side and r['control']==control);bs,cs=float(b['sigma_speed_mps']),float(c['sigma_speed_mps']);increments.append({'quality':qname,'geometry':bname,'case_id':case,'side':side,'control':control,'baseline_sigma_speed_mps':bs,'candidate_sigma_speed_mps':cs,'speed_std_ratio':cs/bs if math.isfinite(bs) and bs>0 else None,'variance_reduction':1-(cs/bs)**2 if math.isfinite(bs) and bs>0 else None,'direction_improved':cs<bs})
 table('FOUR_PARAMETER_DEMO_TABLE.csv',report_rows);table('INFORMATION_SUMMARY.csv',info_rows);table('INFORMATION_INCREMENT_BY_SCENE.csv',increments);table('B3_DIRECTION_BY_CASE.csv',directions);table('RESOURCE_AND_COVERAGE_SUMMARY.csv',resources)
 # Scientific figures use sealed outputs only; no data selection or threshold fitting.
 colors=['#475569','#2563eb','#d97706','#be123c'];separation=[r['node_separation_km'] for r in resources]
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for qname,style in [('Q1','-o'),('Q2','--s')]:
  rr=[next(r for r in report_rows if r['quality']==qname and r['geometry']==b) for b in cfg['geometries']];axes[0].plot(separation,[r['speed_worst_case_P95_percent'] for r in rr],style,label=qname);axes[1].plot(separation,[r['range_worst_case_P95_percent'] for r in rr],style,label=qname)
 for ax in axes:ax.set_xlabel('Node separation (km)');ax.grid(alpha=.25);ax.legend()
 axes[0].set_ylabel('Worst registered-case speed P95 error (%)');axes[1].set_ylabel('Worst registered-case range P95 error (%)');fig.suptitle('Same two-node resource count; full paired test errors, failures = infinity');fig.tight_layout();fig.savefig(OUT/'FIG1_RESOURCE_ERROR.png',dpi=170);fig.savefig(OUT/'FIG1_RESOURCE_ERROR.pdf');plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(10,4),sharey=True)
 for ax,qname in zip(axes,['Q1','Q2']):
  for b,c in zip(cfg['geometries'],colors):
   data=np.array([float(r['speed_rel']) for r in runs if r['config']==qname+'_'+b]);a=np.sort(data[np.isfinite(data)])*100;ax.step(a,np.arange(1,len(a)+1)/len(data),where='post',label=f'{b}; failures {len(data)-len(a)}',color=c)
  ax.set_title(qname);ax.set_xlabel('Absolute relative speed error (%)');ax.grid(alpha=.25);ax.legend(fontsize=8)
 axes[0].set_ylabel('Empirical CDF (denominator includes failures)');fig.suptitle('All6000records per configuration; no tail trimming');fig.tight_layout();fig.savefig(OUT/'FIG2_SPEED_DISTRIBUTIONS.png',dpi=170);fig.savefig(OUT/'FIG2_SPEED_DISTRIBUTIONS.pdf');plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(10,4),sharey=True)
 for ax,qname in zip(axes,['Q1','Q2']):
  rr=[next(r for r in info_rows if r['quality']==qname and r['geometry']==b) for b in cfg['geometries']];ax.plot(separation,[r['free_bias_worst_1p96sigma_speed_percent'] for r in rr],'-o',label='Bias unknown; no prior');ax.plot(separation,[r['known_bias_worst_1p96sigma_speed_percent'] for r in rr],'--s',label='Bias known: information control');ax.set_title(qname);ax.set_xlabel('Node separation (km)');ax.grid(alpha=.25);ax.legend(fontsize=8)
 axes[0].set_ylabel('Worst local 1.96-sigma speed equivalent (%)');fig.suptitle('Local Fisher diagnostics; not empirical P95 or coverage certificate');fig.tight_layout();fig.savefig(OUT/'FIG3_BIAS_INFORMATION_LIMIT.png',dpi=170);fig.savefig(OUT/'FIG3_BIAS_INFORMATION_LIMIT.pdf');plt.close(fig)
 mat=np.load(OUT/'INFORMATION_MATRICES.npz');fig,axes=plt.subplots(1,2,figsize=(10,4),sharex=True,sharey=True)
 for ax,qname in zip(axes,['Q1','Q2']):
  extent=0
  for b,color in zip(cfg['geometries'],colors):
   row=next(r for r in infos if r['quality']==qname and r['geometry']==b and r['case_id']=='H01' and r['side']=='1' and r['control'].startswith('BIAS_FREE'));E=mat['effective_velocity_scaled'][int(row['matrix_index'])];ev,U=np.linalg.eigh(E)
   if np.all(ev>0):
    widths=2/np.sqrt(ev)*2;angle=math.degrees(math.atan2(U[1,0],U[0,0]));patch=Ellipse((0,0),widths[0],widths[1],angle=angle,fill=False,color=color,label=b);ax.add_patch(patch);extent=max(extent,float(max(widths)/2))
  ax.set_title(qname+'; H01; positive side');ax.set_aspect('equal');ax.set_xlim(-extent*1.08,extent*1.08);ax.set_ylim(-extent*1.08,extent*1.08);ax.set_xlabel('Delta vx (m/s)');ax.grid(alpha=.25);ax.legend(fontsize=8)
 axes[0].set_ylabel('Delta vy (m/s)');fig.suptitle('Profiled local1-sigma velocity ellipses; unknown position and bias');fig.tight_layout();fig.savefig(OUT/'FIG4_VELOCITY_WEAK_DIRECTION.png',dpi=170);fig.savefig(OUT/'FIG4_VELOCITY_WEAK_DIRECTION.pdf');plt.close(fig);mat.close()
 figures=[{'file':p.name,'source':'sealed estimator statistics or nominal Fisher matrices; no new computation of observations/estimates','sha256':sha(p)} for p in OUT.glob('FIG*')];write(OUT/'FIGURE_INDEX.json',figures)
 consistency={q:{'median_speed_cases':sum(r['median_speed_improved'] for r in directions if r['quality']==q),'P95_speed_cases':sum(r['P95_speed_improved'] for r in directions if r['quality']==q),'P95_range_cases':sum(r['P95_range_improved'] for r in directions if r['quality']==q),'P95_heading_cases':sum(r['P95_heading_improved'] for r in directions if r['quality']==q),'Fisher_free_side_cases_improved':sum(r['direction_improved'] for r in increments if r['quality']==q and r['geometry']=='B3' and r['control'].startswith('BIAS_FREE'))} for q in ['Q1','Q2']}
 final={'stage':'R4-T3-W1','parent_SHA':cfg['parent_SHA'],'design_SHA':read(OUT/'EXECUTION_STARTED.json')['design_SHA'],'research_result':'CONDITIONAL_WIDE_BASELINE_HORIZONTAL_INFORMATION_AND_DEMO' if cold_ok else 'NOT_ADMITTED_INDEPENDENT_AUDIT_FAILED','claim_level':'RESOURCE_INFORMATION_DESIGN_REFERENCE_ONLY;NO_RECEIVED_MEASUREMENT_CHAIN_CERTIFICATION','consistency':consistency,'primary_geometry':'B3=(3,+/-15)km','monte_carlo_configurations':48000,'independent_base_innovations':6000,'cases':12,'known_development_geometries':True,'two_qualities':True,'cold_PASS':cold_ok,'cold_checks':audit['checks'],'cold_FAIL':audit['FAIL'],'bias_prior':'uniformsupportonly;noaddedGaussianpenalty','novelty':'increasingbaselinealoneisstandardgeometricinformationgain;noalgorithmicnoveltyclaim','route':'STOP_WIDE_BASELINE_PARAMETER_SWEEPS;RETAIN_CONDITIONAL_DEMO_AND_MEASUREMENT_CHAIN_GAP','physical_bearing_navigation_hardware':'UNKNOWN','synchronization_association':'ASSUMED_NOT_VALIDATED','depth':'NOT_EVALUATED','T2_STOP':'UNCHANGED','historical_Gates':'UNCHANGED','R4_percent':0,'research_lead_audit':'PENDING','automatic_next_stage':False,'STOP':True};write(OUT/'FINAL_ADJUDICATION.json',final)
 title='''# T3-W1：科研去留与适用条件

当前阶段完成一次宽基线几何资源—信息Demo，不继续扫描更远布放、不更换估计器寻找增益。独立冷审计通过与否见FINAL_ADJUDICATION；它不把合成特征变成真实接收链能力。

## 机制与新意

两个节点保持同一运动、同一观测数与时长。改变基线使两节点视线导数的差异增强，减轻距离与沿视线速度弱方向的耦合；该解释来自冻结Jacobian和剖面速度信息，不依赖高精度会聚区复场。此次没有验证辅助节点独立Z形或其他新航迹。

增大基线提高三角定位精度本身是标准几何效果，不能作为新方法发明。具有研究价值的问题是：远距离、有限时长、固定系统偏置和有限被动资源同时存在时，哪些弱方向得到改善、哪部分收益依赖校准支持，以及部署/覆盖代价是否可接受。当前证据为后续科学问题定标，不证明创新性已经完成。

## 误差与先验

主估计器显式估计共模和半差，按登记Uniform支持作有界拟合，没有真实偏置输入，也没有额外高斯校准正则。偏置触界率必须与误差表并读，触界的干扰状态不视为已校准。四个物理状态同样继承有界域，边界触碰不事后剔除，也不把有限终点叫作完整候选域。

完整局部Fisher包括方位均值和一阶导航等效方差的导数，对位置和无先验偏置作剖面；已知偏置只是信息上限对照。Uniform边界不是正则Gaussian先验，不能将其方差倒数直接加到Fisher。由有限支持和可能有偏的估计得到的MC误差不受局部无偏CRLB直接约束，局部1.96σ也不是经验P95。

## 资源与物理合同

B3接收节点间距约15.30 km，是5 km的约3.06倍；节点数与242个方位样本不变。它要求维持较宽协同几何及相应通信、导航、同步和同目标关联。45–60 km内的几何距离不能证明第一会聚区探测覆盖，也不能证明每个历元可提供0.05/0.075度的有符号方位。

实际测向、导航及阵列能力继续UNKNOWN。申请阶段允许条件研究，不要求先完成设备验收；但本轮不能用仿真假设替代接收链证据，也不直接承诺工程精度。

## 去留

保留该条件水平Demo及完整原始观测，用于说明资源—信息的弱方向与限制；在缺少可靠测量链证据时，整体按资源—信息设计参照定位。停止扩大同类布放仿真，不自动开启新的宽基线参数扫描。

若负责人后续决定继续，下一问题应围绕已有接收资源能否在相应距离持续提供登记质量的有符号方位、固定偏置支持是否有依据、关联与同步条件是否成立，而不是继续寻找更远或更有利的坐标。这些是研究选项，本轮没有执行授权。

没有新增目标深度观测，深度NOT_EVALUATED，完整五维联合定位未证明。原T2数值修补链保持关闭，历史Gate和R4=0%不变。提交推送后STOP。
''';(OUT/'RESEARCH_ROUTE_AND_LIMITS.md').write_text(title,encoding='utf-8')
 tablelines=['|质量|几何|距离最坏场景P95 %|方位 °|速度最坏场景P95 %|速度绝对误差 m/s|航向 °|输出/总数|','|---|---|---:|---:|---:|---:|---:|---:|']
 for r in report_rows:tablelines.append('|'+ '|'.join([r['quality'],r['geometry'],number(r['range_worst_case_P95_percent']),number(r['bearing_worst_case_P95_deg']),number(r['speed_worst_case_P95_percent']),number(r['speed_abs_worst_case_P95_mps']),number(r['heading_worst_case_P95_deg']),f"{r['accepted']}/{r['n']}"])+'|')
 infotable=['|质量|几何|未知偏置、无先验的最坏局部1.96σ速度 %|已知偏置信息对照 %|','|---|---|---:|---:|']
 for r in info_rows:infotable.append(f"|{r['quality']}|{r['geometry']}|{number(r['free_bias_worst_1p96sigma_speed_percent'])}|{number(r['known_bias_worst_1p96sigma_speed_percent'])}|")
 resource_table=['|几何|节点间距 km|辅助—目标距离 km|超出45–60 km的时空点/2904|','|---|---:|---|---:|']
 for r in resources:resource_table.append(f"|{r['geometry']}|{r['node_separation_km']:.3f}|{r['aux_range_min_km']:.3f}–{r['aux_range_max_km']:.3f}|{r['outside45to60km_epoch_locations']}/2904|")
 benefits=[]
 for q in ['Q1','Q2']:
  c=consistency[q];base=next(r for r in report_rows if r['quality']==q and r['geometry']=='B0');wide=next(r for r in report_rows if r['quality']==q and r['geometry']=='B3');benefits.append(f"{q}：B3相对B0，速度中位改善{c['median_speed_cases']}/12，速度P95改善{c['P95_speed_cases']}/12，距离P95改善{c['P95_range_cases']}/12，航向P95改善{c['P95_heading_cases']}/12；无先验偏置Fisher的速度1σ在{c['Fisher_free_side_cases_improved']}/24个状态—侧别条件中下降。最坏速度P95由{base['speed_worst_case_P95_percent']:.4f}%变为{wide['speed_worst_case_P95_percent']:.4f}%。")
 pairlines=[]
 for q in ['Q1','Q2']:
  for metric in ['range_rel','bearing_deg','speed_abs_mps','speed_rel','heading_deg']:
   pp=[r for r in pairs if r['quality']==q and r['geometry']=='B3' and r['metric']==metric];pairlines.append(f"{q}/{metric}：改善{sum(int(r['better']) for r in pp)}，变差{sum(int(r['worse']) for r in pp)}，相等{sum(int(r['equal']) for r in pp)}；分母6000。")
 sync=f'''# R4-T3-W1 宽基线协同观测Demo同步

父提交：`{cfg['parent_SHA']}`；设计SHA：`{final['design_SHA']}`。执行SHA与remote/main由最终回传给出，避免自引用。

## 阶段结论

本轮完成{raw_decision.get('configurations',0)}个配置估计，独立冷审计PASS={cold_ok}，{audit['checks']}项检查、{audit['FAIL']}失败。结论限于冻结的双节点有符号方位、导航误差与Uniform偏置模型；不是实际装备性能。未用5%/10%作为项目验收线，也没有赋予R4进度信用。

{chr(10).join(benefits)}

## 四参数结果

下表均为**12个场景中最坏的无条件经验P95**，失败按INF；不是把全部样本池化后取最有利结果。完整中位数、P90、P95、P99、成功条件统计和两侧分层见CASE_STATISTICS/GLOBAL_STATISTICS。每配置6000次，12个已见开发状态各500次，两侧各250次。四几何两质量共享6000组标准化基础创新，48000不是48000组独立观测。

{chr(10).join(tablelines)}

## 信息来源与偏置条件

{chr(10).join(infotable)}

局部信息共192个名义几何、384行含已知偏置信息对照。全部解析—有限差分控制通过；完整矩阵、参数尺度、弱方向和剖面速度信息均保存。这里的1.96σ不是经验P95，也不是全局保证。主MC并未获得已知偏置：它联合估计两个干扰参数，并使用Uniform支持的硬边界。没有把Uniform方差倒数当作Gaussian校准信息。

测量质量从Q1降至Q2同时改变方位、导航和固定偏置范围，是完整合同对照，不把三者变化解释成单因素因果。估计器包含一阶导航等效Gaussian似然及其log方差项；这是明确的局部误差近似，不是对真实导航分布、相关误差或接收SNR的认证。

## 原始观测、初始化与失败

所有96个公开输入块在拟合前封存，12个私有创新块隔离保存。Q/B内每个记录的时间、两节点报告轨迹和方位均可直接重用；同质量、同场景下四组主节点输入逐字节相同。估计器只接收这些公开数组，使用射线交会与5轮IRLS初值，偏置初始0，失败时按冻结物理中心回退。没有读取真实偏置、真实节点位置或目标状态作初值。

48,000个有限局部终点不是连续全域候选证书。全部失败、优化状态、边界触碰、偏置触界及弱秩诊断保留；没有通过真值重新选择输出。各方案配对变差也保留：

{chr(10).join(pairlines)}

每组最大场景偏置触界率和物理边界触碰见FOUR_PARAMETER_DEMO_TABLE.csv。触界干扰估计不是已完成校准；偏置和运动状态仍可能耦合。

## 资源、覆盖与研究新意

{chr(10).join(resource_table)}

四组辅助节点均与主平台共速平移，维持固定全局偏移；没有新增独立机动。B3横向15 km并前置3 km，节点间距15.297 km，约为B0的3.06倍。所有方案观测时长、节点数和历元数相同，但部署、通信、相对导航与关联成本不同。

距离落在名义研究范围只表示几何适用性；不能推出会聚区探测覆盖或全程能稳定取得登记方位质量。实际硬件、测向、导航能力仍UNKNOWN，同步与目标关联被假定而未验证。没有新增深度量，深度NOT_EVALUATED，不能填写零误差或宣称五维联合定位。

增大基线的几何信息收益本身不是新算法。值得保留的申请研究问题是有限远距离观测与偏置、资源、覆盖条件共同作用时的弱方向与可辨识性。本轮整体定位为**条件水平Demo及资源—信息设计参照**；缺少接收链证据时不升级工程能力。建议停止扩大同类布放扫描，下一研究需负责人单独选题授权。

## 独立复核、预算与停止

冷审计未导入主生成器、估计器、信息或统计代码；从已存创新和轨迹独立重建H01全部4000个配置记录，复核全部48000条误差和似然、全部分位数及配对计数，并独立重建192个Fisher几何。没有冷审计优化重跑。详细范围与检查数见INDEPENDENT_COLD_AUDIT.json。

主分析一次、冷审计一次、科研重试0次。科学链约{budget['wall_s']:.3f}秒，峰值进程树约{budget['peak_tree_RSS_bytes']/1024**2:.2f} MiB。新增原始数据在D盘；源文件及公开/私有输入均以SHA256索引，最终体积见OUTPUT_MANIFEST。前置阶段修正了历史E1目录别名、Git行尾检出及历史函数导入路径；暂停记录保留，这些均发生于Commit A和正式科研前，未更改历史内容或看到登记结果后调参。

新增KRAKEN/KRAKENC/FIELD/FIELDS/SCOOTER、T2补算、径向/频率特征、波形与深度试验均为0。T2/FUNC2 STOP及全部历史Gate保持，R4=0%。提交推送并核对后STOP，不自动开新一轮。
''';(OUT/'GPT_SYNC_T3_W1.md').write_text(sync,encoding='utf-8')
 write(OUT/'DATA_AND_STATISTICAL_SCOPE.json',{'data':'syntheticdirectedbearingfeaturesandreportedpositions,noreceivedwaveforms','state':'4horizontalparametersplus2fixedbiasnuisances;depthnotobserved','baseinnovations':6000,'configurations':48000,'cases':12,'geometries':'seenpanel,newseed;noexternalgeometricconfirmation','pairednotindependent':'B/Qsameinnovations;allmaininputsbyteequalacrossBwithinQ','P95':'nearest-rankceil(.95*n),500case/250side;unconditionalINFfailures','likelihood':'firstorderGaussianmarginalnavigationapproximation;state-dependentvarianceandlognormalizationincluded','prior':'UniformsupportofcommonandHALFdifferential;noGaussianpenalty','Fisher':'localunbiasedGaussianreference;Uniformboundarynotadded;knownbiasinformationcontrolnotestimator','physicaldomain':'sameboundedr,theta,v,psi forallconfigurations;boundaryratesreported;noglobalbranchguarantee','novelty':'baselinealoneisstandardgeometry,notnewinvention','R4_percent':0})
 logs=OUT/'runtime_logs';logs.mkdir(exist_ok=True)
 for p in (LOCAL/'runtime_logs').glob('*'):shutil.copyfile(p,logs/p.name)
 row=f"R4-T3-W1,../R4_T3_W1_WIDE_BASELINE_DEMO/FINAL_ADJUDICATION.json,CONDITIONAL_WIDE_BASELINE_RESOURCE_INFORMATION_DEMO,PAIRED48000_CONFIGS_6000_BASE_DRAWS;NO_RECEIVED_CHAIN;DEPTH_NOT_EVALUATED;R4_delta=0,COLD{audit['checks']}_CHECKS_{audit['FAIL']}FAIL;HISTORICAL_UNCHANGED;RESEARCH_LEAD_PENDING;STOP,0\r\n"
 ledger.write_bytes(original+row.encode('ascii'));write(OUT/'LEDGER_APPEND_PROVENANCE.json',{'old_sha256':hashlib.sha256(original).hexdigest(),'new_sha256':sha(ledger),'old_bytes':len(original),'preserved_exact_prefix':ledger.read_bytes().startswith(original),'one_row_appended':True,'R4_delta':0})
 write(OUT/'DELIVERY_VALIDATION.json',{'main_seal_preserved':True,'frozen_code_unchanged':True,'historical_source_unchanged_except_zero_credit_ledger_append':True,'main_calls':1,'cold_calls':1,'scientific_retries':0,'cold_PASS':cold_ok,'depth':'NOT_EVALUATED','R4_percent':0,'STOP':True})
 write(OUT/'OUTPUT_MANIFEST.json',{'files':[{'path':p.relative_to(OUT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'],'D_archive_files':len(rows(OUT/'PUBLIC_BLOCK_IDENTITIES.csv'))+len(rows(OUT/'PRIVATE_INNOVATION_IDENTITIES.csv')),'D_raw_data_indices':['PUBLIC_BLOCK_IDENTITIES.csv','PRIVATE_INNOVATION_IDENTITIES.csv'],'scope':'allGitdeliveryfilesexceptmanifest;largedataremainDindexed','R4_percent':0,'STOP':True})
 print('DELIVERY_COMPLETE',json.dumps(final,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
