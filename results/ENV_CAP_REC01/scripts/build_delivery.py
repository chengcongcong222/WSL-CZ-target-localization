"""Delivery-only summaries and figure from saved REC01 outputs. No scientific rerun."""
from pathlib import Path
import csv,json,hashlib,time,threading,os
import numpy as np
import psutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parents[1];OLD=O.parent/'ENV_CAP_FULL_CAMPAIGN';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_REC01');begin=time.time();peaks=[0];done=threading.Event()
def monitor():
 while not done.wait(.1):
  p=psutil.Process();peaks[0]=max(peaks[0],sum(q.memory_info().rss for q in [p]+p.children(recursive=True) if q.is_running()))
threading.Thread(target=monitor,daemon=True).start()
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def text(p,s):p.write_text(s,encoding='utf-8',newline='\n')
def tab(p,a):
 with p.open('w',encoding='utf-8',newline='') as f:
  d=csv.DictWriter(f,list(dict.fromkeys(k for v in a for k in v)));d.writeheader();d.writerows(a)
cold=json.loads((O/'REC01_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'));assert cold['PASS']
case=json.loads((O/'REC01_CZ_PHYSICS_ADJUDICATION.json').read_text(encoding='utf-8'));(O/'REC01_CZ_PHYSICS_ADJUDICATION_PRE_COLD.json').write_bytes((O/'REC01_CZ_PHYSICS_ADJUDICATION.json').read_bytes());oldcases=json.dumps(case['cases'],sort_keys=True)
case['stage']='REC01_FINAL_AFTER_INDEPENDENT_COLD';case['global_energy_reconstruction_pending_cold']=False;case['independent_cold_PASS']=True;case['cold_audit_sha256']=h(O/'REC01_INDEPENDENT_COLD_AUDIT.json');case['exclusive_refractive_energy_identity']='UNRESOLVED_BOUNDARY_CONTRIBUTIONS_NOT_SEPARATED';assert json.dumps(case['cases'],sort_keys=True)==oldcases;js(O/'REC01_CZ_PHYSICS_ADJUDICATION.json',case)
full=list(csv.DictReader((O/'REC01_FULL_RANGE_ENERGY.csv').open(encoding='utf-8-sig')));curves=list(csv.DictReader((O/'REC01_E1_NUMERICAL_RECONSTRUCTION.csv').open(encoding='utf-8-sig')))
legacy=json.loads((OLD/'E1_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'));mapping=[]
for a in legacy['rows']:
 if a['PASS']:continue
 n=a['check'];cid=next(c['case_id'] for c in case['cases'] if n.startswith(c['case_id']))
 matches=[b for b in curves if b['case_id']==cid and b.get('max_legacy_TL_delta_dB') and float(b['max_legacy_TL_delta_dB'])>1e-9]
 mapping.append({'old_check':n,'old_detail':a['detail'],'case_id':cid,'new_case_class':next(c['classification'] for c in case['cases'] if c['case_id']==cid),'attribution':'legacy rolling signed energy / cancellation, amplified by dB floor; original check remains FAIL','affected_energy_curves':';'.join(b['level'] for b in matches),'largest_reproduced_old_TL_delta_dB':max(float(b['max_legacy_TL_delta_dB']) for b in matches),'old_result_changed':False})
assert len(mapping)==12;tab(O/'REC01_LEGACY_FAILURE_ATTRIBUTION.csv',mapping)
# Full-range refinement differences are descriptive only; no post-hoc new Gate.
refinement=[]
for c in case['cases']:
 cid=c['case_id'];ts={lev:np.array([float(a['positive_TL_floor1e30_dB']) for a in full if a['case_id']==cid and a['level']==lev]) for lev in ['I0','I1','I2']}
 if not len(ts['I0']):continue
 for la,lb in [('I0','I1'),('I1','I2')]:
  delta=abs(ts[la]-ts[lb]);refinement.append({'case_id':cid,'comparison':la+'_'+lb,'range_domain_km':'10..120','all_bins':881,'median_TL_delta_dB':float(np.median(delta)),'P95_TL_delta_dB':float(np.quantile(delta,.95)),'max_TL_delta_dB':float(delta.max()),'interpretation':'FULL_DOMAIN_DIAGNOSTIC_ONLY_NOT_GLOBAL_PROPAGATION_ADMISSION'})
tab(O/'REC01_FULL_DOMAIN_REFINEMENT_DIAGNOSTIC.csv',refinement)
# Publication-style diagnostic, saved data only.
fig,ax=plt.subplots(2,2,figsize=(11,6.5),layout='constrained');provenance=[]
for row,(cid,lev) in enumerate([('SCS17_ANNUAL_z200','I0'),('ESTD_IDEAL_z250','I0')]):
 a=[a for a in full if a['case_id']==cid and a['level']==lev];x=np.array([float(v['range_km']) for v in a]);en=np.array([float(v['positive_window_energy']) for v in a]);rolling=np.array([float(v['legacy_rolling_energy_signed']) for v in a]);old=np.array([float(v['old_saved_TL_dB']) for v in a]);new=np.array([float(v['positive_TL_floor1e30_dB']) for v in a]);ax[row,0].plot(x,en,label='Nonnegative reconstruction',lw=1);ax[row,0].plot(x,rolling,label='Legacy rolling energy',lw=.7,alpha=.7);ax[row,0].set_yscale('symlog',linthresh=1e-30);ax[row,0].set_ylabel('Stored-field energy (model units)');ax[row,0].set_title(cid+' / '+lev);ax[row,1].plot(x,old-new,lw=.8);ax[row,1].set_ylabel('Legacy TL minus reconstructed TL (dB)');ax[row,1].set_title('All ranges retained; no weak-bin deletion')
 for a1 in ax[row]:a1.set_xlabel('Range (km)');a1.grid(alpha=.2);a1.set_xlim(10,120)
ax[0,0].legend(fontsize=8);fig.savefig(O/'figures/REC01_energy_roundoff.png',dpi=180);fig.savefig(O/'figures/REC01_energy_roundoff.pdf');plt.close(fig)
tab(O/'FIGURE_DATA_PROVENANCE.csv',[{'figure':'figures/REC01_energy_roundoff.png; figures/REC01_energy_roundoff.pdf','input':'REC01_FULL_RANGE_ENERGY.csv','input_sha256':h(O/'REC01_FULL_RANGE_ENERGY.csv'),'selection':'two previously known legacy-failure examples shown descriptively; all cases and all bins remain in full CSV; no selection for scientific Gate'}])
usage=[json.loads((O/(p+'_PROCESS_USAGE.json')).read_text(encoding='utf-8')) for p in ['ANALYSIS','COLD']];wall=sum(a['elapsed_s'] for a in usage);peak=max(a['tree_peak_RSS_bytes'] for a in usage)
counts={k:sum(a['classification']==k for a in case['cases']) for k in sorted(set(a['classification'] for a in case['cases']))};peaks_km=[a['focus_metrics'][-1]['peak_km'] for a in case['cases'] if a['ray_TL_conditional_supported']]
res={'stage':'ENV_CAP_REC01','parent_SHA':'a34a1466d21e74b0ea2cfbad53957e1d6c3f49bc','design_SHA':'ebb1f45fddf468c04bedfd027e8bce54d1cbec65','status':'STOP','decision':'CONDITIONAL_CZ_NUMERICS_RECOVERED_COHERENT_PROVIDER_NOT_READY','independent_cold_PASS':True,'cold_checks':cold['checks'],'cold_FAIL':0,'R1':{'all_cases':15,'classes':counts,'full_energy_bins':len(full),'distinct_source_families':1,'source_family':'WOA23_17N118E','exclusive_refractive_energy':'UNRESOLVED','negative_controls_preserved':3,'E_STD_physical_identity_unresolved':3},'R2':'FLOAT64_MIRROR_CONTROL_SUPPORTED_OLD_25_OF_26_UNCHANGED','R3':{'current_coherent_provider':'NOT_READY','new_solver_calls':0,'Yang_CZ_ideal_depth':'NOT_EVALUATED','MMAR_CZ_nonoracle_extraction':'NOT_EVALUATED','reference_paths':2},'R4_stage_research_decision':{'A_nominal_CZ_environment':'YES_CONDITIONAL_ONE_SOURCE_FAMILY','B_coherent_provider':'NO_CURRENT_CERTIFIED_PROVIDER','C_Yang_MMAR_performance_prerequisites':'NOT_READY','D_proven_frontmost_obstacle':'COHERENT_FORWARD_AND_NONORACLE_OBSERVATION_CONTRACT_MISSING','E_resource_lower_bounds':'NOT_ESTABLISHED_FOR_NEW_CZ'},'budget':{'new_MC':0,'new_solver_calls':0,'science_wall_s':wall,'science_peak_RSS_bytes':peak},'historical_R4':'0%','historical_ENV_CAP':'NO_GO_CURRENT_CZ_METHOD_PERFORMANCE_DEMO_UNCHANGED','T2_FUNC2':'STOP_UNCHANGED','next_execution_authorized':False}
js(O/'REC01_FINAL_DECISION.json',res)
text(O/'GPT_SYNC_ENV_CAP_REC01.md',f'''# ENV-CAP-REC-01：阶段同步

正式结论：**CONDITIONAL_CZ_NUMERICS_RECOVERED_COHERENT_PROVIDER_NOT_READY**。REC01已执行并恢复STOP；原ENV-CAP不追溯改判，R4仍为0%。

- 父提交：`{res['parent_SHA']}`。
- 设计提交：`{res['design_SHA']}`。
- 最终提交：本报告所在的Commit B；完整SHA见远端科研分支HEAD和Git提交记录。main不合并。
- 科研分支：`research/env-cap-rec01`。

## 真正完成了什么

完整核验819份历史输入，加冻结代码/合同共825项执行身份检查；15项合成实现控制通过。一次主分析、一次独立冷审计，**{cold['checks']}项通过、0失败**。

全部15个旧环境单元保留：**9个南海单元条件性通过、3个E-STD物理身份未闭合、3个常声速负例**。南海细级能量峰在{min(peaks_km):.3f}–{max(peaks_km):.3f} km。三套气候组合只算WOA23_17N118E一个资料来源族，不能宣称跨海域验证。

**数值实现恢复：** 36条存档SHD曲线、31,716个全域窗口用非负求和重建，冷程序以独立struct解析、math.fsum和Decimal哨兵验证。565个精确零窗口、67个旧滑动结果负窗全部保留。旧12项冷审计失败的7条能量曲线/5项派生比较可追溯到极弱窗舍入与dB地板；最大重现TL差65.72 dB。该差值不是定位误差，也不是有能量焦区误差。

**条件物理证据：** 9个单元满足登记的零边界首深回返、角域折叠与局部总TL增强/细化规则。所有这些焦区都有可能重叠的边界路径，未分离各路径能量；C-linear节点附近的折叠也保留限制。这支持名义模型的条件机制环境，不能证明排他折射能量或真实海域可检测性。全域能量重建通过不等于全域传播精度认证。

**几何软件恢复：** 水平8元、正/负倾斜8元和HLA加4垂向元共16组镜像源/阵组合通过float64控制及80位Decimal距离核对。共x-z平面阵型仍存在y镜像距离恒等关系；没有解决真实左右舷观测歧义，旧E2的25/26保持。

## 仍未解决什么

E1保存的是射线和非相干能量，SHD的复数存储字段不构成相干复场。新SCS合同为无体吸收、理想刚性底，不能移植旧E-STD的Thorp修正、抽样根或函数认证。

已经逐项审查模态波数、函数、完整性、复声压和边界误差，并提出完整高精度模态与独立谱积分两条可追溯路径。当前两者均未就绪，未执行新求解。**Yang/MMAR第一CZ模式提取、测深精度与资源可行域仍NOT_EVALUATED。** 没有用恢复的镜像控制或射线证据跳过相干提供器准入。

## 资源与去留

科学主分析{usage[0]['elapsed_s']:.2f}秒，冷审计{usage[1]['elapsed_s']:.2f}秒，累计{wall:.2f}秒；监测进程树峰值{peak/1024**2:.1f} MiB。新增声学求解0、新MC0、新接收记录0，所有新增材料在D盘。

本轮获得了有价值的条件环境和可复核的软件恢复结果，尚未得到新的深度或多参数性能Demo。停止相应性能研究，保留解除条件与资源清单。下一次是否投资窄样本相干提供器认证，由研究负责人另行决定；不自动启动E3、不恢复T2/FUNC2。

详细证据见[最终判定](REC01_FINAL_DECISION.json)、[数值重建](REC01_E1_NUMERICAL_RECONSTRUCTION.csv)、[物理裁定](REC01_CZ_PHYSICS_ADJUDICATION.json)、[独立冷审计](REC01_INDEPENDENT_COLD_AUDIT.json)、[提供器审查](REC01_COHERENT_FORWARD_FEASIBILITY.md)及[资源缺口](REC01_RESOURCE_AND_RESEARCH_GAPS.md)。
''')
text(O/'REC01_NUMERICAL_ATTRIBUTION.md','''# 全域能量重建的解释边界

旧12项FAIL包含7条SHD重建差异和5条派生细化统计差异。SCS17_ANNUAL_z200只涉及I0极弱区；E-STD占其余11项。逐项映射见REC01_LEGACY_FAILURE_ATTRIBUTION.csv，旧审核原件保持不变。

旧滑动滤波在强弱能量并存时可能保留消去残差或产生极小负数。新程序不以其结果为准，直接加和非负原始能量。显示TL需要地板，线性能量本身不需要；所以同一极弱窗可以出现几十dB的显示差，却不是同样大的焦区声场差异。

所有31,716个窗口均交付，包括565个精确零窗。精确零仅指已存复分量在该窗口为零，不推断真实海洋能量严格为零。新线性能量主冷一致性只认证算术重建。接收动态范围40/60/80 dB与三种显示地板均在执行前冻结，未用它们删除任何范围或改判Gate。

全域I0/I1/I2差异另见REC01_FULL_DOMAIN_REFINEMENT_DIAGNOSTIC.csv；它们没有被焦区规则隐藏，也不参与任何新增事后Gate。焦区准入沿用既有几何定义，E-STD缺少登记角域内内点折叠，因而仍为物理身份未闭合，不能因其线性能量重建成功便改称首CZ成立。

图figures/REC01_energy_roundoff.png/.pdf只是两个已知失败实例的全距离展示；所有案例仍在完整CSV，图选例不参与判定。科学代码及阈值在Commit A以后未修改；独立冷审计后仅给物理裁定增加已审核状态，原始主分析版本单独保留。
''')
js(O/'STOP.json',{'stage':'ENV_CAP_REC01','status':'STOP','time_unix':time.time(),'reason':res['decision'],'cancelled_dependent_work':['Yang/MMAR CZ ideal/nonoracle depth performance','E3 sensitivity and performance MC'],'new_solver_calls':0,'new_MC':0,'history_unchanged':True,'next_stage_requires_new_authorization':True})
text(O/'REVIEW_ENTRYPOINT.md','''# REC01审阅入口

先读GPT_SYNC_ENV_CAP_REC01.md，再读REC01_FINAL_DECISION.json。代码、规则和输入身份锁定在Commit A，正式主分析与独立冷审计均各一次；执行记录和局部原始记录路径见RUN_STATE.json与RAW_DATA_INDEX.csv。

关键核查：31,716个全域能量窗是否完整；旧12项失败是否全部映射；9个条件单元是否仍只算一个资料族；边界路径未分离的限制是否明确；R2新通过是否被错误升级为旧E2通过；相干提供器未准入时是否停止Yang/MMAR性能。FIGURE_DATA_PROVENANCE.csv给出科研图的数据来源。

main、旧ENV-CAP、T2/FUNC2及R4台账不改。本分支STOP，不自动合并main。
''')
# Do not copy raw science files; index original paths and exact identities instead.
freeze=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'));idx=[{'scope':a['scope'],'path':a['path'],'sha256':a['sha256'],'bytes':Path(a['path']).stat().st_size} for a in freeze['execution_inputs']]
for p in sorted(W.glob('*.json')):
 if p.name=='ARCHIVE_IDENTITY_PRECHECK.json':continue
 idx.append({'scope':'REC01_local_execution_record','path':str(p),'sha256':h(p),'bytes':p.stat().st_size})
for p in sorted(W.glob('*.log')):idx.append({'scope':'REC01_execution_log','path':str(p),'sha256':h(p),'bytes':p.stat().st_size})
tab(O/'RAW_DATA_INDEX.csv',idx)
# Copy small provenance logs into Git; large archive remains D-only.
for p in sorted(W.glob('*.json')):
 if p.name=='ARCHIVE_IDENTITY_PRECHECK.json':continue
 target=O/'inputs'/p.name;target.write_bytes(p.read_bytes())
for p in W.glob('*.log'):(O/'inputs'/p.name).write_bytes(p.read_bytes())
js(O/'RUN_STATE.json',{'stage':'ENV_CAP_REC01','phase':'STOP','formal_runs':1,'cold_runs':1,'analysis_launch':'inputs/ANALYSIS_LAUNCH.json','cold_launch':'inputs/COLD_LAUNCH.json','scientific_threshold_changes_after_A':0,'new_solver_calls':0,'new_MC':0,'coherent_performance':'NOT_EVALUATED','raw_files_copied':False})
done.set();js(O/'DOCUMENTATION_PROCESS_USAGE.json',{'purpose':'saved-result summaries/figure only; no scientific execution','elapsed_s':time.time()-begin,'sampled_peak_RSS_bytes':peaks[0]})
print(json.dumps({'cold_PASS':cold['PASS'],'science_wall_s':wall,'classes':counts,'science_peak_MiB':peak/1024**2,'documentation_peak_MiB':peaks[0]/1024**2,'legacy_failures_mapped':len(mapping)}))
