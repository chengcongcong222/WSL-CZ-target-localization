"""Saved-result reporting only; no science implementation import, no new ray or quadrature execution."""
from pathlib import Path
import json,csv,hashlib,time,threading
import numpy as np
import psutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZROB02');start=time.time();end=threading.Event();peak=[0]
def monitor():
 while not end.wait(.1):
  p=psutil.Process();peak[0]=max(peak[0],sum(q.memory_info().rss for q in [p]+p.children(recursive=True) if q.is_running()))
threading.Thread(target=monitor,daemon=True).start()
def h(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p,a):p.write_text(json.dumps(a,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def text(p,s):p.write_text(s,encoding='utf-8',newline='\n')
def read(p):return list(csv.DictReader(p.open(encoding='utf-8-sig')))
def tab(p,a):
 with p.open('w',encoding='utf-8',newline='') as f:
  d=csv.DictWriter(f,list(dict.fromkeys(k for v in a for k in v)));d.writeheader();d.writerows(a)
f=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'));p=json.loads((O/'CZROB02_PRIMARY_RESULTS.json').read_text(encoding='utf-8'));cold=json.loads((O/'CZROB02_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8'));rows=read(O/'CZROB02_RETURN_CURVES.csv');top=read(O/'CZROB02_SSP_TOPOLOGY.csv');paths=read(O/'CZROB02_PATH_ATTRIBUTION.csv');qc=read(O/'CZROB02_COLD_QUADRATURE.csv');failed=[a for a in cold['checks_detail'] if not a['PASS']]
# Formal outputs are immutable; reporting adds diagnostics without revising any frozen scientific verdict.
assert len(p['cases'])==9 and cold['checks']==39601 and cold['FAIL_count']==6 and not cold['PASS'];assert all(not a['numerical_PASS'] and not a['geometry_stable'] for a in p['cases'])
summary=[];worst=[];coldfail=[]
for c in p['cases']:
 name=c['case_id'];rec={'case_id':name,'old_primary_node_distance_m':c['old_primary_node_distance_m'],'old_refined_node_distance_m':c['old_refined_node_distance_m'],'first_candidate_shift_m':c['range_shift_km']*1000,'turn_depth_shift_m':c['turn_depth_shift_m'],'formal_geometry_admitted':False,'formal_classification':c['classification']}
 for m in ['C_LINEAR','PCHIP']:
  mm=c['methods'][m];ff=mm['first_fold'];lookup=name+('_C_LINEAR_angle' if m=='C_LINEAR' else '_PCHIP_adaptive_angle')+str(ff['angle_index']);test=next(a for a in cold['checks_detail'] if a['check']==lookup)
  rec[m+'_candidate_range_km']=ff['deep_return_km'];rec[m+'_candidate_angle_deg']=ff['angle_deg'];rec[m+'_fold_count']=mm['folds_count'];rec[m+'_candidate_cold_delta_m']=test['detail']*1000;rec[m+'_candidate_cold_PASS']=test['PASS'];rec[m+'_full_angle_Q32_Q64_max_delta_m']=mm['quadrature_max_delta_km']*1000;rec[m+'_full_angle_quadrature_PASS']=mm['quadrature_PASS'];rec[m+'_ODE_runs']=mm['ODE_cases']
  rr=[a for a in rows if a['case_id']==name and a['SSP_method']==m and a['status']=='FINITE'];w=max(rr,key=lambda a:float(a['max_order_difference_km']));worst.append(dict(w,diagnostic='MAIN_Q32_VS_Q64_FULL_ANGLE_MAXIMUM'))
 summary.append(rec)
for a in failed:
 name=next(c['case_id'] for c in p['cases'] if a['check'].startswith(c['case_id']));m='C_LINEAR' if '_C_LINEAR_angle' in a['check'] else 'PCHIP';idx=int(a['check'].rsplit('angle',1)[1]);v=next(v for v in rows if v['case_id']==name and v['SSP_method']==m and int(v['angle_index'])==idx);coldfail.append({'check':a['check'],'case_id':name,'SSP_method':m,'angle_index':idx,'angle_deg':v['angle_deg'],'independent_reference_delta_m':a['detail']*1000,'deep_turn_z_m':v['deep_turn_z_m'],'nearest_SSP_node_distance_m':v['nearest_node_distance_m'],'main_deep_range_km':v['deep_return_km'],'main_upper_range_km':v['upper_return_km'],'old_result_or_gate_changed':False})
tab(O/'CZROB02_SAVED_CANDIDATE_DIAGNOSTICS.csv',summary);tab(O/'CZROB02_WORST_ORDER_DIFFERENCES.csv',worst);tab(O/'CZROB02_NUMERICAL_FAILURE_ATTRIBUTION.csv',coldfail)
# Figures show both first-candidate observations and the failed all-angle criterion.
fig,ax=plt.subplots(3,3,figsize=(13,10),layout='constrained')
for axis,c in zip(ax.flat,p['cases']):
 name=c['case_id']
 for m,col in [('C_LINEAR','black'),('PCHIP','tab:blue')]:
  a=[a for a in rows if a['case_id']==name and a['SSP_method']==m and a['status']=='FINITE' and a['deep_bottom_bounce']=='False'];axis.plot([float(v['angle_deg']) for v in a],[float(v['deep_return_km']) for v in a],c=col,lw=.8,label=m);first=c['methods'][m]['first_fold'];axis.plot(first['angle_deg'],first['deep_return_km'],'o',c=col,ms=3)
 axis.set_title(name.replace('SCS17_',''));axis.set_xlim(.1,20);axis.set_ylim(10,120);axis.grid(alpha=.2);axis.set_xlabel('Downward launch angle (deg)');axis.set_ylabel('Deep return range (km)')
ax.flat[0].legend(fontsize=8);fig.suptitle('Registered range window; saved ray-integral candidates — global numerical gate NOT PASSED')
for suffix in ['png','pdf']:fig.savefig(O/f'figures/CZROB02_return_curves.{suffix}',dpi=160)
plt.close(fig)
fig,ax=plt.subplots(3,3,figsize=(13,10),layout='constrained')
for axis,c in zip(ax.flat,p['cases']):
 name=c['case_id']
 for m,col in [('C_LINEAR','black'),('PCHIP','tab:blue')]:
  a=[a for a in rows if a['case_id']==name and a['SSP_method']==m and a['status']=='FINITE'];axis.plot([float(v['angle_deg']) for v in a],[max(1e-9,float(v['max_order_difference_km'])*1000) for v in a],c=col,lw=.7,label=m)
 axis.axhline(5,color='tab:red',ls='--',lw=.9,label='Frozen 5 m limit');axis.set_yscale('log');axis.set_ylim(1e-9,100);axis.set_xlim(.1,20);axis.set_title(name.replace('SCS17_',''));axis.set_xlabel('Launch angle (deg)');axis.set_ylabel('Q32/Q64 range difference (m)');axis.grid(alpha=.2)
ax.flat[0].legend(fontsize=7);fig.suptitle('All registered angles retained — display floor 1e-9 m; no exclusion of difficult rays')
for suffix in ['png','pdf']:fig.savefig(O/f'figures/CZROB02_numerical_precision.{suffix}',dpi=160)
plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained');labels=[a['case_id'].replace('SCS17_','').replace('_z','/') for a in summary];xx=np.arange(9)
axes[0].bar(xx,[a['first_candidate_shift_m'] for a in summary]);axes[0].set_ylabel('First candidate range shift (m)');axes[1].bar(xx,[a['turn_depth_shift_m'] for a in summary],color='tab:orange');axes[1].set_ylabel('Turning-depth shift (m)')
for a in axes:a.set_xticks(xx,labels,rotation=45,ha='right');a.grid(axis='y',alpha=.2)
fig.suptitle('C-linear vs C1 PCHIP: candidate diagnostics only; none is globally certified')
for suffix in ['png','pdf']:fig.savefig(O/f'figures/CZROB02_candidate_shifts.{suffix}',dpi=160)
plt.close(fig)
tab(O/'FIGURE_DATA_PROVENANCE.csv',[{'figures':'figures/CZROB02_return_curves.png/.pdf; figures/CZROB02_numerical_precision.png/.pdf','input':'CZROB02_RETURN_CURVES.csv','sha256':h(O/'CZROB02_RETURN_CURVES.csv'),'scope':'All 9 cases; registered launch/range domains; no performance or energy inference'},{'figures':'figures/CZROB02_candidate_shifts.png/.pdf','input':'CZROB02_SAVED_CANDIDATE_DIAGNOSTICS.csv','sha256':h(O/'CZROB02_SAVED_CANDIDATE_DIAGNOSTICS.csv'),'scope':'Descriptive saved candidates; does not revise failed formal Gates'}])
usage=[json.loads((O/(x+'_PROCESS_USAGE.json')).read_text(encoding='utf-8')) for x in ['ANALYSIS','COLD']];wall=sum(a['elapsed_s'] for a in usage);rss=max(a['tree_peak_RSS_bytes'] for a in usage);shift=[a['first_candidate_shift_m'] for a in summary];turnshift=[a['turn_depth_shift_m'] for a in summary];clcount=[a['C_LINEAR_fold_count'] for a in summary];pcount=[a['PCHIP_fold_count'] for a in summary];maxorder=max(float(v['max_order_difference_km'])*1000 for v in worst);minorder=min(float(v['max_order_difference_km'])*1000 for v in worst);cd=[a['independent_reference_delta_m'] for a in coldfail]
text(O/'CZROB02_RESEARCH_DECISION.md',f'''# CZROB-02：科研去留

结论：**CZROB02_NUMERICAL_ADMISSION_NOT_CLOSED_AND_ENERGY_ATTRIBUTION_UNRESOLVED**。执行流程完成，物理稳健性没有获得本轮正式准入；不建议现在开启高精度相干提供器开发。本轮STOP，历史REC01接受结论与R4=0%保持。

## 1. 获得的物理诊断

同一WOA派生节点下，原分段线性与C1 PCHIP的9个单元均仍出现首深回返局部极小值。保存的首候选位置相差{min(shift):.3f}–{max(shift):.3f} m，转向深度相差{min(turnshift):.3f}–{max(turnshift):.3f} m；18个首候选点的独立参考比较均通过登记的5 m比较要求。它们支持一个有限观察：**本次C1对照没有显示首候选整体消失或公里级迁移。** 不能因转向点靠近原节点就断言CZ是假的。

与此同时，登记曲线中局部极小值计数从C-linear的{min(clcount)}–{max(clcount)}个变为PCHIP的{min(pcount)}–{max(pcount)}个，细节拓扑依赖表示。该计数是离散角网格诊断，不能当作真实焦散总数，且仍受下述数值未闭合限制。标准PCHIP无观测后调参，未增加或删除SSP节点，也未改变底深。

**以上均为保留的局部诊断，不构成稳健第一CZ正式通过。** 全域Gate未过，不将局部好结果升级为全域认证，也不把新PCHIP几何配上旧C-linear总能量后宣称新能量增强。

## 2. 本轮失败属于哪个层次

全部9单元、两种表示共18个全角域32/64阶积分比较都超过预冻结5 m要求，各组最大差{minorder:.3f}–{maxorder:.3f} m。故9个单元均为GEOMETRY_NUMERIC_UNCERTIFIED；代表射线ODE按条件门取消，实际执行0，不能称为独立射线验证通过。

冷审计完成{cold['checks']:,}项，{cold['checks']-cold['FAIL_count']:,}项通过、{cold['FAIL_count']}项失败。失败集中在JW_z200的0.40845°、JW_z250的0.72685°、ANNUAL_z250的2.33875°，C-linear和PCHIP各3项；独立参考差{min(cd):.3f}–{max(cd):.3f} m。没有重试、加阶、缩角域或放宽门槛。

失败角的转向点落在节点下方约3.8、8.3、15.9 mm，两种表示都出现。冻结主实现只对**包含转向点的一段**进行奇性变换，转向刚越过节点时，前一段端点仍有尖锐近奇性；该段固定Gauss求积可能不足。C-linear完整曲线的独立解析原函数以及PCHIP独立适应积分在这三处提供了直接差异证据。此为**冻结积分实现的精度限制**，不是介质在C1下失去CZ的物理证据，也不是对Snell积分或射线方法族的否定。

冷审计未通过，因此不能将本轮所有诊断统一称为“已认证”。失败原数据、主输出封存清单与冷审计原件全部交付。没有为给出正结论改写分类。

## 3. 边界路径归因获得了什么

原C-linear R0/R1的全部源深穿越已按**穿越前**的边界接触重新分类，独立冷解析的18组统计全部一致。原焦区中确有零边界回返，也有此前已接触边界的回返；每单元每级射线有77–87个边界类焦区穿越。数量不等于独立声源数、能量比例或污染强度；同一条射线可有多个回返。

REFRACTED_PATH_IDENTITY_SUPPORTED仅限旧C-linear存档零边界路径身份。边界几何重叠得到支持，但若BOUNDARY_PATH_CONTAMINATION_PRESENT意指非零或主要能量贡献，本轮没有足够幅度证据，仍为NOT_ESTABLISHED_AS_ENERGY。全部9单元的EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED保持。

旧RAY只含向下发射扇，不能穷举初始向上及所有高阶路径；I类SHD没有各支能量标签，不能反演出可验证分解。新PCHIP的总TL、反射/折射能量比均NOT_EVALUATED。本轮没有发明替代能量分解、没有混用非相干能量与相干场。

## 4. 对项目申请及后续投入的意义

目前可保留的是：同一个WOA23_17N118E资料族、气候均值、平底刚性、C-linear名义条件下已有第一深回返及总TL增强证据；新C1对照出现相近的首回返候选。它们是**适用条件与机制线索**，不是9个海域的环境验证，也不是Yang/MMAR可提取模态或测深性能。

关口A仍未完成全域数值认证，主要路径的能量贡献也未分开，因此暂缓关口B的高精度相干提供器投入。不能把此次失败写成环境不成立、SSP插值制造了全部CZ，或方法已不可行。

解除条件分别是：给出可独立验证的跨节点近转向积分/射线证据；以及带路径标签、具有幅度与非相干加和合同的能量证据，或其它能区分折射和边界贡献的合格观测。两者都是下一次独立授权才可执行的研究要求，不是本轮自动修补清单。旧T2/FUNC2继续关闭，本轮不提出或启动连续调参链。

## 5. 执行与资源

主分析一次、冷审计一次；累计科学壁钟{wall:.2f} s，采样进程树RSS峰值{rss/1024**2:.1f} MiB。外部声学求解0、MC0、新SHD0、Yang/MMAR性能0。全部新增数据和缓存位于D盘。完成Commit B推送后STOP，不自动开展后续研究。
''')
final={'stage':'ENV_CAP_CZROB02','parent_SHA':f['parent_SHA'],'design_SHA':'3d2a56e1821c3c94949346bb5d9e80170dfc221b','branch':f['branch'],'status':'STOP','execution':'COMPLETED_NEGATIVE_BRANCH_NO_RETRY','decision':'CZROB02_NUMERICAL_ADMISSION_NOT_CLOSED_AND_ENERGY_ATTRIBUTION_UNRESOLVED','science_gate_PASS':False,'cold_audit':{'PASS':False,'checks':cold['checks'],'PASS_count':cold['checks']-cold['FAIL_count'],'FAIL_count':cold['FAIL_count']},'A':{'cases_retained':9,'source_families':1,'full_angle_return_rows':len(rows),'numeric_admitted_cases':0,'formal_robust_CZ_cases':0,'all18_full_angle_integral_gates_failed':True,'SSP_INTERPOLATION_SENSITIVE':'LOCAL_MINIMUM_COUNTS_CHANGE_OBSERVED_ONLY; GLOBAL_ROBUSTNESS_UNCERTIFIED','CZ_GEOMETRY_NOT_ROBUST':'NOT_ESTABLISHED_PHYSICALLY_NUMERICAL_ADMISSION_FAILED','diagnostic_candidate_shift_m':[min(shift),max(shift)],'diagnostic_turn_shift_m':[min(turnshift),max(turnshift)],'independent_ODE_runs':0,'ODE_status':'NOT_EVALUATED_CANCELLED_BY_QUADRATURE_GATE'},'B':{'old_C_linear_refracted_identity':'SUPPORTED_FROM_SAVED_RAYS','boundary_geometry_overlap':'PRESENT','BOUNDARY_PATH_CONTAMINATION_PRESENT':'NOT_ESTABLISHED_AS_ENERGY','exclusive_energy_attribution':'EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED','PCHIP_energy':'NOT_EVALUATED'},'C':{'invest_coherent_provider_now':'DEFER','reason':'GLOBAL_GEOMETRY_NUMERICAL_ADMISSION_AND_PATH_ENERGY_ATTRIBUTION_NOT_CLOSED','Yang_MMAR_depth':'NOT_EVALUATED'},'budget':{'science_wall_s':wall,'sampled_peak_RSS_bytes':rss,'new_external_solver_calls':0,'new_MC':0,'new_acoustic_fields':0},'history':{'ENV_CAP':'UNCHANGED','REC01':'UNCHANGED_ACCEPTED','T2_FUNC2':'STOP_UNCHANGED','R4':'0%','main':'UNCHANGED'},'next_stage_execution_authorized':False}
js(O/'CZROB02_FINAL_DECISION.json',final)
text(O/'GPT_SYNC_CZROB02.md',f'''# ENV-CAP-CZROB-02：最终同步

**结论：CZROB02_NUMERICAL_ADMISSION_NOT_CLOSED_AND_ENERGY_ATTRIBUTION_UNRESOLVED。** 本轮负结果分支交付完成并STOP，不建议现在启动相干提供器开发。

- 父SHA：`{f['parent_SHA']}`。
- 设计SHA：`{final['design_SHA']}`。
- 最终SHA：本报告所在Commit B；完整身份见Git科研分支HEAD。
- 分支：`{f['branch']}`；不合并main。

## 核心科研结果

同一WOA派生节点上的C-linear和C1 PCHIP，在9个单元中都出现首深回返候选；保存的首候选移动{min(shift):.2f}–{max(shift):.2f} m，转向深度移动{min(turnshift):.2f}–{max(turnshift):.2f} m。18个首候选点的独立参考比较通过。局部极小值计数由C-linear的{min(clcount)}–{max(clcount)}个变为PCHIP的{min(pcount)}–{max(pcount)}个，显示细节表示依赖。这些是**局部诊断，非稳健CZ通过**。

18个全角域积分阶次Gate全部失败，最大差{minorder:.2f}–{maxorder:.2f} m；因此0/9取得本轮数值准入，代表射线ODE按条件门取消。**0/9不是CZ消失或物理不可行。** 冷审计{cold['checks']:,}项中{cold['FAIL_count']}项失败，集中在三个近节点转向角，独立差{min(cd):.2f}–{max(cd):.2f} m。冻结实现对跨节点近转向的积分精度不足，没有修补、重试或放宽标准。

18组存档射线前缀分类统计经独立解析核对一致：原焦区同时有零边界与边界接触回返。**RAY计数不能换算能量占比。** 全部9单元能量归因仍EXCLUSIVE_ENERGY_ATTRIBUTION_UNRESOLVED；PCHIP总TL与各支能量均NOT_EVALUATED。

## 研究去留

现有局部诊断没有显示C1对照令首候选整体消失，但全域精度和主要路径能量归因均未闭合。暂缓高精度相干提供器；不启动Yang/MMAR性能，不重开T2/FUNC2。原REC01已接受的条件证据保持，三套气候组合仍只算一个来源族。新研究没有证明真实海域CZ、可检测性或测深精度。

## 执行核对

922项历史输入身份检查通过，冻结后931项执行身份检查；25项合成控制通过。主分析一次，冷审计一次，科学壁钟{wall:.2f} s，采样RSS峰值{rss/1024**2:.1f} MiB。新增外部声学求解0、MC0、新接收/SHD0，所有新增文件在D盘。main、历史STOP及R4=0%不变。

详细结论见[CZROB02_RESEARCH_DECISION.md](CZROB02_RESEARCH_DECISION.md)，所有失败见[CZROB02_INDEPENDENT_COLD_AUDIT.json](CZROB02_INDEPENDENT_COLD_AUDIT.json)，两种表示全曲线、原节点/系数、边界路径穿越和原始哈希索引均一并交付。交付完整性检查与科学Gate必须分别阅读。
''')
js(O/'STOP.json',{'stage':'ENV_CAP_CZROB02','status':'STOP','reason':final['decision'],'cold_audit_FAIL':6,'cancelled':['representative ODE under failed full-angle quadrature Gate','coherent-provider development','Yang/MMAR ideal/nonoracle depth performance'],'no_science_retry':True,'history_unchanged':True,'R4':'0%','next_stage_requires_new_authorization':True,'time_unix':time.time()})
text(O/'REVIEW_ENTRYPOINT.md','''# CZROB02审核入口

先读GPT_SYNC_CZROB02.md及CZROB02_RESEARCH_DECISION.md，再读最终JSON与冷审计。注意：交付完整不等于科学Gate通过；本轮冷审计6项FAIL全部保留。

设计代码与阈值在Commit A固定，主分析与冷审计各一次。PRIMARY_OUTPUT_FREEZE.json封存冷审计前的输入。CZROB02_SSP_TOPOLOGY.csv是正式分类；CZROB02_SAVED_CANDIDATE_DIAGNOSTICS.csv仅为保存结果的诊断汇总，不替代失败Gate。全部72,018个角度—表示—单元记录在CZROB02_RETURN_CURVES.csv。原节点/系数、节点邻域探针、全体局部最小值均保留。

边界路径的几何存在与能量贡献必须分开。原C-linear RAY有前者，存档SHD没有可用的分支能量标签；不能把统计回返次数作为能量权重，也不能把旧TL移植到PCHIP。STOP后不自动重跑、加阶或开启相干提供器。
''')
# Index old inputs and new execution evidence without copying sealed raw data.
idx=[{'scope':a['source'],'path':a['path'],'bytes':Path(a['path']).stat().st_size,'sha256':a['sha256']} for a in f['execution_inputs']]
for pth in sorted(W.iterdir()):
 if pth.is_file() and pth.suffix in ('.json','.log','.npz'):
  idx.append({'scope':'CZROB02_new_execution','path':str(pth),'bytes':pth.stat().st_size,'sha256':h(pth)})
  if pth.suffix in ('.json','.log'):(O/'inputs'/pth.name).write_bytes(pth.read_bytes())
tab(O/'RAW_DATA_INDEX.csv',idx)
js(O/'RUN_STATE.json',{'phase':'STOP','formal_runs':1,'cold_runs':1,'main_execution_exit':0,'cold_execution_exit':2,'cold_exit2_is_reported_audit_failure_not_crash':True,'scientific_code_modified_after_A':False,'solver_calls':0,'new_MC':0,'ODE_runs':0,'new_PCHIP_energy':'NOT_EVALUATED'})
end.set();js(O/'DOCUMENTATION_PROCESS_USAGE.json',{'purpose':'saved-result rendering and summaries only','wall_s':time.time()-start,'sampled_peak_RSS_bytes':peak[0]})
print(json.dumps({'scientific_GATE_PASS':False,'cold_checks':cold['checks'],'cold_FAIL':cold['FAIL_count'],'candidate_shift_m':[min(shift),max(shift)],'max_order_difference_m':maxorder,'science_wall_s':wall,'science_peak_MiB':rss/1024**2,'documentation_peak_MiB':peak[0]/1024**2,'STOP':True}))
