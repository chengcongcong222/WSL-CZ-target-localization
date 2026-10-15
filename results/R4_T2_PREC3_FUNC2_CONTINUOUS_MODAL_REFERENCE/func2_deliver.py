"""FUNC2 post-execution reporting only. No mode, root, ODE, FIM or signal calculations."""
from pathlib import Path
import json,csv,hashlib,shutil,subprocess,datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(name,data):
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:
  fields=list(dict.fromkeys(k for r in data for k in r));w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(data)
def git(*a):return subprocess.check_output(['git','-c','safe.directory='+ROOT.as_posix(),*a],cwd=ROOT).decode().strip()
def main():
 cfg=read(OUT/'FUNC2_DESIGN_FREEZE.json');decision=read(OUT/'FUNC2_DECISION.json');cold=read(OUT/'FUNC2_INDEPENDENT_COLD_AUDIT.json');budget=read(OUT/'EXECUTION_BUDGET.json');seal=read(OUT/'MAIN_OUTPUT_SEAL.json');checks=[]
 assert decision['A']=='FAIL' and decision['B']==decision['C']==decision['D']=='NOT_EVALUATED'
 assert decision['stop_reason']=='SHOOT_ROOT_NOT_CONVERGED'
 for r in cfg['inputs']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
 for r in cfg['vectors']:assert sha(r['path'])==r['sha256'],r['path']
 for r in cfg['frozen_code']:assert sha(OUT/r['path'])==r['sha256'],r['path']
 for r in seal['files']:assert sha(r['path'])==r['sha256'],r['path']
 assert budget['main_calls']==budget['cold_calls']==1 and budget['retries']==0 and not read(OUT/'FUNC2_RAW_DATA_INDEX.json')['files']
 original=(LOCAL/'LEDGER_BEFORE.bin').read_bytes();ledger=ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv';assert ledger.read_bytes()==original
 output_logs=OUT/'runtime_logs';output_logs.mkdir(exist_ok=True)
 for p in (LOCAL/'runtime_logs').glob('*'):shutil.copyfile(p,output_logs/p.name)
 shutil.copyfile(LOCAL/'RUNTIME_DEPENDENCY_PROVENANCE.json',OUT/'RUNTIME_DEPENDENCY_PROVENANCE.json')
 gate_rows=[{'stage':'A','item':'analytic method-level-mode-loss endpoint controls','planned':24,'attempted':1,'completed_result_rows':0,'failure_count':1,'not_evaluated':23,'status':'FUNC2_ANALYTIC_REFERENCE_GATE_FAIL','k_real_imag_tolerance':cfg['k_budget'],'shape_tolerance':cfg['shape_budget'],'max_observed_k_deviation':None,'max_observed_shape_deviation':None,'reason':'DOP853 level0 mode1 loss0 root exhausted10iterations; no converged root or function saved'},
 {'stage':'B','item':'E-STD complex continuation roots','planned':48,'attempted':0,'completed_result_rows':0,'failure_count':0,'not_evaluated':48,'status':'NOT_EVALUATED','reason':'Gate A failed before any E-STD continuous solve'},
 {'stage':'B','item':'independent continuous mode references','planned':3,'attempted':0,'completed_result_rows':0,'failure_count':0,'not_evaluated':3,'status':'NOT_EVALUATED','reason':'Gate A failed'},
 {'stage':'C','item':'original FE/FD function versus continuous reference','planned':10,'attempted':0,'completed_result_rows':0,'failure_count':0,'not_evaluated':10,'status':'NOT_EVALUATED','reason':'no admitted reference; FD961 not authorized'},
 {'stage':'C','item':'200m error bands','planned':250,'attempted':0,'completed_result_rows':0,'failure_count':0,'not_evaluated':250,'status':'NOT_EVALUATED','reason':'no admitted reference'},
 {'stage':'D','item':'original FE/FD Q versus continuous reference','planned':200,'attempted':0,'completed_result_rows':0,'failure_count':0,'not_evaluated':200,'status':'NOT_EVALUATED','reason':'no admitted reference'},
 {'stage':'D','item':'reference-path and own-precision Q comparisons','planned':180,'attempted':0,'completed_result_rows':0,'failure_count':0,'not_evaluated':180,'status':'NOT_EVALUATED','reason':'no admitted reference'},
 {'stage':'COLD','item':'identity checks actually available','planned':160,'attempted':160,'completed_result_rows':160,'failure_count':cold['FAIL'],'not_evaluated':0,'status':'INTEGRITY_ONLY_PASS','reason':'144 parent inputs +10 NPZ +6 frozen Python files'},
 {'stage':'COLD','item':'fresh highmode RK45 reconstruction','planned':1,'attempted':0,'completed_result_rows':0,'failure_count':0,'not_evaluated':1,'status':'NOT_EVALUATED','reason':'no admitted highmode root; cold scope must not imply numeric certification'},
 {'stage':'E','item':'go/no-go document','planned':1,'attempted':1,'completed_result_rows':1,'failure_count':0,'not_evaluated':0,'status':'NO_GO_CURRENT_IMPLEMENTATION','reason':'close current numerical repair chain; no automatic new study'}]
 table('FUNC2_STAGE_COMPLETION.csv',gate_rows)
 register=[]
 for mode in [1,481,961]:
  for stage in ['E_STD_CONTINUOUS_REFERENCE','FE_VS_REFERENCE','FD_VS_REFERENCE','LOCAL_COUPLING','REFERENCE_K_FUNCTION_BOUNDARY_COLD_RECONSTRUCTION']:
   register.append({'mode':mode,'stage':stage,'status':'NOT_EVALUATED_FD_VECTOR_NOT_ARCHIVED' if mode==961 and stage=='FD_VS_REFERENCE' else 'NOT_EVALUATED','reason':'FD961 not archived and recomputation prohibited' if mode==961 and stage=='FD_VS_REFERENCE' else 'analytic reference gate failed before admission'})
 table('FUNC2_NOT_EVALUATED_REGISTER.csv',register)
 final={'stage':'T2-PREC3-FUNC2','parent_SHA':cfg['parent_SHA'],'design_SHA':read(OUT/'EXECUTION_STARTED.json')['design_SHA'],'classifications':decision['classifications'],'scientific_result':'NEW_REFERENCE_NOT_ADMITTED_AT_ANALYTIC_GATE','first_failed_case':{'method':'DOP853','level':0,'mode':1,'environment':'constant1500m/s,5000m,pressure_release_surface,rigid_bottom','loss':0,'reason':'SHOOT_ROOT_NOT_CONVERGED','max_registered_iterations':10},'A':{'planned_endpoint_cases':24,'attempted':1,'completed':0,'failed':1,'deviation_metrics':'UNAVAILABLE_NO_CONVERGED_REFERENCE_SAVED'},'B':'NOT_EVALUATED','C':'NOT_EVALUATED','D':'NOT_EVALUATED','cold_integrity_checks':160,'cold_integrity_FAIL':0,'cold_numeric_reconstructions':0,'cold_scope':'IDENTITY_ONLY_NOT_K_FUNCTION_OR_COUPLING_CERTIFICATION','continuous_reference_supported_modes':[],'FE_vs_FD_accuracy_ranking':'NOT_ESTABLISHED','local_coupling_accuracy':'NOT_EVALUATED','root_iteration_trace_saved':False,'diagnostic_limit':'root nonconvergence traceback saved; no iteration history survived exception; do not attribute detailed numerical cause','route':'NO_GO_CURRENT_IMPLEMENTATION_AND_T2_PERFORMANCE_RESTART_UNDER_CURRENT_CONTRACT','next_priority':'research-lead review of T3 paired-node and broader observation/state constraints; NOT_AUTHORIZED_FOR_EXECUTION','history':'UNCHANGED','PREC3_B':'COMPLEX_K_NOT_CERTIFIED','PREC3_C':'NOT_EVALUATED','T2_A_B_actual_performance':'NOT_EVALUATED','CHIGH_and_real_bottom':'NOT_CERTIFIED','R4_percent':0,'research_lead_review':'PENDING','STOP':True}
 write(OUT/'FUNC2_FINAL_ADJUDICATION.json',final)
 write(OUT/'FUNC2_DATA_DICTIONARY.json',{'empty_science_CSV':'header-only means0 completed numeric rows, not zero-valued error or successful precision','planned_counts':'A24 endpoints=2methods*2precisions*3modes*2reportedlosses; B48 continuation points=2*2*3*4','A_attempted_failure':'first control; neither k nor function deviation quantified before root convergence','identity_checks':'160 checks are144 repository inputs+10 originalNPZ+6 frozen Python files, not160 independent physical samples','cold_PASS':'identity-only because no numeric reference arrays existed; k/function/boundary/highmodeQ cold scopes NOT_EVALUATED','historical_metrics':'all old numbers remain historical diagnostics, no new reference-based error or FE/FD ranking','runtime_admission':'first wrapper import lackedpsutil before execution marker; dependency added only onD; no formal scientific retry','ledger':'physicalCRLF and GitLF differ only by checkout text conversion; original bytes preserved before zero-credit append','root_trace_limit':'frozen solver kept unsuccessful Newton trace only in process memory; no converged snapshot or iteration values archived; no posthoc scientific rerun permitted'})
 stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
 route='''# T2数值路线：本阶段去留判断

结论：**NO-GO——停止当前冻结连续参考实现的后续执行，也不恢复T2-A/B接收或定位性能试验。建议把下一阶段主要研究资源转向T3辅助节点及观测—状态联合约束，但本轮不启动这些工作。**

这是一项阶段选路决定，不是原单HLA、倾斜阵或短垂直子阵物理不可行的证明。FUNC2在解析等声速无损波导的模式1控制中，DOP853未能在冻结的10次迭代内完成求根；因此只得到`FUNC2_ANALYTIC_REFERENCE_GATE_FAIL`。E-STD的481、961均未执行，不能把本轮写成“已证明两种高阶模态无法建立连续参考”。

## 1. 有限模态参考能否建立？

**本轮未建立。** 两条新路径虽已冻结实现，但主路径在第一个解析控制失败，独立Gauss路径尚未进行登记模态验证。没有任何模态获得`CONTINUOUS_COMPLEX_MODE_REFERENCE_SUPPORTED_SUBSET`。不能用预检中的代数控制，或旧FE/FD相互接近，替代新参考准入。

## 2. 200 m附近耦合精度有无可信的新证据？

**没有。** 本轮局部函数、源—阵元乘积Q和完整复场均没有新增数值结果。FUNC1的局部FE/FD差异和弱幅值诊断保持原作用域，不能改写成相对连续真函数的误差。共深度单模态向量容易共线，既有近零形状差异仍不是全场认证。

## 3. FE与FD谁更准确？

**当前无法判别。** 481两种方法之差呈二阶规律的历史观察继续保留；在缺少独立连续参考时，不能把其中任何一种指定为真值。961只保存FE完整向量，没有新增FD961，也没有补算授权。

## 4. 全模态与真实传播还缺什么？

有限模态连续参考与有损函数精度仍未准入；全模态相干叠加的误差控制未建立；N2的CHIGH窗口完整性问题独立存在。名义刚性5000 m海底、CVWT/Thorp和单条SSP只是冻结研究模型，不能视为真实海域已认证合同。提高波数实部精度不自动修复这些缺口。

## 5. 现有甲方资料下继续投入T2是否合理？

**不建议把主要资源继续投入这一数值修补链或硬件优选。** 这一轮事先被定义为有限模态函数精度专项的收口验证；解析准入未完成，且即便数值参考建立，真实环境合同仍是独立障碍。继续加网格、放宽阈值、追加求根方法或自动开启FUNC3，都不属于本轮授权。

这不是要求项目申请阶段先完成设备验收。应先把传播环境假设、可获得资料和不确定性范围写清；再由负责人判断T2的投入是否有相对于T3的科学价值。当前没有依据承诺倾斜或短垂直结构的距离精度。

## 6. 后续资源与解除条件

建议优先审查T3辅助节点已有配对观测及四参数、拒判和尾部风险问题，或研究对传播模型更稳健的观测及多维候选约束；下一任务需独立授权，本轮不运行。

若未来考虑恢复T2，应先明确拟研究海域、距离/深度范围及环境资料来源；需要可追溯的水深与地形、SSP测量/重建及其时空不确定性，以及海底声速、密度、衰减和必要的分层/粗糙度信息或有依据的先验范围。哪些参数需实测、哪些允许采用文献或区域先验，应由研究问题决定，不预设全部真实参数必须精确已知。随后才有理由另行评估有限模态参考、边界截断和全模态误差合同。

历史PREC3-B、N2及更早Gate全部不变；T2-A/B实际性能`NOT_EVALUATED`；R4=0%；提交推送后STOP。
'''
 (OUT/'T2_NUMERICAL_ROUTE_GO_NO_GO.md').write_text(route,encoding='utf-8')
 (OUT/'FUNC2_NEXT_RESEARCH_DECISION.md').write_text('''# 下一研究决定

本阶段关闭当前T2数值修补执行链，不自动开启FUNC3。建议负责人优先审查`T2_NUMERICAL_ROUTE_GO_NO_GO.md`，决定后续主资源是否转向T3辅助节点、传播观测鲁棒性或五维状态候选约束。

该建议未授权新实验。不能在本轮追加求根迭代、修改DOP853设置、替换控制模态、启动Gauss补跑、补FD961或恢复阵型性能计算。

本次失败首先说明冻结的连续求根实现没有通过解析准入；它没有回答原E-STD高阶参考是否可建立，也没有给出新的FE/FD真实误差或局部耦合精度。详细数值原因尚未得到保留下来的逐次迭代证据，不能事后补造。

R4=0%；STOP；等待研究负责人审阅。
''',encoding='utf-8')
 sync=f'''# T2-PREC3-FUNC2 同步报告

父提交：`{cfg['parent_SHA']}`。设计SHA：`{final['design_SHA']}`。执行SHA和remote/main在最终回传给出，避免自引用。

## 正式判定

`FUNC2_ANALYTIC_REFERENCE_GATE_FAIL`

**本轮没有建立新的连续复模态参考。** 正式主分析在第一个解析基准即停止：201 Hz、5000 m、等声速1500 m/s、压力释放海面及刚性海底，模式1、无损、DOP853较宽容差档，冻结求根最多10次迭代后抛出`SHOOT_ROOT_NOT_CONVERGED`。未修改代码或阈值重试；没有更换模态，也没有继续运行独立Gauss登记验证来绕过失败。

这是当前冻结求根实现的解析准入失败，不是原E-STD传播模型、481/961模态、倾斜阵或短垂直子阵的负物理证据。

## 实际完成范围

|环节|预登记|完成结果|失败／未执行|
|---|---:|---:|---|
|A解析端点控制|24|0个收敛结果|首个尝试失败1；其余23未评价|
|B E-STD延续根|48|0|Gate A失败，全部NOT_EVALUATED|
|B独立连续模态参考|3|0|模式1/481/961均未认证|
|C原FE/FD函数对参考|10|0|无参考，不能判FE或FD更准确|
|C 200 m误差带|250|0|未计算新的真实误差分布|
|D旧函数耦合对参考|200|0|未形成新局部精度证据|
|D参考间局部对照|180|0|未执行|
|冷审计身份检查|160|160通过|0 FAIL；作用域限身份|
|冷审计高阶函数/耦合重建|1个预注册高阶模态|0|没有获准进入的高阶根与函数|

A的24项是2方法×2精度×3模态×2个报告损耗端点；B另含4个损耗延续点。失败没有产生已收敛根或函数，因此最大波数偏差、最大函数形状误差、边界数值及局部Q误差均是**未获得**，不能填0。科学CSV保留表头而无数据行；完整完成数与原因见`FUNC2_STAGE_COMPLETION.csv`。

冻结参考预算：实部和虚部分别≤1.6666667e-9 m^-1，函数形状≤5e-4；双线性归一化误差≤1e-6，边界/匹配残差分别≤1e-8/2e-8。失败触发的是迭代收敛硬门，不是已经测得函数误差超过0.05%。旧0.2%门槛没有改变。

## 独立复核及证据边界

一次冷审计重新核对144份仓库输入、10份原始NPZ和6份冻结Python代码，共160项、0失败。由于主程序没有保存新的参考数组，根、函数、边界残差及至少一个高阶模态Q的数值冷重建全部未执行。原始`FUNC2_INDEPENDENT_COLD_AUDIT.json`中的PASS仅表示本次实际可用的身份检查通过；不能解释为数值参考通过或全部冷审计目标完成。

失败堆栈和原主日志已封存。冻结求根器在抛出异常前没有把逐次Newton历史写盘，因此本轮不能据此认定是迭代停滞、导数估计、容差还是另一项细节主导；不得补跑来填充该证据。`FUNC2_ROOT_TRACES.json`为空，明确作为交付限制保留。

## 输入完整性与启动记录

完整预检347项通过。原台账工作区CRLF与父Git blob的LF不同，统一行尾后内容逐字节一致；首次暂停记录和两种哈希保留，台账科学内容未变。所有正式科研代码均在Commit A前冻结。

第一次启动器导入时缺少psutil，尚未创建EXECUTION_STARTED、主分析和冷审计调用均为0。原pip自身导入损坏，随后仅在D盘下载并按发布SHA核对psutil5.9.8轮子，未改C盘环境或冻结科研代码。`FUNC2_RUNTIME_ADMISSION.json`记录这项启动依赖遗漏。正式主分析1次、冷审计1次、科研重试0次。

## 预算与禁止范围

科学链约{budget['wall_s']:.3f}秒，监控进程树峰值{budget['peak_tree_RSS_bytes']/1024**2:.2f} MiB；未触及3600秒/4 GiB硬门。全部新增数据位于D盘，主监控时新增数据{budget['new_D_bytes']}字节，Git结果{budget['git_result_bytes']}字节；最终交付体积见OUTPUT_MANIFEST。

新增KRAKEN、KRAKENC、FIELD、FIELDS、SCOOTER、原二阶FE/FD/Prüfer重算、全模态场、波形、MC、DFT、FIM、定位优化均为0。原10份函数未重归一化、平滑、替换或补算；FD961继续未归档、未评价。

## 去留判断

**NO-GO：停止当前冻结连续参考实现及本阶段T2数值修补链，不恢复T2-A/B性能试验。** 这不证明其他算法无效，也不证明481/961参考在预算内不可能建立。建议把后续主要研究资源优先转向T3及观测—状态联合约束，由负责人另行选题授权。继续T2还需先明确实际SSP与海底资料及其不确定性合同，而不只是再加网格。

完整理由、环境资料需求和解除条件见`T2_NUMERICAL_ROUTE_GO_NO_GO.md`。当前没有正面参考结果可供误差归因，不能沿用FUNC1数字冒充FUNC2新结果。

G0/N1/N2/N3/PREC1/PREC2、PREC3-A/B/C及FUNC1全部保持原判。CHIGH与真实海底未认证；T2-A/B实际性能NOT_EVALUATED；R4=0%。本轮交付推送后STOP，不自动启动FUNC3或T3。
'''
 (OUT/'GPT_SYNC_FUNC2.md').write_text(sync,encoding='utf-8')
 # This is a completion diagram, not a numerical error curve or scientific performance plot.
 fig,ax=plt.subplots(figsize=(10,3.5));ax.axis('off');labels=['A: analytic control\nFAILED: first attempted case','B: E-STD reference\nNOT EVALUATED','C: FE/FD attribution\nNOT EVALUATED','D: local coupling\nNOT EVALUATED'];colors=['#fee2e2','#f1f5f9','#f1f5f9','#f1f5f9']
 for i,(label,color) in enumerate(zip(labels,colors)):
  x=.02+i*.25;ax.text(x,.55,label,transform=ax.transAxes,ha='left',va='center',fontsize=9,bbox={'boxstyle':'round,pad=.7','facecolor':color,'edgecolor':'#94a3b8'})
 ax.text(.02,.1,'Execution-status diagram only. No converged root/function arrays were obtained. R4 = 0%.',transform=ax.transAxes,fontsize=10);ax.set_title('FUNC2: downstream science cancelled at the analytic gate',fontsize=13);fig.tight_layout();fig.savefig(OUT/'FUNC2_EXECUTION_SCOPE.png',dpi=180);fig.savefig(OUT/'FUNC2_EXECUTION_SCOPE.pdf');plt.close(fig)
 # Append exactly one zero-credit row; retain all old physical bytes as prefix.
 append=b'T2-PREC3-FUNC2,../R4_T2_PREC3_FUNC2_CONTINUOUS_MODAL_REFERENCE/FUNC2_FINAL_ADJUDICATION.json,FUNC2_ANALYTIC_REFERENCE_GATE_FAIL,A0_OF24_FIRST_ROOT_FAIL;B_C_D_NOT_EVALUATED;NO_NEW_EXTERNAL_SOLVER;R4_delta=0,COLD160_IDENTITY_PASS_0FAIL;NUMERIC_COLD_NOT_EVALUATED;HISTORICAL_UNCHANGED;RESEARCH_LEAD_PENDING;STOP,0\r\n'
 ledger.write_bytes(original+append);write(OUT/'LEDGER_APPEND_PROVENANCE.json',{'old_physical_sha256':hashlib.sha256(original).hexdigest(),'new_physical_sha256':sha(ledger),'old_bytes':len(original),'preserved_exact_prefix':ledger.read_bytes().startswith(original),'one_appended_row':True,'R4_delta':0})
 write(OUT/'FUNC2_DELIVERY_VALIDATION.json',{'source_identities':len(cfg['inputs']),'original_vectors':len(cfg['vectors']),'frozen_code_files':len(cfg['frozen_code']),'main_output_seal_preserved':True,'source_files_unchanged_before_ledger_append':True,'ledger_exact_prefix_preserved':True,'main_calls':1,'cold_calls':1,'scientific_retries':0,'science_admission':'FAIL','cold_scope':'IDENTITY_ONLY','numeric_reference_arrays':0,'R4_percent':0,'STOP':True,'UTC':stamp})
 (OUT/'FUNC2_FIGURE_SCOPE.md').write_text('执行状态图仅说明硬门与未执行范围，不是函数误差或阵型性能图。本轮没有数值参考数组，不能绘制或制造误差曲线。\n',encoding='utf-8')
 write(OUT/'OUTPUT_MANIFEST.json',{'files':[{'path':p.relative_to(OUT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'],'scope':'all delivery files except this manifest; main scientific outputs and frozen code unchanged','R4_percent':0,'STOP':True})
 print('DELIVERY_PREPARED',len(list(OUT.rglob('*'))),'files/data paths; formal science NOT retried')
if __name__=='__main__':main()
