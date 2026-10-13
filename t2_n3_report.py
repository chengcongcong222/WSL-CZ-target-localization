"""N3 scoped delivery and one append-only R4 entry; no scientific retries."""
from t2_n3_common import *
import shutil,subprocess,io

def main():
 seal=read(OUT/'SCIENCE_OUTPUT_FREEZE.json');audit=read(OUT/'COLD_AUDIT.json');physical=read(OUT/'PHYSICAL_SOURCE_AUDIT.json');calls=rows(OUT/'SOLVER_CALLS.csv');cross=rows(OUT/'ANALYTIC_KRAKEN_MODE_AND_FIELD_CROSSCHECK.csv');mesh=rows(OUT/'ANALYTIC_MESH_COMPARISON.csv');ext=rows(OUT/'E_STD_BOTTOM_EXTREMES.csv');failure=read(OUT/'RUN_FAILURES.json')
 passed=seal['analytic_configs_PASS'];maincalls=[r for r in calls if r['kind']=='ANALYTIC' and r['mesh']=='160001'];coarse=[r for r in calls if r['kind']=='ANALYTIC' and r['mesh']=='80001'];gate='PASS' if passed==16 and audit['PASS'] else 'FAIL' if any(r.get('analytic_gate_PASS')=='False' for r in calls) else 'NOT_EVALUATED'
 full=passed==16 and len(ext)==144 and audit['PASS'];maxfield=max((float(r['spatial_shape_error']) for r in cross),default=None);maxkr=max((float(r['max_kr_relative_error']) for r in cross),default=None)
 decision={'stage':'T2-N3','overall':'T2_N3_PHYSICS_AND_ANALYTIC_DIAGNOSTIC_COMPLETE_SCOPED' if full else 'T2_N3_PARTIAL_NUMERIC_NOT_ADMITTED','PHYSICAL_SOURCE_AND_TURNING_AUDIT':'COMPLETE','ANALYTIC_ISOVELOCITY_GATE':gate,'analytic_main':{'attempted':len(maincalls),'planned':12,'PASS':sum(r.get('analytic_gate_PASS')=='True' for r in maincalls)},'analytic_secondary':{'attempted':len(coarse),'planned':4,'PASS':sum(r.get('analytic_gate_PASS')=='True' for r in coarse)},'analytic_field_vectors':{'completed':len(cross),'planned':320,'max_shape_error':maxfield,'threshold':.002,'max_kr_relative_error':maxkr,'kr_threshold':2e-5},'IDEAL_RIGID_AND_PRESSURE_RELEASE_BOUNDARY_COMPARISON':'COMPLETE','ideal_comparison_scope':'SOLVER_FREE_ANALYTIC_REFERENCE_ONLY;16 count,240window and80radial rows independently checked; not evidence all KRAKEN configs executed','E_STD_BOTTOM_R_V_DIAGNOSTIC':'COMPLETE' if len(ext)==144 else 'PARTIAL' if ext else 'NOT_EVALUATED','E_STD_vectors':{'completed':len(ext),'planned':144,'max_shape_difference':max((float(r['spatial_shape_difference']) for r in ext),default=None)},'failure':failure,'cold_audit':audit,'new_calls':{'KRAKEN':len(calls),'SCOOTER':0,'FIELDS':0,'FIELD':0,'waveforms_DFT_MC_FIM_optimizers':0},'REAL_CZ1_PHYSICAL_PROPAGATION':'NOT_CERTIFIED','T2_A_B_RECEIVER_PERFORMANCE':'NOT_EVALUATED','historical_G0_N1_N2':'UNCHANGED','CHIGH_global_completeness':'NOT_CERTIFIED','real_seabed_parameters':'UNKNOWN','R4_PERCENT':0,'research_lead_audit':'PENDING','STOP':True,'next_stage_authorized':False}
 dump(OUT/'N3_DECISION.json',decision)
 archive=OUT/'modal_archive';archive.mkdir();indices=[]
 for p in sorted(LOCAL.rglob('*')):
  if not p.is_file():continue
  upload=p.is_relative_to(LOCAL/'runs') and p.suffix.lower() in ('.mod','.prt','.env','.txt')
  public=None
  if upload:
   public=archive/p.parent.name/p.name;public.parent.mkdir(exist_ok=True);shutil.copyfile(p,public);assert sha(p)==sha(public)
  indices.append({'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p),'public_copy':str(public) if public else 'LOCAL_ONLY_SHA_INDEXED'})
 table(OUT/'LOCAL_ARCHIVE_INDEX.csv',indices)
 expected={r['path']:r['sha256'] for r in read(OUT/'DESIGN_FREEZE.json')['frozen_files']};unchanged=all(sha(p)==h for p,h in expected.items())
 validation={'PASS':audit['PASS'] and unchanged,'cold_reconstruction_PASS':audit['PASS'],'scientific_gate_PASS':gate=='PASS','all_frozen_code_inputs_unchanged':unchanged,'no_second_solver_attempt':True,'stop_after_first_gate_failure':True,'new_field_or_signal_calls':0,'completed_analytic_configurations':len(maincalls)+len(coarse),'not_evaluated_E_STD_is_not_PASS':True}
 dump(OUT/'VALIDATION.json',validation)
 ledger=ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv';before=(LOCAL/'LEDGER_ORIGINAL.csv').read_bytes();assert ledger.read_bytes()==before
 txt=io.StringIO();w=csv.writer(txt,lineterminator='\n');w.writerow(['T2-N3','../R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE/N3_DECISION.json',decision['overall'],f'OLD_PHYSICS_READ_ONLY;ANALYTIC_{passed}_OF16_PASS;KRAKEN_{len(calls)};E_STD_BOTTOM_{len(ext)}_OF144;NO_RECEIVER_MC;R4_delta=0',f'COLD_{audit["checks"]}_ROWS_{audit["FAIL"]}_FAIL;RESEARCH_LEAD_PENDING;STOP',0]);ledger.write_bytes(before+(b'' if before.endswith(b'\n') else b'\n')+txt.getvalue().encode('utf-8'))
 dump(OUT/'LEDGER_APPEND_PROVENANCE.json',{'parent_SHA':PARENT,'original_sha256':hashlib.sha256(before).hexdigest(),'parent_blob_sha256':read(OUT/'LEDGER_BEFORE.json')['parent_blob_sha256'],'working_prefix_preserved':ledger.read_bytes().startswith(before),'new_rows':1,'after_sha256':sha(ledger),'R4_delta':0})
 lines=['# T2-N3 同步报告','',f'设计SHA：`{seal["design_SHA"]}`。执行SHA由最终回传列示，避免文件自引用。',f'任务SHA256：`{TASKSHA}`。','',f'本轮判定：`{decision["overall"]}`。','', '|阶段|完成范围与结论|','|---|---|',f'|N3-1 原物理合同|COMPLETE；{physical["modal_records"]}条存档模态分类，{physical["contribution_records"]}条复相干带贡献记录|',f'|N3-2 主解析|{sum(r.get("analytic_gate_PASS")=="True" for r in maincalls)}/{len(maincalls)}已执行配置通过；预登记12项|',f'|N3-2 粗网格|{sum(r.get("analytic_gate_PASS")=="True" for r in coarse)}/{len(coarse)}已执行；预登记4项，零执行不称通过|',f'|N3-3 E-STD底界|{len(ext)}/144向量；{decision["E_STD_BOTTOM_R_V_DIAGNOSTIC"]}|','',f'原SSP逐字节核验：水深5000 m，自由海面V、刚性海底R、CVWT和实际PRT的Thorp；c_min={physical["c_min"]} m/s，c_max={physical["c_max"]} m/s，声道轴1300 m。水柱密度1 g/cm³。所有新增1800–2400 m/s模态的相速度高于水柱最大声速，其kz²在全部原SSP节点为正。此为名义水柱无内部转向的判定，不是实测海底回波证书。','', '模态带在复数域完整相加并独立重建；相消和弱场均保留。各带范数之比可能大于1，不能称为独立声线能量比例。','',f'独立解析参考已复现并冷核对：16条计数、240条窗口差异、80条Hankel/远场差异。它们描述两个等声速理想波导，不是原Munk型SSP或真实第一会聚区。','']
 if cross:
  worst=max(cross,key=lambda r:float(r['spatial_shape_error']));lines += [f'已执行KRAKEN解析场最大空间形状差异：{maxfield*100:.8f}%（门槛0.2%），最大水平波数相对偏差{maxkr:.10g}（门槛2e-5）。最坏：{worst["stem"]}、{worst["probe"]}、{worst["resource"]}。逐配置模态数量、原始复场差、全五探针四资源共用一个复比例的诊断全部保存。','']
 if failure:lines += ['冻结停止原因：`'+json.dumps(failure,ensure_ascii=False)+'`。没有重试、修改环境、增加预算或改门槛；后续配置明确标记NOT_EVALUATED。','']
 if ext:lines += [f'R/V数学边界对照最大空间形状差异{decision["E_STD_vectors"]["max_shape_difference"]*100:.8f}%。两种边界只是数学极端，不代表真实海底选择。','']
 lines += [f'冷复核：{audit["checks"]}项，{audit["FAIL"]} FAIL，PASS={audit["PASS"]}。使用独立MOD解析、本征根brentq、J0-iY0径向项及QR投影。冷复核通过只认证结果和Gate记录可重建，不能改变科学Gate失败。','',f'新增KRAKEN {len(calls)}/19；SCOOTER/FIELDS/FIELD、波形、DFT、MC、FIM、优化、N1模板重评分均0。求解期壁钟{seal["budget"]["wall_s"]:.2f} s，进程树峰RSS {seal["peak_process_tree_RSS_bytes"]/1024**2:.2f} MiB；科学全链预算另见PIPELINE_BUDGET.json。','', '原G0/N1/N2全部冻结；N1阵型证据仅限旧1800窗口，不移植到本轮。CHIGH全域完整性未认证。T2-A/B接收和定位性能NOT_EVALUATED，真实CZ1传播NOT_CERTIFIED，实际底质UNKNOWN，R4=0%。','', '后续建议：若解析Gate失败，先研究本次已保存的模型/工具数值合同差异；若全部通过且底界敏感，则须获得有来源的底质与传播合同，再讨论阵型优选。任何后续工作均需独立任务，本轮不自动执行。','', '负责人审计PENDING。STOP。', '', '工具边界定义参考：[官方KRAKEN手册](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node47.html)、[边界条件](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node10.html)。']
 (OUT/'GPT_SYNC.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 (OUT/'NEXT_RESEARCH_DECISION.md').write_text('本轮不授权下一试验。\n\n'+('\n'.join(lines[-9:]))+'\n',encoding='utf-8')
 dump(OUT/'STOP.json',{'STOP':True,'R4':0,'next_stage_authorized':False,'receiver_performance':'NOT_EVALUATED'})
 try:
  import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
  fig,ax=plt.subplots(1,2,figsize=(11,4));ssp=env_data(mod_path('K0',201).with_suffix('.env'))['ssp'];ax[0].plot(ssp[:,1],ssp[:,0]);ax[0].invert_yaxis();ax[0].set(xlabel='Sound speed (m/s)',ylabel='Depth (m)',title='Archived E-STD: rigid bottom, Thorp')
  summaries=[r for r in physical['mode_summary'] if r['variant']=='K2'];bottom=np.zeros(3)
  for j,band in enumerate(('(1500,cmax]','(cmax,1800]','(1800,2000]','(2000,2400]')):
   n=np.array([r['bands'][j] for r in summaries]);ax[1].bar([str(r['frequency_hz']) for r in summaries],n,bottom=bottom,label=band);bottom+=n
  ax[1].set(xlabel='Frequency (Hz)',ylabel='Returned modes',title='K2: phase-speed groups, not ray energy');ax[1].legend(fontsize=8);fig.tight_layout();(OUT/'figures').mkdir();fig.savefig(OUT/'figures/PHYSICAL_CONTRACT_AND_MODAL_GROUPS.png',dpi=160);plt.close(fig)
 except ImportError:(OUT/'FIGURE_NOT_READY.txt').write_text('Matplotlib unavailable; numeric tables complete.\n',encoding='utf-8')
 print('DELIVERY',decision['overall'],'STOP',flush=True)
if __name__=='__main__':main()