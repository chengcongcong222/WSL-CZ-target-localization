"""Post-A report only; no numerical method changes or extra experiments."""
from hla_sp1_core import *
import itertools,subprocess

def main():
 design=read(OUT/'DESIGN_FREEZE.json');validation=read(OUT/'VALIDATION.json');summary=rows(OUT/'RANGE_RESULTS_BY_GEOMETRY.csv');comp=rows(OUT/'METHOD_AND_CALIBRATION_COMPARISON.csv');numeric=rows(OUT/'GRID_AND_NUMERIC_SENSITIVITY.csv');decision_methods={};classifications=[]
 for method in ['P0','P1','P2']:
  bygrid={}
  for step in [.10,.05]:
   a=[r for r in summary if r['method']==method and r['C']=='C0' and int(r['reference_SNR_dB'])==20 and float(r['grid_step_km'])==step and int(r['modal_mesh'])==160001];bygrid[str(step)]=[int(r['geometry']) for r in a if r['meets_signal']=='True']
  good=len(bygrid['0.05'])>=5 and len(bygrid['0.1'])>=5;decision_methods[method]={'signal':good,'geometries_by_grid':bygrid,'judgment':'CONDITIONAL_FINITE_CONTRACTION' if good else 'RANGE_CONTRACTION_NOT_SUPPORTED'}
 limited=validation['direct_boundary_flips']>0 or validation['direct_score_max_difference']>1e-7
 diag={}
 for method in ['P0','P1','P2']:
  a=[r for r in numeric if r['method']==method];gridwidth=[abs(float(r['fine_minus_main_width_km'])) for r in a if r['fine_minus_main_width_km']!='UNRESOLVED_EMPTY'];meshwidth=[abs(float(r['modal_width_difference_km'])) for r in a if r['modal_width_difference_km']!='UNRESOLVED_EMPTY'];diag[method]={'median_absolute_grid_envelope_change_km':float(np.median(gridwidth)),'median_absolute_modal_envelope_change_km':float(np.median(meshwidth)),'maximum_grid_point_change_km':max(abs(float(r['fine_minus_main_point_km'])) for r in a),'maximum_modal_point_change_km':max(abs(float(r['modal_point_difference_km'])) for r in a),'grid_horizontal_label_flips':sum(r['grid_horizontal_retention_flip']=='True' for r in a),'modal_horizontal_label_flips':sum(r['modal_horizontal_retention_flip']=='True' for r in a)};limited|=diag[method]['median_absolute_grid_envelope_change_km']>1 or diag[method]['median_absolute_modal_envelope_change_km']>1
  for step in [.10,.05]:
   hi=[r for r in summary if r['method']==method and r['C']=='C0' and int(r['reference_SNR_dB'])==20 and float(r['grid_step_km'])==step and int(r['modal_mesh'])==160001];lo=[r for r in summary if r['method']==method and r['C']=='C0' and int(r['reference_SNR_dB'])==20 and float(r['grid_step_km'])==step and int(r['modal_mesh'])==80001];limited|=(sum(r['meets_signal']=='True' for r in hi)>=5)!=(sum(r['meets_signal']=='True' for r in lo)>=5)
 # Task section8 independently requires reporting major ranking instability.
 # This post-execution audit flag does not modify the preregistered contraction signal.
 ranking_limited=True  # independent audit: observed14.15km/12.65km ranking jumps under task8; qualitative limitation, no added acceptance threshold
 limited=limited or ranking_limited
 if decision_methods['P1']['signal']:classifications.append('HLA_SP1_FINITE_RANGE_CONTRACTION_SUPPORTED_CONDITIONALLY')
 else:classifications.append('HLA_SP1_RANGE_CONTRACTION_NOT_SUPPORTED_IN_TESTED_CATALOGUE')
 classifications.append('HLA_SP1_PHASE_CALIBRATION_INVARIANCE_ONLY_OR_CONDITIONAL_GAIN')
 if limited:classifications.append('HLA_SP1_NUMERICAL_OR_GRID_LIMITED')
 if not validation['PASS']:classifications.append('HLA_SP1_INPUT_OR_IMPLEMENTATION_INCOMPLETE')
 decision={'classifications':classifications,'P1_preregistered_primary':decision_methods['P1'],'methods':decision_methods,'numeric_sensitivity':diag,'ranking_limit_scope':'best-point ranking only: observed grid jump14.15km and modal jump12.65km; envelope verdict consistent; post-execution audit under task section8, not a new performance threshold','phase_statement':'noise-free registered fixed cross-frequency gain invariance confirmed; no distance contraction or noisy dominance established','relative_to_P0':'per-method independently calibrated filters; not equal-coverage optimal-estimator comparison','scope':['STATIC_CONTROL','FINITE_CATALOGUE','DEPTH_NOT_TARGET','NO_VELOCITY_ESTIMATE'],'v_psi':'NOT_EVALUATED','R4_percent':0,'new_KRAKEN':6,'new_FIELD':3,'calibration_records':1024,'test_records':512,'test_method_grid_modal_configurations':6144,'new_dynamic_signals':0,'CHIGH_UNCERTIFIED':True,'old_E2_status_unchanged':True,'next':'STOP'};dump(OUT/'DECISION.json',decision)
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 # first registered geometry, first realization C0/20dB, full domain both catalogues
 public=rows(LOCAL/'test/PUBLIC_RECORD_INDEX.csv');first=next(i for i,r in enumerate(public) if r['C']=='C0' and int(r['reference_SNR_dB'])==20);tau=rows(OUT/'CALIBRATION_THRESHOLDS.csv');fig,axs=plt.subplots(1,2,figsize=(12,4))
 for ax,step in zip(axs,[.1,.05]):
  rs,ths,cat=catalogue(step);profile=np.load(LOCAL/f'RANGE_PROFILES_{step:.2f}_n160001.npy',mmap_mode='r')[first]
  for method in range(3):
   ax.plot(rs,profile[method],label='P'+str(method));t=float(next(r['tau'] for r in tau if r['C']=='C0' and r['reference_SNR_dB']=='20' and r['method']=='P'+str(method) and float(r['grid_step_km'])==step));ax.axhline(1-t,color=['C0','C1','C2'][method],linestyle='--',alpha=.5)
  ax.set_xlim(45,60);ax.set_ylim(0,1.02);ax.set_title(str(step)+' km catalogue; max over theta,z');ax.set_xlabel('Range km');ax.set_ylabel('Normalized score and empirical threshold');ax.legend()
 fig.suptitle('First registered geometry/realization, C0 20dB; no truth-local crop');fig.tight_layout();fig.savefig(OUT/'FIG_FIRST_FULL_RANGE_PROFILE.png',dpi=140);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(11,4))
 for method in ['P0','P1','P2']:
  a=sorted([r for r in summary if r['method']==method and r['C']=='C0' and r['reference_SNR_dB']=='20' and float(r['grid_step_km'])==.05 and r['modal_mesh']=='160001'],key=lambda r:int(r['geometry']));axs[0].plot(range(8),[float(r['median_envelope_km']) for r in a],'-o',label=method);axs[1].plot(range(8),[float(r['horizontal_label_retention']) for r in a],'-o',label=method)
 axs[0].set_ylim(0,16);axs[0].set_ylabel('Median full cell envelope km (empty=INF)');axs[1].set_ylim(0,1.05);axs[1].set_ylabel('Empirical horizontal label retention,16 repeats');
 for ax in axs:ax.set_xlabel('All geometry IDs');ax.legend()
 fig.suptitle('C0/20dB, fine catalogue, nominal160001 provider');fig.tight_layout();fig.savefig(OUT/'FIG_ALL_GEOMETRY_WIDTH_RETENTION.png',dpi=140);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(10,4))
 for method in ['P0','P1','P2']:
  a=[next(r for r in comp if r['method']==method and r['C']==c and r['reference_SNR_dB']=='20') for c in ['C0','C1']];axs[0].plot(['C0','C1'],[float(r['median_geometry_envelope_km']) for r in a],'-o',label=method);axs[1].plot(['C0','C1'],[float(r['median_geometry_point_error_km']) for r in a],'-o',label=method)
 axs[0].set_ylabel('Median of geometry median envelope km');axs[1].set_ylabel('Median of geometry median point error km');
 for ax in axs:ax.legend()
 fig.suptitle('20dB paired records; separate fixed empirical thresholds; point is not support');fig.tight_layout();fig.savefig(OUT/'FIG_CALIBRATION_METHOD_COMPARISON.png',dpi=140);plt.close(fig)
 lines=['| 方法/校准/参考SNR | 输出 | 水平标签保留 | 距离包络 km | 点误差 km | 全集合最坏 km |','|---|---:|---:|---:|---:|---:|']
 for r in comp:
  worst=float(r['maximum_cell_worst_error_km']);lines.append(f"| {r['method']}/{r['C']}/{r['reference_SNR_dB']}dB | {r['effective_outputs']}/128 | {100*float(r['horizontal_label_retention']):.2f}% | {float(r['median_geometry_envelope_km']):.2f} | {float(r['median_geometry_point_error_km']):.3f} | {worst:.3f} |")
 report=f'''# HLA-SP1同步报告

parent：{PARENT}。设计A：{read(OUT/'EXECUTION_STARTED.json')['design_SHA']}。执行B为包含本报告的提交，最终回报核对remote/main。任务SHA256：{design['task_sha256']}。

主要结论：{'; '.join(classifications)}。P1主条件C0/20dB，两套距离网格均{len(decision_methods['P1']['geometries_by_grid']['0.05'])}/8几何达到登记收紧信号；三方法均未建立距离总包络收紧。不能由最高峰较准推断全候选集合较窄，也不能上升为HLA物理不可观测。

以下细网格/160001结果：包络与点误差列为8个几何各16次中位数的中位，包含空集INF；最坏列包含拒判INF。水平标签保留允许全21深度搜索。联合生成深度标签、距离单元保留及逐例连续真值评分另见CSV，三者不可混用。

{chr(10).join(lines)}

C0/20dB三方法虽均保留128/128水平标签，包络中位仍覆盖完整45–60km域。完整不相连区间、占用长度和全部theta镜像/z标签保留于CSV及packed支持文件；15km总包络不意味着其中每点均被接受。5dB有空集，其宽度不是零，点/集合无条件误差计INF。全部512条测试、1024条校准记录完成，三方法共享数据；1536条波形均实际DFT，K片段和共享频对不是额外独立几何。

阈值由独立16几何×16次的最近标签D取第244顺序统计量，24个门槛在测试前冻结，未按测试调节。不同方法门槛不同，不构成理论95%置信域或等覆盖最优比较。源谱各频/片段未知。C1仅跨三频固定的±1dB/±30°通道复增益；P2无噪声不变性成立，但未形成距离收紧，也不优于所有正确未知增益估计器。

新静态KRAKEN6/6、FIELD3/3；FIELD共同源比例后最大形状残差{max(float(r['shape_relative_error']) for r in rows(OUT/'FIELD_CONTROLS.csv')):.4%}。生成采用160001完整复模态场，副本用去载波复样条，80001及两网格结果均保存。最佳点和区间边界直接复核{validation['all_best_and_boundary_checks']}项，最大得分差{validation['direct_score_max_difference']:.3g}，阈值翻转{validation['direct_boundary_flips']}；冷复核{validation['PASS_count']}PASS/{validation['FAIL_count']}FAIL。最佳点排序仍有大幅变化：距离网格变更最大跳14.15km，模态网格最大跳12.65km，故依任务第8节并列记录NUMERICAL_OR_GRID_LIMITED，仅限最佳点排序；总包络不收紧的方向两网格一致，未改门槛。上述证据不认证CHIGH遗漏模态，也不补发旧E2准入。

这是名义匹配环境的静态三频DFT格点控制：8元14m原HLA、阵深200m，8个独立2s片段，20/5dB为参考复谱SNR而非实际宽带SNR；弱场实际SNR单列未删除。没有盲频率、漂移、动态群延迟或1200s移动定位认证，没有新增VLA或辅助节点。深度仅搜索未知量，v/psi=NOT_EVALUATED。新数据和临时文件均在D。旧结论完整保留，R4=0%；提交推送后STOP，不加网格/频带或切换方法。
'''
 (OUT/'GPT_SYNC.md').write_text(report,encoding='utf-8');dump(OUT/'REPORT_BUDGET.json',{'wall_RSS':guard(),'repo_stage_bytes':sum(p.stat().st_size for p in OUT.iterdir() if p.is_file()),'local_bytes':sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file()),'scope':'single registered execution followed by serial independent review'})
 dump(OUT/'OUTPUT_MANIFEST.json',{'time_utc':stamp(),'parent':PARENT,'design_SHA':read(OUT/'EXECUTION_STARTED.json')['design_SHA'],'sha256':{p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'},'local_manifest':'public indexes and private archive manifests, DFT freeze and inference freeze bind all local big data','cold_audit_source_sha256':sha(ROOT/'hla_sp1_audit.py'),'report_source_sha256':sha(ROOT/'hla_sp1_report.py')});print(json.dumps(decision),flush=True)
if __name__=='__main__':main()
