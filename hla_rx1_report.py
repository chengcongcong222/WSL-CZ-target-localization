"""Delivery only, no alteration of the frozen experiment."""
import datetime,json,csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hla_rx1_common import *
def main():
 cfg=read(OUT/'DESIGN_FREEZE.json');status=read(OUT/'EXECUTION_STATUS.json');man=read(OUT/'RECEIVER_DATA_MANIFEST.json');decision=read(OUT/'DECISION.json');val=read(OUT/'VALIDATION.json');cold=read(OUT/'SIGNAL_PROVENANCE_VALIDATION.json')
 provider=pd.read_csv(OUT/'PROVIDER_GROUP_AND_BANDWIDTH_CHECK.csv')
 provider_failed=not bool(provider.drift_ready.all())
 decision['classifications']=list(dict.fromkeys(decision['classifications']+(['IMPLEMENTATION_OR_PROVIDER_LIMIT'] if provider_failed or not val['PASS'] or not cold['PASS'] else [])))
 decision['independent_arithmetic_review_PASS']=val['PASS'];decision['independent_signal_review_PASS']=cold['PASS'];decision['provider_group_delay_certified']=not provider_failed
 decision['stable_delay_fallback']='1500 m/s for failed group provider; dynamic modal delay not certified' if provider_failed else 'finite-difference modal group delay'
 conditions=[]
 for family in ['STABLE','DRIFT']:
  for snr in [20,5]:
   count=sum(r['source_family']==family and r['SNR_dB']==snr for r in man['records'])
   conditions.append(dict(source_family=family,SNR_dB=snr,registered_records=64,completed_records=count,classification='IMPLEMENTATION_OR_PROVIDER_LIMIT' if provider_failed else 'FEATURE_MAPPING_NOT_ESTABLISHED',mapping='NOT_ESTABLISHED' if count else 'NOT_EVALUATED'))
 decision['per_condition']=conditions;decision['completed_stable_reference_dependence']='REF/RX mapping diagnostic differs strongly under the unvalidated delay fallback; reference dependent observation only, no certified CZ capability'
 decision['singlepath_scope_only']='REFERENCE_DEPENDENT_ONLY: REF retains absolute motion; receiver-estimated carrier removes mean Doppler, not a general No-Go'
 dump(OUT/'DECISION.json',decision);table(OUT/'CONDITION_COMPLETENESS.csv',conditions)
 stats=pd.read_csv(OUT/'FEATURE_ERROR_AND_COVERAGE.csv',float_precision='round_trip');m=pd.read_csv(OUT/'WINDOW_FEATURE_MAPPING.csv',float_precision='round_trip');primary=stats[(stats.source_family=='STABLE')&(stats.SNR_dB==20)&(stats.lag_s==200)]
 scene=[]
 for key,g in m[m.lag_s==200].groupby(['geometry','source_family','SNR_dB','method']):
  valid=g[g.accepted];err=np.sort(valid.qeff_vs_abs_qbar_error_mps.to_numpy());un=np.sort(np.where(g.accepted,g.qeff_vs_abs_qbar_error_mps,np.inf))
  scene.append(dict(geometry=int(key[0]),source_family=key[1],SNR_dB=int(key[2]),method=key[3],registered_records=8,completed_records=g.record_id.nunique(),registered_line_windows=144,accepted_line_windows=len(valid),conditional_median_mps=float(err[int(np.ceil(len(err)*.5))-1]) if len(err) else float('inf'),conditional_P95_mps=float(err[int(np.ceil(len(err)*.95))-1]) if len(err) else float('inf'),reject_INF_P95_mps=float(un[int(np.ceil(len(un)*.95))-1]),fraction_registered_error_le_0_05=float(sum(un<=.05)/144)))
 table(OUT/'PER_GEOMETRY_CONDITIONS.csv',scene)
 summary='；'.join(f"{r['method']}输出{int(r['accepted_line_windows'])}/1152，映射诊断中位误差{r['conditional_median_abs_error_mps']:.4f}、P95 {r['conditional_P95_abs_error_mps']:.4f} m/s，误差≤0.05占登记线窗{r['fraction_registered_with_error_le_0.05_mps']:.1%}" for r in primary.to_dict('records'))
 rec=man['records'][0];x=np.load(rec['path']);s=np.load(rec['correlation_spectra_path']);key='w0_f1_RX_L200'
 fig,ax=plt.subplots(1,3,figsize=(12,3.5));fx=np.fft.fftshift(np.fft.fftfreq(800,.25));px=np.fft.fftshift(abs(np.fft.fft(x[:800,1,0]*np.hanning(800)))**2);ax[0].plot(fx+235,10*np.log10(np.maximum(px/max(px),1e-15)));ax[0].set(title='First registered actual receive spectrum',xlabel='Hz',xlim=(234.5,235.5));ax[1].plot(np.arange(len(s[key+'_correlation'])),s[key+'_correlation']);ax[1].set(title='Actual RX lag correlation',xlabel='lag s');ax[2].plot(s[key+'_frequencies'],s[key+'_spectrum']);ax[2].set(xlim=(0,.15),title='RX spectrum: include DC',xlabel='Hz');fig.tight_layout();fig.savefig(OUT/'FIG1_RECEIVER_AND_CORRELATION.png',dpi=130);plt.close(fig);s.close()
 a=stats[stats.lag_s==200];fig,ax=plt.subplots(figsize=(8,4));ix=np.arange(len(a));ax.bar(ix,a.conditional_median_abs_error_mps,label='median diagnostic error');ax.scatter(ix,a.conditional_P95_abs_error_mps,color='red',label='P95 diagnostic error');ax.set_xticks(ix,[f'{r.source_family} {r.SNR_dB}dB {r.method}' for r in a.itertuples()],rotation=25,ha='right');ax.set(ylabel='|q_eff - |qbar|| m/s',title='Conditional mapping diagnostic; DRIFT not evaluated');ax.legend();fig.tight_layout();fig.savefig(OUT/'FIG2_ERROR_BY_CONDITION.png',dpi=130);plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4))
 for method in ['RX','REF']:
  a=m[(m.source_family=='STABLE')&(m.SNR_dB==20)&(m.lag_s==200)&(m.line==1)&(m.method==method)];ax.scatter(abs(a.qbar_mps),a.q_eff_mps,s=9,alpha=.3,label=method)
 ax.plot([0,.8],[0,.8],'k--');ax.set(xlabel='H2 |net range displacement / 200s|',ylabel='uncertified correlation q_eff',title='All eight geometries, 235 band');ax.legend();fig.tight_layout();fig.savefig(OUT/'FIG3_REFERENCE_AND_MAPPING.png',dpi=130);plt.close(fig)
 bad=m[(m.qeff_vs_abs_qbar_error_mps>.2)&(m.method=='RX')&(m.lag_s==200)]
 if len(bad):
  first=bad.iloc[0];r=next(r for r in man['records'] if r['record_id']==first.record_id);j=int(first.window);np.savez(OUT/'RECOMPUTABLE_FIRST_MAPPING_FAILURE.npz',samples=np.load(r['path'])[j*800:(j+1)*800],record_id=r['record_id'],window=j,note='first registered mapping error >0.2 m/s; not necessarily line rejection')
 sync=f"""# HLA-RX1同步报告

parent：{PARENT}。设计A：{read(OUT/'EXECUTION_STARTED.json')['design_sha']}。执行B和remote/main以本次最终回报及提交历史为准，文件不内嵌自身提交SHA。

本轮只执行一次冻结配置，KRAKEN {status['new_KRAKEN_calls']}/30登记调用（上限36）、FIELD=0。保存{status['receiver_records_completed']}/256条原8元/14m单HLA接收记录：STABLE全部128条完成，DRIFT全部128条PROVIDER_NOT_READY。不是q+Gaussian注入，也不是H直接进入提取器；信号包含非整源载频、随机源幅/初相位、传播相位和带限通道噪声。1200s复基带、完整相关及滞后谱均留在D:/ProjectStorage/WSL-CZ/HLA_RX1/receiver_records，路径/SHA见RECEIVER_DATA_MANIFEST.json。

结论并列：FEATURE_MAPPING_NOT_ESTABLISHED、IMPLEMENTATION_OR_PROVIDER_LIMIT、PARTIAL_EXECUTION。三个频带的双网格全模态群延迟匹配均未通过，未删模态、未缩步长、未追加求解，故没有产生漂移接收数据。STABLE沿冻结代码使用1500m/s模态延迟后备假设，因此也不能认证动态模态物理映射；静态精确中心频率完整复场仍可复核。这是提供器/实现局限，不是漂移物理无效。

主条件STABLE/20dB：{summary}。数字只表示q_eff与H2的|qbar|之间的诊断差异，不是认证径向精度或Gaussian sigma。其余噪声档、100s对照、逐场景分母和全部峰/拒判保留于CSV；六窗三线不是额外独立航迹，每几何仅8个配对源/噪声实现。

五项实际接收控制通过：单路径延迟/相位/单位、非整频率双数字参考、源频—运动同接收记录歧义、4/8Hz采样与实值/基带等价。单路径REF保留已知源参考下的运动；RX估计接收载频时同时减去了未知源频偏和平均Doppler。这个零点限制是本轮处理选择暴露的结果，不能把REF对照借给原被动任务。所有输出保留±方向与多峰，不用真值挑峰。

REF只额外知道源频和接收时刻源相位，不使用真实发射延迟，不是ORACLE_ALIGNMENT。相关多模态峰没有认证单一径向时间权、有效相速度或未知源零点；因此不能直接改名为H2的200s净位移率。中点速度、净位移率与平均绝对速度的差异分别保存。没有方位提取或方位噪声输入，与方位误差依赖未评价。

冷复核PASS={val['PASS']}：{val['checks']}项、{val['FAIL']} FAIL；独立接收/发射时刻重建PASS={cold['PASS']}。检查通过只证明保存信号与算术可重建，不证明动态传播合同和量测映射通过。设计说明的中文曾被PowerShell ASCII管道破坏，原冻结文件/哈希保留，补交METHOD_INPUT_LOCK_READABLE.md；代码、参数与原任务副本未改。

严格单计算进程协议还有一项偏差：两次约1秒的存档模式诊断曾与主生成并行，未新增求解或改变配置，见BUDGET_PROTOCOL_DEVIATIONS.json；不能将执行协议整体标为PASS。

本轮不具备回灌H2条件。不拟合Gaussian替换实际误差，不真值去偏，不追加提供器或提取器修补。原H1/H2/E1/H3结论保持，R4=0%；B提交推送并核对远端后STOP。
"""
 (OUT/'GPT_SYNC.md').write_text(sync,encoding='utf-8')
 # Original frozen evaluator's Chinese table was also sent through ASCII stdin;
 # preserve its output and provide an explicitly readable delivery table.
 lines=['# 条件表（可读交付）','','| 条件 | 输出/登记线窗 | 中位映射误差 m/s | P95 m/s | ≤0.05比例 |','|---|---:|---:|---:|---:|']
 for row in stats[stats.lag_s==200].to_dict('records'):
  lines.append(f"| {row['source_family']}/{row['SNR_dB']}dB/{row['method']} | {int(row['accepted_line_windows'])}/1152 | {row['conditional_median_abs_error_mps']:.4f} | {row['conditional_P95_abs_error_mps']:.4f} | {row['fraction_registered_with_error_le_0.05_mps']:.1%} |")
 lines+=['| DRIFT/两SNR/REF与RX | 0/4608；提供器未就绪 | NOT_EVALUATED | NOT_EVALUATED | NOT_EVALUATED |','','STABLE也使用1500m/s后备延迟，全部表值只为未认证动态模型的映射诊断。源频、源相位、深度、模态与真运动均不交RX；公开LO、精确导航、固定端射合并与接收处声速为条件。无真实UUV源谱、海噪、导航/阵元校准认证。']
 (OUT/'MEASUREMENT_CONDITION_TABLE_READABLE.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 (OUT/'MEASUREMENT_CONDITION_TABLE_INITIAL_ENCODED.md').write_bytes((OUT/'MEASUREMENT_CONDITION_TABLE.md').read_bytes())
 (OUT/'MEASUREMENT_CONDITION_TABLE.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 (OUT/'BUDGET_PROTOCOL_DEVIATIONS.json').write_text(json.dumps({'frozen_compute_process_limit':1,'temporary_saved_modal_diagnostic_concurrent_with_main_worker':True,'observed_python_process_count_during_diagnostic':2,'diagnostic_commands_count':2,'combined_diagnostic_compute_seconds_approximately':2,'additional_solver_calls':0,'effect':'arithmetic inspection only; no change to configuration or receiver generation','interpretation':'strict single-process execution protocol was not fully met; report as a deviation, not full protocol PASS'},indent=2)+'\n',encoding='utf-8')
 total=(datetime.datetime.now(datetime.timezone.utc)-datetime.datetime.fromisoformat(read(OUT/'EXECUTION_STARTED.json')['time_utc'])).total_seconds()
 dump(OUT/'BUDGET_ACCOUNTING.json',{'elapsed_since_execution_start_s':total,'limit_s':14400,'local_data_bytes':sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file()),'repo_new_stage_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),'RSS_limit_bytes':8*1024**3,'note':'single computational worker; execution/cold guards; no claimed continuous peak RSS measurement'})
 dump(OUT/'OUTPUT_MANIFEST.json',{'artifacts_sha256':{p.relative_to(ROOT).as_posix():sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='OUTPUT_MANIFEST.json'},'code_sha256':{p.name:sha(p) for p in ROOT.glob('hla_rx1_*.py')},'large_records':'D only; no copyrighted papers committed'})
 print('DELIVERY',decision['classifications'])
if __name__=='__main__':main()
