from hla_rx1z_common import *
import math,subprocess

def main():
 guard();checks=[]
 def check(name,ok):checks.append({'name':name,'PASS':bool(ok)})
 design=read(OUT/'DESIGN_FREEZE.json');manifest=read(OUT/'FEATURE_OUTPUT_MANIFEST.json');start=read(OUT/'EVALUATION_STARTED.json')
 check('freeze_before_truth',manifest['time_utc']<start['time_utc'])
 for n,s in manifest['feature_sha256'].items():check('feature_immutable:'+n,sha(OUT/n)==s)
 for n,s in design['code_sha256'].items():check('code_immutable:'+n,sha(ROOT/n)==s)
 for n,s in design['input_sha256'].items():check('input_immutable:'+n,sha(OUT/n)==s)
 for r in rows(OUT/'INPUT_PROVENANCE.csv'):check('historical_immutable:'+r['path'],sha(ROOT/r['path'])==r['working_sha256'])
 for r in rows(OUT/'PUBLIC_RECORD_INDEX.csv'):check('raw_immutable:'+r['record_id'],sha(r['path'])==r['sha256'])
 ledger={(r['record_id'],int(r['line']),int(r['window'])):r for r in rows(OUT/'RECEIVER_FREQUENCY_LEDGER.csv')};features=rows(OUT/'INTERWINDOW_FEATURES.csv');ev=rows(OUT/'MAPPING_DIAGNOSTIC_ALL.csv');truth=rows(ROOT/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv');meta=read(OUT/'PUBLIC_RECEIVER_METADATA.json')
 for f,r in zip(features,ev):
  j=int(f['difference_j']);key=(f['record_id'],int(f['line']));a=ledger[key+(j,)];b=ledger[key+(j+1,)];d=float(b['f_hat_hz'])-float(a['f_hat_hz']);mean=sum(float(ledger[key+(k,)]['f_hat_hz']) for k in range(6))/6;v=-meta['c_ref_mps']*d/mean
  check('signed_difference:'+str(key)+':'+str(j),abs(d-float(f['d_frequency_hz']))<1e-12 and abs(v-float(f['apparent_radial_change_mps']))<1e-12)
  h=truth[int(r['geometry'])];angle=math.radians(float(h['theta_deg']));course=math.radians(float(h['psi_deg']));r0=float(h['r_km'])*1000;speed=float(h['v_mps']);positions=[]
  for k in range(7):
   t=k*200;pos=complex(r0*math.cos(angle),r0*math.sin(angle))+speed*t*complex(math.cos(course),math.sin(course));platform=complex(2*min(t,600),0)+2*max(t-600,0)*complex(math.cos(math.pi/12),math.sin(math.pi/12));positions.append(abs(pos-platform))
  dq=(positions[j+2]-2*positions[j+1]+positions[j])/200
  check('cold_geometry:'+str(key)+':'+str(j),abs(dq-float(r['true_delta_qbar_mps']))<1e-10 and abs(v-dq-float(r['error_mps']))<1e-10)
 check('counts',len(ledger)==2304 and len(features)==1920 and len(ev)==1920)
 check('recompute',manifest['max_frequency_recompute_difference_hz']<=1e-8)
 check('digital_reference',manifest['digital_reference_PASS'])
 summary=rows(OUT/'RELATIVE_MAPPING_BY_GEOMETRY.csv')
 for row in summary:
  sub=[r for r in ev if r['geometry']==row['geometry'] and r['SNR_dB']==row['SNR_dB'] and r['line']==row['line'] and (row['scope']=='ALL' or r['group']==row['scope'])]
  er=[float(r['error_mps']) for r in sub];q=[float(r['true_delta_qbar_mps']) for r in sub];rm=math.sqrt(sum(x*x for x in er)/len(er));zr=math.sqrt(sum(x*x for x in q)/len(q))
  check('cold_summary:'+str((row['geometry'],row['SNR_dB'],row['line'],row['scope'])),len(sub)==int(row['planned']) and abs(rm-float(row['RMSE_mps']))<1e-12 and abs(zr-float(row['zero_predictor_RMSE_mps']))<1e-12)
 passing={scope:[] for scope in ['SAME_LEG','TURN']}
 for scope in passing:
  for g in range(8):
   if sum(r['beats_zero']=='True' for r in summary if int(r['geometry'])==g and int(r['SNR_dB'])==20 and r['scope']==scope)>=2:passing[scope].append(g)
 common=set(passing['SAME_LEG'])&set(passing['TURN']);mapping='SUPPORTED_IN_SAVED_MODEL_DIAGNOSTIC' if len(common)>=5 else ('NOT_SUPPORTED_IN_SAVED_MODEL' if not passing['SAME_LEG'] else 'UNRESOLVED')
 valid=all(r['PASS'] for r in checks);decision={'REFERENCE_INFORMATION_ACCOUNTED':valid,'INTERWINDOW_FREQUENCY_FEATURE':'REPRODUCIBLE' if valid else 'TRACKING_UNSTABLE','RELATIVE_MOTION_MAPPING':mapping if valid else 'UNRESOLVED','same_leg_geometry_pass_ids':passing['SAME_LEG'],'turn_geometry_pass_ids':passing['TURN'],'main_SNR_dB':20,'ABSOLUTE_RADIAL_ZERO_NOT_ESTABLISHED':True,'FIRST_CZ_DYNAMIC_PROVIDER_UNCERTIFIED_FROM_RX1':True,'DRIFT_NOT_EVALUATED':True,'H2_FEEDBACK_NOT_AUTHORIZED':True,'R4_percent':0,'new_KRAKEN_FIELD':0,'new_signal_noise_state_optimization':0,'status':'COMPLETED' if valid else 'IMPLEMENTATION_LIMIT','stop':True}
 dump(OUT/'DECISION.json',decision)
 dump(OUT/'VALIDATION.json',{'time_utc':stamp(),'all_PASS':valid,'checks':checks,'PASS_count':sum(r['PASS'] for r in checks),'FAIL_count':sum(not r['PASS'] for r in checks),'scope':'independent scalar geometry and feature arithmetic, byte immutability; not external physical mapping certification','resource_sample':guard()})
 # Two predetermined figures; no fitted curves or success filtering.
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 first=sorted({r['record_id'] for r in rows(OUT/'PUBLIC_RECORD_INDEX.csv') if int(r['SNR_dB'])==20})[0]
 fig,axs=plt.subplots(2,1,figsize=(9,6))
 for line in range(3):
  a=[ledger[(first,line,k)] for k in range(6)];axs[0].plot([float(r['time_center_s']) for r in a],[float(r['center_offset_hz']) for r in a],'-o',label=str(meta['digital_reference_hz'][line])+' Hz');axs[1].plot([float(r['time_center_s']) for r in a],[float(r['old_residual_nu_hz']) for r in a],'-o')
 axs[0].set_ylabel('Receive center minus public LO (Hz)');axs[1].set_ylabel('Old residual peak (Hz)');axs[1].set_xlabel('Window center time (s)');axs[0].legend();fig.suptitle(first+' (preselected first 20 dB record)');fig.tight_layout();fig.savefig(OUT/'FIRST_RECORD_FREQUENCY.png',dpi=140);plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(11,4))
 for ax,scope in zip(axs,['SAME_LEG','TURN']):
  for line in range(3):
   a=sorted([r for r in summary if int(r['SNR_dB'])==20 and int(r['line'])==line and r['scope']==scope],key=lambda r:int(r['geometry']));ax.plot(range(8),[float(r['RMSE_mps']) for r in a],'-o',label=str(meta['digital_reference_hz'][line])+' Hz error')
  ax.plot(range(8),[float(r['zero_predictor_RMSE_mps']) for r in a],'k--',label='True change RMS / zero predictor error');ax.set_title(scope);ax.set_xlabel('All geometry IDs');ax.set_ylabel('m/s (mapping diagnostic)');ax.legend(fontsize=7)
 fig.tight_layout();fig.savefig(OUT/'ALL_GEOMETRY_SIGNAL_ERROR.png',dpi=140);plt.close(fig)
 pooled=rows(OUT/'POOLED_DESCRIPTIVE_ONLY.csv');lines=[]
 for r in pooled:lines.append(f"- {r['SNR_dB']} dB / {r['scope']}: RMSE={float(r['RMSE_mps']):.6f} m/s，零变化RMSE={float(r['zero_predictor_RMSE_mps']):.6f} m/s，比值={float(r['RMSE_to_zero_ratio']):.3f}，有效={r['valid']}/{r['planned']}。")
 text=f'''# HLA-RX1Z 同步报告

新任务SHA256：{design['task_sha256']}。
Parent：{PARENT}。设计SHA：{read(OUT/'EXECUTION_STARTED.json')['design_sha']}。
执行SHA为包含本文件及结果的Commit B，见最终回报；该提交不自引用哈希。

旧RX1已经完成128条稳定源记录，漂移记录0/128。稳定源使用1500 m/s群延迟后备，动态映射未认证。旧20 dB RX绝对径向映射诊断中位差0.4380 m/s，REF仅0.0038 m/s，表明参考依赖；两者均不能赋予未知源主任务绝对速度精度。

本轮没有补跑RX1，只承接128条记录（两档SNR各64），保留2304个接收载频值，产生1920个有符号相邻窗差分。接收侧先独立冻结输出哈希，再由评价进程读取真值；禁止源频率、源相位、REF和真值挑峰。128条原始记录各同设置复算一次，最大载频差{manifest['max_frequency_recompute_difference_hz']:.3g} Hz；两项数字参考坐标控制保持通过。旧逐窗去中心会消去常量运动频移；把接收中心加回仍不能区分未知源频率与运动零点。

主表示仍为Hz。以下m/s仅为-c_ref Delta F/Fbar表观变化，比较六个200 s距离增量的相邻差，绝不是绝对径向幅值或H2独立Gaussian特征：

{chr(10).join(lines)}

跨600 s转向项与同航段四项独立列报。主20 dB同航段至少两频优于零预测的几何为{passing['SAME_LEG']}；转向项为{passing['TURN']}。上述是描述性路线判据，8次重复与共享三频/相邻差分不构成总体覆盖率证明。实际STFT首尾裁剪后的时间支持及辅助比较均在逐例表中保留，没有选择更有利的定义替代主值。

判定：REFERENCE_INFORMATION_ACCOUNTED={valid}；INTERWINDOW_FREQUENCY_FEATURE={decision['INTERWINDOW_FREQUENCY_FEATURE']}；RELATIVE_MOTION_MAPPING={decision['RELATIVE_MOTION_MAPPING']}。能重复接收频率估计并不证明多模态峰身份稳定或动态物理映射正确；无真值支持的峰切换没有被手工修正。只有转向信号时不能宣布目标速度可精确估计；本轮未拟合目标状态或用oracle方位扣除平台。

冷复核{len(checks)}项，失败{sum(not r['PASS'] for r in checks)}，覆盖独立标量几何、差分/尺度、输入与旧文件字节及冻结顺序。外部导数/动态提供器认证不在其作用域。新传播、接收记录、源、噪声、状态优化全部0；DRIFT仍未评价，绝对径向零点未建立，H2回灌未授权。R4=0%。所有几何和两档SNR完整保存；提交推送后STOP，不自动修复提供器、换峰或开启新定位试验。
'''
 (OUT/'GPT_SYNC.md').write_text(text,encoding='utf-8');guard()
 dump(OUT/'OUTPUT_MANIFEST.json',{'time_utc':stamp(),'sha256':{p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'},'design_sha':read(OUT/'EXECUTION_STARTED.json')['design_sha'],'parent':PARENT,'total_result_bytes':sum(p.stat().st_size for p in OUT.iterdir() if p.is_file())})
 print(json.dumps(decision),flush=True)
if __name__=='__main__':main()
