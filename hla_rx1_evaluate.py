"""Truth-only evaluation, independently recomputes kinematic quantities.
No calibration, source correction, sign choice or peak selection is returned to RX.
"""
import csv,json
import numpy as np
import pandas as pd
from hla_rx1_common import *

def radial(row,t):
 p,v=state(row);t=np.asarray(t);d=p+t[:,None]*v-platform(t)
 vp=np.where((t<600)[:,None],np.array([2.,0.]),2*np.array([np.cos(np.pi/12),np.sin(np.pi/12)]))
 return np.linalg.norm(d,axis=1),np.sum(d*(v-vp),axis=1)/np.linalg.norm(d,axis=1)

def quantile(x,p):
 x=sorted(x);return x[max(0,int(np.ceil(len(x)*p))-1)] if x else float('inf')

def main():
 cfg=read(OUT/'DESIGN_FREEZE.json');status=read(OUT/'EXECUTION_STATUS.json');truth=list(csv.DictReader((ROOT/cfg['truth_file']).open(encoding='utf-8-sig')))
 if not (OUT/'EXTRACTED_FEATURES.csv').exists():
  dump(OUT/'DECISION.json',{'classifications':['IMPLEMENTATION_OR_PROVIDER_LIMIT','PARTIAL_EXECUTION'],'reason':status['partial'],'H2_mapping':'NOT_EVALUATED','R4_percent':0});return
 frame=pd.read_csv(OUT/'EXTRACTED_FEATURES.csv',float_precision='round_trip',keep_default_na=False);maps=[]
 for d in frame.to_dict('records'):
  row=truth[int(d['geometry'])];j=int(d['window']);t=np.linspace(j*200,(j+1)*200,801);r,q=radial(row,t);qbar=float((r[-1]-r[0])/200);mid=float(radial(row,np.array([j*200+100.]))[1][0]);meanabs=float(np.trapz(abs(q),t)/200)
  # Bias of a finite-time correlation peak cannot be decomposed as white Gaussian.
  error=abs(float(d['q_eff_mps'])-abs(qbar)) if d['q_eff_mps']!='' else float('inf')
  maps.append(dict(**d,qbar_mps=qbar,midpoint_signed_q_mps=mid,window_mean_abs_q_mps=meanabs,abs_qbar_vs_mean_abs_difference=abs(qbar)-meanabs,
   qeff_vs_abs_qbar_error_mps=error,qeff_vs_abs_midpoint_error_mps=abs(float(d['q_eff_mps'])-abs(mid)),
   uncond_error_mps=error if d['accepted'] else float('inf'),signed_bias='NOT_DEFINED_UNRESOLVED_SIGN',
   correlation_weighted_radial_truth='NOT_CERTIFIED_NONLINEAR_LAG_PEAK',mapping_status='RX_FEATURE_MAPPING_NOT_ESTABLISHED',truth_sign_branch_retained=True,wrong_forced_sign=False))
 m=pd.DataFrame(maps);m.to_csv(OUT/'WINDOW_FEATURE_MAPPING.csv',index=False)
 stats=[]
 for (family,snr,method,lag),g in m.groupby(['source_family','SNR_dB','method','lag_s'],sort=False):
  accepted=g[g.accepted];values=accepted.qeff_vs_abs_qbar_error_mps.tolist();uv=g.uncond_error_mps.tolist()
  stats.append(dict(source_family=family,SNR_dB=int(snr),method=method,lag_s=int(lag),registered_records=64,completed_records=g.record_id.nunique(),registered_line_windows=1152,evaluated_line_windows=len(g),accepted_line_windows=len(accepted),
   conditional_median_abs_error_mps=quantile(values,.5),conditional_P90_abs_error_mps=quantile(values,.9),conditional_P95_abs_error_mps=quantile(values,.95),conditional_magnitude_RMSE=float(np.sqrt(np.mean(np.square(values)))) if values else float('inf'),
   reject_INF_unconditional_median_mps=quantile(uv,.5),reject_INF_unconditional_P95_mps=quantile(uv,.95),
   **{f'fraction_registered_with_error_le_{v}_mps':sum(a<=v for a in uv)/1152 for v in [.02,.05,.10,.20]},
   mapping='not a Gaussian sigma; qbar comparisons diagnostic only'))
 table(OUT/'FEATURE_ERROR_AND_COVERAGE.csv',stats)
 # Fixed equal median fusion of three selected peaks; retains +/- and all raw peaks.
 fusion=[]
 for key,g in m.groupby(['record_id','geometry','source_family','SNR_dB','replicate','window','method','lag_s'],sort=False):
  ok=bool(g.accepted.all());q=float(np.median(g.q_eff_mps));err=abs(q-abs(g.qbar_mps.iloc[0]))
  fusion.append(dict(zip(['record_id','geometry','source_family','SNR_dB','replicate','window','method','lag_s'],key))|dict(q_eff_mps=q,accepted=ok,error_vs_abs_qbar_mps=err if ok else float('inf'),direction='UNRESOLVED_PLUS_MINUS',fusion='equal median of three lines; no truth peak choice',mapping='NOT_ESTABLISHED'))
 table(OUT/'FUSED_FEATURES.csv',fusion)
 pairs=m.pivot(index=['record_id','window','line','lag_s'],columns='method',values='q_eff_mps').reset_index();pairs['REF_minus_RX_mps']=pairs.REF-pairs.RX;pairs.to_csv(OUT/'SOURCE_REFERENCE_COMPARISON.csv',index=False)
 dep=[]
 for family in ['STABLE','DRIFT']:
  for snr in [20,5]:
   for method in ['RX','REF']:
    g=m[(m.source_family==family)&(m.SNR_dB==snr)&(m.method==method)&(m.lag_s==200)]
    if len(g):
     wide=g.pivot(index='record_id',columns=['window','line'],values='qeff_vs_abs_qbar_error_mps');corr=wide.corr();dep.append(dict(source_family=family,SNR_dB=snr,method=method,record_count=len(wide),max_abs_cross_frequency_or_window_correlation=float(np.nanmax(abs(corr.to_numpy()-np.eye(18)))),bearing_error_dependence='NOT_EVALUATED_NO_BEARING_OBSERVATIONS_OR_BEARING_ESTIMATOR_IN_RX1'))
 table(OUT/'EMPIRICAL_DEPENDENCE.csv',dep)
 conditions=[]
 for s in stats:
  if s['lag_s']==200:
   conditions.append(f"| {s['source_family']}/{s['SNR_dB']} dB/{s['method']} | {s['accepted_line_windows']}/{s['registered_line_windows']} | {s['conditional_median_abs_error_mps']:.4f} | {s['conditional_P95_abs_error_mps']:.4f} | {s['fraction_registered_with_error_le_0.05_mps']:.1%} |")
 (OUT/'MEASUREMENT_CONDITION_TABLE.md').write_text('# ???????????\n\n| ?? | ??/???? | ?abs(qbar)???? m/s | P95 m/s | ?????0.05?? |\n|---|---:|---:|---:|---:|\n'+'\n'.join(conditions)+'\n\n??????H2???????????????????Gaussian sigma?RX?????LO????????????c_ref????REF??????????????????????????????????????????????????????????????????????????????????????\n',encoding='utf-8')
 dump(OUT/'DECISION.json',{'classifications':['FEATURE_MAPPING_NOT_ESTABLISHED']+(['PARTIAL_EXECUTION'] if status['partial'] else []),'conditions':'per-source/SNR/method full denominators in FEATURE_ERROR_AND_COVERAGE.csv','receiver_extraction':'actual multichannel sampled records processed, not q+Gaussian','RX_mean_frequency_removal':'receiver peak centering intentionally removes unknown constant source offset and mean Doppler together; absolute speed zero not identifiable by this preprocessing','REF_permission':'source phase at receiver time; not retarded oracle alignment','H2_feedback':'NOT_AUTHORIZED_AND_MAPPING_NOT_ESTABLISHED','no_general_FDSL_NoGo':True,'R4_percent':0})
if __name__=='__main__':main()
