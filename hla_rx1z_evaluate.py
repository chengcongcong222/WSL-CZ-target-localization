import math,itertools
from hla_rx1z_common import *
OLD=ROOT/'results/R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT'
TRUTH=ROOT/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv'
def rho(h,t):
 a=math.radians(float(h['theta_deg']));p=math.radians(float(h['psi_deg']));r=float(h['r_km'])*1000;v=float(h['v_mps']);turn=math.radians(15)
 px=2*min(t,600)+2*max(t-600,0)*math.cos(turn);py=2*max(t-600,0)*math.sin(turn)
 return math.hypot(r*math.cos(a)+v*t*math.cos(p)-px,r*math.sin(a)+v*t*math.sin(p)-py)
def stats(a):
 e=np.array([float(r['error_mps']) if r['valid']=='True' else np.inf for r in a]);ok=e[np.isfinite(e)];q=np.array([float(r['true_delta_qbar_mps']) for r in a]);aux=np.array([float(r['aux_error_mps']) if r['valid']=='True' else np.inf for r in a]);zero=float(np.sqrt(np.mean(q*q)));rmse=float(np.sqrt(np.mean(e*e)))
 return dict(planned=len(a),valid=len(ok),missing=len(a)-len(ok),bias_mps=float(np.mean(ok)) if len(ok) else 'NOT_AVAILABLE',MAE_mps=float(np.mean(abs(e))),RMSE_mps=rmse,median_abs_error_mps=float(np.median(abs(e))),maximum_abs_error_mps=float(np.max(abs(e))),conditional_RMSE_mps=float(np.sqrt(np.mean(ok*ok))) if len(ok) else 'NOT_AVAILABLE',true_change_RMS_mps=zero,true_change_median_abs_mps=float(np.median(abs(q))),true_change_max_abs_mps=float(np.max(abs(q))),zero_predictor_RMSE_mps=zero,zero_predictor_median_abs_error_mps=float(np.median(abs(q))),RMSE_to_zero_ratio=rmse/zero if zero else 'UNDEFINED',aux_support_RMSE_mps=float(np.sqrt(np.mean(aux*aux))),beats_zero=rmse<zero)
def main():
 m=read(OUT/'FEATURE_OUTPUT_MANIFEST.json')
 for n,s in m['feature_sha256'].items():assert sha(OUT/n)==s
 dump(OUT/'EVALUATION_STARTED.json',{'time_utc':stamp(),'receiver_manifest_sha256':sha(OUT/'FEATURE_OUTPUT_MANIFEST.json'),'truth_opened_after_receiver_freeze':True})
 guard();truth=rows(TRUTH);records={r['record_id']:r for r in read(OLD/'RECEIVER_DATA_MANIFEST.json')['records']};a=[]
 for f in rows(OUT/'INTERWINDOW_FEATURES.csv'):
  rec=records[f['record_id']];g=int(rec['geometry']);h=truth[g];j=int(f['difference_j']);q=[(rho(h,200*(k+1))-rho(h,200*k))/200 for k in range(6)];dq=q[j+1]-q[j];low=(rho(h,float(f['from_end_s']))-rho(h,float(f['from_start_s'])))/(float(f['from_end_s'])-float(f['from_start_s']));high=(rho(h,float(f['to_end_s']))-rho(h,float(f['to_start_s'])))/(float(f['to_end_s'])-float(f['to_start_s']));aux=high-low;valid=f['valid']=='True' and bool(f['apparent_radial_change_mps']);v=float(f['apparent_radial_change_mps']) if valid else float('inf')
  a.append(dict(f,geometry=g,SNR_dB=int(rec['SNR_dB']),rep=int(rec['replicate']),valid=valid,true_delta_qbar_mps=dq,aux_support_delta_qbar_mps=aux,error_mps=v-dq,aux_error_mps=v-aux,diagnostic='MAPPING_DIAGNOSTIC_UNDER_RX1_FALLBACK_MODEL'))
 table(OUT/'MAPPING_DIAGNOSTIC_ALL.csv',a);summary=[]
 for g,s,l,scope in itertools.product(range(8),[20,5],range(3),['ALL','TURN','SAME_LEG']):
  sub=[r for r in a if r['geometry']==g and r['SNR_dB']==s and int(r['line'])==l and (scope=='ALL' or r['group']==scope)]
  # normalize in-memory boolean for the same reader-facing statistics
  for r in sub:r['valid']=str(r['valid'])
  summary.append(dict(geometry=g,SNR_dB=s,line=l,scope=scope,**stats(sub)))
 table(OUT/'RELATIVE_MAPPING_BY_GEOMETRY.csv',summary);table(OUT/'TURN_VS_SAME_LEG.csv',[r for r in summary if r['scope']!='ALL']);pooled=[]
 for s,scope in itertools.product([20,5],['ALL','TURN','SAME_LEG']):pooled.append(dict(SNR_dB=s,scope=scope,**stats([r for r in a if r['SNR_dB']==s and (scope=='ALL' or r['group']==scope)])))
 table(OUT/'POOLED_DESCRIPTIVE_ONLY.csv',pooled);guard()
 print('EVALUATED',len(a),'contrasts',flush=True)
if __name__=='__main__':main()
