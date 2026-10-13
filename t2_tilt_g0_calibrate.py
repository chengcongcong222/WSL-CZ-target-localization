"""Independent calibration, no test records/optimizer/full-catalogue inference."""
from t2_tilt_g0_core import *

def main():
 if (OUT/'THRESHOLD_FREEZE.json').exists():raise RuntimeError('thresholds already sealed')
 start=time.monotonic();public=rows(OUT/'CALIBRATION_DFT_INDEX.csv');private=read(LOCAL/'calibration/private/MANIFEST.json');truth={r['base_id']:r['truth'] for r in private['private_records']};result=[];all_d=[]
 for fine in [False,True]:
  name='fine' if fine else 'primary';s=catalogue(fine);axes=[np.unique(s[:,j]) for j in range(3)]
  for mesh in [80001,160001]:
   for res in RES:
    guard(start);h=np.load(LOCAL/'replicas'/f'{res}_{name}_n{mesh}.npy')
    for c in ['C0','C1']:
     for snr in ['20','5']:
      pp=[r for r in public if r['resource']==res and r['calibration']==c and r['reference_SNR_dB']==snr];ds=[]
      for r in pp:
       t=truth[int(r['base_id'])];near=[axes[j][np.argmin(abs(axes[j]-v))] for j,v in enumerate([t['r_km'],abs(t['theta_deg']),t['z_m']])];ci=int(np.flatnonzero(np.all(s==near,axis=1))[0]);y=np.load(r['DFT_path']);d=float(1-score(h[ci:ci+1],y)[0]);ds.append(d);all_d.append({'resource':res,'catalogue':name,'mesh':mesh,'calibration':c,'reference_SNR_dB':snr,'record_id':int(r['record_id']),'base_id':int(r['base_id']),'nearest_candidate_id':ci,'D':d})
      assert len(ds)==256;result.append({'resource':res,'catalogue':name,'mesh':mesh,'calibration':c,'reference_SNR_dB':snr,'n':256,'rank_one_based':244,'tau_D':float(np.sort(ds)[243])})
 table(OUT/'CALIBRATION_THRESHOLDS.csv',result);table(OUT/'CALIBRATION_NEAREST_SCORES.csv',all_d)
 dump(OUT/'THRESHOLD_FREEZE.json',{'utc':utc(),'thresholds_SHA256':sha(OUT/'CALIBRATION_THRESHOLDS.csv'),'scores_SHA256':sha(OUT/'CALIBRATION_NEAREST_SCORES.csv'),'thresholds':len(result),'test_generated_yet':False,'empirical_filter_not_theoretical95pct_domain':True,'science_wall_s':time.monotonic()-start});print('THRESHOLDS SEALED',len(result),flush=True)
if __name__=='__main__':main()