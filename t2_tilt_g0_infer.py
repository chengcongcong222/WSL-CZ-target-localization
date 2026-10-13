"""Complete catalogue inference, reads no test truth or generator values."""
from t2_tilt_g0_core import *

def main():
 if (OUT/'INFERENCE_STARTED.json').exists():raise RuntimeError('inference no rerun')
 assert read(OUT/'THRESHOLD_FREEZE.json')['thresholds_SHA256']==sha(OUT/'CALIBRATION_THRESHOLDS.csv')
 dump(OUT/'INFERENCE_STARTED.json',{'utc':utc()});start=time.monotonic();pub=rows(OUT/'TEST_DFT_INDEX.csv');taus=rows(OUT/'CALIBRATION_THRESHOLDS.csv');inferred=[];archives=[];folder=LOCAL/'inference';folder.mkdir()
 for fine in [False,True]:
  cat='fine' if fine else 'primary';s=catalogue(fine);rs=np.unique(s[:,0]);nt=len(np.unique(s[:,1]));step=rs[1]-rs[0]
  for mesh in [80001,160001]:
   for res in RES:
    guard(start);pp=[r for r in pub if r['resource']==res];y=np.stack([np.load(r['DFT_path']) for r in pp]);h=np.load(LOCAL/'replicas'/f'{res}_{cat}_n{mesh}.npy');d=1-scores_many(h,y);sp=folder/f'SCORES_{res}_{cat}_n{mesh}.npy';np.save(sp,d);bitmap=[]
    for j,r in enumerate(pp):
     t=next(t for t in taus if t['resource']==res and t['catalogue']==cat and int(t['mesh'])==mesh and t['calibration']==r['calibration'] and t['reference_SNR_dB']==r['reference_SNR_dB']);tau=float(t['tau_D']);accepted=d[j]<=tau;bitmap.append(np.packbits(accepted,bitorder='little'));ci=int(np.argmin(d[j]));kept=s[accepted];ids=np.flatnonzero(np.any(accepted.reshape(len(rs),nt,5),axis=(1,2)));ints=[]
     if len(ids):
      for chunk in np.split(ids,np.where(np.diff(ids)>1)[0]+1):ints.append([float(max(45,rs[chunk[0]]-step/2)),float(min(60,rs[chunk[-1]]+step/2))])
     width=float(ints[-1][1]-ints[0][0]) if ints else 'INF';occup=float(sum(b-a for a,b in ints)) if ints else 'INF';rank=np.argsort(d[j],kind='stable')[:2]
     inferred.append({'record_id':int(r['record_id']),'base_id':int(r['base_id']),'resource':res,'calibration':r['calibration'],'reference_SNR_dB':int(r['reference_SNR_dB']),'catalogue':cat,'mesh':mesh,'accepted_output':bool(len(kept)),'accepted_nonnegative_states':int(accepted.sum()),'accepted_mirrored_states':int(2*accepted.sum()-np.sum(accepted&(s[:,1]==0))),'range_labels':json.dumps(rs[ids].tolist()),'range_cell_intervals':json.dumps(ints),'envelope_km':width,'occupied_km':occup,'disconnected_intervals':len(ints),'best_candidate_id':ci,'best_r_km':float(s[ci,0]),'best_abs_theta_deg':float(s[ci,1]),'best_z_m':float(s[ci,2]),'best_D':float(d[j,ci]),'second_D':float(d[j,rank[1]]),'tau_D':tau,'best_hits_boundary':bool(s[ci,0] in [45,60] or s[ci,1] in [0,5] or s[ci,2] in [180,220]),'theta_mirrors_retained':True,'support_is_finite_not_continuous':True})
    bp=folder/f'ACCEPTED_{res}_{cat}_n{mesh}.npy';np.save(bp,np.array(bitmap,dtype=np.uint8));ids_path=folder/f'RECORD_IDS_{res}_{cat}_n{mesh}.npy';np.save(ids_path,np.array([int(r['record_id']) for r in pp]));archives.extend([sp,bp,ids_path]);print('INFER',res,cat,mesh,len(pp),'configs',flush=True)
 table(OUT/'INFERENCE_RESULTS_NO_TRUTH.csv',inferred)
 dump(OUT/'INFERENCE_FREEZE.json',{'utc':utc(),'result_SHA256':sha(OUT/'INFERENCE_RESULTS_NO_TRUTH.csv'),'archives':{str(p):sha(p) for p in archives},'resource_configurations':len(pub),'paired_numeric_outputs':len(inferred),'truth_not_accessed':True,'science_wall_s':time.monotonic()-start})
 # Compact complete support exported, all catalogues and all mirrors implied by frozen exact symmetry.
 (OUT/'complete_support').mkdir()
 for p in archives:
  if p.name.startswith(('ACCEPTED','RECORD_IDS')):shutil_copy(p,OUT/'complete_support'/p.name)
 print('INFERENCE SEALED',len(inferred),flush=True)
def shutil_copy(src,dst):
 import shutil
 shutil.copyfile(src,dst)
if __name__=='__main__':main()