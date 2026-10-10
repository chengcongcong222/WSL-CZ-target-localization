"""Full finite-catalogue inference. No test truth, source or gain access."""
from hla_sp1_core import *
def main():
 if (OUT/'TEST_INFERENCE_STARTED.json').exists():raise RuntimeError('once only inference')
 freeze=read(OUT/'THRESHOLD_FREEZE.json');assert sha(OUT/'CALIBRATION_THRESHOLDS.csv')==freeze['threshold_sha256'];dump(OUT/'TEST_INFERENCE_STARTED.json',{'time_utc':stamp(),'threshold_freeze_sha256':sha(OUT/'THRESHOLD_FREEZE.json'),'private_truth_reads':0})
 m=read(OUT/'TEST_DFT_FREEZE.json');assert sha(m['file'])==m['sha256'];y=np.load(m['file'],mmap_mode='r');public=rows(LOCAL/'test/PUBLIC_RECORD_INDEX.csv');taus=rows(OUT/'CALIBRATION_THRESHOLDS.csv');moments=[np.array([realmoment(moment(a,method)).ravel()/3 for a in y]) for method in range(3)];results=[];boundaries=[];catalogues=[]
 for step in [.10,.05]:
  rs,ths,cat=catalogue(step);np.savez(OUT/f'CATALOGUE_{step:.2f}.npz',ranges_km=rs,nonnegative_theta_deg=ths,depths_m=ZS,theta_signs=np.array([-1,1]),candidate_id_order=np.array(['r','abs_theta','z']));catalogues.append(dict(step_km=step,nonnegative_representatives=len(cat),mirror_rule='theta=0 once, positive theta restored +/-; all z labels preserved'))
  for mesh in [160001,80001]:
   provider=Provider(mesh);path=LOCAL/f'SCORES_step{step:.2f}_n{mesh}.npy';scores=np.lib.format.open_memmap(path,mode='w+',dtype=np.float64,shape=(len(public),3,len(cat)));packed=np.empty((len(public),3,(len(cat)+7)//8),np.uint8);profiles=np.empty((len(public),3,len(rs)),float)
   for a in range(0,len(cat),4096):
    guard();h=provider.field(cat[a:a+4096])
    for method in range(3):
     feat=realfeatures(transform(h,method)).reshape(len(h),-1);b=moments[method]@feat.T
     if not np.isfinite(b).all() or b.min()<-1e-10 or b.max()>1+1e-10:raise ValueError('score bounds/formula invalid')
     scores[:,method,a:a+len(h)]=b
    if a%65536==0:print('CATALOGUE',step,mesh,a,'/',len(cat),'budget',guard(),flush=True)
   scores.flush()
   for i,r in enumerate(public):
    for method in range(3):
     tau=float(next(t['tau'] for t in taus if t['C']==r['C'] and t['reference_SNR_dB']==r['reference_SNR_dB'] and float(t['grid_step_km'])==step and t['method']=='P'+str(method)));b=np.asarray(scores[i,method]);accept=1-b<=tau;packed[i,method]=np.packbits(accept);br=b.reshape(len(rs),len(ths),len(ZS));profile=br.max(axis=(1,2));profiles[i,method]=profile;ir=np.flatnonzero(accept.reshape(br.shape).any(axis=(1,2)));ints=intervals(ir,rs,step);best=int(np.argmax(b));valid=len(ir)>0;width=ints[-1][1]-ints[0][0] if valid else 'EMPTY_NOT_ZERO';length=sum(hi-lo for lo,hi in ints);states=cat[np.flatnonzero(accept)];theta_max=float(states[:,1].max()) if len(states) else None;zvalues=sorted(set(states[:,2])) if len(states) else []
     results.append(dict(record_id=r['record_id'],C=r['C'],reference_SNR_dB=r['reference_SNR_dB'],method='P'+str(method),grid_step_km=step,modal_mesh=mesh,numerical_valid=True,effective_output=valid,tau=tau,best_id=best,best_score=float(b[best]),best_r_km=float(cat[best,0]),best_abs_theta_deg=float(cat[best,1]),best_z_m=float(cat[best,2]),best_r_at_boundary=bool(cat[best,0] in [45,60]),best_z_at_boundary=bool(cat[best,2] in [150,250]),accepted_nonnegative_states=int(accept.sum()),accepted_mirrored_states=int(2*accept.sum()-accept.reshape(br.shape)[:,0,:].sum()),range_labels=json.dumps(rs[ir].tolist()),range_intervals_km=json.dumps(ints),interval_count=len(ints),envelope_width_km=width,occupied_cell_length_km=length,theta_mirror_min_deg=-theta_max if theta_max is not None else None,theta_mirror_max_deg=theta_max,z_values_m=json.dumps(zvalues),full_support_file=f'FULL_RANGE_SUPPORT_{step:.2f}_n{mesh}.npz'))
     ids={best}
     for lo,hi in ints:
      for idx in [nearest(rs,lo+step/2),nearest(rs,hi-step/2)]:
       k=int(np.argmax(br[idx]));ids.add(idx*len(ths)*len(ZS)+k)
     for idx in sorted(ids):boundaries.append(dict(record_id=r['record_id'],record_index=i,method=method,grid_step_km=step,modal_mesh=mesh,candidate_id=idx,interpolated_score=float(b[idx]),tau=tau,role='BEST_OR_PROJECTED_ACCEPTED_INTERVAL_BOUNDARY'))
   np.savez_compressed(OUT/f'FULL_RANGE_SUPPORT_{step:.2f}_n{mesh}.npz',packed_acceptance=packed,record_ids=np.array([r['record_id'] for r in public]),candidate_count=len(cat),methods=np.array(['P0','P1','P2']),bitorder='big',mirrors='restore theta +/- except0; all depth labels included');np.save(LOCAL/f'RANGE_PROFILES_{step:.2f}_n{mesh}.npy',profiles);del scores,provider,packed
 table(OUT/'INFERENCE_RESULTS_NO_TRUTH.csv',results);table(OUT/'DIRECT_REVIEW_CANDIDATES.csv',boundaries);dump(OUT/'INFERENCE_FREEZE.json',{'time_utc':stamp(),'private_truth_source_gain_reads':0,'DFT_sha256':m['sha256'],'results_sha256':sha(OUT/'INFERENCE_RESULTS_NO_TRUTH.csv'),'catalogues':catalogues,'completed_configurations':len(results),'score_files':{p.name:sha(p) for p in LOCAL.glob('SCORES*.npy')},'full_support_files':{p.name:sha(p) for p in OUT.glob('FULL_RANGE_SUPPORT*.npz')}})
if __name__=='__main__':main()
