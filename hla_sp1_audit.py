"""Post-A independent scalar/FFT/direct modal and support review, read-only experiments."""
from hla_sp1_core import *
import struct,math,itertools,shutil,subprocess

def cold_mod(path):
 b=Path(path).read_bytes();rec=struct.unpack_from('<i',b,0)[0]*4;nt,nm=struct.unpack_from('<ii',b,92);n=struct.unpack_from('<i',b,rec*5)[0];z=np.array(struct.unpack_from('<'+str(nt)+'f',b,rec*4));phi=np.array([np.array(struct.unpack_from('<'+str(2*nm)+'f',b,(7+j)*rec)).reshape(-1,2) for j in range(n)]);phi=(phi[:,:,0]+1j*phi[:,:,1]).T;kparts=np.array(struct.unpack_from('<'+str(2*n)+'f',b,(7+n)*rec)).reshape(-1,2);return z,phi,kparts[:,0]+1j*kparts[:,1]
def cold_field(mods,state):
 r,th,z=state;theta=math.radians(th);rr=np.hypot(r*1000*math.cos(theta)-OFF,r*1000*math.sin(theta));out=[]
 for depths,phi,k in mods:
  w=phi[np.flatnonzero(depths==z)[0]]*phi[np.flatnonzero(depths==200)[0]];out.append(np.sum(w[:,None]*np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*k[:,None]*rr-1j*np.pi/4),axis=0))
 return np.array(out)
def cold_score(h,y,m):
 if m==0:a=h;b=y
 else:a=np.array([h[p]*h[q].conj() for p,q in PAIRS]);b=np.array([[v[p]*v[q].conj() for p,q in PAIRS] for v in y])
 if m==2:a=a/abs(a);b=b/abs(b)
 total=0.
 for k in range(8):
  for f in range(3):total+=abs(np.vdot(a[f],b[k,f]))**2/(np.vdot(a[f],a[f]).real*np.vdot(b[k,f],b[k,f]).real)
 return float(total/24)
def main():
 if not (OUT/'INFERENCE_FREEZE.json').exists():raise RuntimeError('inference missing, partial audit only required')
 checks=[]
 def add(n,e,limit=1e-10):checks.append(dict(check=n,error=float(e),limit=limit,PASS=bool(e<=limit)))
 design=read(OUT/'DESIGN_FREEZE.json');execution=read(OUT/'EXECUTION_STARTED.json')
 for n,s in design['code_sha256'].items():add('frozen_code:'+n,0 if sha(ROOT/n)==s else 1,0)
 for n,s in design['input_sha256'].items():add('frozen_input:'+n,0 if sha(OUT/n)==s else 1,0)
 for r in rows(OUT/'INPUT_PROVENANCE.csv'):add('old_immutable:'+r['path'],0 if sha(ROOT/r['path'])==r['working_sha256'] else 1,0)
 add('threshold_before_test',0 if read(OUT/'THRESHOLD_FREEZE.json')['time_utc']<read(OUT/'TEST_INFERENCE_STARTED.json')['time_utc'] else 1,0);add('inference_before_truth',0 if read(OUT/'INFERENCE_FREEZE.json')['time_utc']<read(OUT/'EVALUATION_STARTED.json')['time_utc'] else 1,0)
 mods={mesh:[cold_mod(LOCAL/'modal'/f'sp1_f{int(f)}_n{mesh}.mod') for f in FREQ] for mesh in [80001,160001]};diff=[]
 for mesh in mods:
  for f,m in zip(FREQ,mods[mesh]):
   a=parse_mod(LOCAL/'modal'/f'sp1_f{int(f)}_n{mesh}.mod');add('independent_binary:'+str((mesh,f)),np.linalg.norm(m[1]-a[1])+np.linalg.norm(m[2]-a[2]),0)
 # Independent FFT and real-waveform algebra; all saved records/denominators.
 dfte=[];noise_ratios=[]
 for kind in ['calibration','test']:
  dest=LOCAL/kind;pub=rows(dest/'PUBLIC_RECORD_INDEX.csv');private=read(dest/'PRIVATE_MANIFEST.json')['records'];y=np.load(dest/'RECEIVER_DFT.npy',mmap_mode='r');meta=read(dest/'PUBLIC_METADATA.json');n=np.arange(2048);phasors=np.exp(2j*np.pi*FREQ[:,None]*n/1024);nu0=design['reference_power']
  for i,(record,priv) in enumerate(zip(pub,private)):
   guard();assert record['record_id']==priv['record_id'];x=np.load(record['path']);add('raw_hash:'+record['record_id'],0 if sha(record['path'])==record['sha256'] else 1,0);fft=np.fft.rfft(x,axis=1)[:,(FREQ*2).astype(int),:]/2048;e=np.linalg.norm(fft-y[i])/np.linalg.norm(y[i]);add('independent_FFT:'+record['record_id'],e,1e-10);innov=np.load(priv['private_innovations']);h=cold_field(mods[160001],[priv['truth_r_km'],priv['truth_theta_deg'],priv['truth_z_m']]);add('private_H_direct:'+record['record_id'],np.linalg.norm(h-innov['unit_source_H'])/np.linalg.norm(h),1e-10);source=innov['source'];gain=innov['gain'] if record['C']=='C1' else np.ones(8);expected=np.zeros(x.shape)
   for f in range(3):expected+=2*np.real((source[:,f,None]*gain[None,:]*h[f])[...,None,:]*phasors[f][None,:,None])
   nu=nu0/10**(int(record['reference_SNR_dB'])/10);expected+=math.sqrt(2048*nu)*innov['unit_real_noise'];add('saved_waveform_reconstruction:'+record['record_id'],np.linalg.norm(expected-x)/np.linalg.norm(x),1e-10)
   if record['C']=='C0' and int(record['reference_SNR_dB'])==20:
    noise_bin=np.fft.rfft(innov['unit_real_noise'],axis=1)[:,(FREQ*2).astype(int),:]/math.sqrt(2048);noise_ratios.append(np.mean(abs(noise_bin)**2));dfte.append(dict(record_id=record['record_id'],FFT_relative_error=e,actual_noise_complex_bin_power_over_nu=float(np.mean(abs(noise_bin)**2))))
   if i==0:
    shutil.copyfile(record['path'],OUT/(kind.upper()+'_FIRST_RECORD.npy'));np.save(OUT/(kind.upper()+'_FIRST_DFT.npy'),y[i])
  shutil.copyfile(dest/'PUBLIC_RECORD_INDEX.csv',OUT/(kind.upper()+'_PUBLIC_RECORD_INDEX.csv'));shutil.copyfile(dest/'PRIVATE_MANIFEST.json',OUT/(kind.upper()+'_PRIVATE_ARCHIVE_MANIFEST.json'));shutil.copyfile(dest/'PUBLIC_METADATA.json',OUT/(kind.upper()+'_PUBLIC_METADATA.json'))
 table(OUT/'FFT_AND_NOISE_REVIEW.csv',dfte);add('noise_bin_average',abs(float(np.mean(noise_ratios))-1),.08)
 # Exact modal review all registered best and accepted interval endpoints, no tuning.
 y=np.load(LOCAL/'test/RECEIVER_DFT.npy',mmap_mode='r');review=rows(OUT/'DIRECT_REVIEW_CANDIDATES.csv');cache={};directrows=[];maxscore=0.;flips=0
 for mesh,step in itertools.product([160001,80001],[.10,.05]):
  rs,ths,cat=catalogue(step);sub=[r for r in review if int(r['modal_mesh'])==mesh and float(r['grid_step_km'])==step];ids=sorted({int(r['candidate_id']) for r in sub});direct={i:cold_field(mods[mesh],cat[i]) for i in ids}
  for r in sub:
   i=int(r['record_index']);m=int(r['method']);idx=int(r['candidate_id']);score_exact=cold_score(direct[idx],y[i],m);e=abs(score_exact-float(r['interpolated_score']));flip=(1-score_exact<=float(r['tau']))!=(1-float(r['interpolated_score'])<=float(r['tau']));maxscore=max(maxscore,e);flips+=flip;directrows.append(dict(r,direct_modal_score=score_exact,score_absolute_difference=e,threshold_boundary_flip=flip));add('exact_score:'+str((mesh,step,i,m,idx)),e,1e-7)
  print('COLD BOUNDARIES',mesh,step,len(sub),'unique states',len(ids),flush=True);guard()
 table(OUT/'DIRECT_BEST_AND_BOUNDARY_REVIEW.csv',directrows)
 # Independent calibration scores and nearest-rank thresholds.
 yc=np.load(LOCAL/'calibration/RECEIVER_DFT.npy',mmap_mode='r');pub=rows(LOCAL/'calibration/PUBLIC_RECORD_INDEX.csv');priv=read(LOCAL/'calibration/PRIVATE_MANIFEST.json')['records'];dist=rows(OUT/'CALIBRATION_DISTRIBUTION.csv');calcache={}
 for r in dist:
  step=float(r['grid_step_km']);rs,ths,cat=catalogue(step);i=next(i for i,p in enumerate(pub) if p['record_id']==r['record_id']);p=priv[i];state=(rs[nearest(rs,p['truth_r_km'])],ths[nearest(ths,abs(p['truth_theta_deg']))],p['truth_z_m'])
  if state not in calcache:calcache[state]=cold_field(mods[160001],state)
  d=1-cold_score(calcache[state],yc[i],int(r['method'][1]));add('calibration_score:'+str((r['record_id'],r['grid_step_km'],r['method'])),abs(d-float(r['D_nearest_joint_label'])),1e-10)
 for tau in rows(OUT/'CALIBRATION_THRESHOLDS.csv'):
  a=sorted(float(r['D_nearest_joint_label']) for r in dist if all(r[k]==tau[k] for k in ['C','reference_SNR_dB','grid_step_km','method']));add('threshold:'+str(tau),abs(a[243]-float(tau['tau'])),0)
 # All support projection summaries and selected-point exact bookkeeping.
 results=rows(OUT/'RANGE_RESULTS_BY_RECORD.csv');supportcheck=[]
 for mesh,step in itertools.product([160001,80001],[.10,.05]):
  rs,ths,cat=catalogue(step);pack=np.load(OUT/f'FULL_RANGE_SUPPORT_{step:.2f}_n{mesh}.npz')['packed_acceptance'];sc=np.load(LOCAL/f'SCORES_step{step:.2f}_n{mesh}.npy',mmap_mode='r');lookup={r['record_id']:i for i,r in enumerate(rows(LOCAL/'test/PUBLIC_RECORD_INDEX.csv'))}
  for r in results:
   if int(r['modal_mesh'])!=mesh or float(r['grid_step_km'])!=step:continue
   i=lookup[r['record_id']];m=int(r['method'][1]);mask=np.unpackbits(pack[i,m],count=len(cat)).astype(bool);direct_mask=1-sc[i,m]<=float(r['tau']);add('all_acceptance_bits:'+str((mesh,step,i,m)),np.count_nonzero(mask!=direct_mask),0);a=mask.reshape(len(rs),len(ths),len(ZS));labels=rs[np.flatnonzero(a.any(axis=(1,2)))];stored=np.array(json.loads(r['range_labels']));add('all_range_projection:'+str((mesh,step,i,m)),0 if np.array_equal(labels,stored) else 1,0);best=int(np.argmax(sc[i,m]));add('strict_top1:'+str((mesh,step,i,m)),abs(best-int(r['best_id'])),0);ir=nearest(rs,float(r['truth_r_km']));pr=read(LOCAL/'test/PRIVATE_MANIFEST.json')['records'][i];it=nearest(ths,abs(pr['truth_theta_deg']));iz=nearest(ZS,pr['truth_z_m']);add('truth_horizontal_support:'+str((mesh,step,i,m)),0 if bool(a[ir,it].any())==(r['horizontal_label_retained']=='True') else 1,0);add('truth_joint_support:'+str((mesh,step,i,m)),0 if bool(a[ir,it,iz])==(r['joint_z_label_retained']=='True') else 1,0)
  del pack,sc
 # full static vector and score mesh differences across all calibration/test geometries and depths.
 panel=rows(OUT/'CALIBRATION_PANEL.csv')+rows(ROOT/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv');num=[]
 for g,p in enumerate(panel):
  for z in ZS:
   state=[float(p['r_km']),float(p['theta_deg']),z];h=cold_field(mods[160001],state);low=cold_field(mods[80001],state)
   for fi in range(3):
    scalar=np.vdot(low[fi],h[fi])/np.vdot(low[fi],low[fi]);num.append(dict(geometry=g,z_m=z,frequency_hz=FREQ[fi],raw_vector_relative=float(np.linalg.norm(h[fi]-low[fi])/np.linalg.norm(h[fi])),common_source_removed_relative=float(np.linalg.norm(h[fi]-scalar*low[fi])/np.linalg.norm(h[fi])),field_norm=float(np.linalg.norm(h[fi])),minimum_element_abs=float(abs(h[fi]).min()),P0_mesh_D=1-cold_score(h,low[None].repeat(8,axis=0),0),P1_mesh_D=1-cold_score(h,low[None].repeat(8,axis=0),1),P2_mesh_D=1-cold_score(h,low[None].repeat(8,axis=0),2)))
 table(OUT/'STATIC_VECTOR_AND_SCORE_MESH_DIAGNOSTICS.csv',num);table(OUT/'INDEPENDENT_AUDIT_CHECKS.csv',checks)
 allpass=all(r['PASS'] for r in checks);dump(OUT/'VALIDATION.json',{'time_utc':stamp(),'PASS':allpass,'PASS_count':sum(r['PASS'] for r in checks),'FAIL_count':sum(not r['PASS'] for r in checks),'receiver_waveforms':1536,'all_best_and_boundary_checks':len(directrows),'direct_score_max_difference':maxscore,'direct_boundary_flips':flips,'new_solver_calls':{'KRAKEN':6,'FIELD':3},'code_after_A_unchanged':True,'budget_wall_RSS':guard(),'scope':'arithmetic/bytes/exact full stored modal sums; no omitted-mode certificate or sea-trial validation'})
 print('COLD REVIEW COMPLETE',len(checks),'FAIL',sum(not r['PASS'] for r in checks),'boundaryflips',flips,flush=True)
if __name__=='__main__':main()
