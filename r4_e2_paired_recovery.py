"""Once-only registered paired-frequency recovery. No FIM, RNG, audio, optimizer or new environment."""
from pathlib import Path
import os,time,subprocess,json
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import r4_e2_g0 as old
import r4_e2_forward_diagnostic as base
O=Path('results/R4_E2_PAIRED_FREQUENCY_RECOVERY')
PARENT='e6f1d68f2aec7d905652e5159593b1985f2e2517'
PAIRS=[(151.,2,150.,152.),(152.5,5,150.,155.),(155.,10,150.,160.),(199.,2,198.,200.),(202.5,5,200.,205.),(200.,10,195.,205.),(245.,2,244.,246.),(246.5,5,244.,249.),(245.,10,240.,250.)]
RAW=sorted({f for p in PAIRS for f in p[2:]})
GRIDS=[20001,40001,80001,160001]
TRANS=list(zip(GRIDS[:-1],GRIDS[1:]))
read,dump,write,git,sha=base.read,base.dump,base.write,base.git,base.sha
def path_for(f,g):
 if g<80001:return old.MOD/f'f{int(f*2):04d}_s+0_n{g}.mod'
 if g==80001 and f in [150,200,244,250]:return base.OUT/'third_mesh'/f'fine_f{int(f)}_n80001.mod'
 return O/'modal'/f'f{int(f*2):04d}_n{g}.mod'
def generate(start):
 rows=[];available={};stop=False
 for g in GRIDS:
  for f in RAW:
   path=path_for(f,g);existed=(g<80001 or g==80001 and f in [150,200,244,250])
   if stop and not existed:
    rows.append(dict(mesh=g,frequency_hz=f,new_call=False,elapsed_s=0,status='NOT_ATTEMPTED_AFTER_PROVIDER_FAILURE',modal_file=str(path),timeout_s=90));continue
   if existed:
    m=old.parse_mod(path);available[g,f]=m
    rows.append(dict(mesh=g,frequency_hz=f,new_call=False,elapsed_s=0,status='REUSED',modal_file=str(path),timeout_s=90));continue
   stem=path.stem;path.parent.mkdir(exist_ok=True)
   path.with_suffix('.env').write_text(old.modal_env(f,0,g),encoding='utf-8')
   begin=time.monotonic();result=None
   try:
    remain=1980-(begin-start)
    if remain<=0:raise TimeoutError('frozen solver budget exhausted')
    result=subprocess.run([str(old.AT/'kraken.exe'),stem],cwd=path.parent,capture_output=True,timeout=min(90,remain))
    path.with_suffix('.stdout.txt').write_bytes(result.stdout+result.stderr)
    if result.returncode!=0:raise RuntimeError('KRAKEN_EXIT_'+str(result.returncode))
    m=old.parse_mod(path);available[g,f]=m;status='GENERATED'
   except Exception as e:
    if isinstance(e,subprocess.TimeoutExpired):path.with_suffix('.stdout.txt').write_bytes((e.stdout or b'')+(e.stderr or b''))
    status=type(e).__name__+':'+str(e);stop=True
   rows.append(dict(mesh=g,frequency_hz=f,new_call=True,elapsed_s=time.monotonic()-begin,status=status,modal_file=str(path),timeout_s=90))
   write(O/'PAIRED_MODAL_GENERATION.csv',rows)
   print('MODAL',g,f,status,flush=True)
 write(O/'PAIRED_MODAL_GENERATION.csv',rows)
 return rows,available
def mode_diagnostics(mods):
 out=[]
 for f in RAW:
  for g,h in TRANS:
   row=dict(frequency_hz=f,coarse_mesh=g,fine_mesh=h,status='MISSING',mode_count_coarse='',mode_count_fine='',same_count='',identity_unambiguous='',min_sampled_phi_correlation='',max_phase_delta_60km_rad='',max_k_imag='',max_signed_phi_relative='',normalization='FULL_DEPTH_INTEGRAL_NOT_AVAILABLE_13_OUTPUT_DEPTHS')
   if (g,f) in mods and (h,f) in mods:
    z,p,k=mods[g,f];zz,q,kk=mods[h,f];same=len(k)==len(kk)
    row.update(status='AVAILABLE',mode_count_coarse=len(k),mode_count_fine=len(kk),same_count=same,max_k_imag=float(max(k.imag.max(),kk.imag.max())))
    if same:
     gap=np.minimum(np.r_[np.inf,abs(np.diff(k.real))],np.r_[abs(np.diff(k.real)),np.inf])
     corr=abs(np.sum(p.conj()*q,axis=0))/np.maximum(np.linalg.norm(p,axis=0)*np.linalg.norm(q,axis=0),1e-300)
     sign=np.where(np.sum(p.real*q.real,axis=0)>=0,1,-1)
     row.update(identity_unambiguous=bool(np.all(abs(kk.real-k.real)<.45*gap)),min_sampled_phi_correlation=float(corr.min()),max_phase_delta_60km_rad=float(np.max(abs(kk.real-k.real))*60000),max_signed_phi_relative=float(np.linalg.norm(p-q*sign)/np.linalg.norm(p)))
    else:row['identity_unambiguous']=False
   out.append(row)
 write(O/'MODAL_CONVERGENCE.csv',out)
 return out
def field_rows(p,q,sc,z,raw,g,h):
 rs=base.decomposition(p,q,sc,z,raw)
 return [dict(**r,coarse_mesh=g,fine_mesh=h,status='AVAILABLE') for r in rs]
def ca_one(p,sc,z,g):
 u=base.pair_vectors(p,PAIRS,RAW)
 noise=.01*np.sqrt(4*(1-np.sum(abs(u)**4,axis=-1)))
 return u,noise
def candidate_rows(sc,p,q,mods):
 a=[mods[80001,f] for f in RAW];b=[mods[160001,f] for f in RAW]
 rs=base.candidate_rows(p,q,a,b,sc,PAIRS,RAW)
 u,n=ca_one(q,sc,200,160001)
 for r in rs:
  r['separation_80001']=r.pop('separation_20001');r['separation_160001']=r.pop('separation_40001')
  r['rank_80001']=r.pop('rank_20001');r['rank_160001']=r.pop('rank_40001')
  ids=[i for i,pair in enumerate(PAIRS) if r['delta_hz']=='JOINT' or pair[1]==r['delta_hz']]
  node=r['node'];nn=n[:,:,ids] if node=='DUAL_NONCOHERENT' else n[:,int(node),ids]
  rms=float(np.sqrt(np.mean(nn**2)))
  r.update(ideal_CA_noise_1pct_rms=rms,ideal_CA_noise_5pct_rms=5*rms,marginal_separation_over_1pct=r['separation_160001']/max(rms,1e-300),marginal_separation_over_5pct=r['separation_160001']/max(5*rms,1e-300),status='AVAILABLE')
 return rs
def fixed_ranges(mods):
 out=[];r=np.array([50000.,50007.,55000.,55007.,60000.,60007.])
 for f in RAW:
  for z in [180.,190.,199.,200.,201.,210.,220.]:
   for g,h in TRANS:
    row=dict(frequency_hz=f,source_depth_m=z,coarse_mesh=g,fine_mesh=h,status='MISSING',raw_relative='',common_removed_relative='',coarse_norm='',fine_norm='')
    if (g,f) in mods and (h,f) in mods:
     p=old.pressure(mods[g,f],r,z);q=old.pressure(mods[h,f],r,z);c=np.vdot(p,q)/np.vdot(p,p)
     row.update(status='AVAILABLE',raw_relative=float(np.linalg.norm(q-p)/np.linalg.norm(q)),common_removed_relative=float(np.linalg.norm(q-c*p)/np.linalg.norm(q)),coarse_norm=float(np.linalg.norm(p)),fine_norm=float(np.linalg.norm(q)))
    out.append(row)
 write(O/'REGISTERED_SIX_RANGE_CONVERGENCE.csv',out)
 return out
def execute():
 if (O/'EXECUTION_STARTED.json').exists():raise RuntimeError('single execution only')
 assert git('log','-1','--pretty=%s')=='R4 E2: freeze paired-frequency numerical recovery'
 assert git('ls-remote','origin','refs/heads/main').split()[0]==git('rev-parse','HEAD')
 f=read(O/'DESIGN_FREEZE.json')
 for p,r in f['bindings'].items():assert sha(p,r['kind']=='RAW')==r['sha256'],p
 dump(O/'EXECUTION_STARTED.json',dict(design_SHA=git('rev-parse','HEAD'),epoch=time.time(),new_MC=0,new_audio=0))
 start=time.monotonic();logs,mods=generate(start);md=mode_diagnostics(mods);six=fixed_ranges(mods)
 complete=all((g,ff) in mods for g in GRIDS for ff in RAW)
 fields=[];cas=[];energy=[];candidates=[];by_scene={}
 for sc in read(old.OUT/'E2_SCENE_FREEZE.json')['scenes']:
  for z in [180.,200.,220.]:
   ps={g:old.field_state([mods[g,ff] for ff in RAW],sc['state'],sc['mirror'],z) for g in GRIDS if all((g,ff) in mods for ff in RAW)}
   us={};ns={}
   for g,p in ps.items():
    np.savez_compressed(O/f"PRESSURE_{sc['scene_id']}_z{int(z)}_n{g}.npz",pressure=p)
    us[g],ns[g]=ca_one(p,sc,z,g)
    for ti,t in enumerate(old.TIMES):
     for node in [0,1]:
      for fi,ff in enumerate(RAW):
       norm=float(np.linalg.norm(p[ti,node,:,fi]))
       energy.append(dict(scene_id=sc['scene_id'],source_depth_m=z,time_s=int(t),node=node,mesh=g,frequency_hz=ff,array_norm=norm,array_power=norm**2,min_element_abs=float(abs(p[ti,node,:,fi]).min()),power_percentile=0.,power_over_frequency_median=0.,relative_noise_sigma1=.01,absolute_source_units='MODEL_UNIT_SOURCE_NO_RECEIVED_SL_CALIBRATION'))
   for g,h in TRANS:
    if g not in ps or h not in ps:continue
    fields.extend(field_rows(ps[g],ps[h],sc,z,RAW,g,h))
    e=base.ca_distance(us[g],us[h]);ov=np.sum(us[g].conj()*us[h],axis=-1)
    for ti,t in enumerate(old.TIMES):
     for node in [0,1]:
      for j,(fc,delta,lo,hi) in enumerate(PAIRS):
       l=RAW.index(lo);hh=RAW.index(hi)
       ag=ps[g][ti,node,:,hh]*ps[g][ti,node,:,l].conj();ah=ps[h][ti,node,:,hh]*ps[h][ti,node,:,l].conj()
       gain=max(1/np.linalg.norm(ag),1/np.linalg.norm(ah))
       phase=np.angle(us[h][ti,node,j,1:]/us[h][ti,node,j,0]*np.conj(us[g][ti,node,j,1:]/us[g][ti,node,j,0]))
       cas.append(dict(scene_id=sc['scene_id'],source_depth_m=z,time_s=int(t),node=node,center_hz=fc,delta_hz=delta,low_hz=lo,high_hz=hi,coarse_mesh=g,fine_mesh=h,status='AVAILABLE',CA_Frobenius=float(e[ti,node,j]),principal_angle_sin=float(e[ti,node,j]/np.sqrt(2)),within_node_phase_rms_rad=float(np.sqrt(np.mean(phase**2))),dual_noncoherent_CA_rms=float(np.sqrt(np.mean(e[ti,:,j]**2))),ideal_noise_1pct=float(ns[h][ti,node,j]),ideal_noise_5pct=float(5*ns[h][ti,node,j]),numerical_over_noise1=float(e[ti,node,j]/max(ns[h][ti,node,j],1e-300)),numerical_over_noise5=float(e[ti,node,j]/max(5*ns[h][ti,node,j],1e-300)),max_normalization_gain=float(gain),power_low_hz=float(np.linalg.norm(ps[h][ti,node,:,l])**2),power_high_hz=float(np.linalg.norm(ps[h][ti,node,:,hh])**2),minimum_power_percentile=0.,cross_node_coherent_phase_used=False))
   if z==200 and 80001 in ps and 160001 in ps:candidates.extend(candidate_rows(sc,ps[80001],ps[160001],mods))
  if time.monotonic()-start>3600:raise TimeoutError('frozen total analysis budget exhausted')
  print('SCENE',sc['scene_id'],flush=True)
 # Equal-weight descriptive registered power distribution per SAME frequency and mesh, not a calibrated SNR.
 for g in GRIDS:
  for ff in RAW:
   rs=[r for r in energy if r['mesh']==g and r['frequency_hz']==ff]
   if not rs:continue
   vals=np.array([r['array_power'] for r in rs]);med=float(np.median(vals))
   for r in rs:
    r['power_percentile']=float(100*(np.sum(vals<r['array_power'])+.5*np.sum(vals==r['array_power']))/len(vals));r['power_over_frequency_median']=r['array_power']/max(med,1e-300)
 en={(r['scene_id'],r['source_depth_m'],r['time_s'],r['node'],r['mesh'],r['frequency_hz']):r for r in energy}
 for r in cas:
  r['minimum_power_percentile']=min(en[r['scene_id'],r['source_depth_m'],r['time_s'],r['node'],r['fine_mesh'],ff]['power_percentile'] for ff in [r['low_hz'],r['high_hz']])
  r['low_signal_level_flag']=r['minimum_power_percentile']<=5
  r['sensitivity_interpretation']='LOW_SIGNAL_LEVEL_SENSITIVITY_ASSOCIATION_ONLY' if r['low_signal_level_flag'] else 'SPATIAL_MODAL_SHAPE_SENSITIVITY_NOT_UNIQUE_CAUSE'
 fg={}
 for r in fields:fg.setdefault((r['scene_id'],r['source_depth_m'],r['time_s'],r['node'],r['frequency_hz']),{})[r['coarse_mesh']]=r
 for v in fg.values():
  for g,h in TRANS:
   if g not in v:continue
   previous={40001:20001,80001:40001}.get(g)
   for metric in ['raw_relative','common_removed_relative']:
    v[g][metric+'_convergence_ratio']=float(v[g][metric]/max(v[previous][metric],1e-300)) if previous in v else ''
 if fields:write(O/'THREE_GRID_FIELD_CONVERGENCE.csv',fields)
 else:write(O/'THREE_GRID_FIELD_CONVERGENCE.csv',[dict(status='NOT_EVALUATED_PROVIDER_UNAVAILABLE')])
 if cas:write(O/'FOUR_GRID_CA_CONVERGENCE.csv',cas)
 else:write(O/'FOUR_GRID_CA_CONVERGENCE.csv',[dict(status='NOT_EVALUATED_PROVIDER_UNAVAILABLE')])
 if energy:write(O/'LOW_ENERGY_DIAGNOSTICS.csv',energy)
 else:write(O/'LOW_ENERGY_DIAGNOSTICS.csv',[dict(status='NOT_EVALUATED_PROVIDER_UNAVAILABLE')])
 if candidates:write(O/'CANDIDATE_NUMERICAL_STABILITY.csv',candidates)
 else:write(O/'CANDIDATE_NUMERICAL_STABILITY.csv',[dict(status='NOT_EVALUATED_PROVIDER_UNAVAILABLE',missing_complete_80001_160001=True)])
 grouped={}
 for r in cas:grouped.setdefault((r['scene_id'],r['source_depth_m'],r['time_s'],r['node'],r['center_hz'],r['delta_hz']),{})[r['coarse_mesh']]=r
 for v in grouped.values():
  for g in [20001,40001,80001]:
   if g not in v:continue
   previous={40001:20001,80001:40001}.get(g)
   v[g]['CA_convergence_ratio']=v[g]['CA_Frobenius']/max(v[previous]['CA_Frobenius'],1e-300) if previous in v else ''
 if cas:write(O/'FOUR_GRID_CA_CONVERGENCE.csv',cas)
 monotonic=[v[80001]['CA_Frobenius']<=v[40001]['CA_Frobenius']+1e-12 and v[40001]['CA_Frobenius']<=v[20001]['CA_Frobenius']+1e-12 for v in grouped.values() if all(g in v for g in GRIDS[:-1])]
 finest=[r for r in cas if r['coarse_mesh']==80001]
 anomalies=[r for r in md if r['status']=='AVAILABLE' and (not r['same_count'] or not r['identity_unambiguous'] or r['max_k_imag']>0)]
 candidate_ok=bool(candidates) and all(r['separation_intervals_exclude_zero'] for r in candidates)
 stable=complete and len(monotonic)==3888 and all(monotonic) and candidate_ok and not anomalies
 noise1=bool(finest) and all(r['numerical_over_noise1']<1 for r in finest)
 noise5=bool(finest) and all(r['numerical_over_noise5']<1 for r in finest)
 if not complete:decision='E2_HIGH_FIDELITY_PROVIDER_UNAVAILABLE';nextstep='STOP'
 elif stable and noise1:decision='E2_PAIRWISE_CA_NUMERICAL_RECOVERABILITY_SUPPORTED';nextstep='FULL_BAND_CERTIFICATION_REVIEW'
 elif stable and noise5:decision='E2_PAIRWISE_CA_STABILITY_CONDITIONAL';nextstep='FULL_BAND_CERTIFICATION_REVIEW'
 else:decision='E2_NUMERICAL_RECOVERY_NOT_ESTABLISHED';nextstep='H3'
 result=dict(decision=decision,parent_SHA=PARENT,design_SHA=git('rev-parse','HEAD'),new_KRAKEN_calls=sum(r['new_call'] for r in logs),provider_complete=complete,pair_count=9,raw_frequency_count=13,CA_comparison_rows=len(cas),complete_CA_state_pairs=len(monotonic),monotonically_decreasing_CA_state_pairs=sum(monotonic),finest_CA_worst=max([r['CA_Frobenius'] for r in finest],default=None),finest_noise1_ratio_max=max([r['numerical_over_noise1'] for r in finest],default=None),finest_noise5_ratio_max=max([r['numerical_over_noise5'] for r in finest],default=None),candidate_rows=len(candidates),candidate_rank_changes=sum(r['rank_changed'] for r in candidates),candidate_unresolved=sum(not r['separation_intervals_exclude_zero'] for r in candidates),mode_anomalies=len(anomalies),original_E2_admission='FAIL_UNCHANGED',E2_scientific_information='NOT_EVALUATED',new_Monte_Carlo=0,new_received_audio=0,R4_percent=0,next=nextstep,next_stage_execution='NOT_AUTHORIZED',elapsed_s=time.monotonic()-start,low_energy='RELATIVE_POWER_PERCENTILES_AND_GAIN_RETAINED_NO_FILTERING_NO_SNR_RETUNING',independent_audit='PENDING')
 dump(O/'NUMERICAL_RECOVERABILITY_DECISION.json',result);print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':execute()
