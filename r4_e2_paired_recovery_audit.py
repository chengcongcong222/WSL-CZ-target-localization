"""Independent cold binary/field/CA audit. Read-only forward evidence; zero solver calls."""
from pathlib import Path
import struct,json
import numpy as np
import r4_e2_paired_recovery as p
O=p.O;checks=[]
def ck(name,ok,detail=''):
 checks.append(dict(check=name,PASS=bool(ok),detail=str(detail)))
def parse(path):
 b=Path(path).read_bytes();r=struct.unpack_from('<i',b)[0]*4
 nt,nm=struct.unpack_from('<ii',b,92);M=struct.unpack_from('<i',b,5*r)[0]
 z=np.array(struct.unpack_from('<'+str(nt)+'f',b,4*r),float);phi=np.empty((nm,M),complex)
 for j in range(M):
  v=np.array(struct.unpack_from('<'+str(2*nm)+'f',b,(7+j)*r));phi[:,j]=v[::2]+1j*v[1::2]
 v=np.array(struct.unpack_from('<'+str(2*M)+'f',b,(7+M)*r))
 return z,phi,v[::2]+1j*v[1::2]
def pressure(mod,r,z):
 zz,phi,k=mod;w=phi[list(zz).index(z)]*phi[list(zz).index(200.)]
 rr=np.asarray(r).ravel();term=w[:,None]*np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*(k[:,None]*rr+np.pi/4))
 return np.add.reduce(term,axis=0).reshape(np.shape(r))
def unit(pres):
 lo=[p.RAW.index(x[2]) for x in p.PAIRS];hi=[p.RAW.index(x[3]) for x in p.PAIRS]
 v=np.moveaxis(pres[...,hi]*pres[...,lo].conj(),-1,-2)
 return v/np.sqrt(np.sum(v.real*v.real+v.imag*v.imag,axis=-1))[...,None]
def distance(a,b):
 c=np.einsum('...i,...j->...ij',a,a.conj())-np.einsum('...i,...j->...ij',b,b.conj())
 return np.sqrt(np.sum(abs(c)**2,axis=(-1,-2)))
def geometry(state,mirror):
 tt=np.array([0.,600.,1200.]);dt=np.maximum(tt-600,0)
 main=np.stack([2*np.minimum(tt,600)+2*dt*np.cos(np.pi/12),2*dt*np.sin(np.pi/12)],axis=-1)
 nodes=np.stack([main,main+[0,mirror*5000]],axis=1)
 elems=nodes[:,:,None,:]+np.stack([np.arange(-7.,8.,2.),np.zeros(8)],axis=-1)[None,None,:,:]
 source=np.array(state[:2])+tt[:,None]*state[2:4]
 return np.sqrt(np.sum((source[:,None,None,:]-elems)**2,axis=-1))
def main():
 freeze=p.read(O/'DESIGN_FREEZE.json');dec=p.read(O/'NUMERICAL_RECOVERABILITY_DECISION.json')
 for name,r in freeze['bindings'].items():ck('binding:'+name,p.sha(name,r['kind']=='RAW')==r['sha256'])
 logs=p.base.rows(O/'PAIRED_MODAL_GENERATION.csv');new=[r for r in logs if r['new_call']=='True']
 ck('call_cap',len(new)<=22 and all(float(r['elapsed_s'])<91 for r in new))
 ck('call_count_record',len(new)==dec['new_KRAKEN_calls'])
 ck('no_substitution',all(int(r['mesh']) in [80001,160001] and float(r['frequency_hz']) in p.RAW for r in new))
 failures=[i for i,r in enumerate(logs) if r['status'] not in ['REUSED','GENERATED','NOT_ATTEMPTED_AFTER_PROVIDER_FAILURE']]
 ck('stop_after_first_failure',len(failures)<=1 and (not failures or not any(r['new_call']=='True' for r in logs[failures[0]+1:])))
 mods={}
 for r in logs:
  if r['status'] not in ['REUSED','GENERATED']:continue
  m=parse(r['modal_file']);ref=p.old.parse_mod(r['modal_file'])
  ck('binary:'+r['modal_file'],all(np.array_equal(a,b) for a,b in zip(m,ref)))
  mods[int(r['mesh']),float(r['frequency_hz'])]=m
 grids=[g for g in p.GRIDS if all((g,ff) in mods for ff in p.RAW)]
 ck('provider_completeness',dec['provider_complete']==(len(grids)==4))
 ca=p.base.rows(O/'FOUR_GRID_CA_CONVERGENCE.csv');fields=p.base.rows(O/'THREE_GRID_FIELD_CONVERGENCE.csv');en=p.base.rows(O/'LOW_ENERGY_DIAGNOSTICS.csv');cand=p.base.rows(O/'CANDIDATE_NUMERICAL_STABILITY.csv')
 ca=[r for r in ca if r['status']=='AVAILABLE'];fields=[r for r in fields if r['status']=='AVAILABLE'];cand=[r for r in cand if r['status']=='AVAILABLE']
 expected=24*3*3*2*9*sum(g in grids and h in grids for g,h in p.TRANS)
 ck('all_CA_records',len(ca)==expected);ck('all_field_records',len(fields)==expected//9*13);ck('all_energy_records',len(en)==24*3*3*2*13*len(grids))
 caix={(r['scene_id'],float(r['source_depth_m']),int(r['time_s']),int(r['node']),float(r['center_hz']),int(r['delta_hz']),int(r['coarse_mesh'])):r for r in ca}
 fix={(r['scene_id'],float(r['source_depth_m']),int(r['time_s']),int(r['node']),float(r['frequency_hz']),int(r['coarse_mesh'])):r for r in fields}
 cix={(r['scene_id'],r['candidate_id'],r['node'],r['delta_hz']):r for r in cand}
 coldmax=0.;camax=0.;noisemax=0.;fieldmax=0.;candidate_max=0.
 for sc in p.read(p.old.OUT/'E2_SCENE_FREEZE.json')['scenes']:
  ranges=geometry(sc['state'],sc['mirror']);ck('independent_geometry:'+sc['scene_id'],np.max(abs(ranges-p.old.geometry(sc['state'],sc['mirror'])[0]))<1e-9)
  for z in [180.,200.,220.]:
   ps={}
   for g in grids:
    pres=np.load(O/f"PRESSURE_{sc['scene_id']}_z{int(z)}_n{g}.npz")['pressure'];ps[g]=pres
    for j,ff in enumerate(p.RAW):
     cold=pressure(mods[g,ff],ranges,z);err=np.max(abs(cold-pres[...,j]))/max(np.max(abs(cold)),1e-300);coldmax=max(coldmax,float(err))
   us={g:unit(pr) for g,pr in ps.items()}
   for g,h in p.TRANS:
    if g not in grids or h not in grids:continue
    e=distance(us[g],us[h]);ns=.01*np.sqrt(4*(1-np.sum(abs(us[h])**4,axis=-1)))
    a=ps[g];b=ps[h];c=np.sum(a.conj()*b,axis=2)/np.sum(abs(a)**2,axis=2)
    raw=np.linalg.norm(b-a,axis=2)/np.linalg.norm(b,axis=2);res=np.linalg.norm(b-c[:,:,None,:]*a,axis=2)/np.linalg.norm(b,axis=2)
    for t in range(3):
     for n in range(2):
      for j,pair in enumerate(p.PAIRS):
       rr=caix[sc['scene_id'],z,int(p.old.TIMES[t]),n,pair[0],pair[1],g]
       camax=max(camax,abs(float(rr['CA_Frobenius'])-e[t,n,j]));noisemax=max(noisemax,abs(float(rr['ideal_noise_1pct'])-ns[t,n,j]))
      for j,ff in enumerate(p.RAW):
       rr=fix[sc['scene_id'],z,int(p.old.TIMES[t]),n,ff,g]
       fieldmax=max(fieldmax,abs(float(rr['raw_relative'])-raw[t,n,j]),abs(float(rr['common_removed_relative'])-res[t,n,j]))
   if z==200 and 80001 in grids and 160001 in grids:
    for name,state,cz,scope in p.old.candidates(sc)[-8:]:
     cr=geometry(state,sc['mirror'])
     cp=np.stack([pressure(mods[80001,ff],cr,cz) for ff in p.RAW],axis=-1)
     cq=np.stack([pressure(mods[160001,ff],cr,cz) for ff in p.RAW],axis=-1)
     u,v=us[80001],us[160001];w,x=unit(cp),unit(cq)
     vv=[distance(u,w),distance(v,x),distance(u,v),distance(w,x)]
     for node in ['0','1','DUAL_NONCOHERENT']:
      for delta in ['2','5','10','JOINT']:
       ids=[j for j,pp in enumerate(p.PAIRS) if delta=='JOINT' or str(pp[1])==delta]
       nums=[float(np.sqrt(np.mean((ar[:,:,ids] if node=='DUAL_NONCOHERENT' else ar[:,int(node),ids])**2))) for ar in vv]
       rr=cix[sc['scene_id'],name,node,delta]
       for key,num in zip(['separation_80001','separation_160001','truth_CA_grid_error','candidate_CA_grid_error'],nums):candidate_max=max(candidate_max,abs(float(rr[key])-num))
  print('AUDIT_SCENE',sc['scene_id'],flush=True)
 ck('all_pressure_cold',coldmax<1e-12,coldmax);ck('all_CA_cold',camax<1e-12,camax);ck('all_field_decomposition',fieldmax<1e-12,fieldmax);ck('all_noise_marginals',noisemax<1e-14,noisemax)
 ck('all_candidate_cold',candidate_max<1e-11,candidate_max)
 for g in grids:
  for ff in p.RAW:
   rs=[r for r in en if int(r['mesh'])==g and float(r['frequency_hz'])==ff];vals=np.array([float(r['array_power']) for r in rs])
   err=max(abs(float(r['power_percentile'])-100*(sum(vals<float(r['array_power']))+.5*sum(vals==float(r['array_power'])))/len(vals)) for r in rs)
   ck('power_percentile:'+str((g,ff)),err<1e-12,err)
 groups={}
 for r in cand:groups.setdefault((r['scene_id'],r['node'],r['delta_hz']),[]).append(r)
 for key,rs in groups.items():
  for field,rank in [('separation_80001','rank_80001'),('separation_160001','rank_160001')]:
   ck('candidate_rank:'+str(key)+rank,all(int(r[rank])==i+1 for i,r in enumerate(sorted(rs,key=lambda r:(float(r[field]),r['candidate_id'])))))
  ck('candidate_triangle:'+str(key),all(float(r['absolute_separation_change'])<=float(r['triangle_bound'])+1e-12 for r in rs))
 finest=[r for r in ca if int(r['coarse_mesh'])==80001];cg={}
 for r in ca:cg.setdefault((r['scene_id'],r['source_depth_m'],r['time_s'],r['node'],r['center_hz'],r['delta_hz']),{})[int(r['coarse_mesh'])]=float(r['CA_Frobenius'])
 mono=[v[80001]<=v[40001]+1e-12 and v[40001]<=v[20001]+1e-12 for v in cg.values() if len(v)==3]
 md=p.base.rows(O/'MODAL_CONVERGENCE.csv');anomaly=any(r['status']=='AVAILABLE' and (r['same_count']=='False' or r['identity_unambiguous']=='False' or float(r['max_k_imag'])>0) for r in md)
 stable=len(grids)==4 and len(mono)==3888 and all(mono) and bool(cand) and all(r['separation_intervals_exclude_zero']=='True' for r in cand) and not anomaly
 n1=bool(finest) and all(float(r['numerical_over_noise1'])<1 for r in finest);n5=bool(finest) and all(float(r['numerical_over_noise5'])<1 for r in finest)
 decision='E2_HIGH_FIDELITY_PROVIDER_UNAVAILABLE' if len(grids)<4 else 'E2_PAIRWISE_CA_NUMERICAL_RECOVERABILITY_SUPPORTED' if stable and n1 else 'E2_PAIRWISE_CA_STABILITY_CONDITIONAL' if stable and n5 else 'E2_NUMERICAL_RECOVERY_NOT_ESTABLISHED'
 ck('decision_reconstructed',dec['decision']==decision)
 ck('no_scientific_credit',dec['original_E2_admission']=='FAIL_UNCHANGED' and dec['E2_scientific_information']=='NOT_EVALUATED' and dec['R4_percent']==dec['new_Monte_Carlo']==dec['new_received_audio']==0)
 p.write(O/'INDEPENDENT_AUDIT_CHECKS.csv',checks)
 val=dict(checks=len(checks),PASS=sum(r['PASS'] for r in checks),FAIL=sum(not r['PASS'] for r in checks),cold_pressure_max_relative_difference=coldmax,CA_max_difference=camax,field_decomposition_max_difference=fieldmax,noise_max_difference=noisemax,candidate_max_difference=candidate_max,scope='Independent struct modal parser; all archived states/frequencies/grids cold pressure; all CA and noise; all192 candidates if provider complete; registered power ranks; outcome independently rebuilt.',new_solver_calls=0,scientific_gate='FAIL_UNCHANGED')
 p.dump(O/'VALIDATION.json',val);print(json.dumps(val,indent=2));assert val['FAIL']==0
if __name__=='__main__':main()
