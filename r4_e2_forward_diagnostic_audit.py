"""Read-only independent reconstruction audit; no external solver or noise samples."""
from pathlib import Path
import csv,json,struct,hashlib
import numpy as np
import r4_e2_forward_diagnostic as d
O=d.OUT
checks=[]
def ck(name,ok,detail=''):
 checks.append(dict(check=name,PASS=bool(ok),detail=str(detail)))
def modal(p):
 b=Path(p).read_bytes();r=struct.unpack_from('<i',b)[0]*4
 nt,nm=struct.unpack_from('<ii',b,92)
 M=struct.unpack_from('<i',b,5*r)[0]
 z=np.array(struct.unpack_from('<'+str(nt)+'f',b,4*r),float)
 phi=np.empty((nm,M),complex)
 for j in range(M):
  vals=np.array(struct.unpack_from('<'+str(2*nm)+'f',b,(7+j)*r))
  phi[:,j]=vals[::2]+1j*vals[1::2]
 vals=np.array(struct.unpack_from('<'+str(2*M)+'f',b,(7+M)*r))
 return z,phi,vals[::2]+1j*vals[1::2]
def pressure(mod,r,z):
 zz,phi,k=mod;w=phi[list(zz).index(z)]*phi[list(zz).index(200.)]
 rr=np.asarray(r).ravel()
 terms=w[:,None]*np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*(k[:,None]*rr+np.pi/4))
 return np.add.reduce(terms,axis=0).reshape(np.shape(r))
def pair(p,lo,hi):
 a=np.moveaxis(p[...,hi]*np.conj(p[...,lo]),-1,-2)
 return a/np.sqrt(np.sum(a.real*a.real+a.imag*a.imag,axis=-1))[...,None]
def dist(a,b):
 ac=np.einsum('...i,...j->...ij',a,a.conj());bc=np.einsum('...i,...j->...ij',b,b.conj())
 return np.sqrt(np.sum(abs(ac-bc)**2,axis=(-1,-2)))
f=d.read(O/'FORWARD_DIAGNOSTIC_FREEZE.json')
for p,r in f['bindings'].items():ck('frozen:'+p,d.sha(p,r['kind']=='RAW')==r['sha256'])
freq=d.read(d.PRIOR/'E2_REQUIRED_FREQUENCIES.json');raw=freq['raw_frequencies_hz'];pairs=freq['frequency_pairs']
lo=[raw.index(x['low_hz']) for x in pairs];hi=[raw.index(x['high_hz']) for x in pairs]
mods={}
for p in list(d.MD.glob('f*_s*_n*.mod'))+list((O/'third_mesh').glob('*.mod')):
 a=modal(p);b=d.old.parse_mod(p)
 ck('independent_binary:'+str(p),all(np.array_equal(x,y) for x,y in zip(a,b)))
 if '_s+0_' in p.name:mods[p.name]=a
calls=d.rows(O/'THIRD_MESH_CALLS.csv')
ck('four_bounded_calls',len(calls)==4 and [float(x['frequency_hz']) for x in calls]==[150,200,244,250] and all(float(x['elapsed_s'])<91 for x in calls))
ck('third_frequencies_no_CA_partner',not any(x['low_hz'] in [150,200,244,250] and x['high_hz'] in [150,200,244,250] for x in pairs))
ca=d.rows(O/'CA_MESH_STABILITY.csv');noise=d.rows(O/'CA_NUMERICAL_VS_NOISE.csv');de=d.rows(O/'COMPLEX_FIELD_ERROR_DECOMPOSITION.csv')
ck('CA_row_count',len(ca)==24*3*3*2*569)
ck('noise_row_count',len(noise)==len(ca));ck('decomposition_count',len(de)==24*3*3*2*201)
candidate_table=d.rows(O/'CA_CANDIDATE_STABILITY.csv')
candidate_index={(r['scene_id'],r['candidate_id'],r['node'],r['delta_hz']):r for r in candidate_table}
candidate_max=0.
caidx={(r['scene_id'],int(float(r['source_depth_m'])),int(r['time_s']),int(r['node']),float(r['center_hz']),int(r['delta_hz'])):r for r in ca}
noidx={(r['scene_id'],int(float(r['source_depth_m'])),int(r['time_s']),int(r['node']),float(r['center_hz']),int(r['delta_hz'])):r for r in noise}
coldmax=0.;camax=0.;noisemax=0.
for sc in d.read(d.PRIOR/'E2_SCENE_FREEZE.json')['scenes']:
 for zs in [180,200,220]:
  ar=np.load(O/f"PRESSURE_{sc['scene_id']}_z{zs}.npz");p=ar['coarse'];q=ar['fine']
  ranges,_=d.old.geometry(sc['state'],sc['mirror'])
  for grid,arr in [(20001,p),(40001,q)]:
   for ff in [150.,200.,244.,250.]:
    cold=pressure(mods[f'f{int(ff*2):04d}_s+0_n{grid}.mod'],ranges,zs)
    err=np.max(abs(cold-arr[...,raw.index(ff)]))/max(np.max(abs(cold)),1e-300);coldmax=max(coldmax,float(err))
  u=pair(p,lo,hi);v=pair(q,lo,hi);e=dist(u,v)
  nr=.01*np.sqrt(4*(1-np.sum(abs(u)**4,axis=-1)))
  for t in range(3):
   for node in range(2):
    for k,x in enumerate(pairs):
     key=(sc['scene_id'],zs,int(d.old.TIMES[t]),node,x['center_hz'],x['delta_hz'])
     camax=max(camax,abs(float(caidx[key]['CA_Frobenius'])-e[t,node,k]))
     noisemax=max(noisemax,abs(float(noidx[key]['CA_noise_rms_1pct'])-nr[t,node,k]))
  if zs==200:
   for name,state,cz,scope in d.old.candidates(sc)[-8:]:
    cr,_=d.old.geometry(state,sc['mirror'])
    cp=np.stack([pressure(mods[f'f{int(ff*2):04d}_s+0_n20001.mod'],cr,cz) for ff in raw],axis=-1)
    cq=np.stack([pressure(mods[f'f{int(ff*2):04d}_s+0_n40001.mod'],cr,cz) for ff in raw],axis=-1)
    w=pair(cp,lo,hi);x=pair(cq,lo,hi)
    values=[dist(u,w),dist(v,x),dist(u,v),dist(w,x)]
    for node in ['0','1','DUAL_NONCOHERENT']:
     for delta in ['2','5','10','JOINT']:
      inds=[j for j,pp in enumerate(pairs) if delta=='JOINT' or str(pp['delta_hz'])==delta]
      nums=[float(np.sqrt(np.mean((vv[:,:,inds] if node=='DUAL_NONCOHERENT' else vv[:,int(node),inds])**2))) for vv in values]
      rr=candidate_index[(sc['scene_id'],name,node,delta)]
      keys=['separation_20001','separation_40001','truth_CA_grid_error','candidate_CA_grid_error']
      candidate_max=max(candidate_max,max(abs(float(rr[kk])-vv) for kk,vv in zip(keys,nums)))
  ck('finite_archived_pressure:'+sc['scene_id']+':'+str(zs),np.all(np.isfinite(p)) and np.all(np.isfinite(q)))
 print('AUDIT_SCENE '+sc['scene_id'],flush=True)
ck('all_candidate_pressure_and_CA_reconstruction',candidate_max<1e-11,candidate_max)
ck('cold_pressure_reconstruction',coldmax<1e-12,coldmax)
ck('all_CA_reconstruction',camax<1e-13,camax);ck('all_noise_marginals',noisemax<1e-14,noisemax)
orthmax=max(abs(float(r['raw_relative'])**2-float(r['common_removed_relative'])**2-float(r['common_component_relative'])**2) for r in de)
ck('common_fit_orthogonal_energy',orthmax<1e-12,orthmax)
shared=d.rows(O/'SHARED_FREQUENCY_NOISE_COVARIANCE.csv')
ck('shared_covariance_controls',len(shared)==432 and all(float(r['shared_cross_covariance_frobenius'])>0 and float(r['disjoint_cross_covariance_frobenius'])==0 and float(r['marginal_formula_max_error'])<1e-17 for r in shared))
candidates=d.rows(O/'CA_CANDIDATE_STABILITY.csv');ck('candidate_count',len(candidates)==2304)
tri=max(float(r['absolute_separation_change'])-float(r['triangle_bound']) for r in candidates)
ck('candidate_triangle_bound',tri<1e-12,tri)
groups={}
for r in candidates:groups.setdefault((r['scene_id'],r['node'],r['delta_hz']),[]).append(r)
for key,rs in groups.items():
 for sep,rank in [('separation_20001','rank_20001'),('separation_40001','rank_40001')]:
  ck('candidate_rank:'+str(key)+':'+rank,all(int(x[rank])==i+1 for i,x in enumerate(sorted(rs,key=lambda r:(float(r[sep]),r['candidate_id'])))))
dec=d.read(O/'FORWARD_FIDELITY_DECISION.json')
ck('no_scientific_stage',dec['E2_scientific_information']=='NOT_EVALUATED' and dec['new_MC']==dec['new_received_audio']==dec['R4_percent']==0 and dec['original_E2_G0_admission']=='FAIL_UNCHANGED')
d.write(O/'INDEPENDENT_AUDIT_CHECKS.csv',checks)
validation=dict(checks=len(checks),PASS=sum(r['PASS'] for r in checks),FAIL=sum(not r['PASS'] for r in checks),cold_pressure_max_relative_difference=coldmax,CA_max_reconstruction_difference=camax,noise_max_reconstruction_difference=noisemax,orthogonal_energy_max_difference=orthmax,scope='Independent struct parser all804+4 MOD; cold pressure all72 scenes/depth arrays at4frequencies; all245808 CA/noise rows; all candidate rank/triangle checks. No new solver.',candidate_forward_recompute='ALL192_CANDIDATES_BOTH_GRIDS_FULL201_FREQUENCIES_COLD_MODAL_SUM',candidate_max_difference=candidate_max,scientific_admission='FAIL_UNCHANGED')
d.dump(O/'VALIDATION.json',validation);print(json.dumps(validation,indent=2))
assert validation['FAIL']==0
