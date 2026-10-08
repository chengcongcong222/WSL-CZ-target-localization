"""Frozen forward-fidelity attribution only. No information matrix, optimizer, audio or RNG."""
from pathlib import Path
import json,csv,hashlib,subprocess,time,os,struct,re
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
import r4_e2_g0 as old
ROOT=Path.cwd()
OUT=Path('results/R4_E2_FORWARD_FIDELITY_DIAGNOSTIC')
PRIOR=old.OUT
MD=PRIOR/'modal'
PARENT='4605cf248f3f51238b7807563a8be60dfe66339f'
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p,raw=False):
 b=Path(p).read_bytes();return hashlib.sha256(b if raw else b.replace(b'\r\n',b'\n')).hexdigest()
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rs):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]));w.writeheader();w.writerows(rs)
def git(*a):return subprocess.check_output(['git',*a],text=True,encoding='utf-8').strip()
def canonical_ca(p,lo,hi):
 a=p[...,hi]*p[...,lo].conj();n=np.sqrt(np.sum(abs(a)**2,axis=-1))
 u=a/n[...,None];return u,n
def ca_distance(u,v):
 # Frobenius distance of uu^H-vv^H; stable direct outer product, not a cancellation-prone 1-|overlap|^2.
 c=u[..., :,None]*u[...,None,:].conj();d=v[..., :,None]*v[...,None,:].conj()
 return np.linalg.norm(c-d,axis=(-1,-2))
def pressure_components(a,b,r,z):
 az,ap,ak=a;bz,bp,bk=b
 idx=lambda v,d:int(np.flatnonzero(v==d)[0])
 w0=ap[idx(az,z)]*ap[idx(az,200.)];w1=bp[idx(bz,z)]*bp[idx(bz,200.)]
 def calc(w,k):
  rr=np.asarray(r).reshape(-1)
  return np.sum(w[:,None]*np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*k[:,None]*rr-1j*np.pi/4),axis=0).reshape(np.shape(r))
 p0=calc(w0,ak);p1=calc(w1,bk)
 if len(ak)!=len(bk):return p0,p1,None,None,None
 gap=np.minimum(np.r_[np.inf,abs(np.diff(ak.real))],np.r_[abs(np.diff(ak.real)),np.inf])
 corr=abs(np.sum(ap.conj()*bp,axis=0))/np.maximum(np.linalg.norm(ap,axis=0)*np.linalg.norm(bp,axis=0),1e-300)
 compatible=bool(np.all(abs(bk.real-ak.real)<.45*gap) and np.all(corr>.99))
 if not compatible:return p0,p1,None,None,None
 return p0,p1,calc(w0,bk)-p0,calc(w1,ak)-p0,compatible
def artifact_format_audit():
 out=[]
 for path in sorted(MD.glob('f*_s*_n*.mod')):
  mt=re.fullmatch(r'f(\d+)_s([+-]\d+)_n(\d+)\.mod',path.name)
  if not mt:continue
  f=int(mt[1])/2;mesh=int(mt[3]);prt=path.with_suffix('.prt').read_text()
  fre=re.search(r'Frequency\s*=\s*([\d.E+-]+)',prt)
  requested=re.search(r'# mesh points\s*=\s*(\d+)',prt)
  z,phi,k=old.parse_mod(path)
  assert fre and abs(float(fre[1])-f)<1e-6 and requested and int(requested[1])==mesh
  assert np.array_equal(z,old.DEPTHS) and np.all(np.isfinite(phi)) and np.all(np.isfinite(k))
  out.append(dict(file=str(path),frequency_hz=f,ssp_variant=int(mt[2]),requested_mesh=mesh,printed_input_mesh=int(requested[1]),returned_modes=len(k),depth_count=len(z),k_imag_max=float(k.imag.max()),phase_convention='EXP_MINUS_I_KR_WITH_COMPLEX_K',attenuation='THORP_DB_PER_WAVELENGTH',internal_final_mesh='NOT_EXPLICITLY_CERTIFIED_BY_PRT',format_valid=True))
 assert len(out)==804
 write(OUT/'HISTORICAL_MODAL_FORMAT_AUDIT.csv',out)
def aggregate_history():
 r=rows(PRIOR/'E2_FORWARD_CONVERGENCE.csv');ind=rows(PRIOR/'E2_INDEPENDENT_FORWARD_RECONSTRUCTION.csv')
 assert len(r)==len(ind)==1416
 out=[]
 for low,high in ((150,174.5),(175,199.5),(200,224.5),(225,249.5),(250,250)):
  q=[x for x in r if x['kind']=='mesh40001_vs20001' and low<=float(x['frequency_hz'])<=high]
  out.append(dict(frequency_low_hz=low,frequency_high_hz=high,checks=len(q),FAIL=sum(float(x['error'])>.002 for x in q),max_error=max(float(x['error']) for x in q),median_error=float(np.median([float(x['error']) for x in q]))))
 write(OUT/'HISTORICAL_FREQUENCY_ERROR.csv',out)
def modes(raw):
 m0=[];m1=[];at=[];dist=[]
 r=np.array([50000.,50007.,55000.,55007.,60000.,60007.])
 for f in raw:
  a=old.parse_mod(MD/f'f{int(f*2):04d}_s+0_n20001.mod');b=old.parse_mod(MD/f'f{int(f*2):04d}_s+0_n40001.mod')
  m0.append(a);m1.append(b)
  z,phi,k=a;zz,p2,k2=b
  equal=len(k)==len(k2)
  corr=abs(np.sum(phi.conj()*p2,axis=0))/np.maximum(np.linalg.norm(phi,axis=0)*np.linalg.norm(p2,axis=0),1e-300) if equal else np.array([0.])
  sign=np.sign(np.sum(phi.real*p2.real,axis=0)) if equal else np.array([0.])
  delta=abs(k-k2) if equal else np.array([0.])
  at.append(dict(frequency_hz=f,coarse_modes=len(k),fine_modes=len(k2),same_count=equal,max_abs_delta_k=float(np.max(delta)),max_phase_difference_60km_rad=float(np.max(abs(k.real-k2.real))*60000) if equal else '',min_sampled_mode_correlation=float(np.min(corr)),sign_reversals=int(np.sum(sign<0)),sampled_phi_relative=float(np.linalg.norm(phi-p2*sign)/np.linalg.norm(phi)) if equal else '',k_imag_max=float(max(k.imag.max(),k2.imag.max())),modal_payload='COMPLEX_FLOAT32',full_depth_normalization_check='NOT_AVAILABLE_13_OUTPUT_DEPTHS_ONLY'))
  for zs in (180.,190.,199.,200.,201.,210.,220.):
   p,q,pk,pw,match=pressure_components(a,b,r,zs)
   c=np.vdot(p,q)/np.vdot(p,p);e=float(np.linalg.norm(q-c*p)/np.linalg.norm(q))
   for j,rr in enumerate(r):
    dist.append(dict(frequency_hz=f,source_depth_m=zs,range_m=rr,coarse_abs=float(abs(p[j])),fine_abs=float(abs(q[j])),raw_point_relative=float(abs(q[j]-p[j])/max(abs(q[j]),1e-300)),common_factor_re=float(c.real),common_factor_im=float(c.imag),residual_point_relative=float(abs(q[j]-c*p[j])/max(abs(q[j]),1e-300)),vector_residual=e))
   denom=max(np.linalg.norm(q-p),1e-300)
   at.append(dict(frequency_hz=f,coarse_modes=len(k),fine_modes=len(k2),same_count=equal,max_abs_delta_k=float(np.max(delta)),max_phase_difference_60km_rad='',min_sampled_mode_correlation=float(np.min(corr)),sign_reversals=int(np.sum(sign<0)),sampled_phi_relative='',k_imag_max=float(max(k.imag.max(),k2.imag.max())),modal_payload='source_depth_'+str(zs),full_depth_normalization_check='WEIGHT_SWAP_CONDITIONAL_'+str(match)))
   # Separate hybrid coefficients; no claim of additive fractions or causal percentages.
   if pk is not None:
    hybrids.append(dict(frequency_hz=f,source_depth_m=zs,deltaP_norm=float(denom),k_swap_norm_over_delta=float(np.linalg.norm(pk)/denom),phi_weight_swap_norm_over_delta=float(np.linalg.norm(pw)/denom),interaction_norm_over_delta=float(np.linalg.norm(q-p-pk-pw)/denom),matched=True))
   else:hybrids.append(dict(frequency_hz=f,source_depth_m=zs,deltaP_norm=float(denom),k_swap_norm_over_delta='',phi_weight_swap_norm_over_delta='',interaction_norm_over_delta='',matched=False))
 write(OUT/'MODAL_ERROR_ATTRIBUTION.csv',at);write(OUT/'HISTORICAL_POINTWISE_ERRORS.csv',dist);write(OUT/'MODAL_HYBRID_ATTRIBUTION.csv',hybrids)
 return m0,m1
def decomposition(p,q,sc,z,raw):
 out=[];offset=np.arange(-7.,8.,2.)
 for t in range(3):
  for node in range(2):
   for fi,f in enumerate(raw):
    a=p[t,node,:,fi];b=q[t,node,:,fi]
    c=np.vdot(a,b)/np.vdot(a,a);res=b-c*a
    rawerr=float(np.linalg.norm(b-a)/np.linalg.norm(b));de=float(np.linalg.norm(res)/np.linalg.norm(b))
    common=float(np.linalg.norm((c-1)*a)/np.linalg.norm(b))
    w=abs(b)**2;w=w/np.sum(w)
    amp=np.log(np.maximum(abs(b/a),1e-300));amp-=np.sum(w*amp)
    phase=np.angle(b/(c*a))
    slope=float(np.polyfit(offset,np.unwrap(np.angle(b/a)),1)[0])
    out.append(dict(scene_id=sc['scene_id'],source_depth_m=z,time_s=int(old.TIMES[t]),node=node,frequency_hz=f,raw_relative=rawerr,common_removed_relative=de,common_component_relative=common,common_energy_fraction=common**2/max(rawerr**2,1e-300),relative_amp_log_rms=float(np.sqrt(np.sum(w*amp**2))),relative_phase_rms_rad=float(np.sqrt(np.sum(w*phase**2))),phase_slope_error_rad_per_m=slope,min_pressure_abs=float(min(abs(a).min(),abs(b).min())),min_element_fraction=float(min(abs(a).min()/np.linalg.norm(a),abs(b).min()/np.linalg.norm(b))),near_null_flag=bool(min(abs(a).min()/np.linalg.norm(a),abs(b).min()/np.linalg.norm(b))<1e-6)))
 return out
def ca_rows(p,q,sc,z,pairs,raw):
 ix={f:i for i,f in enumerate(raw)};rr=[];nn=[]
 for f,d,fl,fh in pairs:
  u,n=canonical_ca(p,ix[fl],ix[fh]);v,n2=canonical_ca(q,ix[fl],ix[fh])
  error=ca_distance(u,v);overlap=np.sum(u.conj()*v,axis=-1);sin=np.sqrt(np.maximum(0,1-abs(overlap)**2))
  rms1=.01*np.sqrt(4*(1-np.sum(abs(u)**4,axis=-1)))
  rms5=5*rms1
  # Noncoherent two-node block diagonal norm, scaled sqrt2; shared pairs are not independent trials.
  dual=np.sqrt(np.mean(error**2,axis=1))
  dualnoise=np.sqrt(np.mean(rms1**2,axis=1))
  for ti,t in enumerate(old.TIMES):
   for node in (0,1):
    rr.append(dict(scene_id=sc['scene_id'],source_depth_m=z,time_s=int(t),node=node,center_hz=f,delta_hz=d,CA_Frobenius=float(error[ti,node]),principal_angle_sin=float(error[ti,node]/np.sqrt(2)),dual_CA_rms=float(dual[ti]),min_autoproduct_norm=float(min(n[ti,node],n2[ti,node])),within_node_relative_phase_rms=float(np.sqrt(np.mean(np.angle(v[ti,node,1:]/v[ti,node,0]*np.conj(u[ti,node,1:]/u[ti,node,0]))**2))),cross_node_coherent_phase_used=False))
    nn.append(dict(scene_id=sc['scene_id'],source_depth_m=z,time_s=int(t),node=node,center_hz=f,delta_hz=d,CA_error=float(error[ti,node]),CA_noise_rms_1pct=float(rms1[ti,node]),CA_noise_rms_5pct=float(rms5[ti,node]),numerical_over_1pct_noise=float(error[ti,node]/max(rms1[ti,node],1e-300)),numerical_over_5pct_noise=float(error[ti,node]/max(rms5[ti,node],1e-300)),dual_numerical_over_1pct_noise=float(dual[ti]/max(dualnoise[ti],1e-300)),aggregation='MARGINAL_RMS_NO_INDEPENDENT_PAIR_COUNT'))
 return rr,nn
def pair_vectors(p,pairs,raw):
 ix={f:i for i,f in enumerate(raw)}
 aa=p[..., [ix[x[3]] for x in pairs]]*p[..., [ix[x[2]] for x in pairs]].conj()
 aa=np.moveaxis(aa,-1,-2)
 return aa/np.linalg.norm(aa,axis=-1)[...,None]
def candidate_rows(p,q,a,b,sc,pairs,raw):
 u=pair_vectors(p,pairs,raw);v=pair_vectors(q,pairs,raw);ep=ca_distance(u,v)
 out=[]
 for name,state,z,scope in old.candidates(sc)[-8:]:
  w=pair_vectors(old.field_state(a,state,sc['mirror'],z),pairs,raw)
  x=pair_vectors(old.field_state(b,state,sc['mirror'],z),pairs,raw)
  arrays=(ca_distance(u,w),ca_distance(v,x),ep,ca_distance(w,x))
  for node in (0,1,2):
   for delta in (2,5,10,'JOINT'):
    ids=[i for i,vv in enumerate(pairs) if delta=='JOINT' or vv[1]==delta]
    def rms(ar):
     ar=ar[:,:,ids] if node==2 else ar[:,node,ids]
     return float(np.sqrt(np.mean(ar**2)))
    s0,s1,et,ev=map(rms,arrays);bound=et+ev
    out.append(dict(scene_id=sc['scene_id'],candidate_id=name,scope=scope,node='DUAL_NONCOHERENT' if node==2 else node,delta_hz=delta,separation_20001=s0,separation_40001=s1,absolute_separation_change=abs(s1-s0),truth_CA_grid_error=et,candidate_CA_grid_error=ev,triangle_bound=bound,numerical_bound_over_min_separation=bound/max(min(s0,s1),1e-300),separation_intervals_exclude_zero=bool(min(s0,s1)>bound),rank_20001=0,rank_40001=0,rank_changed=False))
 for node in (0,1,'DUAL_NONCOHERENT'):
  for delta in (2,5,10,'JOINT'):
   rs=[r for r in out if r['node']==node and r['delta_hz']==delta]
   for key,rank in [('separation_20001','rank_20001'),('separation_40001','rank_40001')]:
    for k,r in enumerate(sorted(rs,key=lambda r:(r[key],r['candidate_id']))):r[rank]=k+1
   for r in rs:r['rank_changed']=r['rank_20001']!=r['rank_40001']
 return out
def shared_noise(p,sc,z,raw):
 # Linear propagation of the SAME proper raw-frequency errors; no Fisher/information matrix.
 ix={f:i for i,f in enumerate(raw)};freqs=[150.,152.,154.,156.,158.];pr=[(150.,152.),(152.,154.),(156.,158.)];out=[]
 for ti in range(3):
  for node in range(2):
   jac=[];expected=[]
   for lo,hi in pr:
    aa=p[ti,node,:,ix[hi]]*p[ti,node,:,ix[lo]].conj();norm=np.linalg.norm(aa);u=aa/norm
    j=np.zeros((128,80))
    for fi,ff in enumerate(freqs):
     for el in range(8):
      for im in (0,1):
       da=np.zeros(8,dtype=complex)
       if ff==hi:da[el]=aa[el]*(1j if im else 1)
       if ff==lo:da[el]=aa[el]*(-1j if im else 1)
       du=(da-u*np.vdot(u,da).real)/norm
       dc=np.outer(du,u.conj())+np.outer(u,du.conj())
       j[:,fi*16+el*2+im]=np.r_[dc.real.ravel(),dc.imag.ravel()]
    jac.append(j);expected.append(.01**2*4*(1-np.sum(abs(u)**4)))
   marginal=[float(np.sum(j*j)*.01**2/2) for j in jac]
   cross=jac[0]@jac[1].T*.01**2/2
   disjoint=jac[0]@jac[2].T*.01**2/2
   out.append(dict(scene_id=sc['scene_id'],source_depth_m=z,time_s=int(old.TIMES[ti]),node=node,shared_raw_frequency_hz=152,shared_cross_covariance_frobenius=float(np.linalg.norm(cross)),disjoint_cross_covariance_frobenius=float(np.linalg.norm(disjoint)),marginal_formula_max_error=float(np.max(abs(np.array(marginal)-expected))),noise='PROPER_COMPLEX_RELATIVE_RAW_1PCT',fivepct_covariance_multiplier=25,pair_independence_assumed=False))
 return out
def third_mesh(a,b,raw):
 third=OUT/'third_mesh';third.mkdir(exist_ok=True);log=[];out=[]
 distances=np.array([50000.,50007.,55000.,55007.,60000.,60007.])
 for f in (150.,200.,244.,250.):
  stem=f'fine_f{int(f)}_n80001';path=third/stem
  Path(str(path)+'.env').write_text(old.modal_env(f,0,80001),encoding='utf-8')
  begin=time.monotonic()
  try:
   result=subprocess.run([str(old.AT/'kraken.exe'),stem],cwd=third,capture_output=True,timeout=90)
   Path(str(path)+'.stdout.txt').write_bytes(result.stdout+result.stderr)
   if result.returncode!=0:raise RuntimeError('solver exit '+str(result.returncode))
   c=old.parse_mod(str(path)+'.mod');status='GENERATED'
  except Exception as e:c=None;status=type(e).__name__+':'+str(e)
  log.append(dict(frequency_hz=f,mesh=80001,elapsed_s=time.monotonic()-begin,status=status,timeout_s=90))
  i=raw.index(f)
  for z in (180.,190.,199.,200.,201.,210.,220.):
   p=old.pressure(a[i],distances,z);q=old.pressure(b[i],distances,z)
   for x,y,tag in ((p,q,'20001_TO_40001'),(q,old.pressure(c,distances,z) if c else None,'40001_TO_80001'),(p,old.pressure(c,distances,z) if c else None,'20001_TO_80001')):
    if y is None:error='';rawerror=''
    else:
     cc=np.vdot(x,y)/np.vdot(x,x);error=float(np.linalg.norm(y-cc*x)/np.linalg.norm(y));rawerror=float(np.linalg.norm(y-x)/np.linalg.norm(y))
    out.append(dict(frequency_hz=f,source_depth_m=z,comparison=tag,raw_relative=rawerror,common_removed_relative=error,modes_20001=len(a[i][2]),modes_40001=len(b[i][2]),modes_80001=len(c[2]) if c else '',status=status,homogeneous_80001_CA='NOT_COMPUTABLE_NO_REGISTERED_PAIR_PARTNERS',mesh80001_is_truth=False))
 write(OUT/'THIRD_MESH_CALLS.csv',log);write(OUT/'THIRD_MESH_CONVERGENCE.csv',out)
 return log,out
def field_audit():
 out=[]
 for f in (150,200,250):
  for z in (180,200,220):
   stem=MD/f'field_f{f}_z{z}';buf=Path(str(stem)+'.shd').read_bytes()
   rec=4*struct.unpack('<i',buf[:4])[0];head=np.frombuffer(buf[2*rec:3*rec],'<i4')
   actual=np.frombuffer(buf[9*rec:9*rec+12],'<f4').astype(float);measured=np.frombuffer(buf[10*rec:10*rec+24],'<c8').astype(complex)
   mod=old.parse_mod(MD/f'f{f*2:04d}_s+0_n20001.mod')
   def residual(p):
    c=np.vdot(p,measured)/np.vdot(p,p);return float(np.linalg.norm(measured-c*p)/np.linalg.norm(measured))
   exact=old.pressure(mod,actual,z);desired=np.array([54993.,55000.,55007.]);coord=float(np.max(abs(actual-desired)))
   depths,phi,k=mod;weight=phi[np.flatnonzero(depths==z)[0]]*phi[np.flatnonzero(depths==200)[0]]
   r32=actual.astype('f4');k32=k.astype('c8')
   term=np.sqrt(np.complex64(2*np.pi)/(k32[:,None]*r32))*np.exp(-np.complex64(1j)*k32[:,None]*r32-np.complex64(1j*np.pi/4))
   single=np.sum(weight.astype('c8')[:,None]*term,axis=0,dtype='c8')
   out.append(dict(frequency_hz=f,source_depth_m=z,nsd=int(head[4]),nrd=int(head[5]),nrr=int(head[6]),range_coordinate_max_error_m=coord,output_complex_precision='FLOAT32_REAL_AND_IMAG',double_reconstruction_residual=residual(exact),single_emulation_residual=residual(single),desired_coordinate_residual=residual(old.pressure(mod,desired,z)),attribution='FIELD_OUTPUT_PRECISION_PLAUSIBLE_NOT_BINARY_INTERNAL_PROOF',same_modal_input=True,independent_mesh_convergence_certificate=False))
 write(OUT/'FIELD_FORMAT_AND_PRECISION_AUDIT.csv',out)
 return out
def execute():
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('single execution only')
 freeze=read(OUT/'FORWARD_DIAGNOSTIC_FREEZE.json')
 assert git('log','-1','--pretty=%s')=='R4 E2: freeze forward fidelity attribution diagnostic'
 assert git('ls-remote','origin','refs/heads/main').split()[0]==git('rev-parse','HEAD')
 for p,r in freeze['bindings'].items():assert sha(p,r['kind']=='RAW')==r['sha256'],p
 dump(OUT/'EXECUTION_STARTED.json',dict(design_SHA=git('rev-parse','HEAD'),new_Monte_Carlo=0,new_audio=0,started_epoch=time.time()))
 start=time.monotonic();freq=read(PRIOR/'E2_REQUIRED_FREQUENCIES.json');raw=freq['raw_frequencies_hz'];pairs=[(r['center_hz'],r['delta_hz'],r['low_hz'],r['high_hz']) for r in freq['frequency_pairs']]
 artifact_format_audit();aggregate_history();a,b=modes(raw);fa=field_audit()
 calls,third=third_mesh(a,b,raw)
 decom=[];caout=[];noise=[];candidates=[];shared=[]
 for sc in read(PRIOR/'E2_SCENE_FREEZE.json')['scenes']:
  for z in (180.,200.,220.):
   p=old.field_state(a,sc['state'],sc['mirror'],z);q=old.field_state(b,sc['state'],sc['mirror'],z)
   np.savez_compressed(OUT/f"PRESSURE_{sc['scene_id']}_z{int(z)}.npz",coarse=p,fine=q)
   shared.extend(shared_noise(p,sc,z,raw));decom.extend(decomposition(p,q,sc,z,raw));r,n=ca_rows(p,q,sc,z,pairs,raw);caout.extend(r);noise.extend(n)
   if z==200:candidates.extend(candidate_rows(p,q,a,b,sc,pairs,raw))
  print('SCENE '+sc['scene_id'],flush=True)
 write(OUT/'COMPLEX_FIELD_ERROR_DECOMPOSITION.csv',decom);write(OUT/'CA_MESH_STABILITY.csv',caout);write(OUT/'CA_NUMERICAL_VS_NOISE.csv',noise);write(OUT/'CA_CANDIDATE_STABILITY.csv',candidates)
 write(OUT/'SHARED_FREQUENCY_NOISE_COVARIANCE.csv',shared)
 okthird=all(r['status']=='GENERATED' for r in calls)
 c01={(r['frequency_hz'],r['source_depth_m']):r['common_removed_relative'] for r in third if r['comparison']=='20001_TO_40001'}
 c12=[r['common_removed_relative']/max(c01[(r['frequency_hz'],r['source_depth_m'])],1e-300) for r in third if r['comparison']=='40001_TO_80001' and r['status']=='GENERATED']
 stable_noise=all(r['numerical_over_1pct_noise']<1 for r in noise)
 stable_candidate=all(r['separation_intervals_exclude_zero'] for r in candidates)
 meshtrend=okthird and all(r<1 for r in c12)
 # Relational diagnostic categories only. No new CA scientific admission threshold.
 if not stable_noise or not stable_candidate:decision='E2_CA_AND_RAW_FORWARD_UNSTABLE'
 elif meshtrend and all(r['modes_20001']==r['modes_40001']==r['modes_80001'] for r in third) and all(r['matched'] for r in hybrids):decision='E2_FORWARD_DISCRETIZATION_LIMIT_IDENTIFIED'
 elif stable_noise and stable_candidate:decision='E2_CA_STABLE_BUT_RAW_FORWARD_UNCERTIFIED'
 else:decision='E2_FORWARD_DIAGNOSTIC_INCONCLUSIVE'
 result=dict(decision=decision,original_E2_G0_admission='FAIL_UNCHANGED',E2_scientific_information='NOT_EVALUATED',new_KRAKEN_calls=len(calls),new_MC=0,new_received_audio=0,R4_percent=0,elapsed_s=time.monotonic()-start,raw_worst=max(r['raw_relative'] for r in decom),differential_worst=max(r['common_removed_relative'] for r in decom),common_energy_median=float(np.median([r['common_energy_fraction'] for r in decom])),CA_worst=max(r['CA_Frobenius'] for r in caout),CA_over_1pct_noise_max=max(r['numerical_over_1pct_noise'] for r in noise),CA_over_5pct_noise_max=max(r['numerical_over_5pct_noise'] for r in noise),candidate_rank_changes=sum(r['rank_changed'] for r in candidates),candidate_count=len(candidates),candidate_unresolved=sum(not r['separation_intervals_exclude_zero'] for r in candidates),third_mesh_all_decreases=meshtrend,third_difference_ratio_max=max(c12) if c12 else None,homogeneous_third_CA='NOT_COMPUTABLE_MISSING_PAIR_PARTNERS',FIELD_root_cause='PRECISION_PLAUSIBLE_UNRESOLVED_INTERNALS',provider_root_cause='MODE_SOLVER_DISCRETIZATION_SUPPORTED' if meshtrend else 'UNRESOLVED',next='STOP',lead_audit='PENDING')
 dump(OUT/'FORWARD_FIDELITY_DECISION.json',result);print(json.dumps(result,indent=2),flush=True)
hybrids=[]
if __name__=='__main__':execute()
