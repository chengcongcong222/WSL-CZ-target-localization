"""Frozen E2-G0 physical mechanism screen; no RNG, audio, optimizer or performance trial."""
from pathlib import Path
import json, csv, hashlib, subprocess, time, os, sys
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy import linalg
OUT=Path('results/R4_E2_G0_HLA_DIFFERENCE_INFORMATION')
MOD=OUT/'modal'
AT=Path('tools/acoustics_toolbox/atWin10/at/bin').resolve()
TEMPLATE=Path('results/R3_C2_Yang_SA_depth/R3_C2_1/_kraken_zgrid/zgrid_f201.env')
PARENT='43fa1f5726cd7ecc97c632a7f652c10f83df9eff'
TIMES=np.array([0.,600.,1200.])
SCALE=np.array([50000.,50000.,2.,2.,200.])
STEP=np.array([.05,.05,.00005,.00005,1.])
DEPTHS=np.array([180.,190.,195.,198.,199.,199.5,200.,200.5,201.,202.,205.,210.,220.])
def dump(p,x):
 Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def sha(p,raw=False):
 b=Path(p).read_bytes()
 return hashlib.sha256(b if raw else b.replace(b'\r\n',b'\n')).hexdigest()
def write(p,rows):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def command(*args):
 return subprocess.check_output(list(args),text=True,encoding='utf-8').strip()
def git(*a):return command('git',*a)
def frequencies():
 center=np.arange(150.,250.01,.5)
 pairs=[(float(f),d,float(f-d/2),float(f+d/2)) for d in (2,5,10) for f in center if 150<=f-d/2 and f+d/2<=250]
 raw=sorted(set(x for p in pairs for x in p[2:]))
 return center,pairs,np.array(raw)
def geometry(state,mirror,t=TIMES):
 post=np.maximum(t-600.,0.)
 main=np.column_stack([2*np.minimum(t,600.)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
 nodes=np.stack([main,main+[0.,mirror*5000.]],axis=1)
 # Fixed global x axis. Translation/turn retained; array attitude is an explicit new design assumption.
 offsets=np.array([-7.,-5.,-3.,-1.,1.,3.,5.,7.])
 pos=nodes[:,:,None,:]+np.stack([offsets,np.zeros(8)],axis=-1)[None,None,:,:]
 target=np.asarray(state[:2])+t[:,None]*np.asarray(state[2:4])
 return np.linalg.norm(target[:,None,None,:]-pos,axis=-1),nodes
def modal_env(freq,ssp,mesh):
 lines=TEMPLATE.read_text(encoding='utf-8').splitlines()
 out=["'E2_ESTD'",f'{freq:.6f}','1',"'CVWT'",f'{mesh} 0.0 5000.0']
 for line in lines[5:]:
  parts=line.split()
  try:
   if len(parts)>=6 and 0<=float(parts[0])<=5000:
    parts[1]=f'{float(parts[1])+ssp:.7f}';out.append(' '.join(parts));continue
  except ValueError:pass
  if line.strip().startswith("'R'"):break
 out+=["'R' 0.0",'1500.0 1800.0','65.0','1','200.0',str(len(DEPTHS))]
 out += [f'{z:.4f}' for z in DEPTHS]
 out += ['R']
 return '\n'.join(out)+'\n'
def parse_mod(path):
 b=Path(path).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0])
 h=np.frombuffer(b[84:108],'<i4');nt,nm=int(h[2]),int(h[3])
 z=np.frombuffer(b[4*rec:5*rec],'<f4')[:nt].astype(float)
 n=int(np.frombuffer(b[5*rec:5*rec+4],'<i4')[0])
 phi=np.stack([np.frombuffer(b[(7+j)*rec:(8+j)*rec],'<c8')[:nm] for j in range(n)],axis=1).astype(complex)
 k=np.frombuffer(b[(7+n)*rec:(7+n)*rec+n*8],'<c8').astype(complex)
 if nt!=nm or n<1 or not np.isfinite(phi).all() or not np.isfinite(k).all() or not np.all(k.real>0):raise ValueError('invalid modal binary')
 if not np.array_equal(z,DEPTHS):raise ValueError('exact requested output depths absent')
 return z,phi,k
def pressure(mod,r,z=200.,zr=200.):
 depths,phi,k=mod
 a=np.flatnonzero(depths==z);b=np.flatnonzero(depths==zr)
 if len(a)!=1 or len(b)!=1:raise ValueError('no exact depth; interpolation prohibited')
 weights=phi[a[0]]*phi[b[0]]
 rr=np.asarray(r).ravel()
 # exp(+i omega t), outgoing exp(-i k r). Complex k retained.
 v=np.sqrt(2*np.pi/(k[:,None]*rr[None,:]))*np.exp(-1j*k[:,None]*rr[None,:]-1j*np.pi/4)
 return (weights@v).reshape(np.shape(r))
def ca(p,pairs,raw):
 ix={float(f):i for i,f in enumerate(raw)}
 a=np.stack([p[...,ix[hi]]*p[...,ix[lo]].conj() for _,_,lo,hi in pairs],axis=-2)
 return a[..., :, :,None]*a[..., :,None,:].conj()/np.sum(abs(a)**2,axis=-1)[..., :,None,None]
def field_state(mods,state,mirror,z=200.,zr=200.):
 r,_=geometry(state,mirror)
 return np.stack([pressure(m,r,z,zr) for m in mods],axis=-1)
def plane_controls(pairs,raw):
 rows=[];offset=np.arange(-7.,8.,2.)
 for d in (2,5,10):
  pp=[p for p in pairs if p[1]==d]
  beta=.173
  spatial=np.exp(2j*np.pi*raw[:,None]*np.cos(beta)*offset[None,:]/1500.)
  def field(r,z,spectrum):
   common=spectrum*np.exp(-2j*np.pi*raw*r/1500.)*(1+.001*z)/np.sqrt(r)
   return common[:,None]*spatial
  p=field(55000.,200.,np.ones(len(raw))).T
  ref=ca(p,pp,raw)
  err=0.
  for r,z in ((50000.,180.),(60000.,220.),(55321.,195.)):
   src=np.exp(.7*(raw-200)/100+1j*raw*.081)
   err=max(err,float(np.max(abs(ca(field(r,z,src).T,pp,raw)-ref))))
  rows.append({'control':f'plane_common_source_range_depth_delta{d}','error':err,'limit':1e-10,'PASS':err<1e-10})
 # Single element C_A=1; arbitrary per-element/frequency response exactly absorbs changes.
 p=np.exp(.01j*raw)[None,:]
 rows.append({'control':'single_element_source_unknown','error':float(np.max(abs(ca(p,pairs,raw)-1))),'limit':1e-12,'PASS':bool(np.allclose(ca(p,pairs,raw),1))})
 p1=np.exp(.02j*raw)[None,:]*np.exp(.07*np.arange(8))[:,None]
 p2=np.exp(-.09j*raw)[None,:]*np.exp(.03*np.arange(8))[:,None]
 gain=p1/p2
 rows.append({'control':'C2_arbitrary_response_absorption','error':float(np.max(abs(ca(p2*gain,pairs,raw)-ca(p1,pairs,raw)))),'limit':1e-12,'PASS':bool(np.allclose(ca(p2*gain,pairs,raw),ca(p1,pairs,raw)))})
 dump(OUT/'PLANE_CONTROL_RECORDS.json',rows)
 return all(r['PASS'] for r in rows)
def make_modes(raw,budget):
 MOD.mkdir(exist_ok=True);start=time.monotonic();rows=[];cache={}
 for ssp,mesh in ((0,20001),(0,40001),(-1,20001),(1,20001)):
  values=[]
  for f in raw:
   stem=f'f{int(round(f*2)):04d}_s{ssp:+d}_n{mesh}'
   path=MOD/stem;Path(str(path)+'.env').write_text(modal_env(f,ssp,mesh),encoding='utf-8')
   remain=budget-(time.monotonic()-start)
   if remain<=0:raise TimeoutError('frozen total modal-generation wall budget exhausted')
   began=time.monotonic()
   result=subprocess.run([str(AT/'kraken.exe'),stem],cwd=MOD,capture_output=True,timeout=min(90,remain))
   Path(str(path)+'.stdout.txt').write_bytes(result.stdout+result.stderr)
   if result.returncode!=0 or not Path(str(path)+'.mod').exists():raise RuntimeError('KRAKEN failed '+stem)
   val=parse_mod(str(path)+'.mod');values.append(val)
   rows.append(dict(frequency_hz=float(f),ssp_offset_mps=ssp,mesh=mesh,modes=len(val[2]),elapsed_s=time.monotonic()-began,mod_sha256=sha(str(path)+'.mod',True)))
   write(OUT/'E2_MODAL_GENERATION.csv',rows)
   if len(rows)%20==0:print(f'MODAL {len(rows)}/{len(raw)*4}; elapsed {time.monotonic()-start:.1f}s',flush=True)
  cache[(ssp,mesh)]=values
 return cache
def verify_physics(cache,raw):
 rows=[];worst=0.
 # Mesh convergence at every raw frequency, all registered source depths and spatial endpoint ranges.
 r=np.array([50000.,50007.,55000.,55007.,60000.,60007.])
 for i,f in enumerate(raw):
  a=cache[(0,20001)][i];b=cache[(0,40001)][i]
  for z in (180.,190.,199.,200.,201.,210.,220.):
   p=pressure(a,r,z);q=pressure(b,r,z)
   c=np.vdot(p,q)/np.vdot(p,p)
   err=float(np.linalg.norm(q-c*p)/np.linalg.norm(q))
   rows.append(dict(kind='mesh40001_vs20001',frequency_hz=float(f),source_depth_m=z,error=err,limit=.002,PASS=err<.002));worst=max(worst,err)
 write(OUT/'E2_FORWARD_CONVERGENCE.csv',rows)
 # Independent FIELD execution: nine comparisons, local aperture span, and separate 50--60 km span.
 for f in (150.,200.,250.):
  i=int(np.flatnonzero(raw==f)[0]);base=MOD/f'f{int(f*2):04d}_s+0_n20001'
  for z in (180.,200.,220.):
   for center in (55000.,):
    stem=f'field_f{int(f)}_z{int(z)}'
    dest=MOD/stem
    Path(str(dest)+'.mod').write_bytes(Path(str(base)+'.mod').read_bytes())
    Path(str(dest)+'.env').write_text(modal_env(f,0,20001),encoding='utf-8')
    # The FIELD requested source and receiver depths are supplied by flp; retained exact modal output.
    text="\n".join(["'E2_FIELD'","'RA'","9999","1","0.0","3",f'{(center-7)/1000:.8f} {(center+7)/1000:.8f} /',"1",f'{z:.3f}',"1","200.0","1","0.0 /"])+"\n"
    Path(str(dest)+'.flp').write_text(text,encoding='utf-8')
    res=subprocess.run([str(AT/'field.exe'),stem],cwd=MOD,capture_output=True,timeout=60)
    Path(str(dest)+'.stdout.txt').write_bytes(res.stdout+res.stderr)
    if res.returncode!=0:raise RuntimeError('FIELD failed')
    b=Path(str(dest)+'.shd').read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0])
    h=np.frombuffer(b[2*rec:3*rec],'<i4');ns,nrd,nr=map(int,h[4:7])
    if (ns,nrd,nr)!=(1,1,3):raise ValueError('FIELD layout not admitted')
    rr=np.frombuffer(b[9*rec:10*rec],'<f4')[:nr].astype(float)
    observed=np.frombuffer(b[10*rec:10*rec+nr*8],'<c8').astype(complex)
    calculated=pressure(cache[(0,20001)][i],rr,z)
    c=np.vdot(calculated,observed)/np.vdot(calculated,calculated)
    err=float(np.linalg.norm(observed-c*calculated)/np.linalg.norm(observed))
    rows.append(dict(kind='independent_FIELD_aperture',frequency_hz=f,source_depth_m=z,error=err,limit=.002,PASS=err<.002));worst=max(worst,err)
 write(OUT/'E2_FORWARD_CONVERGENCE.csv',rows)
 return all(r['PASS'] for r in rows),worst
def spatial_basis():
 # Orthonormal complement of a common complex source for eight elements.
 return linalg.null_space(np.ones((1,8))).T
def freq_basis(pairs,raw,selected):
 ix={float(f):i for i,f in enumerate(raw)}
 p=[p for p in pairs if selected=='JOINT' or p[1]==selected]
 plus=np.zeros((len(p),len(raw)));minus=plus.copy()
 for j,(_,d,lo,hi) in enumerate(p):
  plus[j,ix[lo]]=1;plus[j,ix[hi]]=1
  minus[j,ix[lo]]=-1;minus[j,ix[hi]]=1
 def basis(a):
  _,s,v=linalg.svd(a,full_matrices=False);return v[s>1e-10*s[0]]
 return basis(plus),basis(minus),plus,minus
def projected(j,qp,qm,qs,calibration):
 # j[t,node,m,f,param] is complex relative-pressure derivative.
 a=np.einsum('sm,tnmfp->tnsfp',qs,j.real)
 b=np.einsum('sm,tnmfp->tnsfp',qs,j.imag)
 a=np.einsum('qf,tnsfp->tnsqp',qp,a)
 b=np.einsum('qf,tnsfp->tnsqp',qm,b)
 if calibration=='C1':
  # Gain g_m fixed over all times/frequencies. AP cancels its phase exactly.
  v=qp@np.ones(j.shape[3])
  aa=a.transpose(1,2,0,3,4).reshape(j.shape[1],7,-1,j.shape[-1])
  vv=np.tile(v,j.shape[0]);den=vv@vv
  if den>1e-18:aa-=vv[None,None,:,None]*np.einsum('d,nsdp->nsp',vv,aa)[:,:,None,:]/den
  a=aa
 elif calibration=='C2':
  # Snapshot-local arbitrary frequency/element gains: saturated absorption control, not a time-constant calibration claim.
  return np.zeros((1,j.shape[-1]))
 return np.concatenate([a.reshape(-1,j.shape[-1]),b.reshape(-1,j.shape[-1])],axis=0)
def raw_project(j,qs,calibration):
 a=np.einsum('sm,tnmfp->tnsfp',qs,j.real);b=np.einsum('sm,tnmfp->tnsfp',qs,j.imag)
 if calibration=='C1':
  # Fixed unknown complex gains across time and frequency at each node/element.
  a-=a.mean(axis=(0,3),keepdims=True);b-=b.mean(axis=(0,3),keepdims=True)
 return np.concatenate([a.reshape(-1,j.shape[-1]),b.reshape(-1,j.shape[-1])])
def profile_col(j,col):
 other=[i for i in range(j.shape[1]) if i!=col]
 q=j[:,other];u,s,_=linalg.svd(q,full_matrices=False)
 use=s>1e-10*s[0] if len(s) and s[0]>0 else np.zeros(len(s),bool)
 residual=j[:,col]-(u[:,use]@(u[:,use].T@j[:,col]) if use.any() else 0)
 return float(residual@residual)
def metrics(j,sigma):
 w=j*np.sqrt(2)/sigma;f=w.T@w
 eigen=np.linalg.eigvalsh(f)
 # initial radial direction transformed to a scene-local orthonormal horizontal frame by caller.
 ez=profile_col(w,4);er=profile_col(w,0)
 exact=profile_col(w[:,:4],0)
 _,_,vh=linalg.svd(w,full_matrices=False)
 return dict(range_information_depth_profiled=er/SCALE[0]**2,range_information_exact_depth=exact/SCALE[0]**2,depth_information=ez/SCALE[4]**2,scaled_weak_eigenvalue=float(eigen[0]),scaled_weak_direction=json.dumps(vh[-1].tolist()),target_rank=int(np.sum(eigen>max(eigen[-1],1)*1e-10)))
def baseline(scene):
 import r4_e1_g0 as e1
 t=np.arange(0.,1200.1,10.)
 _,_,_,_,bx,_,_=e1.geometry(scene,t)
 j=bx.reshape(-1,4)
 sig=np.deg2rad(.05)
 nav=np.sum(j[:,:2]**2,axis=1)*25**2
 cov=np.diag(sig**2+nav)
 common=np.ones(len(j));diff=np.tile([1.,-1.],len(t))
 cov+=np.deg2rad(.05)**2/3*np.outer(common,common)+np.deg2rad(.025)**2/3*np.outer(diff,diff)
 w=linalg.solve_triangular(linalg.cholesky(cov,lower=True),j,lower=True)
 return w.T@w
def candidates(scene):
 # Nominal observation-side bearing triangulation; no optimizer and no truth-selected favorable range.
 generating=np.array(scene['state']);_,nodes=geometry(generating,scene['mirror'])
 observed=generating[:2]+TIMES[:,None]*generating[2:4]
 d=observed[:,None,:]-nodes
 bearing=np.arctan2(d[:,:,1],d[:,:,0])
 uv=np.stack([np.cos(bearing),np.sin(bearing)],axis=-1)
 cross=lambda a,b:a[...,0]*b[...,1]-a[...,1]*b[...,0]
 radial=cross(nodes[:,1]-nodes[:,0],uv[:,1])/cross(uv[:,0],uv[:,1])
 positions=nodes[:,0]+radial[:,None]*uv[:,0]
 state=np.r_[positions[0],(positions[-1]-positions[0])/1200.]
 b=float(bearing[0,0]);u=np.array([np.cos(b),np.sin(b)])
 rows=[]
 for r in np.arange(50000.,60000.1,1000.):
  for z in (180.,190.,200.,210.,220.):
   # Velocity from exact nominal dual-bearing triangulation is a diagnostic control; not an estimator performance result.
   x=state.copy();x[:2]=r*u
   rows.append((f'R{int(r)}_Z{int(z)}',x,float(z),'ABSOLUTE_FULL_REGISTERED_GRID'))
 for dr,dz in ((-250.,0.),(250.,0.),(0.,-10.),(0.,10.),(-250.,-10.),(-250.,10.),(250.,-10.),(250.,10.)):
  # Local perturbations centered on the observation-side nominal bearing triangulation are finite sensitivity controls, not a coverage certificate.
  r=np.linalg.norm(state[:2])+dr;x=state.copy();x[:2]=r*u
  rows.append((f'LOCAL_{dr:+g}_{dz:+g}',x,200.+dz,'NOMINAL_BEARING_CENTERED_DIAGNOSTIC_NOT_COVERAGE'))
 return rows
def finite_scores(p,q,pairs,raw,cal):
 c=ca(p,pairs,raw);d=ca(q,pairs,raw)
 # exact normalized spatial distance. Calibration is profiled in log-amplitude ratios, phases cancel for fixed gains.
 if cal=='C1':
  la=np.log(abs(q))-np.log(abs(p))
  gain=np.exp(-la.mean(axis=(0,3),keepdims=True))
  q=q*gain
  d=ca(q,pairs,raw)
 return float(np.sqrt(np.mean(np.sum(abs(d-c)**2,axis=(-1,-2)))))
def noise_ca_rms(p,pairs,raw,sigma):
 ix={float(f):i for i,f in enumerate(raw)}
 vals=[]
 for _,_,lo,hi in pairs:
  a=p[...,ix[hi]]*p[...,ix[lo]].conj()
  w=abs(a)**2/np.sum(abs(a)**2,axis=-1,keepdims=True)
  vals.append(4*sigma**2*(1-np.sum(w*w,axis=-1)))
 return float(np.sqrt(np.mean(vals)))
def covariance_controls(pairs,raw):
 # Independent dense covariance reconstruction of a shared-pressure graph, all eight elements.
 tiny=np.arange(150.,160.1,.5)
 pp=[(float(f),d,float(f-d/2),float(f+d/2)) for d in (2,5) for f in np.arange(153.,157.1,.5)]
 qp,qm,bp,bm=freq_basis(pp,tiny,'JOINT');qs=spatial_basis()
 d=np.zeros((7,8))
 for k in range(7):d[k,k+1]=1;d[k,0]=-1
 j=np.empty((1,1,8,len(tiny),5),complex)
 for m in range(8):
  for f in range(len(tiny)):
   for c in range(5):j[0,0,m,f,c]=np.sin((m+1)*(f+2)*(c+1)*.071)+1j*np.cos((m+2)*(f+1)*(c+1)*.093)
 a=np.kron(d,bp);b=np.kron(d,bm)
 jr=j[0,0].real.reshape(-1,5);ji=j[0,0].imag.reshape(-1,5)
 ya=a@jr;yb=b@ji
 fa=ya.T@linalg.pinvh(a@a.T,rtol=1e-10)@ya
 fb=yb.T@linalg.pinvh(b@b.T,rtol=1e-10)@yb
 z=projected(j,qp,qm,qs,'C0');f=z.T@z
 error=float(np.linalg.norm(f-fa-fb)/np.linalg.norm(f))
 offdiag=float(np.max(abs((bp@bp.T)-np.diag(np.diag(bp@bp.T)))))
 result=dict(dense_vs_factored_relative_error=error,shared_frequency_offdiagonal=offdiag,sigma_scaling_5percent_information_ratio=.04,normalized_CA_chart='local amplitude log ratio and phase ratio; invertible rank-one CA chart',raw_relative_noise='proper complex Gaussian; each real component sigma^2/2; covariance fixed at generating field; no covariance score',PASS=error<1e-9 and offdiag>0)
 dump(OUT/'COVARIANCE_CONTROL.json',result)
 if not result['PASS']:raise ValueError('joint autoproduct covariance control failed')
def evaluate(cache,pairs,raw,forward_error):
 scenes=read('results/R4_E1_G0_FREQUENCY_INFORMATION/E1_G0_SCENE_BINDINGS.json')['scenes']
 covariance_controls(pairs,raw)
 qs=spatial_basis();bases={s:freq_basis(pairs,raw,s) for s in (2,5,10,'JOINT')}
 info=[];nuis=[];cand=[];cal=[];valid=[];saved={};started=time.monotonic()
 for sc in scenes:
  if time.monotonic()-started>1800:raise TimeoutError('frozen information/candidate wall budget exhausted')
  x=np.r_[sc['state'],200.];p=field_state(cache[(0,20001)],x,sc['mirror'])
  fine=field_state(cache[(0,40001)],x,sc['mirror'])
  jac=[];half=[];finejac=[]
  for col,h in enumerate(STEP):
   def derivative(mods,step):
    up=x.copy();dn=x.copy();up[col]+=step;dn[col]-=step
    a=field_state(mods,up,sc['mirror'],up[4]);b=field_state(mods,dn,sc['mirror'],dn[4])
    base=p if mods is cache[(0,20001)] else fine
    return (a-b)/(2*step)/base*SCALE[col]
   jac.append(derivative(cache[(0,20001)],h));half.append(derivative(cache[(0,20001)],h/2));finejac.append(derivative(cache[(0,40001)],h))
  j=np.stack(jac,axis=-1);h=np.stack(half,axis=-1);fj=np.stack(finejac,axis=-1)
  angle=np.arctan2(x[1],x[0]);rot=np.eye(5);rot[:2,:2]=[[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
  j=j@rot;h=h@rot;fj=fj@rot
  np.savez_compressed(OUT/f"INPUT_{sc['scene_id']}.npz",relative_pressure_J=j,pressure=p)
  # Errors must be assessed after source and calibration projection, not on common phase.
  for nnode in (1,2):
   for sel in (2,5,10,'JOINT'):
    qp,qm,_,_=bases[sel]
    for c in ('C0','C1','C2'):
     z=projected(j[:,:nnode],qp,qm,qs,c)
     zh=projected(h[:,:nnode],qp,qm,qs,c)
     zf=projected(fj[:,:nnode],qp,qm,qs,c)
     den=max(np.linalg.norm(zh),1e-30)
     fd=float(np.linalg.norm(z-zh)/den) if c!='C2' else 0.
     me=float(np.linalg.norm(z-zf)/den) if c!='C2' else 0.
     er=profile_col(z,0);erh=profile_col(zh,0);erf=profile_col(zf,0)
     reld=max(abs(er-erh),abs(er-erf))/max(er,1e-30) if c!='C2' else 0.
     valid.append(dict(scene_id=sc['scene_id'],nodes=nnode,delta_hz=sel,calibration=c,fd_relative=fd,mesh_relative=me,PASS=fd<.02 and me<.02 and (c=='C2' or reld<.1)))
     for sigma in (.01,.05):
      m=metrics(z,sigma);m.update(scene_id=sc['scene_id'],representation='P1' if nnode==1 else 'P2',condition=c,delta_hz=sel,sigma_relative=sigma,source_spectrum='UNKNOWN_PROFILED_BY_NORMALIZATION',navigation='EXACT_ARRAY_GEOMETRY_MECHANISM_CONDITION',physical_error_bound=forward_error,status='EVALUATED_LOCAL_ONLY')
      info.append(m)
     if sel=='JOINT':saved[(sc['scene_id'],nnode,c)]=z
  for c in ('C0','C1'):
   z=raw_project(j,qs,c)
   for sigma in (.01,.05):
    m=metrics(z,sigma);m.update(scene_id=sc['scene_id'],representation='P0',condition=c,delta_hz='RAW',sigma_relative=sigma,source_spectrum='UNKNOWN_PER_NODE_FREQ_TIME',navigation='EXACT_ARRAY_GEOMETRY_MECHANISM_CONDITION',physical_error_bound=forward_error,status='EVALUATED_LOCAL_ONLY');info.append(m)
   z2=saved[(sc['scene_id'],2,c)]
   mine=float(np.linalg.eigvalsh(z.T@z-z2.T@z2)[0])
   valid.append(dict(scene_id=sc['scene_id'],nodes=2,delta_hz='DPI',calibration=c,fd_relative=0.,mesh_relative=0.,PASS=mine>=-1e-7*max(float(np.linalg.norm(z.T@z)),1.)))
  b0=baseline(sc);dump(OUT/f"B0_{sc['scene_id']}.json",dict(FIM=b0.tolist(),covariance='original bearing random RMS+25m nav+uniform common/differential moments',fused_with_AP=False,reason='AP is exact-navigation physical control; no shared-navigation-independent fusion claim'))
  for condition in ('SSP_MINUS1','SSP_PLUS1','ZR_MINUS5','ZR_PLUS5'):
   mods=cache[(-1,20001)] if condition=='SSP_MINUS1' else cache[(1,20001)] if condition=='SSP_PLUS1' else cache[(0,20001)]
   zr=195. if condition=='ZR_MINUS5' else 205. if condition=='ZR_PLUS5' else 200.
   q=field_state(mods,x,sc['mirror'],200.,zr)
   for c in ('C0','C1'):
    distance=finite_scores(p,q,pairs,raw,c)
    nuis.append(dict(scene_id=sc['scene_id'],condition=condition,calibration=c,normalized_CA_rms=distance,status='SINGLE_FACTOR_MISMATCH_DISTANCE_NOT_PERFORMANCE'))
  for name,state,z,scope in candidates(sc):
   r=np.linalg.norm(state[:2])
   if not 50000<=r<=60000:
    cand.append(dict(scene_id=sc['scene_id'],candidate_id=name,range_m=r,depth_m=z,scope=scope,C0_CA_rms='',C1_CA_rms='',noise_rms_1pct='',noise_rms_5pct='',status='OUTSIDE_REGISTERED_INITIAL_RANGE'));continue
   q=field_state(cache[(0,20001)],state,sc['mirror'],z)
   cand.append(dict(scene_id=sc['scene_id'],candidate_id=name,range_m=r,depth_m=z,scope=scope,C0_CA_rms=finite_scores(p,q,pairs,raw,'C0'),C1_CA_rms=finite_scores(p,q,pairs,raw,'C1'),noise_rms_1pct=noise_ca_rms(p,pairs,raw,.01),noise_rms_5pct=noise_ca_rms(p,pairs,raw,.05),status='FINITE_DIAGNOSTIC_NOT_GLOBAL_EXCLUSION'))
  c0=saved[(sc['scene_id'],2,'C0')];c1=saved[(sc['scene_id'],2,'C1')]
  cal.append(dict(scene_id=sc['scene_id'],C0_trace=float(np.sum(c0*c0)),C1_trace=float(np.sum(c1*c1)),trace_loss_fraction=float(1-np.sum(c1*c1)/max(np.sum(c0*c0),1e-30)),phase_gain='CANCELS_EXACTLY_IF_FREQUENCY_INDEPENDENT',amplitude_gain='UNKNOWN_FIXED_ALL_TIMES_FREQS',C2='ZERO_SNAPSHOT_FREE_RESPONSE_CONTROL'))
  print('SCENE '+sc['scene_id'],flush=True)
 write(OUT/'E2_FULL_SCENE_INFORMATION.csv',info);write(OUT/'E2_NUISANCE_PROFILE_RESULTS.csv',nuis)
 write(OUT/'E2_CANDIDATE_SEPARATION.csv',cand);write(OUT/'E2_CALIBRATION_SENSITIVITY.csv',cal)
 write(OUT/'E2_LOCAL_NUMERICAL_CHECKS.csv',valid)
 dump(OUT/'INFORMATION_SUMMARY.json',dict(checks_pass=all(r['PASS'] for r in valid),information_rows=len(info),candidate_rows=len(cand),all24_scenes=True,exact_navigation_assumption=True,bearing_fusion=False,not_global_uniqueness=True))
 # Navigation/association/timing are not validated by this ideal complex-field mechanism.
 if not all(r['PASS'] for r in valid):return 'E2_G0_PHYSICS_OR_COVARIANCE_INCOMPLETE'
 primary=[r for r in info if r['representation']=='P2' and r['condition']=='C1' and r['delta_hz']=='JOINT' and r['sigma_relative']==.01]
 if all(r['range_information_depth_profiled']>1e-14 and r['depth_information']>1e-10 for r in primary):return 'E2_G0_HLA_NONBEARING_INFORMATION_PRESENT'
 upper=[r for r in info if r['representation']=='P2' and r['condition']=='C0' and r['delta_hz']=='JOINT' and r['sigma_relative']==.01]
 if all(r['range_information_depth_profiled']>1e-14 for r in upper):return 'E2_G0_CALIBRATION_CONDITIONAL'
 return 'E2_G0_CURRENT_HLA_SPATIAL_AUTOPRODUCT_NOT_SUPPORTED'
def blocked_tables(reason):
 scenes=read('results/R4_E1_G0_FREQUENCY_INFORMATION/E1_G0_SCENE_BINDINGS.json')['scenes']
 for name in ('E2_FULL_SCENE_INFORMATION.csv','E2_NUISANCE_PROFILE_RESULTS.csv','E2_CANDIDATE_SEPARATION.csv','E2_CALIBRATION_SENSITIVITY.csv'):
  write(OUT/name,[dict(scene_id=s['scene_id'],status='NOT_EVALUATED',reason=reason) for s in scenes])
def execute():
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('one execution only; no restart')
 freeze=read(OUT/'E2_G0_DESIGN_FREEZE.json')
 if git('log','-1','--pretty=%s')!=freeze['commit_A_title']:raise RuntimeError('execute only immediately after design commit')
 if git('ls-remote','origin','refs/heads/main').split()[0]!=git('rev-parse','HEAD'):raise RuntimeError('design not remote')
 for p,h in freeze['bindings'].items():
  if sha(p,p.endswith(('.exe','.pdf','.mod'))) !=h:raise RuntimeError('frozen binding changed: '+p)
 dump(OUT/'EXECUTION_STARTED.json',dict(design_SHA=git('rev-parse','HEAD'),started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),new_MC=0,new_audio=0))
 _,pairs,raw=frequencies()
 decision='E2_G0_PHYSICS_OR_COVARIANCE_INCOMPLETE';reason='';error=None
 start=time.monotonic()
 try:
  if not plane_controls(pairs,raw):
   decision='E2_G0_IMPLEMENTATION_INVALID';reason='plane-wave null control failed';blocked_tables(reason)
  else:
   cache=make_modes(raw,freeze['budgets']['modal_wall_seconds'])
   good,error=verify_physics(cache,raw)
   dump(OUT/'PHYSICS_ADMISSION.json',dict(PASS=good,max_relative_error=error,threshold=.002))
   if not good:reason='complex-field mesh/FIELD admission not closed at 0.2% tolerance';blocked_tables(reason)
   else:decision=evaluate(cache,pairs,raw,error);reason='Registered physical-field screen; numerical admission required; no extracted performance.'
 except Exception as e:
  reason=type(e).__name__+': '+str(e);blocked_tables(reason)
 result=dict(decision=decision,reason=reason,elapsed_s=time.monotonic()-start,plane_null_PASS=all(x['PASS'] for x in read(OUT/'PLANE_CONTROL_RECORDS.json')) if (OUT/'PLANE_CONTROL_RECORDS.json').exists() else False,physics_error=error,new_Monte_Carlo=0,new_received_audio=0,R4_percent=0,next='STOP',automatic_other_stage=False,pending_independent_research_lead_audit=True)
 dump(OUT/'E2_G0_DECISION.json',result);print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':execute()
