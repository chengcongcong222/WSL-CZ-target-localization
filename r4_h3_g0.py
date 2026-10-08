"""Once-only finite vertical mechanism. No solver, RNG, audio or optimizer."""
from pathlib import Path
import os,time,json,csv
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy import linalg
import r4_e2_g0 as modal
import r4_e2_paired_recovery as paths
import r4_e2_forward_diagnostic as base
O=Path('results/R4_H3_G0_FINITE_VERTICAL_RESPONSE')
PARENT='4868353124c31c64d4283589339e9ae5bad72aa9'
RAW=paths.RAW;GRIDS=[80001,160001];TIMES=np.array([0.,600.,1200.])
LABELS=[180.,190.,200.,210.,220.]
SCALE=np.array([50000.,1.,200.,2.,1.]) # r,theta,z,v,psi
REF=.00010824612978073539
read,write,sha,git=base.read,base.write,base.sha,base.git
def native(x):
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,dict):return {k:native(v) for k,v in x.items()}
 if isinstance(x,(tuple,list)):return [native(v) for v in x]
 return x
def dump(p,x):base.dump(p,native(x))
def rows(p):return base.rows(p)
def geometry(st,resource):
 post=np.maximum(TIMES-600,0)
 node=np.stack([2*np.minimum(TIMES,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)],axis=-1)
 off=np.array([-7.,-5.,-3.,-1.,1.,3.,5.,7.])
 zz=np.full(8,200.)
 if resource in ['B1','B2','B3']:
  off=np.r_[off,np.full(4,-1.)];zz=np.r_[zz,[198.,199.,201.,202.] if resource!='B2' else np.full(4,200.)]
 target=np.array(st[:2])+TIMES[:,None]*np.array(st[2:4])
 delta=target[:,None,:]-node[:,None,:]-np.stack([off,np.zeros(len(off))],axis=-1)[None,:,:]
 r=np.sqrt(np.sum(delta*delta,axis=-1))
 return r,delta/r[...,None],zz
def field(mods,st,resource,z=200.,derivative=True):
 r,unit,zr=geometry(st,resource);ps=[];rad=[];dz=[];dh=[];bounds=[];gain_errors=[]
 for depths,phi,k in mods:
  ix={float(a):i for i,a in enumerate(depths)}
  assert all(v in ix for v in LABELS+[198.,199.,201.,202.,199.5,200.5])
  receiver=phi[[ix[v] for v in zr]]
  rr=r.reshape(-1)
  E=np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*k[:,None]*rr-1j*np.pi/4)
  # Same Cartesian distance array for pressure and every derivative.
  R=np.tile(receiver,(3,1)).T
  T=phi[ix[z],:,None]*R*E
  p=np.sum(T,axis=0).reshape(r.shape);ps.append(p)
  eps=np.finfo(float).eps;nop=64*len(k)+128;gamma=nop*eps/(1-nop*eps)
  bounds.append((2*gamma*np.sum(abs(T),axis=0)).reshape(r.shape))
  if derivative:
   rad.append(np.sum(T*(-1j*k[:,None]-.5/rr),axis=0).reshape(r.shape))
   for step,acc in [(1.,dz),(.5,dh)]:
    dphi=(phi[ix[z+step]]-phi[ix[z-step]])/(2*step)
    acc.append(np.sum(dphi[:,None]*R*E,axis=0).reshape(r.shape))
   src=phi[ix[z]];dsrc=(phi[ix[z+1]]-phi[ix[z-1]])/2
   assert np.all(src!=0),'nonregular modal source node: free multiplicative gain tangent not valid'
   g=dsrc/src
   err=np.sum(T*g[:,None],axis=0).reshape(r.shape)-dz[-1]
   gain_errors.append(float(np.linalg.norm(err)/max(np.linalg.norm(dz[-1]),1e-30)))
 p=np.stack(ps,axis=-1);bound=np.stack(bounds,axis=-1)
 if not derivative:return p,bound
 st=np.array(st);rn=np.linalg.norm(st[:2]);vn=np.linalg.norm(st[2:])
 geo=np.zeros((3,2,5));geo[:,:,0]=st[:2]/rn;geo[:,:,1]=[-st[1],st[0]]
 geo[:,:,3]=TIMES[:,None]*st[2:]/vn
 geo[:,:,4]=TIMES[:,None]*np.array([-st[3],st[2]])
 grad=np.stack(rad,axis=-1)[...,None]*unit[:,:,None,:]
 J=np.einsum('tnfd,tdc->tnfc',grad,geo)*SCALE
 H=J.copy();J[...,2]=np.stack(dz,axis=-1)*SCALE[2];H[...,2]=np.stack(dh,axis=-1)*SCALE[2]
 return p,J/p[...,None],H/p[...,None],bound,max(gain_errors)
def rvec(a):return np.r_[a.real.reshape(-1,*a.shape[2:]),a.imag.reshape(-1,*a.shape[2:])]
def remove(j,n):
 if n.size==0:return j.copy()
 norms=np.linalg.norm(n,axis=0);n=n[:,norms>1e-12]/norms[norms>1e-12]
 if not n.size:return j.copy()
 u,s,_=linalg.svd(n,full_matrices=False);q=u[:,s>s[0]*1e-10]
 return j-q@(q.T@j)
def system(p,j,noise,cal,raw_oracle=False):
 nt,N,F=p.shape;count=2*N*F
 y=[];gain=[];source=[]
 size=N if cal=='C1a' else N*F
 for t in range(nt):
  var=np.ones(count)/2 if noise=='RELATIVE' else np.tile((REF/abs(p[t]).ravel())**2/2,2)
  W=1/np.sqrt(var)
  yy=W[:,None]*rvec(j[t]);y.append(yy)
  src=linalg.block_diag(np.kron(np.ones((N,1)),np.eye(F)),np.kron(np.ones((N,1)),np.eye(F)))
  ns=np.zeros((count,nt*2*F));ns[:,t*2*F:(t+1)*2*F]=W[:,None]*src;source.append(ns)
  if cal=='C1a':ng=linalg.block_diag(np.kron(np.eye(N),np.ones((F,1))),np.kron(np.eye(N),np.ones((F,1))))
  else:ng=np.eye(count)
  gain.append(W[:,None]*ng)
 y=np.concatenate(y);src=np.concatenate(source)
 if raw_oracle:return y,np.zeros((len(y),0))
 n=src
 if cal in ['C1a','C1b']:n=np.column_stack([src,np.concatenate(gain)])
 return y,n
def projected(p,j,noise,cal,oracle=False):
 y,n=system(p,j,noise,cal,oracle)
 return np.zeros_like(y) if cal=='C2' else remove(y,n)
def metrics(z,sigma):
 a=z/sigma;u,s,vh=linalg.svd(a,full_matrices=False)
 rank=int(np.sum(s>s[0]*1e-10)) if s[0]>0 else 0
 null=vh[rank:];null_error=float(np.linalg.norm(a@null.T)/max(np.linalg.norm(a),1e-30)) if len(null) else 0.
 nz=float(np.linalg.norm(null[:,2])) if len(null) else 0.
 def prof(keep):
  rem=remove(a[:,2:3],a[:,keep]);v=float(np.sum(rem*rem))
  return v/SCALE[2]**2 if v>(max(s[0],1e-30)*1e-10)**2 else 0.
 iz=0. if nz>1e-8 else prof([0,1,3,4])
 return dict(rank=rank,structural_null_residual=null_error,condition_identifiable=float(s[0]/s[rank-1]) if rank else 'INF',depth_null_loading=nz,depth_information_horizontal_profiled=iz,depth_information_oracle_horizontal=prof([]),depth10_signal=10*np.sqrt(iz),depth_local_scale_m=1/np.sqrt(iz) if iz>0 else 'INF',singular_values=json.dumps(s.tolist()),weak_direction=json.dumps(vh[-1].tolist()),full_rank_required=False,information_trace=float(np.sum(a*a)))
def response_distance(p,q):
 u=p/np.linalg.norm(p,axis=1,keepdims=True);v=q/np.linalg.norm(q,axis=1,keepdims=True)
 a=np.einsum('tnf,tmf->tfnm',u,u.conj());b=np.einsum('tnf,tmf->tfnm',v,v.conj())
 return float(np.sqrt(np.mean(np.sum(abs(a-b)**2,axis=(-2,-1)))))
def response_noise(p,sigma,noise):
 if noise=='RELATIVE':
  w=abs(p)**2/np.sum(abs(p)**2,axis=1,keepdims=True)
  return float(np.sqrt(np.mean(2*sigma**2*np.sum(w*(1-w),axis=1))))
 return float(np.sqrt(np.mean(2*sigma**2*REF**2*(p.shape[1]-1)/np.sum(abs(p)**2,axis=1))))
CONTRAST_CACHE={}
def contrast(p,q,noise,cal):
 rel=np.log(abs(q/p))+1j*np.unwrap(np.angle(q/p),axis=1)
 J=np.zeros(p.shape+(5,),complex);J[...,0]=rel
 key=(id(p),noise,cal)
 if key not in CONTRAST_CACHE:
  y,n=system(p,J,noise,cal);v=np.linalg.norm(n,axis=0);nn=n[:,v>1e-12]/v[v>1e-12]
  u,s,_=linalg.svd(nn,full_matrices=False);CONTRAST_CACHE[key]=u[:,s>s[0]*1e-10]
 y,_=system(p,J,noise,cal);Q=CONTRAST_CACHE[key];d=y[:,0];d=d-Q@(Q.T@d)
 return float(np.linalg.norm(d)/.01)
def ck(out,name,error,limit,ok=None):
 out.append(dict(control=name,error=float(error),limit=float(limit),PASS=bool(error<=limit if ok is None else ok)))
def controls():
 out=[]
 p=np.ones((3,12,13),complex);j=np.zeros(p.shape+(5,),complex)
 j[...,2]=np.arange(13)[None,None,:]+1j
 for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
  z=projected(p,j,noise,'C0');raw,_=system(p,j,noise,'C0');e=np.linalg.norm(z)/max(np.linalg.norm(raw),1)
  ck(out,'separable_source_depth_'+noise,e,1e-10)
  for c in ['C1a','C1b']:
   jj=np.zeros_like(j);jj[...,2]=np.arange(12)[None,:,None]*(1+1j)
   zz=projected(p,jj,noise,c);raw,_=system(p,jj,noise,c);e=np.linalg.norm(zz)/max(np.linalg.norm(raw),1)
   ck(out,'fixed_gain_depth_absorption_'+noise+c,e,1e-10)
  yy,_=system(p,j,noise,'C0');ck(out,'C2_free_response_'+noise,np.linalg.norm(yy-np.eye(len(yy))@yy),0.)
 # Two independently parameterized group energies: derivative is exactly gain tangent.
 B=np.array([2.,3.]);dz=np.array([.7,-.2])
 ck(out,'group_gain_depth_absorption',np.linalg.norm(dz-np.diag(B)@(dz/B)),1e-14)
 return out
def execute():
 if (O/'EXECUTION_STARTED.json').exists():raise RuntimeError('registered once-only execution')
 f=read(O/'DESIGN_FREEZE.json')
 assert git('log','-1','--pretty=%s')=='R4 H3-G0: freeze finite vertical depth information mechanism'
 assert git('rev-parse','HEAD')==git('ls-remote','origin','refs/heads/main').split()[0]
 for n,r in f['bindings'].items():assert sha(n,r['kind']=='RAW')==r['sha256'],n
 dump(O/'EXECUTION_STARTED.json',dict(design_SHA=git('rev-parse','HEAD'),epoch=time.time(),new_KRAKEN=0,new_FIELD=0,new_MC=0))
 start=time.monotonic();ctrl=controls();
 from r4_h3_g0_audit import cold
 import r4_e2_nine_pair_pilot_audit as independent_binary
 info=[];fid=[];stable=[];finite=[];saved={}
 for sc in f['scenes']:
  for grid in GRIDS:
   mods=[modal.parse_mod(paths.path_for(ff,grid)) for ff in RAW]
   imod=[independent_binary.mod(paths.path_for(ff,grid)) for ff in RAW]
   for ff,m in zip(RAW,mods):
    fid.append(dict(scene_id=sc['scene_id'],mesh=grid,frequency_hz=ff,exact_source_depths=True,exact_receiver_depths=True,modes=len(m[2]),mode_identity_used_by_estimator=False))
   for res in ['B0','B1','B2']:
    p,j,h,bound,ge=field(mods,sc['state'],res)
    ck(ctrl,'free_modal_gain_'+sc['scene_id']+str(grid)+res,ge,1e-10)
    pp,jj,hh=cold(imod,sc['state'],res)
    ck(ctrl,'cold_pressure_term_bound_'+sc['scene_id']+str(grid)+res,float(np.max(abs(pp-p)/np.maximum(bound,1e-300))),1.)
    ck(ctrl,'cold_analytic_J_'+sc['scene_id']+str(grid)+res,float(np.linalg.norm(jj/pp[...,None]-j)/np.linalg.norm(j)),1e-9)
    ck(ctrl,'cold_source_derivative_'+sc['scene_id']+str(grid)+res,float(np.linalg.norm(hh/pp[...,None]-h)/np.linalg.norm(h)),1e-9)
    np.savez_compressed(O/f"INPUT_{sc['scene_id']}_{res}_n{grid}.npz",pressure=p,J=j,J_half=h,roundoff_bound=bound)
    saved[sc['scene_id'],grid,res]=(p,j,h)
  print('PRE_FORMULA',sc['scene_id'],flush=True)
 write(O/'SOURCE_RECEIVER_DEPTH_FIDELITY.csv',fid);write(O/'CONTROL_CHECKS.csv',ctrl)
 if not all(x['PASS'] for x in ctrl):
  dump(O/'H3_G0_DECISION.json',dict(decision='H3_G0_IMPLEMENTATION_INVALID',regression_executed=False,next='STOP'));return
 # Controls passed; one deterministic mechanism regression follows.
 for sc in f['scenes']:
  for res in ['B0','B1','B2']:
   projected_grid={}
   for grid in GRIDS:
    p,j,h=saved[sc['scene_id'],grid,res]
    for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
     for cal in ['C0','C1a','C1b','C2','B3_KNOWN_SOURCE_AND_CALIBRATION']:
      oracle=cal.startswith('B3')
      z=projected(p,j,noise,cal,oracle);zh=projected(p,h,noise,cal,oracle)
      projected_grid[grid,noise,cal]=(z,zh)
      for sig in [.01,.05]:
       info.append(dict(scene_id=sc['scene_id'],mesh=grid,resource=res,noise_model=noise,calibration=cal,sigma=sig,scope='DETERMINISTIC_MECHANISM_UPPER_BOUND',**metrics(z,sig)))
      # Independent source spatial chart and dense covariance including reference correlations.
      N=p.shape[1];D=np.zeros((N-1,N));D[:,0]=-1
      for i in range(N-1):D[i,i+1]=1
      L=linalg.block_diag(np.kron(D,np.eye(13)),np.kron(D,np.eye(13)))
      yy=[];ng=[]
      for t in range(3):
       v=np.ones(2*N*13)/2 if noise=='RELATIVE' else np.tile((REF/abs(p[t]).ravel())**2/2,2)
       A=linalg.solve_triangular(linalg.cholesky((L*v)@L.T,lower=True),L,lower=True)
       yy.append(A@rvec(j[t]))
       gg=linalg.block_diag(np.kron(np.eye(N),np.ones((13,1))),np.kron(np.eye(N),np.ones((13,1)))) if cal=='C1a' else np.eye(2*N*13)
       ng.append(A@gg)
      yy=np.concatenate(yy)
      zz=remove(yy,np.concatenate(ng)) if cal in ['C1a','C1b'] else np.zeros_like(yy) if cal=='C2' else yy
      if not oracle:
       e=np.linalg.norm(zz.T@zz-z.T@z)/max(np.linalg.norm(z.T@z),1e-30)
       ck(ctrl,'cross_chart_'+sc['scene_id']+str(grid)+res+noise+cal,e,1e-7)
      if res=='B1' and not oracle:
       raw=projected(p,j,noise,'B3_KNOWN_SOURCE_AND_CALIBRATION',True)
       df=raw.T@raw-z.T@z;den=max(np.linalg.norm(raw.T@raw),1e-30)
       mineig=np.linalg.eigvalsh((df+df.T)/2)[0]/den
       ck(ctrl,'DPI_'+sc['scene_id']+str(grid)+noise+cal,max(0.,-mineig),1e-7)
   p0=saved[sc['scene_id'],80001,res][0];p1=saved[sc['scene_id'],160001,res][0]
   rawdiff=float(np.linalg.norm(p0-p1)/np.linalg.norm(p1));resp=response_distance(p0,p1)
   for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
    nr=response_noise(p1,.01,noise);response_ratio=resp/max(nr,1e-30)
    for cal in ['C0','C1a','C1b','C2','B3_KNOWN_SOURCE_AND_CALIBRATION']:
     a,ah=projected_grid[80001,noise,cal];b,bh=projected_grid[160001,noise,cal]
     jac=float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-30))
     step=max(float(np.linalg.norm(a-ah)/max(np.linalg.norm(ah),1e-30)),float(np.linalg.norm(b-bh)/max(np.linalg.norm(bh),1e-30)))
     zz=[metrics(x,.01)['depth_information_horizontal_profiled'] for x in [a,ah,b,bh]]
     eff=(max(zz)-min(zz))/max(max(zz),1e-30)
     ranks=[metrics(x,.01)['rank'] for x in [a,ah,b,bh]]
     good=response_ratio<=.5 and jac<=.02 and step<=.02 and eff<=.1 and len(set(ranks))==1
     stable.append(dict(scene_id=sc['scene_id'],resource=res,noise_model=noise,calibration=cal,raw_field_relative=rawdiff,source_free_response_difference=resp,response_to_noise1_ratio=response_ratio,Jacobian_mesh_relative=jac,depth_step_relative=step,depth_information_relative=eff,rank_agreement=len(set(ranks))==1,status='STABLE' if good else 'VERTICAL_NUMERICAL_FIDELITY_INCOMPLETE'))
  # Finite F: exact registered 25 nodes, all5 labels; no continuous support/search claim.
  mods=[modal.parse_mod(paths.path_for(ff,160001)) for ff in RAW]
  st=np.array(sc['state']);ur=st[:2]/np.linalg.norm(st[:2]);uv=st[2:]/np.linalg.norm(st[2:])
  for dr in [-250.,-125.,0.,125.,250.]:
   for dv in [-.05,-.025,0.,.025,.05]:
    candst=st.copy();candst[:2]+=dr*ur;candst[2:]+=dv*uv
    if time.monotonic()-start>1800:raise TimeoutError('registered1800s total budget')
    for label in LABELS:
     for res in ['B0','B1','B2']:
      nom=saved[sc['scene_id'],160001,res][0]
      q,_=field(mods,candst,res,label,False)
      for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
       finite.append(dict(scene_id=sc['scene_id'],resource=res,noise_model=noise,range_offset_m=dr,speed_offset_mps=dv,depth_label_m=label,C0_joint_chart_scale_1pct=contrast(nom,q,noise,'C0'),C1b_joint_chart_scale_1pct=contrast(nom,q,noise,'C1b'),scope='FINITE_F_GRID_ONLY_NOT_COVERAGE'))
  if time.monotonic()-start>1800:raise TimeoutError('registered1800s total budget')
  print('MECHANISM',sc['scene_id'],flush=True)
 write(O/'STRUCTURAL_NULLSPACE_CHECKS.csv',[{k:r[k] for k in ['scene_id','mesh','resource','noise_model','calibration','sigma','rank','structural_null_residual','depth_null_loading','depth_local_scale_m']} for r in info]);write(O/'VERTICAL_INFORMATION_BY_SCENE.csv',info);write(O/'NUMERICAL_STABILITY.csv',stable);write(O/'FINITE_HORIZONTAL_LABEL_DIAGNOSTIC.csv',finite);write(O/'CONTROL_CHECKS.csv',ctrl)
 for r in info:ck(ctrl,'structural_null_'+r['scene_id']+str(r['mesh'])+r['resource']+r['noise_model']+r['calibration']+str(r['sigma']),r['structural_null_residual'],1e-9)
 write(O/'CONTROL_CHECKS.csv',ctrl)
 main=[r for r in info if r['resource']=='B1' and r['calibration']=='C1b' and r['sigma']==.01 and r['mesh']==160001]
 mainst=[r for r in stable if r['resource']=='B1' and r['calibration']=='C1b']
 matching={(r['scene_id'],r['noise_model']):r for r in info if r['resource']=='B2' and r['calibration']=='C1b' and r['sigma']==.01 and r['mesh']==160001}
 resource=[]
 for r in main:
  b=matching[r['scene_id'],r['noise_model']]
  resource.append(dict(scene_id=r['scene_id'],noise_model=r['noise_model'],B1_depth_information=r['depth_information_horizontal_profiled'],B2_depth_information=b['depth_information_horizontal_profiled'],B1_to_B2=r['depth_information_horizontal_profiled']/max(b['depth_information_horizontal_profiled'],1e-30),B1_depth10_signal=r['depth10_signal'],scope='LOCAL_HORIZONTAL_PROFILE_NOT_RC2_COVERAGE'))
 write(O/'RESOURCE_MATCHED_CONTROL.csv',resource)
 write(O/'NOISE_AND_CALIBRATION_SENSITIVITY.csv',[r for r in info if r['mesh']==160001 and r['sigma']==.01])
 formula=all(r['PASS'] for r in ctrl);numeric=all(r['status']=='STABLE' for r in mainst)
 strong=all(r['depth10_signal']>=1 for r in main)
 gain=all(r['B1_to_B2']>=1.1 for r in resource)
 if not formula:dec='H3_G0_IMPLEMENTATION_INVALID'
 elif not numeric:dec='H3_G0_VERTICAL_NUMERICAL_FIDELITY_INCOMPLETE'
 elif strong and gain:dec='H3_G0_VERTICAL_DEPTH_INFORMATION_PRESENT_CONDITIONALLY'
 elif not strong and all(r['depth10_signal']>=1 for r in info if r['resource']=='B1' and r['calibration'] in ['C0','C1a'] and r['sigma']==.01 and r['mesh']==160001):dec='H3_G0_CALIBRATION_CONDITIONAL'
 else:dec='H3_G0_CURRENT_4M_APERTURE_NOT_SUPPORTED'
 dump(O/'H3_G0_DECISION.json',dict(decision=dec,parent_SHA=PARENT,design_SHA=git('rev-parse','HEAD'),formula_gate='PASS' if formula else 'FAIL',numerical_gate='PASS' if numeric else 'FAIL',primary_noise_scene_count=len(mainst),primary_numerically_stable=sum(r['status']=='STABLE' for r in mainst),C1b_depth_min=min(r['depth_information_horizontal_profiled'] for r in main),B1_to_B2_min=min(r['B1_to_B2'] for r in resource),C1g_FREE_MODAL_GAIN='DEPTH_ABSORBED_REGULAR_SOURCE_NODES',C1g_PHYSICALLY_CONSTRAINED_GROUP_GAIN='NOT_ADMITTED_NO_INDEPENDENT_GROUP_CONTRACT',C1e='NOT_EVALUATED',horizontal_support='H3_FULL_HORIZONTAL_SUPPORT_NOT_ESTABLISHED',regression_executed=True,new_KRAKEN=0,new_FIELD=0,new_MC=0,new_audio=0,H3_EXTRACTED_OBSERVABLE='NOT_OPENED',E2='PAUSED_WITH_FORMULA_REPAIR_FAILURE',R4_percent=0,next='H3_EXTRACTION_REVIEW' if dec=='H3_G0_VERTICAL_DEPTH_INFORMATION_PRESENT_CONDITIONALLY' else 'H3_A_REVIEW' if dec in ['H3_G0_CALIBRATION_CONDITIONAL','H3_G0_CURRENT_4M_APERTURE_NOT_SUPPORTED'] else 'STOP',next_stage_execution='NOT_AUTHORIZED',elapsed_s=time.monotonic()-start))
 print(json.dumps(native(read(O/'H3_G0_DECISION.json')),indent=2))
if __name__=='__main__':execute()
