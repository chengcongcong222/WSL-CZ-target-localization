"""Exploratory nine-pair pilot only; no new modal solver, audio, RNG, or performance gate."""
from pathlib import Path
import os,time,json
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy import linalg
import r4_e2_g0 as old
import r4_e2_paired_recovery as paired
import r4_e2_forward_diagnostic as base
O=Path('results/R4_E2_NINE_PAIR_INFORMATION_PILOT')
PARENT='f596f17c3f18d07f4b85daf8f75382fc8789b1d1'
SCALE=np.array([50000.,1.,2.,2.,200.])
STEP=np.array([.05,1e-6,.00005,.00005,1.])
RAW=paired.RAW;PAIRS=paired.PAIRS
read,dump,write,sha,git=base.read,base.dump,base.write,base.sha,base.git
QS=linalg.null_space(np.ones((1,8))).T
BP=np.zeros((9,13));BM=BP.copy()
for j,(_,_,lo,hi) in enumerate(PAIRS):
 BP[j,RAW.index(lo)]=BP[j,RAW.index(hi)]=1
 BM[j,RAW.index(lo)]=-1;BM[j,RAW.index(hi)]=1
L0=linalg.block_diag(np.kron(QS,np.eye(13)),np.kron(QS,np.eye(13)))
LCA=linalg.block_diag(np.kron(QS,BP),np.kron(QS,BM))
GB=linalg.block_diag(np.kron(np.eye(8),np.ones((13,1))),np.kron(np.eye(8),np.ones((13,1))))
def realvec(a):return np.concatenate([a.real.reshape(-1,*a.shape[2:]),a.imag.reshape(-1,*a.shape[2:])],axis=0)
def basis(n):
 if n.size==0:return np.zeros((n.shape[0],0))
 v=np.linalg.norm(n,axis=0);n=n[:,v>1e-12]/v[v>1e-12]
 if not n.size:return np.zeros((n.shape[0],0))
 u,s,_=linalg.svd(n,full_matrices=False);return u[:,s>s[0]*1e-10]
def remove(j,n):
 q=basis(n);return j-q@(q.T@j)
def cart(x):return np.array([x[0]*np.cos(x[1]),x[0]*np.sin(x[1]),x[2],x[3]])
def state(sc):s=sc['state'];return np.array([np.hypot(s[0],s[1]),np.arctan2(s[1],s[0]),s[2],s[3],200.])
def field(mods,x,mirror):return old.field_state(mods,cart(x),mirror,x[4])
def derivatives(mods,sc):
 x=state(sc);p=field(mods,x,sc['mirror']);jj=[];hh=[]
 for fac,acc in [(1,jj),(.5,hh)]:
  for c,h in enumerate(STEP*fac):
   up=x.copy();dn=x.copy();up[c]+=h;dn[c]-=h
   a=field(mods,up,sc['mirror']);b=field(mods,dn,sc['mirror'])
   acc.append((a-b)/(2*h)/p*SCALE[c])
 _,nodes=old.geometry(sc['state'],sc['mirror'])
 target=np.array(sc['state'][:2])+old.TIMES[:,None]*np.array(sc['state'][2:4])
 v=target[:,None,:]-nodes
 elems=nodes[:,:,None,:]+np.stack([np.arange(-7.,8.,2.),np.zeros(8)],axis=-1)[None,None,:,:]
 ang=1e-6;mat=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]])
 rp=np.linalg.norm(nodes[:,:,None,:]+(v@mat.T)[:,:,None,:]-elems,axis=-1)
 rm=np.linalg.norm(nodes[:,:,None,:]+(v@mat)[:,:,None,:]-elems,axis=-1)
 ba=(np.stack([old.pressure(m,rp,200.) for m in mods],axis=-1)-np.stack([old.pressure(m,rm,200.) for m in mods],axis=-1))/(2*ang)/p
 return p,np.stack(jj,axis=-1),np.stack(hh,axis=-1),ba
def system(p,j,b,rep,nnode,noise,ref):
 L=L0 if rep=='P0' else LCA;size=L.shape[0];y=[];gain=[];bear=[];trans=[]
 for t in range(3):
  for node in range(nnode):
   var=np.ones(208)/2 if noise=='RELATIVE' else np.tile((ref/abs(p[t,node]).ravel())**2/2,2)
   cov=(L*var[None,:])@L.T
   ch=linalg.cholesky(cov,lower=True);A=linalg.solve_triangular(ch,L,lower=True)
   y.append(A@realvec(j[t,node]))
   gg=np.zeros((size,nnode*16));gg[:,node*16:(node+1)*16]=A@GB;gain.append(gg)
   bb=np.zeros((size,3*nnode));bb[:,t*nnode+node]=A@realvec(b[t,node]);bear.append(bb)
   trans.append(A)
 return np.concatenate(y),np.concatenate(gain),np.concatenate(bear),trans
def metrics(j,sigma):
 w=j/sigma;s=linalg.svdvals(w);rank=int(np.sum(s>s[0]*1e-10)) if s[0]>0 else 0
 _,_,vh=linalg.svd(w,full_matrices=False)
 def prof(col,keep):
  rem=remove(w[:,col:col+1],w[:,keep]);v=float(np.sum(rem**2))
  return v/SCALE[col]**2 if v>(max(s[0],1e-30)*1e-10)**2 else 0.
 ir=prof(0,[1,2,3,4]);iz=prof(4,[0,1,2,3])
 return dict(range_information_depth_profiled=ir,range_information_exact_depth=prof(0,[1,2,3]),depth_information_horizontal_profiled=iz,depth_information_exact_range=prof(4,[1,2,3]),vx_information_profiled=prof(2,[0,1,3,4]),vy_information_profiled=prof(3,[0,1,2,4]),scaled_rank=rank,scaled_condition=(float((s[0]/s[-1])**2) if rank==5 else 'INF'),scaled_weak_direction=json.dumps(vh[-1].tolist()),range_local_scale_m=(1/np.sqrt(ir) if ir>0 else 'INF'),depth_local_scale_m=(1/np.sqrt(iz) if iz>0 else 'INF'),range250_signal_sqrt_information=250*np.sqrt(ir),depth10_signal_sqrt_information=10*np.sqrt(iz),information_trace=float(np.sum(w*w)),range_depth_F_coupling=float((w.T@w)[0,4]))
def covariance_crosscheck(p,j,noise,ref):
 # Independent reference-element contrasts and dense covariance, different chart from orthonormal QS.
 D=np.zeros((7,8));D[:,0]=-1
 for k in range(7):D[k,k+1]=1
 L=linalg.block_diag(np.kron(D,BP),np.kron(D,BM));Lr=linalg.block_diag(np.kron(D,np.eye(13)),np.kron(D,np.eye(13)))
 rj=realvec(j);v=np.ones(208)/2 if noise=='RELATIVE' else np.tile((ref/abs(p).ravel())**2/2,2)
 def gram(M):
  C=(M*v)@M.T;yy=M@rj;return yy.T@linalg.solve(C,yy,assume_a='pos')
 return gram(L),gram(Lr)
def ca_formula_error(p,j):
 worst=0.
 for t in range(3):
  for node in range(2):
   for _,_,lo,hi in PAIRS:
    il=RAW.index(lo);ih=RAW.index(hi);a=p[t,node,:,ih]*p[t,node,:,il].conj();n=np.linalg.norm(a);u=a/n
    rel=j[t,node,:,ih]+j[t,node,:,il].conj();da=a[:,None]*rel
    du=(da-u[:,None]*np.real(u.conj()@da)[None,:])/n
    analytic=np.r_[QS@rel.real,QS@rel.imag]
    chart=np.r_[QS@(du/u[:,None]).real,QS@(du/u[:,None]).imag]
    worst=max(worst,float(np.linalg.norm(analytic-chart)/max(np.linalg.norm(analytic),1e-30)))
 return worst
def candidate_profile(p,q,trans,nnode,noise_model='RELATIVE',ref=1.):
 dr=np.log(abs(q/p))+1j*np.unwrap(np.angle(q/p),axis=2)
 ds=[];gain=[]
 for t in range(3):
  for node in range(nnode):
   A=trans[t*nnode+node];ds.append(A@realvec(dr[t,node]))
   gg=np.zeros((A.shape[0],nnode*16));gg[:,node*16:(node+1)*16]=A@GB;gain.append(gg)
 d=np.concatenate(ds);G=np.concatenate(gain)
 beta=linalg.lstsq(G,d,cond=1e-10,lapack_driver='gelsd')[0]
 corrected=q.copy()
 for node in range(nnode):corrected[:,node]*=np.exp(-beta[node*16:node*16+8])[None,:,None]
 residual=d-G@beta
 u=base.pair_vectors(p[:,:nnode],PAIRS,RAW);v=base.pair_vectors(q[:,:nnode],PAIRS,RAW);c=base.pair_vectors(corrected[:,:nnode],PAIRS,RAW)
 rms=lambda a,b:float(np.sqrt(np.mean(base.ca_distance(a,b)**2)))
 nv=[]
 for _,_,lo,hi in PAIRS:
  il=RAW.index(lo);ih=RAW.index(hi);a=p[:,:nnode,:,ih]*p[:,:nnode,:,il].conj()
  w=abs(a)**2/np.sum(abs(a)**2,axis=-1,keepdims=True)
  vv=np.ones_like(w)*2 if noise_model=='RELATIVE' else ref**2*(1/abs(p[:,:nnode,:,il])**2+1/abs(p[:,:nnode,:,ih])**2)
  nv.append(2*.01**2*np.sum(w*(1-w)*vv,axis=-1))
 noise=float(np.sqrt(np.mean(nv)))
 return corrected,dict(C0_CA_Frobenius_rms=rms(u,v),C1_CA_Frobenius_rms=rms(u,c),C0_joint_chart_scale_1pct=float(np.linalg.norm(d)/.01),C1_joint_chart_scale_1pct=float(np.linalg.norm(residual)/.01),C1_joint_chart_scale_5pct=float(np.linalg.norm(residual)/.05),CA_marginal_noise1_rms=float(noise),CA_marginal_noise5_rms=float(5*noise),gain_profile='JOINT_FIXED_GAIN_GLS_LOG_CHART_FEASIBLE_NOT_GLOBAL_CA_FROB_MINIMUM',gain_log_amplitudes=json.dumps(beta.reshape(nnode,16)[:,:8].tolist()))
def controls():
 out=[];p=np.exp(.001j*np.arange(3*2*8*13).reshape(3,2,8,13))*(1+.001*np.arange(8)[None,None,:,None])
 j=np.zeros(p.shape+(5,),complex)
 for c in [0,4]:j[...,c]=(.07+.11j)*np.arange(13)[None,None,None,:]
 j[...,1]=1j*np.arange(8)[None,None,:,None]*np.array(RAW)[None,None,None,:]*.01
 b=j[...,1].copy()
 for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
  for rep in ['P0','P1','P2']:
   nn=2 if rep in ['P0','P2'] else 1
   y,g,bb,_=system(p,j,b,rep,nn,noise,.00010824612978073539)
   z=remove(y,bb)
   out.append(dict(control='plane_range_depth_and_bearing_profile',representation=rep,noise=noise,error=float(np.linalg.norm(z)/max(np.linalg.norm(y),1e-30)),PASS=bool(np.linalg.norm(z)<1e-10*max(np.linalg.norm(y),1))))
 # Finite fixed complex gains: phases cancel CA; log amplitudes identical for all times/frequencies.
 gain=np.exp(np.linspace(-.1,.1,16).reshape(2,8)+1j*np.linspace(-1,1,16).reshape(2,8))
 q=p*gain[None,:,:,None]
 _,_,_,tr=system(p,j,b,'P2',2,'RELATIVE',1)
 _,rr=candidate_profile(p,q,tr,2)
 out.append(dict(control='fixed_gain_joint_absorption',representation='P2',noise='RELATIVE',error=rr['C1_CA_Frobenius_rms'],PASS=rr['C1_CA_Frobenius_rms']<1e-10))
 qfree=p*np.exp(.01*np.arange(p.size).reshape(p.shape)+.003j*np.arange(p.size).reshape(p.shape))
 absorbed=qfree*(p/qfree);e=float(np.max(abs(absorbed-p))/np.max(abs(p)))
 out.append(dict(control='C2_free_response_absorption',representation='ALL',noise='BOTH',error=e,PASS=e<1e-12))
 out.append(dict(control='shared_frequency_covariance_offdiagonal',representation='P2',noise='RELATIVE',error=float(np.max(abs(BP@BP.T-np.diag(np.diag(BP@BP.T))))),PASS=bool(np.max(abs(BP@BP.T-np.diag(np.diag(BP@BP.T))))>0)))
 return out
def execute():
 if (O/'EXECUTION_STARTED.json').exists():raise RuntimeError('once only')
 f=read(O/'DESIGN_FREEZE.json');assert git('log','-1','--pretty=%s')=='R4 E2: freeze nine-pair nonbearing information pilot'
 assert git('rev-parse','HEAD')==git('ls-remote','origin','refs/heads/main').split()[0]
 for name,r in f['bindings'].items():assert sha(name,r['kind']=='RAW')==r['sha256'],name
 dump(O/'EXECUTION_STARTED.json',dict(design_SHA=git('rev-parse','HEAD'),epoch=time.time(),new_KRAKEN=0,new_MC=0,new_audio=0))
 ref=f['absolute_noise_reference_amplitude'];ctrl=controls();write(O/'CONTROL_CHECKS.csv',ctrl)
 if not all(r['PASS'] for r in ctrl):
  dump(O/'PILOT_DECISION.json',dict(decision='IMPLEMENTATION_INVALID',reason='PRE_SCIENCE_CONTROLS_FAIL',R4_percent=0,original_E2_admission='FAIL_UNCHANGED',E2_G0_information='NOT_EVALUATED',pilot='EXPLORATORY_PILOT_ONLY',next='STOP'));return
 start=time.monotonic();info=[];stable=[];cal=[];cand=[];saved={};bearing=[]
 reps=[('P0_SINGLE',1,'P0'),('P1',1,'P1'),('P0_DUAL',2,'P0'),('P2',2,'P2')]
 for sc in read(old.OUT/'E2_SCENE_FREEZE.json')['scenes']:
  x=state(sc);T=np.eye(4);T[:2,:2]=[[np.cos(x[1]),-x[0]*np.sin(x[1])],[np.sin(x[1]),x[0]*np.cos(x[1])]]
  fb=T.T@old.baseline(sc)@T
  bearing.append(dict(scene_id=sc['scene_id'],polar_state_FIM=json.dumps(fb.tolist()),epochs=121,depth_information=0,no_fusion=True,role='INHERITED_BEARING_REFERENCE_ONLY_DIFFERENT_RESOURCE_ASSUMPTIONS'))
  inputs={}
  for grid in [80001,160001]:
   mods=[old.parse_mod(paired.path_for(ff,grid)) for ff in RAW]
   p,j,h,b=derivatives(mods,sc);inputs[grid]=(p,j,h,b,mods)
   np.savez_compressed(O/f"INPUT_{sc['scene_id']}_n{grid}.npz",pressure=p,J=j,J_half=h,bearing_tangent=b)
   err=ca_formula_error(p,j);ctrl.append(dict(control='analytic_CA_chart_formula_'+sc['scene_id']+'_'+str(grid),representation='P2',noise='BOTH',error=err,PASS=err<1e-9))
   for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
    for rep,nn,kind in reps:
     y,g,bb,tr=system(p,j,b,kind,nn,noise,ref);yh,_,_,_=system(p,h,b,kind,nn,noise,ref)
     if rep in ['P0_DUAL','P2']:
      fc,fr=covariance_crosscheck(p[0,0],j[0,0],noise,ref)
      L=LCA if kind=='P2' else L0;fixed=(fc if kind=='P2' else fr)
      vv=np.ones(208)/2 if noise=='RELATIVE' else np.tile((ref/abs(p[0,0]).ravel())**2/2,2)
      ch=linalg.cholesky((L*vv)@L.T,lower=True);yy=linalg.solve_triangular(ch,L@realvec(j[0,0]),lower=True)
      er=float(np.linalg.norm(yy.T@yy-fixed)/max(np.linalg.norm(fixed),1e-30))
      ctrl.append(dict(control='dense_chart_covariance_'+sc['scene_id']+str(grid),representation=rep,noise=noise,error=er,PASS=er<1e-9))
     for c in ['C0','C1','C2']:
      for bear in [False,True]:
       N=g if c=='C1' else np.zeros((len(y),0))
       if bear:N=np.column_stack([N,bb])
       z=remove(y,N) if c!='C2' else np.zeros_like(y);zh=remove(yh,N) if c!='C2' else np.zeros_like(yh)
       saved[sc['scene_id'],grid,rep,noise,c,bear]=(z,zh)
       for sigma in [.01,.05]:
        rr=metrics(z,sigma);info.append(dict(scene_id=sc['scene_id'],mesh=grid,representation=rep,noise_model=noise,calibration=c,bearing_profiled=bear,sigma=sigma,pilot_status='EXPLORATORY_PILOT_ONLY',array_norm_min=float(np.linalg.norm(p[:,:nn],axis=2).min()),array_power_mean=float(np.mean(np.sum(abs(p[:,:nn])**2,axis=2))),**rr))
     if kind=='P2':
      _,_,_,tr=system(p,j,b,kind,nn,noise,ref)
      for name,state0,zz,scope in old.candidates(sc)[-8:]:
       q=old.field_state(mods,state0,sc['mirror'],zz);corrected,r=candidate_profile(p,q,tr,nn,noise,ref)
       np.savez_compressed(O/f"CAND_{sc['scene_id']}_{name}_{grid}_{noise}.npz",pressure=q,corrected=corrected)
       cand.append(dict(scene_id=sc['scene_id'],mesh=grid,noise_model=noise,candidate_id=name,scope=scope,**r,numerical_triangle_bound='',numerical_interval_positive='',rank=0))
  for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
   for rep,nn,kind in reps:
    for c in ['C0','C1','C2']:
     for bear in [False,True]:
      z,h=saved[sc['scene_id'],80001,rep,noise,c,bear];w,wh=saved[sc['scene_id'],160001,rep,noise,c,bear]
      ja=max(np.linalg.norm(z-h)/max(np.linalg.norm(h),1e-30),np.linalg.norm(w-wh)/max(np.linalg.norm(wh),1e-30))
      me=np.linalg.norm(z-w)/max(np.linalg.norm(w),1e-30)
      m0=metrics(z,.01);m1=metrics(w,.01);mh=metrics(wh,.01)
      ir=[m0['range_information_depth_profiled'],m1['range_information_depth_profiled'],mh['range_information_depth_profiled']]
      er=(max(ir)-min(ir))/max(max(ir),1e-30)
      iz=[m0['depth_information_horizontal_profiled'],m1['depth_information_horizontal_profiled'],mh['depth_information_horizontal_profiled']]
      ez=(max(iz)-min(iz))/max(max(iz),1e-30)
      ok=ja<.02 and me<.02 and er<.1 and ez<.1
      stable.append(dict(scene_id=sc['scene_id'],representation=rep,noise_model=noise,calibration=c,bearing_profiled=bear,fd_relative=float(ja),mesh_relative=float(me),range_information_relative=float(er),depth_information_relative=float(ez),status='NUMERICALLY_STABLE' if ok else 'LOCAL_INFORMATION_NUMERICALLY_UNRESOLVED'))
      if c!='C2':
       r0,_=saved[sc['scene_id'],160001,'P0_SINGLE' if nn==1 else 'P0_DUAL',noise,c,bear]
       d=r0.T@r0-w.T@w;mn=float(np.linalg.eigvalsh((d+d.T)/2)[0]);den=max(float(np.linalg.norm(r0.T@r0)),1e-30)
       ctrl.append(dict(control='data_processing_'+sc['scene_id']+'_160001',representation=rep,noise=noise+'_'+c+'_'+str(bear),error=mn/den,PASS=mn>=-1e-7*den))
       lowraw,_=saved[sc['scene_id'],80001,'P0_SINGLE' if nn==1 else 'P0_DUAL',noise,c,bear]
       dd=lowraw.T@lowraw-z.T@z;mm=float(np.linalg.eigvalsh((dd+dd.T)/2)[0]);dn=max(float(np.linalg.norm(lowraw.T@lowraw)),1e-30)
       ctrl.append(dict(control='data_processing_'+sc['scene_id']+'_80001',representation=rep,noise=noise+'_'+c+'_'+str(bear),error=mm/dn,PASS=mm>=-1e-7*dn))
   c0,_=saved[sc['scene_id'],160001,'P2',noise,'C0',True];c1,_=saved[sc['scene_id'],160001,'P2',noise,'C1',True]
   cal.append(dict(scene_id=sc['scene_id'],noise_model=noise,C0_trace=float(np.sum(c0*c0)),C1_trace=float(np.sum(c1*c1)),gain_information_retained=float(np.sum(c1*c1)/max(np.sum(c0*c0),1e-30)),C2_trace=0.,C2='SATURATED_SNAPSHOT_FREQUENCY_ELEMENT_RESPONSE'))
   for name,state0,zz,scope in old.candidates(sc)[-8:]:
    a,b=inputs[80001][0],inputs[160001][0]
    q0=np.load(O/f"CAND_{sc['scene_id']}_{name}_80001_{noise}.npz")['corrected']
    q1=np.load(O/f"CAND_{sc['scene_id']}_{name}_160001_{noise}.npz")['corrected']
    rms=lambda aa,bb:float(np.sqrt(np.mean(base.ca_distance(base.pair_vectors(aa,PAIRS,RAW),base.pair_vectors(bb,PAIRS,RAW))**2)))
    bound=rms(a,b)+rms(q0,q1)
    rs=[r for r in cand if r['scene_id']==sc['scene_id'] and r['noise_model']==noise and r['candidate_id']==name]
    sep=min(r['C1_CA_Frobenius_rms'] for r in rs)
    for r in rs:r.update(numerical_triangle_bound=bound,numerical_interval_positive=bool(sep>bound))
   for grid in [80001,160001]:
    rs=sorted([r for r in cand if r['scene_id']==sc['scene_id'] and r['noise_model']==noise and r['mesh']==grid],key=lambda r:(r['C1_CA_Frobenius_rms'],r['candidate_id']))
    for i,r in enumerate(rs):r['rank']=i+1
  np.savez_compressed(O/f"PROFILED_{sc['scene_id']}.npz",**{'_'.join(map(str,k[1:])):v[0] for k,v in saved.items() if k[0]==sc['scene_id']})
  if time.monotonic()-start>1800:raise TimeoutError('frozen execution budget')
  print('SCENE',sc['scene_id'],flush=True)
 si={(r['scene_id'],r['representation'],r['noise_model'],r['calibration'],r['bearing_profiled']):r['status'] for r in stable}
 for r in info:r['local_information_status']=si[r['scene_id'],r['representation'],r['noise_model'],r['calibration'],r['bearing_profiled']]
 write(O/'BEARING_REFERENCE.csv',bearing)
 write(O/'LOCAL_INFORMATION_BY_SCENE.csv',info);write(O/'NUMERICAL_INFORMATION_STABILITY.csv',stable);write(O/'CALIBRATION_PROFILE.csv',cal);write(O/'CANDIDATE_SEPARATION.csv',cand);write(O/'CONTROL_CHECKS.csv',ctrl)
 low=[r for r in info if r['representation']=='P2' and r['calibration']=='C1' and r['bearing_profiled'] and r['mesh']==160001]
 write(O/'LOW_ENERGY_NOISE_SENSITIVITY.csv',low)
 primary=[r for r in info if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C1' and r['bearing_profiled'] and r['mesh']==160001 and r['sigma']==.01]
 upper=[r for r in info if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C0' and r['bearing_profiled'] and r['mesh']==160001 and r['sigma']==.01]
 st=[r for r in stable if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C1' and r['bearing_profiled']]
 localok=all(r['status']=='NUMERICALLY_STABLE' for r in st)
 strong=lambda rs:all(r['range250_signal_sqrt_information']>=1 and r['depth10_signal_sqrt_information']>=1 for r in rs)
 cs=[r for r in cand if r['mesh']==160001 and r['noise_model']=='RELATIVE']
 candidateok=all(r['numerical_interval_positive'] and r['C1_CA_Frobenius_rms']>r['CA_marginal_noise1_rms'] and r['C1_joint_chart_scale_1pct']>1 for r in cs)
 if not all(r['PASS'] for r in ctrl):decision='IMPLEMENTATION_INVALID'
 elif localok and strong(primary) and candidateok:decision='NINE_PAIR_NONBEARING_MECHANISM_PROMISING'
 elif localok and all(r['status']=='NUMERICALLY_STABLE' for r in stable if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C0' and r['bearing_profiled']) and strong(upper) and not strong(primary):decision='NINE_PAIR_CALIBRATION_CONDITIONAL'
 else:decision='NINE_PAIR_INCREMENT_WEAK_OR_NOT_RESOLVED'
 result=dict(decision=decision,parent_SHA=PARENT,design_SHA=git('rev-parse','HEAD'),local_primary_stable_scenes=sum(r['status']=='NUMERICALLY_STABLE' for r in st),strong_C1_scenes=sum(r['range250_signal_sqrt_information']>=1 and r['depth10_signal_sqrt_information']>=1 for r in primary),C1_range_information_min=min(r['range_information_depth_profiled'] for r in primary),C1_depth_information_min=min(r['depth_information_horizontal_profiled'] for r in primary),candidate_primary_resolved=sum(r['numerical_interval_positive'] and r['C1_CA_Frobenius_rms']>r['CA_marginal_noise1_rms'] for r in cs),candidate_primary_count=len(cs),control_checks=len(ctrl),control_failures=sum(not r['PASS'] for r in ctrl),original_E2_admission='FAIL_UNCHANGED',E2_G0_information='NOT_EVALUATED',pilot_status='EXPLORATORY_PILOT_ONLY',R4_percent=0,new_KRAKEN=0,new_MC=0,new_audio=0,next='FULL_BAND_AND_EXTRACTION_RESEARCH_REVIEW' if decision=='NINE_PAIR_NONBEARING_MECHANISM_PROMISING' else 'H3_REVIEW' if decision!='IMPLEMENTATION_INVALID' else 'STOP',next_stage_execution='NOT_AUTHORIZED',elapsed_s=time.monotonic()-start,independent_audit='PENDING')
 dump(O/'PILOT_DECISION.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':execute()
