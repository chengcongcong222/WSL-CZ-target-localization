"""Bounded analytic formula repair. Imported pilot is read-only; all writes use the new directory."""
from pathlib import Path
import os,time,json,csv,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import numpy as np
from scipy import linalg
import r4_e2_nine_pair_pilot as legacy
import r4_e2_g0 as old
import r4_e2_paired_recovery as paired
import r4_e2_forward_diagnostic as base
O=Path('results/R4_E2_NINE_PAIR_FORMULA_REPAIR')
OLD=legacy.O
PARENT='2574fd59cc04ae49d2fb9e7669a79b2d80376caa'
SCALE=legacy.SCALE;STEP=legacy.STEP;RAW=legacy.RAW;PAIRS=legacy.PAIRS
QS=legacy.QS;BP=legacy.BP;BM=legacy.BM;L0=legacy.L0;LCA=legacy.LCA;GB=legacy.GB
read,write,sha,git=base.read,base.write,base.sha,base.git
realvec,basis,remove,cart,state=legacy.realvec,legacy.basis,legacy.remove,legacy.cart,legacy.state
system,covariance_crosscheck,ca_formula_error,candidate_profile,controls=legacy.system,legacy.covariance_crosscheck,legacy.ca_formula_error,legacy.candidate_profile,legacy.controls
REF=.00010824612978073539
REPS=[('P0_SINGLE',1,'P0'),('P1',1,'P1'),('P0_DUAL',2,'P0'),('P2',2,'P2')]
def native(x):
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,dict):return {k:native(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [native(v) for v in x]
 return x
def dump(path,obj):base.dump(path,native(obj))
def scene_geometry(sc):
 x=state(sc);rr,nodes=old.geometry(sc['state'],sc['mirror'])
 target=cart(x)[:2]+old.TIMES[:,None]*cart(x)[2:]
 v=target[:,None,:]-nodes
 offsets=np.stack([np.arange(-7.,8.,2.),np.zeros(8)],axis=-1)
 er=(v[:,:,None,:]-offsets[None,None,:,:])/rr[...,None]
 geo=np.zeros((3,2,4))
 geo[:,:,0]=[np.cos(x[1]),np.sin(x[1])]
 geo[:,:,1]=[-x[0]*np.sin(x[1]),x[0]*np.cos(x[1])]
 geo[:,0,2]=old.TIMES;geo[:,1,3]=old.TIMES
 return rr,v,er,geo
def analytic(mods,sc):
 x=state(sc);rr,v,er,geo=scene_geometry(sc)
 ps=[];rs=[]
 for zz,phi,k in mods:
  w=phi[list(zz).index(200.)]**2;r=rr.ravel()
  terms=w[:,None]*np.sqrt(2*np.pi/(k[:,None]*r))*np.exp(-1j*k[:,None]*r-1j*np.pi/4)
  ps.append(np.sum(terms,axis=0).reshape(rr.shape))
  rs.append(np.sum(terms*(-1j*k[:,None]-1/(2*r)),axis=0).reshape(rr.shape))
 p=np.stack(ps,axis=-1);rad=np.stack(rs,axis=-1)
 grad=rad[...,None]*er[:,:,:,None,:]/p[...,None]
 j=np.zeros(p.shape+(5,),complex)
 j[...,:4]=np.einsum('tnmfd,tdc->tnmfc',grad,geo)*SCALE[:4]
 h=j.copy()
 for step,acc in [(1.,j),(.5,h)]:
  pu=old.field_state(mods,cart(x),sc['mirror'],200.+step)
  pd=old.field_state(mods,cart(x),sc['mirror'],200.-step)
  acc[...,4]=(pu-pd)/(2*step)/p*SCALE[4]
 bv=np.stack([-v[:,:,1],v[:,:,0]],axis=-1)
 b=np.einsum('tnmfd,tnd->tnmf',grad,bv)
 return p,j,h,b
def structural_null(sc):
 rr,v,er,geo=scene_geometry(sc)
 unit=v[:,0,:]/np.linalg.norm(v[:,0,:],axis=1)[:,None]
 D=np.einsum('td,tdc->tc',unit,geo)*SCALE[:4]
 # Cofactors of row-normalized geometric matrix, independent of field/SVD.
 dn=D/np.linalg.norm(D,axis=1)[:,None]
 h=np.array([(-1.)**i*np.linalg.det(np.delete(dn,i,axis=1)) for i in range(4)])
 h/=np.linalg.norm(h)
 h=np.r_[h,0.]
 return h,D
def quotient(j,null=None):
 # Never used to make the structural control pass: only after residual is measured.
 if null is None:return j
 q=linalg.null_space(null[None,:])
 return (j@q)@q.T
def metrics(j,sigma):
 w=j/sigma;_,s,vh=linalg.svd(w,full_matrices=False)
 rank=int(np.sum(s>s[0]*1e-10)) if s[0]>0 else 0
 null=vh[rank:]
 def prof(col,keep):
  # An estimand with a structural null component has infinite uncertainty.
  if rank<5 and np.linalg.norm(null[:,col])>1e-8:return 0.
  rem=remove(w[:,col:col+1],w[:,keep]);v=float(np.sum(rem*rem))
  return v/SCALE[col]**2 if v>(max(s[0],1e-30)*1e-10)**2 else 0.
 ir=prof(0,[1,2,3,4]);iz=prof(4,[0,1,2,3])
 return dict(range_information_depth_profiled=ir,range_information_exact_depth=profile_subset(w,0,[1,2,3]),depth_information_horizontal_profiled=iz,depth_information_exact_range=profile_subset(w,4,[1,2,3]),vx_information_profiled=prof(2,[0,1,3,4]),vy_information_profiled=prof(3,[0,1,2,4]),scaled_rank=rank,scaled_condition=float((s[0]/s[-1])**2) if rank==5 else 'INF',scaled_weak_direction=json.dumps(vh[-1].tolist()),range_local_scale_m=1/np.sqrt(ir) if ir>0 else 'INF',depth_local_scale_m=1/np.sqrt(iz) if iz>0 else 'INF',range250_signal_sqrt_information=250*np.sqrt(ir),depth10_signal_sqrt_information=10*np.sqrt(iz),information_trace=float(np.sum(w*w)),range_depth_F_coupling=float((w.T@w)[0,4]))
def profile_subset(w,col,keep):
 rem=remove(w[:,col:col+1],w[:,keep]);v=float(np.sum(rem*rem));s=linalg.svdvals(w)
 return v/SCALE[col]**2 if v>(max(s[0],1e-30)*1e-10)**2 else 0.
def ck(rows,name,rep,noise,error,passed):
 rows.append(dict(control=name,representation=rep,noise=noise,error=float(error),PASS=bool(passed)))
def precontrols(f):
 rows=controls();struct=[];fd=[];comparison=[]
 for sc in f['scenes']:
  null,D=structural_null(sc);Q=linalg.null_space(null[None,:])
  for grid in [80001,160001]:
   mods=[old.parse_mod(paired.path_for(ff,grid)) for ff in RAW]
   p,j,h,b=analytic(mods,sc)
   np.savez_compressed(O/f"INPUT_{sc['scene_id']}_n{grid}.npz",pressure=p,J=j,J_half=h,bearing_tangent=b,structural_null=null,center_range_geometry=D)
   orig=np.load(OLD/f"INPUT_{sc['scene_id']}_n{grid}.npz")
   ep=np.max(abs(p-orig['pressure']))/np.max(abs(p))
   ck(rows,'unchanged_pressure_'+sc['scene_id']+str(grid),'ALL','BOTH',ep,ep<1e-12)
   geomres=np.linalg.norm(D@null[:4])/np.linalg.norm(D)
   ck(rows,'cofactor_geometric_null_'+sc['scene_id']+str(grid),'P1','BOTH',geomres,geomres<1e-12)
   for rep,nn,kind in REPS:
    for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
     y,g,bb,_=system(p,j,b,kind,nn,noise,REF)
     yh,_,_,_=system(p,h,b,kind,nn,noise,REF)
     yf,gf,bf,_=system(p,orig['J'],orig['bearing_tangent'],kind,nn,noise,REF)
     for cal in ['C0','C1']:
      n=g if cal=='C1' else np.zeros((len(y),0))
      n=np.column_stack([n,bb]);nf=gf if cal=='C1' else np.zeros((len(yf),0));nf=np.column_stack([nf,bf])
      z=remove(y,n);zh=remove(yh,n);zf=remove(yf,nf)
      sv=linalg.svdvals(z);rank=int(np.sum(sv>sv[0]*1e-10))
      nr=np.linalg.norm(z@null)/max(np.linalg.norm(z),1e-30) if nn==1 else 0.
      nhr=np.linalg.norm(zh@null)/max(np.linalg.norm(zh),1e-30) if nn==1 else 0.
      if nn==1:
       passed=rank<=4 and nr<1e-10 and nhr<1e-10
       ck(rows,'structural_rank_null_'+sc['scene_id']+str(grid)+cal,rep,noise,max(nr,nhr),passed)
       struct.append(dict(scene_id=sc['scene_id'],mesh=grid,representation=rep,calibration=cal,noise_model=noise,numeric_rank=rank,rank_limit=4,weak_singular_ratio=float(sv[-1]/sv[0]),null_response=nr,half_depth_null_response=nhr,geometry_null_response=geomres,null_vector=json.dumps(null.tolist()),PASS=passed))
      # Project old FD into the independently derived physically identifiable state subspace.
      e=float(np.linalg.norm((z-zf)@Q)/max(np.linalg.norm(z@Q),1e-30)) if nn==1 else float(np.linalg.norm(z-zf)/np.linalg.norm(z))
      ck(rows,'analytic_fd_identifiable_'+sc['scene_id']+str(grid)+cal,rep,noise,e,e<.02)
      fd.append(dict(scene_id=sc['scene_id'],mesh=grid,representation=rep,calibration=cal,noise_model=noise,identifiable_relative_error=e,original_rank=int(np.sum(linalg.svdvals(zf)>linalg.svdvals(zf)[0]*1e-10)),analytic_rank=rank,old_null_response=float(np.linalg.norm(zf@null)/np.linalg.norm(zf)) if nn==1 else '',depth_step_relative=float(np.linalg.norm(z-zh)/max(np.linalg.norm(zh),1e-30)),PASS=e<.02))
      mr=metrics(z,.01);mf=legacy.metrics(zf,.01)
      comparison.append(dict(scene_id=sc['scene_id'],mesh=grid,representation=rep,calibration=cal,noise_model=noise,original_rank=mf['scaled_rank'],repaired_rank=mr['scaled_rank'],original_range_information=mf['range_information_depth_profiled'],repaired_range_information=mr['range_information_depth_profiled'],original_depth_information=mf['depth_information_horizontal_profiled'],repaired_depth_information=mr['depth_information_horizontal_profiled'],range_relative_change=abs(mr['range_information_depth_profiled']-mf['range_information_depth_profiled'])/max(mr['range_information_depth_profiled'],mf['range_information_depth_profiled'],1e-30),depth_relative_change=abs(mr['depth_information_horizontal_profiled']-mf['depth_information_horizontal_profiled'])/max(mr['depth_information_horizontal_profiled'],mf['depth_information_horizontal_profiled'],1e-30)))
   e=ca_formula_error(p,j);ck(rows,'CA_exact_chart_'+sc['scene_id']+str(grid),'P2','BOTH',e,e<1e-9)
   print('FORMULA_CONTROL',sc['scene_id'],grid,flush=True)
 write(O/'STRUCTURAL_NULLSPACE_CHECKS.csv',struct);write(O/'ANALYTIC_VS_FD_JACOBIAN.csv',fd);write(O/'ORIGINAL_VS_REPAIRED_COMPARISON.csv',comparison);write(O/'PRE_DEVELOPMENT_CONTROLS.csv',rows)
 return rows
def derivatives(mods,sc):
 # Prepared once by analytic pre-controls; no finite horizontal differences.
 grid=ACTIVE_GRID
 a=np.load(O/f"INPUT_{sc['scene_id']}_n{grid}.npz")
 return a['pressure'],a['J'],a['J_half'],a['bearing_tangent']
def weak_stability(a,b,null=None):
 if null is not None:
  q=linalg.null_space(null[None,:]);a=a@q;b=b@q
 _,sa,va=linalg.svd(a,full_matrices=False);_,sb,vb=linalg.svd(b,full_matrices=False)
 ra=int(np.sum(sa>sa[0]*1e-10)) if sa[0]>0 else 0
 rb=int(np.sum(sb>sb[0]*1e-10)) if sb[0]>0 else 0
 if ra==rb==0:return 0.,0.,True
 if ra!=rb:return 1.,1.,False
 err=float(abs(sa[ra-1]-sb[rb-1])/max(sa[ra-1],sb[rb-1],1e-30))
 # Compare weakest identifiable direction, not the known structural null.
 dot=min(1.,abs(float(va[ra-1]@vb[rb-1])))
 return err,float(np.sqrt(max(0.,1-dot*dot))),True

def execute():
 global ACTIVE_GRID
 if (O/'EXECUTION_STARTED.json').exists():raise RuntimeError('once only')
 f=read(O/'DESIGN_FREEZE.json');assert git('log','-1','--pretty=%s')=='R4 E2: freeze analytic Jacobian and physical rank repair'
 assert git('rev-parse','HEAD')==git('ls-remote','origin','refs/heads/main').split()[0]
 for name,r in f['bindings'].items():assert sha(name,r['kind']=='RAW')==r['sha256'],name
 dump(O/'EXECUTION_STARTED.json',dict(design_SHA=git('rev-parse','HEAD'),epoch=time.time(),new_KRAKEN=0,new_MC=0,new_audio=0))
 ref=f['absolute_noise_reference_amplitude'];ctrl=precontrols(f);write(O/'CONTROL_CHECKS.csv',ctrl)
 if not all(r['PASS'] for r in ctrl):
  dump(O/'REPAIR_DECISION.json',dict(decision='E2_FORMULA_REPAIR_FAILED',reason='PRE_SCIENCE_CONTROLS_FAIL',R4_percent=0,original_E2_admission='FAIL_UNCHANGED',E2_G0_information='NOT_EVALUATED',pilot='FORMULA_REPAIRED_DEVELOPMENT_ONLY',next='STOP'));return
 start=time.monotonic();info=[];stable=[];cal=[];cand=[];saved={};bearing=[]
 reps=[('P0_SINGLE',1,'P0'),('P1',1,'P1'),('P0_DUAL',2,'P0'),('P2',2,'P2')]
 for sc in read(old.OUT/'E2_SCENE_FREEZE.json')['scenes']:
  x=state(sc);T=np.eye(4);T[:2,:2]=[[np.cos(x[1]),-x[0]*np.sin(x[1])],[np.sin(x[1]),x[0]*np.cos(x[1])]]
  fb=T.T@old.baseline(sc)@T
  bearing.append(dict(scene_id=sc['scene_id'],polar_state_FIM=json.dumps(fb.tolist()),epochs=121,depth_information=0,no_fusion=True,role='INHERITED_BEARING_REFERENCE_ONLY_DIFFERENT_RESOURCE_ASSUMPTIONS'))
  inputs={}
  for grid in [80001,160001]:
   ACTIVE_GRID=grid
   mods=[old.parse_mod(paired.path_for(ff,grid)) for ff in RAW]
   p,j,h,b=derivatives(mods,sc);inputs[grid]=(p,j,h,b,mods)

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
        rr=metrics(z,sigma);info.append(dict(scene_id=sc['scene_id'],mesh=grid,representation=rep,noise_model=noise,calibration=c,bearing_profiled=bear,sigma=sigma,pilot_status='FORMULA_REPAIRED_DEVELOPMENT_ONLY',array_norm_min=float(np.linalg.norm(p[:,:nn],axis=2).min()),array_power_mean=float(np.mean(np.sum(abs(p[:,:nn])**2,axis=2))),**rr))
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
      null=structural_null(sc)[0] if nn==1 and bear and c!='C2' else None
      we,wa,wr=weak_stability(z,w,null)
      he,ha,hr=weak_stability(w,wh,null)
      ok=ja<.02 and me<.02 and er<.1 and ez<.1 and max(we,he)<.1 and max(wa,ha)<.1 and wr and hr
      stable.append(dict(scene_id=sc['scene_id'],representation=rep,noise_model=noise,calibration=c,bearing_profiled=bear,fd_relative=float(ja),mesh_relative=float(me),weak_singular_relative=max(we,he),weak_direction_sine=max(wa,ha),rank_agreement=bool(wr and hr),range_information_relative=float(er),depth_information_relative=float(ez),status='NUMERICALLY_STABLE' if ok else 'LOCAL_INFORMATION_NUMERICALLY_UNRESOLVED'))
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
 write(O/'REPAIRED_INFORMATION_BY_SCENE.csv',info);write(O/'EFFECTIVE_INFORMATION_STABILITY.csv',stable);write(O/'CALIBRATION_PROFILE.csv',cal);write(O/'CANDIDATE_SEPARATION.csv',cand);write(O/'CONTROL_CHECKS.csv',ctrl)
 low=[r for r in info if r['representation']=='P2' and r['calibration']=='C1' and r['bearing_profiled'] and r['mesh']==160001]
 write(O/'LOW_ENERGY_NOISE_SENSITIVITY.csv',low)
 primary=[r for r in info if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C1' and r['bearing_profiled'] and r['mesh']==160001 and r['sigma']==.01]
 upper=[r for r in info if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C0' and r['bearing_profiled'] and r['mesh']==160001 and r['sigma']==.01]
 st=[r for r in stable if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C1' and r['bearing_profiled']]
 localok=all(r['status']=='NUMERICALLY_STABLE' for r in st)
 strong=lambda rs:all(r['range250_signal_sqrt_information']>=1 and r['depth10_signal_sqrt_information']>=1 for r in rs)
 cs=[r for r in cand if r['mesh']==160001 and r['noise_model']=='RELATIVE']
 candidateok=all(r['numerical_interval_positive'] and r['C1_CA_Frobenius_rms']>r['CA_marginal_noise1_rms'] and r['C1_joint_chart_scale_1pct']>1 for r in cs)
 if not all(r['PASS'] for r in ctrl):decision='E2_FORMULA_REPAIR_FAILED'
 elif localok and strong(primary) and candidateok:decision='E2_NINE_PAIR_ANALYTIC_CHAIN_REPAIRED_CONDITIONALLY_PROMISING'
 elif localok and all(r['status']=='NUMERICALLY_STABLE' for r in stable if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C0' and r['bearing_profiled']) and strong(upper) and not strong(primary):decision='E2_NINE_PAIR_CALIBRATION_CONDITIONAL'
 else:decision='E2_NINE_PAIR_INFORMATION_NOT_ESTABLISHED'
 result=dict(decision=decision,parent_SHA=PARENT,design_SHA=git('rev-parse','HEAD'),local_primary_stable_scenes=sum(r['status']=='NUMERICALLY_STABLE' for r in st),strong_C1_scenes=sum(r['range250_signal_sqrt_information']>=1 and r['depth10_signal_sqrt_information']>=1 for r in primary),C1_range_information_min=min(r['range_information_depth_profiled'] for r in primary),C1_depth_information_min=min(r['depth_information_horizontal_profiled'] for r in primary),candidate_primary_resolved=sum(r['numerical_interval_positive'] and r['C1_CA_Frobenius_rms']>r['CA_marginal_noise1_rms'] for r in cs),candidate_primary_count=len(cs),control_checks=len(ctrl),control_failures=sum(not r['PASS'] for r in ctrl),original_E2_admission='FAIL_UNCHANGED',E2_G0_information='NOT_EVALUATED',pilot_status='FORMULA_REPAIRED_DEVELOPMENT_ONLY',R4_percent=0,new_KRAKEN=0,new_MC=0,new_audio=0,next='E2_RESEARCH_REVIEW' if decision=='E2_NINE_PAIR_ANALYTIC_CHAIN_REPAIRED_CONDITIONALLY_PROMISING' else 'H3_REVIEW' if decision!='E2_FORMULA_REPAIR_FAILED' else 'H3_REVIEW',next_stage_execution='NOT_AUTHORIZED',elapsed_s=time.monotonic()-start,independent_audit='PENDING')
 dump(O/'REPAIR_DECISION.json',result);print(json.dumps(native(result),indent=2))

if __name__=='__main__':execute()
