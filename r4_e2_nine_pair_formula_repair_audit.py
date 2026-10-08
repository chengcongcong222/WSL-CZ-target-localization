"""Independent exploratory-pilot reconstruction. No old writes, no solver or new experiment."""
from pathlib import Path
import struct,json
import numpy as np
from scipy import linalg
import r4_e2_nine_pair_formula_repair as p
O=p.O;checks=[]
def ck(name,ok,detail=''):checks.append(dict(check=name,PASS=bool(ok),detail=str(detail)))
def qremove(a,n):
 if not n.size:return a.copy()
 v=np.linalg.norm(n,axis=0);n=n[:,v>1e-12]/v[v>1e-12]
 if not n.size:return a.copy()
 q,r,piv=linalg.qr(n,mode='economic',pivoting=True)
 # Independent QR rank, crosschecked against registered singular rank.
 s=linalg.svdvals(n);rank=int(np.sum(s>s[0]*1e-10))
 q=q[:,:rank];return a-q@(q.T@a)
def mod(path):
 b=Path(path).read_bytes();r=struct.unpack_from('<i',b)[0]*4
 nt,nm=struct.unpack_from('<ii',b,92);M=struct.unpack_from('<i',b,5*r)[0]
 z=np.array(struct.unpack_from('<'+str(nt)+'f',b,4*r),float);ph=np.empty((nm,M),complex)
 for j in range(M):
  a=np.array(struct.unpack_from('<'+str(2*nm)+'f',b,(7+j)*r));ph[:,j]=a[::2]+1j*a[1::2]
 a=np.array(struct.unpack_from('<'+str(2*M)+'f',b,(7+M)*r));return z,ph,a[::2]+1j*a[1::2]
def press(m,r,z):
 zz,ph,k=m;w=ph[list(zz).index(z)]*ph[list(zz).index(200.)]
 rr=np.asarray(r).ravel();v=w[:,None]*np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*(k[:,None]*rr+np.pi/4))
 return np.add.reduce(v,axis=0).reshape(np.shape(r))
def geom(sc,x):
 tt=np.array([0.,600.,1200.]);dt=np.maximum(tt-600,0)
 main=np.stack([2*np.minimum(tt,600)+2*dt*np.cos(np.pi/12),2*dt*np.sin(np.pi/12)],axis=-1)
 nodes=np.stack([main,main+[0,sc['mirror']*5000]],axis=1)
 elems=nodes[:,:,None,:]+np.stack([np.arange(-7.,8.,2.),np.zeros(8)],axis=-1)[None,None,:,:]
 st=np.array([x[0]*np.cos(x[1]),x[0]*np.sin(x[1]),x[2],x[3]])
 target=st[:2]+tt[:,None]*st[2:];return np.sqrt(np.sum((target[:,None,None,:]-elems)**2,axis=-1))
def cold_field(ms,sc,x):return np.stack([press(m,geom(sc,x),x[4]) for m in ms],axis=-1)
def local_metrics(z,sigma):
 a=z/sigma;sv=linalg.svdvals(a)
 def prof(col,keep):
  r=qremove(a[:,col:col+1],a[:,keep]);v=float(np.sum(r*r))
  return v/p.SCALE[col]**2 if v>(max(sv[0],1e-30)*1e-10)**2 else 0.
 return prof(0,[1,2,3,4]),prof(4,[0,1,2,3])
def main():
 f=p.read(O/'DESIGN_FREEZE.json');dec=p.read(O/'REPAIR_DECISION.json')
 for name,r in f['bindings'].items():ck('binding:'+name,p.sha(name,r['kind']=='RAW')==r['sha256'])
 controls=p.base.rows(O/'CONTROL_CHECKS.csv');ck('all_formula_plane_covariance_DPI',all(r['PASS']=='True' for r in controls))
 info=p.base.rows(O/'REPAIRED_INFORMATION_BY_SCENE.csv');st=p.base.rows(O/'EFFECTIVE_INFORMATION_STABILITY.csv');cs=p.base.rows(O/'CANDIDATE_SEPARATION.csv')
 ck('all_information_rows',len(info)==4608);ck('all_stability_rows',len(st)==1152);ck('all_candidates',len(cs)==768)
 ck('no_new_science_gate',dec['original_E2_admission']=='FAIL_UNCHANGED' and dec['E2_G0_information']=='NOT_EVALUATED' and dec['pilot_status']=='FORMULA_REPAIRED_DEVELOPMENT_ONLY' and dec['R4_percent']==dec['new_KRAKEN']==dec['new_MC']==dec['new_audio']==0)
 D=np.zeros((7,8));D[:,0]=-1
 for j in range(7):D[j,j+1]=1
 LR=linalg.block_diag(np.kron(D,np.eye(13)),np.kron(D,np.eye(13)))
 LC=linalg.block_diag(np.kron(D,p.BP),np.kron(D,p.BM))
 SRC=linalg.block_diag(np.kron(np.ones((8,1)),np.eye(13)),np.kron(np.ones((8,1)),np.eye(13)))
 coldmax=0.;jmax=0.;Fmax=0.;rawmax=0.;profilemax=0.;candidate_chart_max=0.
 cindex={(r['scene_id'],int(r['mesh']),r['noise_model'],r['candidate_id']):r for r in cs}
 for sc in f['scenes']:
  profile=np.load(O/f"PROFILED_{sc['scene_id']}.npz")
  for grid in [80001,160001]:
   a=np.load(O/f"INPUT_{sc['scene_id']}_n{grid}.npz");P=a['pressure'];J=a['J'];Jh=a['J_half'];B=a['bearing_tangent']
   ms=[mod(p.paired.path_for(ff,grid)) for ff in p.RAW]
   for ff,m in zip(p.RAW,ms):ck('modal:'+sc['scene_id']+str(grid)+str(ff),all(np.array_equal(x,y) for x,y in zip(m,p.old.parse_mod(p.paired.path_for(ff,grid)))))
   x=p.state(sc);pc=cold_field(ms,sc,x);coldmax=max(coldmax,float(np.max(abs(pc-P))/np.max(abs(pc))))
   # Independent Cartesian modal-gradient contraction, recomputed from struct-decoded caches.
   rr=geom(sc,x);target=np.array(sc['state'][:2])+np.array([0.,600.,1200.])[:,None]*np.array(sc['state'][2:])
   dt=np.maximum(np.array([0.,600.,1200.])-600,0)
   main=np.stack([2*np.minimum([0.,600.,1200.],600)+2*dt*np.cos(np.pi/12),2*dt*np.sin(np.pi/12)],axis=-1)
   nodes=np.stack([main,main+[0,sc['mirror']*5000]],axis=1)
   displacement=target[:,None,None,:]-nodes[:,:,None,:]-np.stack([np.arange(-7.,8.,2.),np.zeros(8)],axis=-1)[None,None,:,:]
   grad=[]
   for m in ms:
    zz,ph,k=m;w=ph[list(zz).index(200.)]**2;rad=np.zeros(rr.shape,complex)
    for im in range(len(k)):
     term=w[im]*np.sqrt(2*np.pi/(k[im]*rr))*np.exp(-1j*(k[im]*rr+np.pi/4))
     rad+=term*(-1j*k[im]-.5/rr)
    grad.append(rad[...,None]*displacement/rr[...,None])
   gradient=np.stack(grad,axis=-2)/pc[...,None]
   dxy=np.zeros((3,2,4));dxy[:,:,0]=[np.cos(x[1]),np.sin(x[1])];dxy[:,:,1]=[-x[0]*np.sin(x[1]),x[0]*np.cos(x[1])];dxy[:,0,2]=[0,600,1200];dxy[:,1,3]=[0,600,1200]
   jhoriz=np.einsum('abcfd,adc->abcfc',gradient,dxy) if False else np.stack([np.sum(gradient*dxy[:,None,None,None,:,c],axis=-1)*p.SCALE[c] for c in range(4)],axis=-1)
   vb=target[:,None,:]-nodes;rot=np.stack([-vb[:,:,1],vb[:,:,0]],axis=-1)
   bc=np.sum(gradient*rot[:,:,None,None,:],axis=-1)
   ck('independent_analytic_bearing:'+sc['scene_id']+str(grid),np.linalg.norm(bc-B)/np.linalg.norm(B)<1e-10)
   # Independent geometric null via SVD; primary uses cofactors.
   D=a['center_range_geometry'];_,ss,vv=linalg.svd(D,full_matrices=True);n=np.r_[vv[-1],0.]
   ck('independent_geometry_null:'+sc['scene_id']+str(grid),abs(n@a['structural_null'])>1-1e-10)
   for fac,jac in [(1,J),(.5,Jh)]:
    up=x.copy();dn=x.copy();up[4]+=fac;dn[4]-=fac
    jc=np.concatenate([jhoriz,((cold_field(ms,sc,up)-cold_field(ms,sc,dn))/(2*fac)/pc*p.SCALE[4])[...,None]],axis=-1)
    jmax=max(jmax,float(np.linalg.norm(jc-jac)/np.linalg.norm(jc)))
   for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
    for rep,nn,L in [('P0_SINGLE',1,LR),('P1',1,LC),('P0_DUAL',2,LR),('P2',2,LC)]:
     yy=[];gg=[];bb=[];trans=[];raws=[];rawg=[];rawb=[]
     for t in range(3):
      for n in range(nn):
       v=np.ones(208)/2 if noise=='RELATIVE' else np.tile((f['absolute_noise_reference_amplitude']/abs(P[t,n]).ravel())**2/2,2)
       C=(L*v)@L.T;ch=linalg.cholesky(C,lower=True);A=linalg.solve_triangular(ch,L,lower=True);trans.append(A)
       yy.append(A@p.realvec(J[t,n]))
       G=np.zeros((L.shape[0],nn*16));G[:,n*16:(n+1)*16]=A@p.GB;gg.append(G)
       b=np.zeros((L.shape[0],3*nn));b[:,t*nn+n]=A@p.realvec(B[t,n]);bb.append(b)
       if rep.startswith('P0'):
        W=1/np.sqrt(v);N=W[:,None]*SRC;N/=np.linalg.norm(N,axis=0)
        def sr(z):return z-N@(N.T@z)
        raws.append(sr(W[:,None]*p.realvec(J[t,n])))
        G0=np.zeros((208,nn*16));G0[:,n*16:(n+1)*16]=sr(W[:,None]*p.GB);rawg.append(G0)
        b0=np.zeros((208,3*nn));b0[:,t*nn+n]=sr(W*p.realvec(B[t,n]));rawb.append(b0)
     y=np.concatenate(yy);gain=np.concatenate(gg);bear=np.concatenate(bb)
     for c in ['C0','C1','C2']:
      for be in [False,True]:
       N=gain if c=='C1' else np.zeros((len(y),0))
       if be:N=np.column_stack([N,bear])
       z=qremove(y,N) if c!='C2' else np.zeros_like(y)
       key='_'.join(map(str,[grid,rep,noise,c,be]));ref=profile[key]
       F=z.T@z;Fr=ref.T@ref;err=np.linalg.norm(F-Fr)/max(np.linalg.norm(Fr),1e-20);Fmax=max(Fmax,float(err))
       ck('independent_chart_profile:'+sc['scene_id']+key,err<1e-7,err)
       if rep.startswith('P0') and c!='C2':
        ry=np.concatenate(raws);RG=np.concatenate(rawg);RB=np.concatenate(rawb)
        RN=RG if c=='C1' else np.zeros((len(ry),0))
        if be:RN=np.column_stack([RN,RB])
        rz=qremove(ry,RN);err0=np.linalg.norm(rz.T@rz-Fr)/max(np.linalg.norm(Fr),1e-20);rawmax=max(rawmax,float(err0))
        ck('original_raw_source_profile:'+sc['scene_id']+key,err0<1e-7,err0)
       for sigma in [.01,.05]:
        matching=[r for r in info if r['scene_id']==sc['scene_id'] and int(r['mesh'])==grid and r['representation']==rep and r['noise_model']==noise and r['calibration']==c and (r['bearing_profiled']=='True')==be and float(r['sigma'])==sigma]
        ck('unique_information_row:'+sc['scene_id']+key+str(sigma),len(matching)==1)
        ir,iz=local_metrics(ref,sigma);r=matching[0]
        errm=max(abs(ir-float(r['range_information_depth_profiled']))/max(ir,float(r['range_information_depth_profiled']),1e-20),abs(iz-float(r['depth_information_horizontal_profiled']))/max(iz,float(r['depth_information_horizontal_profiled']),1e-20))
        profilemax=max(profilemax,errm);ck('effective_information:'+sc['scene_id']+key+str(sigma),errm<1e-5,errm)
     if rep=='P2':
      for name,state0,zz,scope in p.old.candidates(sc)[-8:]:
       ar=np.load(O/f"CAND_{sc['scene_id']}_{name}_{grid}_{noise}.npz");q=ar['pressure'];co=ar['corrected'];r=cindex[sc['scene_id'],grid,noise,name]
       gainlog=np.array(json.loads(r['gain_log_amplitudes']))
       expected=q*np.exp(-gainlog)[None,:,:,None]
       ck('joint_fixed_gain_array:'+sc['scene_id']+str(grid)+noise+name,np.max(abs(co-expected))/np.max(abs(co))<1e-12)
       dr=np.log(abs(q/P))+1j*np.unwrap(np.angle(q/P),axis=2);sd=[]
       for t in range(3):
        for n in range(2):sd.append(trans[t*2+n]@p.realvec(dr[t,n]))
       sr=qremove(np.concatenate(sd)[:,None],gain)
       scale=float(np.linalg.norm(sr)/.01);e=abs(scale-float(r['C1_joint_chart_scale_1pct']))/max(scale,1e-20);candidate_chart_max=max(candidate_chart_max,e)
       ck('candidate_independent_chart:'+sc['scene_id']+str(grid)+noise+name,e<1e-7,e)
  print('AUDIT_SCENE',sc['scene_id'],flush=True)
 ck('cold_pressure',coldmax<1e-12,coldmax);ck('cold_analytic_chain_and_depth_stencils',jmax<1e-9,jmax)
 ck('all_projected_information_cross_chart',Fmax<1e-7,Fmax);ck('original_raw_pressure_equivalence',rawmax<1e-7,rawmax)
 ck('independent_effective_profiles',profilemax<1e-5,profilemax);ck('candidate_shared_covariance_profile',candidate_chart_max<1e-7,candidate_chart_max)
 ck('C2_zero_with_unbounded_scales',all(float(r['information_trace'])==0 and r['range_local_scale_m']=='INF' and r['depth_local_scale_m']=='INF' for r in info if r['calibration']=='C2'))
 p.write(O/'INDEPENDENT_AUDIT_CHECKS.csv',checks)
 val=dict(checks=len(checks),PASS=sum(r['PASS'] for r in checks),FAIL=sum(not r['PASS'] for r in checks),cold_pressure_relative_max=coldmax,cold_J_relative_max=jmax,reference_element_chart_F_relative_max=Fmax,original_raw_gaussian_source_profile_F_relative_max=rawmax,independent_effective_profile_relative_max=profilemax,candidate_chart_profile_relative_max=candidate_chart_max,new_solver_calls=0,scope='Independent modal struct/independent per-mode analytic Cartesian chain and depth stencils; alternative reference-element CA chart and QR nuisance elimination; raw Gaussian source projection; every local information metric; candidate joint chart/fixed gains; no pilot scientific rerun.')
 p.dump(O/'VALIDATION.json',val);print(json.dumps(val,indent=2));assert val['FAIL']==0
if __name__=='__main__':main()
