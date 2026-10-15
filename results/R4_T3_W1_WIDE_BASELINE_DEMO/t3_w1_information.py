"""Nominal information calculation, distinct Cartesian parameterization from fitter.
No Monte Carlo observations or estimator solutions determine these matrices.
"""
import math
import numpy as np
SCALE=np.array([50000.,50000.,2.,2.,math.pi/180*.1,math.pi/180*.05])

def mean_variance(state,t,P,quality):
    x=np.asarray(state);delta=(x[:2]+t[:,None]*x[2:4])[:,None,:]-P
    r2=np.sum(delta*delta,axis=-1)
    mu=np.arctan2(delta[:,:,1],delta[:,:,0])+x[4]+x[5]*np.array([-1.,1.])
    var=np.deg2rad(quality['bearing_sigma_deg'])**2+quality['nav_sigma_m']**2/r2
    return mu.ravel(),var.ravel()

def design(state,t,P,quality):
    d=(state[:2]+t[:,None]*state[2:4])[:,None,:]-P;r2=np.sum(d*d,axis=-1)
    gx=-d[:,:,1]/r2;gy=d[:,:,0]/r2
    H=np.stack([gx,gy,gx*t[:,None],gy*t[:,None],np.ones_like(gx),np.broadcast_to([-1.,1.],gx.shape)],axis=-1).reshape(-1,6)
    dv=np.stack([d[:,:,0],d[:,:,1],d[:,:,0]*t[:,None],d[:,:,1]*t[:,None],np.zeros_like(gx),np.zeros_like(gx)],axis=-1)*(-2*quality['nav_sigma_m']**2/r2**2)[:,:,None]
    var=np.deg2rad(quality['bearing_sigma_deg'])**2+quality['nav_sigma_m']**2/r2
    # Fisher of independent Gaussian angular observations with state-dependent variance.
    A=np.vstack([H*SCALE/np.sqrt(var.ravel())[:,None],dv.reshape(-1,6)*SCALE/(np.sqrt(2)*var.ravel()[:,None])])
    return H,dv.reshape(-1,6),A

def analyze(state,t,P,quality,bias_known=False):
    H,dv,A=design(state,t,P,quality)
    ix=[0,1,2,3] if bias_known else list(range(6));A=A[:,ix];scale=SCALE[ix]
    U,s,V=np.linalg.svd(A,full_matrices=False);rank=int(np.sum(s/s[0]>1e-10));Fscaled=A.T@A;Fphysical=Fscaled/scale[:,None]/scale[None,:]
    nuisance=[i for i,v in enumerate(ix) if v not in [2,3]];velocity=[ix.index(2),ix.index(3)]
    N=A[:,nuisance];u,ss,_=np.linalg.svd(N,full_matrices=False);nr=int(np.sum(ss/ss[0]>1e-10));B=A[:,velocity]-u[:,:nr]@(u[:,:nr].T@A[:,velocity]);eff=B.T@B;ev,vec=np.linalg.eigh(eff);order=np.argsort(ev);ev=ev[order];vec=vec[:,order];effective_rank=int(np.sum(ev>max(ev[-1]*1e-12,0)))
    C=None;speed=math.hypot(*state[2:4]);r=math.hypot(*state[:2]);gradients=[np.array([state[0]/r,state[1]/r,0,0]+([] if bias_known else [0,0])),np.array([-state[1]/r**2,state[0]/r**2,0,0]+([] if bias_known else [0,0])),np.array([0,0,state[2]/speed,state[3]/speed]+([] if bias_known else [0,0])),np.array([0,0,-state[3]/speed**2,state[2]/speed**2]+([] if bias_known else [0,0]))]
    if rank==len(ix):
        C=(V.T*(1/s**2))@V;C=C*scale[:,None]*scale[None,:];errors=[math.sqrt(max(0,float(g@C@g))) for g in gradients]
    else:errors=[math.inf]*4
    weak=V[-1];row={'control':'BIAS_KNOWN_INFORMATION_UPPER_CONTROL' if bias_known else 'BIAS_FREE_NO_PRIOR_LOCAL_DATA_INFORMATION','rank':rank,'parameters':len(ix),'condition_design':float(s[0]/s[-1]) if s[-1]>0 else math.inf,'singular_min':float(s[-1]),'singular_max':float(s[0]),'effective_velocity_rank':effective_rank,'effective_velocity_eigen_min':float(ev[0]),'effective_velocity_eigen_max':float(ev[1]),'sigma_range_m':errors[0],'sigma_bearing_deg':math.degrees(errors[1]),'sigma_speed_mps':errors[2],'speed_rel_1p96sigma':1.96*errors[2]/speed,'sigma_heading_deg':math.degrees(errors[3]),'weak_velocity_vx_scaled':float(vec[0,0]),'weak_velocity_vy_scaled':float(vec[1,0]),'scope':'LOCAL_LINEARIZED_NAV_GAUSSIAN_FISHER;NO_EMPIRICAL_P95_OR_GLOBAL_COVERAGE'}
    row.update({f'weak_full_scaled_component_{ix[i]}':float(weak[i]) for i in range(len(ix))})
    return row,Fphysical,Fscaled,eff,C

def finite_difference_check(state,t,P,quality):
    H,dv,A=design(state,t,P,quality);numeric=np.zeros_like(H);numvar=np.zeros_like(dv)
    for j,scale in enumerate(SCALE):
        step=1e-5*scale;a=state.copy();b=state.copy();a[j]+=step;b[j]-=step;ma,va=mean_variance(a,t,P,quality);mb,vb=mean_variance(b,t,P,quality);diff=np.arctan2(np.sin(ma-mb),np.cos(ma-mb));numeric[:,j]=diff/(2*step);numvar[:,j]=(va-vb)/(2*step)
    angle=float(np.max(abs((numeric-H)*SCALE)));variance=float(np.max(abs((numvar-dv)*SCALE))/np.max(mean_variance(state,t,P,quality)[1]));return angle,variance
