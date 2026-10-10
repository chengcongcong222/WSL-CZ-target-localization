"""HLA-H2 pure kinematics, observation-only starts, smooth branch LS."""
from pathlib import Path
from dataclasses import dataclass
import hashlib,json,itertools
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_SINGLE_HLA_RADIAL_INCREMENT_PRECISION_PILOT'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/HLA_H2')
PARENT='dd13a0a32a4bd279c03c272b9ee8b86ceea14730'
AXES=('r_km','theta_deg','v_mps','psi_deg')
ERROR_NAMES=('relative_range_error','absolute_bearing_error_deg','relative_speed_error','absolute_heading_error_deg')
LOW=np.array([45.,-5.,1.,-15.])
HIGH=np.array([60.,5.,3.,15.]);WIDTH=HIGH-LOW
CENTER=np.array([52.5,0.,2.,0.])
TIMES=np.arange(0.,1201.,10.)
ENDPOINTS=np.arange(0,121,20)
WINDOW_SECONDS=200.
SIG_B=np.radians(.1)
SIGMAS=(.02,.05,.10,.20)
SIGNS=np.array(list(itertools.product((-1.,1.),repeat=6)))
B_RANGE_CENTERS=45+15*(np.arange(16)+.5)/16
T_B=float(chi2.ppf(.975,121));T_Q=float(chi2.ppf(.975,6))
SVD_RCOND=1e-12
DEDUP=np.array([.005,.001,.005,.05,.005])
def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [clean(x) for x in v]
    if isinstance(v,np.ndarray):return clean(v.tolist())
    if isinstance(v,(np.floating,float)):
        if np.isnan(v):return None
        if np.isinf(v):return 'INF' if v>0 else '-INF'
        return float(v)
    if isinstance(v,np.integer):return int(v)
    if isinstance(v,np.bool_):return bool(v)
    return v
def write_json(path,obj):Path(path).write_text(json.dumps(clean(obj),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def wrap(a):return (np.asarray(a)+np.pi)%(2*np.pi)-np.pi
def platform(t=TIMES):
    t=np.asarray(t);d=np.maximum(t-600,0)
    return np.stack((2*np.minimum(t,600)+2*d*np.cos(np.radians(15)),2*d*np.sin(np.radians(15))),axis=-1)
PLATFORM=platform()
def geometry(states,jac=False):
    s=np.atleast_2d(np.asarray(states,dtype=float));r,th,v,ps=s.T
    th=np.radians(th);ps=np.radians(ps)
    x=1000*r[:,None]*np.cos(th[:,None])+v[:,None]*TIMES*np.cos(ps[:,None])-PLATFORM[:,0]
    y=1000*r[:,None]*np.sin(th[:,None])+v[:,None]*TIMES*np.sin(ps[:,None])-PLATFORM[:,1]
    rho=np.hypot(x,y);beta=np.arctan2(y,x);q=np.diff(rho[:,ENDPOINTS],axis=1)/200
    if not jac:return beta,q,rho
    dx=np.empty((*x.shape,4));dy=np.empty_like(dx);deg=np.pi/180
    dx[:,:,0]=1000*np.cos(th[:,None]);dy[:,:,0]=1000*np.sin(th[:,None])
    dx[:,:,1]=-1000*r[:,None]*np.sin(th[:,None])*deg;dy[:,:,1]=1000*r[:,None]*np.cos(th[:,None])*deg
    dx[:,:,2]=TIMES*np.cos(ps[:,None]);dy[:,:,2]=TIMES*np.sin(ps[:,None])
    dx[:,:,3]=-v[:,None]*TIMES*np.sin(ps[:,None])*deg;dy[:,:,3]=v[:,None]*TIMES*np.cos(ps[:,None])*deg
    jb=(x[:,:,None]*dy-y[:,:,None]*dx)/rho[:,:,None]**2
    jr=(x[:,:,None]*dx+y[:,:,None]*dy)/rho[:,:,None]
    jq=np.diff(jr[:,ENDPOINTS],axis=1)/200
    return beta,q,rho,jb,jq
def error_vector(states,true):
    s=np.asarray(states,dtype=float);true=np.asarray(true,dtype=float)
    e=np.abs(s-true);e[...,0]/=true[0];e[...,2]/=true[2]
    e[...,1]=np.degrees(abs(wrap(np.radians(s[...,1]-true[1]))))
    e[...,3]=np.degrees(abs(wrap(np.radians(s[...,3]-true[3]))))
    return e
def nearest_rank(values,p):
    x=np.sort(np.asarray(values,dtype=float))
    return float(x[max(0,int(np.ceil(p*len(x)))-1)]) if len(x) else np.nan
@dataclass(frozen=True)
class Observation:
    beta:np.ndarray
    radial:np.ndarray|None
    sigma:float|None
    kind:str
    def __post_init__(self):
        assert self.beta.shape==(121,) and np.isfinite(self.beta).all()
        assert self.kind in ('B','S','U','U-mismatch','UB-matched')
        if self.kind!='B':
            assert self.radial.shape==(6,) and np.isfinite(self.radial).all() and self.sigma>0
            if self.kind!='S':assert np.all(self.radial>=0)
def scaled_lstsq(A,y,scales):
    matrix=A*scales[None,:]/52500.;rhs=y/52500.
    u,res,rank,sv=np.linalg.lstsq(matrix,rhs,rcond=SVD_RCOND)
    return u*scales,dict(rank=int(rank),columns=A.shape[1],condition=float(sv[0]/sv[-1]) if sv[-1]>0 else np.inf,min_singular=float(sv[-1]))
def convert_start(cart,details,bias=False):
    p=cart[:2];vel=cart[2:4]
    raw=np.array([np.linalg.norm(p)/1000,np.degrees(np.arctan2(p[1],p[0])),np.linalg.norm(vel),np.degrees(np.arctan2(vel[1],vel[0]))])
    h=np.clip(raw,LOW,HIGH);details=dict(details,clipped=bool(np.any(h!=raw)),raw_state=raw.tolist())
    if bias:
        b=float(np.clip(cart[-1],-.2,.2));details['bias_clipped']=bool(b!=cart[-1])
        return np.r_[h,b],details
    return h,details
def linear_radial_start(beta,qobs,bias=False):
    tt=TIMES[ENDPOINTS];angles=np.asarray(beta)[ENDPOINTS];e=np.column_stack((np.cos(angles),np.sin(angles)))
    d=np.r_[0.,200*np.cumsum(qobs)];n=7
    A=np.zeros((n,2,6 if bias else 5))
    A[:,0,0]=1;A[:,1,1]=1;A[:,0,2]=tt;A[:,1,3]=tt;A[:,:,4]=-e
    if bias:A[:,:,5]=tt[:,None]*e
    rhs=platform(tt)+d[:,None]*e
    scales=np.array([52500.,52500.,2.,2.,52500.]+([.2] if bias else []))
    cart,meta=scaled_lstsq(A.reshape(14,-1),rhs.ravel(),scales)
    meta['rho0_linear_m']=float(cart[4]);meta['nonpositive_rho0']=bool(cart[4]<=0)
    if bias:
        initial,meta=convert_start(np.r_[cart[:4],cart[5]],meta,True)
    else:initial,meta=convert_start(cart[:4],meta)
    return initial,meta
def bearing_starts(beta):
    e=np.asarray(beta);s=np.sin(e);co=np.cos(e)
    A=np.column_stack((s,-co,TIMES*s,-TIMES*co))
    y=PLATFORM[:,0]*s-PLATFORM[:,1]*co
    cart,meta=scaled_lstsq(A,y,np.array([52500.,52500.,2.,2.]))
    linear,meta=convert_start(cart,meta)
    starts=[('LINEAR',linear,meta)]
    for i,r in enumerate(B_RANGE_CENTERS):starts.append((f'RANGE_{i:02}',np.array([r,0.,2.,0.]),dict(rank=4,columns=4,condition=np.nan,min_singular=np.nan,clipped=False,raw_state=[r,0.,2.,0.])))
    return starts
def profile_bias(q,u):
    q=np.asarray(q);u=np.asarray(u)
    breaks=np.unique(np.r_[-.2,.2,np.clip(-q,-.2,.2)])
    candidates=list(breaks)
    for lo,hi in zip(breaks[:-1],breaks[1:]):
        signs=np.sign(q+(lo+hi)/2)
        optimum=float(np.mean(signs*u-q))
        candidates.append(float(np.clip(optimum,lo,hi)))
    candidates=np.array(sorted(set(candidates)));cost=np.square(abs(q[None,:]+candidates[:,None])-u).sum(axis=1)
    best=int(np.argmin(cost));return float(candidates[best]),float(cost[best])
class Evaluator:
    def __init__(self,obs,unfolded):
        self.obs=obs;self.unfolded=unfolded;self.bias=obs.kind=='UB-matched'
        self.last=None;self.r=None;self.j=None;self.geometry_calls=0
    def evaluate(self,u):
        if self.last is not None and np.array_equal(u,self.last):return
        h=LOW+WIDTH*u[:4];b=-.2+.4*u[4] if self.bias else 0.
        beta,q,rho,jb,jq=geometry(h,True)
        residual=wrap(beta[0]-self.obs.beta)/SIG_B
        jac=jb[0]*WIDTH[None,:]/SIG_B
        if self.obs.kind!='B':
            qr=(q[0]+b-self.unfolded)/self.obs.sigma
            residual=np.r_[residual,qr];jac=np.r_[jac,jq[0]*WIDTH[None,:]/self.obs.sigma]
        if self.bias:
            jac=np.column_stack((jac,np.r_[np.zeros(121),np.full(6,.4/self.obs.sigma)]))
        self.last=np.array(u,copy=True);self.r=residual;self.j=jac;self.geometry_calls+=1
    def fun(self,u):self.evaluate(u);return self.r.copy()
    def derivative(self,u):self.evaluate(u);return self.j.copy()
def fit_start(obs,unfolded,start):
    bias=obs.kind=='UB-matched'
    u=(np.asarray(start)[:4]-LOW)/WIDTH
    if bias:u=np.r_[u,(start[4]+.2)/.4]
    u=np.clip(u,1e-12,1-1e-12);ev=Evaluator(obs,unfolded)
    fit=least_squares(ev.fun,u,jac=ev.derivative,bounds=(np.zeros(len(u)),np.ones(len(u))),method='trf',loss='linear',max_nfev=100,ftol=1e-9,xtol=1e-9,gtol=1e-9)
    state=LOW+WIDTH*fit.x[:4];b=-.2+.4*fit.x[4] if bias else 0.
    return state,b,fit,ev.geometry_calls
def physical_scores(obs,state,b=0.,unfolded=None):
    beta,q,_=geometry(state);jb=float(np.square(wrap(beta[0]-obs.beta)/SIG_B).sum())
    if obs.kind=='B':return dict(J_beta=jb,J_branch=0.,J_folded=0.,J_profile=0.,profile_b=0.,accepted_original=bool(jb<=T_B),accepted_profile=bool(jb<=T_B))
    qs=q[0];sigma=obs.sigma
    js=float(np.square((qs+b-unfolded)/sigma).sum())
    jf=float(np.square((abs(qs+b)-obs.radial)/sigma).sum()) if obs.kind!='S' else js
    pb,prof=(profile_bias(qs,obs.radial) if obs.kind=='UB-matched' else (b,jf*sigma**2))
    jp=prof/sigma**2
    return dict(J_beta=jb,J_branch=js,J_folded=jf,J_profile=jp,profile_b=pb,accepted_original=bool(jb<=T_B and jf<=T_Q),accepted_profile=bool(jb<=T_B and jp<=T_Q))
