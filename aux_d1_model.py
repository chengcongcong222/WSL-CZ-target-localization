"""AUX-D1 observation-only bounded multistart processor. No generator imports."""
import time,math
import numpy as np
from scipy.optimize import least_squares
LOW=np.array([45.,-5.,1.,-15.]);SPAN=np.array([15.,10.,2.,30.]);HIGH=LOW+SPAN
SIGMA=np.deg2rad(.05);A=1.5;CENTER=LOW+SPAN/2
FIELDS=['r_km','theta_deg','v_mps','psi_deg','Q','solver_cost','success','status','nfev','njev','fun_calls','jac_calls','geometry_calls','optimality','boundary','elapsed_s']
def wrap(x):return (x+np.pi)%(2*np.pi)-np.pi

def score(residual):
 e=np.ascontiguousarray(residual,dtype=np.float64).ravel()
 return float(np.sum(2*e*e/(np.sqrt(1+(e/A)**2)+1),dtype=np.float64)/e.size)

def polar_state(cart,origin):
 p=np.asarray(cart[:2])-origin;v=np.asarray(cart[2:]);return np.array([np.hypot(*p)/1000,np.rad2deg(np.arctan2(p[1],p[0])),np.hypot(*v),np.rad2deg(np.arctan2(v[1],v[0]))])

def initialize_main(t,P,y):
 n=np.column_stack([-np.sin(y[:,0]),np.cos(y[:,0])]);matrix=np.column_stack([n,n*t[:,None]]);scale=np.array([50000.,50000.,2.,2.])
 try:
  q,_,rank,sv=np.linalg.lstsq(matrix*scale,np.sum(n*P[:,0],axis=1),rcond=1e-12);cart=q*scale;h=polar_state(cart,P[0,0]);valid=rank==4 and np.isfinite(h).all()
  return np.clip(h if valid else CENTER,LOW,HIGH),{'kind':'MAIN_ONLY_LINEAR','rank':int(rank),'raw_cartesian':cart.tolist(),'raw_h':h.tolist(),'fallback':not valid,'clipped':bool(np.any((h<LOW)|(h>HIGH)))}
 except Exception as e:return CENTER.copy(),{'kind':'MAIN_ONLY_LINEAR','fallback':True,'reason':repr(e)}

def initialize_aux(t,P,y):
 diag={'kind':'DIRECTED_INTERSECTION_5_IRLS','fallback':True}
 try:
  u=np.stack([np.cos(y),np.sin(y)],axis=-1);delta=P[:,1]-P[:,0];det=u[:,0,0]*u[:,1,1]-u[:,0,1]*u[:,1,0]
  with np.errstate(divide='ignore',invalid='ignore'):
   l0=(delta[:,0]*u[:,1,1]-delta[:,1]*u[:,1,0])/det;l1=(delta[:,0]*u[:,0,1]-delta[:,1]*u[:,0,0])/det;points=P[:,0]+l0[:,None]*u[:,0]
  valid=(abs(det)>1e-12)&(l0>0)&(l1>0)&np.isfinite(points).all(axis=1);diag.update(valid_intersections=int(valid.sum()),invalid_intersections=int((~valid).sum()))
  if valid.sum()<4 or np.ptp(t[valid])<=0:return CENTER.copy(),diag
  tt=t[valid];xy=points[valid];T=np.column_stack([np.ones(len(tt)),tt/1200]);coef=np.linalg.lstsq(T,xy,rcond=None)[0]
  for _ in range(5):
   norm=np.linalg.norm(xy-T@coef,axis=1);s=max(10.,float(np.median(norm)));weight=np.minimum(1.,1.5*s/np.maximum(norm,1e-15));w=np.sqrt(weight);coef=np.linalg.lstsq(T*w[:,None],xy*w[:,None],rcond=None)[0]
  cart=np.r_[coef[0],coef[1]/1200];h=polar_state(cart,P[0,0]);diag.update(raw_cartesian=cart.tolist(),raw_h=h.tolist(),clipped=bool(np.any((h<LOW)|(h>HIGH))),fallback=not bool(np.isfinite(h).all()),min_weight=float(weight.min()))
  return np.clip(h,LOW,HIGH) if np.isfinite(h).all() else CENTER.copy(),diag
 except Exception as e:diag['reason']=repr(e);return CENTER.copy(),diag

def starts(t,P,y,fixed):
 tic=time.perf_counter();main,md=initialize_main(t,P[:,:1],y[:,:1]);resource,rd=(CENTER.copy(),{'kind':'PHYSICAL_CENTER','fallback':False}) if P.shape[1]==1 else initialize_aux(t,P,y)
 h=np.vstack([main,resource,fixed]);u=np.clip((h-LOW)/SPAN,1e-12,1-1e-12)
 return u,{'main':md,'resource':rd,'initialization_s':time.perf_counter()-tic,'duplicate_starts':128-len(np.unique(u,axis=0))}

class Model:
 def __init__(self,t,P,y,guard=None):
  self.t=np.asarray(t,dtype=float);self.P=np.ascontiguousarray(P,dtype=float);self.y=np.ascontiguousarray(y,dtype=float);self.origin=self.P[0,0].copy();self.guard=guard
  if self.P.shape!=(len(self.t),self.y.shape[1],2) or self.y.shape[0]!=len(self.t) or self.y.shape[1] not in [1,2] or not np.isfinite(self.P).all() or not np.isfinite(self.y).all():raise ValueError('INPUT_INVALID')
  self.key=None;self.f=None;self.J=None;self.fun_calls=0;self.jac_calls=0;self.geometry_calls=0
 def geometry(self,u):
  if self.guard is not None:self.guard()
  if self.key is not None and np.array_equal(u,self.key):return
  h=LOW+SPAN*u;r,th,v,ps=h;th,ps=np.deg2rad([th,ps]);ct,st,cp,sp=np.cos(th),np.sin(th),np.cos(ps),np.sin(ps)
  x=self.origin[0]+1000*r*ct+self.t*v*cp;y=self.origin[1]+1000*r*st+self.t*v*sp
  dx=x[:,None]-self.P[:,:,0];dy=y[:,None]-self.P[:,:,1];den=dx*dx+dy*dy
  if np.any(den<=0):raise ValueError('target at reported receiver')
  self.f=(wrap(np.arctan2(dy,dx)-self.y)/SIGMA).ravel()
  bx=np.column_stack([np.full(len(x),1000*ct),np.full(len(x),-1000*r*st*np.pi/180),self.t*cp,-self.t*v*sp*np.pi/180])
  by=np.column_stack([np.full(len(x),1000*st),np.full(len(x),1000*r*ct*np.pi/180),self.t*sp,self.t*v*cp*np.pi/180])
  self.J=((dx[:,:,None]*by[:,None,:]-dy[:,:,None]*bx[:,None,:])/den[:,:,None]*SPAN/SIGMA).reshape(-1,4)
  self.key=np.array(u,copy=True);self.geometry_calls+=1
 def residual(self,u):self.fun_calls+=1;self.geometry(u);return self.f.copy()
 def jac(self,u):self.jac_calls+=1;self.geometry(u);return self.J.copy()

def fit_one(t,P,y,u,guard=None):
 tic=time.perf_counter();m=Model(t,P,y,guard)
 fit=least_squares(m.residual,u,jac=m.jac,bounds=(np.zeros(4),np.ones(4)),method='trf',loss='soft_l1',f_scale=1.5,x_scale=1.,tr_solver='exact',max_nfev=100,ftol=1e-10,xtol=1e-10,gtol=1e-10)
 # Direct registered scoring is separate from residual evaluation budget.
 m.geometry(fit.x);Q=score(m.f);h=LOW+SPAN*fit.x
 boundary=bool(np.any((fit.x<=1e-8)|(fit.x>=1-1e-8)))
 return np.r_[h,Q,fit.cost,int(fit.success),fit.status,fit.nfev,fit.njev,m.fun_calls,m.jac_calls,m.geometry_calls,fit.optimality,int(boundary),time.perf_counter()-tic]

def accepted(a,threshold):return (a[:,6]==1)&np.isfinite(a[:,:6]).all(axis=1)&(a[:,4]<=threshold)&np.all((a[:,:4]>=LOW)&(a[:,:4]<=HIGH),axis=1)

def choose(a,threshold,shared=False):
 good=accepted(a,threshold)
 if shared:good[1]=False
 ids=np.flatnonzero(good);sel=int(sorted(ids,key=lambda i:(float(a[i,4]),int(i)))[0]) if len(ids) else -1
 reason='ACCEPTED' if sel>=0 else ('ABSTAIN_RESIDUAL' if np.any((a[:,6]==1)&np.isfinite(a[:,4])) else 'NUMERICAL_OR_OPTIMIZER_FAILURE')
 return sel,ids,reason

def envelope(a,ids):
 reps=[];labels=[]
 for i in sorted(ids,key=lambda i:(float(a[i,4]),int(i))):
  h=a[i,:4];cluster=-1
  for k,j in enumerate(reps):
   d=np.abs(h-a[j,:4]);d[[1,3]]=np.abs(np.rad2deg(wrap(np.deg2rad(h[[1,3]]-a[j,[1,3]]))))
   if np.max(d/np.array([.05,.01,.02,.2]))<=1:cluster=k;break
  if cluster<0:cluster=len(reps);reps.append(int(i))
  labels.append([int(i),cluster])
 return {'representative_start_ids':reps,'assignment':labels,'accepted_count':len(ids),'cluster_count':len(reps),'minimum':a[ids,:4].min(axis=0).tolist() if len(ids) else None,'maximum':a[ids,:4].max(axis=0).tolist() if len(ids) else None,'scope':'FINITE_TERMINAL_ENVELOPE_ONLY'}