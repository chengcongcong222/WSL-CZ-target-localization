"""Public-input-only T3-W1 estimator. No truth/generator/archive imports."""
import math
import numpy as np
from scipy.optimize import least_squares

def wrap(a):return (np.asarray(a)+np.pi)%(2*np.pi)-np.pi

def initialize(t,P,y):
    u=np.stack([np.cos(y),np.sin(y)],axis=-1);d=P[:,1]-P[:,0]
    det=u[:,0,0]*u[:,1,1]-u[:,0,1]*u[:,1,0]
    with np.errstate(divide='ignore',invalid='ignore'):
        a=(d[:,0]*u[:,1,1]-d[:,1]*u[:,1,0])/det
        b=(d[:,0]*u[:,0,1]-d[:,1]*u[:,0,0])/det
        xy=P[:,0]+a[:,None]*u[:,0]
    good=(abs(det)>1e-12)&(a>0)&(b>0)&np.isfinite(xy).all(axis=1)
    if np.sum(good)<4:return np.array([52.5,0.,2.,0.]),True,int(np.sum(good))
    T=np.column_stack([np.ones(np.sum(good)),t[good]/1200]);xy=xy[good];c=np.linalg.lstsq(T,xy,rcond=1e-12)[0]
    for _ in range(5):
        norm=np.linalg.norm(xy-T@c,axis=1);sigma=max(10.,float(np.median(norm)))
        w=np.sqrt(np.minimum(1.,1.5*sigma/np.maximum(norm,1e-15)))
        c=np.linalg.lstsq(T*w[:,None],xy*w[:,None],rcond=1e-12)[0]
    dp=c[0]-P[0,0];vel=c[1]/1200
    initial=np.array([np.hypot(*dp)/1000,np.rad2deg(np.arctan2(dp[1],dp[0])),np.hypot(*vel),np.rad2deg(np.arctan2(vel[1],vel[0]))])
    return initial,False,int(np.sum(good))

class Model:
    def __init__(self,t,P,y,quality):
        self.t=np.asarray(t);self.P=np.asarray(P);self.y=np.asarray(y);self.origin=self.P[0,0].copy();self.sigma=np.deg2rad(quality['bearing_sigma_deg']);self.nav=quality['nav_sigma_m'];c,d=quality['common_bound_deg'],quality['half_diff_bound_deg']
        self.low=np.array([45.,-5.,1.,-15.,-c,-d]);self.span=np.array([15.,10.,2.,30.,2*c,2*d]);self.key=None
    def evaluate(self,q):
        if self.key is not None and np.array_equal(q,self.key):return self.f,self.J
        r,theta,v,psi,c,d=self.low+self.span*q;theta,psi=np.deg2rad([theta,psi]);co,si,cp,sp=np.cos(theta),np.sin(theta),np.cos(psi),np.sin(psi)
        x=self.origin[0]+1000*r*co+self.t*v*cp;y=self.origin[1]+1000*r*si+self.t*v*sp
        D=np.stack([x[:,None]-self.P[:,:,0],y[:,None]-self.P[:,:,1]],axis=-1);rho2=np.sum(D*D,axis=-1)
        if np.any(rho2<=1):raise FloatingPointError('TARGET_AT_RECEIVER')
        dx=np.column_stack([np.full(len(x),1000*co),np.full(len(x),-1000*r*si*np.pi/180),self.t*cp,-self.t*v*sp*np.pi/180,np.zeros(len(x)),np.zeros(len(x))])*self.span
        dy=np.column_stack([np.full(len(x),1000*si),np.full(len(x),1000*r*co*np.pi/180),self.t*sp,self.t*v*cp*np.pi/180,np.zeros(len(x)),np.zeros(len(x))])*self.span
        pred=np.arctan2(D[:,:,1],D[:,:,0])+np.deg2rad(c+d*np.array([-1.,1.]))
        e=wrap(pred-self.y);variance=self.sigma**2+self.nav**2/rho2;std=np.sqrt(variance)
        Jmean=(D[:,:,0,None]*dy[:,None,:]-D[:,:,1,None]*dx[:,None,:])/rho2[:,:,None]
        Jmean[:,:,4]=np.deg2rad(self.span[4]);Jmean[:,:,5]=np.array([-1,1])*np.deg2rad(self.span[5])
        radial_dot=D[:,:,0,None]*dx[:,None,:]+D[:,:,1,None]*dy[:,None,:]
        dlog=-self.nav**2*radial_dot/(rho2[:,:,None]**2*variance[:,:,None])
        r1=e/std;r2=np.sqrt(np.log1p(self.nav**2/(rho2*self.sigma**2)))
        j1=Jmean/std[:,:,None]-r1[:,:,None]*dlog;j2=dlog/r2[:,:,None]
        self.f=np.r_[r1.ravel(),r2.ravel()];self.J=np.vstack([j1.reshape(-1,6),j2.reshape(-1,6)]);self.key=np.array(q,copy=True)
        return self.f,self.J
    def fun(self,q):return self.evaluate(q)[0].copy()
    def jac(self,q):return self.evaluate(q)[1].copy()

def estimate(t,P,y,quality):
    model=Model(t,P,y,quality);initial,fallback,count=initialize(t,P,y)
    h0=np.r_[initial,0.,0.];q0=np.clip((h0-model.low)/model.span,1e-10,1-1e-10)
    fit=least_squares(model.fun,q0,jac=model.jac,bounds=(np.zeros(6),np.ones(6)),method='trf',loss='linear',x_scale=1.,tr_solver='exact',max_nfev=160,ftol=1e-10,xtol=1e-10,gtol=1e-10)
    h=model.low+model.span*fit.x;singular=np.linalg.svd(fit.jac,compute_uv=False);ratio=float(singular[-1]/singular[0]);rank=int(np.sum(singular/singular[0]>1e-10));near=(fit.x<=1e-6)|(fit.x>=1-1e-6)
    return h,{'accepted':bool(fit.success and np.isfinite(h).all() and np.isfinite(fit.fun).all()),'solver_status':int(fit.status),'nfev':int(fit.nfev),'cost':float(fit.cost),'optimality':float(fit.optimality),'initializer_fallback':bool(fallback),'initializer_intersections':count,'initial_h':h0.tolist(),'state_at_bound':bool(np.any(near[:4])),'bias_at_bound':bool(np.any(near[4:])),'common_at_bound':bool(near[4]),'diff_at_bound':bool(near[5]),'rank_augmented_likelihood_jac':rank,'smallest_largest_jac_ratio':ratio,'mean_residual_RMS_standard':float(np.sqrt(np.mean(fit.fun[:len(t)*2]**2)))}
