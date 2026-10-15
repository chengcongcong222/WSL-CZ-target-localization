"""Primary continuous reference: segmented DOP853, two-sided shooting, complex continuation.
No historical eigensolver or second reference is imported.
"""
import math
import numpy as np
from scipy.integrate import solve_ivp

S=2*math.pi*201/1500

def potential_data(depth,speed,loss):
    x=(201/1000)**2
    alpha=(.0033+.11*x/(1+x)+44*x/(4100+x)+.0003*x)/8685.8896
    return np.asarray(speed,dtype=complex)+1j*loss*alpha*np.asarray(speed)**2/(2*math.pi*201)

class Shoot:
    def __init__(self,z,c,match,rtol,atol,max_step):
        self.z=np.asarray(z,float);self.c=np.asarray(c,complex);self.match=float(match)
        self.rtol=rtol;self.atol=atol;self.step=max_step;self.last=None
    def side(self,k,start,end,initial,record=False):
        points=self.z[(self.z>min(start,end))&(self.z<max(start,end))]
        if start>end:points=points[::-1]
        edges=np.r_[start,points,end];v=np.array(initial,complex);log=0.;segments=[]
        for a,b in zip(edges[:-1],edges[1:]):
            ca=np.interp(a,self.z,self.c.real)+1j*np.interp(a,self.z,self.c.imag)
            cb=np.interp(b,self.z,self.c.real)+1j*np.interp(b,self.z,self.c.imag)
            def f(t,y):
                cc=ca+(cb-ca)*(t-a)/(b-a)
                return np.array([S*y[1],(k*k-(2*math.pi*201/cc)**2)*y[0]/S])
            sol=solve_ivp(f,(a,b),v,method='DOP853',rtol=self.rtol,atol=self.atol,max_step=self.step,dense_output=record)
            if not sol.success or not np.isfinite(sol.y[:,-1]).all():raise RuntimeError('DOP853_INTEGRATION_FAILED')
            if record:segments.append((float(a),float(b),sol.sol,log))
            w=sol.y[:,-1];amp=float(np.linalg.norm(w))
            if not math.isfinite(amp) or amp==0:raise RuntimeError('SHOOT_CONDITION_FAILED')
            log+=math.log(amp);v=w/amp
        return v,log,segments
    def mismatch(self,k):
        left=self.side(k,0.,self.match,[0,1])[0];right=self.side(k,5000.,self.match,[1,0])[0]
        return left[0]*right[1]-left[1]*right[0]
    def root(self,guess):
        k=complex(guess);history=[]
        for i in range(10):
            value=self.mismatch(k);eps=2e-7
            derivative=(self.mismatch(k+eps)-self.mismatch(k-eps))/(2*eps)
            if abs(derivative)<1e-12:raise RuntimeError('SHOOT_SINGULAR_ROOT')
            change=value/derivative
            history.append([k.real,k.imag,float(abs(value)),float(abs(change))])
            if abs(change)<1e-13 and abs(value)<=2e-8:
                return k,history
            if abs(change)>2e-6:raise RuntimeError('ROOT_LEFT_FROZEN_IDENTITY_NEIGHBORHOOD')
            k-=change
        raise RuntimeError('SHOOT_ROOT_NOT_CONVERGED')
    def function(self,k):
        l,ll,ls=self.side(k,0.,self.match,[0,1],True)
        r,lr,rs=self.side(k,5000.,self.match,[1,0],True)
        factor=np.vdot(r,l)/np.vdot(r,r)
        jump=float(np.linalg.norm(l-factor*r)/np.linalg.norm(l))
        def raw(depth):
            zz=np.atleast_1d(np.asarray(depth,float));out=np.zeros((2,len(zz)),complex)
            for a,b,fn,log in ls:
                ix=(zz>=min(a,b))&(zz<=max(a,b))&(zz<=self.match)
                if ix.any():out[:,ix]=fn(zz[ix])*math.exp(log-ll)
            for a,b,fn,log in rs:
                ix=(zz>=min(a,b))&(zz<=max(a,b))&(zz>self.match)
                if ix.any():out[:,ix]=fn(zz[ix])*math.exp(log-lr)*factor
            return out
        gx,gw=np.polynomial.legendre.leggauss(32);integral=0j
        aa=self.z[:-1,None];bb=self.z[1:,None]
        zz=aa+(gx+1)*(bb-aa)/2;v=raw(zz.ravel())[0].reshape(zz.shape)
        integral=np.sum(gw*v*v*(bb-aa)/2)
        if abs(integral)<1e-20:raise RuntimeError('BILINEAR_NORMALIZATION_SINGULAR')
        norm=np.sqrt(integral)
        def normalized(z):return raw(z)/norm
        self.last={'match_jump':jump,'pre_norm_real':float(integral.real),'pre_norm_imag':float(integral.imag),'normalization':'BILINEAR_INT_PHI_SQUARED'}
        return normalized,self.last

def solve(z,speed,mode,initial,level,losses,observer=None):
    opt=[(2e-10,2e-12,2.),(2e-12,2e-14,1.)][level]
    match=1300. if mode==1 else 2500.
    current=complex(initial);roots=[];last=None
    for loss in losses:
        cp=potential_data(z,speed,loss);solver=Shoot(z,cp,match,*opt)
        current,h=solver.root(current)
        if abs(current.real-complex(initial).real)>3e-7 or not -1e-5<=current.imag<=1e-12:raise RuntimeError('ROOT_OUTSIDE_FROZEN_DOMAIN')
        fn,info=solver.function(current)
        if loss==0:
            pp=fn(np.arange(0.,5000.25,.25))[0].real;vv=pp[abs(pp)>1e-12*max(abs(pp))];info['sturm_internal_nodes']=int(np.sum(vv[:-1]*vv[1:]<0))
        if observer is not None:observer(loss,current,fn,info)
        roots.append({'loss':loss,'k_real':current.real,'k_imag':current.imag,'root_trace':h,**info})
        last=(solver,fn,info)
    return current,last[1],roots,last[2]
