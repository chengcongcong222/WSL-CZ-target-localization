"""Independent high-order Gauss-Legendre collocation transfer BVP.
No DOP853, historical FE/FD solver, or primary roots/functions imported.
"""
import math
import numpy as np
from scipy.special import eval_legendre

S=math.tau*201/1500

def medium(depth,speed,scale):
    f2=.201**2
    attenuation=(.0033+.11*f2/(1+f2)+44*f2/(4100+f2)+.0003*f2)/8685.8896
    return np.asarray(speed)+1j*scale*attenuation*np.asarray(speed)**2/(math.tau*201)

class Collocation:
    def __init__(self,z,c,mode,level):
        self.z=np.asarray(z);self.c=np.asarray(c);self.match=1300. if mode==1 else 2500.
        self.order=[8,10][level];self.ci,self.w=np.polynomial.legendre.leggauss(self.order);self.ci=(self.ci+1)/2;self.w=self.w/2
        # Integrated cardinal basis, represented in shifted Legendre polynomials.
        V=np.polynomial.legendre.legvander(2*self.ci-1,self.order-1)
        self.coeff=np.linalg.solve(V,np.eye(self.order))
        self.anti=np.polynomial.legendre.legint(self.coeff,axis=0)/2
        self.offset=np.polynomial.legendre.legval(-1.,self.anti)
        self.A=(np.polynomial.legendre.legval(2*self.ci-1,self.anti).T-self.offset)
        self.id=np.eye(self.order*2,dtype=complex)
    def transfer(self,k,a,b):
        step=b-a;zz=a+self.ci*step
        cc=np.interp(zz,self.z,self.c.real)+1j*np.interp(zz,self.z,self.c.imag)
        v=(math.tau*201/cc)**2
        As=np.zeros((self.order,2,2),complex);As[:,0,1]=S;As[:,1,0]=(k*k-v)/S
        # Stage i: Y_i=y_a+h sum_j A_ij L_j Y_j.
        blocks=-step*self.A[:,:,None,None]*As[None,:,:,:]
        mat=self.id+blocks.transpose(0,2,1,3).reshape(self.order*2,self.order*2)
        rhs=np.tile(np.eye(2),(self.order,1))
        stages=np.linalg.solve(mat,rhs).reshape(self.order,2,2)
        derivative=np.einsum('sij,sjk->sik',As,stages)
        trans=np.eye(2)+step*np.einsum('s,sij->ij',self.w,derivative)
        return trans,derivative
    def side(self,k,start,end,initial,record=False):
        # 5 m steps always divide the 20 m SSP intervals; no step crosses a knot.
        grid=np.linspace(start,end,int(abs(end-start)/5)+1);value=np.array(initial,complex);log=0.;parts=[]
        for a,b in zip(grid[:-1],grid[1:]):
            T,der=self.transfer(k,a,b)
            if record:parts.append((float(a),float(b),value.copy(),der,log))
            nxt=T@value;amp=float(np.linalg.norm(nxt))
            if not math.isfinite(amp) or amp==0:raise RuntimeError('COLLOCATION_CONDITION_FAILED')
            value=nxt/amp;log+=math.log(amp)
        return value,log,parts
    def determinant(self,k):
        u=self.side(k,0.,self.match,[0,1])[0];v=self.side(k,5000.,self.match,[1,0])[0]
        return np.linalg.det(np.column_stack((u,v)))
    def root(self,guess):
        k=complex(guess);trace=[]
        for _ in range(10):
            f=self.determinant(k);h=1e-7;df=(self.determinant(k+h)-self.determinant(k-h))/(2*h)
            if abs(df)<1e-12:raise RuntimeError('COLLOCATION_SINGULAR_ROOT')
            dk=f/df;trace.append([k.real,k.imag,float(abs(f)),float(abs(dk))])
            if abs(dk)<1e-13 and abs(f)<=2e-8:return k,trace
            if abs(dk)>2e-6:raise RuntimeError('COLLOCATION_ROOT_LEFT_FIXED_NEIGHBORHOOD')
            k-=dk
        raise RuntimeError('COLLOCATION_ROOT_NOT_CONVERGED')
    def function(self,k):
        u,lu,left=self.side(k,0.,self.match,[0,1],True);v,lv,right=self.side(k,5000.,self.match,[1,0],True)
        beta=(v.conj()@u)/(v.conj()@v);jump=float(np.linalg.norm(u-beta*v))
        def raw(z):
            zz=np.atleast_1d(np.asarray(z,float));out=np.zeros((2,len(zz)),complex)
            for pieces,totlog,factor,which in [(left,lu,1,0),(right,lv,beta,1)]:
                for a,b,y,der,log in pieces:
                    ix=(zz>=min(a,b))&(zz<=max(a,b))&((zz<=self.match) if which==0 else (zz>self.match))
                    if not ix.any():continue
                    t=(zz[ix]-a)/(b-a)
                    ints=np.polynomial.legendre.legval(2*t-1,self.anti).T-self.offset
                    stage_derivatives=np.einsum('sij,j->si',der,y)
                    yy=y[:,None]+(b-a)*(ints@stage_derivatives).T
                    out[:,ix]=yy*math.exp(log-totlog)*factor
            return out
        x,w=np.polynomial.legendre.leggauss(16);bilinear=0j
        zz=np.arange(0.,5000.,5.)[:,None]+2.5*(x+1)
        pp=raw(zz.ravel())[0].reshape(zz.shape);bilinear=2.5*np.sum(w*pp*pp)
        if abs(bilinear)<1e-20:raise RuntimeError('COLLOCATION_BILINEAR_SINGULAR')
        den=np.sqrt(bilinear)
        return lambda z:raw(z)/den,{'match_jump':jump,'pre_norm_real':float(bilinear.real),'pre_norm_imag':float(bilinear.imag),'normalization':'BILINEAR_INT_PHI_SQUARED'}

def solve(z,speed,mode,initial,level,losses,observer=None):
    k=complex(initial);records=[];last=None
    for loss in losses:
        obj=Collocation(z,medium(z,speed,loss),mode,level);k,trace=obj.root(k)
        if abs(k.real-complex(initial).real)>3e-7 or not -1e-5<=k.imag<=1e-12:raise RuntimeError('COLLOCATION_ROOT_OUTSIDE_FROZEN_DOMAIN')
        fn,info=obj.function(k)
        if loss==0:
            pp=fn(np.arange(0.,5000.25,.25))[0].real;vv=pp[abs(pp)>1e-12*max(abs(pp))];info['sturm_internal_nodes']=int(np.sum(vv[:-1]*vv[1:]<0))
        if observer is not None:observer(loss,k,fn,info)
        records.append({'loss':loss,'k_real':k.real,'k_imag':k.imag,'root_trace':trace,**info});last=fn
    return k,last,records,info
