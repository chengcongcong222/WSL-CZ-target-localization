"""Independent references, never imports integrator: own PCHIP slopes, local adaptive mapping, arbitrary precision."""
import math
import numpy as np
from scipy.integrate import quad
import mpmath as mp
class Reference:
 def __init__(self,z,c,kind):
  self.z=list(map(float,z));self.c=list(map(float,c));self.kind=kind;n=len(z);h=np.diff(z);v=np.diff(c)/h;d=np.zeros(n)
  if kind=='PCHIP':
   for j in range(1,n-1):
    if v[j-1]*v[j]>0:
     a=2*h[j]+h[j-1];b=h[j]+2*h[j-1];d[j]=(a+b)/(a/v[j-1]+b/v[j])
   for j,hh,ii,vv,vi in [(0,h[0],h[1],v[0],v[1]),(n-1,h[-1],h[-2],v[-1],v[-2])]:
    x=((2*hh+ii)*vv-hh*vi)/(hh+ii);d[j]=0 if x*vv<=0 else (3*vv if vv*vi<0 and abs(x)>3*abs(vv) else x)
  self.coeff=[]
  for j in range(n-1):
   self.coeff.append(((d[j]+d[j+1]-2*v[j])/h[j]**2,(3*v[j]-2*d[j]-d[j+1])/h[j],d[j],c[j]) if kind=='PCHIP' else (0.,0.,v[j],c[j]))
 def index(self,x):return min(len(self.z)-2,max(0,int(np.searchsorted(self.z,x,side='right')-1)))
 def value(self,x):
  j=self.index(x);a,b,c,d=self.coeff[j];s=x-self.z[j];return ((a*s+b)*s+c)*s+d
 def endpoint(self,zs,cp):
  for j in range(len(self.z)-1):
   if self.z[j+1]>zs and self.c[j]<cp<=self.c[j+1] and self.c[j+1]>self.c[j]:
    a,b,c,d=self.coeff[j];lo=max(zs,self.z[j])-self.z[j];hi=self.z[j+1]-self.z[j]
    for k in range(70):
     mid=(lo+hi)/2;val=((a*mid+b)*mid+c)*mid+(d-cp)
     if val<0:lo=mid
     else:hi=mid
    return self.z[j]+(lo+hi)/2,False
  return self.z[-1],True
 def calc(self,zs,angle):
  cp=self.value(zs)/math.cos(math.radians(angle));end,bb=self.endpoint(zs,cp);cuts=[zs]+[z for z in self.z if zs<z<end]+[end];vals=[];errors=[]
  for k,(lo,hi) in enumerate(zip(cuts[:-1],cuts[1:])):
   j=self.index((lo+hi)/2);a,b,c,d=self.coeff[j];cl=self.value(lo);ch=self.value(hi);turned=not bb and k==len(cuts)-2
   anchor=hi if turned or ch>=cl else lo;sgn=-1 if anchor==hi else 1;cv=cp if turned else max(cl,ch);delta=0. if turned else cp-cv;x=anchor-self.z[j];d1=(3*a*x+2*b)*x+c;d2=6*a*x+2*b;d3=6*a
   if turned and d1<=1e-12:raise ArithmeticError('nonsimple turn')
   def fun(v):
    hh=sgn*v*v;factor=-sgn*d1-.5*d2*v*v-sgn*d3/6*v**4
    if turned:
     cc=cp-v*v*factor;return 2*cc/math.sqrt(factor*(cp+cc))
    dd=delta+v*v*factor;cc=cp-dd;return 0. if v==0 else 2*v*cc/math.sqrt(dd*(cp+cc))
   vmax=math.sqrt(hi-lo);points=[]
   if delta>0 and abs(d1)>1e-12:
    vs=math.sqrt(delta/abs(d1))
    if 0<vs<vmax:points=[vs]
   value,error=quad(fun,0,vmax,epsabs=1e-7/max(1,len(cuts)),epsrel=2e-12,points=points or None,limit=150);vals.append(value);errors.append(error)
  return {'range_m':2*math.fsum(vals),'turn_z_m':end,'bottom_encountered':bb,'error_estimate_m':2*math.fsum(errors)}

def mp_reference(z,c,zs,angle,kind,dps=60):
 """Rebuild slopes/roots in arbitrary precision; tanh-sinh quadrature, no float64 coefficient lifting."""
 with mp.workdps(dps):
  z=[mp.mpf(float(v)) for v in z];c=[mp.mpf(float(v)) for v in c];zs=mp.mpf(float(zs));angle=mp.mpf(float(angle));n=len(z);h=[z[j+1]-z[j] for j in range(n-1)];v=[(c[j+1]-c[j])/h[j] for j in range(n-1)];s=[mp.mpf(0)]*n
  if kind=='PCHIP':
   for j in range(1,n-1):
    if v[j-1]*v[j]>0:
     w1=2*h[j]+h[j-1];w2=h[j]+2*h[j-1];s[j]=(w1+w2)/(w1/v[j-1]+w2/v[j])
   for j,h0,h1,d0,d1 in [(0,h[0],h[1],v[0],v[1]),(n-1,h[-1],h[-2],v[-1],v[-2])]:
    x=((2*h0+h1)*d0-h0*d1)/(h0+h1);s[j]=0 if x*d0<=0 else (3*d0 if d0*d1<0 and abs(x)>3*abs(d0) else x)
  coeff=[((s[j]+s[j+1]-2*v[j])/h[j]**2,(3*v[j]-2*s[j]-s[j+1])/h[j],s[j],c[j]) if kind=='PCHIP' else (mp.mpf(0),mp.mpf(0),v[j],c[j]) for j in range(n-1)]
  def index(x):return min(n-2,max(0,sum(zz<=x for zz in z)-1))
  def value(x):
   j=index(x);a,b,d,e=coeff[j];u=x-z[j];return ((a*u+b)*u+d)*u+e
  cp=value(zs)/mp.cos(angle*mp.pi/180);end=z[-1];bb=True
  for j in range(n-1):
   if z[j+1]>zs and c[j]<cp<=c[j+1] and c[j+1]>c[j]:
    lo=max(zs,z[j]);hi=z[j+1]
    for k in range(4*dps):
     mid=(lo+hi)/2
     if value(mid)<cp:lo=mid
     else:hi=mid
    end=(lo+hi)/2;bb=False;break
  cuts=[zs]+[zz for zz in z if zs<zz<end]+[end];terms=[];errsum=mp.mpf(0)
  if kind=='C_LINEAR':
   for a,b in zip(cuts[:-1],cuts[1:]):
    ca=value(a);cb=cp if not bb and b==end else value(b);sa=mp.sqrt((cp-ca)*(cp+ca));sb=mp.sqrt(max(mp.mpf(0),(cp-cb)*(cp+cb)));terms.append((ca+cb)*(b-a)/(sa+sb))
  else:
   for k,(lo,hi) in enumerate(zip(cuts[:-1],cuts[1:])):
    j=index((lo+hi)/2);a,b,d,e=coeff[j];turned=not bb and k==len(cuts)-2;cl=value(lo);ch=value(hi);anchor=hi if turned or ch>=cl else lo;sgn=-1 if anchor==hi else 1;cv=cp if turned else max(cl,ch);delta=0 if turned else cp-cv;x=anchor-z[j];d1=(3*a*x+2*b)*x+d;d2=6*a*x+2*b;d3=6*a
    def fun(u):
     factor=-sgn*d1-d2*u*u/2-sgn*d3*u**4/6
     if turned:
      cc=cp-u*u*factor;return 2*cc/mp.sqrt(factor*(cp+cc))
     dc=delta+u*u*factor;cc=cp-dc;return 0 if u==0 else 2*u*cc/mp.sqrt(dc*(cp+cc))
    vmax=mp.sqrt(hi-lo);points=[mp.mpf(0)]
    if delta>0 and d1!=0:
     sc=mp.sqrt(delta/abs(d1))
     if 0<sc<vmax:points.append(sc)
    points.append(vmax);vv,ee=mp.quad(fun,points,method='tanh-sinh',error=True,maxdegree=7);terms.append(vv);errsum+=ee
  return {'range_m':float(2*mp.fsum(terms)),'turn_z_m':float(end),'bottom_encountered':bool(bb),'quad_error_estimate_m':float(2*errsum),'dps':dps,'range_decimal':mp.nstr(2*mp.fsum(terms),dps-5)}
