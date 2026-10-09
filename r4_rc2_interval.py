"""Directed binary64 interval geometry, with rational Machin/Taylor remainder certificates."""
from fractions import Fraction
from functools import lru_cache
from math import factorial
import numpy as np
DN=lambda x:np.nextafter(x,-np.inf)
UP=lambda x:np.nextafter(x,np.inf)
def rational(q):
    f=float(q)
    lo=f if Fraction.from_float(f)<=q else np.nextafter(f,-np.inf)
    hi=f if Fraction.from_float(f)>=q else np.nextafter(f,np.inf)
    return float(lo),float(hi)
def atan_rational(x,n=35):
    s=sum(((-1)**i*x**(2*i+1)/Fraction(2*i+1) for i in range(n)),Fraction())
    t=(-1)**n*x**(2*n+1)/Fraction(2*n+1)
    return min(s,s+t),max(s,s+t)
a,b=atan_rational(Fraction(1,5));c,d=atan_rational(Fraction(1,239))
PI_EXACT_LO=16*a-4*d;PI_EXACT_HI=16*b-4*c
PI=(rational(PI_EXACT_LO)[0],rational(PI_EXACT_HI)[1])
RAD=(rational(PI_EXACT_LO/180)[0],rational(PI_EXACT_HI/180)[1])
def add(a,b):return DN(a[0]+b[0]),UP(a[1]+b[1])
def sub(a,b):return DN(a[0]-b[1]),UP(a[1]-b[0])
def mul(a,b):
    xs=np.asarray([a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]])
    return DN(xs.min(axis=0)),UP(xs.max(axis=0))
def div_positive(a,b):
    assert np.all(b[0]>0)
    xs=np.asarray([a[0]/b[0],a[0]/b[1],a[1]/b[0],a[1]/b[1]])
    return DN(xs.min(axis=0)),UP(xs.max(axis=0))
SIN=[rational(Fraction((-1)**i,factorial(2*i+1))) for i in range(7)]
COS=[rational(Fraction((-1)**i,factorial(2*i))) for i in range(7)]
ATAN=[rational(Fraction((-1)**i,2*i+1)) for i in range(9)]
REM_S=rational(Fraction(1,3)**15/Fraction(factorial(15)))[1]
REM_C=rational(Fraction(1,3)**14/Fraction(factorial(14)))[1]
REM_A=rational(Fraction(1,3)**19/Fraction(19))[1]
def polynomial(x,coeff,odd,rem):
    sq=mul(x,x);p=coeff[-1]
    for coef in reversed(coeff[:-1]):p=add(mul(p,sq),coef)
    if odd:p=mul(p,x)
    return DN(p[0]-rem),UP(p[1]+rem)
def sinpoint(deg):return polynomial(mul((deg,deg),RAD),SIN,True,REM_S)
def cospoint(deg):return polynomial(mul((deg,deg),RAD),COS,False,REM_C)
@lru_cache(maxsize=65536)
def trig(lo,hi):
    assert -15<=lo<=hi<=15
    sl=sinpoint(lo);sh=sinpoint(hi);cl=cospoint(lo);ch=cospoint(hi)
    return (sl[0],sh[1]),(min(cl[0],ch[0]),1. if lo<=0<=hi else max(cl[1],ch[1]))
def predict(box,t,platform):
    # box columns [r_km,theta_deg,v_mps,psi_deg]; all coordinates remain intervals.
    sn,cs=trig(*box[1]);sv,cv=trig(*box[3])
    rr=mul(tuple(box[0]),(1000.,1000.))
    x=sub(add(mul(rr,cs),mul(mul(tuple(box[2]),cv),(t,t))),(platform[:,0],platform[:,0]))
    y=sub(add(mul(rr,sn),mul(mul(tuple(box[2]),sv),(t,t))),(platform[:,1],platform[:,1]))
    ratio=div_positive(y,x)
    assert np.all(abs(ratio[0])<1/3) and np.all(abs(ratio[1])<1/3)
    low=polynomial((ratio[0],ratio[0]),ATAN,True,REM_A)[0]
    high=polynomial((ratio[1],ratio[1]),ATAN,True,REM_A)[1]
    return low,high
def distance(lo,hi,observations):
    # Physical +x-sector and archived observations stay inside(-1/3,1/3)rad.
    # No +/-pi seam occurs here. Optional reflection is a conservative added ambiguity.
    assert np.max(abs(observations))<1/3
    ds=[]
    for b in [observations,-observations]:
        below=DN(lo-b);above=DN(b-hi)
        ds.append(np.maximum(0,np.maximum(below,above)))
    return np.minimum(ds[0],ds[1])
def bound(box,t,platform,obs,sigma,slack=1e-10):
    lo,hi=predict(box,t,platform);dist=distance(lo,hi,obs)
    dist=np.maximum(0,DN(dist-slack));z=DN(dist/sigma)
    ss=np.maximum(0,DN(z*z));acc=0.
    for value in ss:acc=max(0.,float(DN(acc+value)))
    return acc,float(z.max()),lo,hi
