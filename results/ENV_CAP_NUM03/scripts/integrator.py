"""NUM03 stable real Snell integrals. No acoustic executable calls."""
from pathlib import Path
import math,json,csv,hashlib
import numpy as np
from scipy.interpolate import PchipInterpolator,PPoly
from scipy.optimize import brentq
from scipy.integrate import quad,solve_ivp
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_NUM03');OLD=R/'results/ENV_CAP_FULL_CAMPAIGN';PRIOR=R/'results/ENV_CAP_CZROB02'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,a):Path(p).write_text(json.dumps(a,ensure_ascii=False,indent=2,allow_nan=False,default=lambda x:x.item() if hasattr(x,'item') else str(x))+'\n',encoding='utf-8',newline='\n')
def tab(p,rows):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,list(dict.fromkeys(k for a in rows for k in a)));w.writeheader();w.writerows(rows)
GL={n:np.polynomial.legendre.leggauss(n) for n in [32,64]}
class Profile:
 def __init__(self,z,c,kind):
  self.z=np.asarray(z,float);self.c=np.asarray(c,float);self.kind=kind;self.H=float(self.z[-1]);assert self.z[0]==0 and np.all(np.diff(self.z)>0) and np.isfinite(self.c).all()
  self.f=PchipInterpolator(self.z,self.c,extrapolate=False) if kind=='PCHIP' else PPoly(np.vstack([np.diff(self.c)/np.diff(self.z),self.c[:-1]]),self.z,extrapolate=False)
  self.df=self.f.derivative()
 def turn(self,zs,cp):
  ids=np.flatnonzero((self.z[1:]>zs)&(self.c[:-1]<cp)&(self.c[1:]>=cp)&(np.diff(self.c)>0))
  if not len(ids):return self.H,True
  j=int(ids[0]);a=max(zs,self.z[j])-self.z[j];b=self.z[j+1]-self.z[j];coef=self.f.c[:,j].copy();coef[-1]-=cp
  root=brentq(lambda x:float(np.polyval(coef,x)),a,b,xtol=1e-12,rtol=1e-14)
  return float(self.z[j]+root),False
 def linear_integral(self,zs,end,cp,turned):
  zz=np.r_[zs,self.z[(self.z>zs)&(self.z<end)],end];cc=self.f(zz);cc[-1]=cp if turned else cc[-1]
  da=np.maximum(0,cp-cc);s=np.sqrt(da*(cp+cc));den=s[:-1]+s[1:]
  if not np.all(den>0):raise ArithmeticError('critical C-linear segment')
  return math.fsum(((cc[:-1]+cc[1:])*np.diff(zz)/den).tolist())
 def kernels(self,zs,end,cp,turned):
  cuts=np.r_[zs,self.z[(self.z>zs)&(self.z<end)],end];out=[]
  for k,(a,b) in enumerate(zip(cuts[:-1],cuts[1:])):
   j=int(np.searchsorted(self.z,(a+b)/2,side='right')-1);coef=self.f.c[:,j];aa=float(self.f(a));bb=float(self.f(b));actual_turn=turned and k==len(cuts)-2
   if actual_turn:ref=end;cv=cp;delta=0.
   elif aa>=bb:ref=float(a);cv=aa;delta=cp-cv
   else:ref=float(b);cv=bb;delta=cp-cv
   x=ref-self.z[j];d1=float(np.polyval(np.polyder(coef),x));d2=float(np.polyval(np.polyder(coef,2),x)) if len(coef)>2 else 0.;d3=float(np.polyval(np.polyder(coef,3),x)) if len(coef)>3 else 0.
   if actual_turn and d1<=1e-12:raise ArithmeticError('nonsimple or grazing turn, no finite simple-turn reference')
   if delta < -1e-10:raise ArithmeticError('forbidden segment')
   ulo=math.sqrt(max(0,end-b));uhi=math.sqrt(max(0,end-a));shift=end-ref
   def fun(u,delta=delta,d1=d1,d2=d2,d3=d3,shift=shift,actual_turn=actual_turn):
    u=np.asarray(u);hh=shift-u*u
    if actual_turn:
     factor=d1-.5*d2*u*u+d3/6*u**4;dc=u*u*factor;ccv=cp-dc;v=2*ccv/np.sqrt(factor*(cp+ccv))
    else:
     dc=delta-hh*(d1+.5*d2*hh+d3/6*hh*hh);ccv=cp-dc
     if np.any(dc<=0):raise ArithmeticError('nonpositive stable endpoint difference')
     v=2*u*ccv/np.sqrt(dc*(cp+ccv))
    return v
   points=[]
   if ulo==0 and not actual_turn and delta>0 and abs(d1)>1e-12:
    scale=math.sqrt(delta/abs(d1))
    if ulo<scale<uhi:points=[scale]
   out.append((ulo,uhi,fun,points))
  return out
 def integrate(self,zs,end,cp,turned):
  if self.kind=='C_LINEAR':return self.linear_integral(zs,end,cp,turned),0.,'STABLE_ANALYTIC',0.
  kernels=self.kernels(zs,end,cp,turned);values=[]
  for n in [32,64]:
   xx,ww=GL[n];values.append(math.fsum(float(np.dot(ww,fn((a+b)/2+(b-a)/2*xx))*(b-a)/2) for a,b,fn,points in kernels))
  diff=abs(values[1]-values[0])*2
  if diff<=1e-4:return values[1],2*diff,'GLOBAL_SQRT_G64',diff
  vals=[];errs=[]
  for a,b,fn,points in kernels:
   v,e=quad(lambda u:float(fn(u)),a,b,epsabs=1e-7/max(1,len(kernels)),epsrel=2e-12,points=points or None,limit=150);vals.append(v);errs.append(e)
  err=2*math.fsum(errs)
  if err>.001:raise ArithmeticError('adaptive error estimate exceeds 1 mm')
  return math.fsum(vals),err,'GLOBAL_SQRT_ADAPTIVE',diff
 def calc(self,zs,angle):
  cs=float(self.f(zs));cp=cs/math.cos(math.radians(angle));end,bottom=self.turn(zs,cp);leg,err,how,diff=self.integrate(zs,end,cp,not bottom)
  return {'angle_deg':float(angle),'range_m':2*leg,'turn_z_m':end,'bottom_encountered':bottom,'zero_boundary_deep_return':bool(not bottom and end<self.H-1),'error_estimate_m':err,'algorithm':how,'G32_G64_difference_m':diff,'phase_speed_mps':cp,'nearest_node_m':float(np.min(abs(self.z-end))),'status':'FINITE'}
 def ode(self,zs,angle,step):
  def fun(x,y):
   z=float(np.clip(y[0],0,self.H));return [math.tan(y[1]),-float(self.df(z))/float(self.f(z))]
  def back(x,y):return y[0]-zs
  back.direction=-1;back.terminal=True
  def bed(x,y):return self.H-y[0]
  bed.direction=-1;bed.terminal=True
  def apex(x,y):return y[1]
  apex.direction=-1;apex.terminal=False
  s=solve_ivp(fun,(0.,120000.),[zs,math.radians(angle)],method='DOP853',rtol=1e-11,atol=[1e-6,1e-12],max_step=step,events=(back,bed,apex))
  pp=np.cos(s.y[1])/self.f(np.clip(s.y[0],0,self.H));p0=math.cos(math.radians(angle))/float(self.f(zs));iv=float(np.max(abs(pp/p0-1)))
  return {'success':bool(s.success),'returned':bool(len(s.t_events[0])),'range_m':float(s.t_events[0][0]) if len(s.t_events[0]) else None,'turn_z_m':float(s.y_events[2][0,0]) if len(s.y_events[2]) else None,'invariant_error':iv,'step_max_m':step,'rhs_evaluations':s.nfev},s.t,s.y
