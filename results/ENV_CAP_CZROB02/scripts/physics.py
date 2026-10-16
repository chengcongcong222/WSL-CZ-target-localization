"""CZROB02 geometric physics. No field/energy surrogate is constructed."""
from pathlib import Path
import json,csv,hashlib,math
import numpy as np
from scipy.interpolate import PchipInterpolator,PPoly
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
O=Path(__file__).resolve().parents[1];R=O.parents[1];OLD=R/'results/ENV_CAP_FULL_CAMPAIGN';REC=R/'results/ENV_CAP_REC01';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZROB02');RAW=W.parent/'ENV_CAP/E1'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False,default=lambda v:v.item() if hasattr(v,'item') else str(v))+'\n',encoding='utf-8',newline='\n')
def tab(p,rows,cols=None):
 keys=cols or list(dict.fromkeys(k for a in rows for k in a))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  d=csv.DictWriter(f,keys);d.writeheader();d.writerows(rows)
NODES={n:np.polynomial.legendre.leggauss(n) for n in [32,64]}
class Profile:
 def __init__(self,z,c,kind):
  self.z=np.asarray(z,float);self.c=np.asarray(c,float);self.kind=kind;self.H=float(z[-1]);assert np.all(np.diff(self.z)>0) and np.isfinite(self.c).all() and np.min(self.c)>0
  if kind=='PCHIP':self.f=PchipInterpolator(self.z,self.c,extrapolate=False)
  elif kind=='C_LINEAR':self.f=PPoly(np.vstack((np.diff(self.c)/np.diff(self.z),self.c[:-1])),self.z,extrapolate=False)
  else:raise ValueError(kind)
  self.df=self.f.derivative()
 def v(self,z):return self.f(z)
 def endpoint(self,zs,cp,down):
  if down:
   idx=np.flatnonzero((self.z[1:]>zs)&(self.c[:-1]<cp)&(self.c[1:]>=cp)&(np.diff(self.c)>0))
   if not len(idx):return self.H,True
   j=int(idx[0]);a=max(zs,float(self.z[j]));b=float(self.z[j+1])
  else:
   idx=np.flatnonzero((self.z[:-1]<zs)&(self.c[:-1]>=cp)&(self.c[1:]<cp)&(np.diff(self.c)<0))
   if not len(idx):return 0.,True
   j=int(idx[-1]);a=float(self.z[j]);b=min(zs,float(self.z[j+1]))
  root=brentq(lambda z:float(self.f(z))-cp,a,b,xtol=1e-9,rtol=1e-14)
  return float(root),False
 def integral(self,a,b,cp,q,turn_end=None):
  """dx= c dz/sqrt(cp^2-c^2); endpoint-square substitution removes simple-turn singularity."""
  assert 0<=a<b<=self.H
  cuts=np.r_[a,self.z[(self.z>a)&(self.z<b)],b];lo=cuts[:-1].copy();hi=cuts[1:].copy();xx,ww=NODES[q];ans=0.
  if turn_end is not None:
   at_top=turn_end=='LOW';k=0 if at_top else len(lo)-1;zt=a if at_top else b;bound=hi[k] if at_top else lo[k];L=abs(bound-zt)
   jj=int(np.searchsorted(self.z,(zt+bound)/2,side='right')-1);coef=self.f.c[:,jj];x=zt-self.z[jj]
   deriv=float(np.polyval(np.polyder(coef),x));d2=float(np.polyval(np.polyder(coef,2),x)) if len(coef)>2 else 0.;d3=float(np.polyval(np.polyder(coef,3),x)) if len(coef)>3 else 0.
   sign=1. if at_top else -1.;g=-sign*deriv
   if g<=1e-12:raise ArithmeticError('GRAZING_OR_NONSIMPLE_TURN')
   u=(xx+1)*math.sqrt(L)/2;dz=sign*u*u;cv=np.asarray(self.f(zt+dz));ratio=g-.5*d2*u*u-sign*d3/6*u**4
   if np.any(ratio<=0):raise ArithmeticError('NONPOSITIVE_TURN_FACTOR')
   val=2*cv/np.sqrt(ratio*(cp+cv));ans=float(np.dot(ww,val)*math.sqrt(L)/2);lo=np.delete(lo,k);hi=np.delete(hi,k)
  if len(lo):
   zz=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*xx[None,:]/2;cv=self.f(zz);den=(cp-cv)*(cp+cv)
   if not np.all(den>0):raise ArithmeticError('INTERNAL_FORBIDDEN_REGION_OR_ROUNDOFF')
   ans+=float(np.sum((hi-lo)/2*np.sum(ww[None,:]*cv/np.sqrt(den),axis=1)))
  return ans
 def evaluate(self,zs,angle,q):
  cp=float(self.f(zs))/math.cos(math.radians(angle));bottom,bb=self.endpoint(zs,cp,True);top,sb=self.endpoint(zs,cp,False)
  deep=2*self.integral(zs,bottom,cp,q,None if bb else 'HIGH')/1000;up=2*self.integral(top,zs,cp,q,None if sb else 'LOW')/1000
  return {'angle_deg':float(angle),'deep_return_km':deep,'deep_turn_z_m':bottom,'deep_bottom_bounce':bb,'deep_refracted_clearance_PASS':bool(not bb and bottom<self.H-1),'upper_return_km':up,'upper_surface_bounce':sb,'upper_turn_z_m':top,'upper_then_deep_km':deep+up,'phase_speed_mps':cp,'nearest_node_distance_m':float(np.min(abs(self.z-bottom))),'status':'FINITE'}
 def ode(self,zs,angle,step):
  def fun(r,y):
   z=float(np.clip(y[0],0,self.H));return [math.tan(y[1]),-float(self.df(z))/float(self.f(z))]
  def back(r,y):return y[0]-zs
  back.direction=-1;back.terminal=True
  def bed(r,y):return self.H-y[0]
  bed.direction=-1;bed.terminal=True
  def top(r,y):return y[0]
  top.direction=-1;top.terminal=True
  def apex(r,y):return y[1]
  apex.direction=-1;apex.terminal=False
  sol=solve_ivp(fun,(0.,120000.),[zs,math.radians(angle)],method='DOP853',rtol=1e-10,atol=[1e-5,1e-11],max_step=step,events=(back,bed,top,apex))
  invariant=np.cos(sol.y[1])/self.f(np.clip(sol.y[0],0,self.H));p0=math.cos(math.radians(angle))/float(self.f(zs));err=float(np.max(abs(invariant/p0-1)))
  out={'success':bool(sol.success),'returned':bool(len(sol.t_events[0])),'range_km':float(sol.t_events[0][0]/1000) if len(sol.t_events[0]) else None,'turn_z_m':float(sol.y_events[3][0,0]) if len(sol.y_events[3]) else None,'max_relative_Snell_invariant_error':err,'rhs_evaluations':sol.nfev,'step_max_m':step,'message':sol.message,'n_points':len(sol.t)}
  return out,sol.t,sol.y

def extrema(rows,stride=1):
 sel=rows[::stride];out=[]
 for i in range(1,len(sel)-1):
  a,b,c=sel[i-1:i+2]
  if any(x['status']!='FINITE' or x['deep_bottom_bounce'] or not x['deep_refracted_clearance_PASS'] for x in (a,b,c)):continue
  if b['deep_return_km']<a['deep_return_km'] and b['deep_return_km']<c['deep_return_km'] and 10<=b['deep_return_km']<=120:out.append(b)
 return out

def path_archive(p,zs,H,lo,hi):
 """Saved downward ray fan: crossings labelled by boundary contacts in the prefix, never total counters."""
 lines=p.read_text().splitlines();assert 'rz' in lines[6].lower();i=7;rows=[]
 while i<len(lines):
  if not lines[i].strip():i+=1;continue
  ang=float(lines[i]);n,ns,nb=map(int,lines[i+1].split());i+=2;ar=np.array([[float(v) for v in l.split()] for l in lines[i:i+n]]);i+=n
  up=ar[:,1]<=1;bed=ar[:,1]>=H-1;us=np.cumsum(up & np.r_[True,~up[:-1]]);bs=np.cumsum(bed & np.r_[True,~bed[:-1]])
  for j in np.flatnonzero(((ar[:-1,1]>zs)&(ar[1:,1]<=zs))|((ar[:-1,1]<zs)&(ar[1:,1]>=zs))):
   f=(zs-ar[j,1])/(ar[j+1,1]-ar[j,1]);rr=(ar[j,0]+f*(ar[j+1,0]-ar[j,0]))/1000
   if not 10<=rr<=120:continue
   sn=int(us[j+1]);bn=int(bs[j+1]);kind='ZERO_BOUNDARY' if sn==bn==0 else ('SURFACE_ONLY' if sn>0 and bn==0 else ('BOTTOM_ONLY' if bn>0 and sn==0 else 'MIXED_BOUNDARY'))
   rows.append({'angle_deg':ang,'return_km':float(rr),'direction':'UP' if ar[j+1,1]<ar[j,1] else 'DOWN','prefix_surface_contacts':sn,'prefix_bottom_contacts':bn,'identity':kind,'in_original_focus_band':bool(lo<=rr<=hi),'whole_trace_surface_count':ns,'whole_trace_bottom_count':nb})
 return rows
