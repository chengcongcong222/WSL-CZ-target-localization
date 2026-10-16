"""CZO04 single-frequency, single-depth archived energy model. No solver imports."""
from pathlib import Path
import json,csv,hashlib,struct
import numpy as np
from scipy.linalg import solve_triangular
from scipy.stats import chi2
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZO04')
T=np.arange(0.,1201.,10.);V=np.column_stack([np.ones(121),(T-600)/600])
SCALE=np.array([1.,1.,.2,1.]);STEPS=np.array([.001,.0001,.0001,.0001])

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def table(p,a):
 if not a:return
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(a[0]));w.writeheader();w.writerows(a)
def shd(p):
 b=Path(p).read_bytes();rec=4*struct.unpack_from('<i',b)[0];dims=struct.unpack_from('<7i',b,2*rec);nz,nr=dims[-2:]
 assert dims[:5]==(1,1,1,1,1) and (nz,nr)==(81,881)
 depths=np.frombuffer(b,'<f4',nz,8*rec).astype(float);r=np.frombuffer(b,'<f4',nr,9*rec).astype(float)
 j=int(np.flatnonzero(depths==200)[0]);q=np.frombuffer(b,'<c8',nr,(10+j)*rec).astype(complex);en=q.real*q.real+q.imag*q.imag
 assert np.isfinite(en).all() and np.all(en[(r>=40000)&(r<=66000)]>0)
 return r,en

def curves(level):
 arr=[]
 for e in ['JW','JA','ANNUAL']:
  for z in [150,200,250]:
   r,p=shd(Path(f'D:/ProjectStorage/WSL-CZ/ENV_CAP/E1/SCS17_{e}_z{z}/I{level}.shd'))
   # Zero/out-of-domain values are not floored into measurable observations.
   y=np.full(p.shape,np.nan);good=p>0;y[good]=10*np.log10(p[good]);arr.append(y)
 return r,np.array(arr)
def grid():return np.round(np.array([45.,-5.,1.,-15.])+np.indices((16,21,11,31)).reshape(4,-1).T*np.array([1.,.5,.2,1.]),10)
def states():
 a=[]
 h=load(R/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/DESIGN_FREEZE.json')['geometry']
 with (R/'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv').open(encoding='utf-8-sig') as f:
  k=list(csv.DictReader(f))
 for name,s in [(f'G{i+1:02}',x) for i,x in enumerate(h)]+[(q['case_id'],[float(q[c]) for c in ['r0_km','theta0_deg','v_mps','psi_deg']]) for q in k]:
  for sign in [1,-1]:a.append(dict(id=name+('P' if sign==1 else 'M'),h=[s[0],sign*s[1],s[2],sign*s[3]],historical_grid= name.startswith('G')))
 return a

def geom(h):
 s=np.atleast_2d(h);r,th,v,ps=s.T;th=np.deg2rad(th);ps=np.deg2rad(ps)
 dt=np.maximum(T-600,0);px=2*np.minimum(T,600)+2*dt*np.cos(np.deg2rad(15));py=2*dt*np.sin(np.deg2rad(15))
 x=1000*r[:,None]*np.cos(th[:,None])+v[:,None]*np.cos(ps[:,None])*T-px
 y=1000*r[:,None]*np.sin(th[:,None])+v[:,None]*np.sin(ps[:,None])*T-py
 angle=np.arctan2(y,x)-np.where(T>600,np.deg2rad(15),0)
 # Unsigned forward-sector line-array angle. Never supply the true side.
 return np.abs(angle),np.hypot(x,y)
def feature(h,level,r,ys):return np.interp(geom(h)[1],r,ys[level])
def weight(rho):
 c=rho**np.abs(np.arange(121)[:,None]-np.arange(121)[None,:]);l=np.linalg.cholesky(c)
 return solve_triangular(l,np.eye(121),lower=True)
def project(a,wh,condition):
 b=np.asarray(a)@wh.T
 if condition=='C2':return np.zeros_like(b)
 cols=1 if condition=='C0' else 2;q=np.linalg.qr(wh@V[:,:cols],mode='reduced')[0]
 return b-(b@q)@q.T

def jac(h,r,y,fac=1):
 jb=[];jl=[]
 for j in range(4):
  d=np.zeros(4);d[j]=STEPS[j]*fac
  ba,ra=geom(np.asarray(h)+d);bb,rb=geom(np.asarray(h)-d)
  jb.append((ba[0]-bb[0])/(2*d[j]));jl.append((np.interp(ra[0],r,y)-np.interp(rb[0],r,y))/(2*d[j]))
 return np.array(jb).T*SCALE,np.array(jl).T*SCALE

def info(j):
 u,s,v=np.linalg.svd(j,full_matrices=False);tol=s[0]*1e-10 if s.size and s[0]>0 else 1e-30;rank=int((s>tol).sum())
 other=j[:,1:];uo,so,vo=np.linalg.svd(other,full_matrices=False);qo=uo[:,so>max(so[0]*1e-10,1e-30)]
 vec=j[:,0]-qo@(qo.T@j[:,0]);eff=float(vec@vec)
 return dict(rank=rank,condition=float(s[0]/s[-1]) if rank==4 else None,range_effective_per_km2=eff,weak_direction=v[-1].tolist(),eigenvalues=(s*s).tolist())
