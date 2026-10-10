"""SP1 pure static model and fixed normalized scores; no old imports."""
from pathlib import Path
import json,csv,hashlib,datetime,os,time
import numpy as np
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/HLA_SP1')
PARENT='f21f430fbe18f13b855042f91163d050c4fd9b78'
FREQ=np.array([201.,235.,283.]);ZS=np.arange(150.,251.,5.);OFF=np.arange(-7.,8.,2.)
PAIRS=[(1,0),(2,1),(2,0)]
def stamp():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def table(p,a):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(a[0]));w.writeheader();w.writerows(a)
def guard():
 import psutil
 elapsed=(datetime.datetime.now(datetime.timezone.utc)-datetime.datetime.fromisoformat(read(LOCAL/'RUN_BUDGET_START.json')['time_utc'])).total_seconds();rss=psutil.Process().memory_info().rss
 if elapsed>14400 or rss>8*1024**3:raise RuntimeError('hard compute budget reached')
 return elapsed,rss
def parse_mod(p):
 b=Path(p).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);hdr=np.frombuffer(b[84:108],'<i4');nt,nm=map(int,hdr[2:4]);z=np.frombuffer(b[4*rec:5*rec],'<f4')[:nt].astype(float);count=int(np.frombuffer(b[5*rec:5*rec+4],'<i4')[0]);phi=np.stack([np.frombuffer(b[(7+j)*rec:(8+j)*rec],'<c8')[:nm] for j in range(count)],axis=1).astype(complex);k=np.frombuffer(b[(7+count)*rec:(7+count)*rec+count*8],'<c8').astype(complex)
 if nt!=nm or not np.array_equal(z,ZS) or not np.isfinite(phi).all() or not np.isfinite(k).all() or not np.all(k.real>0) or not np.all(k.imag<=1e-12):raise ValueError('static provider invalid/exact depths absent/attenuation sign')
 return z,phi,k
def pressure(mod,rr,z):
 depths,phi,k=mod;r=np.asarray(rr);w=phi[np.flatnonzero(depths==z)[0]]*phi[np.flatnonzero(depths==200)[0]];v=np.sqrt(2*np.pi/(k[:,None]*r.ravel()[None,:]))*np.exp(-1j*k[:,None]*r.ravel()[None,:]-1j*np.pi/4)
 return (w@v).reshape(r.shape)
def catalogue(step):
 rs=np.round(np.arange(45,60+step/2,step),8);th=np.round(np.arange(0,5+step*1.25/2,step*2.5),8)
 return rs,th,np.stack(np.meshgrid(rs,th,ZS,indexing='ij'),axis=-1).reshape(-1,3)
def nearest(a,x):return int(np.argmin(abs(np.asarray(a)-x)))
class Provider:
 def __init__(self,mesh,cache=True):
  self.mesh=mesh;self.mods=[parse_mod(LOCAL/'modal'/f'sp1_f{int(f)}_n{mesh}.mod') for f in FREQ];self.spl=[];self.carriers=[]
  if cache:
   x=np.arange(44990.,60011.,1.)
   for fi,m in enumerate(self.mods):
    guard();k=m[2];carrier=float((k.real.min()+k.real.max())/2);self.carriers.append(carrier);path=LOCAL/f'spline_n{mesh}_f{int(FREQ[fi])}.npy'
    if path.exists():v=np.load(path)
    else:
     w=m[1]*m[1][np.flatnonzero(m[0]==200)[0]];v=np.empty((len(x),len(ZS)),complex)
     for a in range(0,len(x),512):
      rr=x[a:a+512];kernel=np.sqrt(2*np.pi/(k[:,None]*rr[None,:]))*np.exp(-1j*(k[:,None]-carrier)*rr[None,:]-1j*np.pi/4);v[a:a+512]=(w@kernel).T
     np.save(path,v)
    self.spl.append(CubicSpline(x,v,axis=0,extrapolate=False))
 def field(self,states,exact=False):
  states=np.atleast_2d(states);r=states[:,0]*1000;th=np.radians(states[:,1]);rr=np.sqrt((r[:,None]*np.cos(th[:,None])-OFF)**2+(r[:,None]*np.sin(th[:,None]))**2);out=np.empty((len(states),3,8),complex)
  for fi,m in enumerate(self.mods):
   if exact:
    for z in np.unique(states[:,2]):
     use=states[:,2]==z;out[use,fi]=pressure(m,rr[use],z)
   else:
    iz=np.array([np.flatnonzero(ZS==z)[0] for z in states[:,2]]);a=self.spl[fi](rr);out[:,fi]=np.take_along_axis(a,iz[:,None,None],axis=-1)[...,0]*np.exp(-1j*self.carriers[fi]*rr)
  if not np.isfinite(out).all():raise ValueError('nonfinite field')
  return out
def transform(v,method):
 a=np.asarray(v)
 if method:
  a=np.stack([a[...,p,:]*a[...,q,:].conj() for p,q in PAIRS],axis=-2)
 if method==2:
  if np.any(abs(a)==0):raise ValueError('phase undefined')
  a=a/abs(a)
 n=np.linalg.norm(a,axis=-1,keepdims=True)
 if np.any(n==0) or not np.isfinite(n).all():raise ValueError('undefined normalized vector')
 return a/n
def score(h,y,method):
 a=transform(h,method);b=transform(y,method)
 return np.mean(abs(np.einsum('...fm,kfm->...kf',a.conj(),b))**2,axis=(-1,-2))
def moment(y,method):
 a=transform(y,method)
 return np.einsum('kfm,kfn->fmn',a,a.conj())/len(a)
def realfeatures(a):
 # v^H C v: diagonals and all off-diagonal conjugate cross products, no independence claim.
 v=transform(a,0);diag=abs(v)**2;idx=np.triu_indices(8,1);cross=v[...,idx[0]].conj()*v[...,idx[1]]
 return np.concatenate([diag,2*cross.real,-2*cross.imag],axis=-1)
def realmoment(c):
 idx=np.triu_indices(8,1)
 return np.concatenate([np.diagonal(c,axis1=-2,axis2=-1).real,c[...,idx[0],idx[1]].real,c[...,idx[0],idx[1]].imag],axis=-1)
def intervals(ids,rs,step):
 if not len(ids):return []
 cuts=np.where(np.diff(ids)>1)[0]+1;chunks=np.split(ids,cuts)
 return [(float(max(45,rs[c[0]]-step/2)),float(min(60,rs[c[-1]]+step/2))) for c in chunks]
