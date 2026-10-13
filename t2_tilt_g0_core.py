"""T2-TILT-G0 common IO and exact all-mode static provider; no historical code edits."""
from pathlib import Path
import json,csv,hashlib,datetime,time,os
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_T2_TILT_VS_SHORT_VERTICAL_G0'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_TILT_G0')
PKG=Path('D:/Downloads/T2_TILT_G0_TASK_PACKAGE_20261013_v1')
PARENT='251aa604af546827ab10080861c776d8666fe780'
TASKSHA='7dddc10e1d504d967c79d96e7437c73eaa1a874390c3e0ae6535a8e0b4a1fb94'
RES=['A0_0','A_P05','A_M05','A_P10','A_M10','A_P4M','A_M4M','B1_V4','B2_DUP']
FREQ=np.array([201.,235.,283.]); K=8; N=2048; FS=1024
AT=ROOT/'tools/acoustics_toolbox/atWin10/at/bin'
TRUTH=ROOT/'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv'
def utc(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''): h.update(b)
 return h.hexdigest()
def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x): Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def table(p,a,fields=None):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(a[0]));w.writeheader();w.writerows(a)
def geometries():
 a=rows(OUT/'ARRAY_GEOMETRY_TABLE.csv');return {r:np.array([[float(t[c]) for c in ['x_m','y_m','z_m']] for t in sorted([x for x in a if x['resource']==r],key=lambda x:int(x['channel']))]) for r in RES}
def catalogue(fine=False):
 rs=np.linspace(45,60,121 if fine else 61);ths=np.linspace(0,5,21 if fine else 11);zs=np.arange(180,221,10,dtype=float)
 return np.stack(np.meshgrid(rs,ths,zs,indexing='ij'),axis=-1).reshape(-1,3)
def normalize(x):
 den=np.linalg.norm(x,axis=-1,keepdims=True)
 if not np.all(np.isfinite(x)) or np.any(den==0):raise ValueError('undefined/zero spatial vector')
 return x/den

def score(h,y):
 a=normalize(h);b=normalize(y)
 return np.mean(np.abs(np.einsum('cfm,kfm->ckf',a.conj(),b,optimize=True))**2,axis=(1,2))
def covariance(y):
 b=normalize(y);return np.einsum('kfm,kfn->fmn',b,b.conj(),optimize=True)/K

def scores_many(h,ys):
 a=normalize(h);b=normalize(ys);c=np.einsum('rkfm,rkfn->rfmn',b,b.conj(),optimize=True)/K
 return np.einsum('cfm,rfmn,cfn->rc',a.conj(),c,a,optimize=True).real/3

def shape_error(pred,ref):
 den=np.sum(abs(pred)**2,axis=-1)
 if np.any(den==0) or np.any(np.linalg.norm(ref,axis=-1)==0):raise ValueError('zero field, cannot align')
 c=np.sum(pred.conj()*ref,axis=-1)/den
 return np.linalg.norm(ref-c[...,None]*pred,axis=-1)/np.linalg.norm(ref,axis=-1)

def parse_mod(path,expected=None):
 b=Path(path).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);hdr=np.frombuffer(b[84:108],'<i4');nt,nm=map(int,hdr[2:4]);z=np.frombuffer(b[4*rec:5*rec],'<f4')[:nt].astype(float);nmodes=int(np.frombuffer(b[5*rec:5*rec+4],'<i4')[0]);phi=np.stack([np.frombuffer(b[(7+j)*rec:(8+j)*rec],'<c8')[:nm] for j in range(nmodes)],axis=1).astype(complex);kk=np.frombuffer(b[(7+nmodes)*rec:(7+nmodes)*rec+8*nmodes],'<c8').astype(complex)
 if nt!=nm or not np.isfinite(phi).all() or not np.isfinite(kk).all() or np.any(kk.real<=0) or np.any(kk.imag>1e-12):raise ValueError('invalid MOD header or attenuation')
 if expected is not None and any(np.min(abs(z-d))>1e-5 for d in expected):raise ValueError('exact requested depth unavailable')
 return z,phi,kk

def z_index(z,d):
 idx=int(np.argmin(abs(z-d)))
 if abs(z[idx]-d)>1e-5:raise ValueError('unregistered source/receiver depth')
 return idx

class Provider:
 def __init__(self,mesh,folder=None,expected=None):
  self.mesh=mesh;self.mods=[parse_mod((folder or LOCAL/'modal')/f't2_f{int(f)}_n{mesh}.mod',expected) for f in FREQ]
 def field(self,states,coords):
  coords,coord_inverse=np.unique(np.asarray(coords,float),axis=0,return_inverse=True);s=np.atleast_2d(np.asarray(states,float));hor,inv=np.unique(s[:,:2],axis=0,return_inverse=True);zs,iz=np.unique(s[:,2],return_inverse=True);rr=np.sqrt((hor[:,0,None]*1000*np.cos(np.radians(hor[:,1,None]))-coords[:,0])**2+(hor[:,0,None]*1000*np.sin(np.radians(hor[:,1,None]))-coords[:,1])**2)
  out=np.empty((len(s),3,len(coords)),complex)
  for fi,(depths,phi,kk) in enumerate(self.mods):
   ps=phi[[z_index(depths,x) for x in zs]]
   for depth in np.unique(coords[:,2]):
    use=np.flatnonzero(coords[:,2]==depth);ranges=rr[:,use].ravel();w=ps*phi[z_index(depths,depth)]
    vals=np.empty((len(zs),len(ranges)),complex)
    for a in range(0,len(ranges),512):
     r=ranges[a:a+512];ker=np.sqrt(2*np.pi/(kk[:,None]*r))*np.exp(-1j*kk[:,None]*r-1j*np.pi/4);vals[:,a:a+512]=w@ker
    vals=vals.reshape(len(zs),len(hor),len(use));out[:,fi,use]=vals[iz,inv]
  if not np.isfinite(out).all():raise ValueError('nonfinite full-modal field')
  return out[:,:,coord_inverse]

TONE=np.exp(2j*np.pi*FREQ[:,None]*np.arange(N)[None,:]/FS)
def waveform(h,source,gain,noise,nu):
 return 2*np.real(np.einsum('kfm,fn->kmn',source[:,:,None]*gain[None,None,:]*h[None],TONE,optimize=True))+noise*np.sqrt(N*nu)
def dft(x): return np.einsum('kmn,fn->kfm',x,TONE.conj(),optimize=True)/N

def guard(start=None):
 import psutil,shutil
 rss=psutil.Process().memory_info().rss
 elapsed=time.monotonic()-start if start is not None else 0
 if rss>8*1024**3 or elapsed>21600:raise RuntimeError('hard science budget exceeded')
 if shutil.disk_usage('D:/').free<1024**3:raise RuntimeError('D disk reserve reached')
 return {'elapsed_s':elapsed,'RSS_bytes':rss}