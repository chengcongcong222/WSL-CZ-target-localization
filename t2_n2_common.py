"""T2-N2 common exact geometry/IO and all-returned-mode propagation; no historical mutations."""
from pathlib import Path
import csv,json,hashlib,datetime,time,os
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_N2')
PKG=Path('D:/Downloads/T2_N2_TASK_PACKAGE_20261013_v1')
G0=ROOT/'results/R4_T2_TILT_VS_SHORT_VERTICAL_G0';N1=ROOT/'results/R4_T2_N1_ARCHIVED_FIELD_SCORE_STABILITY'
OLDLOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_TILT_G0')
AT=ROOT/'tools/acoustics_toolbox/atWin10/at/bin'
PARENT='eeb06809b1b37852ef41d5d03b6348dcf95310a8';TASKSHA='e2cc30ce930932f5135bb33e63f916eefa25d6c6f9e4679a5a6f43dd4693df17'
RES=['A0_0','A_P4M','B1_V4','B2_DUP'];FREQ=[201,235,283];CATS=['primary','fine'];KVAR={'K0':(160001,1800),'K1':(160001,2000),'K2':(160001,2400),'K3':(240001,2400)};FVAR={'F1':(16001,2400),'F2':(32001,2400),'F3':(32001,3000)}

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def table(p,a,fields=None):
 fields=fields or list(dict.fromkeys(k for r in a for k in r))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(a)
def geometries():
 rr=rows(OUT/'ARRAY_GEOMETRY_TABLE.csv');return {res:np.array([[float(t[k]) for k in ['x_m','y_m','z_m']] for t in sorted([t for t in rr if t['resource']==res],key=lambda t:int(t['channel']))]) for res in RES}
def states(cat):return np.load(G0/f'CATALOGUE_{cat}.npy',mmap_mode='r',allow_pickle=False)
def sentinels():return rows(OUT/'T2_N2_SENTINELS_20261013_v1.csv')
def sstates():return np.array([[float(t[k]) for k in ['r_km','abs_theta_deg','source_z_m']] for t in sentinels()])
def ranges(s,coords):
 s=np.atleast_2d(s);return np.sqrt((s[:,0,None]*1000*np.cos(np.radians(s[:,1,None]))-coords[:,0])**2+(s[:,0,None]*1000*np.sin(np.radians(s[:,1,None]))-coords[:,1])**2)
def shape(a,b):
 n=np.linalg.norm(b,axis=-1);den=np.sum(abs(a)**2,axis=-1)
 if np.any(n==0) or np.any(den==0) or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('ZERO_OR_NONFINITE_FIELD_NOT_ZERO_ERROR')
 c=np.sum(a.conj()*b,axis=-1)/den;return np.linalg.norm(b-c[...,None]*a,axis=-1)/n

def parse_mod(p):
 b=Path(p).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);hdr=np.frombuffer(b[84:108],'<i4');nt,nm=map(int,hdr[2:4]);z=np.frombuffer(b[4*rec:5*rec],'<f4')[:nt].astype(float);count=int(np.frombuffer(b[5*rec:5*rec+4],'<i4')[0]);phi=np.stack([np.frombuffer(b[(7+j)*rec:(8+j)*rec],'<c8')[:nm] for j in range(count)],axis=1).astype(complex);k=np.frombuffer(b[(7+count)*rec:(7+count)*rec+count*8],'<c8').astype(complex)
 if nt!=nm or not np.isfinite(phi).all() or not np.isfinite(k).all() or np.any(k.real<=0) or np.any(k.imag>1e-12):raise ValueError('invalid MOD')
 return z,phi,k

def zi(z,d):
 i=int(np.argmin(abs(z-d)))
 if abs(z[i]-d)>1e-5:raise ValueError('SOURCE_OR_RECEIVER_DEPTH_NOT_READY')
 return i

def mod_path(variant,f):
 return OLDLOCAL/'modal'/f't2_f{f}_n160001.mod' if variant=='K0' else LOCAL/'K'/variant/f'n2_{variant}_f{f}.mod'

class Provider:
 def __init__(self,variant):self.variant=variant;self.mods=[parse_mod(mod_path(variant,f)) for f in FREQ]
 def field(self,s,coords,band_lower=None):
  coords,rev=np.unique(np.asarray(coords,float),axis=0,return_inverse=True);s=np.atleast_2d(s);hor,inv=np.unique(s[:,:2],axis=0,return_inverse=True);sz,iz=np.unique(s[:,2],return_inverse=True);rr=ranges(hor,coords);out=np.empty((len(s),3,len(coords)),complex);band=np.zeros_like(out)
  for fi,(depths,phi,k) in enumerate(self.mods):
   ps=phi[[zi(depths,x) for x in sz]];mask=2*np.pi*FREQ[fi]/k.real>band_lower if band_lower is not None else np.zeros(len(k),bool)
   for d in np.unique(coords[:,2]):
    cols=np.flatnonzero(coords[:,2]==d);r=rr[:,cols].ravel();w=ps*phi[zi(depths,d)];v=np.empty((len(sz),len(r)),complex);bv=np.zeros_like(v)
    for a in range(0,len(r),512):
     rpart=r[a:a+512];kernel=np.sqrt(2*np.pi/(k[:,None]*rpart))*np.exp(-1j*k[:,None]*rpart-1j*np.pi/4);v[:,a:a+512]=w@kernel
     if np.any(mask):bv[:,a:a+512]=w[:,mask]@kernel[mask]
    v=v.reshape(len(sz),len(hor),len(cols));bv=bv.reshape(len(sz),len(hor),len(cols));out[:,fi,cols]=v[iz,inv];band[:,fi,cols]=bv[iz,inv]
  if not np.isfinite(out).all():raise ValueError('invalid modal field')
  return out[:,:,rev],band[:,:,rev]

def parse_shd(p,green=False):
 b=Path(p).read_bytes();rec=4*int(np.frombuffer(b[:4],'<i4')[0]);head=np.frombuffer(b[2*rec:2*rec+32],'<i4');nf,nt,nx,ny,ns,nz,nr=map(int,head[:7]);atten=float(np.frombuffer(b[2*rec+28:2*rec+32],'<f4')[0]);freq=np.frombuffer(b[3*rec:4*rec],'<f8')[:nf].astype(float);sd=np.frombuffer(b[7*rec:8*rec],'<f4')[:ns].astype(float);zd=np.frombuffer(b[8*rec:9*rec],'<f4')[:nz].astype(float);r=np.frombuffer(b[9*rec:10*rec],'<f4')[:nr].astype(float)
 if nf!=1 or (nt,nx,ny)!=(1,1,1) or nr<2 or ns<1 or nz<1 or len(b)<(10+ns*nz)*rec:raise ValueError('SHD_GRN_FORMAT_NOT_READY')
 value=np.empty((ns,nz,nr),complex)
 for i in range(ns):
  for j in range(nz):
   offset=(10+i*nz+j)*rec;value[i,j]=np.frombuffer(b[offset:offset+8*nr],'<c8').astype(complex)
 if not np.isfinite(value).all() or not np.isfinite(r).all() or not np.isfinite(freq).all():raise ValueError('invalid SHD_GRN values')
 if not green and np.any(np.diff(r)<=0):raise ValueError('SHD_RANGE_NONMONOTONE')
 return {'freq':float(freq[0]),'atten':atten,'source_depths':sd,'receiver_depths':zd,'axis':r,'value':value,'record_bytes':rec,'plot_type':b[rec:rec+10].decode('ascii','replace')}

def interpolate_shd(d,s,coords):
 if abs(d['freq']-float(s['frequency_hz']))>1e-6:raise ValueError('FFP_FREQUENCY_NOT_READY')
 iz=zi(d['source_depths'],float(s['source_z_m']));rr=ranges([[float(s['r_km']),float(s['abs_theta_deg'])]],coords)[0];axis=d['axis'];out=[]
 if np.any(rr<axis[0]) or np.any(rr>axis[-1]):raise ValueError('FFP_NO_EXTRAPOLATION_NOT_READY')
 for i,z in enumerate(coords[:,2]):
  y=d['value'][iz,zi(d['receiver_depths'],z)];out.append(np.interp(rr[i],axis,y.real)+1j*np.interp(rr[i],axis,y.imag))
 return np.array(out)

def quant(v):
 a=np.sort(np.asarray(v));n=len(a);return {'N':n,'min':float(a[0]),'median':float(np.median(a)),'P95':float(a[max(0,int(np.ceil(.95*n))-1)]),'max':float(a[-1])}

def verify_sources():
 for r in read(OUT/'SOURCE_IDENTITY.json')['records']:
  if sha(r['path'])!=r['sha256']:raise RuntimeError('SOURCE_MUTATED '+r['path'])

def storage_guard(start=None):
 import psutil,shutil
 elapsed=0 if start is None else time.monotonic()-start;rss=psutil.Process().memory_info().rss;localbytes=sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file());repo=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
 if elapsed>21600 or rss>8*1024**3 or localbytes+repo>16*1024**3 or repo>512*1024**2 or shutil.disk_usage('D:/').free<1024**3:raise RuntimeError('HARD_BUDGET_STOP')
 return {'elapsed_s':elapsed,'RSS_bytes':rss,'new_D_bytes':localbytes+repo,'new_git_bytes':repo}