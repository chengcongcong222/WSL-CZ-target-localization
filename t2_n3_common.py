"""Frozen T2-N3 IO, exact stored-depth fields, closed-form ideal references."""
from pathlib import Path
import csv,json,hashlib,datetime,time,os,struct
import numpy as np
from scipy.special import hankel2
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_T2_N3_BOUNDARY_PHYSICS_ANALYTIC_REFERENCE'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_N3')
PKG=Path('D:/Downloads/T2_N3_TASK_PACKAGE_20261013_v1')
N2=ROOT/'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET'
G0=ROOT/'results/R4_T2_TILT_VS_SHORT_VERTICAL_G0'
AT=ROOT/'tools/acoustics_toolbox/atWin10/at/bin'
PARENT='1ddc329edec9cbf07ffd6f2f5bbd63805104964e'
TASKSHA='5082da133db2ab54de93c71cd2179f8b2df531cc5dbefbcc86335994c714e127'
RES=('A0_0','A_P4M','B1_V4','B2_DUP');FREQ=(201,235,283)
PROBES=(('P01',45.5,0.,180.),('P02',50.,0.,200.),('P03',55.,2.5,200.),('P04',58.,4.5,210.),('P05',59.5,5.,220.))
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
 fields=fields or list(dict.fromkeys(k for r in a for k in r)) or ['status']
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(a)
def arrays():
 a=rows(OUT/'package/T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv')
 return {n:np.array([[float(r[k]) for k in ('x_m','y_m','z_m')] for r in sorted([r for r in a if r['resource']==n],key=lambda r:int(r['channel']))]) for n in RES}
def sentinels():return rows(OUT/'package/T2_N2_SENTINELS_20261013_v1.csv')
def mod_path(v,f):
 return Path('D:/ProjectStorage/WSL-CZ/T2_TILT_G0/modal')/f't2_f{f}_n160001.mod' if v=='K0' else Path('D:/ProjectStorage/WSL-CZ/T2_N2/K')/v/f'n2_{v}_f{f}.mod'
def parse_mod(p):
 b=Path(p).read_bytes();rec=4*struct.unpack('<i',b[:4])[0];nf,nmedia,nt,nm=struct.unpack('<4i',b[84:100]);f=struct.unpack('<d',b[3*rec:3*rec+8])[0]
 z=np.frombuffer(b[4*rec:5*rec],'<f4')[:nt].astype(float);M=struct.unpack('<i',b[5*rec:5*rec+4])[0]
 if M<1 or nf!=1 or nt!=nm:raise ValueError('INVALID_MODAL_HEADER')
 phi=np.stack([np.frombuffer(b[(7+j)*rec:(8+j)*rec],'<c8')[:nm] for j in range(M)],axis=1).astype(complex)
 k=np.frombuffer(b[(7+M)*rec:(7+M)*rec+M*8],'<c8').astype(complex)
 if len(k)!=M or not np.isfinite(phi).all() or not np.isfinite(k).all() or np.any(k.real<=0) or np.any(k.imag>1e-12):raise ValueError('INVALID_MODAL_DATA')
 return z,phi,k,f

def zi(z,d):
 i=int(np.argmin(abs(z-d)))
 if abs(z[i]-d)>1e-5:raise ValueError('DEPTH_NOT_READY')
 return i

def shape(a,b):
 na=np.linalg.norm(a);nb=np.linalg.norm(b)
 if na==0 or nb==0 or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('UNDEFINED_NORMALIZATION')
 u=a/na;v=b/nb;return float(np.linalg.norm(u-v*np.vdot(v,u)))
def distance(pos,coords):
 r,t,z=pos;xy=1000*r*np.array([np.cos(np.radians(t)),np.sin(np.radians(t))]);return np.linalg.norm(xy-coords[:,:2],axis=1)
def modal_field(mod,pos,coords,mask=None):
 z,phi,k,f=mod
 if mask is None:mask=np.ones(len(k),bool)
 kk=k[mask];weights=phi[zi(z,pos[2]),mask,None]*phi[[zi(z,x) for x in coords[:,2]]][:,mask].T
 rr=distance(pos,coords);kernel=np.sqrt(2*np.pi/(kk[:,None]*rr))*np.exp(-1j*kk[:,None]*rr-1j*np.pi/4)
 return np.sum(weights*kernel,axis=0,dtype=np.complex128)
def spectrum(f,bc,high):
 order=np.arange(int(np.ceil(2*f*5000/1500))+2,dtype=float)+(.5 if bc=='PR_R' else 1.)
 order=order[order<2*f*5000/1500];mu=order*np.pi/5000;kr=np.sqrt((2*np.pi*f/1500)**2-mu*mu);cp=2*np.pi*f/kr;mask=cp<high
 return mu[mask],kr[mask],cp[mask]
def analytic_field(f,bc,high,pos,coords,exact=True):
 mu,kr,cp=spectrum(f,bc,high);rr=distance(pos,coords);x=kr[:,None]*rr
 p=.0004*np.sin(mu*pos[2])[:,None]*np.sin(mu[:,None]*coords[None,:,2])
 radial=-1j*np.pi*hankel2(0,x) if exact else np.sqrt(2*np.pi/x)*np.exp(-1j*x-1j*np.pi/4)
 return np.sum(p*radial,axis=0,dtype=np.complex128)
def env_data(p):
 lines=Path(p).read_text(encoding='utf-8-sig').splitlines();b=next(i for i in range(5,len(lines)) if lines[i].strip().startswith(("'R'","'V'","'A'")))
 ssp=np.array([[float(v) for v in l.split()[:6]] for l in lines[5:b]])
 ns=int(lines[b+3]);sd=np.array([float(x) for x in lines[b+4:b+4+ns]]);j=b+4+ns;nr=int(lines[j]);rd=np.array([float(x) for x in lines[j+1:j+1+nr]])
 return {'lines':lines,'b':b,'ssp':ssp,'option':lines[3].strip().strip("'"),'bottom':lines[b].split()[0].strip("'"),'f':float(lines[1]),'mesh':int(lines[4].split()[0]),'clow':float(lines[b+1].split()[0]),'chigh':float(lines[b+1].split()[1]),'rmax':float(lines[b+2]),'source':sd,'receiver':rd}
def verify_sources():
 for r in read(OUT/'OLD_N2_SOURCE_IDENTITIES.json')['records']:
  if sha(r['path'])!=r['sha256']:raise RuntimeError('SOURCE_CHANGED '+r['path'])
def verify_freeze():
 for r in read(OUT/'DESIGN_FREEZE.json')['frozen_files']:
  if sha(r['path'])!=r['sha256']:raise RuntimeError('FROZEN_CHANGED '+r['path'])
 verify_sources()
def guard(start=None):
 import psutil,shutil
 rss=psutil.Process().memory_info().rss;wall=0 if start is None else time.monotonic()-start
 a=sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file());b=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
 if wall>7200 or rss>8*1024**3 or a+b>4*1024**3 or b>256*1024**2 or shutil.disk_usage('D:/').free<1024**3:raise RuntimeError('BUDGET_LIMIT_STOP')
 return {'wall_s':wall,'RSS_bytes':rss,'D_bytes':a+b,'git_output_bytes':b}