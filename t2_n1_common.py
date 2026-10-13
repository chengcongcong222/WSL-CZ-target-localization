"""T2-N1 archive-only diagnostics: IO, frozen identities, stable projectors, no propagation imports."""
from pathlib import Path
import hashlib,json,csv,datetime,time,os
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_T2_N1_ARCHIVED_FIELD_SCORE_STABILITY'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_N1')
OLD=ROOT/'results/R4_T2_TILT_VS_SHORT_VERTICAL_G0'
ARCHIVE=Path('D:/ProjectStorage/WSL-CZ/T2_TILT_G0')
PKG=Path('D:/Downloads/T2_N1_TASK_PACKAGE_20261013_v1')
PARENT='6b60c39f52b68a361b2c005e74126b9f68401b77'
TASKSHA='36a0d6341046263490a69ae12dff5f64e25d0e4c4d7068dd6327344b3d0f0fc9'
RES=['A0_0','A_P05','A_M05','A_P10','A_M10','A_P4M','A_M4M','B1_V4','B2_DUP']
FREQ=[201,235,283]; CATS=['primary','fine']; PROBES=[('P01',45.5,0.,180.),('P02',50.,0.,200.),('P03',55.,2.5,200.),('P04',58.,4.5,210.),('P05',59.5,5.,220.)]
PAIRS=[('A','A_P4M','A0_0'),('B','B1_V4','B2_DUP')]; TOL=1e-10

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def table(p,a,fields=None):
 columns=fields or list(dict.fromkeys(k for r in a for k in r))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=columns,lineterminator='\n');w.writeheader();w.writerows(a)
def catalog(cat):return np.load(OLD/f'CATALOGUE_{cat}.npy',mmap_mode='r',allow_pickle=False)
def archived(res,cat,mesh):
 name=f'{res}_{cat}_n{mesh}.npy';item=next(r for r in read(OUT/'INPUT_IDENTITY_AND_SCOPE.json')['arrays'] if r['filename']==name)
 return np.load(item['path'],mmap_mode='r',allow_pickle=False)
def normalize(h):
 norm=np.linalg.norm(h,axis=-1)
 bad=(norm==0)|(~np.isfinite(norm))|(~np.isfinite(h).all(axis=-1))
 if np.any(bad):raise ValueError('ARCHIVED_INPUT_INVALID_VECTOR')
 return h/norm[...,None],norm

def residual_error(hi,lo):
 uh,nh=normalize(hi);ul,nl=normalize(lo);corr=np.sum(ul.conj()*uh,axis=-1);e=np.linalg.norm(uh-corr[...,None]*ul,axis=-1)
 return uh,ul,nh,nl,e

def probe_scores(u,w):
 return np.mean(abs(np.einsum('cfm,fm->cf',u.conj(),w,optimize=True))**2,axis=1)
def ranking(v):return np.lexsort((np.arange(len(v)),-v))
def quantiles(v):
 a=np.sort(np.asarray(v,float));n=len(a)
 if not n:return {'N':0,'min':'NOT_EVALUATED','median':'NOT_EVALUATED','P90':'NOT_EVALUATED','P95':'NOT_EVALUATED','P99':'NOT_EVALUATED','max':'NOT_EVALUATED'}
 return {'N':n,'min':float(a[0]),'median':float(np.median(a)),'P90':float(a[math_rank(n,.9)]),'P95':float(a[math_rank(n,.95)]),'P99':float(a[math_rank(n,.99)]),'max':float(a[-1])}
def math_rank(n,p):return max(0,int(np.ceil(n*p))-1)
def shape_stats(v):
 d=quantiles(v);d.update({'above_historical_0p002':int(np.sum(np.asarray(v)>.002)),'above_fraction':float(np.mean(np.asarray(v)>.002)) if len(v) else 'NOT_EVALUATED'});return d

def verify_identities():
 d=read(OUT/'INPUT_IDENTITY_AND_SCOPE.json')
 for item in d['sources']+d['arrays']:
  if not Path(item['path']).is_file() or sha(item['path'])!=item['sha256']:raise RuntimeError('ARCHIVE_INPUT_INCOMPLETE_OR_MUTATED '+item['path'])
 return len(d['sources'])+len(d['arrays'])

def budget(start=None):
 import psutil
 rss=psutil.Process().memory_info().rss;sec=time.monotonic()-start if start is not None else 0
 localbytes=sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file()) if LOCAL.exists() else 0;gitbytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()) if OUT.exists() else 0
 if sec>1800 or rss>4*1024**3 or localbytes+gitbytes>2*1024**3 or gitbytes>64*1024**2:raise RuntimeError('PARTIAL_BUDGET')
 return {'wall_s':sec,'RSS_bytes':rss,'new_D_bytes':localbytes+gitbytes,'new_git_result_bytes':gitbytes}