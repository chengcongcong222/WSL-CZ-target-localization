"""RX1 shared immutable utilities. No import-time writes or historical imports."""
from pathlib import Path
import json,csv,hashlib,time
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_SINGLE_HLA_RADIAL_RECEIVER_EXTRACTION_PILOT'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/HLA_RX1')
PARENT='132cdd22eba8ce7f5965b342a588c7697c6c772e'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(p,rows):
 if not rows:return
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def platform(t):
 t=np.asarray(t);post=np.maximum(t-600.,0.)
 return np.stack([2*np.minimum(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)],axis=-1)
def state(row):
 r,th,v,ps=[float(row[x]) for x in ['r_km','theta_deg','v_mps','psi_deg']]
 return 1000*r*np.array([np.cos(np.deg2rad(th)),np.sin(np.deg2rad(th))]),v*np.array([np.cos(np.deg2rad(ps)),np.sin(np.deg2rad(ps))])
def guard(start,limit=12600):
 if time.monotonic()-start>limit:raise TimeoutError('frozen compute budget')
 import psutil
 if psutil.Process().memory_info().rss>8*1024**3:raise MemoryError('frozen 8GiB RSS budget')
 if sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file())>32*1024**3:raise OSError('frozen 32GiB local budget')
