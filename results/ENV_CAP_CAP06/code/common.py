from pathlib import Path
import json,csv,hashlib, numpy as np
O=Path(__file__).resolve().parents[1]
W=Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_CAP06')
RAW=W/'raw'
def read(n):return json.loads((O/n).read_text(encoding='utf-8-sig'))
def write(n,obj):
 p=O/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def tab(n,rows):
 with (O/n).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['status']);w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def platform(t):
 t=np.asarray(t);u=np.minimum(t,600);v=np.maximum(t-600,0)
 return np.column_stack((2*u+2*v*np.cos(np.deg2rad(15)),2*v*np.sin(np.deg2rad(15))))
def axis(t):return np.where(np.asarray(t)>600,np.deg2rad(15),0.)
def cn(rng,shape):return (rng.normal(size=shape)+1j*rng.normal(size=shape))/np.sqrt(2)
def wilson(k,n):
 if not n:return [None,None]
 z=1.959963984540054;p=k/n;den=1+z*z/n
 return [(p+z*z/(2*n)+sgn*z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/den for sgn in [-1,1]]
def qtiles(x):
 x=np.sort(np.asarray(x,float));n=len(x)
 # inverted empirical CDF, also well-defined with INF failures
 return [float(x[min(n-1,int(np.ceil(p*n))-1)]) if np.isfinite(x[min(n-1,int(np.ceil(p*n))-1)]) else 'INF' for p in [.5,.9,.95,.99]] if n else [None]*4
