"""Shared serialization and frozen constants only. No truth loading at import."""
from pathlib import Path
import json,csv,hashlib,importlib.util,math
import numpy as np
O=Path(__file__).resolve().parent; ROOT=O.parents[2]; D=O.parent
DATA=Path('D:/ProjectStorage/WSL-CZ/CZ_D2/RESUME_01_DATA')
REF=ROOT/'results/R4_T3_W1_WIDE_BASELINE_DEMO'
CFG=json.loads((D/'DESIGN_FREEZE.json').read_text(encoding='utf-8-sig'))
PARENT='eab37fe96383951bbbf57ba6b30726644b423bae';DESIGN='dc91982e67921eeb7a6aac4e12c482a980d5ad89'
V3=D/'HISTORY_REPAIR_02/HISTORICAL_STATE_EXCLUSIONS_V3.csv'
# Authorized V3 identity; no fallback to V1.
V3HASH='d70e441b8c2b743555eed9e8102d2bc75d628ad34aa0640747a493877434ba9e'
METRICS=['range_abs_km','range_rel','bearing_circular_deg','speed_abs_mps','speed_rel','heading_circular_deg']
DT=np.dtype([('executed','i1'),('accepted','i1'),('h','f8',(6,)),('initial','f8',(6,)),('solver_status','i4'),('nfev','i4'),('cost','f8'),('optimality','f8'),('fallback','i1'),('intersections','i4'),('rank','i1'),('ratio','f8'),('rms','f8'),('near','i1',(6,)),('active_bias_distance','f8',(2,)),('original_bias_distance','f8',(2,)),('outside_original','i1',(2,))])
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def js(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def clean(x):
 if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,np.generic):return clean(x.item())
 if isinstance(x,float) and not math.isfinite(x):return 'INF' if x>0 else '-INF' if x<0 else 'NOT_EVALUATED'
 return x
def dump(p,x,exclusive=False):
 with Path(p).open('x' if exclusive else 'w',encoding='utf-8') as f:json.dump(clean(x),f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def csvwrite(p,rr,fields=None):
 fields=fields or list(dict.fromkeys(k for x in rr for k in x))
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows([{k:clean(v) for k,v in x.items()} for x in rr])
def csvread(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def module(name):
 s=importlib.util.spec_from_file_location(name,REF/(name+'.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def arrayhash(*a):
 h=hashlib.sha256()
 for x in a:
  x=np.ascontiguousarray(x);h.update(str(x.dtype).encode());h.update(str(x.shape).encode());h.update(x.tobytes())
 return h.hexdigest()
def quality(q,p):
 a=dict(CFG['quality'][q]);a['common_bound_deg']*=CFG['bias_support_multiplier'][p];a['half_diff_bound_deg']*=CFG['bias_support_multiplier'][p];return a
def verify_binding():
 b=js(O/'EFFECTIVE_INPUT_BINDING.json');assert b['registry']['path']==str(V3) and b['registry']['sha256']==V3HASH and sha(V3)==V3HASH
 assert len(csvread(V3))==762
 for x in b['inputs']:assert sha(ROOT/x['path'])==x['sha256'],x['path']
 return b
def verify_code():
 verify_binding()
 for x in js(O/'PRODUCTION_CODE_FREEZE.json')['files']:assert sha(Path(x['path']))==x['sha256'],x['path']
def index(q,b,c,r,p):return ((((q*2+b)*12+c)*500+r)*2+p)
def recordpath(q,b,c):return DATA/'public'/f'{q}_{b}_D{c+1:02d}.npz'
def quantile(x,p):
 if len(x)==0:return 'NOT_EVALUATED'
 a=np.sort(x);return float(a[math.ceil(p*len(a))-1])
