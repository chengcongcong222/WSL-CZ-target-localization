"""AUX-D1 storage, identity and resource accounting; no private panel access."""
import csv,datetime,hashlib,json,pathlib,platform,time,ctypes,os
ROOT=pathlib.Path(__file__).resolve().parent
OUT=ROOT/'results/R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT';LOCAL=pathlib.Path('D:/ProjectStorage/WSL-CZ/AUX_D1')
PARENT='9060ca4ccd82f572ef099b185c3d742674752099';TASKSHA='d3079beff433bfa07d84a06769b45c878e18d61ad991f5b6cb1a4728025d42f5'
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(path):
 h=hashlib.sha256()
 with pathlib.Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def dump(name,x):
 p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def read(name):return json.loads((OUT/name).read_text(encoding='utf-8-sig'))
def csvwrite(name,data,fields=None):
 with (OUT/name).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(data[0]),lineterminator='\n');w.writeheader();w.writerows(data)
def csvread(name):return list(csv.DictReader((OUT/name).read_text(encoding='utf-8').splitlines()))
def format_float(x):return format(float(x),'.17g')
def filemeta(path):
 p=pathlib.Path(path);return {'path':str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
def verify():
 x=read('DESIGN_FREEZE.json')
 for p,h in x['code_and_input_sha256'].items():
  if sha(ROOT/p)!=h:raise RuntimeError('Frozen identity mismatch '+p)
 return x

def rss():
 class Info(ctypes.Structure):_fields_=[('cb',ctypes.c_ulong),('PageFaultCount',ctypes.c_ulong),('PeakWorkingSetSize',ctypes.c_size_t),('WorkingSetSize',ctypes.c_size_t),('QuotaPeakPagedPoolUsage',ctypes.c_size_t),('QuotaPagedPoolUsage',ctypes.c_size_t),('QuotaPeakNonPagedPoolUsage',ctypes.c_size_t),('QuotaNonPagedPoolUsage',ctypes.c_size_t),('PagefileUsage',ctypes.c_size_t),('PeakPagefileUsage',ctypes.c_size_t)]
 i=Info();i.cb=ctypes.sizeof(i);kernel=ctypes.windll.kernel32;kernel.GetCurrentProcess.restype=ctypes.c_void_p
 if not ctypes.windll.psapi.GetProcessMemoryInfo(ctypes.c_void_p(kernel.GetCurrentProcess()),ctypes.byref(i),i.cb):raise RuntimeError('RSS monitoring unavailable')
 return int(i.WorkingSetSize)
def storage_bytes():return sum(p.stat().st_size for base in [OUT,LOCAL] for p in base.rglob('*') if p.is_file())
def elapsed_completed():
 p=OUT/'SCIENCE_CLOCK.json';x=read('SCIENCE_CLOCK.json') if p.exists() else {'stages':{}}
 return sum(v['elapsed_s'] for v in x['stages'].values())
def record_clock(stage,elapsed,**kw):
 x=read('SCIENCE_CLOCK.json') if (OUT/'SCIENCE_CLOCK.json').exists() else {'definition':'Sum elapsed times of scientific stages, excluding authoring/Git; main includes full worker walltime. No parallel scientific processes.','stages':{}}
 if stage in x['stages']:raise RuntimeError('Stage cannot be counted/executed twice '+stage)
 x['stages'][stage]={'elapsed_s':elapsed,**kw};x['total_scientific_elapsed_s']=sum(v['elapsed_s'] for v in x['stages'].values());dump('SCIENCE_CLOCK.json',x)