"""One-shot guarded launcher. No scientific subprocess retry is permitted."""
from pathlib import Path
import sys,os,json,hashlib,time,subprocess
import psutil
O=Path(__file__).resolve().parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZROB02')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8',newline='\n')
def usage():return sum(p.stat().st_size for root in [W,O] for p in root.rglob('*') if p.is_file())
def main():
 phase=sys.argv[1];assert phase in ('analysis','cold');marker=W/(phase.upper()+'_LAUNCH.json')
 if marker.exists():raise RuntimeError('Refusing second launch of '+phase)
 cfg=json.loads((O/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
 for a in cfg['execution_inputs']:assert sha(a['path'])==a['sha256'],a['path']
 if phase=='cold':
  assert json.loads((O/'ANALYSIS_EXIT.json').read_text(encoding='utf-8'))['completed']
  for a in json.loads((O/'PRIMARY_OUTPUT_FREEZE.json').read_text(encoding='utf-8'))['files']:assert sha(a['path'])==a['sha256'],a['path']
 start=time.perf_counter();before=usage();old=[]
 for p in W.glob('*_PROCESS_USAGE.json'):old.append(json.loads(p.read_text(encoding='utf-8')))
 elapsed_prior=sum(a['elapsed_s'] for a in old);remaining=cfg['science_wall_max_s']-elapsed_prior;assert remaining>0
 write(marker,{'phase':phase,'time_unix':time.time(),'inputs_sha_verified':len(cfg['execution_inputs']),'freeze_sha':sha(O/'DESIGN_FREEZE.json'),'remaining_wall_s':remaining})
 env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env['PYTHONIOENCODING']='utf-8';env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
 script='analysis.py' if phase=='analysis' else 'cold.py';peak=0;reason=None;samples=0
 with (W/(phase+'.log')).open('wb') as f:
  child=subprocess.Popen([sys.executable,str(O/'scripts'/script)],cwd=O.parents[1],env=env,stdout=f,stderr=subprocess.STDOUT)
  while child.poll() is None:
   try:
    root=psutil.Process(child.pid);tree=[root]+root.children(recursive=True);rss=sum(x.memory_info().rss for x in tree if x.is_running())+psutil.Process().memory_info().rss
   except psutil.Error:rss=peak
   peak=max(peak,rss);samples+=1
   if peak>cfg['process_tree_RSS_max_bytes']:reason='MEMORY_LIMIT'
   elif time.perf_counter()-start>remaining:reason='CUMULATIVE_WALL_LIMIT'
   elif usage()>cfg['new_D_max_bytes']:reason='NEW_DATA_LIMIT'
   if reason:
    for p in reversed(tree):
     try:p.kill()
     except psutil.Error:pass
    child.wait();break
   time.sleep(.2)
  code=child.wait()
 result={'phase':phase,'elapsed_s':time.perf_counter()-start,'tree_peak_RSS_bytes':peak,'samples':samples,'new_D_total_bytes':usage(),'returncode':code,'hard_limit_failure':reason,'new_solver_calls':0,'new_MC':0}
 if phase=='analysis' and code==0 and reason is None:
  files=[p for root in [O,W] for p in root.rglob('*') if p.is_file() and p.suffix.lower() in ('.csv','.json','.npz') and p.name!='PRIMARY_OUTPUT_FREEZE.json'];write(O/'PRIMARY_OUTPUT_FREEZE.json',{'files':[{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size} for p in files]})
 write(W/(phase.upper()+'_PROCESS_USAGE.json'),result);write(O/(phase.upper()+'_PROCESS_USAGE.json'),result);print(json.dumps(result));return 0 if code==0 and reason is None else 2
if __name__=='__main__':raise SystemExit(main())

