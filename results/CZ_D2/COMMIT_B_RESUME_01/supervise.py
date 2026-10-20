"""Independent watchdog and exclusive locks. One main and one cold call maximum."""
from common import *
import os,sys,time,subprocess,psutil

def main():
 assert all(os.environ.get(k)=='1' for k in ['OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'])
 verify_code();assert js(O/'PRODUCTION_PREFLIGHT.json')['PASS'];assert not DATA.exists()
 dump(O/'FORMAL_EXECUTION_STARTED.json',dict(parent=PARENT,design=DESIGN,code_sha256=sha(O/'PRODUCTION_CODE_FREEZE.json'),binding_sha256=sha(O/'EFFECTIVE_INPUT_BINDING.json'),main_calls=1,cold_max_calls=1,retries=0,pid=os.getpid()),True)
 start=time.monotonic();peak=0;stages=[];failure=None;logs=O/'runtime_logs';logs.mkdir()
 for name,cap in [('main_run.py',5400),('cold_audit.py',1800)]:
  if name=='cold_audit.py' and not all((O/n).exists() for n in ['MAIN_COMPLETION.json','EVALUATION_ARRAYS.npz','TERMINALS_ARCHIVE.npz']):break
  verify_code();lock='MAIN_STARTED.json' if name=='main_run.py' else 'COLD_STARTED.json';dump(O/lock,dict(program=name,parent=PARENT,code_sha256=sha(O/'PRODUCTION_CODE_FREEZE.json'),seconds_since_start=time.monotonic()-start),True)
  begin=time.monotonic();last=begin;path=logs/(name+'.log');reason=None
  with path.open('xb') as f:
   child=subprocess.Popen([sys.executable,'-B',str(O/name)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
   while child.poll() is None:
    rss=psutil.Process().memory_info().rss
    for proc in psutil.Process().children(recursive=True):
     try:rss+=proc.memory_info().rss
     except psutil.NoSuchProcess:pass
    peak=max(peak,rss);now=time.monotonic();data=sum(p.stat().st_size for p in DATA.rglob('*') if p.is_file());gitbytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file())
    if now-begin>cap or now-start>7200:reason='WALL_BUDGET'
    if rss>4*1024**3:reason='RSS_BUDGET'
    if data>2*1024**3 or gitbytes>64*1024**2:reason='STORAGE_BUDGET'
    if reason:
     for proc in psutil.Process(child.pid).children(recursive=True):
      try:proc.kill()
      except psutil.NoSuchProcess:pass
     child.kill();child.wait();break
    if now-last>=20:
     record=dict(stage=name,wall_s=now-start,stage_wall_s=now-begin,tree_RSS_bytes=rss,D_bytes=data,Git_bytes=gitbytes);dump(O/'EXECUTION_PROGRESS.json',record);print('HEARTBEAT',json.dumps(record),flush=True);last=now
    time.sleep(.2)
  stages.append(dict(program=name,wall_s=time.monotonic()-begin,cap_s=cap,returncode=child.returncode,termination=reason,log=str(path),log_sha256=sha(path)));print('STAGE_COMPLETE',stages[-1],flush=True)
  if child.returncode!=0:failure=failure or reason or 'SCIENTIFIC_PROCESS_FAILED_NO_RETRY'
  # No retry; one independent audit of completed data is still permitted after main failure.
  if time.monotonic()-start>=7200:break
 dump(O/'EXECUTION_BUDGET.json',dict(main_calls=sum(x['program']=='main_run.py' for x in stages),cold_calls=sum(x['program']=='cold_audit.py' for x in stages),retries=0,total_wall_s=time.monotonic()-start,peak_tree_RSS_bytes_sampled=peak,monitor_sleep_s=.2,peak_is_sampled_not_continuous_guarantee=True,D_bytes=sum(p.stat().st_size for p in DATA.rglob('*') if p.is_file()),Git_bytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),stages=stages,failure=failure),True)
 dump(O/'SCIENCE_STOP.json',dict(STOP=True,formal_calls_ended=True,subsequent_science_authorized=False,CZ_D3_authorized=False,depth='NOT_EVALUATED',R4_percent=0,failure=failure),True)
 print('STOP',failure,flush=True)
if __name__=='__main__':main()
