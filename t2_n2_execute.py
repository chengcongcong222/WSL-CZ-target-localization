"""Once-only N2 pipeline, one numeric worker, whole-tree RSS and6h hard budget."""
from t2_n2_common import *
import subprocess,sys,psutil

def main():
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('NO_SECOND_EXECUTION')
 for p,h in read(OUT/'DESIGN_FREEZE.json')['code_and_input_sha256'].items():assert sha(ROOT/p)==h,p
 verify_sources();head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==head and head!=PARENT
 dump(OUT/'EXECUTION_STARTED.json',{'design_SHA':head,'utc':utc(),'once_only':True});total=0.;records=[]
 for stage in ['t2_n2_run.py','t2_n2_evaluate.py','t2_n2_audit.py','t2_n2_report.py']:
  log=LOCAL/(stage+'.log');beg=time.monotonic();peak=0
  with log.open('wb',buffering=0) as stream:
   p=subprocess.Popen([sys.executable,'-B',stage],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
   while p.poll() is None:
    try:proc=psutil.Process(p.pid);rss=sum(r.memory_info().rss for r in [proc]+proc.children(recursive=True) if r.is_running());peak=max(peak,rss)
    except (psutil.NoSuchProcess,psutil.AccessDenied):rss=0
    if total+time.monotonic()-beg>21600 or rss>8*1024**3:
     try:
      for child in psutil.Process(p.pid).children(recursive=True):child.kill()
     except psutil.NoSuchProcess:pass
     p.kill();p.wait();break
    time.sleep(.25)
  seconds=time.monotonic()-beg;total+=seconds;records.append({'stage':stage,'returncode':p.returncode,'wall_s':seconds,'peak_process_tree_RSS_bytes':peak,'log':str(log)});dump(OUT/'EXECUTION_PROCESSES.json',records);print('STAGE',stage,'rc',p.returncode,'s',seconds,flush=True)
  if p.returncode:
   dump(OUT/'EXECUTION_FAILURE.json',{'stage':stage,'returncode':p.returncode,'action':'STOP_NO_RETRY_OR_SCIENCE_REPAIR','wall_s':total});return
 dump(OUT/'BUDGET_AND_DELIVERY.json',{'wall_s_all_stages':total,'peak_process_tree_RSS_bytes':max(r['peak_process_tree_RSS_bytes'] for r in records),'storage':storage_guard(),'one_science_process':True,'BLAS_threads':1,'STOP':True});print('N2 SINGLE EXECUTION COMPLETE; STOP',flush=True)
if __name__=='__main__':main()