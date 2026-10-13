"""One execution pipeline with whole-stage wall/RSS monitoring and D-only logs."""
from t2_n3_common import *
import subprocess,sys,psutil

def main():
 if (LOCAL/'PIPELINE_STARTED.json').exists():raise RuntimeError('PIPELINE_ONCE_ONLY')
 dump(LOCAL/'PIPELINE_STARTED.json',{'utc':utc()});start=time.monotonic();peak=0;stages=[]
 for name in ('t2_n3_run.py','t2_n3_audit.py','t2_n3_report.py'):
  begin=time.monotonic();log=LOCAL/(Path(name).stem+'.pipeline.log')
  with log.open('wb') as f:
   p=subprocess.Popen([sys.executable,'-B',str(ROOT/name)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
   try:
    while p.poll() is None:
     process=psutil.Process();rss=process.memory_info().rss
     for child in process.children(recursive=True):
      try:rss+=child.memory_info().rss
      except psutil.NoSuchProcess:pass
     peak=max(peak,rss)
     if rss>8*1024**3 or time.monotonic()-start>7200:
      for child in process.children(recursive=True):
       try:child.kill()
       except psutil.NoSuchProcess:pass
      p.wait();raise RuntimeError('PIPELINE_HARD_BUDGET_STOP')
     time.sleep(.2)
   except BaseException:
    if p.poll() is None:p.kill();p.wait()
    raise
  stages.append({'script':name,'returncode':p.returncode,'wall_s':time.monotonic()-begin,'log':str(log),'sha256':sha(log)})
  print(name,'rc',p.returncode,'wall_s',stages[-1]['wall_s'],flush=True)
  if p.returncode!=0:break
 dump(OUT/'PIPELINE_BUDGET.json',{'stages':stages,'wall_s':time.monotonic()-start,'peak_tree_RSS_bytes':peak,'storage':guard(),'all_stages_completed':len(stages)==3 and all(r['returncode']==0 for r in stages),'science_retries':0})
if __name__=='__main__':main()