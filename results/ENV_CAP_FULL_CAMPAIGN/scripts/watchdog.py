"""Process-tree time/RSS watchdog; resource guard, not scientific code."""
import os,sys,time,json,subprocess
from pathlib import Path
import psutil
script=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]);limit=float(sys.argv[3]);marker=out.with_suffix('.SCIENCE_STARTED.json')
if marker.exists():raise SystemExit('Repeated campaign process refused')
marker.write_text(json.dumps({'script':str(script),'time_unix':time.time(),'timeout_s':limit,'pid_supervisor':os.getpid()}),encoding='utf-8')
start=time.monotonic();peak=0;reason=None
with out.open('wb') as f:
 p=subprocess.Popen([sys.executable,str(script)],stdout=f,stderr=subprocess.STDOUT);parent=psutil.Process(p.pid)
 while p.poll() is None:
  try:
   tree=[parent]+parent.children(recursive=True);rss=sum(q.memory_info().rss for q in tree if q.is_running());peak=max(peak,rss)
   if rss>8*1024**3:reason='RSS_BUDGET_EXCEEDED'
   if time.monotonic()-start>limit:reason='WALL_BUDGET_EXCEEDED'
   if reason:
    for q in reversed(tree):
     try:q.kill()
     except psutil.Error:pass
    break
  except psutil.Error:pass
  time.sleep(.2)
 rc=p.wait()
v={'process_PID':p.pid,'wall_s':time.monotonic()-start,'peak_tree_RSS_bytes':peak,'returncode':rc,'budget_stop':reason,'script':str(script)}
out.with_suffix('.JOB_EXIT.json').write_text(json.dumps(v,indent=2),encoding='utf-8');print(json.dumps(v));sys.exit(0 if rc==0 and reason is None else 2)
