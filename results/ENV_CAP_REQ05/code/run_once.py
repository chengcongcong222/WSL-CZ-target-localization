"""Bounded single-execution launcher, writes only to D:, no scientific retry."""
import sys, os, subprocess, pathlib, time, json, psutil, hashlib
O=pathlib.Path(__file__).resolve().parents[1];W=pathlib.Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_REQ05')
phase=sys.argv[1];assert phase in ['analysis','audit']
freeze=json.loads((O/'DESIGN_FREEZE.json').read_text())
for row in freeze['code_and_spec']:
 assert hashlib.sha256((O/row['path']).read_bytes()).hexdigest()==row['sha256'],row['path']
marker=W/(phase+'.started');marker.open('x').write(str(time.time()))
t=time.monotonic();peak=0;previous=0
if phase=='audit':
 prior=json.loads((O/'ANALYSIS_RUN.json').read_text());assert prior['exit_code']==0
 previous=prior['wall_s']
child=subprocess.Popen([sys.executable,str(O/'code'/('analyze.py' if phase=='analysis' else 'cold_audit.py'))],env=os.environ.copy())
reason=None
while child.poll() is None:
 try:
  p=psutil.Process(child.pid);rss=sum(q.memory_info().rss for q in [p]+p.children(recursive=True) if q.is_running());peak=max(peak,rss)
 except psutil.Error:pass
 if time.monotonic()-t+previous>600 or peak>4*1024**3:
  reason='BUDGET_EXCEEDED'
  for q in psutil.Process(child.pid).children(recursive=True):q.kill()
  child.kill();break
 time.sleep(.05)
rc=child.wait();rec=dict(phase=phase,wall_s=time.monotonic()-t,peak_sampled_tree_RSS=peak,exit_code=rc,stop_reason=reason,acoustic_calls=0,MC=0)
(O/(phase.upper()+'_RUN.json')).write_text(json.dumps(rec,indent=2)+'\n')
print(json.dumps(rec));sys.exit(rc or bool(reason))
