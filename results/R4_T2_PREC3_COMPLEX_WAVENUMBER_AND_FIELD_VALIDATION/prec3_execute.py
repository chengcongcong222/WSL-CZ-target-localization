from pathlib import Path
import json,hashlib,subprocess,time,sys,psutil,shutil,traceback
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def verify():
 cfg=read(OUT/'PREC3_DESIGN_FREEZE.json')
 for r in cfg['inputs']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
 for r in cfg['frozen_code']:assert sha(OUT/r['path'])==r['sha256'],r['path']
 for r in cfg['external_sources']:assert sha(r['path'])==r['sha256'],r['path']
def main():
 verify();assert not (OUT/'EXECUTION_STARTED.json').exists()
 head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip();remote=subprocess.check_output(['git','-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main'],cwd=ROOT).decode().split()[0]
 assert head==remote and head!='d5c0da07c17b12c6e169713a23caf4f2bcfd317c';assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=ROOT).strip()
 write(OUT/'EXECUTION_STARTED.json',{'design_SHA':head,'main_max_calls':1,'cold_max_calls':1,'scientific_retries':0})
 start=time.monotonic();peak=0;stages=[];fail=None;logs=LOCAL/'logs';logs.mkdir(exist_ok=True)
 for prog,limit in [('prec3_main.py',3000),('prec3_cold.py',600)]:
  verify();t=time.monotonic();last=t;p=logs/(prog+'.log')
  try:
   with p.open('wb') as f:
    child=subprocess.Popen([sys.executable,'-B',str(OUT/prog)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT);stages.append({'program':prog,'wall_limit_s':limit})
    while child.poll() is None:
     rss=psutil.Process().memory_info().rss
     for q in psutil.Process().children(recursive=True):
      try:rss+=q.memory_info().rss
      except psutil.NoSuchProcess:pass
     peak=max(peak,rss)
     local=sum(q.stat().st_size for q in LOCAL.rglob('*') if q.is_file());gitbytes=sum(q.stat().st_size for q in OUT.rglob('*') if q.is_file())
     reason=None
     if time.monotonic()-t>limit or time.monotonic()-start>3600:reason='WALL_HARD_CAP'
     if rss>4*1024**3:reason='RSS_HARD_CAP'
     if local>2*1024**3 or gitbytes>64*1024**2 or shutil.disk_usage('D:/').free<1024**3:reason='STORAGE_HARD_CAP'
     if reason:
      child.kill();child.wait();raise RuntimeError(reason)
     if time.monotonic()-last>30:
      data={'stage':prog,'stage_wall_s':time.monotonic()-t,'total_wall_s':time.monotonic()-start,'RSS_bytes':rss,'local_bytes':local,'git_result_bytes':gitbytes};write(OUT/'EXECUTION_PROGRESS.json',data);print('HEARTBEAT',json.dumps(data),flush=True);last=time.monotonic()
     time.sleep(.3)
   stages[-1].update(returncode=child.returncode,wall_s=time.monotonic()-t,log=str(p),sha256=sha(p));print(prog,'rc',child.returncode,flush=True)
   if child.returncode!=0:fail='PROGRAM_EXCEPTION_NO_RETRY '+prog;break
   if prog=='prec3_main.py':
    write(OUT/'MAIN_SEAL.json',{'files':[{'path':str(q),'sha256':sha(q)} for q in OUT.glob('PREC3_*') if q.is_file()],'main_calls':1})
    # Main scientific Gate failures cancel B/C there; this separate call is only cold audit.
  except Exception as e:fail=str(e);(logs/'exception.txt').write_text(traceback.format_exc(),encoding='utf-8');break
 write(OUT/'EXECUTION_BUDGET.json',{'stages':stages,'failure':fail,'wall_s':time.monotonic()-start,'peak_tree_RSS_bytes':peak,'main_calls':sum(s['program']=='prec3_main.py' for s in stages),'cold_calls':sum(s['program']=='prec3_cold.py' for s in stages),'retries':0,'new_named_acoustic_solver_MC_waveform_DFT_FIM_localization_calls':0,'local_bytes':sum(q.stat().st_size for q in LOCAL.rglob('*') if q.is_file()),'git_result_bytes':sum(q.stat().st_size for q in OUT.rglob('*') if q.is_file())})
 write(OUT/'FINAL_SCIENCE_SEAL.json',{'files':[{'path':str(q),'sha256':sha(q)} for q in OUT.glob('PREC3_*') if q.is_file()],'STOP':True});print('SCIENCE_STOP',fail or 'COMPLETED_FROZEN_EXECUTION',flush=True)
if __name__=='__main__':main()
