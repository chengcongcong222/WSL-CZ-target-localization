"""PREC2 one-shot serial envelope; byte-exact task science programs are not patched."""
from pathlib import Path
import subprocess,sys,json,hashlib,time,datetime,traceback,shutil,psutil
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'results/R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY';LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC2');PARENT='ac43a8702282fc6c80663319cb6757a30b084e09'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def verify():
 for r in read(OUT/'DESIGN_FREEZE.json')['files']:
  if sha(r['path'])!=r['sha256']:raise RuntimeError('FROZEN_IDENTITY_MISMATCH '+r['path'])
def main():
 if (OUT/'PIPELINE_STARTED.json').exists():raise RuntimeError('ONCE_ONLY_GUARD')
 verify();head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip();remote=subprocess.check_output(['git','-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main'],cwd=ROOT).decode().split()[0]
 assert head==remote and head!=PARENT
 assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=ROOT).strip()
 write(OUT/'PIPELINE_STARTED.json',{'design_SHA':head,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main_limit':1,'cold_limit':1,'science_processes':1})
 started=time.monotonic();peak=0;stages=[];failure=None;logs=LOCAL/'science_logs';logs.mkdir(exist_ok=True);config=OUT/'T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json'
 for name,limit in [('t2_prec2_readonly.py',1800),('t2_prec2_cold_audit.py',900)]:
  verify();began=time.monotonic();log=logs/(Path(name).stem+'.log');last_heartbeat=began
  try:
   with log.open('wb') as f:
    proc=subprocess.Popen([sys.executable,'-B',str(ROOT/name),'--repo',str(ROOT),'--config',str(config),'--out',str(OUT)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    stages.append({'program':name,'limit_s':limit,'status':'STARTED','log':str(log)});write(OUT/'EXECUTION_PROGRESS.json',{'stages':stages,'total_wall_s':time.monotonic()-started})
    while proc.poll() is None:
     current=psutil.Process();rss=current.memory_info().rss
     for child in current.children(recursive=True):
      try:rss+=child.memory_info().rss
      except psutil.NoSuchProcess:pass
     peak=max(peak,rss);outbytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
     reason=None
     if time.monotonic()-began>limit:reason='STAGE_WALL_HARD_CAP'
     if time.monotonic()-started>2700:reason='PIPELINE_WALL_HARD_CAP'
     if rss>4*1024**3:reason='RSS_HARD_CAP'
     if outbytes>64*1024**2:reason='GIT_RESULT_BYTES_HARD_CAP'
     if shutil.disk_usage('D:/').free<=1024**3:reason='D_RESERVE_HARD_CAP'
     if reason:
      for child in current.children(recursive=True):
       try:child.kill()
       except psutil.NoSuchProcess:pass
      proc.wait();raise RuntimeError(reason)
     if time.monotonic()-last_heartbeat>=30:
      progress={'stage':name,'stage_wall_s':time.monotonic()-began,'total_wall_s':time.monotonic()-started,'RSS_bytes':rss,'log_bytes':log.stat().st_size};write(OUT/'EXECUTION_PROGRESS.json',progress);print('HEARTBEAT',json.dumps(progress),flush=True);last_heartbeat=time.monotonic()
     time.sleep(.2)
   stages[-1].update(status='SUCCESS' if proc.returncode==0 else 'FAILED',returncode=proc.returncode,wall_s=time.monotonic()-began,sha256=sha(log));print(name,'rc',proc.returncode,'wall',stages[-1]['wall_s'],flush=True)
   if proc.returncode!=0:failure='PROGRAM_FAILED_NO_RETRY '+name;break
   if name=='t2_prec2_readonly.py':
    write(OUT/'MAIN_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in OUT.glob('PREC2_*') if p.is_file()],'main_calls':1})
    # A false preregistered FDM phase-budget flag stops further science per lead instruction.
    if not read(OUT/'PREC2_DECISION.json')['FDM_R4_gap_all_within_target']:
     failure='FDM_REGISTERED_INTERNAL_CONVERGENCE_NOT_ESTABLISHED';break
  except Exception as e:
   failure=str(e);(logs/'exception.txt').write_text(traceback.format_exc(),encoding='utf-8');
   if stages:stages[-1].update(status='STOPPED',reason=failure,wall_s=time.monotonic()-began,sha256=sha(log))
   break
 budget={'stages':stages,'failure':failure,'wall_s':time.monotonic()-started,'peak_tree_RSS_bytes':peak,'new_result_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),'main_calls':sum(r['program']=='t2_prec2_readonly.py' for r in stages),'cold_calls':sum(r['program']=='t2_prec2_cold_audit.py' for r in stages),'science_retries':0,'new_acoustic_solver_MC_waveform_DFT_FIM_inference_calls':0}
 write(OUT/'EXECUTION_BUDGET.json',budget)
 write(OUT/'SCIENTIFIC_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in OUT.glob('PREC2_*') if p.is_file()],'failure':failure,'final_science_STOP':True})
 print('SCIENCE_STOP',failure or 'BOUNDED_MAIN_AND_COLD_COMPLETE',flush=True)
if __name__=='__main__':main()