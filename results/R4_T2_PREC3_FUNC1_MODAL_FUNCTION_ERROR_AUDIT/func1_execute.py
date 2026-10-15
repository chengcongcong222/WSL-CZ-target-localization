"""One guarded archived-data analysis and one independent audit."""
from pathlib import Path
import json,hashlib,sys,subprocess,time,traceback,shutil,psutil
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC1')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def verify():
 cfg=read(OUT/'FUNC1_DESIGN_FREEZE.json')
 for r in cfg['inputs']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
 for r in cfg['vectors']:assert sha(r['path'])==r['sha256'],r['path']
 for r in cfg['frozen_code']:assert sha(OUT/r['path'])==r['sha256'],r['path']
def main():
 verify();assert not (OUT/'EXECUTION_STARTED.json').exists();head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip();remote=subprocess.check_output(['git','-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main'],cwd=ROOT).decode().split()[0];assert head==remote and head!='6b4bd489e7fd292adb64e25b72632948f20b2c4c';assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=ROOT).strip()
 write(OUT/'EXECUTION_STARTED.json',{'design_SHA':head,'main_max_calls':1,'cold_max_calls':1,'retries':0})
 start=time.monotonic();peak=0;stages=[];failure=None;logs=LOCAL/'runtime_logs';logs.mkdir(exist_ok=True)
 for prog,limit in [('func1_main.py',1200),('func1_cold.py',600)]:
  verify();t=time.monotonic();last=t;log=logs/(prog+'.log')
  try:
   with log.open('wb') as f:
    child=subprocess.Popen([sys.executable,'-B',str(OUT/prog)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT);stages.append({'program':prog,'wall_limit_s':limit})
    while child.poll() is None:
     rss=psutil.Process().memory_info().rss
     for p in psutil.Process().children(recursive=True):
      try:rss+=p.memory_info().rss
      except psutil.NoSuchProcess:pass
     peak=max(peak,rss);localbytes=sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file());gitbytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
     reason=None
     if time.monotonic()-t>limit or time.monotonic()-start>1800:reason='WALL_HARD_CAP'
     if rss>2*1024**3:reason='RSS_HARD_CAP'
     if localbytes>2*1024**3 or gitbytes>64*1024**2 or shutil.disk_usage('D:/').free<1024**3:reason='STORAGE_HARD_CAP'
     if reason:child.kill();child.wait();raise RuntimeError(reason)
     if time.monotonic()-last>30:
      progress={'stage':prog,'stage_wall_s':time.monotonic()-t,'total_wall_s':time.monotonic()-start,'RSS_bytes':rss,'new_D_bytes':localbytes,'git_result_bytes':gitbytes};write(OUT/'EXECUTION_PROGRESS.json',progress);print('HEARTBEAT',json.dumps(progress),flush=True);last=time.monotonic()
     time.sleep(.2)
   stages[-1].update(returncode=child.returncode,wall_s=time.monotonic()-t,log=str(log),sha256=sha(log));print(prog,'rc',child.returncode,flush=True)
   if child.returncode!=0:failure='PROGRAM_FAILED_NO_RETRY '+prog;break
   if prog=='func1_main.py':write(OUT/'MAIN_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in OUT.glob('FUNC1_*') if p.is_file()],'main_calls':1})
  except Exception as e:failure=str(e);(logs/'exception.txt').write_text(traceback.format_exc(),encoding='utf-8');break
 write(OUT/'EXECUTION_BUDGET.json',{'stages':stages,'failure':failure,'wall_s':time.monotonic()-start,'peak_tree_RSS_bytes':peak,'main_calls':sum(r['program']=='func1_main.py' for r in stages),'cold_calls':sum(r['program']=='func1_cold.py' for r in stages),'science_retries':0,'new_KRAKEN_KRAKENC_FIELD_FIELDS_SCOOTER_FE_FD_PRUFER_fullmode_MC_waveform_DFT_FIM_localization_calls':0,'new_D_bytes':sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file()),'git_result_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())})
 write(OUT/'SCIENCE_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in OUT.glob('FUNC1_*') if p.is_file()],'STOP':True});print('SCIENCE_STOP',failure or 'COMPLETED_ONCE',flush=True)
if __name__=='__main__':main()
