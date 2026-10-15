"""Hard-cap wrapper: one main and one independent cold, no retries."""
from pathlib import Path
import json,hashlib,subprocess,sys,time,psutil,shutil,traceback
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def git(*args):return subprocess.check_output(['git','-c','safe.directory='+ROOT.as_posix(),*args],cwd=ROOT).decode().strip()
def verify():
 cfg=read(OUT/'FUNC2_DESIGN_FREEZE.json')
 for r in cfg['inputs']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
 for r in cfg['vectors']:assert sha(r['path'])==r['sha256'],r['path']
 for r in cfg['frozen_code']:assert sha(OUT/r['path'])==r['sha256'],r['path']
def main():
 verify();assert not (OUT/'EXECUTION_STARTED.json').exists();assert not git('diff','--name-only');assert not git('diff','--cached','--name-only');head=git('rev-parse','HEAD');remote=git('-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main').split()[0];assert head==remote and head!=read(OUT/'FUNC2_DESIGN_FREEZE.json')['parent_SHA']
 save(OUT/'EXECUTION_STARTED.json',{'design_SHA':head,'main_calls_max':1,'cold_calls_max':1,'retries':0});logs=LOCAL/'runtime_logs';logs.mkdir();start=time.monotonic();peak=0;stages=[];fail=None
 for name,cap in [('func2_main.py',2900),('func2_cold.py',700)]:
  if time.monotonic()-start>=3600:fail='TOTAL_WALL_BUDGET';break
  verify();t=time.monotonic();last=t;log=logs/(name+'.log');stage={'program':name,'start_since_s':t-start,'wall_limit_s':cap};stages.append(stage)
  with log.open('wb') as handle:
   child=subprocess.Popen([sys.executable,'-B',str(OUT/name)],cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT)
   while child.poll() is None:
    rss=psutil.Process().memory_info().rss
    for p in psutil.Process().children(recursive=True):
     try:rss+=p.memory_info().rss
     except psutil.NoSuchProcess:pass
    peak=max(peak,rss);data=sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file());results=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file());reason=None
    if time.monotonic()-start>3600 or time.monotonic()-t>cap:reason='WALL_HARD_CAP'
    if rss>4*1024**3:reason='MEMORY_HARD_CAP'
    if data>2*1024**3 or results>64*1024**2:reason='STORAGE_HARD_CAP'
    if reason:
     child.kill();child.wait();fail=reason;break
    if time.monotonic()-last>25:
     progress={'program':name,'science_wall_s':time.monotonic()-start,'stage_wall_s':time.monotonic()-t,'tree_RSS_bytes':rss,'new_D_data_bytes':data,'git_bytes':results};save(OUT/'EXECUTION_PROGRESS.json',progress);print('HEARTBEAT',json.dumps(progress),flush=True);last=time.monotonic()
    time.sleep(.2)
  stage.update(returncode=child.returncode,wall_s=time.monotonic()-t,log=str(log),sha256=sha(log))
  print('FINISHED',name,child.returncode,flush=True)
  if child.returncode!=0:fail=fail or 'PROCESS_EXCEPTION_NO_RETRY'
  if name=='func2_main.py':save(OUT/'MAIN_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in OUT.glob('FUNC2_*') if p.is_file()],'raw_index':read(OUT/'FUNC2_RAW_DATA_INDEX.json') if (OUT/'FUNC2_RAW_DATA_INDEX.json').exists() else None})
  if fail:break
 save(OUT/'EXECUTION_BUDGET.json',{'wall_s':time.monotonic()-start,'peak_tree_RSS_bytes':peak,'main_calls':sum(s['program']=='func2_main.py' for s in stages),'cold_calls':sum(s['program']=='func2_cold.py' for s in stages),'retries':0,'failure':fail,'stages':stages,'new_D_bytes':sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file()),'git_result_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),'new_external_acoustic_solvers':0,'new_old_FE_FD_Pruefer_fullmode_MC_waveform_DFT_FIM_localization':0})
 save(OUT/'STOP.json',{'stage':'T2-PREC3-FUNC2','science_execution_ended':True,'automatic_next_stage':False,'R4_percent':0,'STOP':True});print('STOP',fail or 'EXECUTED_ONCE',flush=True)
if __name__=='__main__':main()
