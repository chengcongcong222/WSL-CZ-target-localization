"""PREC1 operational envelope: original supplied scientific programs stay byte-exact."""
from pathlib import Path
import json,hashlib,subprocess,sys,time,os,datetime,traceback,psutil
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_T2_N3_PREC1_MOD_PRECISION_ATTRIBUTION'
LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_N3_PREC1')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def frozen():
 for r in json.loads((OUT/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))['files']:
  if sha(r['path'])!=r['sha256']:raise RuntimeError('FROZEN_IDENTITY_MISMATCH '+r['path'])
def main():
 if (OUT/'FORMAL_PIPELINE_STARTED.json').exists():raise RuntimeError('PIPELINE_ONCE_ONLY')
 frozen();head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip();remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],cwd=ROOT).decode().split()[0]
 assert head==remote and head!='8ee2bd53fe1b7a50725833e6c7b5c5dea35fbe23'
 assert not subprocess.check_output(['git','status','--porcelain=v1'],cwd=ROOT).strip()
 write(OUT/'FORMAL_PIPELINE_STARTED.json',{'design_SHA':head,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main_execution_limit':1,'cold_execution_limit':1})
 start=time.monotonic();peak=0;stages=[];failure=None;logs=LOCAL/'logs';logs.mkdir(exist_ok=True)
 for name,args in [('t2_n3_prec1.py',['--repo',str(ROOT),'--out',str(OUT)]),('t2_n3_prec1_cold_audit.py',['--repo',str(ROOT),'--out',str(OUT)])]:
  frozen();log=logs/(Path(name).stem+'.log');begin=time.monotonic()
  try:
   with log.open('wb') as f:
    p=subprocess.Popen([sys.executable,'-B',str(ROOT/name)]+args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    while p.poll() is None:
     process=psutil.Process();rss=process.memory_info().rss
     for child in process.children(recursive=True):
      try:rss+=child.memory_info().rss
      except psutil.NoSuchProcess:pass
     peak=max(peak,rss);resultbytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())+sum(p.stat().st_size for p in logs.rglob('*') if p.is_file())
     if time.monotonic()-start>600 or rss>2*1024**3 or resultbytes>1024**3:
      for child in process.children(recursive=True):
       try:child.kill()
       except psutil.NoSuchProcess:pass
      p.wait();raise RuntimeError('HARD_BUDGET_STOP')
     time.sleep(.1)
   stages.append({'program':name,'returncode':p.returncode,'wall_s':time.monotonic()-begin,'log':str(log),'sha256':sha(log)})
   print(name,'rc',p.returncode,'wall_s',stages[-1]['wall_s'],flush=True)
   if p.returncode!=0:failure='PROGRAM_FAILED_NO_RETRY '+name;break
   if name=='t2_n3_prec1.py':
    files=[p for p in OUT.glob('PREC1_*') if p.is_file()]
    write(OUT/'MAIN_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in files],'main_runs':1,'no_cold_run_yet':True})
  except Exception as e:
   failure=str(e);(logs/'exception.txt').write_text(traceback.format_exc(),encoding='utf-8');break
 write(OUT/'EXECUTION_BUDGET.json',{'stages':stages,'failure':failure,'wall_s':time.monotonic()-start,'peak_tree_RSS_bytes':peak,'result_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),'formal_main_calls':sum(r['program']=='t2_n3_prec1.py' for r in stages),'formal_cold_calls':sum(r['program']=='t2_n3_prec1_cold_audit.py' for r in stages),'new_KRAKEN_SCOOTER_FIELD_FIELDS_MC_waveforms_DFT_inference':0,'retries':0})
 files=[p for p in OUT.glob('PREC1_*') if p.is_file()]
 write(OUT/'SCIENTIFIC_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in files],'failure':failure,'final':True})
 if failure:print('STOP',failure,flush=True)
 else:print('ONE MAIN + ONE COLD COMPLETE; scientific outputs sealed',flush=True)
if __name__=='__main__':main()