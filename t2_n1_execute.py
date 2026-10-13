"""Once-only archived analysis and read-only cold audit, shared hard budget, single numeric worker."""
from t2_n1_common import *
import subprocess,sys,psutil

def main():
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('EXECUTION_ONCE_GUARD')
 f=read(OUT/'N1_DESIGN_FREEZE.json')
 for path,h in f['code_and_input_SHA256'].items():assert sha(ROOT/path)==h,path
 verify_identities();head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0];assert head==remote and head!=PARENT
 dump(OUT/'EXECUTION_STARTED.json',{'utc':utc(),'design_SHA':head,'only_one_registered_analysis':True,'new_solver_or_observation':0});processes=[];total=0.
 for script in ['t2_n1_analyze.py','t2_n1_audit.py','t2_n1_report.py']:
  started=time.monotonic();peak=0;log=LOCAL/(script+'.log')
  with log.open('wb',buffering=0) as file:
   p=subprocess.Popen([sys.executable,'-B',script],cwd=ROOT,stdout=file,stderr=subprocess.STDOUT)
   while p.poll() is None:
    try:rss=psutil.Process(p.pid).memory_info().rss;peak=max(peak,rss)
    except psutil.NoSuchProcess:rss=0
    if total+time.monotonic()-started>1800 or rss>4*1024**3:p.kill();p.wait();break
    time.sleep(.1)
  sec=time.monotonic()-started;total+=sec;processes.append({'stage':script,'returncode':p.returncode,'wall_s':sec,'sampled_peak_RSS_bytes':peak,'log':str(log)});dump(OUT/'EXECUTION_PROCESSES.json',processes);print('STAGE',script,'rc',p.returncode,'seconds',sec,flush=True)
  if p.returncode:
   dump(OUT/'EXECUTION_FAILURE.json',{'stage':script,'rc':p.returncode,'science_wall_s':total,'action':'STOP_NO_SCIENCE_RERUN_OR_CORE_CHANGE'})
   if not (OUT/'DECISION.json').exists():dump(OUT/'DECISION.json',{'stage':'T2-N1','status':'PARTIAL_BUDGET_OR_IMPLEMENTATION_FAILURE','failed_stage':script,'HISTORICAL_G0_FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED':True,'A_NOT_EVALUATED':True,'B_NOT_EVALUATED':True,'CHIGH_UNCERTIFIED':True,'R4_percent':0,'state':'STOP'})
   return
 dump(OUT/'BUDGET_AND_STOP.json',{'all_stages_wall_s':total,'sampled_peak_RSS_bytes':max(r['sampled_peak_RSS_bytes'] for r in processes),'storage':budget(),'new_KRAKEN_FIELD_propagation':0,'new_MC_waveform_DFT_spectra_optimization_FIM':0,'numeric_worker_processes':1,'BLAS_threads':1,'STOP':True});print('READ-ONLY N1 COMPLETE; STOP',flush=True)
if __name__=='__main__':main()