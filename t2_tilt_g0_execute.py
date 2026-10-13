"""Frozen single-compute-process once-only orchestration; enforce budget and scientific stop."""
from t2_tilt_g0_core import *
import subprocess,sys,psutil

def main():
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('once-only execution already started')
 freeze=read(OUT/'DESIGN_FREEZE.json')
 for name,s in freeze['code_and_input_SHA256'].items():assert sha(ROOT/name)==s,name
 head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==head and head!=PARENT
 dump(OUT/'EXECUTION_STARTED.json',{'utc':utc(),'design_SHA':head,'once_only':True});elapsed=0.;proc=[]
 stages=[('t2_tilt_g0_forward.py',),('t2_tilt_g0_generate.py','calibration'),('t2_tilt_g0_receiver.py','calibration'),('t2_tilt_g0_calibrate.py',),('t2_tilt_g0_generate.py','test'),('t2_tilt_g0_receiver.py','test'),('t2_tilt_g0_infer.py',),('t2_tilt_g0_evaluate.py',)]
 for stage in stages:
  began=time.monotonic();log=LOCAL/(stage[0]+('_'+stage[1] if len(stage)>1 else '')+'.log');peak=0
  with log.open('wb',buffering=0) as file:
   p=subprocess.Popen([sys.executable,'-B',*stage],cwd=ROOT,stdout=file,stderr=subprocess.STDOUT)
   while p.poll() is None:
    try:rss=psutil.Process(p.pid).memory_info().rss;peak=max(peak,rss)
    except psutil.NoSuchProcess:rss=0
    if rss>8*1024**3 or elapsed+time.monotonic()-began>21600:p.kill();p.wait();break
    time.sleep(.2)
  sec=time.monotonic()-began;elapsed+=sec;proc.append({'stage':list(stage),'returncode':p.returncode,'wall_s':sec,'peak_RSS_bytes':peak,'log':str(log)});dump(OUT/'EXECUTION_PROCESSES.json',proc);print('STAGE',stage,'rc',p.returncode,'seconds',sec,flush=True)
  if p.returncode:
   dump(OUT/'EXECUTION_FAILURE.json',{'stage':list(stage),'returncode':p.returncode,'action':'STOP_NO_SCIENTIFIC_RERUN','science_wall_s':elapsed});return
  if stage[0]=='t2_tilt_g0_forward.py' and not read(OUT/'FORWARD_DECISION.json')['admitted']:
   dump(OUT/'DECISION.json',{'stage':'T2-TILT-G0','status':'FORWARD_NUMERIC_NOT_ADMITTED','A':'NOT_EVALUATED','B':'NOT_EVALUATED','primary_geometry_counts':{'A':'NOT_EVALUATED','B':'NOT_EVALUATED','denominator':3},'calibration_records_executed':0,'test_resource_configs_executed':0,'planned_test_resource_configs':5184,'main_inference':'NOT_EVALUATED','propagation_status':'CONDITIONAL_MODAL_WINDOW_LIMIT','R4_percent':0,'research_lead_audit':'PENDING','state':'STOP','next_stage_authorized':False});print('FORWARD GATE BLOCKED; MAIN NOT_EVALUATED',flush=True);return
 print('REGISTERED TRIAL COMPLETE; cold review only next',flush=True)
if __name__=='__main__':main()