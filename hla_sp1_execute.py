"""Single once-only execution launcher; no concurrent numeric workers."""
import subprocess,sys,time,os,json
from pathlib import Path
import psutil
from hla_sp1_core import ROOT,OUT,LOCAL,read,dump,sha,stamp,guard

def main():
 if (OUT/'EXECUTION_STARTED.json').exists():raise RuntimeError('no rerun')
 cfg=read(OUT/'DESIGN_FREEZE.json')
 for n,s in cfg['code_sha256'].items():assert sha(ROOT/n)==s,n
 for n,s in cfg['input_sha256'].items():assert sha(OUT/n)==s,n
 head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip();assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==head
 dump(OUT/'EXECUTION_STARTED.json',{'time_utc':stamp(),'design_SHA':head,'once_only':True});log=LOCAL/'EXECUTION.log';processes=[]
 stages=[('hla_sp1_generate.py','calibration'),('hla_sp1_receiver.py','calibration'),('hla_sp1_calibrate.py',),('hla_sp1_generate.py','test'),('hla_sp1_receiver.py','test'),('hla_sp1_infer.py',),('hla_sp1_evaluate.py',)]
 with log.open('ab',buffering=0) as file:
  for stage in stages:
   began=time.monotonic();elapsed,_=guard();p=subprocess.Popen([sys.executable,'-B',*stage],cwd=ROOT,stdout=file,stderr=subprocess.STDOUT);peak=0
   while p.poll() is None:
    try:rss=psutil.Process(p.pid).memory_info().rss;peak=max(peak,rss)
    except psutil.NoSuchProcess:rss=0
    if rss>8*1024**3 or time.monotonic()-began>14400-elapsed:p.kill();p.wait();raise RuntimeError('hard budget stop '+str(stage))
    time.sleep(.25)
   processes.append({'stage':list(stage),'returncode':p.returncode,'wall_s':time.monotonic()-began,'sampled_peak_RSS':peak});dump(OUT/'EXECUTION_PROCESSES.json',processes);print('STAGE',stage,'exit',p.returncode,'s',time.monotonic()-began,flush=True)
   if p.returncode:dump(OUT/'EXECUTION_FAILURE.json',{'stage':list(stage),'returncode':p.returncode,'action':'STOP; no repair/rerun same panel'});return
 print('REGISTERED EXECUTION COMPLETE',flush=True)
if __name__=='__main__':main()
