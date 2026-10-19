from common import *
import sys,subprocess,time,psutil,datetime
C=Path(__file__).resolve().parent

def run(unit,script,args=[]):
 manifest=read('CODE_FREEZE.json')
 for rel,h in manifest.items():
  if sha(O/rel)!=h:raise RuntimeError('FROZEN_CODE_IDENTITY:'+rel)
 folder=O/'units';folder.mkdir(exist_ok=True);mark=folder/(unit+'_SCIENCE_STARTED.json')
 if mark.exists():raise RuntimeError('NO_RETRY:'+unit)
 spec=read('CAP06_FROZEN_SPEC.json');budget=spec['budgets'];done=sum(json.loads(q.read_text()).get('elapsed_s',0) for q in folder.glob('*_ENDED.json'))
 mark.write_text(json.dumps(dict(unit=unit,code_SHA=manifest,config_SHA=sha(O/'CAP06_FROZEN_SPEC.json'),input_SHA=sha(O/'P0_SOURCE_IDENTITY.csv'),seed_contract=spec['seed_contract'],start_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),budget=budget,model_level=spec['evidence_scope']),indent=2))
 start=time.monotonic();peak=0;reason=None;log=W/'logs'/(unit+'.log')
 with log.open('w') as f:
  process=subprocess.Popen([sys.executable,str(C/script)]+args,stdout=f,stderr=subprocess.STDOUT,cwd=str(O.parents[1]))
  lastscan=0
  while process.poll() is None:
   try:
    tree=[psutil.Process(process.pid)]+psutil.Process(process.pid).children(recursive=True);memory=sum(q.memory_info().rss for q in tree);peak=max(peak,memory)
   except psutil.Error:pass
   elapsed=time.monotonic()-start
   if done+elapsed>budget['scientific_wall_s']:reason='WALL_BUDGET'
   if peak>budget['memory_bytes']:reason='MEMORY_BUDGET'
   if elapsed-lastscan>20:
    lastscan=elapsed;size=sum(q.stat().st_size for q in W.rglob('*') if q.is_file())
    if size>budget['D_bytes']:reason='D_STORAGE_BUDGET'
   if reason:
    for q in reversed(tree):
     try:q.kill()
     except psutil.Error:pass
    break
   time.sleep(.5)
  rc=process.wait()
 outcome=dict(unit=unit,elapsed_s=time.monotonic()-start,peak_process_tree_RSS=peak,exit_code=rc,status='COMPLETE' if rc==0 and reason is None else 'FAILED_NO_RETRY',reason=reason,log=str(log),log_SHA=sha(log))
 (folder/(unit+'_ENDED.json')).write_text(json.dumps(outcome,indent=2)+'\n');print(outcome,flush=True)
 if outcome['status']!='COMPLETE':return False
 return True

def main():
 phase=sys.argv[1];state=read('RUN_STATE.json');state['current_phase']=phase;state['status']='RUNNING';write('RUN_STATE.json',state)
 if phase=='P1':
  actions=[('P1_GENERATE','p1_toy.py',['generate']),('P1_EXTRACT','p1_toy.py',['extract']),('P1_EVALUATE','p1_toy.py',['evaluate']),('P1_COLD','audit_p1.py',[])]
 elif phase=='P2':actions=[('P2_GENERATE','generate_p2.py',[]),('P2_EXTRACT','extract_p2.py',[]),('P2_EVALUATE','evaluate_p2.py',[]),('P2_INFORMATION','information_p2.py',[]),('P2_COLD','audit_p2.py',[])]
 elif phase=='P3':
  audit=read('P2_INDEPENDENT_COLD_AUDIT.json');gate=read('P2_GATE.json')
  if not audit['PASS']:gate['P3_eligible']=[];gate['audit_rejection']=True;write('P2_GATE.json',gate)
  actions=[('P3_SEARCH','support_p3.py',[]),('P3_EVALUATE','evaluate_p3.py',[]),('P3_COLD','audit_p3.py',[])]
 else:raise ValueError(phase)
 ok=True
 for unit,script,args in actions:
  if not run(unit,script,args):ok=False;break
 state['phases'][phase]='EXECUTED' if ok else 'FAILED_NO_RETRY';state['status']='CHECKPOINT_PENDING';write('RUN_STATE.json',state)
 (O/'SESSION_HANDOFF.md').write_text('# CAP06 recovery\nRead RUN_STATE.json and all units/* markers. Never rerun a started unit. Current phase '+phase+'; executed='+str(ok)+'. Preserve failures. Independent other branch may continue. No acoustic solvers.\n')
if __name__=='__main__':main()
