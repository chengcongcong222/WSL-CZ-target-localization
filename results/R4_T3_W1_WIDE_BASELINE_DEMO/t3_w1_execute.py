"""T3-W1 guarded one-main one-cold wrapper; no scientific retry."""
from pathlib import Path
import hashlib,json,subprocess,sys,time,psutil
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T3_W1')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def git(*a):return subprocess.check_output(['git','-c','safe.directory='+ROOT.as_posix(),*a],cwd=ROOT).decode().strip()
def verify():
    cfg=read(OUT/'T3_W1_DESIGN_FREEZE.json')
    for r in cfg['inputs']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    for r in cfg['frozen_code']:assert sha(OUT/r['path'])==r['sha256'],r['path']
    return cfg

def main():
    cfg=verify();assert not (OUT/'EXECUTION_STARTED.json').exists();assert not git('diff','--name-only') and not git('diff','--cached','--name-only');head=git('rev-parse','HEAD');remote=git('-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main').split()[0];assert head==remote and git('rev-parse','HEAD^')==cfg['parent_SHA']
    save(OUT/'EXECUTION_STARTED.json',{'design_SHA':head,'main_max_calls':1,'cold_max_calls':1,'scientific_retries':0});logs=LOCAL/'runtime_logs';logs.mkdir();start=time.monotonic();peak=0;failure=None;stages=[]
    for program,cap in [('t3_w1_main.py',5400),('t3_w1_cold.py',1800)]:
        verify();t=time.monotonic();last=t;log=logs/(program+'.log');stage={'program':program,'wall_limit_s':cap,'start_since_s':t-start};stages.append(stage)
        with log.open('wb') as handle:
            child=subprocess.Popen([sys.executable,'-B',str(OUT/program)],cwd=ROOT,stdout=handle,stderr=subprocess.STDOUT)
            while child.poll() is None:
                rss=psutil.Process().memory_info().rss
                for p in psutil.Process().children(recursive=True):
                    try:rss+=p.memory_info().rss
                    except psutil.NoSuchProcess:pass
                peak=max(peak,rss);now=time.monotonic();data=sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file());gitbytes=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file());reason=None
                if now-start>7200 or now-t>cap:reason='WALL_HARD_CAP'
                if rss>4*1024**3:reason='MEMORY_HARD_CAP'
                if data>2*1024**3 or gitbytes>64*1024**2:reason='STORAGE_HARD_CAP'
                if reason:child.kill();child.wait();failure=reason;break
                if now-last>25:
                    progress={'program':program,'total_wall_s':now-start,'stage_wall_s':now-t,'tree_RSS_bytes':rss,'local_D_bytes':data,'Git_result_bytes':gitbytes};save(OUT/'EXECUTION_PROGRESS.json',progress);print('HEARTBEAT',json.dumps(progress),flush=True);last=now
                time.sleep(.3)
        stage.update(returncode=child.returncode,wall_s=time.monotonic()-t,log=str(log),sha256=sha(log));print('COMPLETE',program,child.returncode,flush=True)
        if child.returncode!=0:failure=failure or 'SCIENTIFIC_PROCESS_FAILED_NO_RETRY'
        if program=='t3_w1_main.py':save(OUT/'MAIN_OUTPUT_SEAL.json',{'files':[{'path':str(p),'sha256':sha(p)} for p in OUT.iterdir() if p.is_file() and p.name not in ['MAIN_OUTPUT_SEAL.json','EXECUTION_PROGRESS.json']]})
        if failure:break
    save(OUT/'EXECUTION_BUDGET.json',{'main_calls':sum(r['program']=='t3_w1_main.py' for r in stages),'cold_calls':sum(r['program']=='t3_w1_cold.py' for r in stages),'scientific_retries':0,'wall_s':time.monotonic()-start,'peak_tree_RSS_bytes':peak,'new_D_bytes':sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file()),'Git_result_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),'failure':failure,'stages':stages,'acoustic_solver_T2_radial_feature_signal_MC_otherthan_W1':0})
    save(OUT/'STOP.json',{'stage':'R4-T3-W1','science_ended':True,'R4_percent':0,'depth':'NOT_EVALUATED','next_stage_authorized':False,'STOP':True});print('SCIENCE_STOP',failure or 'EXECUTED_ONCE',flush=True)
if __name__=='__main__':main()
