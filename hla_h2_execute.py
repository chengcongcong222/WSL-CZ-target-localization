"""Single frozen observation-only kinematic estimation execution."""
import time,json,gzip,csv,traceback,subprocess
import numpy as np,pandas as pd,psutil
from threadpoolctl import threadpool_limits
import hla_h2_core as c
FIELDS=['terminal_id','config_id','geometry','replicate','method','sigma','branch_id','start_id','start_source',*c.AXES,'b_original','profile_b','J_beta','J_branch','J_folded','J_profile','accepted_original','accepted_profile','solver_success','solver_status','nfev','njev','geometry_calls','optimality','bound_touch','exception','init_rank','init_columns','init_condition','init_clipped','init_nonpositive_rho0','init_rho0_linear_m','init_bias_clipped','init_r_km','init_theta_deg','init_v_mps','init_psi_deg','init_b']
GLOBAL_START=None
PEAK=0
def guard():
    global PEAK
    rss=psutil.Process().memory_info().rss;PEAK=max(PEAK,rss)
    if time.time()-GLOBAL_START>12600:raise RuntimeError('HARD_WALL_BUDGET')
    if rss>8*1024**3:raise RuntimeError('HARD_MEMORY_BUDGET')
def plan():
    rows=[]
    for g in range(8):
        for rep in range(32):
            variants=[('B',0.,-1)]
            for si,sigma in enumerate(c.SIGMAS):variants.extend([('S',sigma,si),('U',sigma,si)])
            variants.extend([('U-mismatch',.05,1),('UB-matched',.05,1)])
            for method,sigma,si in variants:
                rows.append(dict(config_id=len(rows),geometry=g,replicate=rep,method=method,sigma=sigma,sigma_index=si))
    return pd.DataFrame(rows)
def load_observation(row,visible):
    g,rep,si=int(row.geometry),int(row.replicate),int(row.sigma_index);method=row.method
    beta=visible['beta'][g,rep].copy()
    if method=='B':return c.Observation(beta,None,None,'B')
    if method=='S':
        # Dedicated control file only accessed for S; never supplied to unsigned API.
        with np.load(c.OUT/'SIGNED_CONTROL_OBSERVATIONS.npz') as data:radial=data['signed'][g,rep,si].copy()
    elif method=='U':radial=visible['unsigned'][g,rep,si].copy()
    else:radial=visible['unsigned_biased'][g,rep].copy()
    return c.Observation(beta,radial,float(row.sigma),method)
def choose(terminals):
    eligible=[v for v in terminals if not v.get('exception') and v['accepted_profile']]
    if not eligible:return None
    return min(eligible,key=lambda v:(v['J_beta']+v['J_profile'],v['terminal_id']))
def fit_configuration(row,obs,start_terminal_id,writer):
    records=[];terminal_id=start_terminal_id
    if obs.kind=='B':
        branches=[(-1,None,c.bearing_starts(obs.beta))]
    else:
        signs=[(0,None)] if obs.kind=='S' else list(enumerate(c.SIGNS))
        branches=[]
        for branch_id,sign in signs:
            unfolded=obs.radial.copy() if obs.kind=='S' else sign*obs.radial
            linear,details=c.linear_radial_start(obs.beta,unfolded,obs.kind=='UB-matched')
            center=np.r_[c.CENTER,0.] if obs.kind=='UB-matched' else c.CENTER.copy()
            center_meta=dict(rank=-1,columns=6 if obs.kind=='UB-matched' else 5,condition=np.nan,clipped=False)
            branches.append((branch_id,unfolded,[('LINEAR_GEOMETRY',linear,details),('PHYSICAL_CENTER',center,center_meta)]))
    for branch_id,unfolded,starts in branches:
        for start_id,(source,initial,meta) in enumerate(starts):
            guard();record={k:np.nan for k in FIELDS}
            record.update(terminal_id=terminal_id,config_id=int(row.config_id),geometry=int(row.geometry),replicate=int(row.replicate),method=row.method,sigma=float(row.sigma),branch_id=branch_id,start_id=start_id,start_source=source,exception='',init_rank=meta.get('rank',-1),init_columns=meta.get('columns',-1),init_condition=meta.get('condition',np.nan),init_clipped=meta.get('clipped',False),init_nonpositive_rho0=meta.get('nonpositive_rho0',False),init_rho0_linear_m=meta.get('rho0_linear_m',np.nan),init_bias_clipped=meta.get('bias_clipped',False),init_b=float(initial[4]) if len(initial)==5 else 0.)
            record.update(dict(zip(['init_'+ax for ax in c.AXES],initial[:4])))
            try:
                state,b,fit,calls=c.fit_start(obs,unfolded,initial)
                if not np.isfinite(state).all() or not np.isfinite(fit.fun).all():raise ValueError('Nonfinite solver terminal')
                scores=c.physical_scores(obs,state,b,unfolded)
                record.update(dict(zip(c.AXES,state)))
                record.update(scores,b_original=b,solver_success=bool(fit.success),solver_status=int(fit.status),nfev=int(fit.nfev),njev=int(fit.njev or 0),geometry_calls=calls,optimality=float(fit.optimality),bound_touch=bool(np.any((state-c.LOW)<1e-6*c.WIDTH)|np.any((c.HIGH-state)<1e-6*c.WIDTH)|(obs.kind=='UB-matched' and min(b+.2,.2-b)<4e-7)))
            except Exception as error:
                record.update(exception=type(error).__name__+': '+str(error),accepted_original=False,accepted_profile=False,solver_success=False,solver_status=-99)
            writer.writerow(record);records.append(record);terminal_id+=1
    selected=choose(records)
    base=dict(config_id=int(row.config_id),geometry=int(row.geometry),replicate=int(row.replicate),method=row.method,sigma=float(row.sigma),status='ACCEPTED' if selected else 'ABSTAIN',terminal_runs=len(records),distinct_expansion_branches=len(set(v['branch_id'] for v in records)),accepted_terminal_records=sum(bool(v['accepted_profile']) for v in records),solver_exceptions=sum(bool(v['exception']) for v in records),solver_nonconverged=sum(not v['solver_success'] and not v['exception'] for v in records))
    if selected:
        base.update(dict(zip(c.AXES,[selected[ax] for ax in c.AXES])),b=selected['profile_b'],selected_terminal_id=selected['terminal_id'],J_beta=selected['J_beta'],J_radial=selected['J_profile'])
    else:
        base.update({ax:np.nan for ax in c.AXES});base.update(b=np.nan,selected_terminal_id=-1,J_beta=np.nan,J_radial=np.nan)
    return base,terminal_id
def execute():
    global GLOBAL_START
    freeze=json.loads((c.OUT/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
    for name,expected in freeze['sha256'].items():
        if c.sha(c.ROOT/name)!=expected:raise RuntimeError('Frozen bytes changed: '+name)
    head=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()
    remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]
    if head!=remote:raise RuntimeError('Design local/remote mismatch')
    if (c.OUT/'EXECUTION_MANIFEST.json').exists() or (c.OUT/'ALL_TERMINALS.csv.gz').exists():raise RuntimeError('One execution only')
    GLOBAL_START=time.time()
    c.write_json(c.OUT/'RUN_CLOCK.json',dict(epoch=GLOBAL_START,design_sha=head,remote_main_at_start=remote))
    import hla_h2_generate
    hla_h2_generate.generate()
    configs=pd.read_csv(c.OUT/'CONFIGURATION_PLAN.csv');estimates=[];terminal_id=0
    with np.load(c.OUT/'OBSERVATIONS.npz') as visible,gzip.open(c.OUT/'ALL_TERMINALS.csv.gz','wt',encoding='utf-8',newline='') as stream,(c.OUT/'ESTIMATES.csv').open('w',encoding='utf-8',newline='') as estimates_stream:
        writer=csv.DictWriter(stream,fieldnames=FIELDS);writer.writeheader()
        estimate_writer=None
        for _,row in configs.iterrows():
            guard();obs=load_observation(row,visible)
            estimate,terminal_id=fit_configuration(row,obs,terminal_id,writer);estimates.append(estimate)
            if estimate_writer is None:
                estimate_writer=csv.DictWriter(estimates_stream,fieldnames=list(estimate));estimate_writer.writeheader()
            estimate_writer.writerow(estimate);estimates_stream.flush()
            if int(row.config_id)%11==10:
                stream.flush()
                c.write_json(c.OUT/'EXECUTION_CHECKPOINT.json',dict(completed_configurations=len(estimates),terminal_runs=terminal_id,elapsed_seconds=time.time()-GLOBAL_START,peak_rss_bytes=PEAK))
                if int(row.replicate)%4==3:print('execution',int(row.geometry)+1,int(row.replicate)+1,'configs',len(estimates),'terminals',terminal_id,'seconds',round(time.time()-GLOBAL_START,1),flush=True)
    c.write_json(c.OUT/'EXECUTION_MANIFEST.json',dict(status='COMPLETE',configurations=len(estimates),terminal_runs=terminal_id,expected_terminal_runs=203008,elapsed_seconds=time.time()-GLOBAL_START,peak_rss_bytes=PEAK,design_sha=head,workers=1,BLAS_threads=4,new_KRAKEN_FIELD_BELLHOP=0,new_audio_or_complex_spectra=0,new_feature_MC=True))
if __name__=='__main__':
    try:
        with threadpool_limits(limits=4):execute()
    except Exception as error:
        c.write_json(c.OUT/'EXECUTION_MANIFEST.json',dict(status='PARTIAL_EXECUTION',design_sha=json.loads((c.OUT/'RUN_CLOCK.json').read_text())['design_sha'] if (c.OUT/'RUN_CLOCK.json').exists() else None,configurations=int(pd.read_csv(c.OUT/'ESTIMATES.csv').shape[0]) if (c.OUT/'ESTIMATES.csv').exists() else 0,terminal_runs=0,reason=str(error),traceback=traceback.format_exc(),elapsed_seconds=time.time()-GLOBAL_START if GLOBAL_START else 0,peak_rss_bytes=PEAK))
        raise
