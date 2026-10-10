"""Independent saved-terminal cold audit and post-estimation truth evaluation."""
import time,json
import numpy as np,pandas as pd,psutil
from threadpoolctl import threadpool_limits
import hla_h2_core as c
import hla_h2_independent as independent
def check_budget():
    clock=json.loads((c.OUT/'RUN_CLOCK.json').read_text())
    if time.time()-clock['epoch']>14400:raise RuntimeError('TOTAL_WALL_BUDGET')
    if psutil.Process().memory_info().rss>8*1024**3:raise RuntimeError('MEMORY_BUDGET')
def observations(row,visible,signed):
    g,rep=int(row.geometry),int(row.replicate);method=row.method
    beta=visible['beta'][g,rep]
    if method=='B':return beta,None
    si=list(c.SIGMAS).index(float(row.sigma))
    if method=='S':return beta,signed['signed'][g,rep,si]
    if method=='U':return beta,visible['unsigned'][g,rep,si]
    return beta,visible['unsigned_biased'][g,rep]
def audit():
    check_budget()
    frozen=json.loads((c.OUT/'DESIGN_FREEZE.json').read_text(encoding='utf-8'))
    input_pass=all(c.sha(c.ROOT/name)==digest for name,digest in frozen['sha256'].items())
    manifest=json.loads((c.OUT/'EXECUTION_MANIFEST.json').read_text());configs=pd.read_csv(c.OUT/'CONFIGURATION_PLAN.csv').set_index('config_id')
    visible=np.load(c.OUT/'OBSERVATIONS.npz');signed=np.load(c.OUT/'SIGNED_CONTROL_OBSERVATIONS.npz')
    count=0;fail=0;maxdiff=0.;branch_sources={};rebuild=[]
    jacmax=0.;jacchecks=0;finite=0;exceptions=0;forward_calls=0
    for chunk in pd.read_csv(c.OUT/'ALL_TERMINALS.csv.gz',chunksize=4096):
        check_budget();chunk['exception']=chunk['exception'].fillna('')
        count+=len(chunk);exceptions+=int((chunk['exception']!='').sum())
        forward_calls+=int(chunk.geometry_calls.fillna(0).sum())
        for config_id,part in chunk.groupby('config_id',sort=False):
            cfg=configs.loc[int(config_id)]
            seen=branch_sources.setdefault(int(config_id),[])
            seen.extend(zip(part.branch_id.astype(int),part.start_id.astype(int)))
            good=part[part['exception']=='']
            if not len(good):continue
            beta,u=observations(cfg,visible,signed)
            h=good[list(c.AXES)].to_numpy(float);bb,qq,_=independent.geometry(h)
            jb=np.square(independent.wrap(bb-beta)/np.radians(.1)).sum(axis=1)
            b=good.b_original.to_numpy(float)
            if cfg.method=='B':
                js=jf=jp=np.zeros(len(good));bp=np.zeros(len(good))
            elif cfg.method=='S':
                js=np.square((qq-u)/float(cfg.sigma)).sum(axis=1);jf=js.copy();jp=js.copy();bp=b.copy()
            else:
                signs=np.array([np.array(c.SIGNS[int(branch)]) for branch in good.branch_id])
                js=np.square((qq+b[:,None]-signs*u)/float(cfg.sigma)).sum(axis=1)
                jf=np.square((abs(qq+b[:,None])-u)/float(cfg.sigma)).sum(axis=1)
                if cfg.method=='UB-matched':
                    bp,raw=independent.profile(qq,np.tile(u,(len(good),1)));jp=raw/float(cfg.sigma)**2
                    # Validate recorded profiled b, not only value of another tied minimum.
                    at_saved=np.square((abs(qq+good.profile_b.to_numpy()[:,None])-u)/float(cfg.sigma)).sum(axis=1)
                else:jp=jf.copy();bp=b.copy()
            original_accept=(jb<=c.T_B)&(jf<=c.T_Q if cfg.method!='B' else True)
            profiled_accept=(jb<=c.T_B)&(jp<=c.T_Q if cfg.method!='B' else True)
            differences=np.column_stack([abs(jb-good.J_beta.to_numpy())/np.maximum(1,abs(good.J_beta.to_numpy())),abs(js-good.J_branch.to_numpy())/np.maximum(1,abs(good.J_branch.to_numpy())),abs(jf-good.J_folded.to_numpy())/np.maximum(1,abs(good.J_folded.to_numpy())),abs(jp-good.J_profile.to_numpy())/np.maximum(1,abs(good.J_profile.to_numpy()))])
            error=differences.max(axis=1)
            if cfg.method=='UB-matched':error=np.maximum(error,abs(at_saved-jp)/np.maximum(1,jp))
            ok=(error<=1e-7)&(original_accept==good.accepted_original.to_numpy())&(profiled_accept==good.accepted_profile.to_numpy())
            fail+=int((~ok).sum());maxdiff=max(maxdiff,float(error.max()));finite+=len(good)
            rebuild.extend(dict(terminal_id=int(i),max_relative_score_difference=float(e),PASS=bool(p)) for i,e,p in zip(good.terminal_id,error,ok))
            if cfg.method=='U' and float(cfg.sigma)==.05 and int(cfg.replicate)==0:
                for _,terminal in good.iterrows():
                    state=terminal[list(c.AXES)].to_numpy(float);sign=c.SIGNS[int(terminal.branch_id)]
                    # Independent centered physical residual differences, no optimization.
                    beta0,q0,_,jba,jqa=c.geometry(state,True)
                    analytic=np.r_[jba[0]/c.SIG_B,jqa[0]/.05]*c.WIDTH[None,:]
                    fd=np.empty_like(analytic)
                    for axis in range(4):
                        delta=np.zeros(4);delta[axis]=c.WIDTH[axis]*1e-6
                        bpv,qpv,_=independent.geometry(state+delta);bmv,qmv,_=independent.geometry(state-delta)
                        fd[:,axis]=np.r_[independent.wrap(bpv[0]-bmv[0])/np.radians(.1),(qpv[0]-qmv[0])/.05]/2e-6
                    je=np.max(np.linalg.norm(analytic-fd,axis=0)/np.maximum(1,np.linalg.norm(fd,axis=0)))
                    jacmax=max(jacmax,float(je));jacchecks+=1
                    if je>1e-5:fail+=1
        if count%32768<len(chunk):print('cold terminals',count,'fail',fail,flush=True)
    estimates=pd.read_csv(c.OUT/'ESTIMATES.csv')
    completed_ids=set(estimates.config_id.astype(int))
    branch_checks=[]
    for cid,cfg in configs.iterrows():
        method=cfg.method
        expected={(i,j) for i in range(64) for j in (0,1)} if method in ('U','U-mismatch','UB-matched') else ({(-1,j) for j in range(17)} if method=='B' else {(0,0),(0,1)})
        actual=branch_sources.get(int(cid),[])
        complete=(set(actual)==expected and len(actual)==len(expected)) if int(cid) in completed_ids else (set(actual).issubset(expected) and len(actual)==len(set(actual)))
        branch_checks.append(dict(config_id=int(cid),method=method,runs=len(actual),expected_runs=len(expected),PASS=complete,completed_configuration=int(cid) in completed_ids))
    pd.DataFrame(branch_checks).to_csv(c.OUT/'BRANCH_EXECUTION_AUDIT.csv',index=False)
    pd.DataFrame(rebuild).to_csv(c.OUT/'TERMINAL_COLD_REBUILD.csv.gz',index=False,compression='gzip')
    # Input truths are first loaded only after all estimator constructions finished.
    hidden=np.load(c.OUT/'GENERATOR_TRUTH_AND_INNOVATIONS.npz')
    observation_fail=0;observation_checks=0
    expected_beta,expected_q,_=independent.geometry(hidden['horizontal_truth'])
    for g in range(8):
        for rep in range(32):
            eb=np.random.default_rng(np.random.SeedSequence([2026101002,g,rep,0])).standard_normal(121)
            eq=np.random.default_rng(np.random.SeedSequence([2026101002,g,rep,1])).standard_normal(6)
            for sigma_index,sigma in enumerate(c.SIGMAS):
                fold=abs(expected_q[g]+sigma*eq)
                observation_fail+=not np.allclose(visible['unsigned'][g,rep,sigma_index],fold,rtol=0,atol=1e-10);observation_checks+=1
                observation_fail+=not np.allclose(signed['signed'][g,rep,sigma_index],expected_q[g]+sigma*eq,rtol=0,atol=1e-10);observation_checks+=1
            observation_fail+=not np.allclose(independent.wrap(visible['beta'][g,rep]-independent.wrap(expected_beta[g]+np.radians(.1)*eb)),0,rtol=0,atol=1e-12);observation_checks+=1
            observation_fail+=not np.allclose(visible['unsigned_biased'][g,rep],abs(expected_q[g]+float(hidden['bias_truth'])+.05*eq),rtol=0,atol=1e-10);observation_checks+=1
    evaluation=[];search=[]
    for _,cfg in configs.reset_index().iterrows():
        cid=int(cfg.config_id);found=estimates[estimates.config_id==cid]
        if not len(found):
            evaluation.append(dict(cfg,status='NOT_EVALUATED'));continue
        est=found.iloc[0];g,rep=int(cfg.geometry),int(cfg.replicate)
        true=hidden['horizontal_truth'][g];beta,u=observations(cfg,visible,signed)
        bb,qq,_=independent.geometry(true);tb=float(np.square(independent.wrap(bb[0]-beta)/c.SIG_B).sum())
        if cfg.method=='B':tq=0.;tp=0.
        elif cfg.method=='S':tq=float(np.square((qq[0]-u)/float(cfg.sigma)).sum());tp=tq
        elif cfg.method=='UB-matched':
            true_b=float(hidden['bias_truth'])
            tq=float(np.square((abs(qq[0]+true_b)-u)/float(cfg.sigma)).sum())
            _,pr=independent.profile(qq,np.atleast_2d(u));tp=float(pr[0]/float(cfg.sigma)**2)
        else:tq=float(np.square((abs(qq[0])-u)/float(cfg.sigma)).sum());tp=tq
        compatible=tb<=c.T_B and (tp<=c.T_Q or cfg.method=='B')
        accepted=est.status=='ACCEPTED'
        record=dict(cfg,status=est.status,true_J_beta=tb,true_J_radial_at_true_bias=tq,true_J_radial_profiled=tp,true_compatible=compatible,search_failure=bool(compatible and not accepted),solver_exceptions=int(est.solver_exceptions),solver_nonconverged=int(est.solver_nonconverged),accepted_terminals=int(est.accepted_terminal_records),b_estimated=est.b,selected_terminal_id=int(est.selected_terminal_id))
        if accepted:
            recovered=est[list(c.AXES)].to_numpy(float);errors=np.abs(recovered-true)
            errors[[0,2]]/=true[[0,2]]
            errors[[1,3]]=np.abs(np.degrees(np.arctan2(np.sin(np.radians(recovered[[1,3]]-true[[1,3]])),np.cos(np.radians(recovered[[1,3]]-true[[1,3]])))))
        else:errors=np.full(4,np.inf)
        for name,error in zip(c.ERROR_NAMES,errors):record[name]=float(error)
        evaluation.append(record)
    results=pd.DataFrame(evaluation);results.to_csv(c.OUT/'CONFIGURATION_RESULTS.csv',index=False,float_format='%.17g')
    results[results.status!='ACCEPTED'].to_csv(c.OUT/'ABSTAIN_AND_NOT_EVALUATED.csv',index=False)
    # All saved terminals available to diagnose rejected best; never primary accuracy.
    terminals=pd.read_csv(c.OUT/'ALL_TERMINALS.csv.gz')
    summaries=[]
    for cid,group in terminals.groupby('config_id',sort=False):
        good=group[group['exception'].isna()]
        if not len(good):continue
        best=good.loc[(good.J_beta+good.J_profile).idxmin()]
        cfg=configs.loc[int(cid)];true=hidden['horizontal_truth'][int(cfg.geometry)]
        error=c.error_vector(best[list(c.AXES)].to_numpy(float),true)
        summaries.append(dict(config_id=int(cid),method=cfg.method,geometry=int(cfg.geometry),replicate=int(cfg.replicate),sigma=float(cfg.sigma),best_diagnostic_terminal=int(best.terminal_id),accepted=bool(best.accepted_profile),best_J_beta=float(best.J_beta),best_J_profile=float(best.J_profile),**{'diagnostic_'+name:float(e) for name,e in zip(c.ERROR_NAMES,error)}))
    pd.DataFrame(summaries).to_csv(c.OUT/'SEARCH_FAILURE_DIAGNOSTICS.csv',index=False)
    selection_fail=0;selection_checks=0
    for _,est in estimates.iterrows():
        group=terminals[(terminals.config_id==est.config_id)&terminals['exception'].isna()]
        eligible=group[group.accepted_profile].copy()
        if len(eligible):
            eligible['total']=eligible.J_beta+eligible.J_profile
            chosen=eligible.sort_values(['total','terminal_id']).iloc[0]
            ok=est.status=='ACCEPTED' and int(est.selected_terminal_id)==int(chosen.terminal_id)
            if ok:ok=bool(np.allclose(est[list(c.AXES)].to_numpy(float),chosen[list(c.AXES)].to_numpy(float),rtol=0,atol=1e-12) and abs(est.b-chosen.profile_b)<=1e-12)
        else:ok=est.status=='ABSTAIN' and int(est.selected_terminal_id)==-1
        selection_checks+=1;selection_fail+=not ok
    full_counts=len(estimates)==2816 and count==203008
    counts_valid=full_counts if manifest['status']=='COMPLETE' else True
    branch_pass=all(r['PASS'] for r in branch_checks)
    result=dict(PASS=fail==0 and observation_fail==0 and selection_fail==0 and input_pass and branch_pass and counts_valid and jacmax<=1e-5,execution_complete=full_counts,observation_rebuild_checks=observation_checks,observation_rebuild_FAIL=observation_fail,selection_rebuild_checks=selection_checks,selection_rebuild_FAIL=selection_fail,terminal_score_checks=finite,terminal_score_FAIL=fail,max_relative_score_difference=maxdiff,frozen_input_hashes_PASS=input_pass,branch_execution_PASS=branch_pass,configuration_rows=len(estimates),terminal_rows=count,terminal_exceptions=exceptions,cold_Jacobian_checks=jacchecks,max_normalized_Jacobian_column_difference=jacmax,geometry_forward_evaluations_recorded=forward_calls,source_depth_used=False,truth_used_for_initialization=False,new_solver_calls=0,independent_research_lead_audit='PENDING',scope='Independent Cartesian + full-sign bias minima + residual/acceptance rebuild for every terminal; 8 registered primary configurations gradient check. Reporting summaries validated separately.')
    c.write_json(c.OUT/'VALIDATION.json',result);print(json.dumps(result,indent=2))
    if not result['PASS']:raise RuntimeError('Independent audit failed')
if __name__=='__main__':
    with threadpool_limits(limits=4):audit()
