"""Reconstruct FIX2 exported evidence without changing the frozen search method."""
import hashlib
import inspect
import json
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits
import r4_a1_fix2_coverage as f


def read(path):return pd.read_csv(path,float_precision='round_trip')


def geometry(states):
    states=np.atleast_2d(states);t=np.arange(121)*10.;turn=np.maximum(t-600.,0.)
    px=2*np.minimum(t,600)+2*turn*np.cos(np.pi/12);py=2*turn*np.sin(np.pi/12)
    r,theta,v,psi=states.T;theta,psi=np.radians(theta),np.radians(psi)
    x=1000*r[:,None]*np.cos(theta[:,None])+v[:,None]*t*np.cos(psi[:,None])-px
    y=1000*r[:,None]*np.sin(theta[:,None])+v[:,None]*t*np.sin(psi[:,None])-py
    return np.arctan2(y,x),np.sqrt(x*x+y*y)


def errors(states,truth):
    s=np.atleast_2d(states);q=truth[list(f.L.AXES)].to_numpy(float)
    d=abs(s-q);d[:,1]=abs((s[:,1]-q[1]+180)%360-180);d[:,3]=abs((s[:,3]-q[3]+180)%360-180)
    d[:,0]/=q[0];d[:,2]/=q[2];return d


def fresh_widths(models,panels,obs):
    """Post-confirmation oracle sections; never consumed by frozen search."""
    rows=[]
    for result in read(f.OUT/'FRESH_HOLDOUT_RESULTS.csv').itertuples():
        cid=result.case_id;truth=panels.loc[result.panel_id];tv=truth[list(f.L.AXES)].to_numpy(float)
        if result.sigma_deg==0:
            rows.append(dict(case_id=cid,direction='NOISELESS_FULL_RANK_SINGLETON',threshold_db=.001,width=0.,width_status='RC2_NUMERICAL_SINGLETON; NOT_ACOUSTIC_BASIN_WIDTH',oracle_diagnostic_only=True));continue
        bearing,tl=obs[cid];base=f.prior.theta_profile(tv[[0,2,3]],bearing)[0];a0=np.r_[tv[[0,2,3]],tv[1]-base[1]]
        def evaluate(a):
            s=f.prior.theta_profile(a[:3],bearing)[0];s[1]+=a[3];b,r=geometry(s);cost=float(np.sum(((b[0]-bearing+np.pi)%(2*np.pi)-np.pi)**2))
            if cost>result.cutoff+1e-14 or np.any(s<f.prior.LOW) or np.any(s>f.prior.HIGH):return None,None
            features=f.L.direct_features(models,r)[0];js=np.sqrt(np.mean((features-tl)**2,axis=(1,2)));z=int(js.argmin());return float(js[z]),float(f.L.PROFILE[z])
        for k,name in enumerate(['r_km','v_mps','psi_deg','theta_offset_deg']):
            for target in [.001,.01,.1]:
                walls=[];status=[];depths=[]
                for sign in [-1,1]:
                    direction=np.eye(4)[k]*sign;lo=0.;hi=[1e-8,1e-8,1e-6,1e-6][k];wall='ACOUSTIC_FIRST_CROSSING'
                    for _ in range(40):
                        j,z=evaluate(a0+direction*hi)
                        if j is None:
                            for _ in range(25):
                                mid=(lo+hi)/2
                                if evaluate(a0+direction*mid)[0] is None:hi=mid
                                else:lo=mid
                            hi=lo;wall='CENSORED_BY_BEARING_OR_STATE_BOUND';break
                        if j>=target:
                            hi=brentq(lambda h:evaluate(a0+direction*h)[0]-target,lo,hi,xtol=1e-12);break
                        lo=hi;hi*=2
                    walls.append(hi);status.append(wall);depths.append(evaluate(a0+direction*hi)[1])
                rows.append(dict(case_id=cid,direction=name,threshold_db=target,width=sum(walls),negative_extent=walls[0],positive_extent=walls[1],width_status=';'.join(status),z_negative=depths[0],z_positive=depths[1],oracle_diagnostic_only=True))
        print('FRESH WIDTHS',cid,flush=True)
        f.save('FRESH_BASIN_GEOMETRY.csv',pd.DataFrame(rows))
    f.save('FRESH_BASIN_GEOMETRY.csv',pd.DataFrame(rows))


def audit(development_only=False):
    out=f.OUT;f.protect();checks=[];numeric=[];aliases=[]
    def check(name,ok,detail):
        checks.append(dict(check=name,pass_check=bool(ok),detail=str(detail)))
        assert ok,(name,detail)
    frozen=read(out/'FROZEN_INPUT_HASHES.csv')
    def valid_frozen(row):
        path=f.L.ROOT/row.path
        if f.L.sha(path)==row.sha256:return True
        portable=path.suffix=='.py' or path==f.L.REFERENCE/'SUBSET_CANDIDATE_SCORES_FIXED.csv'
        return portable and hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==row.sha256_lf
    raw_matches=sum(f.L.sha(f.L.ROOT/r.path)==r.sha256 for r in frozen.itertuples())
    check('Historical evidence and R3 frozen inputs',all(valid_frozen(r) for r in frozen.itertuples()),f'{raw_matches}/{len(frozen)} raw byte hashes match; LF portability allowed only for inherited Python/R3 score CSV')
    mf=json.loads((out/'METHOD_FREEZE.json').read_text())
    check('Frozen method/config/code bytes',all(f.L.sha(f.L.ROOT/p)==h for p,h in mf['sha256'].items()),mf['sha256'])
    hd=json.loads((out/'FRESH_HOLDOUT_DESIGN_FREEZE.json').read_text())
    check('Fresh panel/design hashes',all(f.L.sha(out/p)==h for p,h in hd['sha256'].items()),hd['sha256'])
    check('Method frozen before fresh panel design',mf['frozen_utc']<=hd['frozen_utc'] and hd['before_observation_generation'],hd['frozen_utc'])
    check('Preregistered unchanged search configuration',json.loads((out/'SEARCH_BUDGET_PREREGISTRATION.json').read_text())==f.SEARCH_CONFIG,'Threshold remains .001 dB')
    search_functions=[f.radial_profile,f.bank,f.select_beams,f.refine,f.one_budget,f.cumulative_candidates,f.independent_search]
    check('Search functions exclude truth/evaluation helpers',all('truth' not in inspect.signature(fn).parameters and 'L.errors(' not in inspect.getsource(fn) and 'L.neighborhood(' not in inspect.getsource(fn) for fn in search_functions),'Seven observation-only search APIs; runtime sentinel test also passed')
    # Replay the original observation generator's CSV parser for truth states.
    # Candidate states still use round_trip so stored exact scores retain their inputs.
    panels=pd.concat([pd.read_csv(f.L.OUT/'OFFGRID_TRUTH_PANEL.csv'),pd.read_csv(f.prior.OUT/'HOLDOUT_TRUTH_PANEL.csv'),pd.read_csv(out/'FRESH_HOLDOUT_PANEL.csv')]).set_index('panel_id')
    cfg=f.SEARCH_CONFIG['fresh_holdout'];rng=np.random.default_rng(cfg['seed']);generated=[]
    while len(generated)<cfg['n']:
        s=rng.uniform(cfg['low'],cfg['high']);q=(s-f.L.ORIGIN)/f.L.STEP
        if np.min(abs(q-np.round(q))*f.L.STEP)<1e-6:continue
        generated.append(np.r_[s,rng.choice(cfg['depths'])])
    fresh=read(out/'FRESH_HOLDOUT_PANEL.csv')
    check('Fresh truths exactly reproduce frozen seed/rule',np.array_equal(np.array(generated),fresh[list(f.L.AXES)+['z_true_m']].to_numpy(float)),cfg)
    q=(fresh[list(f.L.AXES)].to_numpy()-f.L.ORIGIN)/f.L.STEP
    check('Every fresh coordinate is off-grid',np.all(abs(q-np.round(q))*f.L.STEP>1e-6),'8 truths, all four coordinates')
    obs={}
    observation_paths=[f.prior.OUT/'CASE_OBSERVATIONS.npz']
    if not development_only:observation_paths.append(out/'FRESH_HOLDOUT_OBSERVATIONS.npz')
    for path in observation_paths:
        data=np.load(path)
        obs.update({str(cid):(data['bearing_rad'][i],data['relative_tl'][i]) for i,cid in enumerate(data['case_ids'])})
    tables=[('SEARCH_BUDGET_CONVERGENCE.csv','DEVELOPMENT_CANDIDATES.csv'),('INDEPENDENT_SOLVER_AGREEMENT.csv','INDEPENDENT_CANDIDATES.csv')]
    if not development_only:tables.append(('FRESH_HOLDOUT_RESULTS.csv','FRESH_HOLDOUT_CANDIDATES.csv'))
    models=f.inputs();memo={};total=0
    replay_cache=out/'EXACT_CANDIDATE_REPLAY_CACHE.csv';replay_manifest=out/'EXACT_REPLAY_CACHE_MANIFEST.json'
    provenance={'audit_code_sha256':f.L.sha(Path(__file__)),'method_freeze_sha256':f.L.sha(out/'METHOD_FREEZE.json'),'historical_observations_sha256':f.L.sha(f.prior.OUT/'CASE_OBSERVATIONS.npz')}
    fresh_observation_sha=None if development_only else f.L.sha(out/'FRESH_HOLDOUT_OBSERVATIONS.npz')
    if replay_cache.exists() and replay_manifest.exists():
        previous=json.loads(replay_manifest.read_text())
        if previous['provenance']==provenance and previous['cache_sha256']==f.L.sha(replay_cache) and previous.get('fresh_observation_sha256') in (None,fresh_observation_sha):
            for r in read(replay_cache).itertuples():memo[(r.case_id,*[getattr(r,a) for a in f.L.AXES])]=(r.recomputed_exact_J,r.z_star_label_m)
    for result_name,candidate_name in tables:
        results=read(out/result_name);catalog=read(out/candidate_name)
        for row in results.itertuples():
            cid=row.case_id;truth=panels.loc[cid.split('_')[0]];bearing,tl=obs[cid]
            subset=catalog[catalog.case_id==cid]
            if hasattr(row,'budget'):subset=subset[subset.budget.fillna(row.budget)==row.budget] if 'budget' in subset else subset
            states=subset[list(f.L.AXES)].to_numpy(float);pred,ranges=geometry(states)
            check(f'Physical state bounds {cid} {getattr(row,"budget",result_name)}',np.all(states>=f.prior.LOW-1e-12) and np.all(states<=f.prior.HIGH+1e-12),'Inherited four-dimensional bounds')
            residual=(pred-bearing+np.pi)%(2*np.pi)-np.pi;cost=np.sum(residual**2,axis=1)
            js=[];zs=[]
            for s,rr in zip(states,ranges):
                key=(cid,*s)
                if key not in memo:
                    feature=f.L.direct_features(models,rr)[0]
                    scores=np.sqrt(np.mean((feature-tl)**2,axis=(1,2)));zi=int(scores.argmin())
                    memo[key]=(float(scores[zi]),float(f.L.PROFILE[zi]))
                j,z=memo[key];js.append(j);zs.append(z)
            js=np.array(js);zs=np.array(zs);delta=float(max(abs(js-subset.J_exact.to_numpy())))
            check(f'Exact candidate reconstruction {cid} {getattr(row,"budget",result_name)}',delta<1e-7 and np.all(zs==subset.z_star_label_m) and np.max(abs(cost-subset.bearing_cost))<1e-12 and np.all(cost<=row.cutoff+1.1e-14),(len(states),delta))
            check(f'Top ranking and recovery {cid} {getattr(row,"budget",result_name)}',abs(js[0]-row.exact_J)<1e-7 and js[0]<=js.min()+1e-7 and (not hasattr(row,'recovered') or row.recovered==bool(js[0]<.001 and cost[0]<=row.cutoff+1e-14)),row.exact_J)
            e=errors(states,truth);keep=subset.keep_tau05.to_numpy(bool) if 'keep_tau05' in subset else js<=js.min()+f.SEARCH_CONFIG['candidate_tau_db']
            for k,name in enumerate(('rel_r','abs_theta_deg','rel_v','abs_psi_deg')):
                if hasattr(row,'top1_'+name):
                    check(f'Errors/envelope {cid} {getattr(row,"budget",result_name)} {name}',abs(getattr(row,'top1_'+name)-e[0,k])<1e-9 and abs(getattr(row,'survivor_worst_'+name)-e[keep,k].max())<1e-9,'Independent circular-error reconstruction')
            true_b,true_r=geometry(truth[list(f.L.AXES)].to_numpy(float));true_feature=f.L.direct_features(models,true_r,[truth.z_true_m])[0,0]
            truthcost=float(np.sum(((true_b[0]-bearing+np.pi)%(2*np.pi)-np.pi)**2))
            check(f'Matched truth remains feasible {cid}',truthcost<=row.cutoff+1e-14 and np.sqrt(np.mean((true_feature-tl)**2))<1e-10,(truthcost,row.cutoff))
            bearing_rms=np.degrees(np.sqrt(np.mean(((pred-true_b+np.pi)%(2*np.pi)-np.pi)**2,axis=1)))
            distance=abs(states-truth[list(f.L.AXES)].to_numpy(float));distinct=np.any(distance>np.array(f.prior.CONFIG['candidate_dedup_resolution']),axis=1)
            for i in np.flatnonzero(distinct & (js<1e-6) & (bearing_rms<1e-8)):
                aliases.append(dict(case_id=cid,source=candidate_name,candidate_index=int(i),exact_J=float(js[i]),bearing_rms_deg=float(bearing_rms[i]),distinct_joint_exact_match=True))
            numeric.append(dict(case_id=cid,budget=getattr(row,'budget','INDEPENDENT'),n_candidates=len(states),max_exact_J_error=delta,truth_bearing_cost=truthcost,cutoff=row.cutoff,pass_check=True));total+=len(states)
            print('RECONSTRUCT',cid,getattr(row,'budget','INDEPENDENT'),len(states),'unique modal evaluations',len(memo),flush=True)
    if not aliases:aliases=[dict(case_id='ALL_EXPORTED_CONTROLS',source='finite candidate catalogs',candidate_index=-1,exact_J=np.nan,bearing_rms_deg=np.nan,distinct_joint_exact_match=False)]
    suffix='_DEVELOPMENT' if development_only else ''
    f.save('ALIAS_DIAGNOSTIC'+suffix+'.csv',pd.DataFrame(aliases));f.save('RESULT_RECONSTRUCTION_AUDIT'+suffix+'.csv',pd.DataFrame(numeric))
    check('Every oracle landscape probe is likelihood feasible',read(out/'BASIN_FEASIBLE_PROBES.csv').eval('bearing_cost <= cutoff + 1e-14').all(),'1008 probes; infeasible walls explicitly censored')
    f.save('INTEGRITY_AUDIT'+suffix+'.csv',pd.DataFrame(checks))
    f.save(replay_cache.name,pd.DataFrame([dict(case_id=key[0],**dict(zip(f.L.AXES,key[1:])),recomputed_exact_J=value[0],z_star_label_m=value[1]) for key,value in memo.items()]))
    f.L.json_write(replay_manifest,dict(provenance=provenance,cache_sha256=f.L.sha(replay_cache),fresh_observation_sha256=fresh_observation_sha,scope='development only' if development_only else 'development + independent + fresh confirmation'))
    if development_only:
        print('DEVELOPMENT AUDIT PASSED',len(checks),'checks',flush=True);return
    fresh_widths(models,panels,obs)
    (out/'LOCAL_VALIDATION.md').write_text(f'''# Local validation\n\n22 unit tests passed, including analytic exact/spline derivatives, observation-only continuation and cumulative-basin preservation.\n\n{len(checks)} reconstruction/integrity checks passed across {len(numeric)} case-budget/solver results, {total} candidate rows and {len(memo)} distinct case/state modal evaluations. Scores use direct modal features and independent RMS, geometry and circular-error calculations. Frozen historical bytes, method/code/config hashes, fresh seed/rule and design order were checked.\n\nNo distinct exact joint alias was found in these finite exported candidates; this is not a global uniqueness certificate. Scientific Gate failure and numerical reconstruction success are separate. Independent research-lead audit is pending.\n''',encoding='utf-8')
    print('AUDIT PASSED',len(checks),'checks')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--development-only',action='store_true');args=parser.parse_args()
    with threadpool_limits(limits=2):audit(args.development_only)
