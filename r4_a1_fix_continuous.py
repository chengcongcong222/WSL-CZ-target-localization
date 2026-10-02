"""R4 repair: observation-generated continuous candidates; frozen legacy untouched."""
from __future__ import annotations
import argparse
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares, differential_evolution
from scipy.stats import qmc
from threadpoolctl import threadpool_limits
import r4_a1_offgrid_bearing as legacy

OUT = legacy.ROOT/'results/R4_A1_FIX_CONTINUOUS_SEARCH'
LOW = np.array([45.,-5.,1.,-15.])
HIGH = np.array([60.,5.,3.,15.])
WIDTH = HIGH-LOW
CONFIG = {
    'implementation_version':2,
    'initialization_correction':'Initialize nondegenerate theta-offset diversity within the continuous likelihood region; V1 evidence retained in INITIAL_IMPLEMENTATION.',
    'baseline_commit':'4633ff091874bc27d54017a9aff3b87cd4eed97d',
    'method':'bounded bearing multistart + observation-profiled theta + two global acoustic DE islands + multibasin local refinement + exact modal verification',
    'RC2_LS_starts':32,'RC2_LS_max_nfev':300,'RC2_sobol_power':12,
    'RC2_cutoff_coefficient':13.3,'noiseless_cost_tolerance':1e-20,
    'forward_acoustic_tolerance_db':1e-10,
    'parameter_bounds':{'low':LOW.tolist(),'high':HIGH.tolist()},
    'pressure_spline_step_m':1.,'pressure_spline_range_m':[39000.,66000.],
    'pressure_spline_validation_seed':20261003,'pressure_spline_validation_points':256,
    'pressure_spline_TL_tolerance_db':.001,
    'acoustic_DE_population':128,'acoustic_DE_generations':160,'acoustic_DE_islands':2,
    'acoustic_DE_seeds':[771001,771002],'DE_mutation':[.5,1.],'DE_recombination':.7,
    'local_basins_per_island':8,'local_max_nfev':150,'exact_polish_basins':4,'exact_polish_max_nfev':80,
    'candidate_tau_db':.5,'acoustic_basin_recovery_db':.001,
    'candidate_dedup_resolution':[.005,.001,.005,.05],
    'development_nominal_seeds':[410001],
    'holdout_rule':{'seed':2026100203,'panel_prefix':'Q','n':6,'r_km':[47.,55.5],'theta_deg':[-1.2,1.2],
        'v_mps':[1.5,2.5],'psi_deg':[-9.,9.],'z_m':[180.,200.,220.],
        'reject_distance_to_any_coarse_axis_less_than':1e-6,'nominal_seeds':[420001,420002]},
    'gate_scope':'noiseless + nominal repair only; no engineering tolerance claim; R4 credit remains zero',
    'all_hyperparameters':'Frozen before development search and holdout observation generation; no retuning to panel success.'}


def save(name,frame):
    frame.to_csv(OUT/name,index=False,float_format='%.17g',lineterminator='\n')


def md_table(frame):
    def value(v):return f'{v:.7g}' if isinstance(v,float) else str(v)
    header='| '+' | '.join(frame.columns)+' |\n| '+' | '.join(['---']*len(frame.columns))+' |'
    return header+'\n'+'\n'.join('| '+' | '.join(value(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None))


def freeze():
    OUT.mkdir(exist_ok=True)
    p=OUT/'R4_A1_FIX_CONFIG.json'
    if p.exists():
        if json.loads(p.read_text())!=CONFIG: raise ValueError('Method configuration changed after freeze')
    else:
        legacy.json_write(p,CONFIG)
    method=OUT/'METHOD_FREEZE.md'
    if not method.exists():
        method.write_text('''# Continuous-search method freeze

The bearing geometry is linear in Cartesian initial position and velocity for noiseless observations. Its full-rank solution is an observation-derived start, verified through bounded nonlinear circular residual least squares with 32 independent Sobol starts. For noisy data, local minima alone do not represent the likelihood region: 4096 global Sobol (r,v,psi) points are continued by optimizing theta using bearing observations, retaining those within the unchanged continuous 13.3 sigma-squared likelihood cutoff.

No grid-bias allowance is added. Search states are continuous, so the representation problem is removed rather than compensated by a fitted threshold. RC2 extent and individual optimizer convergence are exported; candidate sampling is finite and is not a certified exhaustive confidence region.

Implementation version 2 corrects a discovered initialization defect: V1 initialized all theta offsets to zero, so differential evolution could not vary that coordinate globally. V2 initializes theta offsets with independent random diversity, accepting only continuously feasible points. Search budgets and scientific thresholds remain unchanged. V1 evidence and its frozen design are retained in INITIAL_IMPLEMENTATION. Since partial V1 holdout results were seen, V2 uses a new deterministic holdout seed 2026100203, frozen before its observations. This is a programming correction, not tuning to pass the nine regression cases.

Acoustic search uses two independent global differential-evolution islands seeded exclusively by the continuous bearing region, followed by local refinement of eight diverse basins per island. It searches r,v,psi and a theta offset from observation-profiled theta; feasibility uses exact bearing cost. A conservative offset bound follows from the minimum bearing theta Jacobian over the full geometry bounds. Shared z is profiled over the inherited 21 nuisance labels. Final candidates are verified with exact modal evaluation, and four distinct best basins are polished with the exact model.

For affordable global evaluation, interpolate complex pressure after removing the common modal carrier, never coarse-grid TL. All modal wavenumbers remain in the forward model. The 1 m spline resolves the maximum demodulated wavenumber below 0.1 rad/m; it is checked at 256 fixed random ranges across all nuisance depths/frequencies, with a predeclared 0.001 dB maximum TL discrepancy. Exact forward scoring determines all final outputs. Range bounds 39--66 km follow the triangle inequality over the full state and platform trajectory bounds.

Benchmark: one 283 Hz cache takes about 1.54 s; a 115-state/121-time spline batch about 0.0065 s. Its fixed random-range maximum TL discrepancy is about 9.45e-6 dB. Exact per-trajectory coarse scoring previously took about 0.015 s/node. These forward-cost/geometry measurements determine the method and budget; they do not use panel recovery to tune parameters.

Finite-budget global search may fail. Failure is reported as search failure, not converted into continuous non-identifiability. Mode counts mean separated sampled local basins at the configured coordinate resolution, not a proof of all physical modes. Best-candidate error is evaluation-only, coordinatewise minima, not a selectable joint oracle state. Noiseless geometry has an independently checked rank-four control; near-zero bearing candidates may collapse to one hypothesis without acoustic grid initialization.

Hyperparameters and deterministic holdout generation rule are in R4_A1_FIX_CONFIG.json. Freeze its hash together with the generated holdout panel before creating holdout observations. The nine legacy truths are a development/regression panel. No full multi-sigma Monte Carlo is run here.
''',encoding='utf-8')
    hp=OUT/'HOLDOUT_TRUTH_PANEL.csv'
    if not hp.exists():
        rule=CONFIG['holdout_rule']; rng=np.random.default_rng(rule['seed']); rows=[]
        while len(rows)<rule['n']:
            s=np.array([rng.uniform(*rule[a]) for a in legacy.AXES])
            q=(s-legacy.ORIGIN)/legacy.STEP
            if np.min(np.abs(q-np.round(q))*legacy.STEP)<1e-6: continue
            rows.append(dict(panel_id=f"{rule['panel_prefix']}{len(rows)+1:02d}",**dict(zip(legacy.AXES,s)),z_true_m=float(rng.choice(rule['z_m']))))
        save(hp.name,pd.DataFrame(rows))
    manifest=OUT/'HOLDOUT_DESIGN_FREEZE.json'
    digests={p.name:legacy.sha(p) for p in (p,method,hp)}
    if manifest.exists():
        assert json.loads(manifest.read_text())['sha256']==digests
    else:
        from datetime import datetime,timezone
        legacy.json_write(manifest,{'frozen_utc':datetime.now(timezone.utc).isoformat(),'sha256':digests,'rule':CONFIG['holdout_rule'],'before_holdout_observations':True})


def angular_residual(states,observed):
    return legacy.wrap_rad(legacy.geometry(states)[0]-observed)


def theta_profile(rvpsi,observed):
    a=np.atleast_2d(rvpsi)
    s=np.column_stack((a[:,0],np.full(len(a),np.degrees(observed[0])),a[:,1],a[:,2]))
    xp,yp=legacy.platform_xy(); t=legacy.TIMES
    for _ in range(4):
        th,psi=np.radians(s[:,1]),np.radians(s[:,3]); r=s[:,0]*1000
        dx=r[:,None]*np.cos(th)[:,None]+s[:,2,None]*t*np.cos(psi)[:,None]-xp
        dy=r[:,None]*np.sin(th)[:,None]+s[:,2,None]*t*np.sin(psi)[:,None]-yp
        jac=r[:,None]*(dx*np.cos(th)[:,None]+dy*np.sin(th)[:,None])/(dx*dx+dy*dy)
        residual=legacy.wrap_rad(np.arctan2(dy,dx)-observed)
        th-=np.sum(residual*jac,axis=1)/np.sum(jac*jac,axis=1)
        s[:,1]=np.clip(np.degrees(th),LOW[1],HIGH[1])
    return s


def unique_indices(states):
    kept=[]; resolution=np.array(CONFIG['candidate_dedup_resolution'])
    for i,s in enumerate(states):
        if not kept or all(np.max(np.abs(s-states[j])/resolution)>1 for j in kept): kept.append(i)
    return np.array(kept,dtype=int)


def rc2(observed,sigma):
    b=observed; xp,yp=legacy.platform_xy()
    A=np.column_stack((-1000*np.sin(b),1000*np.cos(b),-legacy.TIMES*np.sin(b),legacy.TIMES*np.cos(b)))
    rhs=-xp*np.sin(b)+yp*np.cos(b)
    cart=np.linalg.lstsq(A,rhs,rcond=None)[0]
    linear=np.array([np.hypot(*cart[:2]),np.degrees(np.arctan2(cart[1],cart[0])),np.hypot(*cart[2:]),np.degrees(np.arctan2(cart[3],cart[2]))])
    starts=qmc.Sobol(4,scramble=True,seed=770001).random_base2(5)
    starts[0]=np.clip((linear-LOW)/WIDTH,1e-9,1-1e-9)
    logs=[]; solutions=[]
    for i,x in enumerate(starts):
        fit=least_squares(lambda u:angular_residual(LOW+WIDTH*u,observed)[0],x,bounds=(np.zeros(4),np.ones(4)),
            max_nfev=CONFIG['RC2_LS_max_nfev'],xtol=1e-12,ftol=1e-12,gtol=1e-12)
        s=LOW+WIDTH*fit.x; cost=float(np.square(fit.fun).sum())
        logs.append(dict(start_id=i,success=fit.success,status=fit.status,nfev=fit.nfev,bearing_cost=cost,optimality=fit.optimality,**dict(zip(legacy.AXES,s))))
        solutions.append(s)
    solutions=np.array(solutions); costs=np.square(angular_residual(solutions,observed)).sum(axis=1); cmin=float(costs.min())
    cutoff=cmin+(CONFIG['noiseless_cost_tolerance'] if sigma==0 else 13.3*np.radians(sigma)**2)
    if sigma==0:
        eligible=solutions[costs<=cutoff]; cloud=eligible[unique_indices(eligible)]
    else:
        samples=qmc.Sobol(3,scramble=True,seed=770002).random_base2(CONFIG['RC2_sobol_power'])
        rv=LOW[[0,2,3]]+samples*WIDTH[[0,2,3]]
        candidate=theta_profile(rv,observed)
        keep=np.square(angular_residual(candidate,observed)).sum(axis=1)<=cutoff
        cloud=np.concatenate((solutions[costs<=cutoff],candidate[keep]))
        cloud=cloud[unique_indices(cloud)]
    return cloud,cmin,cutoff,pd.DataFrame(logs),np.linalg.svd(A,compute_uv=False)


class ContinuousAcoustic:
    def __init__(self,models):
        self.models=models; self.splines={}; records=[]
        x=np.arange(39000.,66001.,CONFIG['pressure_spline_step_m'])
        points=np.random.default_rng(CONFIG['pressure_spline_validation_seed']).uniform(39000,66000,256)
        for freq,m in models.items():
            iz,_=legacy.effective_depths(m,legacy.PROFILE); receiver,_=legacy.effective_depths(m,[200.])
            w=m['phi'][iz]*m['phi'][receiver[0]]; k=m['k'].real.astype(float); alpha=-m['k'].imag.astype(float); carrier=(k.min()+k.max())/2
            def pressure(rr):
                rr=np.atleast_1d(rr)[None,:]
                factor=np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*(k[:,None]-carrier)*rr-alpha[:,None]*rr-1j*np.pi/4)
                return (w@factor).T
            p=np.empty((len(x),21),complex)
            for start in range(0,len(x),2048):p[start:start+2048]=pressure(x[start:start+2048])
            spline=CubicSpline(x,p,axis=0,extrapolate=False); self.splines[freq]=spline
            delta=np.abs(20*np.log10(np.abs(spline(points)))-20*np.log10(np.abs(pressure(points))))
            records.append(dict(frequency_hz=freq,max_TL_error_db=float(delta.max()),P95_TL_error_db=float(np.quantile(delta,.95)),tolerance_db=CONFIG['pressure_spline_TL_tolerance_db'],pass_validation=delta.max()<CONFIG['pressure_spline_TL_tolerance_db']))
        save('CONTINUOUS_PRESSURE_INTERPOLATION_CONTROL.csv',pd.DataFrame(records))
        if not all(x['pass_validation'] for x in records):raise ValueError('Continuous pressure interpolation control failed')

    def features(self,states):
        ranges=legacy.geometry(states)[1]; n=len(ranges)
        result=np.empty((n,21,3,121))
        for fi,freq in enumerate(legacy.FREQS):
            p=self.splines[freq](ranges)
            result[:,:,fi,:]=(20*np.log10(np.maximum(np.abs(p),1e-30))).transpose(0,2,1)
        return legacy.demean_windows(result)

    def score(self,states,tl):return legacy.score_features(self.features(states),tl)

    def exact(self,states,tl):return legacy.score_features(legacy.direct_features(self.models,legacy.geometry(states)[1]),tl)


def acoustic_population(cloud,bearing,sigma,cutoff,rng):
    dims=[0,2,3]
    jac_lower=45000*39000/66000**2
    delta_bound=sigma*np.sqrt(13.3/len(bearing))/jac_lower
    ids=rng.choice(len(cloud),size=CONFIG['acoustic_DE_population'],replace=len(cloud)<CONFIG['acoustic_DE_population'])
    base=theta_profile(cloud[ids][:,dims],bearing)
    offsets=np.zeros(len(base))
    pending=np.arange(len(base))
    for _ in range(32):
        if not len(pending):break
        delta=rng.uniform(-1,1,len(pending)); trial=base[pending].copy(); trial[:,1]+=delta*delta_bound
        costs=np.square(angular_residual(trial,bearing)).sum(axis=1)
        valid=(costs<=cutoff)&(trial[:,1]>=LOW[1])&(trial[:,1]<=HIGH[1])
        offsets[pending[valid]]=delta[valid]; pending=pending[~valid]
    population=np.column_stack(((base[:,dims]-LOW[dims])/WIDTH[dims],offsets))
    if np.ptp(offsets)<1e-8:raise RuntimeError('Degenerate theta-offset DE initialization')
    return population


def acoustic_search(model,bearing,tl,sigma,cloud,cutoff):
    logs=[]; candidates=[]
    if sigma==0:
        js,zs=model.exact(cloud,tl)
        return cloud,js,zs,pd.DataFrame([dict(stage='NOISELESS_RC2_HYPOTHESES',island=0,success=True,nfev=len(cloud),Jmin=float(js.min()))])
    # Global lower bound on d(bearing)/d(theta) from triangle inequalities.
    jac_lower=45000*39000/66000**2
    delta_bound=sigma*np.sqrt(13.3/len(bearing))/jac_lower
    dims=[0,2,3]
    def unpack(u):
        u=np.atleast_2d(u)
        rv=LOW[dims]+u[:,:3]*WIDTH[dims]
        s=theta_profile(rv,bearing); s[:,1]+=u[:,3]*delta_bound
        return s
    def objective(u):
        s=unpack(np.atleast_2d(u).T)
        costs=np.square(angular_residual(s,bearing)).sum(axis=1)
        feasible=(costs<=cutoff)&(s[:,1]>=LOW[1])&(s[:,1]<=HIGH[1])
        value=500+np.maximum((costs-cutoff)/np.radians(sigma)**2,0)
        if feasible.any():value[feasible]=model.score(s[feasible],tl)[0]
        return value
    for island,seed in enumerate(CONFIG['acoustic_DE_seeds']):
        rng=np.random.default_rng(seed)
        init=acoustic_population(cloud,bearing,sigma,cutoff,rng)
        fit=differential_evolution(objective,[(0,1)]*3+[(-1,1)],init=init,maxiter=CONFIG['acoustic_DE_generations'],
            mutation=tuple(CONFIG['DE_mutation']),recombination=CONFIG['DE_recombination'],seed=seed,polish=False,
            vectorized=True,updating='deferred',tol=0,atol=0)
        ranked=np.argsort(fit.population_energies); pop=unpack(fit.population[ranked]); selected=unique_indices(pop)[:8]
        candidates.extend(pop[selected])
        logs.append(dict(stage='GLOBAL_DE',island=island,success=fit.success,nfev=fit.nfev,nit=fit.nit,Jmin=float(fit.fun),message=str(fit.message)))
        for start in pop[selected]:
            _,zz=model.score(start,tl); depth=int(np.flatnonzero(legacy.PROFILE==zz[0])[0])
            def residual(u):
                s=LOW+WIDTH*u
                acoustic=(model.features(s)[0,depth]-tl).ravel()/np.sqrt(tl.size)
                excess=max(float(np.square(angular_residual(s,bearing)).sum())-cutoff,0)
                return np.r_[acoustic,10*np.sqrt(excess)/np.radians(sigma)]
            local=least_squares(residual,np.clip((start-LOW)/WIDTH,1e-10,1-1e-10),bounds=(np.zeros(4),np.ones(4)),
                max_nfev=CONFIG['local_max_nfev'],xtol=1e-10,ftol=1e-10,gtol=1e-10,diff_step=1e-6)
            s=LOW+WIDTH*local.x; candidates.append(s)
            j,z=model.score(s,tl)
            logs.append(dict(stage='LOCAL_SPLINE',island=island,success=local.success,nfev=local.nfev,Jmin=float(j[0]),bearing_cost=float(np.square(angular_residual(s,bearing)).sum())))
    candidates=np.array(candidates)
    # Preserve unpolished global hypotheses as well as all local basins.
    candidates=np.concatenate((candidates,cloud[np.argsort(model.score(cloud,tl)[0])[:8]]))
    j,z=model.score(candidates,tl); order=np.argsort(j); candidates=candidates[order]; candidates=candidates[unique_indices(candidates)]
    j,z=model.exact(candidates,tl)
    for i in np.argsort(j)[:CONFIG['exact_polish_basins']]:
        start=candidates[i]; depth=z[i]
        def exact_residual(u):
            s=LOW+WIDTH*u
            f=legacy.direct_features(model.models,legacy.geometry(s)[1],[depth])[0,0]
            excess=max(float(np.square(angular_residual(s,bearing)).sum())-cutoff,0)
            return np.r_[(f-tl).ravel()/np.sqrt(tl.size),10*np.sqrt(excess)/np.radians(sigma)]
        fit=least_squares(exact_residual,np.clip((start-LOW)/WIDTH,1e-10,1-1e-10),bounds=(np.zeros(4),np.ones(4)),
            max_nfev=CONFIG['exact_polish_max_nfev'],xtol=1e-11,ftol=1e-11,gtol=1e-11,diff_step=1e-7)
        s=LOW+WIDTH*fit.x; candidates=np.vstack((candidates,s))
        js,zs=model.exact(s,tl)
        logs.append(dict(stage='LOCAL_EXACT',island=-1,success=fit.success,nfev=fit.nfev,Jmin=float(js[0]),bearing_cost=float(np.square(angular_residual(s,bearing)).sum())))
    costs=np.square(angular_residual(candidates,bearing)).sum(axis=1)
    keep=costs<=cutoff+1e-14
    candidates=candidates[keep]
    if not len(candidates):raise RuntimeError('No acoustic candidate satisfies continuous RC2')
    j,z=model.exact(candidates,tl); order=np.argsort(j); candidates,j,z=candidates[order],j[order],z[order]
    keep=unique_indices(candidates)
    return candidates[keep],j[keep],z[keep],pd.DataFrame(logs)


def controls(models):
    panel=pd.read_csv(legacy.OUT/'OFFGRID_TRUTH_PANEL.csv'); obs=np.load(legacy.OUT/'PILOT_OBSERVATIONS.npz'); rows=[]
    for pi,truth in panel.iterrows():
        perfect=legacy.generate_observation(truth,0,0,pi,models)
        s=truth[list(legacy.AXES)].to_numpy(float)
        j,z=legacy.score_features(legacy.direct_features(models,legacy.geometry(s)[1]),obs['relative_tl'][pi])
        cost=float(np.square(angular_residual(s,perfect.bearing_rad)).sum())
        rows.append(dict(panel_id=truth.panel_id,bearing_cost_noiseless=cost,J_TRIPLE=float(j[0]),z_star_label_m=float(z[0]),bearing_tolerance=1e-20,acoustic_tolerance_db=1e-10,pass_self_match=cost<1e-20 and j[0]<1e-10))
    save('CONTINUOUS_FORWARD_SELF_MATCH.csv',pd.DataFrame(rows))
    if not all(x['pass_self_match'] for x in rows):raise RuntimeError('A1_FIX_BLOCKED_BY_CONTINUOUS_FORWARD_INCONSISTENCY')
    # Reproduce legacy R3 IDs and scores using read-only models / reference CSV.
    estimator=legacy.FrozenEstimator(legacy.grid())
    control=pd.Series(dict(r_km=50.,theta_deg=0.,v_mps=2.,psi_deg=5.,z_true_m=200.))
    b,rg=legacy.geometry(control[list(legacy.AXES)].to_numpy(float))
    b=b[0]+np.random.default_rng(20260912).normal(0,np.radians(.1),121)
    ids,_=estimator.select_cloud(estimator.bearing_costs(b)[:,0],.1)
    stored=pd.read_csv(legacy.REFERENCE/'SUBSET_CANDIDATE_SCORES_FIXED.csv'); stored=stored[stored['subset']=='201+235+283']
    feats=np.concatenate([legacy.direct_features(models,legacy.geometry(estimator.states[ids[i:i+16]])[1]) for i in range(0,len(ids),16)])
    checks=[dict(check='exact_385_RC2_IDs',pass_check=set(ids)==set(stored.node_id),max_error=0.)]
    for z in (180.,200.,220.):
        tl=legacy.direct_features(models,rg,[z])[0,0]; js,_=legacy.score_features(feats,tl)
        ref=stored[stored.z_true_m==z].set_index('node_id').loc[ids]; _,kept,_=estimator.acoustic_estimate(ids,js,.5)
        error=float(np.max(abs(js-ref.J.to_numpy())))
        checks.append(dict(check=f'exact_TRIPLE_scores_survivors_z{z:g}',pass_check=error<1e-8 and set(kept)==set(ref[ref.keep_tau05].index),max_error=error))
    save('R3_LEGACY_CONTROL.csv',pd.DataFrame(checks)); assert all(x['pass_check'] for x in checks)
    catalog=pd.read_csv(legacy.OUT/'CANDIDATE_SCORE_CATALOG.csv.gz'); old=pd.read_csv(legacy.OUT/'CASE_LEVEL_RESULTS.csv'); clouds=np.load(legacy.OUT/'RC2_ACCEPTED_CLOUDS.npz'); replay=[]
    for ci,cid in enumerate(clouds['case_ids']):
        row=old[old.case_id==cid].iloc[0]; costs=estimator.bearing_costs(obs['bearing_rad'][ci])[:,0]; accepted,_=estimator.select_cloud(costs,.1)
        lo,hi=clouds['offsets'][ci:ci+2]; exact_ids=np.array_equal(accepted,clouds['node_ids'][lo:hi]); js=catalog[catalog.panel_id==row.panel_id].set_index('node_id').loc[accepted,'J'].to_numpy()
        top,kept,j=estimator.acoustic_estimate(accepted,js,.5)
        replay.append(dict(case_id=cid,pass_check=exact_ids and top==row.top1_node_id and abs(j-row.Jmin)<1e-12,top1_cell_hit=bool(row.top1_in_bracketing_cell)))
    save('LEGACY_A1_FAILURE_REPLAY.csv',pd.DataFrame(replay)); assert all(x['pass_check'] for x in replay)


def run_case(truth,pi,sigma,seed,models,model,group):
    started=time.perf_counter(); cid=f'{truth.panel_id}_s{sigma:g}_seed{seed}'
    obs=legacy.generate_observation(truth,sigma,seed,pi,models)
    cloud,cmin,cutoff,logs,sv=rc2(obs.bearing_rad,sigma)
    logs.insert(0,'case_id',cid)
    s,j,z,aclogs=acoustic_search(model,obs.bearing_rad,obs.relative_tl,sigma,cloud,cutoff); aclogs.insert(0,'case_id',cid)
    tv=truth[list(legacy.AXES)].to_numpy(float); err=legacy.errors(s,truth); region=legacy.errors(cloud,truth)
    kept=j<=j.min()+.5; floor=legacy.errors(legacy.grid(),truth).min(axis=0)
    truthcost=float(np.square(angular_residual(tv,obs.bearing_rad)).sum())
    row=dict(case_id=cid,panel_id=truth.panel_id,group=group,sigma_deg=sigma,seed=seed,Jmin=float(j[0]),bearing_cost_top1=float(np.square(angular_residual(s[0],obs.bearing_rad)).sum()),
        rc2_cmin=cmin,rc2_cutoff=cutoff,truth_bearing_cost_eval_only=truthcost,truth_inside_continuous_likelihood=bool(truthcost<=cutoff+1e-14),
        n_rc2_hypotheses=len(cloud),n_continuous_candidates=len(s),n_survivors=int(kept.sum()),n_sampled_acoustic_modes=len(s),
        rc2_starts=len(logs),rc2_convergence_fraction=float(logs.success.mean()),geometry_rank=int((sv>1e-8).sum()),geometry_smallest_singular=float(sv[-1]),
        recovered_matched_acoustic_basin=bool(j[0]<CONFIG['acoustic_basin_recovery_db']),elapsed_seconds=time.perf_counter()-started)
    for a,name in enumerate(('rel_r','abs_theta_deg','rel_v','abs_psi_deg')):
        row['top1_'+name]=float(err[0,a]); row['best_candidate_'+name]=float(err[:,a].min()); row['survivor_worst_'+name]=float(err[kept,a].max()); row['floor_'+name]=float(floor[a]); row['rc2_best_'+name]=float(region[:,a].min())
    for a,name in enumerate(legacy.AXES):
        row['top1_'+name]=float(s[0,a]); row['survivor_'+name+'_min']=float(s[kept,a].min()); row['survivor_'+name+'_max']=float(s[kept,a].max()); row['rc2_'+name+'_min']=float(cloud[:,a].min()); row['rc2_'+name+'_max']=float(cloud[:,a].max())
    cr=pd.DataFrame(s,columns=legacy.AXES); cr.insert(0,'case_id',cid); cr['J_exact']=j; cr['z_star_label_m']=z; cr['bearing_cost']=np.square(angular_residual(s,obs.bearing_rad)).sum(axis=1); cr['keep_tau05']=kept
    rc=pd.DataFrame(cloud,columns=legacy.AXES); rc.insert(0,'case_id',cid); rc['bearing_cost']=np.square(angular_residual(cloud,obs.bearing_rad)).sum(axis=1)
    print(cid,'RC2',len(cloud),'J',f'{j[0]:.6g}','errors',err[0].tolist(),'seconds',round(row['elapsed_seconds'],2),flush=True)
    return row,logs,aclogs,cr,rc


def run():
    freeze(); models={f:legacy.parse_mod(legacy.MODES/f'zgrid_f{f}.mod') for f in legacy.FREQS}
    controls(models); print('CONTINUOUS_AND_LEGACY_CONTROLS_PASS',flush=True)
    model=ContinuousAcoustic(models)
    dev=pd.read_csv(legacy.OUT/'OFFGRID_TRUTH_PANEL.csv'); hold=pd.read_csv(OUT/'HOLDOUT_TRUTH_PANEL.csv')
    allrows=[]; allrc=[]; allac=[]; allc=[]; regions=[]
    for group,panel,seeds in [('REGRESSION',dev,CONFIG['development_nominal_seeds']),('HOLDOUT',hold,CONFIG['holdout_rule']['nominal_seeds'])]:
        for pi,truth in panel.iterrows():
            for sigma,seed in [(0.,0)]+[(.1,s) for s in seeds]:
                checkpoint=OUT/'case_cache'/f'{truth.panel_id}_s{sigma:g}_seed{seed}'
                if (checkpoint/'result.json').exists():
                    row=json.loads((checkpoint/'result.json').read_text()); logs=pd.read_csv(checkpoint/'rc2.csv'); ac=pd.read_csv(checkpoint/'ac.csv'); c=pd.read_csv(checkpoint/'candidates.csv'); rc=pd.read_csv(checkpoint/'region.csv')
                    print('RESUME',row['case_id'],flush=True)
                else:
                    row,logs,ac,c,rc=run_case(truth,pi,sigma,seed,models,model,group)
                    checkpoint.mkdir(parents=True,exist_ok=True)
                    for name,df in [('rc2.csv',logs),('ac.csv',ac),('candidates.csv',c),('region.csv',rc)]:df.to_csv(checkpoint/name,index=False,float_format='%.17g')
                    legacy.json_write(checkpoint/'result.json',row)
                allrows.append(row); allrc.append(logs); allac.append(ac); allc.append(c); regions.append(rc)
                save('CASE_LEVEL_RESULTS.csv',pd.DataFrame(allrows))
    save('RC2_OPTIMIZATION_STARTS.csv',pd.concat(allrc)); save('ACOUSTIC_OPTIMIZATION_LOG.csv',pd.concat(allac)); save('CONTINUOUS_ACOUSTIC_SEARCH.csv',pd.concat(allc)); save('RC2_CONTINUOUS_HYPOTHESES.csv',pd.concat(regions))
    report()


def report():
    cases=pd.read_csv(OUT/'CASE_LEVEL_RESULTS.csv'); dev=cases[cases.group=='REGRESSION']; hold=cases[cases.group=='HOLDOUT']
    save('REGRESSION_PANEL_RESULTS.csv',dev); save('HOLDOUT_RESULTS.csv',hold)
    save('RC2_CONTINUOUS_VALIDATION.csv',cases[[c for c in cases if c.startswith(('rc2_','truth_','geometry_')) or c in ('case_id','sigma_deg','n_rc2_hypotheses','group')]])
    old=pd.read_csv(legacy.OUT/'CASE_LEVEL_RESULTS.csv'); paired=dev[dev.sigma_deg==.1].merge(old,on='case_id',suffixes=('_continuous','_legacy'))
    save('LEGACY_VS_CONTINUOUS.csv',paired)
    categories=[]
    for r in cases.itertuples():
        if r.sigma_deg==0 and r.rc2_cmin>=1e-20: mechanism='RC2_REPRESENTATION_FAILURE'
        elif not r.truth_inside_continuous_likelihood: mechanism='MIXED'
        elif not r.recovered_matched_acoustic_basin: mechanism='SEARCH_OPTIMIZATION_FAILURE'
        else: mechanism='LEGACY_REPRESENTATION_FAILURE_REPAIRED_FOR_THIS_CASE'
        categories.append(dict(case_id=r.case_id,classification=mechanism,
            legacy_regression_mechanisms='RC2_REPRESENTATION_FAILURE; ACOUSTIC_COARSE_REPRESENTATION_FAILURE' if r.group=='REGRESSION' else 'NOT_A_FROZEN_LEGACY_CASE',exact_alias_demonstrated=False,
            detail='Finite global-search result; no indistinguishable distinct exact continuous state has been demonstrated.'))
    save('FAILURE_MECHANISM_CLASSIFICATION.csv',pd.DataFrame(categories))
    noiseless=cases[cases.sigma_deg==0]
    gates={'CONTINUOUS_FORWARD_SELF_MATCH_VALIDATED':bool(pd.read_csv(OUT/'CONTINUOUS_FORWARD_SELF_MATCH.csv').pass_self_match.all()),
        'TRUTH_INDEPENDENT_RC2_CONTINUATION_VALIDATED':bool((noiseless.rc2_cmin<1e-20).all() and (cases.geometry_rank==4).all()),
        'CONTINUOUS_MULTIFREQUENCY_SEARCH_VALIDATED':bool(cases.recovered_matched_acoustic_basin.all()),
        'LEGACY_R3_CONTROL_PRESERVED':bool(pd.read_csv(OUT/'R3_LEGACY_CONTROL.csv').pass_check.all()),
        'FROZEN_OFFGRID_REGRESSION_IMPROVED':bool(dev.recovered_matched_acoustic_basin.all()),
        'HOLDOUT_OFFGRID_STRUCTURAL_GENERALIZATION_CONFIRMED':bool(hold.recovered_matched_acoustic_basin.all())}
    passed=all(gates.values()); decision='A1_FIX_PASS' if passed else 'A1_FIX_BLOCKED_BY_FINITE_CONTINUOUS_SEARCH_FAILURE'
    leakage=[dict(check='truth_free_'+name,pass_check=not any(x in __import__('inspect').signature(fn).parameters for x in ('truth','panel','oracle')),detail=str(__import__('inspect').signature(fn))) for name,fn in [('rc2',rc2),('acoustic_search',acoustic_search),('theta_profile',theta_profile)]]
    frozen=json.loads((OUT/'HOLDOUT_DESIGN_FREEZE.json').read_text())
    for name,digest in frozen['sha256'].items():leakage.append(dict(check='frozen_'+name,pass_check=legacy.sha(OUT/name)==digest,detail=digest))
    for directory in ('R3_C2_Yang_SA_depth','R3_RC23_CLOSEDLOOP','R4_A1_OFFGRID_BEARING_BOUNDARY'):
        paths=__import__('subprocess').check_output(['git','diff','--name-only','4633ff0','--','results/'+directory],text=True).splitlines()
        leakage.append(dict(check='legacy_unchanged_'+directory,pass_check=not paths,detail=str(paths)))
    save('LEAKAGE_AND_INTEGRITY_AUDIT.csv',pd.DataFrame(leakage)); assert all(x['pass_check'] for x in leakage)
    legacy.json_write(OUT/'R4_A1_FIX_DECISION.json',dict(decision=decision,gates=gates,overall_progress_percent=0,A1_status='READY_TO_RESUME_PENDING_INDEPENDENT_AUDIT' if passed else 'BLOCKED',no_engineering_tolerance_claim=True,no_exact_continuous_alias_demonstrated=True,
        scope={'V2_cases':len(cases),'regression_noiseless':9,'regression_nominal':9,'holdout_noiseless':6,'holdout_nominal':12,'V1_excluded_from_confirmation':True},
        frozen_basin_recovery_threshold_db=CONFIG['acoustic_basin_recovery_db'],nominal_cases_unrecovered=cases[(cases.sigma_deg==.1)&~cases.recovered_matched_acoustic_basin].case_id.tolist()))
    master=legacy.MASTER/'R4_PROGRESS.json'; progress=json.loads(master.read_text(encoding='utf-8')); progress.update(current_stage='A1-FIX',A1_FIX=decision,A1_FIX_gates=gates,A1_status='READY_TO_RESUME_PENDING_INDEPENDENT_AUDIT' if passed else 'BLOCKED_BY_CONTINUOUS_SEARCH',overall_progress_percent=0,independent_audit_status='A1_1_ACCEPTED; A1_FIX_COMMIT_PENDING_AUDIT')
    progress['A1_gates']['A1-1']='OFFGRID_PIPELINE_INDEPENDENT_AUDIT_PASSED'
    progress['blockers']=[] if passed else ['Finite continuous acoustic search has not recovered every matched truth basin; no physical non-identifiability conclusion is warranted.']
    progress['next_recommended_stage']='Independent A1-FIX audit before resuming full frozen A1 noise axis' if passed else 'A1-FIX search coverage/convergence review; do not open A2/B/P5'
    legacy.json_write(master,progress)
    def table(df):
        return md_table(df[['case_id','rc2_cmin','Jmin','top1_rel_r','top1_abs_theta_deg','top1_rel_v','top1_abs_psi_deg','recovered_matched_acoustic_basin']])
    text=f'''# R4-A1-FIX continuous search

**{decision}; R4 progress 0%.** This is a repair Gate, not the full sensor-accuracy experiment. R3 and structural-stop evidence remain frozen. A1-1 was independently accepted at baseline 4633ff0.

## Method and frozen scope

See METHOD_FREEZE.md and R4_A1_FIX_CONFIG.json. Starts are generated from bearing observations, the full bounded state space and independent Sobol sequences. Two global acoustic DE islands search the continuous bearing likelihood region; up to sixteen diverse spline-refined basins and up to four exact-polished basins are preserved. Collapsed populations can produce fewer separated starts; actual counts are exported. Exact modal scores determine final outputs. Neither truth coordinates nor truth-cell vertices enter either estimator API. Continuous source depth is not estimated: the shared inherited nuisance profile is unchanged.

The six holdout states follow seed 2026100203, fixed interior ranges, rejection of on-grid coordinates and the frozen generation rule. Their panel and method hash were frozen before any holdout observation. Development: 9 truths x (noiseless + one nominal seed). Holdout: 6 truths x (noiseless + two nominal seeds). A V1 theta-offset initialization defect was corrected and a fresh holdout frozen; V1 evidence remains in INITIAL_IMPLEMENTATION. Search budgets and scientific thresholds were not changed to improve panel success. No full multi-sigma MC is run.

## Continuous self-match and legacy controls

{md_table(pd.read_csv(OUT/'CONTINUOUS_FORWARD_SELF_MATCH.csv'))}

Tolerance: bearing cost <1e-20 rad squared; exact profiled acoustic J <1e-10 dB. Truth evaluation is an oracle forward-integrity control only. R3's exact 385-node IDs and all TRIPLE scores/survivors at three depths are checked against frozen evidence. All 27 legacy A1 failure cases are replayed from frozen observations/catalogs; legacy scoring remains BASELINE_LEGACY.

## Regression results

{table(dev)}

## Holdout results

{table(hold)}

## Interpretation and limits

Gate flags: `{json.dumps(gates)}`.

Recovered exact acoustic basins: regression noiseless {int(dev[dev.sigma_deg==0].recovered_matched_acoustic_basin.sum())}/9, nominal {int(dev[dev.sigma_deg==.1].recovered_matched_acoustic_basin.sum())}/9; holdout noiseless {int(hold[hold.sigma_deg==0].recovered_matched_acoustic_basin.sum())}/6, nominal {int(hold[hold.sigma_deg==.1].recovered_matched_acoustic_basin.sum())}/12. Continuous-truth bearing-likelihood coverage is {int(cases.truth_inside_continuous_likelihood.sum())}/{len(cases)}. These fixed small panels are structural tests, not Monte Carlo accuracy boundaries.

Noiseless bearing solutions use independently full-rank Cartesian geometry and bounded nonlinear multistart refinement; distance/region extent and all convergence starts are exported. For noisy cases the cutoff is unchanged, but acts on exact continuous states. RC2 sampling and acoustic global search remain finite and do not certify exhaustive modes or a confidence envelope. Survivor envelopes describe exported candidates only. DE reaching its fixed generation budget is explicitly recorded, not counted as optimizer convergence. Local convergence fractions are in raw logs.

`n_sampled_acoustic_modes` in raw case data is a finite separated-candidate proxy at the frozen deduplication resolution, not a certified count of continuous local or physical modes. MODE_INVENTORY.csv distinguishes the certified unique noiseless bearing solution from non-exhaustive noisy hypothesis clusters. DE `nfev` denotes vectorized objective batches, not individual state evaluations.

An exact synthetic truth has near-zero acoustic score, so a positive optimized minimum, while truth remains inside the bearing region, is evidence of finite-search failure. It does not establish physical aliasing. No two distinct continuous states with indistinguishable observations have been demonstrated. A recovered acoustic basin requires exact J <0.001 dB, not just an improvement in top1 error.

CASE_LEVEL_RESULTS.csv reports angular circular errors, relative r/v errors, coordinatewise best-candidate errors, candidate envelopes, bearing/acoustic costs, sampled basin counts and convergence. Coarse quantization floors are reference only. HOLDOUT_RESULTS.csv is structurally confirmatory for this fixed interior design, not an engineering bearing tolerance or general ocean guarantee.

## Reproduce

`python r4_a1_fix_continuous.py run` (checkpoints allow deterministic resume); `python r4_a1_fix_continuous.py report` regenerates summaries. Forward interpolation is independently checked at fixed random ranges and final candidates use exact scoring. API signatures, frozen hashes and legacy isolation are locally checked; independent scientific audit remains pending.

`python r4_a1_fix_audit.py` reconstructs all final acoustic scores, bearing costs, errors and candidate envelopes from separately regenerated observations, verifies legacy challenge bytes, exports optimizer starts/convergence and checks exact-alias diagnostics. `python r4_a1_fix_figures.py` creates scientific figures from saved data only. See LOCAL_VALIDATION.md for final check counts and limitations.

Recommended next Gate if blocked: strengthen truth-independent acoustic basin coverage and verify global-search convergence, for example with physically resolved hierarchical range/profile searches or independent global solvers. Preserve V2 regression/holdout failures as development evidence and preregister a new confirmation panel after any method change. Do not inflate the bearing threshold or begin full A1 statistics before repair validation.
'''
    (OUT/'R4_A1_FIX_REPORT.md').write_text(text,encoding='utf-8'); (OUT/'GPT_SYNC.md').write_text(text,encoding='utf-8')
    print(decision,gates,flush=True)


def main():
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['freeze','run','report']); args=p.parse_args()
    with threadpool_limits(limits=2):
        {'freeze':freeze,'run':run,'report':report}[args.action]()


if __name__=='__main__':main()
