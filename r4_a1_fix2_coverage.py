"""FIX2: development basin diagnostics and observation-driven acoustic coverage."""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq, least_squares, minimize
from scipy.stats import qmc
from threadpoolctl import threadpool_limits
import r4_a1_fix_continuous as prior
L=prior.legacy
OUT=L.ROOT/'results/R4_A1_FIX2_ACOUSTIC_COVERAGE'
SEARCH_CONFIG={
    'baseline_commit':'956f2dfba8e719561641fd135f246a650fe76dca',
    'architecture':'deterministic nested range/radial-velocity bank, bearing continuation, depth-branch spatial beams, analytic-Jacobian local refinement, exact verification',
    'budgets':[{'name':'B1','r_step_m':40.,'u_step_mps':.04,'starts_per_depth':3},
               {'name':'B2','r_step_m':20.,'u_step_mps':.02,'starts_per_depth':6},
               {'name':'B3','r_step_m':10.,'u_step_mps':.01,'starts_per_depth':12}],
    'range_bounds_m':[45000.,60000.],'radial_velocity_bounds_mps':[.94,3.],
    'bearing_profile_GN_iterations':5,'batch_size':512,'beam_pool_per_depth':1024,
    'beam_max_range_separation_m':16.,'local_max_nfev':160,'exact_polish_starts':8,'exact_polish_max_nfev':120,
    'local_tolerances':{'xtol':1e-12,'ftol':1e-12,'gtol':1e-12},
    'basin_gate_db':.001,'candidate_tau_db':.5,'dedup_resolution':prior.CONFIG['candidate_dedup_resolution'],
    'independent_solver':{'family':'independent Sobol radial bank + derivative-free Nelder-Mead',
        'sobol_power':18,'seed':773003,'starts_per_depth':4,'maxiter':2500,'xatol':1e-9,'fatol':1e-14},
    'hardest_cases':['P09_s0.1_seed410001','P06_s0.1_seed410001','P04_s0.1_seed410001','Q01_s0.1_seed420001','Q04_s0.1_seed420001'],
    'convergence_state_tolerances':[.00001,.001,.0001,.01],
    'independent_agreement_state_tolerances':[.00001,.001,.0001,.01],
    'fresh_holdout':{'seed':2026100204,'n':8,'prefix':'F','nominal_seeds':[430001,430002],
        'low':[47.,-1.2,1.5,-9.],'high':[55.5,1.2,2.5,9.],'depths':[180.,200.,220.]}}


def preregister_budgets():
    OUT.mkdir(exist_ok=True)
    p=OUT/'SEARCH_BUDGET_PREREGISTRATION.json'
    if p.exists():assert json.loads(p.read_text())==SEARCH_CONFIG
    else:L.json_write(p,SEARCH_CONFIG)


def radial_profile(r_m,u,bearing):
    """Observation-only theta/tangential continuation at fixed range/radial speed."""
    r_m=np.atleast_1d(r_m).astype(float); u=np.broadcast_to(u,r_m.shape).astype(float)
    theta=np.full(len(r_m),np.mean(bearing)); w=np.zeros(len(r_m)); t=L.TIMES
    xp,yp=L.platform_xy()
    for _ in range(SEARCH_CONFIG['bearing_profile_GN_iterations']):
        ct,st=np.cos(theta)[:,None],np.sin(theta)[:,None]
        longitudinal=r_m[:,None]+u[:,None]*t; transverse=w[:,None]*t
        dx=longitudinal*ct-transverse*st-xp; dy=longitudinal*st+transverse*ct-yp
        rr=dx*dx+dy*dy; residual=L.wrap_rad(np.arctan2(dy,dx)-bearing)
        jt=(dx*(longitudinal*ct-transverse*st)-dy*(-longitudinal*st-transverse*ct))/rr
        jw=t*(dx*ct+dy*st)/rr
        aa=np.sum(jt*jt,axis=1); bb=np.sum(jt*jw,axis=1); cc=np.sum(jw*jw,axis=1)
        g1=np.sum(jt*residual,axis=1); g2=np.sum(jw*residual,axis=1); det=aa*cc-bb*bb
        theta-= (cc*g1-bb*g2)/det; w-= (aa*g2-bb*g1)/det
        theta=np.clip(theta,np.radians(-5),np.radians(5))
        cap=np.sqrt(np.maximum(9-u*u,0)); lo=np.maximum(u*np.tan(np.radians(-15)-theta),-cap); hi=np.minimum(u*np.tan(np.radians(15)-theta),cap)
        w=np.clip(w,lo,hi)
        minimum=np.sqrt(np.maximum(1-u*u,0)); w=np.where(abs(w)<minimum,np.where(w>=0,minimum,-minimum),w)
    v=np.hypot(u,w); psi=theta+np.arctan2(w,u)
    return np.column_stack((r_m/1000,np.degrees(theta),v,np.degrees(psi)))


def trajectory_jacobian(state):
    r,theta,v,psi=state; theta,psi=np.radians(theta),np.radians(psi); t=L.TIMES
    xp,yp=L.platform_xy(); dx=1000*r*np.cos(theta)+v*t*np.cos(psi)-xp; dy=1000*r*np.sin(theta)+v*t*np.sin(psi)-yp
    ranges=np.hypot(dx,dy)
    dxds=np.column_stack((np.full(121,1000*np.cos(theta)),np.full(121,-1000*r*np.sin(theta)*np.pi/180),t*np.cos(psi),-v*t*np.sin(psi)*np.pi/180))
    dyds=np.column_stack((np.full(121,1000*np.sin(theta)),np.full(121,1000*r*np.cos(theta)*np.pi/180),t*np.sin(psi),v*t*np.cos(psi)*np.pi/180))
    range_jac=(dx[:,None]*dxds+dy[:,None]*dyds)/ranges[:,None]
    bearing_jac=(dx[:,None]*dyds-dy[:,None]*dxds)/(ranges[:,None]**2)
    return np.arctan2(dy,dx),ranges,range_jac,bearing_jac


def feature_jacobian(model,state,depth,exact=False):
    _,ranges,range_jac,_=trajectory_jacobian(state); values=[]; jac=[]
    for frequency in L.FREQS:
        if exact:
            mod=model.models[frequency]; iz,_=L.effective_depths(mod,[depth]); izr,_=L.effective_depths(mod,[200.])
            k=mod['k'].real.astype(float); alpha=-mod['k'].imag.astype(float)
            weights=mod['phi'][iz[0]]*mod['phi'][izr[0]]
            factors=np.sqrt(2*np.pi/(k[:,None]*ranges))*np.exp(-1j*k[:,None]*ranges-alpha[:,None]*ranges-1j*np.pi/4)
            p=weights@factors; dp=weights@(factors*(-.5/ranges-1j*k[:,None]-alpha[:,None]))
        else:
            iz=int(np.flatnonzero(L.PROFILE==depth)[0]); spline=model.depth_splines[frequency,iz]
            p=spline(ranges); dp=spline(ranges,1)
        tl=20*np.log10(np.maximum(abs(p),1e-30)); derivative=20/np.log(10)*np.real(dp/p)
        j=derivative[:,None]*range_jac
        for sl in (slice(0,61),slice(61,None)):
            tl[sl]-=tl[sl].mean(); j[sl]-=j[sl].mean(axis=0)
        values.append(tl); jac.append(j)
    return np.array(values),np.array(jac)


def refine(model,start,depth,bearing,tl,sigma,cutoff,exact=False):
    cached={}
    def both(u):
        key=u.tobytes()
        if key!=cached.get('key'):
            s=prior.LOW+prior.WIDTH*u; f,j=feature_jacobian(model,s,depth,exact)
            predicted,_,_,bj=trajectory_jacobian(s); res=L.wrap_rad(predicted-bearing); cost=float(res@res)
            excess=max(cost-cutoff-1e-14,0); penalty=10*np.sqrt(excess)/np.radians(sigma)
            pg=10*(res@bj)/np.radians(sigma)/np.sqrt(excess) if excess else np.zeros(4)
            residual=np.r_[(f-tl).ravel()/np.sqrt(tl.size),penalty]
            jac=np.vstack((j.reshape(-1,4)/np.sqrt(tl.size),pg))*prior.WIDTH
            cached.update(key=key,residual=residual,jac=jac)
        return cached['residual'],cached['jac']
    fit=least_squares(lambda u:both(u)[0],np.clip((start-prior.LOW)/prior.WIDTH,1e-12,1-1e-12),jac=lambda u:both(u)[1],
        bounds=(np.zeros(4),np.ones(4)),x_scale='jac',max_nfev=SEARCH_CONFIG['exact_polish_max_nfev'] if exact else SEARCH_CONFIG['local_max_nfev'],**SEARCH_CONFIG['local_tolerances'])
    return prior.LOW+prior.WIDTH*fit.x,dict(success=bool(fit.success),nfev=int(fit.nfev),optimality=float(fit.optimality),status=int(fit.status))


def bank(model,bearing,tl,cutoff,kind='GRID',power=None):
    if kind=='GRID':
        r=np.arange(45000,60001,10.); u=np.round(np.arange(.94,3.001,.01),10)
        ri,ui=np.meshgrid(np.arange(len(r)),np.arange(len(u)),indexing='ij'); ri,ui=ri.ravel(),ui.ravel()
        rr,uu=r[ri],u[ui]
    else:
        points=qmc.Sobol(2,scramble=True,seed=SEARCH_CONFIG['independent_solver']['seed']).random_base2(power)
        rr=45000+15000*points[:,0]; uu=.94+2.06*points[:,1]; ri=ui=np.arange(len(rr))
    states=[]; scores=[]; cost_out=[]; ids=[]; last=time.perf_counter(); started=last
    for st in range(0,len(rr),SEARCH_CONFIG['batch_size']):
        s=radial_profile(rr[st:st+512],uu[st:st+512],bearing)
        costs=np.square(prior.angular_residual(s,bearing)).sum(axis=1)
        valid=(costs<=cutoff)&np.all(s>=prior.LOW-1e-12,axis=1)&np.all(s<=prior.HIGH+1e-12,axis=1)
        if valid.any():
            selected=s[valid]; features=model.features(selected); j=np.sqrt(np.square(features-tl[None,None,:,:]).mean(axis=(2,3)))
            states.append(selected); scores.append(j); cost_out.append(costs[valid]); ids.append(st+np.flatnonzero(valid))
        if time.perf_counter()-last>35:
            print('BANK',kind,f'{st+512}/{len(rr)}','elapsed',round(time.perf_counter()-started,1),flush=True); last=time.perf_counter()
    ids=np.concatenate(ids)
    return np.concatenate(states),np.concatenate(scores),np.concatenate(cost_out),ri[ids],ui[ids],dict(total_bank_nodes=len(rr),feasible_bank_nodes=len(ids),elapsed_seconds=time.perf_counter()-started)


def select_beams(states,scores,mask,per_depth):
    result=[]; indices=np.flatnonzero(mask)
    t=L.TIMES[[0,60,120]]; xp,yp=L.platform_xy(); xp=xp[[0,60,120]];yp=yp[[0,60,120]]
    r,theta,v,psi=states.T; theta=np.radians(theta);psi=np.radians(psi)
    sample=np.hypot(1000*r[:,None]*np.cos(theta)[:,None]+v[:,None]*t*np.cos(psi)[:,None]-xp,
                    1000*r[:,None]*np.sin(theta)[:,None]+v[:,None]*t*np.sin(psi)[:,None]-yp)
    for zi,z in enumerate(L.PROFILE):
        subset=indices[np.argsort(scores[indices,zi])[:SEARCH_CONFIG['beam_pool_per_depth']]]; chosen=[]
        for i in subset:
            if not chosen or all(np.max(abs(sample[i]-sample[j]))>=SEARCH_CONFIG['beam_max_range_separation_m'] for j in chosen):chosen.append(i)
            if len(chosen)>=per_depth:break
        result.extend((int(i),float(z),float(scores[i,zi])) for i in chosen)
    return result


def one_budget(model,bearing,tl,sigma,cutoff,states,scores,mask,budget,cid):
    beams=select_beams(states,scores,mask,budget['starts_per_depth']); refined=[]; logs=[]; starts=[]
    for i,(index,depth,initial_j) in enumerate(beams):
        start=states[index]; s,info=refine(model,start,depth,bearing,tl,sigma,cutoff)
        refined.append(s); logs.append(dict(case_id=cid,budget=budget['name'],stage='LOCAL_SURROGATE',depth_label_m=depth,**info))
        starts.append(dict(case_id=cid,budget=budget['name'],start_id=i,depth_label_m=depth,proposal_J=initial_j,**dict(zip(L.AXES,start))))
    refined=np.concatenate((np.array(refined),states[[i for i,_,_ in beams]]))
    scores_exact,depths=model.exact(refined,tl)
    # Polish spatially distinct basins, not merely the best candidate.
    ranked=np.argsort(scores_exact); unique=prior.unique_indices(refined[ranked]); chosen=ranked[unique][:SEARCH_CONFIG['exact_polish_starts']]
    for i in chosen:
        s,info=refine(model,refined[i],depths[i],bearing,tl,sigma,cutoff,exact=True)
        refined=np.vstack((refined,s)); logs.append(dict(case_id=cid,budget=budget['name'],stage='LOCAL_EXACT',depth_label_m=depths[i],**info))
    costs=np.square(prior.angular_residual(refined,bearing)).sum(axis=1); refined=refined[costs<=cutoff+1e-14]
    j,z=model.exact(refined,tl); order=np.argsort(j); refined,j,z=refined[order],j[order],z[order]
    keep=prior.unique_indices(refined); refined,j,z=refined[keep],j[keep],z[keep]
    frame=pd.DataFrame(refined,columns=L.AXES); frame['J_exact']=j; frame['z_star_label_m']=z; frame['bearing_cost']=np.square(prior.angular_residual(refined,bearing)).sum(axis=1)
    frame.insert(0,'budget',budget['name']); frame.insert(0,'case_id',cid); frame['keep_tau05']=j<=j.min()+.5
    return frame,pd.DataFrame(logs),pd.DataFrame(starts)


def cumulative_candidates(frames,budget_names):
    """Truth-free nested mode bank: never discard a basin found at lower budgets."""
    result=[]
    for k,name in enumerate(budget_names):
        pooled=frames[frames.budget.isin(budget_names[:k+1])].copy().sort_values('J_exact',kind='stable')
        states=pooled[list(L.AXES)].to_numpy();pooled=pooled.iloc[prior.unique_indices(states)].copy()
        pooled['origin_budget']=pooled.budget;pooled['budget']=name;pooled['keep_tau05']=pooled.J_exact<=pooled.J_exact.min()+.5
        result.append(pooled)
    return pd.concat(result,ignore_index=True)


def close_development_pools():
    raw=[];frames=[]
    for cache in sorted((OUT/'development_cache').iterdir()):
        if (cache/'complete.json').exists():raw.append(pd.read_csv(cache/'results.csv'));frames.append(pd.read_csv(cache/'candidates.csv'))
    raw=pd.concat(raw,ignore_index=True);frames=pd.concat(frames,ignore_index=True)
    panels=pd.concat([pd.read_csv(L.OUT/'OFFGRID_TRUTH_PANEL.csv'),pd.read_csv(prior.OUT/'HOLDOUT_TRUTH_PANEL.csv')]).set_index('panel_id')
    rows=[];pooled=[]
    for cid,group in frames.groupby('case_id',sort=False):
        p=cumulative_candidates(group,[x['name'] for x in SEARCH_CONFIG['budgets']]);pooled.append(p)
        for budget,bundle in p.groupby('budget',sort=False):
            r=raw[(raw.case_id==cid)&(raw.budget==budget)].iloc[0].to_dict();r['raw_component_exact_J']=r['exact_J'];r['raw_component_recovered']=r['recovered']
            top=bundle.iloc[0];truth=panels.loc[r['panel_id']];e=L.errors(bundle[list(L.AXES)].to_numpy(),truth)
            r.update(exact_J=float(top.J_exact),recovered=bool(top.J_exact<.001),bearing_cost=float(top.bearing_cost),z_star_label_m=float(top.z_star_label_m),n_separated_candidates=len(bundle),origin_budget_of_best=top.origin_budget)
            for a in L.AXES:r[a]=float(top[a])
            for k,n in enumerate(('rel_r','abs_theta_deg','rel_v','abs_psi_deg')):r['top1_'+n]=float(e[0,k]);r['survivor_worst_'+n]=float(e[bundle.keep_tau05.to_numpy(),k].max())
            rows.append(r)
    save('RAW_BUDGET_COMPONENT_RESULTS.csv',raw);save('RAW_DEVELOPMENT_CANDIDATES.csv',frames)
    save('SEARCH_BUDGET_CONVERGENCE.csv',pd.DataFrame(rows));save('DEVELOPMENT_RECOVERY.csv',pd.DataFrame(rows).query("budget=='B3'"));save('DEVELOPMENT_CANDIDATES.csv',pd.concat(pooled))


def development():
    preregister_budgets(); protect(); model=Acoustic(inputs())
    cases=pd.read_csv(prior.OUT/'CASE_LEVEL_RESULTS.csv'); obs=np.load(prior.OUT/'CASE_OBSERVATIONS.npz'); panels=pd.concat([pd.read_csv(L.OUT/'OFFGRID_TRUTH_PANEL.csv'),pd.read_csv(prior.OUT/'HOLDOUT_TRUTH_PANEL.csv')]).set_index('panel_id')
    allresults=[]; allcandidates=[]; alllogs=[]; allstarts=[]; coverage=[]
    for ci,row in cases.iterrows():
        cache=OUT/'development_cache'/row.case_id; bearing=obs['bearing_rad'][ci]; tl=obs['relative_tl'][ci]
        if (cache/'complete.json').exists():
            allresults.append(pd.read_csv(cache/'results.csv')); allcandidates.append(pd.read_csv(cache/'candidates.csv')); alllogs.append(pd.read_csv(cache/'logs.csv')); allstarts.append(pd.read_csv(cache/'starts.csv')); coverage.append(json.loads((cache/'complete.json').read_text()))
            print('RESUME_DEVELOPMENT',row.case_id,flush=True); continue
        start_time=time.perf_counter(); result=[]; cf=[]; lf=[]; sf=[]
        if row.sigma_deg==0:
            cloud,cmin,cutoff,logs,_=prior.rc2(bearing,0); j,z=model.exact(cloud,tl)
            states=cloud; scores=None; cov=dict(total_bank_nodes=len(cloud),feasible_bank_nodes=len(cloud),elapsed_seconds=0.)
        else:
            states,scores,_,ri,ui,cov=bank(model,bearing,tl,row.rc2_cutoff)
        for budget in SEARCH_CONFIG['budgets']:
            tick=time.perf_counter()
            if row.sigma_deg==0:
                c=pd.DataFrame(states,columns=L.AXES); c['J_exact']=j; c['z_star_label_m']=z; c['bearing_cost']=np.square(prior.angular_residual(states,bearing)).sum(axis=1); c['case_id']=row.case_id;c['budget']=budget['name'];c['keep_tau05']=j<=j.min()+.5
                logs=pd.DataFrame([dict(case_id=row.case_id,budget=budget['name'],stage='NOISELESS_RC2_UNIQUE',success=True,nfev=len(states))]); starts=pd.DataFrame([dict(case_id=row.case_id,budget=budget['name'],start_id=0,depth_label_m=z[0],proposal_J=j[0],**dict(zip(L.AXES,states[0])))])
                n_nodes=len(states)
            else:
                stride=int(budget['r_step_m']/10); mask=((ri%stride==0)|(ri==1500))&((ui%stride==0)|(ui==206))
                c,logs,starts=one_budget(model,bearing,tl,row.sigma_deg,row.rc2_cutoff,states,scores,mask,budget,row.case_id); n_nodes=int(mask.sum())
            truth=panels.loc[row.panel_id]; error=L.errors(c[list(L.AXES)].to_numpy(),truth); top=c.iloc[0]
            r=dict(case_id=row.case_id,panel_id=row.panel_id,origin_group=row.group,budget=budget['name'],sigma_deg=row.sigma_deg,seed=row.seed,
                exact_J=float(top.J_exact),bearing_cost=float(top.bearing_cost),cutoff=float(row.rc2_cutoff),feasible=True,recovered=bool(top.J_exact<.001),
                n_bank_feasible_nodes=n_nodes,n_separated_candidates=len(c),n_local_starts=len(starts),local_convergence_fraction=float(logs.success.mean()),
                budget_r_step_m=budget['r_step_m'],budget_u_step_mps=budget['u_step_mps'],elapsed_local_seconds=time.perf_counter()-tick,
                z_star_label_m=float(top.z_star_label_m),**{a:float(top[a]) for a in L.AXES})
            for k,name in enumerate(('rel_r','abs_theta_deg','rel_v','abs_psi_deg')):r['top1_'+name]=float(error[0,k]); r['survivor_worst_'+name]=float(error[c.keep_tau05.to_numpy(),k].max())
            result.append(r); cf.append(c); lf.append(logs); sf.append(starts)
            print('DEVELOPMENT',row.case_id,budget['name'],'J',r['exact_J'],'recovered',r['recovered'],'local seconds',round(r['elapsed_local_seconds'],1),flush=True)
        cache.mkdir(parents=True,exist_ok=True)
        for name,df in [('results.csv',pd.DataFrame(result)),('candidates.csv',pd.concat(cf)),('logs.csv',pd.concat(lf)),('starts.csv',pd.concat(sf))]:df.to_csv(cache/name,index=False,float_format='%.17g')
        cov.update(case_id=row.case_id,elapsed_total_seconds=time.perf_counter()-start_time); L.json_write(cache/'complete.json',cov)
        allresults.append(pd.DataFrame(result)); allcandidates.append(pd.concat(cf)); alllogs.append(pd.concat(lf)); allstarts.append(pd.concat(sf)); coverage.append(cov)
        save('SEARCH_BUDGET_CONVERGENCE.csv',pd.concat(allresults)); save('DEVELOPMENT_RECOVERY.csv',pd.concat(allresults).query("budget=='B3'")); save('GLOBAL_BANK_COVERAGE.csv',pd.DataFrame(coverage))
    save('DEVELOPMENT_CANDIDATES.csv',pd.concat(allcandidates)); save('LOCAL_CONVERGENCE_LOG.csv',pd.concat(alllogs)); save('OBSERVATION_GENERATED_STARTS.csv',pd.concat(allstarts))
    close_development_pools()


def independent_search(model,bearing,tl,sigma,cutoff):
    """Separate low-discrepancy global initialization and derivative-free family."""
    cfg=SEARCH_CONFIG['independent_solver']; states,scores,_,_,_,coverage=bank(model,bearing,tl,cutoff,'SOBOL',cfg['sobol_power'])
    beams=select_beams(states,scores,np.ones(len(states),bool),cfg['starts_per_depth']); proposals=[]; solutions=[]; logs=[]
    # Local radial coordinates remove the dominant r/v phase coupling.
    scale=np.array([10.,.01,.1,.1])
    def radial(s):
        theta=np.radians(s[1]); angle=np.radians(s[3])-theta
        return np.array([1000*s[0],s[2]*np.cos(angle),s[2]*np.sin(angle),s[1]])
    def original(a):return np.array([a[0]/1000,a[3],np.hypot(a[1],a[2]),a[3]+np.degrees(np.arctan2(a[2],a[1]))])
    for ii,(index,z,initial_j) in enumerate(beams):
        start=states[index]; a0=radial(start)
        def objective(x):
            s=original(a0+x*scale)
            if (s<prior.LOW).any() or (s>prior.HIGH).any():return 1e6+float(np.square(np.maximum(prior.LOW-s,0)+np.maximum(s-prior.HIGH,0)).sum())
            f=feature_jacobian(model,s,z)[0]; acoustic=float(np.square(f-tl).mean())
            excess=max(float(np.square(prior.angular_residual(s,bearing)).sum())-cutoff-1e-14,0)
            return acoustic+100*excess/np.radians(sigma)**2
        fit=minimize(objective,np.zeros(4),method='Nelder-Mead',options={'maxiter':cfg['maxiter'],'xatol':cfg['xatol'],'fatol':cfg['fatol'],'initial_simplex':np.vstack((np.zeros(4),np.eye(4)))})
        s=original(a0+fit.x*scale); proposals.append(dict(start_id=ii,depth_label_m=z,proposal_J=initial_j,**dict(zip(L.AXES,start))))
        solutions.extend((s,start)); logs.append(dict(start_id=ii,success=bool(fit.success),nfev=int(fit.nfev),nit=int(fit.nit),objective=float(fit.fun)))
    solutions=np.array(solutions); cost=np.square(prior.angular_residual(solutions,bearing)).sum(axis=1)
    solutions=solutions[(cost<=cutoff+1e-14)&np.all(solutions>=prior.LOW-1e-12,axis=1)&np.all(solutions<=prior.HIGH+1e-12,axis=1)]
    j,z=model.exact(solutions,tl); order=np.argsort(j); solutions,j,z=solutions[order],j[order],z[order]
    keep=prior.unique_indices(solutions); solutions,j,z=solutions[keep],j[keep],z[keep]
    frame=pd.DataFrame(solutions,columns=L.AXES);frame['J_exact']=j;frame['z_star_label_m']=z;frame['bearing_cost']=np.square(prior.angular_residual(solutions,bearing)).sum(axis=1)
    return frame,pd.DataFrame(proposals),pd.DataFrame(logs),coverage


def independent():
    preregister_budgets(); protect(); model=Acoustic(inputs()); observations=np.load(prior.OUT/'CASE_OBSERVATIONS.npz'); cases=pd.read_csv(prior.OUT/'CASE_LEVEL_RESULTS.csv').set_index('case_id')
    ids=list(observations['case_ids']); records=[]; finals=[]; logs=[]; proposals=[]
    supplement=OUT/'INDEPENDENT_SUPPLEMENT_DESIGN.json'
    failed=pd.read_csv(OUT/'DEVELOPMENT_RECOVERY.csv').query('not recovered').case_id.tolist()
    extra=[cid for cid in failed if cid not in SEARCH_CONFIG['hardest_cases']]
    design=dict(reason='Additional failures exposed by completed development search; diagnostic only, unchanged independent family/budget',case_ids=extra,solver=SEARCH_CONFIG['independent_solver'])
    if supplement.exists():assert json.loads(supplement.read_text())==design
    else:L.json_write(supplement,design)
    for cid in SEARCH_CONFIG['hardest_cases']+extra:
        cache=OUT/'independent_cache'/cid
        if (cache/'result.json').exists():
            info=json.loads((cache/'result.json').read_text());frame=pd.read_csv(cache/'candidates.csv');starts=pd.read_csv(cache/'starts.csv');log=pd.read_csv(cache/'logs.csv')
        else:
            ci=ids.index(cid); row=cases.loc[cid]
            frame,starts,log,coverage=independent_search(model,observations['bearing_rad'][ci],observations['relative_tl'][ci],row.sigma_deg,row.rc2_cutoff)
            main=pd.read_csv(OUT/'DEVELOPMENT_RECOVERY.csv').set_index('case_id').loc[cid];top=frame.iloc[0]
            diff=abs(top[list(L.AXES)].to_numpy(float)-main[list(L.AXES)].to_numpy(float))
            for a in (1,3):diff[a]=float(L.angle_error_deg(top[L.AXES[a]],main[L.AXES[a]]))
            agreement=bool(top.J_exact<.001 and main.exact_J<.001 and np.all(diff<=SEARCH_CONFIG['independent_agreement_state_tolerances']) and top.z_star_label_m==main.z_star_label_m and top.bearing_cost<=main.cutoff+1e-14)
            info=dict(case_id=cid,exact_J=float(top.J_exact),primary_exact_J=float(main.exact_J),bearing_cost=float(top.bearing_cost),cutoff=float(main.cutoff),z_star_label_m=float(top.z_star_label_m),agreement=agreement,
                differences=json.dumps(diff.tolist()),**{a:float(top[a]) for a in L.AXES},**coverage)
            cache.mkdir(parents=True,exist_ok=True)
            for name,df in [('candidates.csv',frame),('starts.csv',starts),('logs.csv',log)]:df.to_csv(cache/name,index=False,float_format='%.17g')
            L.json_write(cache/'result.json',info)
        for df in (frame,starts,log):df.insert(0,'case_id',cid)
        records.append(info);finals.append(frame);proposals.append(starts);logs.append(log)
        save('INDEPENDENT_SOLVER_AGREEMENT.csv',pd.DataFrame(records))
        print('INDEPENDENT',cid,'J',info['exact_J'],'agreement',info['agreement'],flush=True)
    save('INDEPENDENT_CANDIDATES.csv',pd.concat(finals));save('INDEPENDENT_STARTS.csv',pd.concat(proposals));save('INDEPENDENT_OPTIMIZATION_LOG.csv',pd.concat(logs))


def freeze_method():
    protect(); preregister_budgets()
    # Confirmation starts only after the development and independent comparison finish.
    assert len(pd.read_csv(OUT/'DEVELOPMENT_RECOVERY.csv'))==36
    extra=json.loads((OUT/'INDEPENDENT_SUPPLEMENT_DESIGN.json').read_text())['case_ids']
    assert len(pd.read_csv(OUT/'INDEPENDENT_SOLVER_AGREEMENT.csv'))==5+len(extra)
    method=OUT/'FIX2_METHOD.md'
    if not method.exists():
        method.write_text('''# FIX2 method frozen before fresh confirmation

The accepted continuous RC2 is reused unchanged. Global proposals partition initial range and target radial speed with deterministic nested 40/20/10 m by 0.04/0.02/0.01 m/s meshes. Bearing observations continue theta and tangential speed through bounded Gauss-Newton profiling; only bearing-feasible, physically bounded proposals are scored. All distinct verified candidates from lower budgets are preserved in a cumulative mode bank. Raw component results are also reported separately, so retention alone is not mistaken for independent fine-resolution convergence. This is a finite coverage-oriented search, not a certified exhaustive four-dimensional branch-and-bound.

Every inherited nuisance-depth branch retains 3/6/12 spatially separated starts. Distinct proposals differ by at least 16 m in range at one of t=0,600,1200 s, approximately half the fastest modal-interference spatial cycle. Depth is a shared profiled nuisance, not a new estimated state. Continuous pressure splines generate proposals, with analytic derivatives and Jacobian scaling for refinement. The final scores and eight final-polish branches use the exact modal forward model. No oracle state, error or basin-diagnostic coordinate enters search.

Three increasing budgets are preregistered in SEARCH_BUDGET_PREREGISTRATION.json. Raw components are evaluated separately and their verified candidates retained cumulatively. B3 confirmation executes the complete B1/B2/B3 hierarchy and returns its cumulative mode pool. The finest bank is calculated once and its nested subsets reused without allowing excluded points into lower-budget selection. Node counts therefore describe effective logical resolutions; shared-cache elapsed time is not a separate B1 runtime measurement. Known regression and former Q holdout cases are development evidence. Oracle landscape measurements characterize feasible basin widths and coupled sensitivity; they never supply proposal coordinates or initialization.

Independent validation uses a separate Sobol range/radial bank with 2^18 points, and a derivative-free Nelder-Mead family in radial coordinates. It does not initialize from the primary recovered state. Five preselected difficult old failures plus additional failed development cases recorded in INDEPENDENT_SUPPLEMENT_DESIGN.json are compared on exact score, all coordinates, shared nuisance depth and bearing feasibility. The supplementary checks are explicitly selected after development outcomes, before fresh confirmation, with the same independent solver/budget.

The direct-modal matched-basin threshold stays J<0.001 dB. Finite candidate envelopes and basin counts are not confidence guarantees or exhaustive alias certificates. No engineering bearing tolerance is claimed. R4 progress remains zero even if this repair Gate passes.

Fresh confirmation: eight deterministic interior off-grid truths, seed 2026100204, all four coordinates off-grid; noiseless plus two nominal 0.1 degree seeds 430001/430002. The method/config/code hashes are frozen before generating this panel; its design/panel hash is separately frozen before observations. No manual editing or post-confirmation tuning is permitted. Independent cases may execute in four process shards (panel index modulo four), each with two BLAS threads; aggregation replays the same case caches in the original panel order. Sharding changes scheduling only, never observations, budgets, selection or search results.
''',encoding='utf-8')
    code=Path(__file__); paths=[code,method,OUT/'SEARCH_BUDGET_PREREGISTRATION.json',OUT/'INDEPENDENT_SUPPLEMENT_DESIGN.json',OUT/'FROZEN_INPUT_HASHES.csv']
    manifest=OUT/'METHOD_FREEZE.json'; hashes={p.relative_to(L.ROOT).as_posix():L.sha(p) for p in paths}
    if manifest.exists():assert json.loads(manifest.read_text())['sha256']==hashes,'Frozen method/code changed'
    else:
        from datetime import datetime,timezone
        L.json_write(manifest,dict(frozen_utc=datetime.now(timezone.utc).isoformat(),sha256=hashes,search_budget='B3',recovery_threshold_db=.001,before_new_holdout_panel_generation=True))
    hold=OUT/'FRESH_HOLDOUT_PANEL.csv'; cfg=SEARCH_CONFIG['fresh_holdout']
    if not hold.exists():
        rng=np.random.default_rng(cfg['seed']); rows=[]
        while len(rows)<cfg['n']:
            s=rng.uniform(cfg['low'],cfg['high']);q=(s-L.ORIGIN)/L.STEP
            if np.min(abs(q-np.round(q))*L.STEP)<1e-6:continue
            rows.append(dict(panel_id=f"{cfg['prefix']}{len(rows)+1:02d}",**dict(zip(L.AXES,s)),z_true_m=float(rng.choice(cfg['depths']))))
        save(hold.name,pd.DataFrame(rows))
    hd=OUT/'FRESH_HOLDOUT_DESIGN_FREEZE.json'; hashes={p.name:L.sha(p) for p in (hold,manifest)}
    if hd.exists():assert json.loads(hd.read_text())['sha256']==hashes
    else:
        from datetime import datetime,timezone
        L.json_write(hd,dict(frozen_utc=datetime.now(timezone.utc).isoformat(),generation_rule=cfg,sha256=hashes,before_observation_generation=True))


def holdout(shard=None):
    freeze_method();model=Acoustic(inputs(),write_control=shard is None);panel=pd.read_csv(OUT/'FRESH_HOLDOUT_PANEL.csv');rows=[];finals=[];starts_all=[];logs_all=[];observations=[]
    prefix='FRESH_HOLDOUT' if shard is None else f'FRESH_HOLDOUT_SHARD_{shard}'
    for pi,truth in panel.iterrows():
        if shard is not None and pi%4!=shard:continue
        for sigma,seed in [(0.,0)]+[(.1,x) for x in SEARCH_CONFIG['fresh_holdout']['nominal_seeds']]:
            cid=f'{truth.panel_id}_s{sigma:g}_seed{seed}';cache=OUT/'holdout_cache'/cid
            obs=L.generate_observation(truth,sigma,seed,pi,model.models);observations.append((cid,obs.bearing_rad,obs.relative_tl))
            if (cache/'result.json').exists():
                row=json.loads((cache/'result.json').read_text());frame=pd.read_csv(cache/'candidates.csv');starts=pd.read_csv(cache/'starts.csv');logs=pd.read_csv(cache/'logs.csv')
            else:
                started=time.perf_counter();cloud,cmin,cutoff,rc2logs,_=prior.rc2(obs.bearing_rad,sigma)
                if sigma==0:
                    j,z=model.exact(cloud,obs.relative_tl);frame=pd.DataFrame(cloud,columns=L.AXES);frame['J_exact']=j;frame['z_star_label_m']=z;frame['bearing_cost']=np.square(prior.angular_residual(cloud,obs.bearing_rad)).sum(axis=1);frame['keep_tau05']=j<=j.min()+.5
                    logs=pd.DataFrame([dict(stage='NOISELESS_RC2',success=True,nfev=len(cloud))]);starts=pd.DataFrame(cloud,columns=L.AXES);cov=dict(total_bank_nodes=len(cloud),feasible_bank_nodes=len(cloud))
                else:
                    states,scores,_,ri,ui,cov=bank(model,obs.bearing_rad,obs.relative_tl,cutoff)
                    frames=[];log_parts=[];start_parts=[]
                    for budget in SEARCH_CONFIG['budgets']:
                        stride=int(budget['r_step_m']/10);mask=((ri%stride==0)|(ri==1500))&((ui%stride==0)|(ui==206))
                        part,log,start=one_budget(model,obs.bearing_rad,obs.relative_tl,sigma,cutoff,states,scores,mask,budget,cid)
                        frames.append(part);log_parts.append(log);start_parts.append(start)
                    frame=cumulative_candidates(pd.concat(frames),[b['name'] for b in SEARCH_CONFIG['budgets']]).query("budget=='B3'")
                    logs=pd.concat(log_parts);starts=pd.concat(start_parts)
                if 'elapsed_seconds' in cov:cov['bank_elapsed_seconds']=cov.pop('elapsed_seconds')
                top=frame.iloc[0];err=L.errors(frame[list(L.AXES)].to_numpy(),truth);row=dict(case_id=cid,panel_id=truth.panel_id,sigma_deg=sigma,seed=seed,budget='B3',exact_J=float(top.J_exact),bearing_cost=float(top.bearing_cost),cutoff=cutoff,cmin=cmin,recovered=bool(top.J_exact<.001 and top.bearing_cost<=cutoff+1e-14),
                    truth_bearing_cost_eval_only=float(np.square(prior.angular_residual(truth[list(L.AXES)].to_numpy(float),obs.bearing_rad)).sum()),z_star_label_m=float(top.z_star_label_m),n_separated_candidates=len(frame),elapsed_seconds=time.perf_counter()-started,**cov,**{a:float(top[a]) for a in L.AXES})
                for k,name in enumerate(('rel_r','abs_theta_deg','rel_v','abs_psi_deg')):row['top1_'+name]=float(err[0,k]);row['survivor_worst_'+name]=float(err[frame.keep_tau05.to_numpy(),k].max())
                cache.mkdir(parents=True,exist_ok=True)
                for name,df in [('candidates.csv',frame),('starts.csv',starts),('logs.csv',logs),('rc2_logs.csv',rc2logs)]:df.to_csv(cache/name,index=False,float_format='%.17g')
                L.json_write(cache/'result.json',row)
            for df in (frame,starts,logs):
                if 'case_id' not in df:df.insert(0,'case_id',cid)
            rows.append(row);finals.append(frame);starts_all.append(starts);logs_all.append(logs)
            save(prefix+'_RESULTS.csv',pd.DataFrame(rows));print('FRESH_HOLDOUT',cid,'J',row['exact_J'],'recovered',row['recovered'],flush=True)
    save(prefix+'_CANDIDATES.csv',pd.concat(finals));save(prefix+'_STARTS.csv',pd.concat(starts_all));save(prefix+'_OPTIMIZATION_LOG.csv',pd.concat(logs_all))
    np.savez_compressed(OUT/(prefix+'_OBSERVATIONS.npz'),case_ids=np.array([x[0] for x in observations]),bearing_rad=np.array([x[1] for x in observations]),relative_tl=np.array([x[2] for x in observations]))


def save(name,df):
    OUT.mkdir(exist_ok=True)
    df.to_csv(OUT/name,index=False,float_format='%.17g',lineterminator='\n')


class Acoustic(prior.ContinuousAcoustic):
    """Same pressure surrogate; constructor writes only this Gate's directory."""
    def __init__(self,models,write_control=True):
        self.models=models; self.splines={}; self.depth_splines={}; rows=[]
        x=np.arange(39000.,66001.,1.)
        points=np.random.default_rng(20261003).uniform(39000,66000,256)
        for freq,m in models.items():
            iz,_=L.effective_depths(m,L.PROFILE); receiver,_=L.effective_depths(m,[200.])
            weights=m['phi'][iz]*m['phi'][receiver[0]]; k=m['k'].real.astype(float); alpha=-m['k'].imag.astype(float); carrier=(k.min()+k.max())/2
            def pressure(rr):
                rr=np.atleast_1d(rr)[None,:]
                factor=np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*(k[:,None]-carrier)*rr-alpha[:,None]*rr-1j*np.pi/4)
                return (weights@factor).T
            values=np.empty((len(x),21),complex)
            for st in range(0,len(x),2048):values[st:st+2048]=pressure(x[st:st+2048])
            spline=CubicSpline(x,values,axis=0,extrapolate=False); self.splines[freq]=spline
            for zi in range(21):self.depth_splines[freq,zi]=CubicSpline.construct_fast(spline.c[:,:,zi].copy(),spline.x,extrapolate=False,axis=0)
            delta=np.abs(20*np.log10(abs(spline(points)))-20*np.log10(abs(pressure(points))))
            rows.append(dict(frequency_hz=freq,max_error_db=float(delta.max()),tolerance_db=.001,pass_check=bool(delta.max()<.001)))
        if write_control:save('PRESSURE_SURROGATE_CONTROL.csv',pd.DataFrame(rows))
        assert all(r['pass_check'] for r in rows)

    def exact(self,states,tl):
        states=np.atleast_2d(states); js=[]; zs=[]
        for start in range(0,len(states),16):
            j,z=super().exact(states[start:start+16],tl); js.extend(j); zs.extend(z)
        return np.array(js),np.array(zs)


def inputs():return {f:L.parse_mod(L.MODES/f'zgrid_f{f}.mod') for f in L.FREQS}


def protect():
    files=list(prior.OUT.rglob('*'))+list(L.OUT.rglob('*'))+[L.ROOT/'r4_a1_offgrid_bearing.py',L.ROOT/'r4_a1_fix_continuous.py']
    files=[f for f in files if f.is_file()]+[L.MODES/f'zgrid_f{freq}.mod' for freq in L.FREQS]+[L.REFERENCE/'SUBSET_CANDIDATE_SCORES_FIXED.csv']
    path=OUT/'FROZEN_INPUT_HASHES.csv'
    if not path.exists():save(path.name,pd.DataFrame([dict(path=f.relative_to(L.ROOT).as_posix(),sha256=L.sha(f),sha256_lf=hashlib.sha256(f.read_bytes().replace(b'\r\n',b'\n')).hexdigest() if f.suffix in ('.py','.csv') else L.sha(f)) for f in sorted(files)]))
    else:
        frame=pd.read_csv(path)
        for row in frame.itertuples():
            f=L.ROOT/row.path
            portable_source=f.suffix=='.py' or f==L.REFERENCE/'SUBSET_CANDIDATE_SCORES_FIXED.csv'
            if L.sha(f)!=row.sha256:
                assert portable_source and hasattr(row,'sha256_lf') and hashlib.sha256(f.read_bytes().replace(b'\r\n',b'\n')).hexdigest()==row.sha256_lf,row.path
        if 'sha256_lf' not in frame:
            frame['sha256_lf']=[hashlib.sha256((L.ROOT/p).read_bytes().replace(b'\r\n',b'\n')).hexdigest() if (L.ROOT/p).suffix in ('.py','.csv') else L.sha(L.ROOT/p) for p in frame.path]
            save(path.name,frame)


def basin_geometry():
    OUT.mkdir(exist_ok=True); protect(); models=inputs()
    config={'purpose':'ORACLE_DIAGNOSTIC_ONLY; never estimator starts','bearing_feasible_parameterization':'r,v,psi with observation-profiled theta plus fixed truth-profile offset; theta direction varies only the offset',
        'thresholds_db':[.001,.01,.1],'finite_difference_steps':[.0001,.0001,.001,.001],
        'coordinates':['r_km','v_mps','psi_deg','theta_offset_deg'],'normalized_scales':[15.,2.,30.,.16]}
    p=OUT/'BASIN_DIAGNOSTIC_CONFIG.json'
    if p.exists():assert json.loads(p.read_text())==config
    else:L.json_write(p,config)
    cases=pd.read_csv(prior.OUT/'CASE_LEVEL_RESULTS.csv'); observations=np.load(prior.OUT/'CASE_OBSERVATIONS.npz')
    panels=pd.concat([pd.read_csv(L.OUT/'OFFGRID_TRUTH_PANEL.csv'),pd.read_csv(prior.OUT/'HOLDOUT_TRUTH_PANEL.csv')]).set_index('panel_id')
    rows=[]; sensitivities=[]; probes=[]
    for ci,row in cases.iterrows():
        truth=panels.loc[row.panel_id]; tv=truth[list(L.AXES)].to_numpy(float)
        if row.sigma_deg==0:
            rows.append(dict(case_id=row.case_id,direction='NOISELESS_FULL_RANK_SINGLETON',threshold_db=0.001,width=0.,width_status='RC2_NUMERICAL_SINGLETON; NOT_ACOUSTIC_BASIN_WIDTH',oracle_diagnostic_only=True))
            continue
        bearing=observations['bearing_rad'][ci]; tl=observations['relative_tl'][ci]; cutoff=row.rc2_cutoff
        base=prior.theta_profile(tv[[0,2,3]],bearing)[0]
        a0=np.r_[tv[[0,2,3]],tv[1]-base[1]]
        def state(a):
            s=prior.theta_profile(a[:3],bearing)[0]; s[1]+=a[3]; return s
        def evaluate(a):
            s=state(a); cost=float(np.square(prior.angular_residual(s,bearing)).sum())
            if cost>cutoff+1e-14 or (s<prior.LOW).any() or (s>prior.HIGH).any():return None,cost,None
            feature=L.direct_features(models,L.geometry(s)[1]); js,z=L.score_features(feature,tl)
            return float(js[0]),cost,float(z[0])
        steps=np.array(config['finite_difference_steps']); center=L.direct_features(models,L.geometry(tv)[1],[truth.z_true_m])[0,0]
        jac=[]
        for k,h in enumerate(steps):
            ap=a0.copy(); am=a0.copy(); ap[k]+=h; am[k]-=h
            # Sensitivity directions are measured only if both probes are RC2-feasible.
            fp,_,_=evaluate(ap); fm,_,_=evaluate(am)
            assert fp is not None and fm is not None,'Chosen diagnostic derivative step violates bearing likelihood'
            plus=L.direct_features(models,L.geometry(state(ap))[1],[truth.z_true_m])[0,0]
            minus=L.direct_features(models,L.geometry(state(am))[1],[truth.z_true_m])[0,0]
            jac.append(((plus-minus)/(2*h)).ravel()/np.sqrt(tl.size))
        jac=np.array(jac).T; scales=np.array(config['normalized_scales']); _,singular,vt=np.linalg.svd(jac*scales,full_matrices=False)
        for k in range(4):sensitivities.append(dict(case_id=row.case_id,sensitivity_rank=k,singular_value_db=singular[k],
            normalized_direction=json.dumps(vt[k].tolist()),coordinate_units=json.dumps(config['coordinates']),
            Hessian_eigenvalue_for_half_J_squared=float(singular[k]**2),oracle_diagnostic_only=True))
        directions=[(name,np.eye(4)[k],steps[k]) for k,name in enumerate(config['coordinates'])]
        directions += [(f'coupled_svd_{k}',vt[k]*scales,1e-6) for k in range(4)]
        for name,direction,step in directions:
            for target in config['thresholds_db']:
                walls=[]; statuses=[]; z_walls=[]
                derivative=np.linalg.norm(jac@direction)
                for sign in (-1,1):
                    h=max(step/100,target/max(derivative,1e-12)/2)
                    prev=0.; value=0.; status='ACOUSTIC_FIRST_CROSSING'
                    for _ in range(32):
                        value,cost,z=evaluate(a0+sign*h*direction)
                        if value is None:
                            # Locate bearing/bounds wall; do not score infeasible landscape points.
                            lo,hi=prev,h
                            for _ in range(24):
                                mid=(lo+hi)/2; v,_,_=evaluate(a0+sign*mid*direction)
                                if v is None:hi=mid
                                else:lo=mid
                            h=lo; status='CENSORED_BY_BEARING_OR_STATE_BOUND'; break
                        if value>=target:
                            h=brentq(lambda x:evaluate(a0+sign*x*direction)[0]-target,prev,h,xtol=1e-12); break
                        prev=h; h*=2
                    walls.append(h); statuses.append(status)
                    val,cost,z=evaluate(a0+sign*h*direction); z_walls.append(z)
                    probes.append(dict(case_id=row.case_id,direction=name,threshold_db=target,sign=sign,offset=h,J=val,bearing_cost=cost,cutoff=cutoff,z_star_label_m=z,feasible=True,oracle_diagnostic_only=True))
                rows.append(dict(case_id=row.case_id,direction=name,threshold_db=target,negative_extent=walls[0],positive_extent=walls[1],width=sum(walls),
                    width_status=';'.join(statuses),z_negative=z_walls[0],z_positive=z_walls[1],nuisance_depth_changes=any(z!=truth.z_true_m for z in z_walls),
                    direction_vector=json.dumps(direction.tolist()),oracle_diagnostic_only=True))
        save('BASIN_GEOMETRY.csv',pd.DataFrame(rows)); save('BASIN_SENSITIVITY.csv',pd.DataFrame(sensitivities)); save('BASIN_FEASIBLE_PROBES.csv',pd.DataFrame(probes))
        print('BASIN_GEOMETRY',row.case_id,'normalized sensitivity',singular.tolist(),flush=True)


def main():
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['geometry','development','pool','independent','freeze','holdout']);p.add_argument('--shard',type=int,choices=range(4));args=p.parse_args()
    assert args.shard is None or args.action=='holdout','Only independent confirmation cases may be sharded'
    with threadpool_limits(limits=2):{'geometry':basin_geometry,'development':development,'pool':close_development_pools,'independent':independent,'freeze':freeze_method,'holdout':lambda:holdout(args.shard)}[args.action]()


if __name__=='__main__':main()
