"""Local reconstruction and leakage checks for the continuous repair Gate."""
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import qmc
from threadpoolctl import threadpool_limits
import r4_a1_fix_continuous as fix


def audit():
    out=fix.OUT; legacy=fix.legacy
    cases=pd.read_csv(out/'CASE_LEVEL_RESULTS.csv')
    candidates=pd.read_csv(out/'CONTINUOUS_ACOUSTIC_SEARCH.csv')
    regions=pd.read_csv(out/'RC2_CONTINUOUS_HYPOTHESES.csv')
    panel=pd.concat([pd.read_csv(legacy.OUT/'OFFGRID_TRUTH_PANEL.csv').assign(group='REGRESSION'),
                     pd.read_csv(out/'HOLDOUT_TRUTH_PANEL.csv').assign(group='HOLDOUT')])
    models={f:legacy.parse_mod(legacy.MODES/f'zgrid_f{f}.mod') for f in legacy.FREQS}
    records=[]; obs_b=[]; obs_tl=[]; convergence=[]; initial=[]; alias=[]; populations=[]; modes=[]
    def check(name,passed,detail=''):
        records.append(dict(check=name,pass_check=bool(passed),detail=detail))
        if not passed:raise AssertionError(name+': '+detail)
    t=np.arange(121)*10.; dt=np.maximum(t-600,0)
    px=2*np.minimum(t,600)+2*dt*np.cos(np.pi/12); py=2*dt*np.sin(np.pi/12)
    def bearing(states):
        r,theta,v,psi=np.atleast_2d(states).T
        x=1000*r[:,None]*np.cos(np.radians(theta))[:,None]+v[:,None]*t*np.cos(np.radians(psi))[:,None]-px
        y=1000*r[:,None]*np.sin(np.radians(theta))[:,None]+v[:,None]*t*np.sin(np.radians(psi))[:,None]-py
        return np.arctan2(y,x)
    def cost(states,observed):
        delta=bearing(states)-observed
        return np.square(np.arctan2(np.sin(delta),np.cos(delta))).sum(axis=1)
    for row in cases.itertuples():
        group=panel[panel.group==row.group].reset_index(drop=True)
        pi=int(group.index[group.panel_id==row.panel_id][0]); truth=group.loc[pi]
        tv=truth[list(legacy.AXES)].to_numpy(float)
        noisy=bearing(tv)[0]
        if row.sigma_deg:
            noisy+=np.radians(row.sigma_deg)*np.random.default_rng(np.random.SeedSequence([int(row.seed),pi])).standard_normal(121)
        observation=legacy.generate_observation(truth,row.sigma_deg,int(row.seed),pi,models)
        check(row.case_id+'_observations',np.allclose(noisy,observation.bearing_rad,rtol=0,atol=1e-14))
        obs_b.append(noisy); obs_tl.append(observation.relative_tl)
        sub=candidates[candidates.case_id==row.case_id]; states=sub[list(legacy.AXES)].to_numpy()
        c=cost(states,noisy)
        check(row.case_id+'_bearing_costs',np.allclose(c,sub.bearing_cost,rtol=0,atol=1e-13))
        check(row.case_id+'_likelihood_feasibility',bool((c<=row.rc2_cutoff+1e-13).all()))
        js=[]; depths=[]
        for start in range(0,len(states),16):
            f=legacy.direct_features(models,legacy.geometry(states[start:start+16])[1])
            # Independent RMS reduction across frequencies/times and shared z.
            diff=f-observation.relative_tl[None,None,:,:]
            rms=np.linalg.norm(diff.reshape(len(f),21,-1),axis=2)/np.sqrt(363)
            iz=rms.argmin(axis=1); js.extend(rms[np.arange(len(f)),iz]); depths.extend(legacy.PROFILE[iz])
        js=np.array(js); depths=np.array(depths)
        check(row.case_id+'_exact_TRIPLE_scores',np.allclose(js,sub.J_exact,rtol=0,atol=1e-7),str(float(np.max(abs(js-sub.J_exact.to_numpy())))))
        check(row.case_id+'_shared_depth_profile',np.array_equal(depths,sub.z_star_label_m))
        top=int(np.argmin(sub.J_exact)); keep=sub.J_exact.to_numpy()<=sub.J_exact.min()+.5
        def errors(s):
            delta=s-tv; e=np.abs(delta)
            for a in (1,3):e[:,a]=np.degrees(np.abs(np.arctan2(np.sin(np.radians(delta[:,a])),np.cos(np.radians(delta[:,a])))))
            e[:,0]/=tv[0]; e[:,2]/=tv[2]
            return e
        e=errors(states); metricdiff=[]
        for a,name in enumerate(('rel_r','abs_theta_deg','rel_v','abs_psi_deg')):
            metricdiff.extend([abs(getattr(row,'top1_'+name)-e[top,a]),abs(getattr(row,'best_candidate_'+name)-e[:,a].min()),abs(getattr(row,'survivor_worst_'+name)-e[keep,a].max())])
        check(row.case_id+'_metrics_reconstructed',max(metricdiff)<1e-9)
        check(row.case_id+'_survivor_counts',int(keep.sum())==row.n_survivors)
        rc=regions[regions.case_id==row.case_id]; rcstates=rc[list(legacy.AXES)].to_numpy()
        check(row.case_id+'_RC2_region',np.allclose(cost(rcstates,noisy),rc.bearing_cost,atol=1e-13,rtol=0) and len(rc)==row.n_rc2_hypotheses)
        if row.sigma_deg:
            for island,seed in enumerate(fix.CONFIG['acoustic_DE_seeds']):
                pop=fix.acoustic_population(rcstates,noisy,row.sigma_deg,row.rc2_cutoff,np.random.default_rng(seed))
                rv=fix.LOW[[0,2,3]]+pop[:,:3]*fix.WIDTH[[0,2,3]]
                init=fix.theta_profile(rv,noisy)
                bound=row.sigma_deg*np.sqrt(13.3/121)/(45000*39000/66000**2)
                init[:,1]+=pop[:,3]*bound
                check(row.case_id+f'_DE_initialization_{island}',np.ptp(pop[:,3])>1e-8 and (cost(init,noisy)<=row.rc2_cutoff+1e-13).all())
                for i,s in enumerate(init):populations.append(dict(case_id=row.case_id,island=island,start_id=i,**dict(zip(legacy.AXES,s))))
        starts=pd.read_csv(out/'case_cache'/row.case_id/'rc2.csv'); ac=pd.read_csv(out/'case_cache'/row.case_id/'ac.csv')
        convergence.append(dict(case_id=row.case_id,RC2_starts=len(starts),RC2_converged=int(starts.success.sum()),RC2_fraction=starts.success.mean(),
            DE_islands=int((ac.stage=='GLOBAL_DE').sum()),DE_budget_reached=int(((ac.stage=='GLOBAL_DE')&~ac.success).sum()),
            DE_objective_batches=int(ac[ac.stage=='GLOBAL_DE'].nfev.sum()),DE_bearing_state_evaluations_upper_bound=int(ac[ac.stage=='GLOBAL_DE'].nfev.sum())*fix.CONFIG['acoustic_DE_population'],
            local_starts=int(ac.stage.str.startswith('LOCAL_').sum()),local_converged=int(ac[ac.stage.str.startswith('LOCAL_')].success.sum()),
            local_convergence_fraction=float(ac[ac.stage.str.startswith('LOCAL_')].success.mean()) if ac.stage.str.startswith('LOCAL_').any() else 1.))
        modes.append(dict(case_id=row.case_id,noiseless_bearing_modes_certified=1 if row.sigma_deg==0 and row.geometry_rank==4 else -1,
            separated_exported_acoustic_candidate_clusters=len(states),joint_exact_aliases_demonstrated=0,
            mode_count_status='UNIQUE_NOISELESS_BEARING_SOLUTION_BY_FULL_RANK' if row.sigma_deg==0 else 'FINITE_HYPOTHESIS_PROXY_NOT_EXHAUSTIVE_MODE_COUNT'))
        # Reconstruct every optimizer start exclusively from its observation.
        A=np.column_stack((-1000*np.sin(noisy),1000*np.cos(noisy),-t*np.sin(noisy),t*np.cos(noisy)))
        cart=np.linalg.lstsq(A,-px*np.sin(noisy)+py*np.cos(noisy),rcond=None)[0]
        linear=np.array([np.hypot(*cart[:2]),np.degrees(np.arctan2(cart[1],cart[0])),np.hypot(*cart[2:]),np.degrees(np.arctan2(cart[3],cart[2]))])
        st=qmc.Sobol(4,scramble=True,seed=770001).random_base2(5); st[0]=np.clip((linear-fix.LOW)/fix.WIDTH,1e-9,1-1e-9)
        for i,s in enumerate(fix.LOW+st*fix.WIDTH):initial.append(dict(case_id=row.case_id,start_id=i,**dict(zip(legacy.AXES,s))))
        # Demonstrating aliases requires two distinct exact observation matches.
        distinct=np.max(np.abs(states-tv)/np.array(fix.CONFIG['candidate_dedup_resolution']),axis=1)>1
        bearing_rms_deg=np.degrees(np.sqrt(cost(states,bearing(tv)[0])/121))
        matched=distinct&(bearing_rms_deg<1e-8)&(js<1e-6)
        alias.append(dict(case_id=row.case_id,n_distinct_joint_exact_matches=int(matched.sum()),bearing_RMS_tolerance_deg=1e-8,acoustic_RMS_tolerance_db=1e-6,
            max_bearing_RMS_deg=float(bearing_rms_deg.max()),status='ALIAS_DEMONSTRATED' if matched.any() else 'NO_EXACT_ALIAS_DEMONSTRATED_IN_EXPORTED_CANDIDATES'))
        modes[-1]['joint_exact_aliases_demonstrated']=int(matched.sum())
    manifest=json.loads((out/'HOLDOUT_DESIGN_FREEZE.json').read_text())
    for name,digest in manifest['sha256'].items():check('freeze_'+name,legacy.sha(out/name)==digest)
    # Physical regression challenge bytes must remain exactly as committed.
    import subprocess
    for name in ('OFFGRID_TRUTH_PANEL.csv','R4_A1_CONFIG.json','CASE_LEVEL_RESULTS.csv','PILOT_OBSERVATIONS.npz','DESIGN_FREEZE_MANIFEST.json'):
        path=(legacy.OUT/name).relative_to(legacy.ROOT).as_posix()
        raw=subprocess.check_output(['git','show','4633ff0:'+path])
        check('legacy_bytes_'+name,hashlib.sha256(raw).hexdigest()==legacy.sha(legacy.OUT/name))
    np.savez_compressed(out/'CASE_OBSERVATIONS.npz',case_ids=cases.case_id.to_numpy(dtype=str),bearing_rad=np.array(obs_b),relative_tl=np.array(obs_tl),times_s=t)
    fix.save('RESULT_RECONSTRUCTION_AUDIT.csv',pd.DataFrame(records)); fix.save('CONVERGENCE_AUDIT.csv',pd.DataFrame(convergence)); fix.save('RC2_INITIAL_STARTS.csv',pd.DataFrame(initial)); fix.save('CONTINUOUS_ALIAS_DIAGNOSTIC.csv',pd.DataFrame(alias))
    fix.save('DE_INITIAL_POPULATION_RECONSTRUCTION.csv',pd.DataFrame(populations))
    fix.save('MODE_INVENTORY.csv',pd.DataFrame(modes))
    print(len(records),'independent local reconstruction checks passed; scientific audit pending')


if __name__=='__main__':
    with threadpool_limits(limits=2):audit()
