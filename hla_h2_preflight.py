"""Pre-freeze fixed mathematical controls; no eight-scene performance runs."""
import time,json
import numpy as np,pandas as pd
from numpy.polynomial.legendre import leggauss
from threadpoolctl import threadpool_limits
import hla_h2_core as c
import hla_h2_independent as independent
def run():
    c.OUT.mkdir(parents=True,exist_ok=True);c.LOCAL.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter();checks=[]
    def check(name,error,tol):
        checks.append(dict(name=name,error=float(error),tolerance=tol,PASS=bool(error<=tol)))
    rng=np.random.default_rng(2026101012)
    states=rng.uniform(c.LOW+.001*c.WIDTH,c.HIGH-.001*c.WIDTH,(16,4))
    states=np.r_[states,[c.LOW+1e-6*c.WIDTH,c.HIGH-1e-6*c.WIDTH]]
    check('independent_cartesian_geometry',max(np.max(abs(a-b)) for a,b in zip(c.geometry(states),independent.geometry(states))),1e-8)
    maxjac=0.
    for h in states:
        beta,q,rho,jb,jq=c.geometry(h,True);analytic=np.r_[jb[0]/c.SIG_B,jq[0]/.05]*c.WIDTH[None,:]
        fd=[]
        for axis in range(4):
            d=np.zeros(4);d[axis]=c.WIDTH[axis]*1e-6
            bp,qp,_=independent.geometry(h+d);bm,qm,_=independent.geometry(h-d)
            fd.append(np.r_[independent.wrap(bp[0]-bm[0])/c.SIG_B,(qp[0]-qm[0])/.05]/2e-6)
        fd=np.array(fd).T
        error=np.linalg.norm(analytic-fd,axis=0)/np.maximum(np.linalg.norm(fd,axis=0),1.)
        maxjac=max(maxjac,float(error.max()))
    check('analytic_normalized_Jacobian_central_difference',maxjac,1e-5)
    ubobs=c.Observation(c.geometry(states[0])[0][0],abs(c.geometry(states[0])[1][0]),.05,'UB-matched')
    point=np.r_[(states[0]-c.LOW)/c.WIDTH,.65]
    evaluator=c.Evaluator(ubobs,-ubobs.radial)
    analytic=evaluator.derivative(point)
    step=np.zeros(5);step[4]=1e-6
    fd=(evaluator.fun(point+step)-evaluator.fun(point-step))/2e-6
    check('UB_normalized_bias_Jacobian',np.linalg.norm(analytic[:,4]-fd)/max(1,np.linalg.norm(fd)),1e-5)
    import hla_h2_execute as execute
    check('fixed_configuration_plan',abs(len(execute.plan())-2816),0)
    check('empty_and_rejected_selection',0 if execute.choose([]) is None and execute.choose([dict(exception='',accepted_profile=False)]) is None else 1,0)
    toy=[dict(exception='',accepted_profile=True,J_beta=1,J_profile=2,terminal_id=i) for i in (2,1)]
    check('accepted_tie_selection',abs(execute.choose(toy)['terminal_id']-1),0)
    import hla_h2_report as report
    check('Wilson_zero_denominator',0 if np.isnan(report.wilson(0,0)).all() else 1,0)
    toyresults=pd.DataFrame([dict(config_id=0,geometry=0,replicate=0,method='B',sigma=0,status='ACCEPTED',true_compatible=True,search_failure=False,solver_exceptions=0,solver_nonconverged=0,**dict(zip(c.ERROR_NAMES,[.1,.2,.3,.4]))),
        dict(config_id=1,geometry=0,replicate=1,method='B',sigma=0,status='ABSTAIN',true_compatible=True,search_failure=True,solver_exceptions=0,solver_nonconverged=0,**dict(zip(c.ERROR_NAMES,[np.inf]*4)))])
    toystats=report.statistics(toyresults).iloc[0]
    check('conditional_unconditional_abstain_separation',0 if toystats.accepted==1 and toystats.output_rate==.5 and np.isinf(toystats.range_relative_unconditional_P95) and toystats.range_relative_conditional_P95==.1 else 1,0)

    # Exact endpoint displacement equals high-order integral of instantaneous rate.
    nodes,weights=leggauss(64);h=np.array([51.123,.7,1.73,-6.123])
    _,q,_=c.geometry(h);integral=[]
    for j in range(6):
        t=200*j+100*(nodes+1)
        th,ps=np.radians(h[[1,3]])
        rel=np.column_stack((1000*h[0]*np.cos(th)+h[2]*t*np.cos(ps),1000*h[0]*np.sin(th)+h[2]*t*np.sin(ps)))-c.platform(t)
        vo=np.tile([2.,0.],(len(t),1));vo[t>600]=2*np.array([np.cos(np.radians(15)),np.sin(np.radians(15))])
        vt=h[2]*np.array([np.cos(ps),np.sin(ps)])
        rate=np.einsum('ij,ij->i',rel,vt-vo)/np.linalg.norm(rel,axis=1)
        integral.append(float(np.dot(weights,rate)/2))
    check('endpoint_increment_vs_Gauss_integrated_signed_rate',np.max(abs(q[0]-integral)),1e-10)
    # Straight co-speed and turn boundary.
    _,co,rho=c.geometry([50.,0.,2.,0.])
    check('co_speed_first_three_windows',np.max(abs(co[0,:3])),0)
    check('turn_remaining_nonzero',0 if np.all(co[0,3:]>0) else 1,0)
    check('wrap_179_minus179_2deg',abs(abs(np.degrees(c.wrap(np.radians(179-(-179)))))-2),1e-10)
    # Toy crossing, not a selected scientific trajectory.
    endpoint_rho=np.hypot(np.array([-100.,100.]),1000.)
    absrate=np.dot(weights,abs(100*nodes)/np.hypot(100*nodes,1000))/2
    check('within_window_sign_change_net_not_mean_absolute',0 if abs(np.diff(endpoint_rho)[0])==0 and absrate>0 else 1,0)
    beta,q,_=c.geometry(h)
    init,meta=c.linear_radial_start(beta[0],q[0])
    check('noiseless_linear_geometry_identity',np.max(abs((init-h)/c.WIDTH)),1e-8)
    initb,metab=c.linear_radial_start(beta[0],q[0]+.1,True)
    check('noiseless_bias_linear_identity',max(np.max(abs((initb[:4]-h)/c.WIDTH)),abs(initb[4]-.1)),1e-7)
    obs=c.Observation(beta[0],q[0],.05,'S')
    terminal,b,fit,calls=c.fit_start(obs,q[0],init)
    check('noiseless_final_self_match',np.square(c.Evaluator(obs,q[0]).fun((terminal-c.LOW)/c.WIDTH)).sum(),1e-12)
    score=c.physical_scores(obs,terminal,0,q[0])
    independent_r=independent.residual(terminal,beta[0],q[0],.05)
    check('direct_residual_independent_score',abs(score['J_beta']+score['J_branch']-np.square(independent_r).sum()),1e-7)
    for i in range(16):
        qtoy=rng.uniform(-.5,.5,6);utoy=abs(rng.normal(0,.3,6))
        sse=np.square(qtoy[None,:]-c.SIGNS*utoy).sum(axis=1).min()
        check(f'folded_min_sign_identity_{i}',abs(sse-np.square(abs(qtoy)-utoy).sum()),1e-12)
        b,cs=c.profile_bias(qtoy,utoy);bi,ci=independent.profile(qtoy,utoy)
        check(f'bias_all_sign_vs_piecewise_{i}',abs(cs-ci[0]),1e-12)
        dense=np.linspace(-.2,.2,20001);dc=np.square(abs(qtoy[None,:]+dense[:,None])-utoy).sum(axis=1).min()
        check(f'bias_dense_floor_{i}',max(cs-dc,0.),1e-12)
    # Noisy observation sign may disagree with underlying q; not filtered.
    trueq=np.array([.001]*6);noise=np.array([-.002]*6);observed=trueq+noise
    check('noise_sign_flip_preserved',0 if np.all(np.sign(observed)!=np.sign(trueq)) and len(c.SIGNS)==64 else 1,0)
    check('nearest_rank_32_P95_is_31st',abs(c.nearest_rank(np.arange(1,33),.95)-31),0)
    check('INF_quantile_retained',0 if np.isinf(c.nearest_rank([1,2,np.inf,np.inf],.95)) else 1,0)
    check('empty_accepted_quantile_undefined',0 if np.isnan(c.nearest_rank([],.5)) else 1,0)
    terminal=dict(config_id=0,exception=np.nan,accepted_profile=True,accepted_original=True,profile_b=0.,b_original=0.,J_beta=0.,J_profile=0.,J_folded=0.,terminal_id=0,branch_id=0,**dict(zip(c.AXES,h)))
    duplicate=dict(terminal,terminal_id=1,branch_id=1)
    envelope=report.terminal_envelopes(pd.DataFrame([terminal,duplicate]),toyresults.iloc[:1],dict(horizontal_truth=np.array([h]))).iloc[0]
    check('duplicate_terminal_cluster_retains_sources',0 if envelope.clusters==1 and envelope.accepted_records==2 and envelope.distinct_expansion_branches==2 else 1,0)
    # Fixed hypothetical numerical noise control for accounting, not panel tuning.
    dummy=c.Observation(beta[0]+np.radians(.1)*rng.normal(size=121),q[0]+.05*rng.normal(size=6),.05,'S')
    cost_start=time.perf_counter()
    for signs in c.SIGNS[:4]:
        un=signs*abs(dummy.radial);initial,details=c.linear_radial_start(dummy.beta,un)
        for start in (initial,c.CENTER):c.fit_start(c.Observation(dummy.beta,abs(dummy.radial),.05,'U'),un,start)
    benchmark=(time.perf_counter()-cost_start)/8
    check('observation_API_no_hidden_generator_fields',0 if set(c.Observation.__dataclass_fields__)=={'beta','radial','sigma','kind'} else 1,0)
    check('single_component_q_omega_gauge',abs((1-.8*.12)-((1+.8*(.32-.12))-.8*.32)),1e-14)
    pd.DataFrame(checks).to_csv(c.OUT/'PREFLIGHT_CHECKS.csv',index=False)
    result=dict(PASS=all(v['PASS'] for v in checks),count=len(checks),max_Jacobian_column_normalized_difference=maxjac,scientific_panel_observations_generated=False,benchmark_fixed_nonpanel_control=dict(fits=8,seconds_per_fit=benchmark,estimated_203008_fit_seconds=benchmark*203008),elapsed_seconds=time.perf_counter()-started,implementation_control_only=True)
    c.write_json(c.OUT/'PREFLIGHT.json',result);print(json.dumps(result,indent=2))
    if not result['PASS']:raise RuntimeError('Preflight controls failed')
if __name__=='__main__':
    with threadpool_limits(limits=4):run()
