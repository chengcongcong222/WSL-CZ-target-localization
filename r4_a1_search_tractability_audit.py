"""Cost-only calibration and independent pre-run evidence reconstruction."""
import argparse
import hashlib
import inspect
import json
import math
import platform
import subprocess
import time
from datetime import datetime,timezone
from pathlib import Path
import importlib.metadata as metadata
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
import r4_a1_search_tractability as t


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_frozen_input(path,raw_sha256,lf_sha256,baseline_git_blob_oid,allow_line_endings):
    data=path.read_bytes()
    if hashlib.sha256(data).hexdigest()==raw_sha256:return 'RAW_BYTES_IDENTICAL'
    normalized=data.replace(b'\r\n',b'\n')
    blob_oid=hashlib.sha1(b'blob '+str(len(normalized)).encode()+b'\0'+normalized).hexdigest()
    assert allow_line_endings and hashlib.sha256(normalized).hexdigest()==lf_sha256 and blob_oid==baseline_git_blob_oid,path
    return 'GIT_TEXT_CRLF_LF_ONLY'


def input_snapshot_row(p,blob_oid,eol):
    data=(t.ROOT/p).read_bytes();normalized=data.replace(b'\r\n',b'\n')
    normalized_oid=hashlib.sha1(b'blob '+str(len(normalized)).encode()+b'\0'+normalized).hexdigest()
    allow=eol.startswith('i/lf') and 'attr/-text' not in eol and normalized_oid==blob_oid
    return dict(path=p,sha256=hashlib.sha256(data).hexdigest(),canonical_lf_sha256=hashlib.sha256(normalized).hexdigest(),baseline_git_blob_oid=blob_oid,allow_line_endings=allow)


def frozen_inputs():
    path=t.OUT/'FROZEN_INPUT_HASHES.csv'
    if path.exists():
        modes=[]
        for row in pd.read_csv(path).itertuples():
            modes.append(verify_frozen_input(t.ROOT/row.path,row.sha256,row.canonical_lf_sha256,row.baseline_git_blob_oid,row.allow_line_endings))
        return modes
    listing=subprocess.check_output(['git','ls-tree','-r',t.BASELINE],text=True).splitlines()
    blobs={line.split('\t',1)[1]:line.split('\t',1)[0].split()[-1] for line in listing};files=list(blobs)
    eols={line.split('\t',1)[1]:line.split('\t',1)[0] for line in subprocess.check_output(['git','ls-files','--eol'],text=True).splitlines()}
    protected=[p for p in files if p.startswith(('results/R3','results/R4_A1_OFFGRID_BEARING_BOUNDARY/','results/R4_A1_FIX_CONTINUOUS_SEARCH/','results/R4_A1_FIX2_ACOUSTIC_COVERAGE/')) or (p.lower().startswith(('r3','r4_a1')) and p.endswith('.py'))]
    t.save(path.name,pd.DataFrame([input_snapshot_row(p,blobs[p],eols.get(p,'')) for p in protected]))
    return ['RAW_BYTES_IDENTICAL']*len(protected)


def fixture():
    models=t.inherited.inputs()
    bearing=t.L.geometry([51.37,.32,2.11,3.42])[0][0]
    return models,bearing,np.zeros((3,121)),13.3*np.radians(.1)**2


class APIQuadratic:
    """Synthetic API-only callback; never claims a direct-modal evaluation."""
    def prepare(self,q):return t.Prepared(np.r_[q,0.,0.],np.zeros((1,121)),0.,False)
    def key(self,p):return (200.,p.state.astype('<f8').tobytes())
    def evaluate(self,p):return t.Evaluation(float(np.sum((p.state[:2]-.37)**2)),None,p.state,0.,200.,False,0)


def calibration():
    frozen_inputs()
    if (t.OUT/'TRACTABILITY_BUDGET_FREEZE.json').exists():raise RuntimeError('Calibration cannot replace frozen budget evidence')
    models,bearing,tl,cut=fixture();rows=[];cost_counters=[]
    points=np.array([[.31,.62],[.42,.53],[.57,.66],[.72,.44]])
    for z in t.L.PROFILE:
        objective=t.ExactObjective(models,bearing,tl,cut,float(z))
        wrapper=t.BudgetObjective(objective,len(points),cache_enabled=False)
        for q in points:
            start=time.perf_counter();wrapper(q);elapsed=time.perf_counter()-start
            row=wrapper.ledger[-1].copy();row.update(depth_label_m=float(z),elapsed_seconds=elapsed,fixture='FIXED_BEARING_ZERO_RELATIVE_TL_COST_ONLY')
            rows.append(row)
        cost_counters.append(dict(depth_label_m=float(z),**wrapper.counters()))
    t.save('OBJECTIVE_RUNTIME_CALIBRATION.csv',pd.DataFrame(rows))
    t.save('OBJECTIVE_RUNTIME_COUNTERS.csv',pd.DataFrame(cost_counters))
    records=[]
    for name in ['SHGO','DIRECT']:
        for cap in [1,3,7,20,64,128]:
            controlled=t.BudgetObjective(APIQuadratic(),cap)
            options={'n':256,'iters':1} if name=='SHGO' else {'maxfun':4096,'maxiter':100}
            _,record=t.run_solver(name,controlled,options)
            assert record['hard_budget_enforcement_supported'] and record['n_objective_requests']==cap
            record['fixture']='API_QUADRATIC_NO_ACOUSTIC_RECOVERY';records.append(record)
            t.save(f'API_REQUEST_LOGS/{name}_CAP_{cap}.csv',pd.DataFrame(controlled.ledger))
    t.save('HARD_BUDGET_TEST_RESULTS.csv',pd.DataFrame(records))
    smokes=[]
    for name in ['SHGO','DIRECT']:
        pair=[]
        for repeat in range(2):
            controlled=t.BudgetObjective(t.ExactObjective(models,bearing,tl,cut,200.),32)
            options={'n':16,'iters':1} if name=='SHGO' else {'maxfun':256,'maxiter':100}
            _,record=t.run_solver(name,controlled,options);catalog=t.candidate_catalog(controlled)
            assert len(catalog)>0
            verified=t.verify_catalog(models,bearing,tl,cut,catalog)
            record.update(repeat=repeat+1,fixture='FIXED_BEARING_ZERO_RELATIVE_TL_COMPATIBILITY_ONLY',n_validation_exact_forward_evaluations=len(verified))
            smokes.append(record);pair.append(catalog)
            t.save(f'API_REQUEST_LOGS/{name}_EXACT_REPEAT_{repeat+1}.csv',pd.DataFrame(controlled.ledger))
            t.save(f'API_REQUEST_LOGS/{name}_EXACT_CANDIDATES_{repeat+1}.csv',catalog)
        pd.testing.assert_frame_equal(pair[0],pair[1],check_exact=True)
    t.save('EXACT_API_SMOKE_RESULTS.csv',pd.DataFrame(smokes))
    environment=dict(recorded_utc=datetime.now(timezone.utc).isoformat(),python=platform.python_version(),platform=platform.platform(),
        packages={p:metadata.version(p) for p in ['numpy','scipy','pandas','pytest','threadpoolctl']},blas_threads=2,
        shgo_signature=str(inspect.signature(t.shgo)),direct_signature=str(inspect.signature(t.direct)),
        direct_exception_semantics='SystemError with ObjectiveBudgetExceeded as cause; first blocked callback ends run; never counted as success',
        noisy_development_runs=0)
    t.json_write('ENVIRONMENT_AND_API_CALIBRATION.json',environment)
    print('Cost/API calibration complete; zero noisy development cases')


def freeze():
    frozen_inputs()
    path=t.OUT/'TRACTABILITY_BUDGET_FREEZE.json'
    if path.exists():raise RuntimeError('Existing budget freeze cannot be overwritten')
    data=pd.read_csv(t.OUT/'OBJECTIVE_RUNTIME_CALIBRATION.csv');exact=data[data.exact_forward_calls==1]
    assert len(exact)>0
    p95=float(exact.elapsed_seconds.quantile(.95));median=float(exact.elapsed_seconds.median())
    target=6*3600;factor=2.;weights=(1,4,16);branches=21;cases=21;families=2
    affordable=target/(factor*p95*branches*cases*families*sum(weights))
    base=min(128,2**int(math.floor(math.log2(affordable))))
    if base<16:raise RuntimeError('Calibration cost exceeds declared feasible budget policy')
    caps=[int(base*x) for x in weights]
    shgo={};direct={}
    for i,(level,cap) in enumerate(zip(['T1','T2','T3'],caps)):
        shgo[level]=dict(hard_budget_limit=cap,n=max(8,base//4)*2**i,iters=i+1,sampling_method='simplicial',
            options={'maxev':cap,'minimize_every_iter':True,'local_iter':1},
            minimizer_kwargs={'method':'SLSQP','options':{'maxiter':64,'ftol':1e-12,'eps':1e-8}})
        direct[level]=dict(hard_budget_limit=cap,options={'maxfun':cap,'maxiter':1000,'eps':1e-4,'locally_biased':True,'vol_tol':1e-16,'len_tol':1e-6})
    budget=dict(baseline_commit=t.BASELINE,frozen_utc=datetime.now(timezone.utc).isoformat(),stage='PRE_RUN_CHECKPOINT_ONLY',
        development_authorized=False,noisy_development_runs=0,recovery_threshold_db=.001,shgo=shgo,direct=direct,
        hard_budget_definition='n_objective_requests admitted before cache/forward; N+1 is logged as one refused attempt and raises',
        hard_budget_scope='per solver, case, raw budget and inherited depth branch; 21 branches uniformly enumerated',
        objective_domain={'normalized_bounds':[[0.,1.],[0.,1.]],'range_m':list(t.R_BOUNDS),'radial_velocity_mps':list(t.U_BOUNDS)},
        dedup_tolerances=list(t.DEDUP),state_agreement_tolerances=list(t.STATE_AGREEMENT),depth_agreement_rule='identical inherited profile label',
        exact_alias_tolerances={'bearing_rms_deg':1e-8,'acoustic_rms_db':1e-6},nuisance_depth_labels_m=t.L.PROFILE.tolist(),
        internal_seeds=None,independence='Fresh objective cache and optimizer state per raw budget/branch; no lower-budget or other-solver candidates',
        callback_exception_policy='preserve native exception class and causal chain; accept causal budget exhaustion only with admitted count N and exactly one refused callback; unrelated SystemError re-raised',
        budget_choice_basis='fixed algorithmic doubling sampling/quadrupling request hierarchy; cost/API calibration only; no recovery result read',
        cost_policy={'target_total_seconds':target,'exact_and_validation_cost_factor':factor,'request_weights':list(weights),'base_upper_bound':128,
            'base_rule':'min(128, largest power of two <= target/(2*p95*21*21*2*21)); minimum feasible base=16',
            'median_exact_request_seconds':median,'p95_exact_request_seconds':p95,'computed_affordable_base':affordable,'selected_base':base,
            'maximum_admitted_development_requests':int(sum(caps)*branches*cases*families),
            'estimated_exact_plus_validation_seconds':float(sum(caps)*branches*cases*families*p95*factor),
            'runtime_estimate_limit':'Estimate from callback cost, not a wall-clock guarantee; excludes catalog/I/O overhead'},
        prerequisite_evidence_sha256={name:sha(t.OUT/name) for name in ['OBJECTIVE_RUNTIME_CALIBRATION.csv','OBJECTIVE_RUNTIME_COUNTERS.csv','HARD_BUDGET_TEST_RESULTS.csv','EXACT_API_SMOKE_RESULTS.csv','ENVIRONMENT_AND_API_CALIBRATION.json']})
    t.json_write(path.name,budget)
    paths=[t.ROOT/'r4_a1_search_tractability.py',Path(__file__),t.ROOT/'tests/test_r4_a1_search_tractability.py']
    paths += [t.ROOT/'.gitattributes',t.OUT/'.gitattributes']
    paths += [t.OUT/name for name in ['A1_SEARCH_TRACTABILITY_AUDIT_DESIGN.md','TRACTABILITY_METHOD.md','TRACTABILITY_BUDGET_FREEZE.json','FROZEN_INPUT_HASHES.csv']]
    paths += list((t.OUT/'API_REQUEST_LOGS').glob('*.csv'))
    t.json_write('METHOD_FREEZE.json',dict(frozen_utc=datetime.now(timezone.utc).isoformat(),checkpoint='PRE_RUN_ONLY_REQUIRES_RESEARCH_LEAD_RELEASE',noisy_development_runs=0,
        sha256={p.relative_to(t.ROOT).as_posix():sha(p) for p in paths},development_enabled_in_code=False))
    print('Frozen per-depth hard request caps',caps,'estimated exact+verification hours',round(budget['cost_policy']['estimated_exact_plus_validation_seconds']/3600,3))


def independent_geometry(s):
    r,theta,v,psi=np.asarray(s,float);theta,psi=np.radians(theta),np.radians(psi)
    times=np.arange(121)*10.;post=np.maximum(times-600,0.)
    x=1000*r*np.cos(theta)+v*times*np.cos(psi)-2*np.minimum(times,600)-2*post*np.cos(np.pi/12)
    y=1000*r*np.sin(theta)+v*times*np.sin(psi)-2*post*np.sin(np.pi/12)
    return np.arctan2(y,x),np.sqrt(x*x+y*y)


def audit():
    input_modes=frozen_inputs();checks=[]
    def check(name,ok,detail):
        checks.append(dict(check=name,pass_check=bool(ok),detail=str(detail)))
        assert ok,(name,detail)
    check('Protected baseline identities',len(input_modes)==len(pd.read_csv(t.OUT/'FROZEN_INPUT_HASHES.csv')),f'{len(input_modes)} files; {input_modes.count("GIT_TEXT_CRLF_LF_ONLY")} Git-declared text checkout conversions; binary/-text files require raw bytes')
    manifest=json.loads((t.OUT/'METHOD_FREEZE.json').read_text())
    check('Frozen code/design/tests/budget bytes',all(sha(t.ROOT/p)==h for p,h in manifest['sha256'].items()),'before any noisy development')
    budget=json.loads((t.OUT/'TRACTABILITY_BUDGET_FREEZE.json').read_text())
    for name in ['shgo','direct']:
        caps=[budget[name][level]['hard_budget_limit'] for level in ['T1','T2','T3']]
        check('Increasing hard caps '+name,0<caps[0]<caps[1]<caps[2],caps)
    rows=pd.read_csv(t.OUT/'OBJECTIVE_RUNTIME_CALIBRATION.csv');p95=float(rows[rows.exact_forward_calls==1].elapsed_seconds.quantile(.95))
    policy=budget['cost_policy'];affordable=policy['target_total_seconds']/(2*p95*21*21*2*21)
    check('Cost-only deterministic budget rule',policy['selected_base']==min(128,2**int(math.floor(math.log2(affordable)))),'No development/recovery input')
    check('Calibration hashes',all(sha(t.OUT/p)==h for p,h in budget['prerequisite_evidence_sha256'].items()),len(budget['prerequisite_evidence_sha256']))
    for path in sorted((t.OUT/'API_REQUEST_LOGS').glob('*.csv')):
        if 'CANDIDATES' in path.name:continue
        frame=pd.read_csv(path);accepted=frame[frame.admitted]
        check('Request sequence '+path.name,list(accepted.request)==list(range(1,len(accepted)+1)),len(accepted))
        check('Refused attempt '+path.name,len(frame[~frame.admitted])==1 and frame.iloc[-1].admitted==False,'exactly one blocked final callback')
        check('Cache/forward accounting '+path.name,not (accepted.cache_hit & (accepted.exact_forward_calls>0)).any(),'cache hit counted as request, no hidden forward')
    for table in ['HARD_BUDGET_TEST_RESULTS.csv','EXACT_API_SMOKE_RESULTS.csv']:
        for row in pd.read_csv(t.OUT/table).itertuples():
            suffix=f'CAP_{row.hard_budget_limit}' if table.startswith('HARD') else f'EXACT_REPEAT_{row.repeat}'
            label=f'{row.solver}_{suffix}'
            log=pd.read_csv(t.OUT/f'API_REQUEST_LOGS/{label}.csv');accepted=log[log.admitted]
            unique=len(accepted[['depth_label_m','state_sha256']].drop_duplicates())
            check('Hard limit '+label,row.n_objective_requests==row.hard_budget_limit and row.n_blocked_requests==1 and row.n_objective_attempts==row.hard_budget_limit+1 and row.hard_budget_triggered and not row.solver_success, row.exception_chain)
            check('All four ledger counters '+label,int(log.admitted.sum())==row.n_objective_requests and int(log.cache_hit.sum())==row.n_cache_hits and int(log.exact_forward_calls.sum())==row.n_exact_forward_evaluations and unique==row.n_unique_states_evaluated,'request / cache / exact forward / unique physical states')
            check('Exception semantics '+label,json.loads(row.exception_chain)==(['ObjectiveBudgetExceeded'] if row.solver=='SHGO' else ['SystemError','ObjectiveBudgetExceeded']),'native class and preserved cause')
    for row in pd.read_csv(t.OUT/'OBJECTIVE_RUNTIME_COUNTERS.csv').itertuples():
        branch=rows[rows.depth_label_m==row.depth_label_m]
        check('Uncached calibration counters '+str(row.depth_label_m),len(branch)==row.n_objective_requests==4 and int(branch.exact_forward_calls.sum())==row.n_exact_forward_evaluations and len(branch.state_sha256.unique())==row.n_unique_states_evaluated and row.n_cache_hits==0 and not row.hard_budget_triggered,'fixed cost fixture; no recovery metric')
    models,bearing,tl,cut=fixture()
    for path in sorted((t.OUT/'API_REQUEST_LOGS').glob('*CANDIDATES*.csv')):
        frame=pd.read_csv(path,float_precision='round_trip')
        for i,row in enumerate(frame.itertuples()):
            s=np.array([getattr(row,a) for a in t.L.AXES]);b,r=independent_geometry(s)
            feature=t.L.direct_features(models,r,[row.z_star_label_m])[0,0]
            j=float(np.sqrt(np.mean((feature-tl)**2)));cost=float(np.sum(((b-bearing+np.pi)%(2*np.pi)-np.pi)**2))
            check(f'Exact saved score {path.name} {i}',abs(j-row.J_exact)<1e-7 and abs(cost-row.bearing_cost)<1e-12 and cost<=cut+1e-14,'independent geometry and RMS')
        check('Exact ranking '+path.name,np.all(np.diff(frame.J_exact)>=0),'sorted by direct score')
    check('Development not reached',budget['development_authorized']==False and not any((t.OUT/name).exists() for name in ['DEVELOPMENT_SHGO_RESULTS.csv','DEVELOPMENT_DIRECT_RESULTS.csv','DEVELOPMENT_EXECUTION_DECISION.json']) and not list(t.OUT.glob('FRESH_*')),'metadata-only harness freeze is allowed; no scientific results or fresh panels')
    harness_active=(t.OUT/'DEVELOPMENT_EXECUTION_FREEZE.json').exists()
    if harness_active:
        from r4_a1_search_tractability_harness_audit import audit_harness
        audit_harness(check)
    t.save('INTEGRITY_AUDIT.csv',pd.DataFrame(checks))
    t.json_write('EXECUTION_HARNESS_CHECKPOINT_DECISION.json' if harness_active else 'PRE_RUN_CHECKPOINT_DECISION.json',dict(decision='DEVELOPMENT_EXECUTION_HARNESS_FREEZE_READY_FOR_INDEPENDENT_AUDIT' if harness_active else 'PRE_RUN_FREEZE_READY_FOR_INDEPENDENT_AUDIT',baseline_commit=t.BASELINE,
        scientific_tractability_decision='NOT_EVALUATED',development_authorized=False,noisy_development_runs=0,
        hard_budget_enforcement='VALIDATED_FOR_BOTH_INSTALLED_APIS_WITH_DISCLOSED_DIRECT_CAUSAL_WRAPPING',
        integrity_checks_passed=len(checks),R4_progress_percent=0,A2_depth_SSP_P5='UNOPENED'))
    print('Tractability integrity checks passed',len(checks))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['calibrate','freeze','audit']);args=parser.parse_args()
    with threadpool_limits(limits=2):{'calibrate':calibration,'freeze':freeze,'audit':audit}[args.action]()
