"""Predeclared fixture calibration only; no development observation access."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,subprocess,time
from decimal import Decimal,ROUND_CEILING
from pathlib import Path
import numpy as np
from flint import arb
from r4_cz_envelope_representation import feature,score
from r4_cz_envelope_support import Cell,Observation,HardBudget,COUNTERS,DOMAIN_LOW,DOMAIN_HIGH
from r4_cz_envelope_model import ModalModel,precision,geometry,point_geometry,project_balls

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_ROUTE_REDESIGN_GATE_3A_CZ_ENVELOPE_PRE_RUN'
SCRATCH=Path('D:/ProjectStorage/WSL-CZ/gate3a-20261004')
POLICY_SHA='a393c7b45c92c42fd148d9c162faeb5ea35160c7'

def load(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def midpoint(values):return np.asarray([float(x.mid()) for x in values])
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    policy_path=OUT/'CZ_NUMERICAL_CALIBRATION_POLICY.json';policy=load(policy_path)
    # Fail before any fixture unless the policy commit is both local and remote main.
    remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],cwd=ROOT,text=True).split()[0]
    if remote!=POLICY_SHA or subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()!=POLICY_SHA:
        raise RuntimeError('Policy commit must be pushed before calibration')
    if policy_path.read_bytes().replace(b'\r\n',b'\n')!=subprocess.check_output(['git','show',POLICY_SHA+':'+policy_path.relative_to(ROOT).as_posix()],cwd=ROOT).replace(b'\r\n',b'\n'):
        raise RuntimeError('Policy changed')
    for rel,sha in policy['input_sha256'].items():
        if digest(ROOT/rel)!=sha:raise RuntimeError('Input identity changed: '+rel)
    budget=HardBudget(dict.fromkeys(COUNTERS,100000),600000)
    model=ModalModel();rows=[];first=None
    offsets=np.zeros((3,121))
    for f,(a,b) in enumerate(policy['analytic_tests']['offsets_line_window']):offsets[f,:61]=a;offsets[f,61:]=b
    for i,state in enumerate(policy['fixed_states_r_km_theta_deg_v_mps_psi_deg']):
        for z in policy['fixed_calibration_depth_labels_m']:
            start=time.perf_counter();point=Cell(tuple(state),tuple(state))
            budget.charge(n_validation_evaluations=1);base=model.levels_numpy(state,z,budget);fbase=feature(base)
            again=model.levels_numpy(state,z,budget);repeat=score(fbase,feature(again))
            with precision(160):balls160=project_balls(model.levels_arb(point,z,budget));f160=midpoint(balls160)
            with precision(256):balls256=project_balls(model.levels_arb(point,z,budget));f256=midpoint(balls256)
            budget.charge(n_validation_evaluations=4)
            metrics={'repeat_dB':repeat,'numpy_vs_arb256_dB':score(fbase,f256),'arb160_vs_arb256_dB':score(f160,f256),'offset_invariance_dB':score(fbase,feature(base+offsets))}
            delta=max(metrics.values())
            if not math.isfinite(delta) or delta>policy['numerical_allowance']['maximum_allowed_delta_num_dB']:
                raise RuntimeError('Predeclared numerical validity limit exceeded; no tau retuning')
            rows.append({'fixture_id':f'DOMAIN_CORNER_{i+1}_Z{z}','state':json.dumps(state),'depth_label_m':z,**metrics,'elapsed_seconds':time.perf_counter()-start})
            if first is None:first=(state,base,fbase,balls256)
            print(json.dumps({'calibration_fixture_completed':rows[-1]['fixture_id'],'elapsed_seconds':rows[-1]['elapsed_seconds']}),flush=True)
    # Analytic feature fixtures exercise constant, linear and quadratic shape plus fixed offsets.
    for order in policy['analytic_tests']['orders']:
        t=np.linspace(-1,1,121);levels=np.tile(t**order,(3,1));f=feature(levels)
        budget.charge(n_validation_evaluations=2)
        rows.append({'fixture_id':f'ANALYTIC_ORDER_{order}','state':'ANALYTIC','depth_label_m':'NONE','repeat_dB':score(f,feature(levels)),
                     'numpy_vs_arb256_dB':0.,'arb160_vs_arb256_dB':0.,'offset_invariance_dB':score(f,feature(levels+offsets)),'elapsed_seconds':0.})
    deltas=[max(r[k] for k in ['repeat_dB','numpy_vs_arb256_dB','arb160_vs_arb256_dB','offset_invariance_dB']) for r in rows]
    maximum=max(deltas);allow=policy['numerical_allowance']
    quantum=Decimal(str(allow['round_up_quantum_dB']))
    raw=max(Decimal(str(allow['minimum_floor_dB'])),Decimal(allow['safety_factor'])*Decimal.from_float(maximum))
    tau=float((raw/quantum).to_integral_value(rounding=ROUND_CEILING)*quantum)
    state,levels,coefficients,reference=first
    angles,_=point_geometry(state);observation=Observation(angles,coefficients,1e9)
    root=Cell(DOMAIN_LOW,DOMAIN_HIGH)
    eps=policy['non_development_cell_bound_tests']['tiny_width_each_axis']
    tiny=Cell(tuple(state),tuple(min(b,a+eps) for a,b in zip(state,DOMAIN_HIGH)))
    runtimes=[];enclosures=[]
    for label,cell in [('ROOT',root),('TINY',tiny)]:
        start=time.perf_counter();budget.charge(n_validation_evaluations=1)
        cert=model.certificate(cell,observation,budget)
        elapsed=time.perf_counter()-start;runtimes.append(elapsed)
        if not cert.bearing.valid() or cert.feature is None or not cert.feature.valid():raise RuntimeError('Invalid whole-cell bounds')
        enclosures.append({'fixture':label,'cell_low':cell.low,'cell_high':cell.high,'elapsed_seconds':elapsed,
                           'bearing_lower':str(cert.bearing.lower),'bearing_upper':str(cert.bearing.upper),
                           'feature_lower':str(cert.feature.lower),'feature_upper':str(cert.feature.upper),'depth_profiles':len(cert.depth_labels)})
        print(json.dumps({'enclosure_fixture_completed':label,'elapsed_seconds':elapsed}),flush=True)
    costs=max(max(runtimes),1e-3);rt=policy['runtime_calibration']
    caps={label:max(1,min(rt['structural_cell_ceiling'][label],math.floor(rt['safety_time_fraction']*rt[label.lower()+'_runtime_envelope_seconds']/costs))) for label in ['COARSE','FINE']}
    with (OUT/'CZ_NUMERICAL_CALIBRATION.csv').open('w',encoding='utf-8',newline='') as handle:
        w=csv.DictWriter(handle,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write(OUT/'CZ_NUMERICAL_CALIBRATION_SUMMARY.json',{'policy_commit':POLICY_SHA,'policy_raw_sha256':digest(policy_path),'fixture_count':len(rows),'delta_num_dB':maximum,'tau_F_dB':tau,
        'tau_formula':allow['tau_formula'],'claim_scope':allow['claim'],'runtime_cell_cost_seconds':costs,'raw_cell_caps':caps,
        'calibration_accounting':budget.report(),'enclosure_fixtures':enclosures,'depth_mapping':model.mapping,
        'dependency':load(SCRATCH/'dependency-install.json'),'development_numerical_runs':0,'development_observations_loaded':False,'fresh':'NOT_GENERATED','scientific_result':False,'source_sha256':{name:digest(ROOT/name) for name in ['r4_cz_envelope_representation.py','r4_cz_envelope_model.py','r4_cz_envelope_support.py','r4_cz_envelope_prerun_audit.py']}})
    write(SCRATCH/'fixed-fixture-reference.json',{'state':state,'depth_label_m':150,'numpy_levels':levels.tolist(),'numpy_coefficients':coefficients.tolist(),'arb256_coefficients_midpoint':midpoint(reference).tolist()})
    print(json.dumps({'tau_F_dB':tau,'raw_cell_caps':caps,'calibration':'PASS'}),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--calibrate',action='store_true');args=parser.parse_args()
    if not args.calibrate:parser.error('Explicit --calibrate required; no development mode exists')
    main()
