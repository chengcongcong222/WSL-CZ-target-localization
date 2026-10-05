"""Frozen application pre-research dynamic acquisition; no acoustics/depth."""
from pathlib import Path
import argparse,csv,json,hashlib,math,subprocess
import numpy as np
from r4_a1_new_estimator import estimate,wrap
OUT=Path('results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC')
BASE='333ead6461a7574a469495033d48fd89682aa496'
COMMIT_PARENT='0eef9ff7d0f40cf05400707b791db8a80e074fac'
TITLE='R4 A1-new: freeze augmented off-grid dynamic baseline'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):
 def clean(v):
  if isinstance(v,float) and not math.isfinite(v):return 'INF'
  if isinstance(v,dict):return {k:clean(a) for k,a in v.items()}
  if isinstance(v,list):return [clean(a) for a in v]
  return v
 Path(p).write_text(json.dumps(clean(x),indent=2,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,data):
 with Path(p).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def git(*a):return subprocess.check_output(['git','-c','core.quotepath=false',*a],text=True,encoding='utf-8').strip()
def trajectories():
 t=np.arange(121)*10.;post=np.maximum(t-600,0.);main=np.column_stack([2*np.minimum(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)]);return t,main

def scene(truth,anchor,u,z,side):
 t,main=trajectories();c=(2*u[0]-1)*anchor['common_deg'];d=(2*u[1]-1)*anchor['half_diff_deg'];beta=np.deg2rad((2*u[2]-1)*anchor['beta_deg']);offset=anchor['baseline_m']*np.array([-np.sin(beta),side*np.cos(beta)]);true_nodes=np.stack([main,main+offset],axis=1)
 r=float(truth['r0_km'])*1000;theta=np.deg2rad(float(truth['theta0_deg']));v=float(truth['v_mps']);psi=np.deg2rad(float(truth['psi_deg']));state=np.array([r*np.cos(theta),r*np.sin(theta),v*np.cos(psi),v*np.sin(psi)]);target=state[:2]+t[:,None]*state[2:];delta=target[:,None,:]-true_nodes
 measured=np.arctan2(delta[:,:,1],delta[:,:,0])+np.deg2rad(np.array([c-d,c+d]))+np.deg2rad(anchor['sigma_deg'])*z[:,:2]
 estimated=true_nodes+anchor['nav_sigma_m']*z[:,2:].reshape(121,2,2)
 return t,estimated,measured,dict(b_common_deg=c,b_diff_HALF_deg=d,beta_deg=float(np.rad2deg(beta)))

def evaluate(state,truth,main_est0):
 rel=state[:2]-main_est0;r=np.linalg.norm(rel);theta=np.arctan2(rel[1],rel[0]);v=np.linalg.norm(state[2:]);psi=np.arctan2(state[3],state[2])
 return dict(r_hat_km=float(r/1000),theta_hat_deg=float(np.rad2deg(theta)),v_hat_mps=float(v),psi_hat_deg=float(np.rad2deg(psi)),range_error=abs(float(r)-float(truth['r0_km'])*1000)/(float(truth['r0_km'])*1000),bearing_error_deg=float(abs(np.rad2deg(wrap(theta-np.deg2rad(float(truth['theta0_deg'])))))),speed_error=abs(float(v)-float(truth['v_mps']))/float(truth['v_mps']),heading_error_deg=float(abs(np.rad2deg(wrap(psi-np.deg2rad(float(truth['psi_deg'])))))))
METRICS=['range_error','bearing_error_deg','speed_error','heading_error_deg']
def quant(values,q):
 value=float(np.sort(np.asarray(values,float))[math.ceil(len(values)*q)-1]);return value if math.isfinite(value) else 'INF'
def summarize(data):
 out=[]
 for anchor in ['A','B']:
  for case in sorted({x['case_id'] for x in data if x['anchor']==anchor}):
   g=[x for x in data if x['anchor']==anchor and x['case_id']==case];row=dict(anchor=anchor,case_id=case,n_runs=len(g),failure_rate=sum(x['failed'] for x in g)/len(g))
   for m in METRICS:
    for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:row[m+'_'+tag]=quant([x[m] for x in g],q)
   row['STRONG_pass']=all(float(row[m+'_P95'])<=v for m,v in zip(METRICS,[.05,.5,.05,2]));row['PROJECT_pass']=all(float(row[m+'_P95'])<=v for m,v in zip(METRICS,[.1,1.,.1,5.]));out.append(row)
 return out

def execute():
 p=load(OUT/'A1_NEW_DESIGN_FREEZE.json')
 for name,h in {**p['bindings'],**p['historical_bindings']}.items():assert sha(name)==h,name
 head=git('rev-parse','HEAD');assert head==git('ls-remote','origin','refs/heads/main').split()[0] and git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==COMMIT_PARENT
 assert not (OUT/'A1_NEW_EXECUTION_START.json').exists(),'one primary execution only'
 dump(OUT/'A1_NEW_EXECUTION_START.json',dict(design_SHA=head,remote_verified_before_MC=True,policy_sha256=sha(OUT/'A1_NEW_DESIGN_FREEZE.json')))
 truth_panel=rows(OUT/'A1_NEW_TRUTH_PANEL.csv');data=[];init=[];diagnostics=[]
 for aname in ['A','B']:
  a=p['anchors'][aname];rng=np.random.Generator(np.random.PCG64(a['seed']));u=rng.random((6000,3));z=rng.standard_normal((6000,121,6));np.savez_compressed(OUT/f'A1_NEW_DRAWS_{aname}.npz',uniforms=u,normals=z)
  group=[]
  for ix in range(6000):
   truth=truth_panel[ix//500];side=-1 if ix%2==0 else 1;t,positions,bearings,errors=scene(truth,a,u[ix],z[ix],side);key=dict(anchor=aname,case_id=truth['case_id'],realization=ix%500,draw_index=ix,mirror=side)
   diag=dict(valid_intersections=0,invalid_intersections=121,min_abs_det=0.,initializer_weight_min=0.,initial_x_m='NA',initial_y_m='NA',initial_vx_mps='NA',initial_vy_mps='NA',solver_success=False,solver_status=0,nfev=0,cost='NA',residual_RMS_deg='NA',jacobian_condition_scaled='NA',approx_se_x_m='NA',approx_se_y_m='NA',approx_se_vx_mps='NA',approx_se_vy_mps='NA',covariance_scope='UNAVAILABLE_FAILED')
   try:
    state,diag=estimate(t,positions,bearings,a['sigma_deg']);metrics=evaluate(state,truth,positions[0,0]);failed=False;reason='NONE'
   except (ValueError,np.linalg.LinAlgError,FloatingPointError) as ex:
    state=np.full(4,np.nan);metrics={k:'INF' for k in METRICS};metrics.update(r_hat_km='NA',theta_hat_deg='NA',v_hat_mps='NA',psi_hat_deg='NA');failed=True;reason=str(ex)
   row=dict(**key,**errors,failed=failed,failure_reason=reason,x_hat_m=float(state[0]) if not failed else 'NA',y_hat_m=float(state[1]) if not failed else 'NA',vx_hat_mps=float(state[2]) if not failed else 'NA',vy_hat_mps=float(state[3]) if not failed else 'NA',**metrics);group.append(row);init.append(dict(**key,**{k:diag[k] for k in ['valid_intersections','invalid_intersections','min_abs_det','initializer_weight_min','initial_x_m','initial_y_m','initial_vx_mps','initial_vy_mps']}));diagnostics.append(dict(**key,**diag))
   if (ix+1)%500==0:print(f'Anchor {aname} {ix+1}/6000',flush=True)
  write(OUT/f'A1_NEW_RUN_RESULTS_{aname}.csv',group)
  if aname=='A':
   write(OUT/'A1_NEW_INITIALIZER_DIAGNOSTICS_A.csv',init);write(OUT/'A1_NEW_ESTIMATOR_DIAGNOSTICS_A.csv',diagnostics);paths=[OUT/'A1_NEW_DRAWS_A.npz',OUT/'A1_NEW_RUN_RESULTS_A.csv',OUT/'A1_NEW_INITIALIZER_DIAGNOSTICS_A.csv',OUT/'A1_NEW_ESTIMATOR_DIAGNOSTICS_A.csv'];dump(OUT/'A1_NEW_PRIMARY_SEAL.json',dict(bindings={x.as_posix():sha(x) for x in paths},before_secondary=True))
  data+=group
 write(OUT/'A1_NEW_RUN_RESULTS.csv',data);write(OUT/'A1_NEW_INITIALIZER_DIAGNOSTICS.csv',init);write(OUT/'A1_NEW_ESTIMATOR_DIAGNOSTICS.csv',diagnostics);summary=summarize(data);write(OUT/'A1_NEW_CASE_SUMMARY.csv',summary)
 global_rows=[]
 for aname in ['A','B']:
  g=[x for x in data if x['anchor']==aname];s=[x for x in summary if x['anchor']==aname];row=dict(anchor=aname,n_runs=len(g),failure_rate=sum(x['failed'] for x in g)/len(g),STRONG_pass=all(x['STRONG_pass'] for x in s),PROJECT_pass=all(x['PROJECT_pass'] for x in s))
  for m in METRICS:
   row[m+'_worst_case_P95']=max(float(x[m+'_P95']) for x in s)
   for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:row[m+'_pooled_'+tag]=quant([x[m] for x in g],q)
  global_rows.append(row)
 write(OUT/'A1_NEW_GLOBAL_SUMMARY.csv',global_rows);primary=global_rows[0];outcome='R4_A1_NEW_STRONG_BASELINE_ESTABLISHED' if primary['STRONG_pass'] else 'R4_A1_NEW_PROJECT_BASELINE_ESTABLISHED' if primary['PROJECT_pass'] else 'R4_A1_NEW_BASELINE_NOT_ESTABLISHED';credit=15 if primary['PROJECT_pass'] else 0
 dump(OUT/'A1_NEW_DECISION.json',dict(scientific_parent_SHA=BASE,commit_parent_SHA=COMMIT_PARENT,design_SHA=head,project_stage='APPLICATION_PRE_RESEARCH',primary_application_scenario='Anchor A',primary_truth_cases=12,primary_runs=6000,secondary_runs=6000,primary_summary=primary,secondary_summary=global_rows[1],scientific_decision=outcome,R4_A1_progress_percent=credit,R4_overall_percent=credit,progress_scope='PRE_REGISTERED_SCIENTIFIC_COMPLETION; PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT',hardware_status='UNKNOWN; NONBLOCKING_DESIGN_ASSUMPTION',depth='NOT_OPENED',next='R4-A2' if credit else 'STOP',automatic_next_experiment=False,stop_after_commit_B=True,new_acoustic_propagation=0,new_depth_score=0))
 print(outcome,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if args.execute:execute()
 else:ap.error('--execute required')
