"""Gate1: frozen cross-track geometry only, no acoustic evaluations."""
from pathlib import Path
import argparse,csv,hashlib,json,math,subprocess
import numpy as np
from r4_aux_gate0_geometry import triangulate
OUT=Path('results/R4_AUX_GATE1_CROSSTRACK_REQUIREMENT')
PARENT='79bc6b5c04e08d6074719e3c045947c1325253c8'
TITLE='R4 AUX Gate1: freeze cross-track requirement boundary'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read_csv(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write_csv(p,rows):
 with Path(p).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def git(*a):return subprocess.run(['git',*a],check=True,capture_output=True,text=True,encoding='utf-8').stdout.strip()
def quantiles(values):
 v=np.sort(values);n=len(v)
 return {tag:float(v[math.ceil(n*q)-1]) for tag,q in [('median',.5),('p90',.9),('p95',.95),('p99',.99)]}
def number(v):return 'INF' if not math.isfinite(v) else v

def evaluate(row,z,policy):
 r,b,s=[float(row[k]) for k in ['range_km','baseline_km','sigma_deg']];side=int(row['side']);p2=np.array([0.,side*b]);alpha=math.atan2(b,r)
 noise=math.radians(s)*side*z
 xy,failed,parallel,behind,numeric,det=triangulate(p2,noise[:,0],-side*alpha+noise[:,1],policy['parallel_det_tolerance'])
 er=np.abs(np.linalg.norm(xy,axis=1)-r)/r;pos=np.linalg.norm(xy-np.array([r,0.]),axis=1);cross=np.abs(xy[:,1])
 for x in [er,pos,cross]:x[failed]=np.inf
 result=dict(row);result.update(alpha_deg=math.degrees(alpha),aux_true_bearing_deg=-side*math.degrees(alpha),n_samples=len(z),execution_valid=True,
  failure_rate=float(failed.mean()),behind_sensor_rate=float(behind.mean()),parallel_rate=float(parallel.mean()),nonfinite_rate=float(numeric.mean()),ill_conditioned_rate=float((1/np.maximum(np.abs(det),np.finfo(float).tiny)>policy['ill_condition_proxy_limit']).mean()))
 for name,arr in [('relative_range',er),('range_abs_km',er*r),('position_2d_km',pos),('cross_range_abs_km',cross)]:
  for tag,v in quantiles(arr).items():result[name+'_'+tag]=number(v)
 # Independent first-order matrix covariance. Truth used only for geometry and scoring.
 q=np.array([r,-side*b]);H=np.array([[0.,1/r],np.array([-q[1],q[0]])/(q@q)]);F=H.T@H/math.radians(s)**2;C=np.linalg.inv(F)
 sig=math.sqrt(C[0,0]);linear=1.959963984540054*sig/r
 result.update(jacobian_condition=float(np.linalg.cond(H)),fim_condition=float(np.linalg.cond(F)),sigma_r_linear_km=sig,sigma_cross_linear_km=math.sqrt(C[1,1]),p95_linear_relative=linear,
  nonlinear_to_linear_p95_ratio=float(result['relative_range_p95'])/linear,nonlinear_deviation_relative=float(result['relative_range_p95'])-linear,
  sigma_rad_over_sin_alpha=math.radians(s)/math.sin(alpha),p95_over_sigma_rad_over_sin_alpha=float(result['relative_range_p95'])/(math.radians(s)/math.sin(alpha)),
  fraction_le5pct=float((er<=.05).mean()),fraction_le10pct=float((er<=.1).mean()))
 n=len(z);zz=1.959963984540054;ordered=np.sort(er);delta=zz*math.sqrt(n*.95*.05)
 result['p95_mc_rank_ci_low']=number(float(ordered[max(0,math.floor(n*.95-delta)-1)]));result['p95_mc_rank_ci_high']=number(float(ordered[min(n-1,math.ceil(n*.95+delta)-1)]))
 for tag,threshold in [('5pct',.05),('10pct',.1)]:
  prop=float((er<=threshold).mean());den=1+zz*zz/n;c=(prop+zz*zz/(2*n))/den;h=zz*math.sqrt(prop*(1-prop)/n+zz*zz/(4*n*n))/den
  result['fraction_'+tag+'_wilson_low']=c-h;result['fraction_'+tag+'_wilson_high']=c+h
 return result

def summarize(results,p):
 summaries=[]
 for s in p['sigma_deg']:
  for b in p['baseline_km']:
   cells=[x for x in results if float(x['baseline_km'])==b and float(x['sigma_deg'])==s]
   finite=all(math.isfinite(float(x['relative_range_p95'])) for x in cells)
   valid=all(x['execution_valid'] is True for x in cells) and len(cells)==2*len(p['range_km']) and sorted({float(x['range_km']) for x in cells})==p['range_km']
   worst=max(cells,key=lambda x:float(x['relative_range_p95']));value=float(worst['relative_range_p95'])
   summaries.append(dict(sigma_deg=s,baseline_km=b,n_ranges=len({x['range_km'] for x in cells}),n_mirror_cells=len(cells),all_ranges_execution_valid=valid,all_p95_finite=finite,
    worst_p95_relative=number(value),worst_range_km=worst['range_km'],worst_side=worst['side'],full_range_5pct_pass=valid and finite and value<=.05,full_range_10pct_pass=valid and finite and value<=.1,
    alpha_min_deg=min(float(x['alpha_deg']) for x in cells),alpha_max_deg=max(float(x['alpha_deg']) for x in cells),max_failure_rate=max(float(x['failure_rate']) for x in cells),
    max_nonlinear_to_linear_ratio=max(float(x['nonlinear_to_linear_p95_ratio']) for x in cells),min_nonlinear_to_linear_ratio=min(float(x['nonlinear_to_linear_p95_ratio']) for x in cells)))
 return summaries

def frontiers(summaries,p):
 frontier=[];fixed=[]
 for s in p['sigma_deg']:
  row=dict(sigma_deg=s)
  for tag in ['5pct','10pct']:
   good=[x for x in summaries if x['sigma_deg']==s and x['full_range_'+tag+'_pass']]
   best=min(good,key=lambda x:x['baseline_km']) if good else None
   row['minimum_tested_B_for_'+tag+'_km']=best['baseline_km'] if best else 'NOT_REACHED_WITHIN_B<=10KM'
   row['worst_range_at_'+tag+'_boundary']=best['worst_range_km'] if best else ''
   row['worst_p95_at_'+tag+'_boundary']=best['worst_p95_relative'] if best else ''
  frontier.append(row)
 for b in p['fixed_baseline_km']:
  row=dict(baseline_km=b)
  for tag in ['5pct','10pct']:
   cells=sorted([x for x in summaries if x['baseline_km']==b],key=lambda x:x['sigma_deg']);good=[x for x in cells if x['full_range_'+tag+'_pass']]
   best=max(good,key=lambda x:x['sigma_deg']) if good else None
   above=[x for x in cells if best and x['sigma_deg']>best['sigma_deg'] and not x['full_range_'+tag+'_pass']]
   row['largest_tested_sigma_for_'+tag+'_deg']=best['sigma_deg'] if best else 'NOT_REACHED_IN_TESTED_SIGMA_GRID'
   row['next_tested_failing_sigma_for_'+tag+'_deg']=above[0]['sigma_deg'] if above else 'NOT_BRACKETED_ABOVE_TESTED_MAX' if best else ''
   row['worst_p95_at_'+tag+'_boundary']=best['worst_p95_relative'] if best else ''
  fixed.append(row)
 return frontier,fixed

def decide(summaries,frontier,fixed,p):
 a=any(x['sigma_deg']==.1 and x['baseline_km']<=5 and x['full_range_5pct_pass'] for x in summaries)
 b=any(x['sigma_deg']>=.05 and x['baseline_km']<=10 and x['full_range_5pct_pass'] for x in summaries)
 outcome='AUX1_CROSSTRACK_GEOMETRY_PRACTICAL_AT_ORIGINAL_SCALE' if a else 'AUX1_BASELINE_BEARING_TRADEOFF_ESTABLISHED' if b else 'AUX1_FULL_RANGE_5PCT_NOT_ESTABLISHED'
 primary=next(x for x in frontier if x['sigma_deg']==.1);selected=[x for x in summaries if x['sigma_deg']==.1 and x['baseline_km']==primary['minimum_tested_B_for_5pct_km']]
 return dict(stage=p['stage'],parent_sha=PARENT,Gate0_independent_audit='ACCEPTED_CONDITIONAL',primary_architecture_decision=outcome,requirement_frontier=frontier,bearing_requirement_fixed_baseline=fixed,primary_sigma_0p10_requirement=primary,
  primary_selected_crossing_angle_range=({k:selected[0][k] for k in ['baseline_km','alpha_min_deg','alpha_max_deg','worst_range_km','worst_p95_relative']} if selected else None),
  requirement_scope='IDEAL_KNOWN_POSITION_DIRECTED_BEARING_GEOMETRY; DISCRETE_11_RANGES; FINITE_MC; NOT_HARDWARE_OR_CONTINUOUS_GUARANTEE',actual_bearing_hardware_capability='UNKNOWN',
  single_array_cold_start='NOT_ESTABLISHED',single_array_depth_engineering='CLOSED',AUX1='CONDITIONAL_ARCHITECTURE_CANDIDATE',R4_A1_NEW='NOT_OPENED',depth='CLOSED_PENDING_INDEPENDENT_HORIZONTAL_ACQUISITION',R4_progress_percent=0,
  A2_A3_A4_B2_B3_B4_C='PAUSED / NOT_OPENED',P5='NOT_OPENED',new_acoustic_propagation=0,new_depth_scores=0,position_bias_time_error_added=False,tracking_filter_implemented=False,audit_status='PENDING_RESEARCH_LEAD_AUDIT',next_stage='STOP_FOR_RESEARCH_LEAD_AUDIT; GATE2_NOT_OPENED')

def execute():
 p=load(OUT/'AUX_GATE1_DESIGN_FREEZE.json')
 for name,d in {**p['bindings'],**p['historical_bindings'],**p['external_bindings']}.items():assert sha(name)==d,name
 head=git('rev-parse','HEAD');remote=git('ls-remote','origin','refs/heads/main').split()[0]
 assert head==remote and git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==PARENT
 assert not (OUT/'AUX_GATE1_EXECUTION_START.json').exists(),'Frozen one-shot execution only'
 dump(OUT/'AUX_GATE1_EXECUTION_START.json',dict(design_sha=head,remote_verified_before_mc=True,policy_sha256=sha(OUT/'AUX_GATE1_DESIGN_FREEZE.json')))
 z=np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_samples'],2));np.savez_compressed(OUT/'AUX_FROZEN_STANDARD_NORMAL_DRAWS.npz',z=z)
 results=[]
 for i,row in enumerate(read_csv(OUT/'AUX_CROSSTRACK_GRID.csv'),1):
  results.append(evaluate(row,z,p))
  if i%220==0:print(f'geometry cells completed {i}/1870',flush=True)
 write_csv(OUT/'AUX_CROSSTRACK_RESULTS.csv',results);summaries=summarize(results,p);write_csv(OUT/'AUX_FULL_RANGE_SUMMARY.csv',summaries)
 frontier,fixed=frontiers(summaries,p);write_csv(OUT/'AUX_BASELINE_REQUIREMENT_FRONTIER.csv',frontier);write_csv(OUT/'AUX_BEARING_REQUIREMENT_AT_FIXED_BASELINE.csv',fixed)
 d=decide(summaries,frontier,fixed,p);d['design_sha']=head;dump(OUT/'AUX_GATE1_DECISION.json',d)
 print(json.dumps(dict(decision=d['primary_architecture_decision'],frontier=frontier)),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if args.execute:execute()
 else:ap.error('--execute required')
