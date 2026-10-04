"""Preregistered paired-bearing geometry only; no acoustic model imports."""
from pathlib import Path
import argparse,csv,hashlib,json,math,subprocess
import numpy as np
OUT=Path('results/R4_ROUTE_RESET_AUXILIARY_NODE_GATE0')
PARENT='6b8d57ef6082f24720c9ece654b9cd710a687e13'
TITLE='R4 route reset: freeze auxiliary-node geometry Gate0'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read_json(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def write_json(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read_csv(p):
 with Path(p).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def write_csv(p,rows):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def git(*a):return subprocess.run(['git',*a],check=True,capture_output=True,text=True,encoding='utf-8').stdout.strip()
def geometry(r,b,alpha,branch,side):
 a=math.radians(alpha);disc=b*b-r*r*math.sin(a)**2
 if disc < -1e-12:return None
 d=r*math.cos(a)+(1 if branch=='FAR' else -1)*math.sqrt(max(0,disc))
 if d<=0:return None
 return np.array([r-d*math.cos(a),side*d*math.sin(a)]),d

def triangulate(p2,b1,b2,tol=1e-12):
 u1=np.column_stack((np.cos(b1),np.sin(b1)));u2=np.column_stack((np.cos(b2),np.sin(b2)))
 det=u1[:,0]*u2[:,1]-u1[:,1]*u2[:,0];parallel=np.abs(det)<=tol;safe=np.where(parallel,np.nan,det)
 t1=(p2[0]*u2[:,1]-p2[1]*u2[:,0])/safe;t2=(p2[0]*u1[:,1]-p2[1]*u1[:,0])/safe
 xy=t1[:,None]*u1;behind=(t1<=0)|(t2<=0);numeric=~np.isfinite(xy).all(axis=1)
 return xy,parallel|behind|numeric,parallel,behind,numeric,det

def covariance(r,p2,s):
 q=np.array([r,0.])-p2;H=np.array([[0.,1/r],np.array([-q[1],q[0]])/(q@q)])
 F=H.T@H/math.radians(s)**2;return H,F,np.linalg.inv(F)
def quantile(v,p):return float(np.sort(v)[math.ceil(len(v)*p)-1])
def serial(v):return 'INF' if not math.isfinite(v) else v

def evaluate(row,z,policy):
 r,b,a,s=[float(row[k]) for k in ['range_km','baseline_km','alpha_deg','sigma_deg']];side=int(row['side']);g=geometry(r,b,a,row['branch'],side)
 basekeys=['aux_x_km','aux_y_km','aux_range_km','unique_noiseless_finite_forward_intersection','jacobian_condition','fim_condition','linear_range_sigma_km','linear_cross_sigma_km','linear_range_p95_relative','n_samples','failure_rate','behind_sensor_rate','parallel_rate','ill_conditioned_rate','numerical_failure_rate','fraction_range_le1pct','fraction_range_le5pct','fraction_range_le10pct','fraction_le5pct_wilson_low','fraction_le5pct_wilson_high','relative_range_p95_mc_ci_low','relative_range_p95_mc_ci_high']
 result=dict(row);result.update({k:'' for k in basekeys});result.update({f'{m}_{p}':'' for m in ['range_abs_km','relative_range','cross_range_abs_km','position_2d_km'] for p in ['median','p90','p95','p99']})
 result['alpha_max_deg']=math.degrees(math.asin(b/r));result['geometry_status']='INFEASIBLE_TRIANGLE' if g is None else 'FEASIBLE'
 if g is None:return result
 p2,d=g;H,F,C=covariance(r,p2,s)
 result.update(aux_x_km=float(p2[0]),aux_y_km=float(p2[1]),aux_range_km=d,unique_noiseless_finite_forward_intersection=True,jacobian_condition=float(np.linalg.cond(H)),fim_condition=float(np.linalg.cond(F)),linear_range_sigma_km=math.sqrt(C[0,0]),linear_cross_sigma_km=math.sqrt(C[1,1]),linear_range_p95_relative=1.959963984540054*math.sqrt(C[0,0])/r)
 noise=z*math.radians(s)*side
 xy,failed,parallel,behind,numeric,det=triangulate(p2,noise[:,0],-side*math.radians(a)+noise[:,1],policy['parallel_det_tolerance'])
 er=np.abs(np.linalg.norm(xy,axis=1)-r)/r;cross=np.abs(xy[:,1]);pos=np.linalg.norm(xy-np.array([r,0.]),axis=1)
 for arr in [er,cross,pos]:arr[failed]=np.inf
 for m,arr in {'range_abs_km':er*r,'relative_range':er,'cross_range_abs_km':cross,'position_2d_km':pos}.items():
  for tag,q in [('median',.5),('p90',.9),('p95',.95),('p99',.99)]:result[f'{m}_{tag}']=serial(quantile(arr,q))
 n=len(z);phat=float((er<=.05).mean());zz=1.959963984540054;den=1+zz*zz/n;center=(phat+zz*zz/(2*n))/den;half=zz*math.sqrt(phat*(1-phat)/n+zz*zz/(4*n*n))/den
 ordered=np.sort(er);delta=zz*math.sqrt(n*.95*.05);lo=max(0,math.floor(n*.95-delta)-1);hi=min(n-1,math.ceil(n*.95+delta)-1)
 result.update(n_samples=n,failure_rate=float(failed.mean()),behind_sensor_rate=float(behind.mean()),parallel_rate=float(parallel.mean()),ill_conditioned_rate=float((1/np.maximum(np.abs(det),np.finfo(float).tiny)>policy['ill_condition_proxy_limit']).mean()),numerical_failure_rate=float(numeric.mean()),fraction_range_le1pct=float((er<=.01).mean()),fraction_range_le5pct=phat,fraction_range_le10pct=float((er<=.1).mean()),fraction_le5pct_wilson_low=center-half,fraction_le5pct_wilson_high=center+half,relative_range_p95_mc_ci_low=serial(float(ordered[lo])),relative_range_p95_mc_ci_high=serial(float(ordered[hi])))
 return result

def summarize(rows,p):
 summaries=[];mapping=[]
 for b in p['baseline_km']:
  for a in p['alpha_deg']:
   for s in p['sigma_deg']:
    good=[x for x in rows if float(x['baseline_km'])==b and float(x['alpha_deg'])==a and float(x['sigma_deg'])==s and x['geometry_status']=='FEASIBLE'];ranges=sorted({float(x['range_km']) for x in good});worst=max((float(x['relative_range_p95']) for x in good),default=None)
    summaries.append(dict(baseline_km=b,alpha_deg=a,sigma_deg=s,n_feasible_cells=len(good),feasible_ranges_km=';'.join(f'{r:g}' for r in ranges),all_three_ranges_feasible=ranges==p['range_km'],worst_feasible_p95_relative_range=worst if worst is not None else '',worst_failure_rate=max((float(x['failure_rate']) for x in good),default=''),all_feasible_cells_strong=worst is not None and worst<=.01,all_feasible_cells_usable=worst is not None and worst<=.05,all_feasible_cells_project_level=worst is not None and worst<=.1,full_range_usable=ranges==p['range_km'] and worst is not None and worst<=.05,full_range_project_level=ranges==p['range_km'] and worst is not None and worst<=.1))
    for r in p['range_km']:
     group=[x for x in good if float(x['range_km'])==r];worst_r=max((float(x['relative_range_p95']) for x in group),default=None)
     mapping.append(dict(range_km=r,baseline_km=b,alpha_deg=a,sigma_deg=s,n_feasible=len(group),worst_p95_relative_range=worst_r if worst_r is not None else '',worst_failure_rate=max((float(x['failure_rate']) for x in group),default=''),status='FEASIBLE' if group else 'INFEASIBLE_TRIANGLE',le5pct=worst_r is not None and worst_r<=.05,le10pct=worst_r is not None and worst_r<=.1))
 return summaries,mapping

def decide(summaries,mapping,p):
 primary=[x for x in summaries if float(x['sigma_deg'])==.1]
 admitted=[x for x in primary if x['full_range_usable'] and float(x['baseline_km'])<=p['practical_baseline_ceiling_km']]
 if admitted:outcome='AUXILIARY_PASSIVE_NODE_GEOMETRY_ADMITTED'
 elif any(x['all_feasible_cells_usable'] or x['full_range_project_level'] for x in primary):outcome='AUXILIARY_NODE_GEOMETRY_CONDITIONALLY_PRACTICAL'
 else:outcome='SINGLE_AUXILIARY_BEARING_NODE_NOT_SUFFICIENT'
 best=[]
 for b in p['baseline_km']:
  options=[x for x in primary if x['baseline_km']==b and x['n_feasible_cells']];full=[x for x in options if x['all_three_ranges_feasible']]
  best.append(dict(baseline_km=b,best_registered_feasible_region=min(options,key=lambda x:float(x['worst_feasible_p95_relative_range'])) if options else None,best_full_50_60km_region=min(full,key=lambda x:float(x['worst_feasible_p95_relative_range'])) if full else None))
 return dict(stage=p['stage'],parent_sha=PARENT,B1B_independent_audit='ACCEPTED_NEGATIVE_ROUTE_CLOSED',CURRENT_SINGLE_ARRAY_EVIDENCE_SUPPORTS_TARGET=False,auxiliary_node_decision=outcome,recommended_architecture='AUX1_CONDITIONAL_CANDIDATE; NO_SYSTEM_OR_DEPTH_ADMISSION',sigma_0p1_baseline_best=best,primary_feasible_regions=[x for x in mapping if x['sigma_deg']==.1 and x['status']=='FEASIBLE'],full_50_60km_usable_regions=[x for x in primary if x['full_range_usable']],full_50_60km_project_level_regions=[x for x in primary if x['full_range_project_level']],practical_baseline_ceiling_km=p['practical_baseline_ceiling_km'],practical_scale_status='PLANNING_ASSUMPTION_NOT_VERIFIED_HARDWARE_CAPABILITY',depth_route='CLOSED_UNDER_CURRENT_SINGLE_ARRAY_ARCHITECTURE; MAY_REOPEN_AFTER_INDEPENDENT_HORIZONTAL_ACQUISITION',R4_original_A_B_sequence='PAUSED_PENDING_ARCHITECTURE_RESET',R4_progress_percent=0,A2_A3_A4_B2_B3_B4_C='PAUSED / NOT_OPENED',P5='NOT_OPENED',next_stage='RESEARCH_LEAD_ARCHITECTURE_DECISION_REQUIRED; NO_AUTOMATIC_RESTART',new_acoustic_propagation_evaluations=0,new_depth_scores=0,full_tracker_implemented=False,audit_status='PENDING_RESEARCH_LEAD_AUDIT')

def execute():
 p=read_json(OUT/'AUX_GATE0_DESIGN_FREEZE.json')
 for path,d in {**p['bindings'],**p['historical_bindings']}.items():assert sha(path)==d,path
 head=git('rev-parse','HEAD');remote=git('ls-remote','origin','refs/heads/main').split()[0]
 assert head==remote and git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==PARENT
 assert not (OUT/'AUX_EXECUTION_START.json').exists(),'One-shot execution; no retuning or rerun'
 write_json(OUT/'AUX_EXECUTION_START.json',dict(design_sha=head,remote_verified_before_geometry=True,policy_sha256=sha(OUT/'AUX_GATE0_DESIGN_FREEZE.json')))
 z=np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_samples'],2));np.savez_compressed(OUT/'AUX_FROZEN_STANDARD_NORMAL_DRAWS.npz',z=z)
 rows=[evaluate(row,z,p) for row in read_csv(OUT/'AUX_GEOMETRY_GRID.csv')];write_csv(OUT/'AUX_GEOMETRY_RESULTS.csv',rows)
 summaries,mapping=summarize(rows,p);write_csv(OUT/'AUX_GEOMETRY_SUMMARY.csv',summaries);write_csv(OUT/'AUX_P95_RANGE_MAP.csv',mapping)
 d=decide(summaries,mapping,p);d['design_sha']=head;write_json(OUT/'AUX_GATE0_DECISION.json',d)
 print(json.dumps(dict(decision=d['auxiliary_node_decision'],requested=len(rows),feasible=sum(x['geometry_status']=='FEASIBLE' for x in rows),samples_per_cell=len(z))))

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--execute',action='store_true');args=parser.parse_args()
 if args.execute:execute()
 else:parser.error('Explicit --execute required')
