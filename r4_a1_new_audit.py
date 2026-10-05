"""Independent scene, slope initializer, cold metrics and finite-difference optimizer audit."""
from pathlib import Path
import argparse,ast,csv,json,hashlib,math,subprocess,xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
OUT=Path('results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC')
METRICS=['range_error','bearing_error_deg','speed_error','heading_error_deg']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,data):
 with Path(p).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def flag(v):return v is True or v=='True'
def w(v):return np.arctan2(np.sin(v),np.cos(v))
def audit():
 p=load(OUT/'A1_NEW_DESIGN_FREEZE.json');checks=[]
 def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
 def near(name,a,b,rtol=2e-7,atol=1e-7):check(name,math.isclose(float(a),float(b),rel_tol=rtol,abs_tol=atol))
 for name,h in {**p['bindings'],**p['historical_bindings']}.items():check('hash:'+name,sha(name)==h)
 panel=rows(OUT/'A1_NEW_TRUTH_PANEL.csv');check('12truths',len(panel)==12)
 for x in panel:
  for name,step in [('r0_km',1),('theta0_deg',.5),('v_mps',.2),('psi_deg',1)]:check('offgrid:'+x['case_id']+name,abs(float(x[name])/step-round(float(x[name])/step))>1e-8)
 tree=ast.parse(Path('r4_a1_new_estimator.py').read_text(encoding='utf-8'));fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='estimate');check('truth_free_API',[x.arg for x in fn.args.args]==['times','node_positions','bearings','sigma_deg']);check('no_generator_or_truth_symbol',not any(isinstance(x,ast.Name) and ('truth' in x.id.lower() or x.id in ['scene','evaluate']) for x in ast.walk(tree)))
 start=load(OUT/'A1_NEW_EXECUTION_START.json');check('pushed_before_MC',start['remote_verified_before_MC'] and start['policy_sha256']==sha(OUT/'A1_NEW_DESIGN_FREEZE.json'))
 check('design_git_parent',subprocess.check_output(['git','rev-parse',start['design_SHA']+'^'],text=True).strip()==p['commit_parent_SHA'])
 check('design_title',subprocess.check_output(['git','log','-1','--format=%s',start['design_SHA']],text=True).strip()=='R4 A1-new: freeze augmented off-grid dynamic baseline')
 seal=load(OUT/'A1_NEW_PRIMARY_SEAL.json');check('primary_sealed_before_secondary',seal['before_secondary'] is True)
 for name,h in seal['bindings'].items():check('primaryseal:'+name,sha(name)==h)
 data=rows(OUT/'A1_NEW_RUN_RESULTS.csv');diagnostics=rows(OUT/'A1_NEW_ESTIMATOR_DIAGNOSTICS.csv');initial=rows(OUT/'A1_NEW_INITIALIZER_DIAGNOSTICS.csv');check('12000runs',len(data)==len(diagnostics)==len(initial)==12000)
 t=np.arange(121)*10.;main=np.array([[2*tt,0.] if tt<=600 else [1200+2*(tt-600)*math.cos(math.pi/12),2*(tt-600)*math.sin(math.pi/12)] for tt in t]);recomputed=[];cold=[]
 for aname in ['A','B']:
  a=p['anchors'][aname];draws=np.load(OUT/f'A1_NEW_DRAWS_{aname}.npz');u=draws['uniforms'];z=draws['normals'];rng=np.random.Generator(np.random.PCG64(a['seed']));check(aname+'uniform_regeneration',np.array_equal(u,rng.random((6000,3))));check(aname+'normal_regeneration',np.array_equal(z,rng.standard_normal((6000,121,6))))
  group=[x for x in data if x['anchor']==aname];saved=rows(OUT/f'A1_NEW_RUN_RESULTS_{aname}.csv');check(aname+'combined_identity',group==saved);check(aname+'6000unique',len({(x['case_id'],int(x['realization'])) for x in group})==6000)
  for ix,x in enumerate(group):
   truth=panel[ix//500];ixall=ix+(0 if aname=='A' else 6000);diag=diagnostics[ixall];ini=initial[ixall];key=aname+':'+str(ix)
   check(key+'case_draw',x['case_id']==truth['case_id'] and int(x['draw_index'])==ix and int(x['realization'])==ix%500);side=-1 if ix%2==0 else 1;check(key+'mirror',int(x['mirror'])==side)
   c=(u[ix,0]*2-1)*a['common_deg'];diff=(u[ix,1]*2-1)*a['half_diff_deg'];beta=(u[ix,2]*2-1)*a['beta_deg'];angle=math.radians(90+beta);offset=np.array([a['baseline_m']*math.cos(angle),side*a['baseline_m']*math.sin(angle)]);nodes=np.stack([main,main+offset],axis=1);estimated=nodes+a['nav_sigma_m']*z[ix,:,2:].reshape(121,2,2)
   r=float(truth['r0_km'])*1000;th=math.radians(float(truth['theta0_deg']));v=float(truth['v_mps']);psi=math.radians(float(truth['psi_deg']));truth_state=np.array([r*math.cos(th),r*math.sin(th),v*math.cos(psi),v*math.sin(psi)]);target=truth_state[:2]+t[:,None]*truth_state[2:];delta=target[:,None,:]-nodes;angles=np.arctan2(delta[:,:,1],delta[:,:,0])+np.deg2rad([c-diff,c+diff])+math.radians(a['sigma_deg'])*z[ix,:,:2]
   for field,val in [('b_common_deg',c),('b_diff_HALF_deg',diff),('beta_deg',beta)]:near(key+field,x[field],val)
   for k in ['anchor','case_id','realization','draw_index','mirror']:check(key+'diag_key'+k,diag[k]==ini[k]==x[k])
   # Independent tangent-slope intersection. No primary initializer import.
   m1=np.tan(angles[:,0]);m2=np.tan(angles[:,1]);dd=estimated[:,1]-estimated[:,0]
   with np.errstate(divide='ignore',invalid='ignore'):
    xx=(dd[:,1]-m2*dd[:,0])/(m1-m2);yy=m1*xx;xy=estimated[:,0]+np.column_stack([xx,yy]);l1=xx*np.cos(angles[:,0])+yy*np.sin(angles[:,0]);l2=(xx-dd[:,0])*np.cos(angles[:,1])+(yy-dd[:,1])*np.sin(angles[:,1])
   det=np.sin(angles[:,1]-angles[:,0]);valid=(np.abs(det)>1e-12)&(l1>0)&(l2>0)&np.isfinite(xy).all(axis=1)
   if flag(x['failed']):
    for m in METRICS:check(key+'failed_inf'+m,float(x[m])==float('inf'))
    recomputed.append(dict(anchor=aname,case_id=truth['case_id'],failed=True,**{m:float('inf') for m in METRICS}));continue
   check(key+'valid_count',int(diag['valid_intersections'])==int(valid.sum()) and int(diag['invalid_intersections'])==int((~valid).sum()))
   T=np.column_stack([np.ones(valid.sum()),t[valid]/1200]);points=xy[valid];coef=np.linalg.lstsq(T,points,rcond=None)[0]
   for _ in range(5):
    norms=np.sqrt(np.sum((points-T@coef)**2,axis=1));scale=max(10.,float(np.median(norms)));weight=np.minimum(1.,1.5*scale/np.maximum(norms,1e-15));coef=np.linalg.lstsq(T*np.sqrt(weight)[:,None],points*np.sqrt(weight)[:,None],rcond=None)[0]
   init_state=np.r_[coef[0],coef[1]/1200]
   for k,value in zip(['initial_x_m','initial_y_m','initial_vx_mps','initial_vy_mps'],init_state):near(key+k,diag[k],value);near(key+'init_table'+k,ini[k],value)
   state=np.array([float(x[k]) for k in ['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps']]);rel=state[:2]-estimated[0,0];rh=np.hypot(*rel);thh=np.arctan2(rel[1],rel[0]);vh=np.hypot(*state[2:]);ph=np.arctan2(state[3],state[2]);values=dict(range_error=abs(rh-r)/r,bearing_error_deg=abs(math.degrees(float(w(thh-th)))),speed_error=abs(vh-v)/v,heading_error_deg=abs(math.degrees(float(w(ph-psi)))))
   for k,value in values.items():near(key+'metric'+k,x[k],value)
   for k,value in [('r_hat_km',rh/1000),('theta_hat_deg',math.degrees(thh)),('v_hat_mps',vh),('psi_hat_deg',math.degrees(ph))]:near(key+'estimate'+k,x[k],value)
   delta_hat=(state[:2]+t[:,None]*state[2:])[:,None,:]-estimated;rr=w(np.arctan2(delta_hat[:,:,1],delta_hat[:,:,0])-angles);near(key+'residual_RMS',diag['residual_RMS_deg'],math.degrees(float(np.sqrt(np.mean(rr**2)))))
   check(key+'solver_valid',flag(diag['solver_success']) and int(diag['nfev'])<=100 and np.isfinite(state).all())
   if ix%500 in [0,499]:
    scaling=np.array([50000.,50000.,2.,2.])
    def objective(q):
     ss=q*scaling;dlt=(ss[:2]+t[:,None]*ss[2:])[:,None,:]-estimated;return (w(np.arctan2(dlt[:,:,1],dlt[:,:,0])-angles)/math.radians(a['sigma_deg'])).ravel()
    fit=least_squares(objective,init_state/scaling,jac='3-point',loss='soft_l1',f_scale=1.5,method='trf',max_nfev=100,ftol=1e-10,xtol=1e-10,gtol=1e-10);difference=np.abs(fit.x*scaling-state);check(key+'independent_fd_solver',fit.success and np.max(difference[:2])<=1 and np.max(difference[2:])<=.001);cold.append(dict(anchor=aname,case_id=truth['case_id'],realization=ix%500,max_position_difference_m=float(max(difference[:2])),max_velocity_difference_mps=float(max(difference[2:])),passed=True))
   recomputed.append(dict(anchor=aname,case_id=truth['case_id'],failed=False,**values))
   if (ix+1)%500==0:print(f'independent {aname} {ix+1}/6000',flush=True)
 summary=rows(OUT/'A1_NEW_CASE_SUMMARY.csv');check('24case_summaries',len(summary)==24)
 for row in summary:
  group=[x for x in recomputed if x['anchor']==row['anchor'] and x['case_id']==row['case_id']];check('500percase:'+str(row),len(group)==int(row['n_runs'])==500)
  near('failure_rate:'+str(row),row['failure_rate'],sum(x['failed'] for x in group)/500)
  for m in METRICS:
   ordered=sorted(x[m] for x in group)
   for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:near('quantile:'+row['anchor']+row['case_id']+m+tag,row[m+'_'+tag],ordered[math.ceil(500*q)-1])
  for tier,limits in [('STRONG',[.05,.5,.05,2]),('PROJECT',[.1,1,.1,5])]:check('gate:'+row['anchor']+row['case_id']+tier,flag(row[tier+'_pass'])==all(float(row[m+'_P95'])<=limit for m,limit in zip(METRICS,limits)))
 glob=rows(OUT/'A1_NEW_GLOBAL_SUMMARY.csv')
 for row in glob:
  group=[x for x in recomputed if x['anchor']==row['anchor']];case=[x for x in summary if x['anchor']==row['anchor']]
  for m in METRICS:
   near('globalworst:'+row['anchor']+m,row[m+'_worst_case_P95'],max(float(x[m+'_P95']) for x in case));ordered=sorted(x[m] for x in group)
   for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:near('pooled:'+row['anchor']+m+tag,row[m+'_pooled_'+tag],ordered[math.ceil(6000*q)-1])
  for tier in ['STRONG','PROJECT']:check('globalgate:'+row['anchor']+tier,flag(row[tier+'_pass'])==all(flag(x[tier+'_pass']) for x in case))
 d=load(OUT/'A1_NEW_DECISION.json');A=glob[0];outcome='R4_A1_NEW_STRONG_BASELINE_ESTABLISHED' if flag(A['STRONG_pass']) else 'R4_A1_NEW_PROJECT_BASELINE_ESTABLISHED' if flag(A['PROJECT_pass']) else 'R4_A1_NEW_BASELINE_NOT_ESTABLISHED';check('primary_only_decision',d['scientific_decision']==outcome and d['primary_application_scenario']=='Anchor A');credit=15 if flag(A['PROJECT_pass']) else 0;check('progress_and_stop',d['R4_A1_progress_percent']==d['R4_overall_percent']==credit and d['depth']=='NOT_OPENED' and d['stop_after_commit_B'] is True and not d['automatic_next_experiment'])
 xml=ET.parse(OUT/'A1_NEW_TESTS.xml').getroot();check('5tests_pass',len(xml.findall('.//testcase'))==5 and not xml.findall('.//failure') and not xml.findall('.//error'))
 for f in ['results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv']:
  before=subprocess.check_output(['git','show',p['commit_parent_SHA']+':'+f]).decode('utf-8').replace('\r\n','\n');check('append_only:'+f,Path(f).read_text(encoding='utf-8').startswith(before))
 return checks,cold

def present():
 d=load(OUT/'A1_NEW_DECISION.json');summary=rows(OUT/'A1_NEW_CASE_SUMMARY.csv');runs=rows(OUT/'A1_NEW_RUN_RESULTS.csv');panel=rows(OUT/'A1_NEW_TRUTH_PANEL.csv');truth={x['case_id']:x for x in panel}
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figs=OUT/'figures';figs.mkdir(exist_ok=True);plt.rcParams['svg.hashsalt']='A1_NEW'
 def save(fig,name):
  fig.savefig(figs/(name+'.png'),dpi=160,bbox_inches='tight');fig.savefig(figs/(name+'.svg'),bbox_inches='tight',metadata={'Date':None});plt.close(fig);f=figs/(name+'.svg');f.write_text('\n'.join(x.rstrip() for x in f.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 fig,axes=plt.subplots(2,2,figsize=(12,7))
 for ax,m,scale,strong,project,label in zip(axes.flat,METRICS,[100,1,100,1],[5,.5,5,2],[10,1,10,5],['Range P95 (%)','Bearing P95 (deg)','Speed P95 (%)','Heading P95 (deg)']):
  for a in ['A','B']:
   g=[x for x in summary if x['anchor']==a];ax.plot(range(1,13),[scale*float(x[m+'_P95']) for x in g],'o-',label=a)
  ax.axhline(strong,color='green',linestyle='--',label='STRONG');ax.axhline(project,color='black',linestyle=':',label='PROJECT');ax.set(xlabel='Off-grid case',ylabel=label);ax.legend();ax.grid(alpha=.2)
 fig.tight_layout();save(fig,'FIG1_FOUR_PARAMETER_P95')
 fig,axes=plt.subplots(1,2,figsize=(11,4))
 for ax,f,tf,label in zip(axes,['r_hat_km','v_hat_mps'],['r0_km','v_mps'],['Range km','Speed m/s']):
  for a in ['A','B']:
   g=[x for x in runs if x['anchor']==a and not flag(x['failed'])];ax.scatter([float(truth[x['case_id']][tf]) for x in g],[float(x[f]) for x in g],s=2,alpha=.1,label=a)
  ax.set(xlabel='True '+label,ylabel='Estimated '+label);ax.legend()
 fig.tight_layout();save(fig,'FIG2_ESTIMATED_VS_TRUE')
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for ax,m in zip(axes,['bearing_error_deg','heading_error_deg']):
  for a in ['A','B']:
   vals=sorted(float(x[m]) for x in runs if x['anchor']==a);ax.plot(vals,np.arange(1,len(vals)+1)/len(vals),label=a)
  ax.set(xlabel=m,ylabel='Unconditional empirical CDF');ax.legend();ax.grid(alpha=.2)
 fig.tight_layout();save(fig,'FIG3_ANGLE_ERROR_DISTRIBUTIONS')
 fig,axes=plt.subplots(1,4,figsize=(13,4))
 for ax,m,scale in zip(axes,METRICS,[100,1,100,1]):ax.bar(['A','B'],[scale*float(d[k][m+'_worst_case_P95']) for k in ['primary_summary','secondary_summary']]);ax.set_title(m);ax.set_ylabel('%' if scale==100 else 'deg')
 fig.tight_layout();save(fig,'FIG4_ANCHOR_COMPARISON')
 lines=['# R4-A1-new application pre-research dynamic baseline','',d['scientific_decision'],'','Project stage: APPLICATION_PRE_RESEARCH. Hardware UNKNOWN is a nonblocking design-assumption boundary, not an actual capability claim. This new authorization supersedes the previous client-confirmation numerical STOP for this stage only. Scientific baseline333ead6 retained; chronological parent0eef9ff retained.','', 'Primary A is the sole verdict source; B is SECONDARY_CONFIRMATION_ONLY.12off-grid truths x500realizations per anchor;121epochs;0..1200s;dt10s;MAIN2m/s,15deg turn after600s. AUX has identical displacement/velocity to MAIN and a fixed global lateral offset relative to nominal initial LoB0deg plus sampled beta, mirrored by side. No target-guided formation control.','', 'Per realization common/HALF-diff/beta sampled independently Uniform over frozen bounds. This APPLICATION_PRE_RESEARCH_BOUNDED_SYSTEMATIC_ENSEMBLE is an assumed statistical design model, not an observed physical distribution or a claim that Gate2B certified the continuous interior. Navigation: independent epoch-wise Gaussian per node per axis; not a specific INS/GNSS temporal model. Bearing random errors independent epoch/node Gaussian. All variates/seeds frozen before MC; no resampling failures.','', 'Observation-only estimator: positive directed-ray intersection initializer, mathematically invalid intersections excluded, fivefixed IRLS line fits with all valid points retained, then four Cartesian initial-position/velocity parameters optimized using wrapped dual-node bearing residual, soft_l1loss1.5sigma, frozen analytic Jacobian and convergence settings. No truth warm start, parameter-grid starts, acoustic model, depth or nuisance truth. r0/theta0 reference estimated MAIN position at t0; v/psi from fitted global velocity. Truth only in generator/evaluator.','', 'Nearest-rank unconditional errors; optimizer/initializer failures count infinite error in all metrics. Per-case and pooled median/P90/P95/P99 retained. Primary Gate requires all12cases and uses worst case-level P95; no pooled substitution. Local covariance/condition/residual diagnostics are approximate and do not account for systematic/navigation model mismatch; Monte Carlo truth errors determine performance.','', '| Anchor | Worst range P95 % | Bearing P95 deg | Speed P95 % | Heading P95 deg | STRONG | PROJECT |','|---|---:|---:|---:|---:|---|---|']
 for k in ['primary_summary','secondary_summary']:
  x=d[k];lines.append(f"| {x['anchor']} | {100*float(x['range_error_worst_case_P95']):.6f} | {float(x['bearing_error_deg_worst_case_P95']):.6f} | {100*float(x['speed_error_worst_case_P95']):.6f} | {float(x['heading_error_deg_worst_case_P95']):.6f} | {x['STRONG_pass']} | {x['PROJECT_pass']} |")
 lines+=['','STRONG:range/speed5%,bearing0.5deg,heading2deg. PROJECT:range/speed10%,bearing1deg,heading5deg. All12cases required. Credit15only if primary PROJECT passes, pending research-lead independent audit; otherwise0. B cannot rescue primary. No application claim of five-parameter/multi-parameter<10% from these four horizontal submetrics.','', 'Independent audit regenerates all draws, reconstructs all scenes with independent geometry, slope intersections and IRLS, cold-recomputes all12000errors/quantiles and decision.48frozen diagnostic runs independently optimized using3-point finite-difference Jacobian; this is saved-scene verification, not extra performance MC.','', 'Single-array conditional depth engineering stays CLOSED; exact-horizontal depth mechanism SUPPORTED_ORACLE_ONLY. Depth may be reconsidered only as a separately authorized augmented-posterior stage after audit. No KRAKEN/BELLHOP/TL/CZ/SSP or depth score. Stop after execution commit/push; A2/A3/A4 are next roadmap stages only, no automatic experiment.','', '```json',json.dumps(d,indent=2),'```']
 text='\n'.join(lines)+'\n'
 for name in ['A1_NEW_REPORT.md','GPT_SYNC.md']:(OUT/name).write_text(text,encoding='utf-8')
 progress=dict(project_stage='PROJECT_APPLICATION_PRE_RESEARCH',active_architecture='MAIN_TOWED_ARRAY_PLUS_MOBILE_PASSIVE_AUXILIARY',active_stage='R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC_LOCALIZATION',hardware_UNKNOWN_blocks_numerical_research=False,scientific_decision=d['scientific_decision'],A1_NEW_progress_percent=d['R4_A1_progress_percent'],overall_progress_percent=d['R4_overall_percent'],independent_audit='PENDING_RESEARCH_LEAD_AUDIT',next=d['next'],automatic_next_stage=False,depth='NOT_OPENED',single_array_depth_engineering='CLOSED',exact_horizontal_depth_mechanism='SUPPORTED_ORACLE_ONLY',legacy_progress_file='R4_PROGRESS.json IS HISTORICAL_FROZEN_SINGLE_ARRAY_SNAPSHOT',scientific_parent_SHA=d['scientific_parent_SHA'],commit_parent_SHA=d['commit_parent_SHA'])
 dump('results/R4_MASTER/R4_APPLICATION_PROGRESS.json',progress)
 with Path('results/R4_MASTER/R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## R4-A1-new execution\n\n'+d['scientific_decision']+'; PENDING_RESEARCH_LEAD_AUDIT. APPLICATION_PRE_RESEARCH design assumptions; hardware UNKNOWN nonblocking. [Active application progress](R4_APPLICATION_PROGRESS.json). A1-new/overall scientific credit='+str(d['R4_overall_percent'])+'%. [Report](../R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_REPORT.md). B secondary only; depth not opened. Commit/push/verify then STOP for audit; no automatic A2.\n')
 with Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow(['A1-NEW','../R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_DECISION.json',d['scientific_decision'],'Frozen application assumptions; primary12off-grid x500; four horizontal metrics; Bsecondary','PENDING_RESEARCH_LEAD_AUDIT; SCIENTIFIC_CREDIT_PER_NEW_AUTHORIZATION',d['R4_overall_percent']])
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--finalize',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 if args.finalize:
  assert not (OUT/'VALIDATION.json').exists(),'one finalize'
  checks,cold=audit();present();write(OUT/'A1_NEW_INDEPENDENT_CHECKS.csv',checks);write(OUT/'A1_NEW_INDEPENDENT_SOLVER_CHECKS.csv',cold);dump(OUT/'VALIDATION.json',dict(status='PASS',checks=len(checks),failed=0,tests_passed=5,primary_runs=6000,secondary_runs=6000,cold_metrics=12000,independent_FD_solver_runs=48,new_acoustic_propagation=0,new_depth_score=0))
  paths=[x for x in OUT.rglob('*') if x.is_file()]+[Path('results/R4_MASTER/R4_PLAN.md'),Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'),Path('results/R4_MASTER/R4_APPLICATION_PROGRESS.json')];dump(OUT/'A1_NEW_RESULTS_FREEZE.json',dict(bindings={x.as_posix():sha(x) for x in sorted(paths)}));print(json.dumps(dict(finalized='PASS',checks=len(checks))),flush=True)
 elif args.verify:
  for name,h in load(OUT/'A1_NEW_RESULTS_FREEZE.json')['bindings'].items():assert sha(name)==h,name
  checks,cold=audit();print(json.dumps(dict(sealed='PASS',checks=len(checks),failed=0)),flush=True)
 else:ap.error('--finalize or --verify required')
