"""Gate2B independent slope reconstruction and preregistered result presentation."""
from pathlib import Path
import argparse,csv,hashlib,json,math,subprocess,itertools,xml.etree.ElementTree as ET
import numpy as np
OUT=Path('results/R4_AUX_GATE2B_JOINT_ERROR_BUDGET')
FAMILIES=['T5_BUDGET_FAMILY','T10_BUDGET_FAMILY']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,data):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def flag(x):return x is True or x=='True'
def key(x):return (x['anchor'],x['family'],float(x['lambda_value']),int(x['sign_common']),int(x['sign_diff']),int(x['sign_beta']),float(x['range_km']),int(x['side']))
def zero_errors(z,a,r,side):
 b=a['baseline_km'];s=math.radians(a['sigma_random_deg']);m1=np.tan(side*s*z[:,0]);m2=np.tan(math.atan2(-side*b,r)+side*s*z[:,1]);xx=side*b/(m1-m2);yy=m1*xx
 return np.abs(np.hypot(xx,yy)-r)/r
def audit():
 p=load(OUT/'AUX_GATE2B_DESIGN_FREEZE.json');data=rows(OUT/'AUX_T5_JOINT_RESULTS.csv')+rows(OUT/'AUX_T10_JOINT_RESULTS.csv');grid=rows(OUT/'AUX_T5_JOINT_GRID.csv')+rows(OUT/'AUX_T10_JOINT_GRID.csv');checks=[]
 def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
 def near(name,x,y):check(name,math.isclose(float(x),float(y),rel_tol=4e-10,abs_tol=4e-10))
 for name,h in {**p['bindings'],**p['historical_bindings'],**p['external_bindings']}.items():check('hash:'+name,sha(name)==h)
 draws=np.load(OUT/'AUX_FROZEN_DRAWS.npz');z=draws['z'];check('draw_regeneration',np.array_equal(z,np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_trials'],6))))
 for name,part in [('bearing',z[:,:2]),('main_position',z[:,2:4]),('aux_position',z[:,4:6])]:check('saved_draws:'+name,np.array_equal(draws[name],part))
 check('2904cells',len(data)==len(grid)==2904);lookup={key(x):x for x in data};check('all_unique',len(lookup)==2904)
 expected=set()
 for family in FAMILIES:
  for anchor in ['A','B']:
   for lam in [0,.25,.5,.75,1]:
    v=p['vectors'][family][anchor]
    for sc,sd,sb in ([(0,0,0)] if lam==0 else itertools.product([-1,1],repeat=3)):
     for r in range(50,61):
      for side in [-1,1]:
       k=(anchor,family,lam,sc,sd,sb,float(r),side);expected.add(k);x=lookup[k]
       for field,value in [('b_common_deg',lam*v[0]*sc),('b_diff_deg',lam*v[1]*sd),('sigma_pos_m',lam*v[2]),('beta_deg',lam*v[3]*sb)]:near('grid_algebra:'+str(k)+field,x[field],value)
 check('independent_cartesian_grid',set(lookup)==expected and {key(x) for x in grid}==expected)
 for i,g in enumerate(grid,1):
  x=lookup[key(g)];name=x['cell_id']
  for k,val in g.items():check(name+':grid:'+k,x[k]==val)
  a=p['anchors'][x['anchor']];r=float(x['range_km']);b=a['baseline_km'];side=int(x['side']);beta=math.radians(float(x['beta_deg']));s=math.radians(a['sigma_random_deg']);c=math.radians(float(x['b_common_deg']));dd=math.radians(float(x['b_diff_deg']));pos=float(x['sigma_pos_m'])/1000
  # Independently construct actual deployment using cos/sin(90+beta).
  aux=np.array([b*math.cos(math.pi/2+beta),side*b*math.sin(math.pi/2+beta)])
  p1=pos*z[:,2:4]*np.array([1.,side]);p2=aux+pos*z[:,4:6]*np.array([1.,side]);delta=p2-p1
  true2=math.atan2(-aux[1],r-aux[0]);t1=side*s*z[:,0]+c-dd;t2=true2+side*s*z[:,1]+c+dd
  m1=np.tan(t1);m2=np.tan(t2);xx=(delta[:,1]-m2*delta[:,0])/(m1-m2);yy=m1*xx
  hatx=p1[:,0]+xx;haty=p1[:,1]+yy;det=np.sin(t2-t1);parallel=np.abs(det)<=p['parallel_det_tolerance'];behind=(xx*np.cos(t1)+yy*np.sin(t1)<=0)|((xx-delta[:,0])*np.cos(t2)+(yy-delta[:,1])*np.sin(t2)<=0);nonfinite=~np.isfinite(hatx)|~np.isfinite(haty);failed=parallel|behind|nonfinite
  er=np.abs(np.hypot(xx,yy)-r)/r;ep=np.hypot(hatx-r,haty)
  for v in [er,ep]:v[failed]=np.inf
  for metric,v in [('relative_range',er),('position_2d_km',ep)]:
   q=np.sort(v)
   for tag,level in [('median',.5),('p90',.9),('p95',.95),('p99',.99)]:
    value=q[math.ceil(len(v)*level)-1]
    if math.isfinite(value):near(name+':'+metric+tag,x[metric+'_'+tag],value)
    else:check(name+':'+metric+tag,x[metric+'_'+tag]=='INF')
  for field,val in [('failure_rate',failed.mean()),('behind_sensor_rate',behind.mean()),('parallel_rate',parallel.mean()),('nonfinite_rate',nonfinite.mean()),('ill_conditioned_rate',(1/np.maximum(np.abs(det),np.finfo(float).tiny)>p['ill_condition_proxy_limit']).mean())]:near(name+':'+field,x[field],val)
  for field,val in [('baseline_km',b),('sigma_random_deg',a['sigma_random_deg']),('true_aux_x_km',aux[0]),('true_aux_y_km',aux[1]),('crossing_angle_deg',abs(math.degrees(true2)))]:near(name+':'+field,x[field],val)
  n=len(z);zz=1.959963984540054;q=np.sort(er);delta_rank=zz*math.sqrt(n*.95*.05)
  for field,ix in [('p95_rank_ci_low',max(0,math.floor(n*.95-delta_rank)-1)),('p95_rank_ci_high',min(n-1,math.ceil(n*.95+delta_rank)-1))]:near(name+':'+field,x[field],q[ix])
  for tag,t in [('5pct',.05),('10pct',.1)]:
   prop=float((er<=t).mean());den=1+zz*zz/n;center=(prop+zz*zz/(2*n))/den;half=zz*math.sqrt(prop*(1-prop)/n+zz*zz/(4*n*n))/den
   near(name+':fraction'+tag,x['fraction_le'+tag],prop);near(name+':wilson_low'+tag,x['fraction_'+tag+'_wilson_low'],center-half);near(name+':wilson_high'+tag,x['fraction_'+tag+'_wilson_high'],center+half)
  check(name+':valid_trials',flag(x['execution_valid']) and int(x['n_trials'])==n)
  if i%220==0:print(f'independent joint cells {i}/2904',flush=True)

 for x in data:
  if int(x['side'])==1:
   mirror=(x['anchor'],x['family'],float(x['lambda_value']),-int(x['sign_common']),-int(x['sign_diff']),int(x['sign_beta']),float(x['range_km']),-1)
   for field in ['relative_range_p95','relative_range_p99','position_2d_km_p95','failure_rate']:near('mirror:'+x['cell_id']+field,x[field],lookup[mirror][field])
 fronts=rows(OUT/'AUX_JOINT_LAMBDA_FRONTIER.csv')+rows(OUT/'AUX_JOINT_T10_LAMBDA_FRONTIER.csv');check('20frontiers',len(fronts)==20)
 for f in fronts:
  group=[x for x in data if (x['anchor'],x['family'],float(x['lambda_value']))==(f['anchor'],f['family'],float(f['lambda_value']))];worst=max(group,key=lambda x:float(x['relative_range_p95']));v=float(worst['relative_range_p95'])
  check('complete_frontier:'+str(f),len(group)==(22 if float(f['lambda_value'])==0 else 176) and all(flag(x['execution_valid']) for x in group))
  for field,value in [('worst_joint_P95',v),('worst_joint_P99',max(float(x['relative_range_p99']) for x in group)),('worst_position_P95_km',max(float(x['position_2d_km_p95']) for x in group)),('max_failure_rate',max(float(x['failure_rate']) for x in group))]:near('aggregate:'+f['anchor']+f['family']+f['lambda_value']+field,f[field],value)
  for outfield,infield in [('worst_range','range_km'),('worst_sign_common','sign_common'),('worst_sign_diff','sign_diff'),('worst_sign_beta','sign_beta'),('worst_mirror','side')]:near('worst_location:'+str(f)+outfield,f[outfield],worst[infield])
  for tag,threshold in [('T5',.05),('T10',.1)]:check('gate:'+str(f)+tag,flag(f['full_range_'+tag+'_pass'])==(math.isfinite(v) and v<=threshold))
 # New seed: freeze a two-sample DKW control, familywise alpha=.001 over44 unique controls.
 prior=Path('results/R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET');oldz=np.load(prior/'AUX_FROZEN_DRAWS.npz')['z'];oldlookup={(x['anchor'],float(x['range_km']),int(x['side'])):x for x in rows(prior/'AUX_POSITION_RESULTS.csv') if float(x['sigma_pos_m'])==0}
 eps=p['zero_control_dkw_bound'];diagnostics=[]
 for anchor in ['A','B']:
  for r in range(50,61):
   for side in [-1,1]:
    a=p['anchors'][anchor];new=np.sort(zero_errors(z,a,r,side));previous=np.sort(zero_errors(oldz,a,r,side));oldrow=oldlookup[(anchor,float(r),side)]
    for tag,q in [('median',.5),('p90',.9),('p95',.95),('p99',.99)]:near('prior_zero_reproduction:'+str((anchor,r,side))+tag,oldrow['relative_range_'+tag],previous[math.ceil(len(previous)*q)-1])
    points=np.sort(np.concatenate([new,previous]));distance=float(np.max(np.abs(np.searchsorted(new,points,side='right')/len(new)-np.searchsorted(previous,points,side='right')/len(previous))))
    check('zero_sampling_DKW:'+str((anchor,r,side)),distance<=eps)
    for family in FAMILIES:
     row=lookup[(anchor,family,0.,0,0,0,float(r),side)];near('zero_family:'+str((anchor,r,side))+family,row['relative_range_p95'],new[math.ceil(len(new)*.95)-1])
    diagnostics.append(dict(anchor=anchor,range_km=r,side=side,old_p95=float(oldrow['relative_range_p95']),new_p95=float(new[math.ceil(len(new)*.95)-1]),ecdf_sup_distance=distance,frozen_dkw_limit=eps,execution_valid=distance<=eps))
 d=load(OUT/'AUX_GATE2B_DECISION.json');half={family:{a:flag(next(x for x in fronts if x['family']==family and x['anchor']==a and float(x['lambda_value'])==.5)['full_range_'+('T5' if family==FAMILIES[0] else 'T10')+'_pass']) for a in ['A','B']} for family in FAMILIES}
 check('half_point_decision',d['half_point_pass']==half)
 primary=any(half[FAMILIES[0]].values());secondary=any(half[FAMILIES[1]].values());outcome='AUX1_JOINT_STATIC_REQUIREMENT_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED' if primary else 'AUX1_JOINT_SUPPORTS_10PCT_NOT_5PCT' if secondary else 'AUX1_STATIC_GEOMETRY_TOO_FRAGILE_FOR_PROJECT_TARGET'
 check('route_decision',d['Gate2B_decision']==outcome);check('primary_fixed',d['primary_decision']==('AUX1_HALF_T5_JOINT_BUDGET_ESTABLISHED' if primary else 'AUX1_HALF_T5_JOINT_BUDGET_NOT_ESTABLISHED'))
 for m in d['maximum_tested_lambdas']:
  group=[x for x in fronts if x['family']==m['family'] and x['anchor']==m['anchor']]
  for tag in ['T5','T10']:near('maximum:'+str(m)+tag,m['maximum_tested_lambda_'+tag],max(float(x['lambda_value']) for x in group if flag(x['full_range_'+tag+'_pass'])))
 xml=ET.parse(OUT/'AUX_GATE2B_TESTS.xml').getroot();check('tests_pass',len(xml.findall('.//testcase'))==4 and not xml.findall('.//failure') and not xml.findall('.//error'))
 check('stop_and_unknown',d['automatic_further_numerical_work']=='STOP' and d['R4_A1_NEW']=='NOT_OPENED' and d['R4_progress_percent']==0 and d['actual_hardware_capability']=='UNKNOWN' and d['time_synchronization']==d['target_association']=='NOT_NUMERICALLY_VALIDATED')
 start=load(OUT/'AUX_GATE2B_EXECUTION_START.json');check('pushed_before_mc',start['design_sha']==d['design_sha'] and start['remote_verified_before_mc'] and start['policy_sha256']==sha(OUT/'AUX_GATE2B_DESIGN_FREEZE.json'))
 title=subprocess.check_output(['git','log','-1','--format=%s',d['design_sha']],text=True).strip();check('design_commit_title',title=='R4 AUX Gate2B: freeze joint static error budget')
 parent=subprocess.check_output(['git','rev-parse',d['design_sha']+'^'],text=True).strip();check('design_parent',parent==p['parent_sha'])
 for filename in ['results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv']:
  baseline=subprocess.check_output(['git','show',p['parent_sha']+':'+filename]).decode('utf-8').replace('\r\n','\n');check('append_only:'+filename,Path(filename).read_text(encoding='utf-8').startswith(baseline))
 return checks,diagnostics

def present():
 p=load(OUT/'AUX_GATE2B_DESIGN_FREEZE.json');d=load(OUT/'AUX_GATE2B_DECISION.json');front=rows(OUT/'AUX_JOINT_LAMBDA_FRONTIER.csv')+rows(OUT/'AUX_JOINT_T10_LAMBDA_FRONTIER.csv');data=rows(OUT/'AUX_T5_JOINT_RESULTS.csv')+rows(OUT/'AUX_T10_JOINT_RESULTS.csv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figdir=OUT/'figures';figdir.mkdir(exist_ok=True);plt.rcParams.update({'svg.hashsalt':'AUX_GATE2B','font.size':10})
 def save(fig,name):
  fig.savefig(figdir/(name+'.png'),dpi=160,bbox_inches='tight');fig.savefig(figdir/(name+'.svg'),bbox_inches='tight',metadata={'Date':None});plt.close(fig)
  path=figdir/(name+'.svg');path.write_text('\n'.join(s.rstrip() for s in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 for family,name,threshold in [(FAMILIES[0],'FIG1_T5_LADDER',5),(FAMILIES[1],'FIG2_T10_LADDER',10)]:
  fig,ax=plt.subplots(figsize=(8,5))
  for a in ['A','B']:
   group=[x for x in front if x['anchor']==a and x['family']==family];ax.plot([float(x['lambda_value']) for x in group],[100*float(x['worst_joint_P95']) for x in group],'o-',label='Anchor '+a)
  ax.axhline(threshold,color='black',linestyle='--',label=f'{threshold}% Gate');ax.axvline(.5,color='grey',linestyle=':',label='Frozen half point');ax.set(xlabel='Budget scaling lambda',ylabel='Worst joint P95 range error (%)',title=family);ax.legend();ax.grid(alpha=.25);save(fig,name)
 signs=list(itertools.product([-1,1],repeat=3));fig,axes=plt.subplots(2,2,figsize=(12,8))
 for ax,(family,a) in zip(axes.flat,itertools.product(FAMILIES,['A','B'])):
  values=np.array([[100*max(float(x['relative_range_p95']) for x in data if x['family']==family and x['anchor']==a and float(x['lambda_value'])==lam and tuple(int(x[k]) for k in ['sign_common','sign_diff','sign_beta'])==sign) for sign in signs] for lam in [.25,.5,.75,1.]])
  im=ax.imshow(values,aspect='auto',vmin=0,vmax=15)
  for i in range(4):
   for j in range(8):ax.text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=8,color='white')
  ax.set_xticks(range(8),[''.join('+' if v==1 else '-' for v in s) for s in signs]);ax.set_yticks(range(4),[.25,.5,.75,1]);ax.set(xlabel='Signs: common / half-diff / beta',ylabel='lambda',title=family+' / '+a)
 fig.tight_layout();save(fig,'FIG3_SIGN_MAP')
 fig,axes=plt.subplots(1,2,figsize=(11,4))
 for ax,family in zip(axes,FAMILIES):
  group=[next(x for x in front if x['family']==family and x['anchor']==a and float(x['lambda_value'])==.5) for a in ['A','B']];ax.bar(['A','B'],[100*float(x['worst_joint_P95']) for x in group]);ax.axhline(5 if family==FAMILIES[0] else 10,color='black',linestyle='--');ax.set(title=family+' / half point',ylabel='Worst joint P95 (%)')
 fig.tight_layout();save(fig,'FIG4_ANCHOR_COMPARISON')
 lines=['# AUX Gate2B joint static error budget','',d['Gate2B_decision'],'','Primary target remains T5-family lambda=0.50; T10 uses its own distinct vector and ladder.','', '| Family | Anchor | lambda | Worst P95 (%) | T5 | T10 |','|---|---|---:|---:|---|---|']
 for x in front:lines.append(f"| {x['family']} | {x['anchor']} | {x['lambda_value']} | {100*float(x['worst_joint_P95']):.6f} | {x['full_range_T5_pass']} | {x['full_range_T10_pass']} |")
 lines+=['','## Interpretation and limitations','','True MAIN=(0,0), AUX=B[-sin(beta),side*cos(beta)], target=(R,0). Bearings use actual geometry, independent Gaussian random terms plus fixed b1=common-half_diff, b2=common+half_diff. Estimated node positions have independent per-axis Gaussian offsets. Deployment and navigation are separate factors, both active. Random error is not scaled by lambda. Range=norm(target_hat-MAIN_est), compared with true R; absolute2D position error is reported separately.','','Each cell has30000 trials; one new frozen PCG64 seed and common six-column standard-normal table across both families. Mirrors reflect bearing noise and navigation y. Every deterministic bias/deployment sign corner is retained. Shared draws/mirrors are not independent repetitions. All parallel, nonfinite and behind-sensor failures remain infinite errors in nearest-rank quantiles. No outlier rejection. The max over11 ranges, both mirrors and8 signed corners is the joint Gate metric, not pooled statistics.','','Lambda0 uses one sign identity per family; independent reproduction of Gate2A saved controls and a preregistered two-sample DKW criterion (familywise alpha=.001 over44 unique controls) check the new seed. Failure means EXECUTION_INVALID. This sampling diagnostic is not hardware assurance.','','Numerical joint requirements apply to frozen static model, tested signed corners, scaling values and11 discrete50:1:60km ranges. No continuous hyperbox, intermediate magnitude/range, correlated-error, dynamic/time-skew or association-failure certificate. A maximum tested lambda is a discrete observed frontier, not an interpolated threshold. The half-point remains primary even when smaller lambda passes.','','Actual auxiliary array, DOA, navigation/attitude, calibration and clocks remain UNKNOWN. Time/association NOT_NUMERICALLY_VALIDATED. No automatic Gate2C or further numerical architecture work. Next action is client requirement confirmation or a research-lead architecture decision. R4=0%; R4-A1-NEW NOT_OPENED; depth closed. Historical Gate0/1/2A artifacts preserved; client table supersedes the operational presentation without rewriting frozen Gate2A files.','','```json',json.dumps(d,indent=2),'```']
 report='\n'.join(lines)+'\n'
 for name in ['AUX_GATE2B_REPORT.md','GPT_SYNC.md']:(OUT/name).write_text(report,encoding='utf-8')
 client=['# Auxiliary measurement requirement confirmation','', 'Gate2A maxima are SINGLE FACTOR ONLY. Joint requirements below come from Gate2B simultaneous error models; the columns cannot be interchanged. Client values remain UNKNOWN.','', '| Quantity | Single-factor tested max | Jointly validated requirement | Client value | Status |','|---|---|---|---|---|']
 for family in FAMILIES:
  tag='T5' if family==FAMILIES[0] else 'T10'
  for a in ['A','B']:
   m=next(x for x in d['maximum_tested_lambdas'] if x['family']==family and x['anchor']==a);lam=m['maximum_tested_lambda_'+tag];v=p['vectors'][family][a];anchor=p['anchors'][a]
   for quantity,maximum,joint in [('baseline km',anchor['baseline_km'],anchor['baseline_km']),('per-node random bearing sigma deg',anchor['sigma_random_deg'],anchor['sigma_random_deg']),('common bias magnitude deg',v[0],lam*v[0]),('differential HALF-bias magnitude deg',v[1],lam*v[1]),('per-node per-axis position sigma m',v[2],lam*v[2]),('known actual deployment |beta| deg',v[3],lam*v[3])]:client.append(f'| {tag} Anchor {a}: {quantity} | {maximum} | {joint} at tested lambda={lam} | UNKNOWN | CLIENT_CONFIRMATION_REQUIRED |')
 client+=['','Each joint vector must be considered as a complete vector with its anchor random sigma. Signed corners tested; no continuous interior guarantee. Differential HALF-bias is (b2-b1)/2, total differential twice this value. Position sigma is Cartesian per-axis1sigma, not radial RMS. Deployment changes true geometry while navigation perturbs estimator positions. Gate2A20deg is a tested grid edge, not a proven tolerance limit.','','Primary half point requirements:']
 for a in ['A','B']:
  client.append(f"- Anchor {a}: T5 half point {'PASS' if d['half_point_pass'][FAMILIES[0]][a] else 'FAIL'}; common magnitude0.05deg, half-diff0.025deg, position25m/axis, deployment10deg, frozen random sigma{p['anchors'][a]['sigma_random_deg']}deg.")
 client+=['','All synchronization, target association/front-back resolution and array/attitude calibration requirements require client evidence. Synchronization/association NOT_NUMERICALLY_VALIDATED, no timing tolerance is claimed. Separate hardware measured random and systematic errors, covariance, node/aperture configuration, directed-bearing accuracy statistics, navigation position covariance, heading calibration, deployable baseline, position updates and clock provenance are required.','','No actual equipment admission; R4-A1-NEW NOT_OPENED; R4=0%. STOP numerical architecture work and obtain client hardware confirmation.']
 (OUT/'AUX_CLIENT_REQUIREMENTS_UPDATED.md').write_text('\n'.join(client)+'\n',encoding='utf-8')
 with Path('results/R4_MASTER/R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## AUX Gate2B execution\n\n'+d['Gate2B_decision']+'; PENDING_RESEARCH_LEAD_AUDIT. Joint static ladders complete. [Report](../R4_AUX_GATE2B_JOINT_ERROR_BUDGET/AUX_GATE2B_REPORT.md). Hardware UNKNOWN; time/association not validated; R4=0%; A1-NEW not opened. STOP numerical architecture work; client confirmation next.\n')
 with Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow(['AUX-GATE2B','../R4_AUX_GATE2B_JOINT_ERROR_BUDGET/AUX_GATE2B_DECISION.json',d['Gate2B_decision'],'Two frozen joint static families; all signed corners;11 ranges; double mirrors','PENDING_RESEARCH_LEAD_AUDIT; HARDWARE_UNKNOWN; STOP; NO_PROGRESS_CREDIT',0])
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--finalize',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 if args.finalize:
  assert not (OUT/'VALIDATION.json').exists(),'One finalize only'
  checks,zero=audit();present();write(OUT/'AUX_INDEPENDENT_CHECKS.csv',checks);write(OUT/'AUX_ZERO_CONTROL_VALIDATION.csv',zero)
  dump(OUT/'VALIDATION.json',dict(status='PASS',independent_checks=len(checks),failed=0,primary_runs=1,tests_passed=4,registered_cells=2904,trials_per_cell=30000,zero_controls=44,new_acoustic_propagation=0,new_depth_scores=0))
  paths=[x for x in OUT.rglob('*') if x.is_file()]+[Path('results/R4_MASTER/R4_PLAN.md'),Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv')]
  dump(OUT/'AUX_RESULTS_FREEZE.json',dict(bindings={x.as_posix():sha(x) for x in sorted(paths)}));print(json.dumps(dict(finalized='PASS',checks=len(checks))),flush=True)
 elif args.verify:
  for name,h in load(OUT/'AUX_RESULTS_FREEZE.json')['bindings'].items():assert sha(name)==h,name
  checks,zero=audit();print(json.dumps(dict(sealed='PASS',checks=len(checks),failed=0)),flush=True)
 else:ap.error('--finalize or --verify required')
