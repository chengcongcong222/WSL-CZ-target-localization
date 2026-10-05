"""Independent saved-scene intervention, metrics, Jacobian and summary audit."""
from pathlib import Path
import argparse,ast,csv,json,hashlib,math,subprocess,xml.etree.ElementTree as ET
import numpy as np
OUT=Path('results/R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC');OLD=Path('results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC');VARIANTS=['D0','D1','D2','D3','D4'];METRICS=['range_error','bearing_error_deg','speed_error','heading_error_deg']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,data):
 with Path(p).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def flag(x):return x is True or x=='True'
def wrap(x):return np.arctan2(np.sin(x),np.cos(x))
def audit():
 p=load(OUT/'SPEED_DIAGNOSTIC_POLICY.json');checks=[]
 def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
 def near(name,a,b,rtol=3e-7,atol=2e-7):check(name,math.isclose(float(a),float(b),rel_tol=rtol,abs_tol=atol))
 for name,h in {**p['bindings'],**p['historical_bindings']}.items():check('hash:'+name,sha(name)==h)
 check('D0_byte_identity',(OUT/'D0_BASELINE_IDENTITY.csv').read_bytes()==(OLD/'A1_NEW_RUN_RESULTS.csv').read_bytes())
 old=rows(OLD/'A1_NEW_RUN_RESULTS.csv');panel=rows(OLD/'A1_NEW_TRUTH_PANEL.csv');oldpolicy=load(OLD/'A1_NEW_DESIGN_FREEZE.json');data=rows(OUT/'SPEED_VARIANT_RUN_RESULTS.csv');diag=rows(OUT/'BIAS_AWARE_ESTIMATOR_DIAGNOSTICS.csv');jdata=rows(OUT/'BIAS_AWARE_JACOBIAN_RUN_RESULTS.csv')
 def key(x):return x['anchor'],int(x['draw_index']),x['variant']
 lookup={key(x):x for x in data};di={(x['anchor'],int(x['draw_index'])):x for x in diag};ji={(x['anchor'],int(x['draw_index']),x['segment']):x for x in jdata}
 check('60000_variant_runs',len(lookup)==len(data)==60000);check('12000_bias_diagnostics',len(di)==len(diag)==12000);check('36000_Jacobians',len(ji)==len(jdata)==36000)
 tree=ast.parse(Path('r4_speed_bias_estimator.py').read_text(encoding='utf-8'));fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='estimate_bias_aware');check('candidate_truth_free_API',[x.arg for x in fn.args.args]==['times','node_positions','bearings','sigma_deg','common_bound_deg','diff_bound_deg']);check('candidate_no_truth_symbols',not any(isinstance(x,ast.Name) and ('truth' in x.id.lower() or x.id in ['scene','evaluate']) for x in ast.walk(tree)))
 runner=ast.parse(Path('r4_speed_diagnostic.py').read_text(encoding='utf-8'));check('no_random_draw_calls',not any(isinstance(x,ast.Attribute) and x.attr in ['random','standard_normal','PCG64','default_rng','Generator'] for x in ast.walk(runner)))
 t=np.arange(121)*10;main=np.array([[2*tt,0] if tt<=600 else [1200+2*(tt-600)*math.cos(math.pi/12),2*(tt-600)*math.sin(math.pi/12)] for tt in t]);recomputed=[];fd=[]
 for aname in ['A','B']:
  a=oldpolicy['anchors'][aname];draw=np.load(OLD/f'A1_NEW_DRAWS_{aname}.npz');base=[x for x in old if x['anchor']==aname]
  for ix,b in enumerate(base):
   truth=panel[ix//500];side=int(b['mirror']);u=draw['uniforms'][ix];z=draw['normals'][ix];c=(2*u[0]-1)*a['common_deg'];diff=(2*u[1]-1)*a['half_diff_deg'];beta=(2*u[2]-1)*a['beta_deg'];angle=math.radians(90+beta);offset=a['baseline_m']*np.array([math.cos(angle),side*math.sin(angle)]);nodes=np.stack([main,main+offset],axis=1);navnodes=nodes+a['nav_sigma_m']*z[:,2:].reshape(121,2,2)
   r=float(truth['r0_km'])*1000;theta=math.radians(float(truth['theta0_deg']));v=float(truth['v_mps']);psi=math.radians(float(truth['psi_deg']));target=np.array([r*math.cos(theta),r*math.sin(theta)])+t[:,None]*np.array([v*math.cos(psi),v*math.sin(psi)]);delta=target[:,None,:]-nodes;raw=np.arctan2(delta[:,:,1],delta[:,:,0])+np.deg2rad([c-diff,c+diff])+math.radians(a['sigma_deg'])*z[:,:2]
   for variant in VARIANTS:
    x=lookup[(aname,ix,variant)];name=f'{aname}:{ix}:{variant}';check(name+'metadata',x['case_id']==b['case_id'] and int(x['realization'])==int(b['realization']) and int(x['mirror'])==side)
    scope='IDENTITY_CONTROL' if variant=='D0' else 'CANDIDATE_DEVELOPMENT_ONLY' if variant=='D4' else 'ORACLE_DIAGNOSTIC_ONLY';check(name+'scope',x['scope']==scope)
    if variant=='D0':
     for m in METRICS+['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps','r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg']:near(name+'identity'+m,x[m],b[m])
     check(name+'identity_failure',flag(x['failed'])==flag(b['failed']))
    if flag(x['failed']):
     for m in METRICS+['absolute_speed_error_mps']:check(name+'infinite'+m,float(x[m])==float('inf'))
     recomputed.append(dict(anchor=aname,case_id=x['case_id'],variant=variant,failed=True,**{m:float('inf') for m in METRICS+['absolute_speed_error_mps']},common_bias_error_deg=float('inf') if variant=='D4' else 'NA',diff_HALF_bias_error_deg=float('inf') if variant=='D4' else 'NA'));continue
    estimator_nodes=navnodes if variant in ['D0','D1','D4'] else nodes;measurements=raw-np.deg2rad([c-diff,c+diff]) if variant in ['D1','D3'] else raw;state=np.array([float(x[k]) for k in ['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps']]);ref=state[:2]-estimator_nodes[0,0];rh=np.hypot(*ref);thh=math.atan2(ref[1],ref[0]);vh=np.hypot(*state[2:]);ph=math.atan2(state[3],state[2]);errors=dict(range_error=abs(rh-r)/r,bearing_error_deg=abs(math.degrees(float(wrap(thh-theta)))),speed_error=abs(vh-v)/v,heading_error_deg=abs(math.degrees(float(wrap(ph-psi)))),absolute_speed_error_mps=abs(vh-v))
    for m,val in errors.items():near(name+'coldmetric'+m,x[m],val)
    if variant=='D4':
     d=di[(aname,ix)];check(name+'biaszero_init',float(d['initial_common_bias_deg'])==float(d['initial_diff_bias_deg'])==0.);check(name+'solver_metadata',flag(d['solver_success']) and int(d['nfev'])<=300)
     # Cold slope E1, no historical-output warm start.
     m1=np.tan(raw[:,0]);m2=np.tan(raw[:,1]);dd=navnodes[:,1]-navnodes[:,0]
     with np.errstate(divide='ignore',invalid='ignore'):
      xx=(dd[:,1]-m2*dd[:,0])/(m1-m2);yy=m1*xx;points=navnodes[:,0]+np.column_stack([xx,yy]);det=np.sin(raw[:,1]-raw[:,0]);l1=xx*np.cos(raw[:,0])+yy*np.sin(raw[:,0]);l2=(xx-dd[:,0])*np.cos(raw[:,1])+(yy-dd[:,1])*np.sin(raw[:,1])
     valid=(np.abs(det)>1e-12)&(l1>0)&(l2>0)&np.isfinite(points).all(axis=1);T=np.column_stack([np.ones(valid.sum()),t[valid]/1200]);points=points[valid];coef=np.linalg.lstsq(T,points,rcond=None)[0]
     for _ in range(5):
      norms=np.linalg.norm(points-T@coef,axis=1);weights=np.minimum(1,1.5*max(10,float(np.median(norms)))/np.maximum(norms,1e-15));coef=np.linalg.lstsq(T*np.sqrt(weights)[:,None],points*np.sqrt(weights)[:,None],rcond=None)[0]
     ini=np.r_[coef[0],coef[1]/1200]
     for field,val in zip(['initial_x_m','initial_y_m','initial_vx_mps','initial_vy_mps'],ini):near(name+'E1'+field,d[field],val)
     ch=float(x['common_hat_deg']);dh=float(x['diff_HALF_hat_deg']);check(name+'bias_support',abs(ch)<=a['common_deg']+1e-12 and abs(dh)<=a['half_diff_deg']+1e-12)
     near(name+'common_error',x['common_bias_error_deg'],abs(ch-c));near(name+'diff_error',x['diff_HALF_bias_error_deg'],abs(dh-diff));errors.update(common_bias_error_deg=abs(ch-c),diff_HALF_bias_error_deg=abs(dh-diff))
     check(name+'boundflagcommon',flag(x['common_at_bound'])==(abs(ch/a['common_deg'])>=1-1e-6));check(name+'boundflagdiff',flag(x['diff_at_bound'])==(abs(dh/a['half_diff_deg'])>=1-1e-6))
     delta_hat=(state[:2]+t[:,None]*state[2:])[:,None,:]-navnodes;pred=np.arctan2(delta_hat[:,:,1],delta_hat[:,:,0])+np.deg2rad([ch-dh,ch+dh]);rr=wrap(pred-raw)/math.radians(a['sigma_deg']);near(name+'cost',d['cost'],2.25*np.sum(np.sqrt(1+(rr/1.5)**2)-1),rtol=1e-6,atol=1e-6);near(name+'rms',d['residual_RMS_deg'],np.sqrt(np.mean(rr**2))*a['sigma_deg'])
     rho=np.sum(delta_hat**2,axis=2);gx=-delta_hat[:,:,1]/rho;gy=delta_hat[:,:,0]/rho;J=np.stack([gx*50000,gy*50000,gx*t[:,None]*2,gy*t[:,None]*2],axis=-1).reshape(242,4)/math.radians(a['sigma_deg']);J=np.column_stack([J,np.full(242,a['common_deg']/a['sigma_deg']),np.tile([-a['half_diff_deg']/a['sigma_deg'],a['half_diff_deg']/a['sigma_deg']],121)])
     for segment,mask in [('STRAIGHT',t<=600),('POST_TURN',t>600),('FULL',np.ones(121,dtype=bool))]:
      saved=ji[(aname,ix,segment)];sv=np.linalg.svd(J.reshape(121,2,6)[mask].reshape(-1,6),compute_uv=False);rank=int(np.sum(sv/sv[0]>1e-10));check(name+'rank'+segment,int(saved['numerical_rank'])==rank and flag(saved['full_rank'])==(rank==6))
      for j,val in enumerate(sv):near(name+'SVD'+segment+str(j),saved[f'singular_value_{j+1}'],val,rtol=1e-6,atol=1e-9)
      near(name+'ratio'+segment,saved['smallest_largest_ratio'],sv[-1]/sv[0],rtol=1e-6,atol=1e-13)
     if ix%500 in [0,499]:
      q=np.r_[state/np.array([50000,50000,2,2]),ch/a['common_deg'],dh/a['half_diff_deg']]
      def fun(q):
       ss=q[:4]*np.array([50000,50000,2,2]);dlt=(ss[:2]+t[:,None]*ss[2:])[:,None,:]-navnodes;pp=np.arctan2(dlt[:,:,1],dlt[:,:,0])+np.deg2rad(q[4]*a['common_deg']+q[5]*a['half_diff_deg']*np.array([-1,1]));return (wrap(pp-raw)/math.radians(a['sigma_deg'])).ravel()
      h=1e-6;finite=np.column_stack([(fun(q+np.eye(6)[j]*h)-fun(q-np.eye(6)[j]*h))/(2*h) for j in range(6)]);ok=np.allclose(finite,J,rtol=1e-5,atol=1e-5);check(name+'independent_FD_Jacobian',ok);fd.append(dict(anchor=aname,case_id=x['case_id'],realization=int(x['realization']),max_absolute_Jacobian_difference=float(np.max(np.abs(finite-J))),passed=bool(ok)))
    else:errors.update(common_bias_error_deg='NA',diff_HALF_bias_error_deg='NA')
    recomputed.append(dict(anchor=aname,case_id=x['case_id'],variant=variant,failed=False,**errors))
   if (ix+1)%500==0:print(f'Independent saved-scene audit {aname} {ix+1}/6000',flush=True)
 summaries=rows(OUT/'SPEED_VARIANT_CASE_SUMMARY.csv');check('120_case_variant_summaries',len(summaries)==120)
 for row in summaries:
  group=[x for x in recomputed if (x['anchor'],x['case_id'],x['variant'])==(row['anchor'],row['case_id'],row['variant'])];name=row['anchor']+row['case_id']+row['variant'];check(name+'500samples',len(group)==int(row['n_runs'])==500);near(name+'failure',row['failure_rate'],sum(x['failed'] for x in group)/500)
  for m in METRICS+['absolute_speed_error_mps']:
   vals=sorted(x[m] for x in group)
   for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:near(name+m+tag,row[m+'_'+tag],vals[math.ceil(500*q)-1])
  for field,m in [('speed_rel_P95','speed_error'),('speed_abs_P95_mps','absolute_speed_error_mps'),('range_P95','range_error'),('bearing_P95','bearing_error_deg'),('heading_P95','heading_error_deg')]:near(name+'alias'+field,row[field],row[m+'_P95'])
  check(name+'speed_PROJECT',flag(row['PROJECT_speed_pass'])==(float(row['speed_rel_P95'])<=.1));check(name+'speed_STRONG',flag(row['STRONG_speed_pass'])==(float(row['speed_rel_P95'])<=.05));check(name+'all_PROJECT',flag(row['PROJECT_all_metrics_pass'])==all(float(row[m+'_P95'])<=lim for m,lim in zip(METRICS,[.1,1,.1,5])));check(name+'failure_gate',flag(row['route_failure_rate_pass'])==(float(row['failure_rate'])<=.01))
  if row['variant']=='D4':
   for m in ['common_bias_error_deg','diff_HALF_bias_error_deg']:near(name+m,row[m+'_P95'],sorted(x[m] for x in group)[474])
   valid=[x for x in data if (x['anchor'],x['case_id'],x['variant'])==(row['anchor'],row['case_id'],'D4') and not flag(x['failed'])]
   for field,source in [('common_at_bound_fraction_valid_fits','common_at_bound'),('diff_at_bound_fraction_valid_fits','diff_at_bound')]:
    if valid:near(name+field,row[field],sum(flag(x[source]) for x in valid)/len(valid))
    else:check(name+field,row[field]=='NA')
 for row in rows(OUT/'SPEED_ATTRIBUTION_SUMMARY.csv'):
  pair=[x for x in summaries if (x['anchor'],x['case_id'])==(row['anchor'],row['case_id'])];b=float(next(x for x in pair if x['variant']=='D0')['speed_rel_P95']);v=float(next(x for x in pair if x['variant']==row['variant'])['speed_rel_P95']);near('attribution:'+str(row),row['P95_reduction_relative_to_D0'],(b-v)/b)
 js=rows(OUT/'BIAS_AWARE_JACOBIAN_SUMMARY.csv');check('72Jacobian_summaries',len(js)==72)
 for row in js:
  group=[x for x in jdata if (x['anchor'],x['case_id'],x['segment'])==(row['anchor'],row['case_id'],row['segment'])];near('rank_rate:'+str(row),row['full_rank_rate'],sum(flag(x['full_rank']) for x in group)/500);vals=sorted(float(x['smallest_largest_ratio']) for x in group)
  for field,value in [('ratio_min',vals[0]),('ratio_median',vals[249]),('ratio_P95',vals[474])]:near('jac_summary:'+str(row)+field,row[field],value)
 d=load(OUT/'SPEED_DIAGNOSTIC_DECISION.json');A={v:[x for x in summaries if x['anchor']=='A' and x['variant']==v] for v in VARIANTS};passes={v:sum(flag(x['PROJECT_speed_pass']) for x in group) for v,group in A.items()};route=all(flag(x['PROJECT_all_metrics_pass']) and float(x['failure_rate'])<=.01 for x in A['D4']);verdict='BIAS_AWARE_SPEED_ROUTE_WORTH_FRESH_VALIDATION' if route else 'SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO' if passes['D3']<12 else 'BIAS_INFORMATION_PRESENT_BUT_ESTIMATOR_NOT_YET_ESTABLISHED' if passes['D1']==12 else 'SPEED_ROUTE_NOT_ESTABLISHED_UNDER_CURRENT_SAVED_SCENE_DIAGNOSTIC';check('frozen_route_decision',d['speed_route_decision']==verdict and d['primary_PROJECT_speed_cases']==passes and d['D4_route_admitted']==route);check('no_credit_no_new_draws',all(d[k]==0 for k in ['new_random_seeds','new_truth_cases','new_noise_draws','new_Monte_Carlo_realizations','new_application_performance_claim','R4_A1_percent','R4_percent']) and not d['automatic_fresh_validation'] and d['depth']=='NOT_OPENED')
 start=load(OUT/'SPEED_EXECUTION_START.json');check('pushed_before_resolves',start['remote_verified_before_resolving'] and start['policy_SHA']==d['policy_SHA'] and start['policy_sha256']==sha(OUT/'SPEED_DIAGNOSTIC_POLICY.json'));check('policy_parent',subprocess.check_output(['git','rev-parse',d['policy_SHA']+'^'],text=True).strip()==p['parent_SHA'])
 xml=ET.parse(OUT/'SPEED_TESTS.xml').getroot();check('4tests_pass',len(xml.findall('.//testcase'))==4 and not xml.findall('.//failure') and not xml.findall('.//error'))
 for name in ['results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv']:
  before=subprocess.check_output(['git','show',p['parent_SHA']+':'+name]).decode('utf-8').replace('\r\n','\n');check('append_only:'+name,Path(name).read_text(encoding='utf-8').startswith(before))
 return checks,fd

def present():
 d=load(OUT/'SPEED_DIAGNOSTIC_DECISION.json');s=rows(OUT/'SPEED_VARIANT_CASE_SUMMARY.csv');data=rows(OUT/'SPEED_VARIANT_RUN_RESULTS.csv');js=rows(OUT/'BIAS_AWARE_JACOBIAN_SUMMARY.csv');jdata=rows(OUT/'BIAS_AWARE_JACOBIAN_RUN_RESULTS.csv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figs=OUT/'figures';figs.mkdir(exist_ok=True);plt.rcParams['svg.hashsalt']='SPEED_DIAG'
 def save(fig,name):
  fig.savefig(figs/(name+'.png'),dpi=160,bbox_inches='tight');fig.savefig(figs/(name+'.svg'),bbox_inches='tight',metadata={'Date':None});plt.close(fig);p=figs/(name+'.svg');p.write_text('\n'.join(x.rstrip() for x in p.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 fig,axes=plt.subplots(1,2,figsize=(13,5))
 for ax,a in zip(axes,['A','B']):
  for v in VARIANTS:
   g=[x for x in s if x['anchor']==a and x['variant']==v];ax.plot(range(1,13),[100*float(x['speed_rel_P95']) for x in g],'o-',label=v)
  ax.axhline(10,color='black',linestyle='--');ax.set(title='Anchor '+a+' saved-scene development',xlabel='Truth case',ylabel='Speed P95 (%)');ax.legend();ax.grid(alpha=.2)
 fig.tight_layout();save(fig,'FIG1_SPEED_VARIANTS')
 fig,axes=plt.subplots(1,2,figsize=(12,5))
 for ax,field,scale,label in zip(axes,['speed_rel_P95','speed_abs_P95_mps'],[100,1],['Relative speed P95 (%)','Absolute speed P95 (m/s)']):
  for v in VARIANTS:
   g=[x for x in s if x['anchor']=='A' and x['variant']==v];ax.plot(range(1,13),[scale*float(x[field]) for x in g],'o-',label=v)
  ax.set(xlabel='Anchor A case',ylabel=label);ax.legend();ax.grid(alpha=.2)
 fig.tight_layout();save(fig,'FIG2_RELATIVE_ABSOLUTE_SPEED')
 fig,axes=plt.subplots(1,2,figsize=(11,5))
 for ax,truth,estimated in zip(axes,['common_true_deg','diff_HALF_true_deg'],['common_hat_deg','diff_HALF_hat_deg']):
  for a in ['A','B']:
   g=[x for x in data if x['anchor']==a and x['variant']=='D4' and not flag(x['failed'])];ax.scatter([float(x[truth]) for x in g],[float(x[estimated]) for x in g],s=2,alpha=.12,label=a)
  ax.set(xlabel=truth,ylabel=estimated);ax.legend()
 fig.tight_layout();save(fig,'FIG3_BIAS_ESTIMATION')
 fig,axes=plt.subplots(1,2,figsize=(12,5))
 for ax,a in zip(axes,['A','B']):
  for seg in ['STRAIGHT','POST_TURN','FULL']:
   g=[x for x in js if x['anchor']==a and x['segment']==seg];ax.semilogy(range(1,13),[max(1e-18,float(x['ratio_median'])) for x in g],'o-',label=seg)
  ax.axhline(1e-10,color='black',linestyle='--',label='Rank threshold');ax.set(title='Anchor '+a,xlabel='Truth case',ylabel='Median smallest/largest singular ratio');ax.legend();ax.grid(alpha=.2)
 fig.tight_layout();save(fig,'FIG4_JACOBIAN_RATIO')
 lines=['# A1-new speed bottleneck causal diagnostic','',d['speed_route_decision'],'','DEVELOPMENT / CAUSAL DIAGNOSTIC. All12000savedrealizations reused, no new seeds/truth/draws/performanceclaim. D0 is byte-identical saved output; D1 removes generatedbias only; D2 uses true positions only; D3 does both; D4 is the sole observation-only candidate. D1-D3 are ORACLE_DIAGNOSTIC_ONLY and cannot be algorithms or project claims.','', '| Anchor | Branch | Worst speed P95 % | Absolute speed P95 m/s (case min..max) | PROJECT speed cases | Worst range P95 % | Bearing P95 deg | Heading P95 deg | Max case failure rate |','|---|---|---:|---|---:|---:|---:|---:|---:|']
 for a in ['A','B']:
  for v in VARIANTS:
   g=[x for x in s if x['anchor']==a and x['variant']==v];absolute=[float(x['speed_abs_P95_mps']) for x in g];lines.append(f"| {a} | {v} | {100*max(float(x['speed_rel_P95']) for x in g):.6f} | {min(absolute):.6f}..{max(absolute):.6f} | {sum(flag(x['PROJECT_speed_pass']) for x in g)}/12 | {100*max(float(x['range_P95']) for x in g):.6f} | {max(float(x['bearing_P95']) for x in g):.6f} | {max(float(x['heading_P95']) for x in g):.6f} | {max(float(x['failure_rate']) for x in g):.6f} |")
 lines+=['','## Attribution scope','','The generator includes fixed realization-wise common/HALF-differential biases while historical4stateestimator has no bias nuisance: MISSPECIFIED_4_STATE_MODEL_UNDER_STATIC_BEARING_BIAS. The mismatch was frozen as a hypothesis, not a predeclared cause. Interventions compare exactly paired savednoise/navigation/deployment; P95 reductions are nonadditive and not a variance decomposition or unique-cause proof. Failure counts remain infinite errors, no sample rejection/replacement. Absolute speed error and relative error share fixed casev, their P95 relation is exact within a case; compare acrosscases without asserting a universal floor.','','D4 uses E1 from original bearings/noisy positions and starts both biases at0. Bounded unregularizedTRF, analytic6stateJacobian, soft_l1f_scale1.5, scales50km/2m/s/designbiasbound, max300evaluations. No truthbias/targetinitialization, historicalrunwarmstart, multistart or added penalty. Bounds are design assumptions, not calibrated hardware knowledge. D4 bias estimates are nuisance diagnostics, not project outputs.','','The full-rank threshold is scaled s/smax>1e-10, evaluated in61straightepochs/60post-turnepochs/full121. Full-rank is not an accuracy certificate; extremely small ratios and bound-hitting rates may expose weak bias/state separation. Turning comparison reuses the same15degturn, no excitation sweep. Fullwindow rank improvement alone does not prove PROJECTprecision.','','D3 failures support the requested SPEED_INFORMATION_INSUFFICIENT_AT_CURRENT_1200S_APPLICATION_SCENARIO label only under this frozen geometry/randomnoise and current4stateestimator; they do not establish a lower bound for allpossibleestimators or physical nonidentifiability. D1passwould demonstrate calibrationrescue onthispanel, not a realizableoraclealgorithm. D4passonlyadmitsfreshvalidation, noA1credit.','','New MC=0; R4-A1=0%; R4=0%; depthNOT_OPENED; noacoustic/depthscore. Existing A1submetrics preserved:range/bearingSTRONG,headingPROJECT, speedbottleneck. Stop after executionpush; noFIX2, nofreshconfirmation automatically.','', '```json',json.dumps(d,indent=2),'```']
 for a in ['A','B']:
  group=[x for x in jdata if x['anchor']==a];lines+=['', '## Jacobian '+a,'']
  for seg in ['STRAIGHT','POST_TURN','FULL']:
   g=[x for x in group if x['segment']==seg];ratios=sorted(float(x['smallest_largest_ratio']) for x in g);lines.append(f'{seg}: full-rank rate={sum(flag(x["full_rank"]) for x in g)/6000:.6f}; ratio min/median/max={ratios[0]:.6e}/{ratios[2999]:.6e}/{ratios[-1]:.6e}.')
 text='\n'.join(lines)+'\n'
 for name in ['SPEED_DIAGNOSTIC_REPORT.md','GPT_SYNC.md']:(OUT/name).write_text(text,encoding='utf-8')
 progress={k:d[k] for k in ['stage','nature','range','bearing','heading','speed','R4_A1_percent','R4_percent','speed_route_decision','next_recommended_stage','automatic_fresh_validation','depth']};progress.update(project_stage='PROJECT_APPLICATION_PRE_RESEARCH',A1_new_audit='R4_A1_NEW_BASELINE_NOT_ESTABLISHED_ACCEPTED',diagnostic_audit='PENDING_RESEARCH_LEAD_AUDIT',hardware_UNKNOWN_blocks_research=False);dump('results/R4_MASTER/R4_SPEED_DIAGNOSTIC_PROGRESS.json',progress)
 with Path('results/R4_MASTER/R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## Speed bottleneck diagnostic execution\n\n'+d['speed_route_decision']+'; DEVELOPMENT_ONLY, PENDING_RESEARCH_LEAD_AUDIT. [Report](../R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC/SPEED_DIAGNOSTIC_REPORT.md). [Active speed status](R4_SPEED_DIAGNOSTIC_PROGRESS.json). AcceptedA1range/bearingSTRONG,headingPROJECT; speedunresolved. Reused12000scenes; newMC0; R4-A1/R4=0%. Depthnotopened. PushthenSTOP; freshonlyifseparatelyauthorized.\n')
 with Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow(['A1-NEW-SPEED-DIAG','../R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC/SPEED_DIAGNOSTIC_DECISION.json',d['speed_route_decision'],'Pairedsavedscenecausaldiagnostic;D1-D3oracle;D4development;newdraws0','PENDING_RESEARCH_LEAD_AUDIT;NO_SCIENTIFIC_CREDIT',0])
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--finalize',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 if args.finalize:
  assert not (OUT/'VALIDATION.json').exists(),'one finalize'
  checks,fd=audit();present();write(OUT/'SPEED_INDEPENDENT_CHECKS.csv',checks)
  if fd:write(OUT/'BIAS_AWARE_FD_JACOBIAN_CHECKS.csv',fd)
  dump(OUT/'VALIDATION.json',dict(status='PASS',checks=len(checks),failed=0,tests_passed=4,saved_realizations=12000,variant_records=60000,cold_metrics=60000,Jacobian_records=36000,new_Monte_Carlo_realizations=0,new_noise_draws=0,new_acoustic_propagation=0,new_depth_score=0))
  paths=[x for x in OUT.rglob('*') if x.is_file()]+[Path('results/R4_MASTER/R4_PLAN.md'),Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'),Path('results/R4_MASTER/R4_SPEED_DIAGNOSTIC_PROGRESS.json')];dump(OUT/'SPEED_RESULTS_FREEZE.json',dict(bindings={x.as_posix():sha(x) for x in sorted(paths)}));print(json.dumps(dict(finalized='PASS',checks=len(checks))),flush=True)
 elif args.verify:
  for name,h in load(OUT/'SPEED_RESULTS_FREEZE.json')['bindings'].items():assert sha(name)==h,name
  checks,fd=audit();print(json.dumps(dict(sealed='PASS',checks=len(checks),failed=0)),flush=True)
 else:ap.error('--finalize or --verify required')
