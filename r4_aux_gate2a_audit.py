"""Independent Gate2A slope audit and saved-result requirements, no acoustic imports."""
from pathlib import Path
import argparse,csv,hashlib,json,math,subprocess,xml.etree.ElementTree as ET
import numpy as np
OUT=Path('results/R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET')
KEYS=['anchor','family','b_common_deg','b_diff_deg','sigma_pos_m','beta_deg']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def flag(x):return x is True or x=='True'
def key(x):return (x['anchor'],x['family'],*[float(x[k]) for k in ['b_common_deg','b_diff_deg','sigma_pos_m','beta_deg','range_km']],int(x['side']))
def audit():
 p=load(OUT/'AUX_GATE2A_DESIGN_FREEZE.json');data=sum([rows(OUT/name) for name in ['AUX_BIAS_RESULTS.csv','AUX_POSITION_RESULTS.csv','AUX_DEPLOYMENT_RESULTS.csv']],[]);grid=rows(OUT/'AUX_GATE2A_ERROR_GRIDS.csv');checks=[]
 def check(name,ok):checks.append(dict(check=name,passed=bool(ok)));assert ok,name
 def near(name,x,y):check(name,math.isclose(float(x),float(y),rel_tol=4e-10,abs_tol=4e-10))
 for name,d in {**p['bindings'],**p['historical_bindings'],**p['external_bindings']}.items():check('hash:'+name,sha(name)==d)
 z=np.load(OUT/'AUX_FROZEN_DRAWS.npz')['z'];check('draw_regeneration',np.array_equal(z,np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_trials'],6))))
 check('4268cells',len(data)==len(grid)==4268);lookup={key(x):x for x in data};check('all_unique',len(lookup)==4268)
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
  if i%440==0:print(f'independent nonideal cells {i}/4268',flush=True)
 for x in data:
  if int(x['side'])==1:
   mirror=(x['anchor'],x['family'],-float(x['b_common_deg']),-float(x['b_diff_deg']),float(x['sigma_pos_m']),float(x['beta_deg']),float(x['range_km']),-1)
   for field in ['relative_range_median','relative_range_p95','relative_range_p99','position_2d_km_p95','failure_rate']:near('mirror:'+x['cell_id']+field,x[field],lookup[mirror][field])
 groups={}
 for x in data:groups.setdefault(tuple(x[k] for k in KEYS),[]).append(x)
 summary=rows(OUT/'AUX_FULL_RANGE_SUMMARY.csv');check('194complete_fullrange_conditions',len(summary)==len(groups)==194)
 for x in summary:
  group=groups[tuple(x[k] for k in KEYS)];worst=max(group,key=lambda y:float(y['relative_range_p95']));v=float(worst['relative_range_p95']);valid=len(group)==22 and {float(y['range_km']) for y in group}==set(p['range_km']) and all(flag(y['execution_valid']) for y in group);finite=all(math.isfinite(float(y['relative_range_p95'])) for y in group)
  near('aggregate:'+str(tuple(x[k] for k in KEYS)),x['worst_p95_relative'],v);near('worst_range:'+str(x),x['worst_range_km'],worst['range_km'])
  for field,val in [('all_execution_valid',valid),('all_p95_finite',finite),('full_range_5pct_pass',valid and finite and v<=.05),('full_range_10pct_pass',valid and finite and v<=.1)]:check(field+str(tuple(x[k] for k in KEYS)),flag(x[field])==val)
  for field,value in [('worst_p99_relative',max(float(y['relative_range_p99']) for y in group)),('worst_p95_position_km',max(float(y['position_2d_km_p95']) for y in group)),('max_failure_rate',max(float(y['failure_rate']) for y in group))]:near(field+str(tuple(x[k] for k in KEYS)),x[field],value)
 def independent_limit(group,axis,tag):
  good=[]
  for cap in sorted({abs(float(x[axis])) for x in group}):
   if all(flag(x['full_range_'+tag+'_pass']) for x in group if abs(float(x[axis]))<=cap):good.append(cap)
  return max(good) if good else 'NOT_REACHED_ON_TESTED_AXIS'
 fronts=[rows(OUT/name) for name in ['AUX_BIAS_REQUIREMENT_FRONTIER.csv','AUX_POSITION_REQUIREMENT_FRONTIER.csv','AUX_DEPLOYMENT_AZIMUTH_REQUIREMENT.csv']]
 for anchor in p['anchors']:
  group=[x for x in summary if x['anchor']==anchor];b,n,g=[next(x for x in f if x['anchor']==anchor) for f in fronts]
  subsets=[([x for x in group if x['family']=='BIAS' and float(x['b_diff_deg'])==0],'b_common_deg',b,'max_tested_abs_common_', '_deg'),([x for x in group if x['family']=='BIAS' and float(x['b_common_deg'])==0],'b_diff_deg',b,'max_tested_abs_diff_','_deg'),([x for x in group if x['family']=='POSITION'],'sigma_pos_m',n,'max_tested_sigma_pos_','_m'),([x for x in group if x['family']=='DEPLOYMENT'],'beta_deg',g,'max_tested_abs_beta_','_deg')]
  for subset,axis,table,prefix,suffix in subsets:
   for tag in ['5pct','10pct']:
    value=independent_limit(subset,axis,tag);field=prefix+tag+suffix
    if isinstance(value,float):near('frontier:'+anchor+field,table[field],value)
    else:check('frontier:'+anchor+field,table[field]==value)
  controls=[next(x for x in group if x['family']==family and all(float(x[k])==0 for k in ['b_common_deg','b_diff_deg','sigma_pos_m','beta_deg'])) for family in ['BIAS','POSITION','DEPLOYMENT']]
  for control in controls:near('zero_controls:'+anchor+control['family'],control['worst_p95_relative'],controls[0]['worst_p95_relative'])
 d=load(OUT/'AUX_GATE2A_DECISION.json');numerical=any(all(float(t[field])>0 for t,field in zip([next(x for x in f if x['anchor']==anchor) for f in [fronts[0],fronts[0],fronts[1],fronts[2]]],['max_tested_abs_common_5pct_deg','max_tested_abs_diff_5pct_deg','max_tested_sigma_pos_5pct_m','max_tested_abs_beta_5pct_deg'])) for anchor in p['anchors'])
 outcome='AUX1_STATIC_MEASUREMENT_CHAIN_ADMITTED' if numerical and p['hardware_all_relevant_requirements_confirmed'] else 'AUX1_NONIDEAL_REQUIREMENTS_ESTABLISHED_CLIENT_CONFIRMATION_REQUIRED' if numerical else 'AUX1_5PCT_ARCHITECTURE_TOO_FRAGILE_UNDER_NONIDEAL_ERRORS'
 check('decision',d['Gate2A_decision']==outcome);check('zero_credit',d['R4_progress_percent']==0 and d['R4_A1_NEW']=='NOT_OPENED' and d['joint_stress_run'] is False)
 start=load(OUT/'AUX_GATE2A_EXECUTION_START.json');check('pushed_before_mc',start['design_sha']==d['design_sha'] and start['remote_verified_before_mc'] and start['policy_sha256']==sha(OUT/'AUX_GATE2A_DESIGN_FREEZE.json'))
 xml=ET.parse(OUT/'AUX_GATE2A_TESTS.xml').getroot();check('tests_pass',len(xml.findall('.//testcase'))==6 and not xml.findall('.//failure') and not xml.findall('.//error'))
 for filename in ['results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv']:
  baseline=subprocess.check_output(['git','show',p['parent_sha']+':'+filename]).decode('utf-8').replace('\r\n','\n');check('append_only:'+filename,Path(filename).read_text(encoding='utf-8').startswith(baseline))
 return checks

def present():
 p=load(OUT/'AUX_GATE2A_DESIGN_FREEZE.json');d=load(OUT/'AUX_GATE2A_DECISION.json');summary=rows(OUT/'AUX_FULL_RANGE_SUMMARY.csv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figdir=OUT/'figures';figdir.mkdir(exist_ok=True);plt.rcParams.update({'svg.hashsalt':'AUX_GATE2A','font.size':10})
 def save(fig,name):
  fig.savefig(figdir/(name+'.png'),dpi=160);fig.savefig(figdir/(name+'.svg'),metadata={'Date':None});plt.close(fig);path=figdir/(name+'.svg');path.write_text('\n'.join(s.rstrip() for s in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 levels=p['bias_levels_deg'];fig,axes=plt.subplots(1,2,figsize=(13,6));fig.subplots_adjust(left=.07,bottom=.19,top=.83,right=.86,wspace=.26)
 for ax,anchor in zip(axes,p['anchors']):
  values=np.array([[100*float(next(x['worst_p95_relative'] for x in summary if x['anchor']==anchor and x['family']=='BIAS' and float(x['b_common_deg'])==c and float(x['b_diff_deg'])==dd)) for dd in levels] for c in levels])
  im=ax.imshow(values,origin='lower',aspect='auto',vmin=0,vmax=10,cmap='viridis')
  for i in range(9):
   for j in range(9):ax.text(j,i,f'{values[i,j]:.1f}',ha='center',va='center',fontsize=8,color='white' if values[i,j]<5 else 'black')
  ax.set_xticks(range(9),levels,rotation=45);ax.set_yticks(range(9),levels);ax.set_xlabel('Differential half-bias (deg)');ax.set_ylabel('Common bias (deg)');ax.set_title('Anchor '+anchor+'; worst11 ranges / mirrors')
 cax=fig.add_axes([.91,.2,.015,.57]);fig.colorbar(im,cax=cax,label='Worst P95 range error (%)');fig.suptitle('Signed bias Cartesian grid; navigation / deployment error zero');save(fig,'FIG1_BIAS_HEATMAP')
 for family,axis,xlabel,name in [('POSITION','sigma_pos_m','Per-axis node-position sigma (m)','FIG2_POSITION_REQUIREMENT'),('DEPLOYMENT','beta_deg','Known actual deployment offset beta (deg)','FIG3_DEPLOYMENT_REQUIREMENT')]:
  fig,ax=plt.subplots(figsize=(8,5));fig.subplots_adjust(left=.12,bottom=.16,top=.86)
  for anchor in p['anchors']:
   group=sorted([x for x in summary if x['anchor']==anchor and x['family']==family],key=lambda x:float(x[axis]));ax.plot([float(x[axis]) for x in group],[100*float(x['worst_p95_relative']) for x in group],'o-',label='Anchor '+anchor)
  ax.axhline(5,color='green',linestyle='--',label='T5');ax.axhline(10,color='black',linestyle=':',label='T10');ax.set_xlabel(xlabel);ax.set_ylabel('Worst tested P95 range error (%)');ax.set_title(f'Single-factor {family.lower()}; signed directions retained');ax.grid(alpha=.25);ax.legend();save(fig,name)
 fig,ax=plt.subplots(figsize=(13,3.8));ax.axis('off');table=[]
 for b,n,g in zip(d['bias_requirements'],d['position_requirements'],d['deployment_requirements']):table.append([b['anchor'],b['baseline_km'],b['sigma_random_deg'],b['max_tested_abs_common_5pct_deg'],b['max_tested_abs_diff_5pct_deg'],n['max_tested_sigma_pos_5pct_m'],g['max_tested_abs_beta_5pct_deg']])
 t=ax.table(cellText=table,colLabels=['Anchor','B km','Random RMS deg','Common bias deg','Diff half-bias deg','Position sigma m','Known beta deg'],loc='center',cellLoc='center');t.auto_set_font_size(False);t.set_fontsize(10);t.scale(1,2);ax.set_title('T5 maximum tested single-factor tolerances; NOT a jointly validated engineering point',pad=25);save(fig,'FIG4_ARCHITECTURE_REQUIREMENT_SUMMARY')
 lines=['# AUX Gate2A nonideal measurement requirements','',d['Gate2A_decision'],'','All reported budgets are SINGLE FACTOR, not a simultaneous deployable error budget. No joint bias/navigation/deployment stress was run.','', '```json',json.dumps({k:d[k] for k in ['bias_requirements','position_requirements','deployment_requirements']},indent=2),'```','',
  'The two anchors are frozen by the research lead: A B5km/random0.05deg, B B7km/random0.075deg. B6km/random0.10deg is preserved only as a historical threshold reference; no new nonideal sweep for it. New zero-control quantiles use a new frozen30k sample table and may differ from Gate1 finite-MC results; they do not replace accepted Gate1 evidence.',
  '', '## Model and interpretation','',
  'Random bearing terms independent zero-mean Gaussian at each anchor sigma. Fixed b1=b_common-b_diff, b2=b_common+b_diff; b_diff is half the difference, so a0.05deg bound corresponds to0.10deg b2-b1. Every signed81 combination kept. Common bias is a shared reference-like parametrization; differential bias is relative mismatch-like, not a verified physical decomposition.',
  'Position error per-axis1-sigma in metres: both estimated node positions receive independent2D Gaussian offsets; observations use true geometry. Range output is norm(p_hat-MAIN_est), compared with true MAIN-to-target R; position error is norm(p_hat-target_true). Navigation-error common translations can cancel in relative range while remaining in absolute2D position; both metrics reported.',
  'Deployment beta modifies true AUX=B[cos(90+beta),sin(90+beta)] and reflected placement. Estimator knows actual position exactly, so this is geometry tolerance, not navigation error. Beta is not converted into a bearing bias. Other nonideal factors zero in each single-factor family.',
  'Exactly30000 trials per registered cell; one frozen30000x6 Gaussian table shared across cells. First2 coordinates bearing, next2 MAIN position, last2 AUX position; position arrays only used in POSITION family. Mirrors reflect noise and navigation y components. Bias mirror compares(c,d,side+) with(-c,-d,side-). Mirrors are controls, not independent trials.',
  'All11 discrete50:1:60km ranges and mirrors included. Full-range P95=max over22 cells, never pooled or averaged across signs. Nearest-rank unconditional quantiles; parallel/nonfinite/behind errors assigned infinity, no outlier removal. Approximate P95 rank/Wilson intervals diagnose finite-MC uncertainty, not simultaneous assurance.',
  'Axis frontiers require every registered inner signed node through the reported absolute cap to pass; common axis diff=0, differential axis common=0. Full bias map must be consulted for combinations: separate maximum axis limits cannot be applied simultaneously. Position and beta caps likewise preserve all lower tested magnitudes/signs. No interpolation or continuous tolerance certificate. A cap at the largest tested magnitude is right-censored; larger values are untested.',
  '', '## Hardware and next decision','',
  'Actual auxiliary bearing array/capability, navigation/attitude, clock and calibration specifications remain UNKNOWN / CLIENT_CONFIRMATION_REQUIRED. Simulated arrays, synthetic0.1deg noise and ideal navigation assumptions are not hardware evidence. See source-bound AUX_HARDWARE_EVIDENCE_AUDIT.md and AUX_CLIENT_REQUIREMENTS.md.',
  'Independent zero-mean random components may enter variance addition only if their independence is substantiated. Systematic biases stay separate; correlated random terms require covariance. No bias is silently folded into Gaussian RMS.',
  'Numerical nonzero single-factor budgets do not admit an actual measurement chain. No Gate2B joint budget, time-skew/association MC, R4-A1-NEW, acoustic/depth/SSP/TDOA/tracker/5D. Single-array depth stays closed; R4=0%. Commit/push/verify then stop for research-lead audit and client confirmation.']
 report='\n'.join(lines)+'\n'
 for name in ['AUX_GATE2A_REPORT.md','GPT_SYNC.md']:(OUT/name).write_text(report,encoding='utf-8')
 client=(OUT/'AUX_CLIENT_REQUIREMENTS_TEMPLATE.md').read_text(encoding='utf-8-sig')
 for b,n,g in zip(d['bias_requirements'],d['position_requirements'],d['deployment_requirements']):
  anchor=b['anchor']
  for tag in ['5pct','10pct']:
   values={'COMMON':b['max_tested_abs_common_'+tag+'_deg'],'DIFF':b['max_tested_abs_diff_'+tag+'_deg'],'POS':n['max_tested_sigma_pos_'+tag+'_m'],'BETA':g['max_tested_abs_beta_'+tag+'_deg']}
   for key,value in values.items():client=client.replace('{{'+anchor+'_'+tag+'_'+key+'}}',str(value))
 assert '{{' not in client
 (OUT/'AUX_CLIENT_REQUIREMENTS.md').write_text(client,encoding='utf-8')
 with Path('results/R4_MASTER/R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## AUX Gate2A execution\n\n'+d['Gate2A_decision']+'; PENDING_RESEARCH_LEAD_AUDIT. Single-factor signed-bias/navigation/deployment budgets only, not a jointly validated point. [Report](../R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET/AUX_GATE2A_REPORT.md). Actual hardware unknown; R4=0%; R4-A1-NEW/Gate2B not opened, depth closed. Commit/push/verify then stop.\n')
 with Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow(['AUX-GATE2A','../R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET/AUX_GATE2A_DECISION.json',d['Gate2A_decision'],'Single-factor signed error budgets; two anchors;11 ranges; no joint stress','PENDING_RESEARCH_LEAD_AUDIT; HARDWARE_UNKNOWN; NO_PROGRESS_CREDIT',0])

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--finalize',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 if args.finalize:
  assert not (OUT/'VALIDATION.json').exists(),'No finalize rerun'
  checks=audit();present()
  with (OUT/'AUX_INDEPENDENT_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['check','passed']);w.writeheader();w.writerows(checks)
  dump(OUT/'VALIDATION.json',dict(status='PASS',independent_checks=len(checks),failed=0,tests_passed=6,primary_runs=1,registered_cells=4268,trials_per_cell=30000,joint_stress_run=False,new_acoustic_propagation=0,new_depth_scores=0))
  paths=[x for x in OUT.rglob('*') if x.is_file()]+[Path('results/R4_MASTER/R4_PLAN.md'),Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv')];dump(OUT/'AUX_RESULTS_FREEZE.json',dict(bindings={x.as_posix():sha(x) for x in sorted(paths)}));print(json.dumps(dict(finalized='PASS',checks=len(checks))),flush=True)
 elif args.verify:
  seal=load(OUT/'AUX_RESULTS_FREEZE.json')
  for name,d in seal['bindings'].items():assert sha(name)==d,name
  checks=audit();print(json.dumps(dict(sealed='PASS',checks=len(checks),failed=0)),flush=True)
 else:ap.error('--finalize or --verify required')
