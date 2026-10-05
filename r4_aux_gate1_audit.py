"""Independent Gate1 slope-intersection audit and saved-data requirement plots."""
from pathlib import Path
import argparse,csv,hashlib,json,math,subprocess,xml.etree.ElementTree as ET
import numpy as np
OUT=Path('results/R4_AUX_GATE1_CROSSTRACK_REQUIREMENT')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def flag(x):return x is True or x=='True'

def audit():
 p=load(OUT/'AUX_GATE1_DESIGN_FREEZE.json');data=rows(OUT/'AUX_CROSSTRACK_RESULTS.csv');grid=rows(OUT/'AUX_CROSSTRACK_GRID.csv');checks=[]
 def check(name,ok):
  checks.append(dict(check=name,passed=bool(ok)));assert ok,name
 def near(name,a,b,atol=3e-10):check(name,math.isclose(float(a),float(b),rel_tol=3e-10,abs_tol=atol))
 for name,d in {**p['bindings'],**p['historical_bindings'],**p['external_bindings']}.items():check('hash:'+name,sha(name)==d)
 z=np.load(OUT/'AUX_FROZEN_STANDARD_NORMAL_DRAWS.npz')['z'];expected=np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_samples'],2));check('frozen_draws',np.array_equal(z,expected))
 check('complete1870cells',len(data)==len(grid)==1870)
 lookup={};zz=1.959963984540054;n=len(z)
 for i,(x,g) in enumerate(zip(data,grid),1):
  name=x['cell_id']
  for key,val in g.items():check(name+':grid:'+key,x[key]==val)
  r,b,s=[float(x[k]) for k in ['range_km','baseline_km','sigma_deg']];side=int(x['side']);alpha=math.atan(b/r);sig=math.radians(s)
  near(name+':crossing_angle',x['alpha_deg'],math.degrees(alpha));near(name+':aux_bearing',x['aux_true_bearing_deg'],-side*math.degrees(alpha))
  n1=side*z[:,0]*sig;n2=-side*alpha+side*z[:,1]*sig
  # Independent straight-line slopes: y=m1*x and y-side*B=m2*x.
  m1=np.tan(n1);m2=np.tan(n2);xx=side*b/(m1-m2);yy=m1*xx;det=np.sin(n2-n1)
  parallel=np.abs(det)<=p['parallel_det_tolerance'];behind=(xx*np.cos(n1)+yy*np.sin(n1)<=0)|(xx*np.cos(n2)+(yy-side*b)*np.sin(n2)<=0)
  nonfinite=~np.isfinite(xx)|~np.isfinite(yy);failed=parallel|behind|nonfinite
  er=np.abs(np.hypot(xx,yy)-r)/r;pos=np.hypot(xx-r,yy);cross=np.abs(yy)
  for v in [er,pos,cross]:v[failed]=np.inf
  sorted_er=np.sort(er)
  for metric,v in [('relative_range',er),('range_abs_km',er*r),('position_2d_km',pos),('cross_range_abs_km',cross)]:
   ordered=np.sort(v)
   for tag,q in [('median',.5),('p90',.9),('p95',.95),('p99',.99)]:
    expected_v=ordered[math.ceil(n*q)-1]
    if math.isfinite(expected_v):near(name+':'+metric+tag,x[metric+'_'+tag],expected_v)
    else:check(name+':'+metric+tag,x[metric+'_'+tag]=='INF')
  for key,value in [('failure_rate',failed.mean()),('behind_sensor_rate',behind.mean()),('parallel_rate',parallel.mean()),('nonfinite_rate',nonfinite.mean()),('ill_conditioned_rate',(1/np.maximum(np.abs(det),np.finfo(float).tiny)>p['ill_condition_proxy_limit']).mean()),('fraction_le5pct',(er<=.05).mean()),('fraction_le10pct',(er<=.1).mean())]:near(name+':'+key,x[key],value)
  # Closed form independently derived from linear bearing sensitivities.
  sx=sig*math.sqrt(r**4+(r*r+b*b)**2)/b;sy=sig*r;lin=zz*sx/r
  near(name+':linear_range_sigma',x['sigma_r_linear_km'],sx);near(name+':linear_cross_sigma',x['sigma_cross_linear_km'],sy);near(name+':linear_p95',x['p95_linear_relative'],lin)
  H=np.array([[0,1/r],[side*b/(r*r+b*b),r/(r*r+b*b)]]);singular=np.linalg.svd(H,compute_uv=False)
  near(name+':H_condition',x['jacobian_condition'],singular[0]/singular[1]);near(name+':F_condition',x['fim_condition'],(singular[0]/singular[1])**2,1e-7)
  nonlinear=float(x['relative_range_p95']);near(name+':nonlinear_ratio',x['nonlinear_to_linear_p95_ratio'],nonlinear/lin);near(name+':nonlinear_deviation',x['nonlinear_deviation_relative'],nonlinear-lin)
  scaling=sig/math.sin(alpha);near(name+':scaling_parameter',x['sigma_rad_over_sin_alpha'],scaling);near(name+':scaling_ratio',x['p95_over_sigma_rad_over_sin_alpha'],nonlinear/scaling)
  delta=zz*math.sqrt(n*.95*.05)
  for key,ix in [('p95_mc_rank_ci_low',max(0,math.floor(n*.95-delta)-1)),('p95_mc_rank_ci_high',min(n-1,math.ceil(n*.95+delta)-1))]:near(name+':'+key,x[key],sorted_er[ix])
  for tag,threshold in [('5pct',.05),('10pct',.1)]:
   prop=float((er<=threshold).mean());den=1+zz*zz/n;c=(prop+zz*zz/(2*n))/den;h=zz*math.sqrt(prop*(1-prop)/n+zz*zz/(4*n*n))/den
   near(name+':wilson_low'+tag,x['fraction_'+tag+'_wilson_low'],c-h);near(name+':wilson_high'+tag,x['fraction_'+tag+'_wilson_high'],c+h)
  check(name+':valid_and_count',flag(x['execution_valid']) and int(x['n_samples'])==n);lookup[(r,b,s,side)]=x
  if i%220==0:print(f'independent geometry checks {i}/1870',flush=True)
 for (r,b,s,side),x in lookup.items():
  if side==-1:
   y=lookup[(r,b,s,1)]
   for key in ['relative_range_median','relative_range_p95','position_2d_km_p95','failure_rate','sigma_r_linear_km']:near('mirror:'+str((r,b,s))+key,x[key],y[key])
 summary=rows(OUT/'AUX_FULL_RANGE_SUMMARY.csv');check('85_full_range_rows',len(summary)==85)
 for x in summary:
  s,b=float(x['sigma_deg']),float(x['baseline_km']);group=[lookup[(r,b,s,side)] for r in p['range_km'] for side in p['sides']];worst=max(group,key=lambda y:float(y['relative_range_p95']));value=float(worst['relative_range_p95']);valid=all(flag(y['execution_valid']) for y in group);finite=all(math.isfinite(float(y['relative_range_p95'])) for y in group)
  near('summary_value'+str((s,b)),x['worst_p95_relative'],value);check('summary_count'+str((s,b)),int(x['n_ranges'])==11 and int(x['n_mirror_cells'])==22)
  check('summary_worst_range'+str((s,b)),float(x['worst_range_km'])==float(worst['range_km']))
  for key,v in [('all_ranges_execution_valid',valid),('all_p95_finite',finite),('full_range_5pct_pass',valid and finite and value<=.05),('full_range_10pct_pass',valid and finite and value<=.1)]:check(key+str((s,b)),flag(x[key])==v)
  for key,value in [('alpha_min_deg',min(float(y['alpha_deg']) for y in group)),('alpha_max_deg',max(float(y['alpha_deg']) for y in group)),('max_failure_rate',max(float(y['failure_rate']) for y in group)),('max_nonlinear_to_linear_ratio',max(float(y['nonlinear_to_linear_p95_ratio']) for y in group)),('min_nonlinear_to_linear_ratio',min(float(y['nonlinear_to_linear_p95_ratio']) for y in group))]:near(key+str((s,b)),x[key],value)
 frontier=rows(OUT/'AUX_BASELINE_REQUIREMENT_FRONTIER.csv');check('all_five_sigma_frontiers',len(frontier)==5)
 for x in frontier:
  s=float(x['sigma_deg'])
  for tag in ['5pct','10pct']:
   good=[y for y in summary if float(y['sigma_deg'])==s and flag(y['full_range_'+tag+'_pass'])];key='minimum_tested_B_for_'+tag+'_km';best=min(good,key=lambda y:float(y['baseline_km'])) if good else None
   if best:
    near(key+str(s),x[key],best['baseline_km']);near('frontier_worst_range'+tag+str(s),x['worst_range_at_'+tag+'_boundary'],best['worst_range_km']);near('frontier_worst_value'+tag+str(s),x['worst_p95_at_'+tag+'_boundary'],best['worst_p95_relative'])
   else:check(key+str(s),x[key]=='NOT_REACHED_WITHIN_B<=10KM')
 fixed=rows(OUT/'AUX_BEARING_REQUIREMENT_AT_FIXED_BASELINE.csv');check('fixed_baselines_complete',[float(x['baseline_km']) for x in fixed]==p['fixed_baseline_km'])
 for x in fixed:
  b=float(x['baseline_km'])
  for tag in ['5pct','10pct']:
   group=sorted([y for y in summary if float(y['baseline_km'])==b],key=lambda y:float(y['sigma_deg']));good=[y for y in group if flag(y['full_range_'+tag+'_pass'])];best=max(good,key=lambda y:float(y['sigma_deg'])) if good else None
   key='largest_tested_sigma_for_'+tag+'_deg'
   if best:
    near(key+str(b),x[key],best['sigma_deg']);above=[y for y in group if float(y['sigma_deg'])>float(best['sigma_deg']) and not flag(y['full_range_'+tag+'_pass'])]
    if above:near('next_failure'+tag+str(b),x['next_tested_failing_sigma_for_'+tag+'_deg'],above[0]['sigma_deg'])
    else:check('next_unbracketed'+tag+str(b),x['next_tested_failing_sigma_for_'+tag+'_deg']=='NOT_BRACKETED_ABOVE_TESTED_MAX')
   else:check(key+str(b),x[key]=='NOT_REACHED_IN_TESTED_SIGMA_GRID')
 d=load(OUT/'AUX_GATE1_DECISION.json');a=any(float(x['sigma_deg'])==.1 and float(x['baseline_km'])<=5 and flag(x['full_range_5pct_pass']) for x in summary);b=any(float(x['sigma_deg'])>=.05 and float(x['baseline_km'])<=10 and flag(x['full_range_5pct_pass']) for x in summary)
 decision='AUX1_CROSSTRACK_GEOMETRY_PRACTICAL_AT_ORIGINAL_SCALE' if a else 'AUX1_BASELINE_BEARING_TRADEOFF_ESTABLISHED' if b else 'AUX1_FULL_RANGE_5PCT_NOT_ESTABLISHED'
 check('frozen_outcome',d['primary_architecture_decision']==decision);check('zero_progress',d['R4_progress_percent']==0);check('A1_not_open',d['R4_A1_NEW']=='NOT_OPENED');check('hardware_unknown',d['actual_bearing_hardware_capability']=='UNKNOWN')
 start=load(OUT/'AUX_GATE1_EXECUTION_START.json');check('policy_push_before_mc',start['design_sha']==d['design_sha'] and start['remote_verified_before_mc'] and start['policy_sha256']==sha(OUT/'AUX_GATE1_DESIGN_FREEZE.json'))
 xml=ET.parse(OUT/'AUX_GATE1_TESTS.xml').getroot();check('unit_tests_pass',len(xml.findall('.//testcase'))==8 and not xml.findall('.//failure') and not xml.findall('.//error'))
 for filename in ['results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv']:
  parent=subprocess.check_output(['git','show',p['parent_sha']+':'+filename]).decode('utf-8').replace('\r\n','\n');check('append_only:'+filename,Path(filename).read_text(encoding='utf-8').startswith(parent))
 return checks

def present():
 p=load(OUT/'AUX_GATE1_DESIGN_FREEZE.json');d=load(OUT/'AUX_GATE1_DECISION.json');data=rows(OUT/'AUX_CROSSTRACK_RESULTS.csv');summary=rows(OUT/'AUX_FULL_RANGE_SUMMARY.csv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figdir=OUT/'figures';figdir.mkdir(exist_ok=True);plt.rcParams.update({'svg.hashsalt':'AUX_GATE1_CROSSTRACK','font.size':10})
 def save(fig,name):
  fig.tight_layout();fig.savefig(figdir/(name+'.png'),dpi=170);fig.savefig(figdir/(name+'.svg'),metadata={'Date':None});plt.close(fig);path=figdir/(name+'.svg');path.write_text('\n'.join(x.rstrip() for x in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 fig,ax=plt.subplots(figsize=(9,5.5))
 for s in p['sigma_deg']:
  group=[x for x in summary if float(x['sigma_deg'])==s];ax.plot([float(x['baseline_km']) for x in group],[100*float(x['worst_p95_relative']) for x in group],'o-',label=f'{s:g} deg')
 ax.axhline(5,color='green',linestyle='--',label='T5');ax.axhline(10,color='black',linestyle=':',label='T10');ax.set_xlabel('Tested cross-track baseline (km)');ax.set_ylabel('Worst tested P95 relative range error (%)');ax.set_title('11 tested ranges (50:1:60 km), mirrors; effective bearing RMS');ax.grid(alpha=.25);ax.legend(title='RMS per node');save(fig,'FIG1_BASELINE_VS_WORST_P95')
 fig,ax=plt.subplots(figsize=(8,5))
 for tag,marker in [('5pct','o'),('10pct','s')]:
  vals=[x for x in d['requirement_frontier'] if isinstance(x['minimum_tested_B_for_'+tag+'_km'],(int,float))];ax.plot([x['sigma_deg'] for x in vals],[x['minimum_tested_B_for_'+tag+'_km'] for x in vals],marker+'--',label=tag+' minimum tested B')
 for x in d['requirement_frontier']:
  if isinstance(x['minimum_tested_B_for_5pct_km'],str):ax.annotate('T5 not reached at B<=10 km',(x['sigma_deg'],10),xytext=(-165,-30),textcoords='offset points',arrowprops={'arrowstyle':'->'})
 ax.set_xlabel('Effective directed-bearing RMS (deg)');ax.set_ylabel('Minimum tested cross-track baseline (km)');ax.set_title('Discrete requirement frontier; connectors do not interpolate thresholds');ax.set_ylim(1,11);ax.grid(alpha=.25);ax.legend();save(fig,'FIG2_BASELINE_REQUIREMENT_FRONTIER')
 selected=sorted({2.,5.,10.}|{float(x['minimum_tested_B_for_5pct_km']) for x in d['requirement_frontier'] if isinstance(x['minimum_tested_B_for_5pct_km'],(int,float))})
 fig,ax=plt.subplots(figsize=(8,5))
 for b in selected:
  group=[x for x in data if float(x['baseline_km'])==b and float(x['sigma_deg'])==.1 and int(x['side'])==1];ax.plot([float(x['range_km']) for x in group],[float(x['alpha_deg']) for x in group],'o-',label=f'B={b:g} km')
 ax.set_xlabel('Tested target range (km)');ax.set_ylabel('True crossing angle (deg)');ax.set_title('Fixed cross-track deployment: alpha=atan(B/R)');ax.grid(alpha=.25);ax.legend();save(fig,'FIG3_CROSSING_ANGLE_VS_RANGE')
 fig,ax=plt.subplots(figsize=(7,5))
 unique=[x for x in data if int(x['side'])==1];xx=[100*float(x['p95_linear_relative']) for x in unique];yy=[100*float(x['relative_range_p95']) for x in unique];colors=[float(x['sigma_deg']) for x in unique];scatter=ax.scatter(xx,yy,c=colors,s=14,cmap='viridis',alpha=.6);lim=max(max(xx),max(yy));ax.plot([0,lim],[0,lim],'k--');fig.colorbar(scatter,ax=ax,label='Effective bearing RMS (deg)');ax.set_xlabel('Local linearized P95 relative range (%)');ax.set_ylabel('Nonlinear MC P95 relative range (%)');ax.set_title('Local approximation is diagnostic; MC governs Gate');ax.grid(alpha=.25);save(fig,'FIG4_LINEAR_VS_NONLINEAR')
 ratios=[float(x['nonlinear_to_linear_p95_ratio']) for x in unique];scale=[float(x['p95_over_sigma_rad_over_sin_alpha']) for x in unique]
 violations=[]
 for s in p['sigma_deg']:
  group=sorted([x for x in summary if float(x['sigma_deg'])==s],key=lambda x:float(x['baseline_km']))
  for a,b in zip(group,group[1:]):
   if float(b['worst_p95_relative'])>float(a['worst_p95_relative']):violations.append([s,a['baseline_km'],b['baseline_km']])
 dump(OUT/'AUX_SCALING_DIAGNOSTIC.json',dict(n_distinct_positive_side_cells=len(unique),mirrors_not_independent=True,nonlinear_over_linear_min=min(ratios),nonlinear_over_linear_max=max(ratios),p95_over_sigma_rad_over_sin_alpha_min=min(scale),p95_over_sigma_rad_over_sin_alpha_max=max(scale),worst_p95_baseline_monotonicity_violations=violations,diagnostic_only=True,formula='sigma_r= sigma_rad*sqrt(R^4+(R^2+B^2)^2)/B; far-field P95 approximately1.96*sqrt(2)*sigma_rad/sin(alpha)'))
 lines=['# AUX Gate1 cross-track baseline and bearing requirement','',d['primary_architecture_decision'],'',
  'Parent79bc6b5c04e08d6074719e3c045947c1325253c8 accepted conditional Gate0. Gate1 uses a new frozen deployment, not a retuned Gate0 triangle. MAIN=(0,0), target=(R,0), AUX=(0,+/-B); alpha=atan(B/R). AUX angle is never adjusted per target range.','',
  '## Requirement frontier','', '| Effective bearing RMS deg | Minimum tested B for5% km | Minimum tested B for10% km |','|---|---|---|']
 for row in d['requirement_frontier']:lines.append(f'| {row["sigma_deg"]:g} | {row["minimum_tested_B_for_5pct_km"]} | {row["minimum_tested_B_for_10pct_km"]} |')
 lines+=['','All11 tested ranges50:1:60km and both mirrors must be valid/finite; worst nonlinear MC P95 governs T5/T10. Finite sample requirement on a discrete grid, no continuous range certificate, no interpolated B/sigma threshold.','',
  'Fixed-baseline bearing requirements:','```json',json.dumps(d['bearing_requirement_fixed_baseline'],indent=2),'```','', 'Primary sigma0.10-deg selection:','```json',json.dumps(d['primary_selected_crossing_angle_range'],indent=2),'```','',
  '## Sampling and propagation of error','',
  'Each of1870 registered mirror cells evaluates50000 trials. The shared single50000x2 PCG64 standard-normal table compares cells; cell-trials are not independent trials across geometry. Mirrors reverse both Gaussian perturbations and are symmetry controls. No repeat seed, new scenario or post-result grid addition. Independent nodes have zero-mean Gaussian effective directed-bearing error of equal RMS; truth enters scene/evaluation only, never the intersection estimator.',
  'Parallel/nonfinite/behind-sensor intersections receive infinite errors and remain in unconditional median/P90/P95/P99. Finite poorly conditioned intersections remain. Failure-rate metrics may overlap. Empirical nearest-rank quantiles, no outlier trimming. P95 rank intervals are approximate binomial-normal diagnostics; Wilson intervals describe proportions<=5/10%. No simultaneous selected-grid or deployment guarantee.',
  'The local linearized covariance check uses an independently derived sensitivity. sigma_r= sigma_rad*sqrt(R^4+(R^2+B^2)^2)/B. At B/R small this gives P95 approximately1.96*sqrt(2)*sigma_rad/sin(alpha). Ratio/scaling diagnostics use actual MC, not surrogate Gate results. Local inverse FIM is a linearized lower-bound/approximation, not estimator guarantee.',
  f'Nonlinear/linear P95 ratio range: {min(ratios):.6f}--{max(ratios):.6f}. MC P95 divided by sigma_rad/sin(alpha): {min(scale):.6f}--{max(scale):.6f}. Baseline monotonicity violations: {len(violations)}. See AUX_SCALING_DIAGNOSTIC.json.',
  '', '## Ideal geometry and measurement requirement','',
  'MAIN/AUX positions exact, timestamps simultaneous, association correct, directed bearing resolved. Perfect cross-track alignment to the reference main LoB is a Gate1 ideal deployment assumption; uncertainty of initial bearing-guided placement/relative alignment is not simulated here. Node positioning, attitude, time offset, systematic bearing bias/correlation are Gate2 topics and NOT_OPENED.',
  'Effective directed-bearing RMS is a total input requirement, not beamformer-only CRLB or known receiver capability. Future measurement budget must include DOA random error, calibration residuals, heading/attitude and relative alignment. Gaussian zero-mean independence excludes systematic/common biases in this Gate. No numerical allocation to components is inferred, no SNR/aperture/snapshot guess.',
  'Existing project materials leave actual element count/aperture/navigation/heading/bearing accuracy UNKNOWN / CLIENT_CONFIRMATION_REQUIRED. Historical model arrays and0.1-degree synthetic noise are not device acceptance. See AUX_BEARING_MEASUREMENT_REQUIREMENT.md and source-bound original excerpts.',
  '', '## Architecture boundary','',
  'Single-array cold start NOT ESTABLISHED; single-array depth engineering CLOSED; exact-horizontal mechanism oracle evidence retained. AUX1 remains CONDITIONAL_ARCHITECTURE_CANDIDATE, even if the requirement tradeoff is established geometrically. No actual auxiliary bearing chain, node pose/time/association/deployment evidence supplied. No velocity/heading/depth precision claim.',
  'R4-A1-NEW NOT_OPENED; depth CLOSED_PENDING_INDEPENDENT_HORIZONTAL_ACQUISITION; AUX Gate2 NOT_OPENED. R4=0%, original A/B sequence paused; A2/A3/A4/B2/B3/B4/C paused/not opened; P5 not opened. No acoustic TL/CZ/depth/SSP/TDOA/5D estimator/tracker. Commit/push/remote verify and STOP.']
 report='\n'.join(lines)+'\n'
 for name in ['AUX_GATE1_REPORT.md','GPT_SYNC.md']:(OUT/name).write_text(report,encoding='utf-8')
 with Path('results/R4_MASTER/R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## AUX Gate1 requirement execution\n\n'+d['primary_architecture_decision']+'; PENDING_RESEARCH_LEAD_AUDIT. Fixed cross-track deployment, discrete11-range requirement frontier; no hardware or continuous-range guarantee. [Report](../R4_AUX_GATE1_CROSSTRACK_REQUIREMENT/AUX_GATE1_REPORT.md). Single-array cold start not established, depth closed, AUX1 conditional, R4=0%; no R4-A1-NEW/Gate2 restart. Commit/push/verify then stop.\n')
 with Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow(['AUX-GATE1','../R4_AUX_GATE1_CROSSTRACK_REQUIREMENT/AUX_GATE1_DECISION.json',d['primary_architecture_decision'],'Cross-track baseline vs effective bearing RMS requirement; discrete11 ranges; ideal geometry','PENDING_RESEARCH_LEAD_AUDIT; NO_HARDWARE_OR_SYSTEM_CREDIT',0])

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--finalize',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 if args.finalize:
  assert not (OUT/'VALIDATION.json').exists(),'No finalization rerun'
  checks=audit();present()
  with (OUT/'AUX_INDEPENDENT_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['check','passed']);w.writeheader();w.writerows(checks)
  dump(OUT/'VALIDATION.json',dict(status='PASS',independent_checks=len(checks),failed=0,tests_passed=8,tests_failed=0,primary_mc_runs=1,registered_cells=1870,trials_per_cell=50000,shared_draws_not_independent_repetitions=True,new_acoustic_propagation=0,new_depth_scores=0,independent_method='No primary implementation import; slope intersection, closed-form linear sensitivity, raw frozen draws'))
  paths=[x for x in OUT.rglob('*') if x.is_file()]+[Path('results/R4_MASTER/R4_PLAN.md'),Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv')]
  dump(OUT/'AUX_RESULTS_FREEZE.json',dict(bindings={x.as_posix():sha(x) for x in sorted(paths)}));print(json.dumps(dict(finalized='PASS',checks=len(checks))),flush=True)
 elif args.verify:
  seal=load(OUT/'AUX_RESULTS_FREEZE.json')
  for name,digest in seal['bindings'].items():assert sha(name)==digest,name
  checks=audit();print(json.dumps(dict(sealed='PASS',checks=len(checks),failed=0)),flush=True)
 else:ap.error('--finalize or --verify required')
