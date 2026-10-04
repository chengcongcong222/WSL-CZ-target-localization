"""Independent cold geometry audit and saved-result presentation; no acoustic imports."""
from pathlib import Path
import argparse,csv,hashlib,json,math,subprocess
import numpy as np
OUT=Path('results/R4_ROUTE_RESET_AUXILIARY_NODE_GATE0')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def flag(v):return v is True or v=='True'

def audit():
 p=load(OUT/'AUX_GATE0_DESIGN_FREEZE.json');data=rows(OUT/'AUX_GEOMETRY_RESULTS.csv');grid=rows(OUT/'AUX_GEOMETRY_GRID.csv');checks=[]
 def check(name,passed):
  checks.append(dict(check=name,passed=bool(passed)))
  assert passed,name
 def near(name,x,y,atol=2e-10):check(name,math.isclose(float(x),float(y),rel_tol=2e-10,abs_tol=atol))
 for path,d in {**p['bindings'],**p['historical_bindings']}.items():check('hash:'+path,sha(path)==d)
 z=np.load(OUT/'AUX_FROZEN_STANDARD_NORMAL_DRAWS.npz')['z'];expected=np.random.Generator(np.random.PCG64(p['seed'])).standard_normal((p['n_samples'],2));check('frozen_gaussian_draws',np.array_equal(z,expected))
 check('complete_grid',len(data)==len(grid)==1920);check('ordered_grid_ids',[x['cell_id'] for x in data]==[x['cell_id'] for x in grid])
 feasible=0
 for row,design in zip(data,grid):
  ident=row['cell_id']
  for key,val in design.items():check(ident+':design:'+key,row[key]==val)
  r,b,a,s=[float(row[k]) for k in ['range_km','baseline_km','alpha_deg','sigma_deg']];side=int(row['side']);ar=math.radians(a)
  delta=b*b-r*r*math.sin(ar)**2;valid=delta>=-1e-12
  check(ident+':feasible',valid==(row['geometry_status']=='FEASIBLE'))
  if not valid:
   check(ident+':no_fictional_metrics',row['relative_range_p95']=='' and row['n_samples']=='');continue
  feasible+=1
  d=r*math.cos(ar)+(1 if row['branch']=='FAR' else -1)*math.sqrt(max(0,delta));p2=np.array([r-d*math.cos(ar),side*d*math.sin(ar)])
  near(ident+':baseline',np.linalg.norm(p2),b);near(ident+':aux_distance',np.linalg.norm(np.array([r,0])-p2),d)
  for k,v in [('aux_x_km',p2[0]),('aux_y_km',p2[1]),('aux_range_km',d)]:near(ident+':'+k,row[k],v)
  n1=z[:,0]*math.radians(s)*side;n2=-side*ar+z[:,1]*math.radians(s)*side
  # Independent slope intersection, not the runner's cross-product formula.
  m1=np.tan(n1);m2=np.tan(n2);x=(p2[1]-m2*p2[0])/(m1-m2);y=m1*x
  det=np.sin(n2-n1);parallel=np.abs(det)<=p['parallel_det_tolerance']
  behind=(x*np.cos(n1)+y*np.sin(n1)<=0)|((x-p2[0])*np.cos(n2)+(y-p2[1])*np.sin(n2)<=0)
  numerical=~np.isfinite(x)|~np.isfinite(y);failed=parallel|behind|numerical
  er=np.abs(np.hypot(x,y)-r)/r;cross=np.abs(y);pos=np.hypot(x-r,y)
  for v in [er,cross,pos]:v[failed]=np.inf
  for metric,v in [('range_abs_km',er*r),('relative_range',er),('cross_range_abs_km',cross),('position_2d_km',pos)]:
   ordered=np.sort(v)
   for tag,q in [('median',.5),('p90',.9),('p95',.95),('p99',.99)]:
    expected_v=float(ordered[math.ceil(len(v)*q)-1]);saved=row[metric+'_'+tag]
    if math.isfinite(expected_v):near(ident+':'+metric+':'+tag,saved,expected_v)
    else:check(ident+':infinite_quantile:'+metric+tag,saved=='INF')
  for k,v in [('failure_rate',failed.mean()),('behind_sensor_rate',behind.mean()),('parallel_rate',parallel.mean()),('numerical_failure_rate',numerical.mean()),('ill_conditioned_rate',(1/np.maximum(np.abs(det),np.finfo(float).tiny)>p['ill_condition_proxy_limit']).mean()),('fraction_range_le1pct',(er<=.01).mean()),('fraction_range_le5pct',(er<=.05).mean()),('fraction_range_le10pct',(er<=.1).mean())]:near(ident+':'+k,row[k],v)
  # Closed-form sensitivity: dy=R*dtheta1; dx=(R*cos(alpha)*dtheta1-d*dtheta2)/(-side*sin(alpha)).
  sig=math.radians(s);sx=sig*math.sqrt((r*math.cos(ar))**2+d*d)/math.sin(ar);sy=sig*r
  near(ident+':linear_range',row['linear_range_sigma_km'],sx);near(ident+':linear_cross',row['linear_cross_sigma_km'],sy);near(ident+':linear_p95',row['linear_range_p95_relative'],1.959963984540054*sx/r)
  H=np.array([[0,1/r],[-side*math.sin(ar)/d,math.cos(ar)/d]]);sv=np.linalg.svd(H,compute_uv=False);near(ident+':jacobian_condition',row['jacobian_condition'],sv[0]/sv[1]);near(ident+':fim_condition',row['fim_condition'],(sv[0]/sv[1])**2,1e-7)
  n=len(z);q=float((er<=.05).mean());zz=1.959963984540054;den=1+zz*zz/n;c=(q+zz*zz/(2*n))/den;h=zz*math.sqrt(q*(1-q)/n+zz*zz/(4*n*n))/den
  near(ident+':wilson_low',row['fraction_le5pct_wilson_low'],c-h);near(ident+':wilson_high',row['fraction_le5pct_wilson_high'],c+h)
  v=np.sort(er);dd=zz*math.sqrt(n*.95*.05)
  for key,ix in [('relative_range_p95_mc_ci_low',max(0,math.floor(n*.95-dd)-1)),('relative_range_p95_mc_ci_high',min(n-1,math.ceil(n*.95+dd)-1))]:near(ident+':'+key,row[key],v[ix])
 check('feasible_count',feasible==240)
 for s in rows(OUT/'AUX_GEOMETRY_SUMMARY.csv'):
  good=[x for x in data if x['geometry_status']=='FEASIBLE' and all(float(x[k])==float(s[k]) for k in ['baseline_km','alpha_deg','sigma_deg'])]
  check('summary_count:'+str(s),len(good)==int(s['n_feasible_cells']))
  if good:
   worst=max(float(x['relative_range_p95']) for x in good);rr=sorted({float(x['range_km']) for x in good});near('summary_worst:'+str(s),s['worst_feasible_p95_relative_range'],worst)
   check('summary_ranges:'+str(s),s['feasible_ranges_km']==';'.join(f'{x:g}' for x in rr));full=rr==p['range_km']
   for k,v in [('all_three_ranges_feasible',full),('all_feasible_cells_strong',worst<=.01),('all_feasible_cells_usable',worst<=.05),('all_feasible_cells_project_level',worst<=.1),('full_range_usable',full and worst<=.05),('full_range_project_level',full and worst<=.1)]:check(k+str(s),flag(s[k])==v)
 for m in rows(OUT/'AUX_P95_RANGE_MAP.csv'):
  good=[x for x in data if x['geometry_status']=='FEASIBLE' and all(float(x[k])==float(m[k]) for k in ['range_km','baseline_km','alpha_deg','sigma_deg'])]
  check('map_count:'+str(m),len(good)==int(m['n_feasible']))
  if good:
   worst=max(float(x['relative_range_p95']) for x in good);near('map_worst:'+str(m),m['worst_p95_relative_range'],worst)
   check('map_5:'+str(m),flag(m['le5pct'])==(worst<=.05));check('map_10:'+str(m),flag(m['le10pct'])==(worst<=.1))
 d=load(OUT/'AUX_GATE0_DECISION.json');summaries=rows(OUT/'AUX_GEOMETRY_SUMMARY.csv');primary=[s for s in summaries if float(s['sigma_deg'])==.1]
 admitted=any(flag(s['full_range_usable']) and float(s['baseline_km'])<=p['practical_baseline_ceiling_km'] for s in primary)
 conditional=any(flag(s['all_feasible_cells_usable']) or flag(s['full_range_project_level']) for s in primary)
 outcome='AUXILIARY_PASSIVE_NODE_GEOMETRY_ADMITTED' if admitted else 'AUXILIARY_NODE_GEOMETRY_CONDITIONALLY_PRACTICAL' if conditional else 'SINGLE_AUXILIARY_BEARING_NODE_NOT_SUFFICIENT'
 check('route_rule',d['auxiliary_node_decision']==outcome);check('zero_credit',d['R4_progress_percent']==0);check('single_array_target_false',d['CURRENT_SINGLE_ARRAY_EVIDENCE_SUPPORTS_TARGET'] is False)
 start=load(OUT/'AUX_EXECUTION_START.json');check('freeze_before_execution',start['design_sha']==d['design_sha'] and start['remote_verified_before_geometry'] and start['policy_sha256']==sha(OUT/'AUX_GATE0_DESIGN_FREEZE.json'))
 return checks

def present():
 p=load(OUT/'AUX_GATE0_DESIGN_FREEZE.json');d=load(OUT/'AUX_GATE0_DECISION.json');data=rows(OUT/'AUX_GEOMETRY_RESULTS.csv');mapping=rows(OUT/'AUX_P95_RANGE_MAP.csv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figdir=OUT/'figures';figdir.mkdir(exist_ok=True);plt.rcParams.update({'svg.hashsalt':'AUX_GATE0_20261005','font.size':9})
 def save(fig,name):
  fig.savefig(figdir/(name+'.png'),dpi=150);fig.savefig(figdir/(name+'.svg'),metadata={'Date':None});plt.close(fig)
  svg=figdir/(name+'.svg');svg.write_text('\n'.join(s.rstrip() for s in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 for sigma in p['sigma_deg']:
  fig,axes=plt.subplots(1,3,figsize=(15,4.6));fig.subplots_adjust(top=.82,bottom=.19,left=.05,right=.87,wspace=.3)
  for ax,r in zip(axes,p['range_km']):
   vals=np.full((5,8),np.nan)
   for cell in mapping:
    if float(cell['range_km'])==r and float(cell['sigma_deg'])==sigma and cell['status']=='FEASIBLE':vals[p['baseline_km'].index(float(cell['baseline_km'])),p['alpha_deg'].index(float(cell['alpha_deg']))]=100*float(cell['worst_p95_relative_range'])
   cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#eeeeee');im=ax.imshow(vals,aspect='auto',cmap=cmap,vmin=0,vmax=30)
   for i in range(5):
    for j in range(8):ax.text(j,i,'X' if np.isnan(vals[i,j]) else f'{vals[i,j]:.1f}',ha='center',va='center',fontsize=8,color='black' if np.isnan(vals[i,j]) or vals[i,j]>15 else 'white')
   ax.set_xticks(range(8),p['alpha_deg']);ax.set_yticks(range(5),p['baseline_km']);ax.set_xlabel('Target crossing angle (deg)');ax.set_ylabel('Baseline (km)');ax.set_title(f'Range {r:g} km; worst root / mirror')
  cax=fig.add_axes([.91,.22,.015,.55]);fig.colorbar(im,cax=cax,label='P95 relative range error (%)');fig.suptitle(f'Geometry only; sigma={sigma:g} deg; X = infeasible triangle')
  save(fig,f'FIG_RANGE_MAP_SIGMA_{sigma:g}')
 fig,axes=plt.subplots(1,2,figsize=(12,4.5));fig.subplots_adjust(top=.80,bottom=.18,wspace=.3)
 for b in p['baseline_km']:
  rr=[x for x in data if x['geometry_status']=='FEASIBLE' and float(x['range_km'])==50 and float(x['baseline_km'])==b and float(x['sigma_deg'])==.1 and x['branch']=='FAR' and int(x['side'])==1];rr.sort(key=lambda x:float(x['alpha_deg']))
  if rr:
   angles=[float(x['alpha_deg']) for x in rr];axes[0].plot(angles,[100*float(x['relative_range_p95']) for x in rr],'o-',label=f'B={b:g} km');axes[1].plot(angles,[float(x['fim_condition']) for x in rr],'o-',label=f'B={b:g} km')
 axes[0].axhline(5,color='green',linestyle='--',label='5%');axes[0].axhline(10,color='orange',linestyle='--',label='10%');axes[0].set_ylabel('P95 relative range error (%)');axes[1].set_ylabel('Local FIM condition number');axes[1].set_yscale('log')
 for ax in axes:ax.set_xlabel('Crossing angle (deg)');ax.legend();ax.grid(alpha=.3)
 fig.suptitle('Near-parallel degradation; R=50 km; sigma=0.1 deg; FAR / upper; lines connect tested nodes only');save(fig,'FIG_NEAR_PARALLEL_CONDITIONING')
 fig,ax=plt.subplots(figsize=(7,5));fig.subplots_adjust(left=.16,bottom=.15)
 valid=[x for x in data if x['geometry_status']=='FEASIBLE'];xx=[100*float(x['linear_range_p95_relative']) for x in valid];yy=[100*float(x['relative_range_p95']) for x in valid]
 ax.scatter(xx,yy,s=12,alpha=.4);lim=max(max(xx),max(yy));ax.plot([0,lim],[0,lim],'k--');ax.set_xlabel('Local linearized 1.96 sigma / range (%)');ax.set_ylabel('Nonlinear empirical P95 range error (%)');ax.set_title('Local approximation versus forward-ray estimator');ax.grid(alpha=.3);save(fig,'FIG_LINEAR_VS_NONLINEAR')
 lines=['# Auxiliary-node geometry Gate0','',d['auxiliary_node_decision'],'','Geometry-only architecture diagnostic; no five-parameter or depth claim.','',f'Requested cells: {len(data)}; feasible: 240; impossible: 1680. Each feasible cell has {p["n_samples"]} Gaussian trials. Every physical branch and mirror retained; summaries take worst branch/mirror.','',
  '## Physical feasibility','', 'With MAIN=(0,0), target=(R,0), auxiliary-target distance d obeys d^2-2R cos(alpha)d+R^2-B^2=0. Feasible iff B>=R sin(alpha), d>0. Both roots and reflected placements are frozen. Impossible baseline/crossing combinations are not simulated. All 20-90 degree cells are impossible at the registered far ranges/baselines. No geometry added after freeze.','',
  '## Sigma 0.1 degree','', '| Baseline km | Best feasible-subset angle | Feasible ranges km | Worst P95 % | Best full 50/55/60 km P95 % |','|---|---|---|---|---|']
 for best in d['sigma_0p1_baseline_best']:
  b=best['baseline_km'];x=best['best_registered_feasible_region'];full=best['best_full_50_60km_region']
  lines.append(f'| {b:g} | '+(f'{x["alpha_deg"]:g} | {x["feasible_ranges_km"]} | {100*float(x["worst_feasible_p95_relative_range"]):.6f} | '+(f'{100*float(full["worst_feasible_p95_relative_range"]):.6f}' if full else 'N/A') if x else 'NONE | NONE | N/A | N/A')+' |')
 lines+=['','Full per-range regions (worst both roots/mirrors):','```json',json.dumps(d['primary_feasible_regions'],indent=2),'```','', 'Full-range <=5% regions: '+json.dumps(d['full_50_60km_usable_regions']), 'Full-range <=10% regions: '+json.dumps(d['full_50_60km_project_level_regions']),'',
 '## Statistical and engineering scope','',
 'Independent zero-mean Gaussian errors between sensors; shared frozen standard-normal draws across cells, mirror sign reversal for exact reflection checks. Cells are not independent trials. Directed target bearings, correct target association and common time are assumed; sensor positions are known exactly. Median/P90/P95/P99 are unconditional nearest-rank order statistics. Parallel, nonfinite and behind-sensor estimates receive infinite errors; finite ill-conditioned cases are retained. Ill-conditioning proxy: 1/abs(sin(observed crossing))>1000. Rates are separate and may overlap.',
 'Range is norm(position-MAIN); cross-range is absolute perpendicular position error to the true MAIN LoB; 2D position error is Euclidean target error. P95 empirical rank uncertainty uses a frozen approximate binomial-normal 95% order-statistic interval; Wilson intervals report fraction <=5%. Primary Gate is empirical P95, as frozen; these intervals do not give selected-grid simultaneous guarantees.',
 'H/FIM inverse is a local linearized lower-bound/approximation, not actual estimator guarantee. At near-parallel geometry nonlinear ratios can have long tails; Monte Carlo quantiles, including failures, govern the Gate. Known positions and unbiased bearings exclude positioning, attitude, shared bias, clock offset, latency, detection and signal-association error.',
 'Baseline 0.5/1 km has no feasible registered angle >=2 deg; N/A is no tested triangle, not proof those baselines are universally impossible. The frozen grid omits their smaller physically possible angles and does not cover every placement. No post-result angle additions. Similarly an infeasible 60 km cell cannot be replaced by a nearer target when claiming full-range performance.',
 'The 5 km practical ceiling is a planning assumption, not verified deployment capability. Hardware, independent directed-bearing accuracy, front/back resolution, self-localization/attitude, time/target association, communication, maneuver limits and acoustic detection remain future gates. An omnidirectional hydrophone alone cannot provide the assumed bearing; no HLA/VLA type is preselected.',
 '', '## Route boundaries','',
 'Two nonparallel directed bearings remove instantaneous horizontal scale ambiguity with known separated nodes. A temporal position sequence could support velocity/heading, but a single epoch cannot observe v/psi, and no tracker or motion accuracy is evaluated. Depth is absent from this model.',
 'Independent horizontal acquisition may justify future conditional-depth validation under a new architecture. Few-percent range is not a validated B1 depth-safe tolerance: B1A/B1B supplied no universally robust envelope. Current single-array depth stays closed, exact-horizontal mechanism retained as oracle evidence.',
 'Outcome A requires full-range <=5% at <=5 km, worst all roots/mirrors. Outcome B admits only conditional consideration (partial-range 5% or full-range 10%); architecture acceptance and hardware scale are research-lead decisions. No automatic R4-A1-NEW or depth restart.',
 '', 'No acoustic propagation/depth scores/global optimizer/SSP/5D estimator/tracker. R4=0%; A2/A3/A4/B2/B3/B4/C paused/not opened; P5 not opened. Commit/push/verify and stop.']
 report='\n'.join(lines)+'\n'
 for name in ['AUX_GATE0_REPORT.md','GPT_SYNC.md']:(OUT/name).write_text(report,encoding='utf-8')
 route='''# R4 route reset

Accepted B1B baseline 6b8d57ef6082f24720c9ece654b9cd710a687e13: CURRENT_CONDITIONAL_DEPTH_ENGINEERING_ROUTE_CLOSED. CZ_ENVELOPE_CERTIFIED_GLOBAL_SEARCH_ROUTE_CLOSED remains frozen. Cold-start single-array acquisition NOT ESTABLISHED. Exact-horizontal depth mechanism RETAINED_AS_ORACLE_EVIDENCE.

CURRENT_SINGLE_ARRAY_EVIDENCE_SUPPORTS_TARGET = false. Continuous bearing remains the most mature information source. It constrains theta and motion; noiseless full-rank singleton is retained as ideal/control evidence. Nominal r-v-psi compensation and wide range ridge remain. The original five-parameter joint <10% norm is not fully defined (relative theta at zero is undefined); no end-to-end range or engineering depth Gate has closed.

AUX1 adds one mobile passive node delivering known position and correctly associated directed target bearing, with no source cooperation or known waveform assumption. No auxiliary array type is preselected. Real bearing extraction, positioning/attitude, common-time registration, target association, detection and communication remain unverified.

AUX2 coarse external range prior must be independent, observation-derived and carry auditable uncertainty/truth-retention. No validated required width is known: 15 km noisy RC2 support did not admit capture; millimetre likelihood sections are not basin attraction widths; no tested B1A rectangle was universally robust. A truth-centered band or a prior produced by the failed global matcher is a hidden oracle/circular bootstrap.

AUX3 active/cooperative ranging: OUTSIDE_CURRENT_PASSIVE_SCOPE. Absolute range is an engineering reference, with changed mission/hardware and possibly source cooperation. AUX4 second passive node + TDOA: CONDITIONAL_ON_SYNCHRONIZATION_AND_SIGNAL_STRUCTURE. Clock alignment, source structure, multipath association and propagation mapping remain unverified; no TDOA simulation here.

Paired bearings over time may support motion structurally; no v/psi accuracy claim. Independent horizontal acquisition under a new architecture may justify separately authorized depth validation; B1 stays closed in the current architecture.

R4 original A/B sequence PAUSED_PENDING_ARCHITECTURE_RESET. R4=0%; A2/A3/A4/B2/B3/B4/C paused/not opened; P5 not opened. No automatic restart. Historical R4_PROGRESS.json remains a dated unchanged snapshot; this stage supplies current route status.

Geometry outcome: '''+d['auxiliary_node_decision']+'\nSee AUX_GATE0_REPORT.md for per-range limits.\n'
 (OUT/'R4_ROUTE_RESET_REPORT.md').write_text(route,encoding='utf-8')
 with Path('results/R4_MASTER/R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## Auxiliary-node geometry Gate0 execution\n\n'+d['auxiliary_node_decision']+'; PENDING_RESEARCH_LEAD_AUDIT. Known-position paired-directed-bearing geometry only. Physically impossible grid cells retained as infeasible. [Report](../R4_ROUTE_RESET_AUXILIARY_NODE_GATE0/AUX_GATE0_REPORT.md). R4=0%; no automatic architecture, motion or depth restart. Commit/push/verify then stop.\n')
 with Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow(['R4-AUX-GATE0','../R4_ROUTE_RESET_AUXILIARY_NODE_GATE0/AUX_GATE0_DECISION.json',d['auxiliary_node_decision'],'Registered paired-directed-bearing geometry only; no acoustic evaluations','PENDING_RESEARCH_LEAD_AUDIT; NO_SYSTEM_OR_DEPTH_CREDIT',0])

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--finalize',action='store_true');ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 if args.finalize:
  assert not (OUT/'VALIDATION.json').exists(),'No finalization rerun'
  checks=audit();present()
  with (OUT/'AUX_INDEPENDENT_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=['check','passed']);w.writeheader();w.writerows(checks)
  dump(OUT/'VALIDATION.json',dict(status='PASS',independent_checks=len(checks),failed=0,new_acoustic_propagation=0,new_depth_scores=0,feasible_cells=240,total_cells=1920,trials_per_feasible_cell=20000,tests='SEE_AUX_TESTS_XML',independent_method='Slope intersection and closed-form covariance; raw frozen Gaussian draws'))
  paths=[x for x in OUT.rglob('*') if x.is_file()]+[Path('results/R4_MASTER/R4_PLAN.md'),Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv')]
  dump(OUT/'AUX_RESULTS_FREEZE.json',dict(bindings={x.as_posix():sha(x) for x in sorted(paths)}))
  print(json.dumps(dict(finalized='PASS',checks=len(checks))))
 elif args.verify:
  seal=load(OUT/'AUX_RESULTS_FREEZE.json')
  for path,digest in seal['bindings'].items():assert sha(path)==digest,path
  checks=audit();print(json.dumps(dict(sealed='PASS',checks=len(checks),failed=0)))
 else:ap.error('--finalize or --verify required')
