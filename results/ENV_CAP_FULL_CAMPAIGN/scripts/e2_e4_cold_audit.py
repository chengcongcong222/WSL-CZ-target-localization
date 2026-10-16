"""Independent E2/E4 cold audit; no main-science imports, no optimizer/solver."""
from pathlib import Path
import json,csv,struct,math,hashlib,numpy as np
R=Path(__file__).resolve().parents[3];O=R/'results/ENV_CAP_FULL_CAMPAIGN';W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP')
def main():
 out=O/'E2_E4_INDEPENDENT_COLD_AUDIT.json'
 if out.exists():raise RuntimeError('Independent cold audit already exists')
 checks=[]
 def ck(name,ok,detail='',scope='identity_or_readonly_reconstruction'):
  checks.append({'check':name,'PASS':bool(ok),'detail':str(detail),'scope':scope})
 def rows(name):return list(csv.DictReader((O/name).open(encoding='utf-8-sig')))
 for f in ['E2_UNIT_FREEZE.json','E4_UNIT_FREEZE.json']:
  for a in json.loads((O/f).read_text())['inputs']:ck(f+'_'+a['path'],hashlib.sha256((R/a['path']).read_bytes()).hexdigest()==a['sha256'])
 p=R/'results/R3_C2_Yang_SA_depth/R3_C2_2C';buf=(p/'_kraken_paper/yang2014_f350.mod').read_bytes();rec=4*struct.unpack_from('<i',buf)[0];nz=struct.unpack_from('<i',buf,92)[0];M=struct.unpack_from('<i',buf,5*rec)[0]
 phis=[]
 for m in range(M):
  vals=struct.unpack_from('<'+str(nz*2)+'f',buf,(7+m)*rec);phis.append([complex(vals[2*j],vals[2*j+1]) for j in range(nz)])
 gtable=list(csv.DictReader((p/'paper_mode_identification.csv').open()));curves=rows('E2_PAPER_DEPTH_CONTROL_CURVES.csv')
 for zs in [4,50]:
  for zr in [18,70]:
   gt=sorted([a for a in gtable if float(a['zs'])==zs and float(a['zr'])==zr and a['method']=='ORACLE_MODE_ID'],key=lambda a:int(a['mode_id']));gm=[complex(float(a['g_real']),float(a['g_imag'])) for a in gt]
   for reg in [.05,.1,.2]:
    delta=reg*max(abs(v) for ph in phis for v in ph);score=[]
    for z in range(nz):score.append(abs(sum(phis[m][z]*gm[m]*phis[m][zr]/(phis[m][zr]**2+delta**2) for m in range(M)))**2)
    total=sum(score);score=[v/total for v in score];saved=sorted([a for a in curves if float(a['zs'])==zs and float(a['zr'])==zr and float(a['regularization_ratio'])==reg],key=lambda a:float(a['z_candidate_m']));err=max(abs(a-float(b['score'])) for a,b in zip(score,saved));ck(f'shallow_depth_{zs}_{zr}_{reg}',len(saved)==nz and err<1e-10,err,'stored_shallow_oracle_only')
 # Python integer identity shows source of NaN rather than pretending main zero control passed.
 observed=rows('E2_ANALYTIC_AND_PAPER_CONTROLS.csv');mirror=next(a for a in observed if a['check']=='HLA_left_right_mirror_zero');overflow=[(50000-(2*i-7))**2 for i in range(8)]
 ck('int32_overflow_failure_preserved',mirror['PASS']=='False' and mirror['error_m']=='nan' and all(v>2**31-1 for v in overflow),'Exact integer square exceeds2147483647; mirror analytic equality still exact','negative_implementation_attribution_not_scientific_gate_repair')
 ck('E2_25_PASS_1_FAIL_accounting',sum(a['PASS']=='True' for a in observed)==25 and sum(a['PASS']=='False' for a in observed)==1)
 # Shared-error differencing: full precision cannot exceed raw information.
 A=np.array([[-1. if j==i else (1. if j==i+1 else 0.) for j in range(6)] for i in range(5)]);j=np.arange(6,dtype=float)**2;cov=A@A.T;direct=float(j@j);trans=float((A@j)@np.linalg.solve(cov,A@j));saved=next(a for a in observed if a['check']=='shared_window_noise_and_DPI');ck('DPI_and_shared_covariance',np.all(np.diag(cov)==2) and np.all(np.diag(cov,1)==-1) and trans<=direct and abs(trans-float(saved['difference_information']))<1e-10)
 for a in rows('E4_ANALYTIC_DESIGN_SCALES.csv'):
  k=a['kind']
  if k=='FOURIER_RAYLEIGH_SCALE':err=abs(float(a['radial_aperture_m'])*float(a['delta_k_rad_per_m'])-2*math.pi);ok=err<1e-10 and abs(float(a['required_time_s'])*float(a['relative_radial_speed_mps'])-float(a['radial_aperture_m']))<1e-8
  elif k=='SINGLE_COMPONENT_REFERENCE_FREQUENCY_SCALE':err=abs(float(a['max_reference_error_hz'])*1500/float(a['frequency_hz'])-float(a['desired_q_error_mps']));ok=err<1e-12
  elif k=='CW_PHASE_COHERENCE_DESIGN_SCALE':err=abs(2*math.pi*float(a['max_uncorrected_frequency_error_hz'])*float(a['observation_s'])-.1);ok=err<1e-12
  else:err=abs(float(a['max_delta_error_m'])*float(a['delta_k_rad_per_m'])-.1);ok=err<1e-12
  ck('analytic_scale_'+k,ok,err,'L0_scale_not_achieved_precision')
 hist=R/'results/R4_T3_W1_WIDE_BASELINE_DEMO';cstats=list(csv.DictReader((hist/'CASE_STATISTICS.csv').open()));pstats=list(csv.DictReader((hist/'FOUR_PARAMETER_DEMO_TABLE.csv').open()));mapping={'range_worst_case_P95_percent':('range_rel_unconditional_P95',100),'bearing_worst_case_P95_deg':('bearing_deg_unconditional_P95',1),'speed_abs_worst_case_P95_mps':('speed_abs_mps_unconditional_P95',1),'heading_worst_case_P95_deg':('heading_deg_unconditional_P95',1)}
 for a in pstats:
  sub=[s for s in cstats if s['quality']==a['quality'] and s['geometry']==a['geometry'] and s['side_group']=='BOTH']
  ck('historical_denominators_'+a['quality']+a['geometry'],len(sub)==12 and sum(int(b['n']) for b in sub)==6000 and sum(int(b['n_failed']) for b in sub)==0,'Historical summaries checked; raw48000 estimations not rerun','historical_readonly_statistics')
  for key,(col,mult) in mapping.items():
   maximum=max(float(b[col])*mult for b in sub);ck('historical_quantile_max_'+key+a['quality']+a['geometry'],abs(maximum-float(a[key]))<1e-10,maximum,'historical_readonly_statistics')
 p=rows('E4_PARETO_RESOURCE_FEASIBILITY.csv')
 for a in p:
  def vec(b):return [float(b['node_separation_km']),1/float(b['bearing_sigma_deg']),1/float(b['nav_sigma_axis_m']),1/float(b['common_bias_bound_deg'])]+[float(b[k+'_worst_case_P95']) for k in ['range','speed','bearing','heading']]
  v=vec(a);dom=[]
  for b in p:
   if b['config']==a['config']:continue
   w=vec(b)
   if all(x<=y for x,y in zip(w,v)) and any(x<y for x,y in zip(w,v)):dom.append(b['config'])
  ck('finite_resource_Pareto_'+a['config'],(not dom)==(a['pareto_nondominated']=='True') and ';'.join(dom)==a['dominated_by'],'No global optimal hardware claim','historical_finite8_resource_map')
 ck('no_new_CZ_feasibility_claims',all(a['status'] not in ['TESTED_CONDITIONAL_FEASIBLE','TESTED_INFEASIBLE_FOR_METHOD'] for a in rows('E4_RESOURCE_INVERSION.csv')))
 ck('all30_performance_cells_masked',len(rows('E3_PERFORMANCE_BY_CASE.csv'))==30 and all(a['status']=='NOT_EVALUATED' and int(a['n_attempted'])==0 for a in rows('E3_PERFORMANCE_BY_CASE.csv')))
 result={'readonly_reconstruction_PASS':all(a['PASS'] for a in checks),'checks':len(checks),'FAIL_count':sum(not a['PASS'] for a in checks),'E2_scientific_implementation_PASS':False,'E1_cold_PASS':False,'new_CZ_scientific_admission':'NOT_ESTABLISHED','scope':'Confirms preserved negative gates and independent archived-paper/analytic/historical-resource calculations. PASS here never repairs scientific failures.','rows':checks}
 out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k!='rows'}));return 0 if result['readonly_reconstruction_PASS'] else 2
if __name__=='__main__':raise SystemExit(main())
