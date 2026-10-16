"""Registered E4 analytic design scales and read-only historical resource inversion."""
from env_provider import *
import struct

def main():
 D=W/'E4';D.mkdir(exist_ok=True);marker=D/'SCIENCE_STARTED.json'
 if marker.exists():raise RuntimeError('E4 already executed')
 js(marker,{'code_sha256':sha(__file__),'config_sha256':sha(O/'E4_UNIT_FREEZE.json'),'seed':'DETERMINISTIC_NO_NEW_MC','scope':'Analytic design coordinates + historical bearing-only conditional evidence; no CZ estimator performance'})
 f=json.loads((O/'E4_UNIT_FREEZE.json').read_text());assert all(sha(R/a['path'])==a['sha256'] for a in f['inputs'])
 scales=[]
 for gap in [1e-5,3e-5,1e-4,3e-4,.001,.003,.01]:
  L=2*np.pi/gap
  for q in [.1,.25,.5,1.,2.]:
   scales.append({'kind':'FOURIER_RAYLEIGH_SCALE','delta_k_rad_per_m':gap,'radial_aperture_m':L,'relative_radial_speed_mps':q,'required_time_s':L/q,'condition':'monotone radial displacement; known coordinates and phase-coherent source','evidence_level':'L0_ANALYTIC','status':'NECESSARY_ANALYTIC_BOUND','scope':'ordinary Fourier adjacent-line Rayleigh separation criterion ONLY, not every high-resolution estimator, not actual current CZ mode spacing'})
 for freq in [150.,200.,350.]:
  for sigq in [.1,.25,.5]:
   scales.append({'kind':'SINGLE_COMPONENT_REFERENCE_FREQUENCY_SCALE','frequency_hz':freq,'desired_q_error_mps':sigq,'max_reference_error_hz':freq*sigq/1500,'condition':'one effective propagation speed1500m/s and signed instantaneous radial rate; multimode/window mapping unvalidated','evidence_level':'L0_ANALYTIC','status':'NECESSARY_ANALYTIC_BOUND','scope':'first-order kinematic frequency-reference sensitivity only'})
 for T in [200.,1200.,2400.,4800.]:
  scales.append({'kind':'CW_PHASE_COHERENCE_DESIGN_SCALE','observation_s':T,'phase_budget_rad':.1,'max_uncorrected_frequency_error_hz':.1/(2*np.pi*T),'evidence_level':'L0_ANALYTIC','status':'NECESSARY_ANALYTIC_BOUND','scope':'illustrative coherent phase budget, not measured source stability or universal method precision'})
 for gap in [1e-5,1e-4,.001,.01]:
  scales.append({'kind':'RELATIVE_MODAL_OFFSET_PHASE_SCALE','delta_k_rad_per_m':gap,'relative_phase_budget_rad':.1,'max_delta_error_m':.1/gap,'evidence_level':'L0_ANALYTIC','status':'NECESSARY_ANALYTIC_BOUND','scope':'relative modal phase at frozen illustrative0.1rad tolerance; common phase may be absorbed by source'})
 table(O/'E4_ANALYTIC_DESIGN_SCALES.csv',scales)
 # Existing shallow modes are a NON-CZ resolution reference only.
 p=R/'results/R3_C2_Yang_SA_depth/R3_C2_2C/_kraken_paper/yang2014_f350.mod';buf=p.read_bytes();rec=4*struct.unpack_from('<i',buf)[0];M=struct.unpack_from('<i',buf,5*rec)[0];k=np.sort(np.frombuffer(buf,dtype='<c8',count=M,offset=(7+M)*rec).real.astype(float));gap=np.diff(k);near=np.minimum(np.r_[np.inf,gap],np.r_[gap,np.inf]);modes=[]
 for L in [14.,600.,1200.,2400.,4990.,10000.]:
  modes.append({'env_id':'YANG_SHALLOW_NON_CZ_H88','frequency_hz':350,'radial_aperture_m':L,'total_modes':M,'energetic_modes':'NOT_EVALUATED','adjacent_gap_min':gap.min(),'adjacent_gap_median':np.median(gap),'rayleigh_unique_all_modes':int(sum(near>=2*np.pi/L)),'actually_extracted_modes':'NOT_EVALUATED','status':'ARCHIVED_MODAL_GEOMETRY_DIAGNOSTIC_ONLY','source_sha256':sha(p)})
 for env in ['SCS17_JW','SCS17_JA','SCS17_ANNUAL','ESTD_IDEAL']:
  modes.append({'env_id':env,'frequency_hz':200,'total_modes':'NOT_EVALUATED','energetic_modes':'NOT_EVALUATED','adjacent_gap_min':'NOT_EVALUATED','rayleigh_unique_all_modes':'NOT_EVALUATED','actually_extracted_modes':'NOT_EVALUATED','status':'NO_CERTIFIED_COMPLETE_MODAL_PROVIDER_FOR_THIS_ENV_IDENTITY'})
 table(O/'E2_MODE_DENSITY_AND_APERTURE.csv',modes)
 hist=R/'results/R4_T3_W1_WIDE_BASELINE_DEMO';stats=list(csv.DictReader((hist/'FOUR_PARAMETER_DEMO_TABLE.csv').open()));resource={a['geometry']:a for a in csv.DictReader((hist/'RESOURCE_AND_COVERAGE_SUMMARY.csv').open())};inv=[];pareto=[]
 grade={'range':([1,5,10],'range_worst_case_P95_percent','percent'),'speed':([.1,.25,.5],'speed_abs_worst_case_P95_mps','m/s'),'bearing':([.2,.5,1],'bearing_worst_case_P95_deg','degree'),'heading':([.2,.5,1],'heading_worst_case_P95_deg','degree')}
 for s in stats:
  r=resource[s['geometry']];q=s['quality'];row={'config':q+'_'+s['geometry'],'scope':'HISTORICAL_T3_SIGNED_BEARING_ONLY_NOT_CZ_ACOUSTIC_EXTRACTION','receivers':2,'element_count':'UNKNOWN','node_separation_km':float(r['node_separation_km']),'observation_s':1200,'bearing_sigma_deg':.05 if q=='Q1' else .075,'nav_sigma_axis_m':25 if q=='Q1' else 37.5,'common_bias_bound_deg':.05 if q=='Q1' else .075,'differential_half_bound_deg':.025 if q=='Q1' else .0375,'independent_phase_sync':'NOT_REQUIRED_BEARING_LEVEL','sync_and_association':'ASSUMED_NOT_VALIDATED','bias_at_bound_max_rate':float(s['max_case_bias_boundary_rate']),'actual_hardware':'UNKNOWN','point_errors_not_candidate_coverage':True}
  for param,(_,key,unit) in grade.items():row[param+'_worst_case_P95']=float(s[key])
  pareto.append(row)
  for param,(levels,key,unit) in grade.items():
   for level in levels:
    val=float(s[key]);inv.append({'cz_env_id':'NO_ACOUSTIC_ENV_T3_KINEMATIC_MODEL','desired_ability':param,'illustrative_error_level':level,'unit':unit,'required_percentile':.95,'observed_worst_case_P95':val,'condition':'HISTORICAL_SIGNED_BEARING_Uniform_bias_support_kinematic_model_and_bounded_estimator','status':'TESTED_CONDITIONAL_FEASIBLE' if val<=level else 'TESTED_INFEASIBLE_FOR_METHOD','tested_config_reference':'T3_W1_'+q+'_'+s['geometry'],'additional_nodes':1,'receiver_elements':'UNKNOWN','baseline_km':r['node_separation_km'],'time_s':1200,'why_not_new_CZ_precision':'No new signal extraction or CZ detectability validation; empirical error quantile is not candidate-set coverage or performance guarantee','engineering_unknowns':'DOA extraction; calibration; navigation; side; synchronization; target association; communication','evidence':'HISTORICAL_CONDITIONAL_MC_REPORTED_NOT_RERUN'})
 # Non-dominance includes performance AND unweighted physical resource coordinates; quality requirement is a cost.
 for a in pareto:
  av=np.array([a['node_separation_km'],1/a['bearing_sigma_deg'],1/a['nav_sigma_axis_m'],1/a['common_bias_bound_deg']]+[a[p+'_worst_case_P95'] for p in grade]);dominators=[]
  for b in pareto:
   if b is a:continue
   bv=np.array([b['node_separation_km'],1/b['bearing_sigma_deg'],1/b['nav_sigma_axis_m'],1/b['common_bias_bound_deg']]+[b[p+'_worst_case_P95'] for p in grade])
   if np.all(bv<=av) and np.any(bv<av):dominators.append(b['config'])
  a['pareto_nondominated']=not dominators;a['dominated_by']=';'.join(dominators);a['pareto_scope']='Finite8 historical configurations, distinct quality contracts costed via inverse sigma/bias support; not global hardware minimum'
 for env in ['SCS17_JW','SCS17_JA','SCS17_ANNUAL','ESTD_IDEAL','ENV_INDEPENDENT_REGIONAL_25N140E','ENV_SCS_WOA_18N114E']:
  for param,levels in [('depth',[5,10,20]),('range',[1,5,10]),('speed',[.1,.25,.5]),('bearing',[.2,.5,1]),('heading',[.2,.5,1])]:
   for level in levels:
    inv.append({'cz_env_id':env,'desired_ability':param,'illustrative_error_level':level,'status':'PHYSICAL_ENV_NOT_ADMITTED' if env.startswith('ENV_') else 'UNKNOWN_OR_UNTESTED','observed_worst_case_P95':'NOT_EVALUATED','tested_config_reference':'NONE','why_not_new_CZ_precision':'No admitted coherent modal-pressure extraction and holdout estimator chain; true coordinate, mode ID, f0 and delta cannot be assumed as observations','engineering_unknowns':'source coherence/frequency; observability of relative radial coordinates; SSP and bottom uncertainty; calibration; branch coverage','evidence':'NO_FEASIBILITY_CLAIM'})
 table(O/'E4_RESOURCE_INVERSION.csv',inv);table(O/'E4_PARETO_RESOURCE_FEASIBILITY.csv',pareto)
 js(O/'E4_GATE.json',{'status':'NECESSARY_SCALES_AND_HISTORICAL_CONDITIONAL_RESOURCE_MAP_ONLY','new_CZ_tested_feasible_rows':0,'new_CZ_precision':'NOT_EVALUATED','minimum_equipment_proven':False,'no_shared_noise_double_count':True,'Pareto_scope':'eight historical kinematic configurations, not new CZ acoustic feasibility','new_MC':0,'science_runs':1})
 print('E4 analytic coordinates and historical resource map completed; no new CZ performance claimed.')
if __name__=='__main__':main()
