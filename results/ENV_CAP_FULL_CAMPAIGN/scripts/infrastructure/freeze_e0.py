"""Prepare immutable scientific inputs and exact per-unit E1 plan before science."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'results/ENV_CAP_FULL_CAMPAIGN/scripts'))
from env_provider import *
import shutil,ast
I=O/'inputs';I.mkdir(exist_ok=True)
for name in ['SCS17_JW','SCS17_JA','SCS17_ANNUAL']:
 z,c=loadprofile(W/'prepared'/f'{name}.csv');z,c=cut(z,c,3909.75);table(I/(name+'.csv'),[{'depth_m':a,'sound_speed_mps':b} for a,b in zip(z,c)])
lines=(R/'results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/environments/n2_K2_f201.env').read_text().splitlines();dat=[]
for l in lines[5:]:
 if "'" in l:break
 a=l.split();dat.append((float(a[0]),float(a[1])))
assert dat[-1][0]==5000
for name,dat in [('ESTD_IDEAL',dat),('NO_CZ_CONSTANT',[(0.,1500.),(5000.,1500.)])]:table(I/(name+'.csv'),[{'depth_m':a,'sound_speed_mps':b} for a,b in dat])
case=[]
for env in ['SCS17_JW','SCS17_JA','SCS17_ANNUAL','ESTD_IDEAL','NO_CZ_CONSTANT']:
 profile=I/(env+'.csv');z,c=loadprofile(profile);H=z[-1]
 for zs in [150,200,250]:
  name=f'{env}_z{zs}';D=W/'E1'/name;ds=[]
  for lev,kind,b,step,angles in [('R0','R',201,25,(.1,20)),('R1','R',201,12.5,(.1,20)),('I0','I',2001,50,(-25,25)),('I1','I',4001,25,(-25,25)),('I2','I',8001,12.5,(-25,25))]:
   p=write_env(D,lev,z,c,H,zs,kind,b,step,angles);ds.append({'path':str(p),'sha256':sha(p)})
  case.append({'case_id':name,'env_id':env,'source_family':'WOA23_17N118E' if env.startswith('SCS17') else env,'profile':profile.relative_to(O).as_posix(),'profile_sha256':sha(profile),'H_m':H,'source_z_m':zs,'source_kind':'RECONSTRUCTED_CLIMATOLOGY_CONDITIONAL' if env.startswith('SCS17') else 'SYNTHETIC_CONTROL','input_sha':ds})
# Inputs and source copies preserve exact bytes.
S=I/'woa_raw';S.mkdir(exist_ok=True)
for p in (W/'sources/subsets').glob('*source17N118E.nc'):shutil.copy2(p,S/p.name)
shutil.copy2(W/'prepared/PROVENANCE.json',I/'WOA_RECONSTRUCTION_PROVENANCE.json')
shutil.copy2(W/'sources/M45_benchmark_selection/C1_ETOPO_POINT_AUDIT.json',I/'ETOPO_POINT_AUDIT.json')
shutil.copy2(W/'sources/independent_25N140E/ACCESS.json',I/'INDEPENDENT_SITE_ACCESS.json')
shutil.copy2(W/'sources/NETCDF_METADATA_REVIEW.json',I/'NETCDF_METADATA_REVIEW.json')
freeze={'parent_SHA':'a9a80191ee4d0f0cecd9db296198606b63a68080','branch':'research/env-cap-20261016','package_master_SHA256':sha(O/'task_package/ENV_CAP_FULL_CAMPAIGN_20261016_v1/01_CODEX_MASTER_PROMPT.md'),'source_selection_SHA256':sha(O/'E0_ACQUISITION_SELECTION_FREEZE.json'),'provider':'ENV_CAP_BELLHOP_CLINEAR_RAY_AND_INCOHERENT_ENERGY_V1','E1_frequency_hz':200,'range_domain_km':[10,120],'receiver_depths_m':{'start':100,'stop':300,'step':2.5},'range_step_km':.125,'TL_smoothing':'9-point moving mean LINEAR ENERGY; source depth +/-25m receiver mean; not coherent modal pressure','surface':'pressure release','bottom':'rigid NOMINAL MATHEMATICAL conditional; only source-site ETOPO depth held flat, no real transect claim','SSP_interpolation':'C-linear, no extrapolation, all4-corner finite required','new_solver_max_calls':78,'new_MC_max_configs':0,'per_call_timeout_s':120,'E1_process_timeout_s':5400,'total_science_max_s':259200,'tree_RSS_max_GiB':8,'new_D_max_GiB':30,'new_Git_max_MiB':120,'retry_after_science_start':0,'E1_cases':case,'E1_gates':{'A':'First interior local minimum of first bottom/surface-free deep return curve in 0.1..20deg, 2001-angle analytic scan; earliest range candidate within10..120km; no fold -> noCZ2','B':'Compare 201 Bellhop rays at25m/12.5m integration step to Snell exact C-linear return; >=10 returns, max<=0.1km','C':'I0/I1/I2 as saved envs; both adjacent comparisons median<=1dB,P95<=2dB over fold+/-5km, peakshift<=1km,widthshift<=2km','D':'Fine TL local minimum in fold+/-5km interior, prominence over5..10km flank >max(0.5dB,2*largestP95numericdifference). Not a universalCZcontrast rule.'},'not_evaluated_sources':{'ENV_SCS_WOA_18N114E':'Cross-study diagnostic deep completion and old NCSS vertical remapping; no new certified full-depth original source SSP; retain diagnostic only','ENV_INDEPENDENT_REGIONAL_25N140E':'Fixed raw T/S requests HTTP400, ETOPO obtained; no site substitution or invented SSP'},'E2_preregistered_units':['Archived shallow Yang 4-case control re-evaluation, single-mode integration/Erratum phase/unknown f0/mirror/free modal-gain controls. No new pressure truth from uncertified modes.','CZ transfer admission requires CZ2 AND L3 certified complete coherent modal provider. This campaign E1 provider is ray/incoherent intensity only: if no existing L3 provider with current source/window/bottom identity, all CZ extraction/depth performance NOT_EVALUATED, complete analytic resource branch.','MMAR original AR k/Hankel amplitude separation audit; archived p7,20dB peak-identity partial result retained; no replacing missing paper choices with tuned algorithm.'],'E3_preregistered_units':['Only if E2 L5 exists: separately freeze finite sensitivity matrix BEFORE execution. If no L3/L5, actual performance rows all NOT_EVALUATED with reason; L0 zero controls and analytic coherence/aperture resource sensitivities remain permitted.'],'E4_preregistered_units':['L0 resource curves for ordinary Fourier delta_k=2pi/Lr; required source/clock frequency error <= f*sigma_q/c; offset phase differential <=0.1rad illustrative across dk. No algorithm-independent bound claim.','Read-only inversion of T3_W1 complete frozen statistics and resources, labelled historical signed-bearing/Uniform-bias controls; no re-run or depth credit.','No new feasible CZ accuracy rows without forward/extraction/holdout gates; all unsupported grade/resource entries UNKNOWN_OR_UNTESTED.'],'precision_ladders_NOT_client_requirement':{'depth_m':[5,10,20],'range_percent':[1,5,10],'speed_mps':[.1,.25,.5],'bearing_deg':[.2,.5,1],'heading_deg':[.2,.5,1],'coverage':[.9,.95]},'E5':'Independent cold audit once from saved raw inputs/outputs; never import science generator in cold implementation; audit failures preserve scientific status and block affirmative claims. Complete all final outputs and STOP.'}
js(O/'E0_FREEZE.json',freeze)
for p in (O/'scripts').glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
js(O/'E0_CODE_FREEZE.json',{str(p.relative_to(O)):sha(p) for p in (O/'scripts').glob('*.py')})
print('Prepared',len(case),'units / max75 E1 calls. No E1 science run.')
