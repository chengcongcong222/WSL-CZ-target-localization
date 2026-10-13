"""N2 input/model whitelist preparation and eight algebra controls, no solver before Commit A."""
from t2_n2_common import *
import shutil,subprocess,sys,ast,platform,re

def prepare():
 if OUT.exists() or (LOCAL/'PREPARED.json').exists():raise RuntimeError('N2 directory already prepared')
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==PARENT
 assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==PARENT
 OUT.mkdir();LOCAL.mkdir(exist_ok=True);(OUT/'environments').mkdir();(OUT/'fields_parameters').mkdir();manifest=read(PKG/'T2_N2_FILE_MANIFEST_20261013_v1.json')
 for r in manifest['files']:
  p=PKG/r['name'];assert sha(p)==r['sha256'] and p.stat().st_size==r['bytes'];shutil.copyfile(p,OUT/p.name)
 shutil.copyfile(PKG/'T2_N2_FILE_MANIFEST_20261013_v1.json',OUT/'T2_N2_FILE_MANIFEST_20261013_v1.json');shutil.copyfile(G0/'ARRAY_GEOMETRY_TABLE.csv',OUT/'ARRAY_GEOMETRY_TABLE.csv');(OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8')
 source_names=[(G0,n) for n in ['GPT_SYNC.md','DESIGN_FREEZE.json','PROVIDER_PROVENANCE.json','PROVIDER_MODE_WINDOW_AND_MESH.csv','FORWARD_DECISION.json','FORWARD_CATALOGUE_NUMERICAL_STABILITY.csv','ARRAY_GEOMETRY_TABLE.csv','CATALOGUE_primary.npy','CATALOGUE_fine.npy','LOCAL_ARCHIVE_INDEX.csv','REPLICA_FREEZE.json']]+[(N1,n) for n in ['GPT_SYNC.md','NEXT_RESEARCH_DECISION.md','TOP_ERROR_CANDIDATES.csv','PAIR_ADVANTAGE_NUMERIC_INTERVALS.csv','ANALYSIS_FREEZE.json']]+[(ROOT/'results/R4_H3_P1_UNIFIED_WINDOW_SSP',n) for n in ['GPT_SYNC.md','FAILURE_SCOPE.md']];sources=[]
 for folder,name in source_names:
  p=folder/name;rel=str(p.relative_to(ROOT)).replace('\\','/');blob=subprocess.check_output(['git','show',PARENT+':'+rel]);sources.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'parent_blob_sha256':hashlib.sha256(blob).hexdigest(),'role':'historical_read_only_source'})
 oldmodes=rows(G0/'PROVIDER_MODE_WINDOW_AND_MESH.csv');index={r['path']:r for r in rows(G0/'LOCAL_ARCHIVE_INDEX.csv')};freeze=read(G0/'REPLICA_FREEZE.json');binary=[]
 for name in ['kraken','scooter','fields','field']:
  p=AT/(name+'.exe');binary.append({'name':name,'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'executable_access':os.access(p,os.X_OK)})
 for r in binary:sources.append(dict(r,role='frozen_binary'))
 for f in FREQ:
  for suffix in ['mod','env']:
   p=OLDLOCAL/'modal'/f't2_f{f}_n160001.{suffix}';h=sha(p);assert h==index[str(p)]['SHA256'];m=next(r for r in oldmodes if r['mesh']=='160001' and int(float(r['frequency_hz']))==f);assert h==m['MOD_SHA256'] if suffix=='mod' else h==m['ENV_SHA256'];sources.append({'path':str(p),'sha256':h,'bytes':p.stat().st_size,'role':'old_K0_'+suffix})
 for cat in CATS:
  for res in RES:
   for mesh in [80001,160001]:
    name=f'{res}_{cat}_n{mesh}.npy';p=OLDLOCAL/'replicas'/name;h=sha(p);assert h==index[str(p)]['SHA256']==freeze['files'][name];a=np.load(p,mmap_mode='r',allow_pickle=False);assert a.shape==(3355 if cat=='primary' else 12705,3,12 if res.startswith('B') else 8);sources.append({'path':str(p),'sha256':h,'bytes':p.stat().st_size,'role':'old_immutable_complex_archive','shape':list(a.shape),'dtype':str(a.dtype)})
 for p in (LOCAL/'source_review').iterdir():
  if p.is_file():sources.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'role':'official_source_static_reference_not_binary_build_proof'})
 dump(OUT/'SOURCE_IDENTITY.json',{'parent_SHA':PARENT,'task_SHA256':TASKSHA,'records':sources,'old_sources_unchanged':True})
 g=geometries();receiver=np.unique(np.concatenate([a[:,2] for a in g.values()]));membership=[]
 for s in sentinels():
  coord=[float(s[k]) for k in ['r_km','abs_theta_deg','source_z_m']];entry={'sentinel_id':s['sentinel_id'],'coordinates':coord}
  for cat in CATS:
   ids=np.flatnonzero(np.all(states(cat)==coord,axis=1));entry[cat+'_candidate_id']=int(ids[0]) if len(ids)==1 else 'NOT_SAMPLED';assert len(ids)==1 if cat=='fine' else len(ids)<=1
  membership.append(entry)
 dump(OUT/'SENTINEL_CATALOGUE_MEMBERSHIP.json',membership)
 plans=[];changes=[]
 for family,variants in [('KRAKEN',{k:v for k,v in KVAR.items() if k!='K0'}),('SCOOTER',FVAR)]:
  for variant,(mesh,high) in variants.items():
   for f in FREQ:
    base=OLDLOCAL/'modal'/f't2_f{f}_n160001.env';original=base.read_text(encoding='utf-8-sig').splitlines();b=next(i for i,l in enumerate(original) if l.strip().startswith("'R'"));new=original.copy();new[0]=f"'T2_N2_{variant}'";new[4]=f'{mesh} 0.0 5000.0';new[b+1]=f'1500.0 {high}.0'
    if family=='SCOOTER':new=new[:b+2]+['65.0','5']+[f'{z:.12f}' for z in [180,190,200,210,220]]+[str(len(receiver))]+[f'{z:.12f}' for z in receiver]+['R']
    p=OUT/'environments'/f'n2_{variant}_f{f}.env';p.write_text('\n'.join(new)+'\n',encoding='utf-8',newline='\n');assert new[1:4]==original[1:4] and new[5:b+1]==original[5:b+1];allowed={0,4,b+1} if family=='KRAKEN' else {0,4,*range(b+1,max(len(new),len(original)))}
    diff=[i for i in range(max(len(new),len(original))) if (new[i] if i<len(new) else None)!=(original[i] if i<len(original) else None)];assert all(i in allowed for i in diff)
    plans.append({'family':family,'variant':variant,'frequency_hz':f,'mesh':mesh,'CLOW':1500,'CHIGH':high,'RMAX':60 if family=='KRAKEN' else 65,'stem':p.stem,'env_path':str(p),'env_sha256':sha(p),'timeout_s':180 if family=='KRAKEN' else 900})
    changes.append({'variant':variant,'frequency_hz':f,'original_env':str(base),'changed_zero_based_lines':json.dumps(diff),'physical_SSP_and_boundary_bytes_unchanged':True,'source_receiver_change':'NONE_K_retains_all_original_depths' if family=='KRAKEN' else 'F_source5_depths_receiver_exact4resource_union'})
 fieldplans=[]
 for variant in FVAR:
  for f in FREQ:
   for s in sentinels():
    for step in ([.05,.025] if variant=='F2' else [.05]):
     tag='d050' if step==.05 else 'd025';stem=f"n2_{variant}_f{f}_{s['sentinel_id']}_{tag}";p=OUT/'fields_parameters'/(stem+'.flp');nr=int(round(40/step))+1;lo=float(s['r_km'])-.020;hi=float(s['r_km'])+.020;content=f"'RPO'\n{nr}\n{lo:.12f} {hi:.12f} /\n";p.write_text(content,encoding='utf-8',newline='\n');fieldplans.append({'variant':variant,'frequency_hz':f,'sentinel_id':s['sentinel_id'],'range_step_m':step,'request_count':nr,'range_min_km':lo,'range_max_km':hi,'stem':stem,'flp_path':str(p),'flp_sha256':sha(p),'timeout_s':120})
 table(OUT/'ENVIRONMENT_WHITELIST_AUDIT.csv',changes);dump(OUT/'SOLVER_PLAN.json',{'primary':plans,'FIELDS':fieldplans,'counts':{'KRAKEN':9,'SCOOTER':9,'FIELDS':144}})
 strings=(AT/'fields.exe').read_bytes();fields_source=(LOCAL/'source_review/fields.f90').read_text(encoding='utf-8-sig');source_positions=(LOCAL/'source_review/SourceReceiverPositions.f90').read_text(encoding='utf-8-sig');ready=all(r['executable_access'] for r in binary if r['name'] in ['kraken','scooter','fields']) and 'CALL ReadRcvrRanges' in fields_source and 'ReadVector( Pos%NRr' in source_positions and b'Option(2:2) should select positive' in strings
 dump(OUT/'SOURCE_RECORDS.json',{'official_manuals':read(OUT/'T2_N2_EXPERIMENT_SPEC_20261013_v1.json')['sources'],'official_source_downloads':read(LOCAL/'source_review/SOURCE_DOWNLOADS.json'),'FIELDS_format':'RPO; then integer NR; then RMIN RMAX slash; ReadRcvrRanges/ReadVector current interface, not1997 RDB/RMIN RMAX NR','FFT_convention':'positive-wavenumber cylindrical, exp(-ikR), frozen polynomial G interpolation; no sourcebeam','same_TopOpt_CVWT_boundary_R_SSP_density':True,'binary_source_exact_build_not_claimed':True,'static_interface_ready_to_attempt':bool(ready),'runtime_ready_not_pretested':True})
 dump(OUT/'EXECUTION_CONTRACT.json',{'static_FFP_ready':bool(ready),'FIELDS_options':'RPO','local_interpolation':'np.interp real and imaginary separately, actual SHD float32 saved axis only, no extrapolation','depth_absolute_tolerance_m':1e-5,'frequency_absolute_tolerance_hz':1e-6,'shape_threshold_K_and_FFP':.002,'FFP_interpolation_threshold':.0005,'cold_relative_error_tolerance':1e-10,'raw_reference':'B field in shape(a,b); raw norm(b-a)/norm(b); crossmethod global scalar/f/source-depth acrossall12sentinelsand4resources, each pressure entry retained','K0_field_source':'qualified original high-grid full-catalogue complex archives, no recomputation of old solver','modal_band_diagnostic':'K1 cp>1800, K2 cp>2000, auxiliary only; fullreturnedmode arrays unchanged','K_failure':'stop subsequentK calls at first failure; preserve completed; FFP still attempts frozen table ifready/budget','FFP_failure':'stop independentFFP path at first nonzero solver or invalidGRN/SHD; remainingregisteredcalls NOT_EVALUATED, no parameterchangeor retry','source_receiver_K_depth_lists':'all originalG0 depths retained per Kwhitelist, includesrequired4resourceunion; unusedregistereddepthsnot newgeometry','new_waveform_MC_DFT_FIM_optimizer_or_template_score':0,'budget':read(OUT/'T2_N2_EXPERIMENT_SPEC_20261013_v1.json')['budget'],'historical_R4_percent':0,'no_next_stage_authorization':True})
 dump(LOCAL/'PREPARED.json',{'utc':utc(),'no_solver_executed':True});print('PREPARED 9K/9S/144FIELDS, staticFFPready',ready,flush=True)

def preflight():
 checks=[]
 def add(name,e,t=1e-10):checks.append({'control':name,'error':float(e),'limit':t,'PASS':bool(e<=t)})
 a=np.exp(.4j*np.arange(8))*(1+np.arange(8)/8);b=np.exp(.33j*np.arange(8))*(1.4-np.arange(8)/20)
 add('common_complex_source_scalar_null',float(shape(a[None],(3+2j)*a[None])[0]));u=a/np.linalg.norm(a);v=b/np.linalg.norm(b);e=shape(a[None],b[None])[0];P=np.outer(u,u.conj())-np.outer(v,v.conj());add('projector_identity',abs(np.linalg.norm(P,2)-e));g=geometries();add('physical_mirror_yzero',max(np.max(abs(v[:,1])) for v in g.values()),0)
 r=ranges([[55,2.5]],g['A_P4M']);mirror=ranges([[55,-2.5]],g['A_P4M']);add('exact_horizontal_mirror',float(np.max(abs(r-mirror))),0);add('same_channel_duplicates',float(np.max(abs(g['B2_DUP'][8:]-g['A0_0'][3]))),0)
 x=np.array([1.,1.2,1.4]);y=np.array([2+1j,3+2j,4+3j]);q=1.1;val=np.interp(q,x,y.real)+1j*np.interp(q,x,y.imag);add('actualaxis_linear_complex_interp',abs(val-(2.5+1.5j)))
 p=read(OUT/'SOLVER_PLAN.json');add('locked_plan_counts_and_N801_1601',0 if len(p['primary'])==18 and len(p['FIELDS'])==144 and all(r['request_count'] in [801,1601] for r in p['FIELDS']) else 1,0);add('frozen_physics_whitelist',0 if all(r['physical_SSP_and_boundary_bytes_unchanged']=='True' for r in rows(OUT/'ENVIRONMENT_WHITELIST_AUDIT.csv')) else 1,0)
 assert len(checks)==8;table(OUT/'PREFLIGHT_CONTROLS.csv',checks);dump(OUT/'PREFLIGHT.json',{'PASS':all(r['PASS'] for r in checks),'algebra_control_cases':8,'new_acoustic_calls':0,'new_samples':0,'FFP':'READY_TO_ATTEMPT_RUNTIME_UNCERTIFIED' if read(OUT/'EXECUTION_CONTRACT.json')['static_FFP_ready'] else 'INDEPENDENT_SOLVER_NOT_READY','budget':storage_guard(),'utc':utc()});print('PREFLIGHT',all(r['PASS'] for r in checks),flush=True)

def freeze_design():
 assert read(OUT/'PREFLIGHT.json')['PASS'];verify_sources();codes=list(ROOT.glob('t2_n2_*.py'))
 for p in codes:ast.parse(p.read_text(encoding='utf-8-sig'))
 inputs=[p for p in OUT.rglob('*') if p.is_file() and p.name!='DESIGN_FREEZE.json'];dump(OUT/'DESIGN_FREEZE.json',{'parent_SHA':PARENT,'task_SHA256':TASKSHA,'utc':utc(),'code_and_input_sha256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in codes+inputs},'runtime':{'python':sys.version,'numpy':np.__version__,'BLAS_threads':1,'PYTHONPATH':'D:/ProjectStorage/WSL-CZ/HLA_H1/runtime_packages'},'solver_budget':{'KRAKEN':9,'SCOOTER':9,'FIELDS':144},'new_acoustics_before_A':0,'historic_G0':'FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED','A_B_performance':'NOT_EVALUATED','CHIGH_global':'UNCERTIFIED','R4_percent':0,'STOP_after_B':True});print('DESIGN FROZEN',len(codes),len(inputs),flush=True)
if __name__=='__main__':{'prepare':prepare,'preflight':preflight,'freeze':freeze_design}[sys.argv[1]]()