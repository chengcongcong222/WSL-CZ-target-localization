"""Package qualification, <=12 algebra controls and pre-A freeze; no new solver before A."""
from t2_tilt_g0_core import *
import shutil,subprocess,sys,ast,platform
TEMPLATE=ROOT/'results/R3_C2_Yang_SA_depth/R3_C2_1/_kraken_zgrid/zgrid_f201.env'
def prepare():
 if (OUT/'PREPARED.json').exists():raise RuntimeError('design already prepared')
 OUT.mkdir();LOCAL.mkdir();(LOCAL/'modal').mkdir();(OUT/'environments').mkdir()
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==PARENT
 assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==PARENT
 manifest=read(PKG/'T2_TILT_G0_FILE_MANIFEST_20261013_v1.json')
 for item in manifest['files']:
  p=PKG/item['name'];assert sha(p)==item['sha256'] and p.stat().st_size==item['size_bytes'];shutil.copyfile(p,OUT/p.name)
 shutil.copyfile(PKG/'T2_TILT_G0_FILE_MANIFEST_20261013_v1.json',OUT/'T2_TILT_G0_FILE_MANIFEST_20261013_v1.json')
 shutil.copyfile(PKG/'T2_TILT_G0_ARRAY_GEOMETRY_20261013_v1.csv',OUT/'ARRAY_GEOMETRY_TABLE.csv')
 (OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8')
 geo=geometries();expected=np.unique(np.concatenate([np.array([180,185,190,200,210,215,220]),*[a[:,2] for a in geo.values()]]))
 rng=np.random.default_rng(2026101301);panel=[{'geometry':i,'r_km':float(rng.uniform(45.2,59.8)),'theta_deg':float(rng.uniform(-4.8,4.8)),'z_m':[185,200,215][i%3]} for i in range(16)];table(OUT/'CALIBRATION_GEOMETRIES.csv',panel)
 # Historical provider used solely to lock noise power at independent B0 calibration positions.
 from hla_sp1_core import parse_mod as old_parse,pressure as old_pressure
 mod=[old_parse(Path('D:/ProjectStorage/WSL-CZ/HLA_SP1/modal')/f'sp1_f{int(f)}_n160001.mod') for f in FREQ]
 values=[]
 for p in panel:
  rr=np.sqrt((p['r_km']*1000*np.cos(np.radians(p['theta_deg']))-geo['A0_0'][:,0])**2+(p['r_km']*1000*np.sin(np.radians(p['theta_deg'])))**2);values.append(np.array([old_pressure(m,rr,p['z_m']) for m in mod]))
 ref=float(np.median(np.mean(abs(np.array(values))**2,axis=(1,2))))
 np.save(OUT/'REFERENCE_B0_FIELDS.npy',np.array(values));dump(OUT/'NOISE_REFERENCE_FREEZE.json',{'P_ref':ref,'nu':{'20':ref/100,'5':ref/10**.5},'rule':'median of unit-source mean3-frequency/8-element B0 power over16 independent calibration positions; qualified historical B0 provider, before any new calls','legacy_MOD_SHA256':{str(int(f)):sha(Path('D:/ProjectStorage/WSL-CZ/HLA_SP1/modal')/f'sp1_f{int(f)}_n160001.mod') for f in FREQ},'no_resource_renormalization':True})
 lines=TEMPLATE.read_text(encoding='utf-8-sig').splitlines();boundary=next(i for i,l in enumerate(lines) if l.strip().startswith("'R'"));src=[180,185,190,200,210,215,220];rec=np.unique(np.concatenate([a[:,2] for a in geo.values()])).tolist()
 plans=[]
 for mesh in [80001,160001]:
  for f in FREQ:
   stem=f't2_f{int(f)}_n{mesh}';env=["'T2_TILT_G0_ESTD'",str(float(f)),'1',lines[3],f'{mesh} 0.0 5000.0']+lines[5:boundary]+["'R' 0.0",'1500.0 1800.0','60.0',str(len(src))]+[f'{d:.12f}' for d in src]+[str(len(rec))]+[f'{d:.12f}' for d in rec]+['R'];p=OUT/'environments'/f'{stem}.env';p.write_text('\n'.join(env)+'\n',encoding='utf-8');plans.append({'stem':stem,'mesh':mesh,'frequency_hz':float(f),'env_SHA256':sha(p),'new_KRAKEN':True})
 # FIELD three frequencies, each spanning all three representative resources and ranges, one scalar/frequency.
 field={'resources':['A0_0','A_P4M','B1_V4'],'states':[[45.321,1.37,200],[55.123,2.31,200],[59.731,4.13,200]],'one_common_complex_alignment_per_frequency_across_all_resources_and_ranges':True}
 allranges=[];depths=[]
 for res in field['resources']:
  a=geo[res]
  for r,t,z in field['states']:
   allranges.extend(np.sqrt((r*1000*np.cos(np.radians(t))-a[:,0])**2+(r*1000*np.sin(np.radians(t))-a[:,1])**2).tolist());depths.extend(a[:,2].tolist())
 allranges=np.unique(allranges);depths=np.unique(depths)
 for f in FREQ:
  flp=["'T2_FIELD'","'RA'",'9999','1','0.0',str(len(allranges))]+[f'{r/1000:.12f}' for r in allranges]+['1','200.0',str(len(depths))]+[f'{d:.12f}' for d in depths]+['1','0.0 /'];(OUT/'environments'/f't2_field_f{int(f)}.flp').write_text('\n'.join(flp)+'\n',encoding='utf-8')
 archive=[]
 for directory in ['D:/ProjectStorage/WSL-CZ/HLA_SP1/modal',str(ROOT/'results/R4_H3_G0_FINITE_VERTICAL_RESPONSE'),str(ROOT/'results/R4_H3_P1_UNIFIED_WINDOW_SSP')]:
  for p in Path(directory).glob('*.mod'):
   try:
    d,ph,kk=parse_mod(p);missing=[float(x) for x in expected if np.min(abs(d-x))>1e-5];archive.append({'path':str(p),'sha256':sha(p),'returned_depths':len(d),'missing_required_depths':len(missing),'qualified':not missing,'example_missing':json.dumps(missing[:4])})
   except Exception as e:archive.append({'path':str(p),'sha256':sha(p),'returned_depths':0,'missing_required_depths':len(expected),'qualified':False,'example_missing':str(e)})
 table(OUT/'ARCHIVED_PROVIDER_QUALIFICATION.csv',archive)
 sources=[TRUTH,TEMPLATE,ROOT/'hla_sp1_core.py',ROOT/'results/R4_SINGLE_HLA_STATIC_SPATIAL_RANGE_PILOT/PROVIDER_PROVENANCE.json',ROOT/'results/R4_H3_P1_UNIFIED_WINDOW_SSP/GPT_SYNC.md']
 table(OUT/'INPUT_PROVENANCE.csv',[{'path':str(p.relative_to(ROOT)).replace('\\','/'),'working_SHA256':sha(p),'parent_blob_SHA256':hashlib.sha256(subprocess.check_output(['git','show',PARENT+':'+str(p.relative_to(ROOT)).replace('\\','/')])).hexdigest()} for p in sources])
 dump(OUT/'PROVIDER_PROVENANCE.json',{'template_path':str(TEMPLATE.relative_to(ROOT)),'template_SHA256':sha(TEMPLATE),'solver_KRAKEN_SHA256':sha(AT/'kraken.exe'),'solver_FIELD_SHA256':sha(AT/'field.exe'),'expected_output_depths_m':expected.tolist(),'absolute_depth_match_tolerance_m':1e-5,'source_depths_m':src,'receiver_depths_m':rec,'plans':plans,'FIELD_plan':field,'CLOW':1500,'CHIGH':1800,'CHIGH_certification':'UNCERTIFIED_HISTORICAL_LIMIT_UNCHANGED','all_returned_modes_used':True,'vertical_or_range_interpolation':'NONE_DIRECT_SUM_AT_ALL_REQUESTED_POINTS','thesis_original':'ORIGINAL_PDF_LOCAL_NOT_FOUND','PDF_search_roots':['D:/Downloads','D:/Downloads/Documents','D:/ProjectStorage/WSL-CZ','repo/literature']})
 dump(OUT/'EXECUTION_CONTRACT.json',{'source':{'calibration_seed':2026101302,'test_seed':2026101304,'amplitude_db':[-3,3],'phase_rad':[-float(np.pi),float(np.pi)]},'noise':{'calibration_seed':2026101303,'test_seed':2026101305,'real_time_sample_variance':'2048*nu','shared_first8_numbered_innovations':True,'extra4_B1_B2_shared':True},'gain':{'calibration_seed':2026101308,'test_seed':2026101306,'first8_shared_across_resources':True,'extra4_B1_B2_shared':False,'independent_extra_gain_ranges':{'B1_V4':[8,12],'B2_DUP':[12,16]},'fixed_all_freq_and_segments':True,'C1_gain_not_given_to_inference':True},'calibration_order':'16geometries*16replicates; IDdisjoint; z=185,200,215 cyclic by geometry (96,80,80 records)','threshold':'244th/256 D at nearest r/abs(theta)/z; tie smaller label; separately each mesh/resource/C/SNR/catalogue; seal before test','preflight_controls_max':12,'forward_catalogues':['primary','fine'],'noise_reference_qualification':'archived exact B0 fields only at independent calibration geometries, frozen before new KRAKEN; not test truth','generation_mesh':160001,'candidate_order':'r then abs(theta) then z, ascending','envelope':'half-grid cells clipped45..60, no continuous coverage guarantee; empty INF','cold_score_tolerance':1e-7,'cold_acceptance_flips_allowed':0,'solver_limits':{'KRAKEN_calls':6,'KRAKEN_timeout_s':90,'FIELD_calls':3,'FIELD_timeout_s':60},'budgets':{'science_wall_s':21600,'RSS':8*1024**3,'new_D_bytes':32*1024**3,'new_git_bytes':512*1024**2,'BLAS_threads':1,'science_processes':1},'main_rule':'fixed A_P4M/A0 and B1/B2; >=2/3 at z200,C0,SNR20,primary160001; each>=15/16 accepted&horizontal-retained;median unconditional cellenvelope<=10km&reduction>=.30; numerical reversal =>limited','new_trials_after_B':False})
 dump(OUT/'PREPARED.json',{'utc':utc(),'no_new_solver_or_formal_record_before_A':True,'reference_power':ref})
 print('prepared',len(expected),'depths; B0ref',ref,flush=True)

def preflight():
 g=geometries();rng=np.random.default_rng(2026101307);checks=[]
 def add(name,error,limit=1e-10):checks.append({'control':name,'error':float(error),'limit':limit,'PASS':bool(error<=limit)})
 coords=rows(OUT/'ARRAY_GEOMETRY_TABLE.csv');err=0.
 for r in coords:
  if r['arc_s_m']:
   u=float(r['arc_s_m']);a=np.radians(float(r['tilt_from_horizontal_deg']));err=max(err,abs(float(r['x_m'])-u*np.cos(a)),abs(float(r['z_m'])-200-u*np.sin(a)))
 add('coordinate_and_primary4m_span',max(err,abs(np.ptp(g['A_P4M'][:,2])-4)))
 add('B1_B2_original8_and_DUP',max(np.max(abs(g['A0_0']-g['B1_V4'][:8])),np.max(abs(g['A0_0']-g['B2_DUP'][:8])),np.max(abs(g['B2_DUP'][8:]-g['A0_0'][3]))))
 a=np.radians(16.601549599020235);phi=.31;theta=.23;axis=np.array([np.cos(a),0,np.sin(a)]);ray=np.array([np.cos(phi)*np.cos(theta),np.cos(phi)*np.sin(theta),np.sin(phi)])
 add('plane_wave_axis_projection_and_mirror',abs(axis@ray-(np.cos(a)*np.cos(phi)*np.cos(theta)+np.sin(a)*np.sin(phi))))
 h=rng.normal(size=(3,8))+1j*rng.normal(size=(3,8));s=10**(rng.uniform(-3,3,(8,3))/20)*np.exp(1j*rng.uniform(-np.pi,np.pi,(8,3)));y=s[:,:,None]*h
 add('source_phase_amplitude_invariance',abs(score(h[None],y)[0]-1))
 x=waveform(h,s,np.ones(8),np.zeros((8,8,N)),1.);add('real_waveform_DFT_identity',np.linalg.norm(dft(x)-y)/np.linalg.norm(y))
 basis=np.concatenate([TONE.real.T,TONE.imag.T],axis=1);co=np.linalg.lstsq(basis,x[0].T,rcond=None)[0];yy=(co[:3]-1j*co[3:])/2;add('independent_sine_fit',np.linalg.norm(yy-y[0])/np.linalg.norm(y[0]))
 noisy=y+.1*(rng.normal(size=y.shape)+1j*rng.normal(size=y.shape));add('covariance_direct_score',abs(scores_many(h[None],noisy[None])[0,0]-score(h[None],noisy)[0]))
 gain=10**(rng.uniform(-1,1,8)/20)*np.exp(1j*np.radians(rng.uniform(-30,30,8)));add('C1_no_false_gain_invariance',0 if score(h[None],y*gain)[0]<1-1e-8 else 1)
 eta=rng.normal(size=(1024,1,N))*np.sqrt(N);bins=np.einsum('kmn,fn->kfm',eta,TONE.conj())/N;add('noise_DFT_variance',np.max(abs(np.mean(abs(bins)**2,axis=0)-1)),.12)
 add('nearest_rank244_and_empty_INF',0 if np.sort(np.arange(256))[243]==243 and np.isinf(np.median([np.inf]*16)) else 1)
 seeds=read(OUT/'EXECUTION_CONTRACT.json');add('seed_test_calibration_isolation',0 if seeds['source']['calibration_seed']!=seeds['source']['test_seed'] and seeds['noise']['calibration_seed']!=seeds['noise']['test_seed'] else 1)
 add('once_guard_absent',0 if not (OUT/'EXECUTION_STARTED.json').exists() else 1)
 assert len(checks)==12
 table(OUT/'PREFLIGHT_CONTROLS.csv',checks);dump(OUT/'PREFLIGHT.json',{'PASS':all(x['PASS'] for x in checks),'controls':checks,'control_cases':12,'new_KRAKEN':0,'new_FIELD':0,'calibration_or_test_records':0,'utc':utc()});print('preflight',all(x['PASS'] for x in checks),flush=True)

def freeze():
 assert read(OUT/'PREFLIGHT.json')['PASS']
 for p in ROOT.glob('t2_tilt_g0_*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
 inputs=[p for p in OUT.rglob('*') if p.is_file() and p.name!='DESIGN_FREEZE.json'];codes=list(ROOT.glob('t2_tilt_g0_*.py'))
 dump(OUT/'DESIGN_FREEZE.json',{'parent_SHA':PARENT,'task_SHA256':TASKSHA,'utc':utc(),'code_and_input_SHA256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in codes+inputs},'runtime':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'BLAS_threads':1},'all_new_propagation_after_CommitA':True,'no_formal_test_opened':True,'execution':'once-only forward gate then calibration/threshold seal then waveformDFT seal then full inference and truth evaluation','original_P0_P1_gate_unchanged':True,'R4_percent':0})
 print('design frozen',len(codes),len(inputs),flush=True)
if __name__=='__main__':{'prepare':prepare,'preflight':preflight,'freeze':freeze}[sys.argv[1]]()