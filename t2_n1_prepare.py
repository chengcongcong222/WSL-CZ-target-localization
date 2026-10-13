"""A-stage identity/header checks and deterministic formula controls, no archived scores."""
from t2_n1_common import *
import subprocess,shutil,sys,ast,platform

def prepare():
 assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==PARENT
 assert subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],text=True).split()[0]==PARENT
 if OUT.exists() or LOCAL.exists():raise RuntimeError('new-stage directories already exist; no overwrite')
 OUT.mkdir();LOCAL.mkdir();m=read(PKG/'T2_N1_FILE_MANIFEST_20261013_v1.json')
 for r in m['files']:
  p=PKG/r['filename'];assert sha(p)==r['sha256'] and p.stat().st_size==r['size_bytes'];shutil.copyfile(p,OUT/p.name)
 shutil.copyfile(PKG/'T2_N1_FILE_MANIFEST_20261013_v1.json',OUT/'T2_N1_FILE_MANIFEST_20261013_v1.json');(OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8')
 names=['GPT_SYNC.md','DECISION.json','FORWARD_DECISION.json','FORWARD_CATALOGUE_NUMERICAL_STABILITY.csv','FORWARD_WEAK_FIELD_AND_FAILURE_COUNTS.csv','OUTPUT_RECOVERY_INPUT_SEAL.json','LOCAL_ARCHIVE_INDEX.csv','REPLICA_FREEZE.json','PROVIDER_MODE_WINDOW_AND_MESH.csv','FIELD_STORAGE_QUANTIZATION_COLD_DIAGNOSTIC.csv','CATALOGUE_primary.npy','CATALOGUE_fine.npy','ARRAY_GEOMETRY_TABLE.csv','T2_TILT_G0_EXPERIMENT_SPEC_20261013_v1.json'];sources=[]
 for name in names:
  p=OLD/name;rel=str(p.relative_to(ROOT)).replace('\\','/');blob=subprocess.check_output(['git','show',PARENT+':'+rel]);sources.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'parent_blob_sha256':hashlib.sha256(blob).hexdigest(),'usage':'read-only historical source, not a new result'})
 p=ROOT/'t2_tilt_g0_core.py';sources.append({'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'parent_blob_sha256':hashlib.sha256(subprocess.check_output(['git','show',PARENT+':t2_tilt_g0_core.py'])).hexdigest(),'usage':'read-only formula provenance; not imported or executed'})
 text=p.read_text(encoding='utf-8-sig');assert 'def shape_error(' in text and 'def score(' in text and 'def scores_many(' in text
 freeze=read(OLD/'REPLICA_FREEZE.json');index={r['path']:r for r in rows(OLD/'LOCAL_ARCHIVE_INDEX.csv')};seal=read(OLD/'OUTPUT_RECOVERY_INPUT_SEAL.json')['files'];arrays=[]
 for cat,count in [('primary',3355),('fine',12705)]:
  c=catalog(cat);assert c.shape==(count,3) and c.dtype==np.float64
  expected=np.stack(np.meshgrid(np.linspace(45,60,61 if cat=='primary' else 121),np.linspace(0,5,11 if cat=='primary' else 21),np.arange(180,221,10),indexing='ij'),axis=-1).reshape(-1,3);assert np.array_equal(c,expected)
  for probe,r,t,z in PROBES:assert len(np.flatnonzero(np.all(c==[r,t,z],axis=1)))==1
  for res in RES:
   for mesh in [80001,160001]:
    name=f'{res}_{cat}_n{mesh}.npy';paths=[k for k in index if Path(k).name==name];assert len(paths)==1;p=Path(paths[0]);h=sha(p);assert h==index[str(p)]['SHA256']==freeze['files'][name]==seal[str(p)] and p.stat().st_size==int(index[str(p)]['bytes']);a=np.load(p,mmap_mode='r',allow_pickle=False);assert a.shape==(count,3,12 if res.startswith('B') else 8) and a.dtype==np.complex128
    arrays.append({'filename':name,'path':str(p),'sha256':h,'bytes':p.stat().st_size,'shape':list(a.shape),'dtype':str(a.dtype),'resource':res,'catalogue':cat,'mesh':mesh,'role':'archived complex field'})
   name=f'SHAPE_ERROR_{res}_{cat}.npy';paths=[k for k in index if Path(k).name==name];assert len(paths)==1;p=Path(paths[0]);h=sha(p);assert h==index[str(p)]['SHA256']==freeze['complete_raw_shape_differences'][name]==seal[str(p)];a=np.load(p,mmap_mode='r',allow_pickle=False);assert a.shape==(count,3) and a.dtype==np.float64
   arrays.append({'filename':name,'path':str(p),'sha256':h,'bytes':p.stat().st_size,'shape':list(a.shape),'dtype':str(a.dtype),'resource':res,'catalogue':cat,'mesh':'comparison','role':'historical per-candidate shape errors'})
 dump(OUT/'INPUT_IDENTITY_AND_SCOPE.json',{'parent_SHA':PARENT,'sources':sources,'arrays':arrays,'36_complex_and18_error_arrays':True,'pre_A_only_hashes_headers_coordinate_identity':True,'new_scores_or_ranks_before_A':False,'historical_G0':'FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED','A':'NOT_EVALUATED','B':'NOT_EVALUATED','CHIGH':'UNCERTIFIED','R4_percent':0})
 dump(OUT/'TASK_SOURCE_LOCK.json',{'task_SHA256':TASKSHA,'files':m['files'],'parent_SHA':PARENT,'utc':utc()})
 print('IDENTITY/HDR qualified',len(sources),len(arrays),'no data metrics computed',flush=True)

def preflight():
 # 10 explicitly constructed small deterministic vectors; no random source/noise.
 a=np.exp(.37j*np.arange(8))*(1+np.arange(8)/10);b=np.cos(.73*np.arange(8))+1j*np.sin(.43*np.arange(8));orth=b-a*np.vdot(a,b)/np.vdot(a,a)
 v={'a':a,'b':b,'scaled_phase_a':a*(.27+1.3j),'near_a':a+1e-11*b,'tiny_a':a*1e-120,'orth_b':orth,'dup_a':np.r_[a,np.repeat(a[3],4)],'dup_b':np.r_[b,np.repeat(b[3],4)],'zero':np.zeros(8,complex),'nan':np.r_[np.nan+0j,a[1:]]};cases=[('general','a','b'),('common_scale_phase','a','scaled_phase_a'),('tiny_angle','a','near_a'),('tiny_nonzero_norm','tiny_a','b'),('orthogonal','a','orth_b'),('independent_duplicate_channels','dup_a','dup_b')];checks=[]
 def add(n,e,t=TOL):checks.append({'control':n,'error':float(e),'limit':t,'PASS':bool(e<=t)})
 for name,ai,bi in cases:
  aa=v[ai]/np.linalg.norm(v[ai]);bb=v[bi]/np.linalg.norm(v[bi]);e=np.linalg.norm(aa-bb*np.vdot(bb,aa));difference=np.outer(aa,aa.conj())-np.outer(bb,bb.conj());opnorm=np.max(abs(np.linalg.eigvalsh(difference)));add(name+'_projector_identity',abs(opnorm-e));add(name+'_trace_zero',abs(np.trace(difference)))
  for wi in [aa,bb]:add(name+'_fixed_observation_bound',max(0.,abs(abs(np.vdot(aa,wi))**2-abs(np.vdot(bb,wi))**2)-e))
  qr=np.linalg.qr(np.column_stack([bb,aa]),mode='reduced')[1];add(name+'_stable_QR_residual',abs(abs(qr[1,1])-e))
 for name in ['zero','nan']:
  rejected=False
  try:normalize(v[name])
  except ValueError:rejected=True
  add(name+'_explicit_input_rejection',0 if rejected else 1,0)
 # K repeats preserve the deterministic worst bound; no noise averaging factor.
 add('K8_fixed_template_repeat',abs(np.mean(np.repeat(.6,8))-.6));add('stable_candidateID_ties',0 if np.array_equal(ranking(np.array([1.,1.,.7])),[0,1,2]) else 1,0)
 table(OUT/'ANALYTIC_PROJECTOR_IDENTITY_CHECKS.csv',checks);dump(OUT/'PREFLIGHT.json',{'PASS':all(r['PASS'] for r in checks),'small_constructed_vectors':10,'random_draws':0,'controls':len(checks),'real_archived_score_computations':0,'budget':budget(),'utc':utc()});print('PREFLIGHT',all(r['PASS'] for r in checks),len(checks),flush=True)

def freeze_design():
 assert read(OUT/'PREFLIGHT.json')['PASS'];verify_identities()
 spec=read(OUT/'T2_N1_EXPERIMENT_SPEC_20261013_v1.json')
 config={'task_parameters':spec,'normalization':'complex128 Euclidean norm; zero/nonfinite vector => input anomaly STOP preserving IDs; never silently zero','e_main':'norm(u_H-u_L*vdot(u_L,u_H)), stable equivalent of sqrt(1-|inner|^2)','e_cold':'Householder thin QR of [u_L,u_H], abs(R22), full vectors; projector eigvalsh separate deterministiccontrols','main_scores':'mean3 abs(einsum(u_candidate.conj(),fixed_H_probe))^2, K identical algebraiccopies; no synthetic spectra','cold_scores':'full Hermitian projector quadratic forms using independently normalized raw archives, chunked; no main scoring imports','ranking':'strict raw float64 descending score then lower candidateID; no epsilon retie, no truth','roundoff_diagnostic':'also flag abs(G/Adv)<=1e-10; do not overwrite raw signs/ranks; independent near-tied argmax ID changes separately reported, score tolerance1e-10 does not certify IDs','statistics':'median conventional, P90/P95/P99 nearest rank ceil(q*N); Spearman average ranks; equal normdeciles floor(10*ascending-stable-position/N)','strata':'marginal range/angle/depth bins plus joint3x3x5, no deletion; complete nonnegative reps, mirrors remain implicit','probe_ids':PROBES,'wrong_distance_radius_km':1.,'historical_error_gate':.002,'historical_gate_interpretation':'UNCHANGED_OLD_G0_ONLY; no new admission threshold','inequality_roundoff_tolerance':TOL,'archive_e_discrepancy_tolerance':TOL,'cold_score_tolerance':TOL,'robust_pair':'all10 (bothcatalogues*5) strict abs(AdvH)>Bsum and consistent direction; report A/B separately','sensitive':'any most-dangerous wrong candidateID change, raw G sign flip or main-pair Adv sign flip; report raw nearzero flags too','unresolved':'pair not all10 robust; may coexist with sensitive','budgets':spec['budgets'],'BLAS_threads':1,'single_science_process':True,'new_propagation_MC_receiver_spectrum_optimization_FIM':0,'R4_percent':0,'historical_G0':'FORWARD_NUMERIC_NOT_ADMITTED_UNCHANGED','A_B':'NOT_EVALUATED','CHIGH':'UNCERTIFIED','STOP_after_B':True}
 dump(OUT/'N1_CONFIG_FREEZE.json',config)
 codes=list(ROOT.glob('t2_n1_*.py'))
 for p in codes:ast.parse(p.read_text(encoding='utf-8-sig'))
 inputs=[p for p in OUT.rglob('*') if p.is_file() and p.name!='N1_DESIGN_FREEZE.json']
 dump(OUT/'N1_DESIGN_FREEZE.json',{'parent_SHA':PARENT,'task_SHA256':TASKSHA,'utc':utc(),'code_and_input_SHA256':{str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in codes+inputs},'runtime':{'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),'PYTHONPATH':'D:/ProjectStorage/WSL-CZ/HLA_H1/runtime_packages','BLAS_threads':1},'no_archive_science_before_A':True,'state':'DESIGN_FROZEN_FOR_ONE_EXECUTION','R4_percent':0})
 print('DESIGN FROZEN',len(codes),'codes',len(inputs),'inputs',flush=True)
if __name__=='__main__':{'prepare':prepare,'preflight':preflight,'freeze':freeze_design}[sys.argv[1]]()