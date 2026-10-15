"""Full preflight. Checks identity and algebra only; does not solve registered modes."""
from pathlib import Path
import ast,csv,json,hashlib,subprocess,sys,shutil,math
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2');PARENT='35ce80055c641711ace8639f528679917bfdbef2'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def git(*a):return subprocess.check_output(['git','-c','safe.directory='+ROOT.as_posix(),*a],cwd=ROOT).decode().strip()
checks=[]
def ck(name,condition):
 checks.append({'check':name,'PASS':bool(condition)})
 if not condition:
  save(LOCAL/'FUNC2_PREFLIGHT_STOP.json',{'failed':name,'checks':checks,'science_executions':0,'STOP':True});raise RuntimeError('PRECHECK_STOP '+name)
def main():
 ck('parent_HEAD',git('rev-parse','HEAD')==PARENT);ck('remote_parent',git('-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main').split()[0]==PARENT);ck('no_tracked_or_staged_changes',not git('diff','--name-only') and not git('diff','--cached','--name-only'))
 initial=read(LOCAL/'INITIAL_STATE.json');ck('ledger_unchanged',sha(ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv')==initial['ledger_sha256']);ck('no_formal_execution',not (OUT/'EXECUTION_STARTED.json').exists());ck('task_bytes_preserved',sha(OUT/'RESEARCH_LEAD_TASK.txt')==initial['task_sha256'])
 needed=[]
 for sub in ['R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY','R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION','R4_T2_PREC3_FUNC1_MODAL_FUNCTION_ERROR_AUDIT']:
  folder=ROOT/'results'/sub
  for p in sorted(folder.rglob('*')):
   if not p.is_file() or 'ORIGINAL_TASK_PACKAGE' in p.parts or p.suffix.lower() in ['.png','.pdf']:continue
   # Read every text/numeric record, including all rows, and bind original blob bytes.
   b=p.read_bytes();rel=p.relative_to(ROOT).as_posix()
   if p.suffix in ['.py','.json','.csv','.md','.txt','.log']:b.decode('utf-8-sig')
   needed.append(rel)
 for rel in ['results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/PROPAGATION_CONTRACT_AND_SCOPE.md','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/SOURCE_IDENTITY.json','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/SCIENCE_OUTPUT_FREEZE.json','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/SOLVER_PLAN.json','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/ARRAY_GEOMETRY_TABLE.csv','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/T2_N2_SENTINELS_20261013_v1.csv','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/environments/n2_K2_f201.env','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv']:
  needed.append(rel)
 needed=list(dict.fromkeys(needed));query=''.join(PARENT+':'+p+'\n' for p in needed).encode('utf-8');batch=subprocess.run(['git','-c','safe.directory='+ROOT.as_posix(),'cat-file','--batch'],cwd=ROOT,input=query,capture_output=True,check=True).stdout;offset=0;inputs=[]
 for path in needed:
  end=batch.index(b'\n',offset);header=batch[offset:end].split();ck('parent_blob_exists/'+path,len(header)==3 and header[1]==b'blob');n=int(header[2]);b=batch[end+1:end+1+n];offset=end+1+n+1;
  physical=(ROOT/path).read_bytes()
  if path=='results/R4_MASTER/R4_EVIDENCE_LEDGER.csv':
   ck('ledger_Git_text_identity',physical.replace(b'\r\n',b'\n')==b)
   ck('ledger_original_physical_bytes_preserved',hashlib.sha256(physical).hexdigest()==initial['ledger_sha256'])
   save(OUT/'LEDGER_REPRESENTATION_AUDIT.json',{'Git_parent_SHA256':hashlib.sha256(b).hexdigest(),'working_SHA256':hashlib.sha256(physical).hexdigest(),'working_CRLF':physical.count(b'\r\n'),'difference':'CRLF_checkout_vs_LF_git_blob_only','content_changed':False,'science_source_identity_failure':False,'original_preflight_preserved':'FUNC2_FIRST_PREFLIGHT_PAUSE.json'})
   b=physical;n=len(b)
  else:ck('exact_parent_bytes/'+path,physical==b)
  inputs.append({'path':path,'sha256':hashlib.sha256(b).hexdigest(),'bytes':n})
 prior=read(ROOT/'results/R4_T2_PREC3_FUNC1_MODAL_FUNCTION_ERROR_AUDIT/FUNC1_DESIGN_FREEZE.json');vectors=prior['vectors']
 for r in vectors:
  ck('vector_SHA/'+r['path'],sha(r['path'])==r['sha256'])
  with np.load(r['path'],allow_pickle=False) as a:
   p=a['phi'];h=float(a['h']);ck('vector_contract/'+r['path'],set(a.files)=={'phi','h'} and p.dtype==np.dtype('complex128') and p.shape==(r['N'],) and h==5000/r['N'] and np.isfinite(p).all());v=np.r_[0j,p];norm=h/3*np.sum(v[:-1]**2+v[:-1]*v[1:]+v[1:]**2) if r['method']=='FE' else h*(np.sum(p[:-1]**2)+p[-1]**2/2);ck('bilinear_norm/'+r['path'],abs(norm-1)<=1e-8)
 ck('exact_vector_set_no_FD961',len(vectors)==10 and not any(v['method']=='FD' and v['mode']==961 for v in vectors))
 env=ROOT/prior['env_path'];lines=env.read_text(encoding='utf-8-sig').splitlines();end=next(i for i,l in enumerate(lines) if l.startswith("'R'"));ssp=np.array([list(map(float,line.split()[:6])) for line in lines[5:end]])
 ck('physical_contract',ssp.shape==(251,6) and np.array_equal(ssp[:,0],np.arange(251)*20) and lines[3].strip()=="'CVWT'" and float(lines[1])==201 and np.all(ssp[:,3]==1) and np.all(ssp[:,2]==0) and np.all(ssp[:,4:]==0))
 import func2_dop853 as primary
 import func2_gauss as secondary
 for level in [0,1]:
  obj=secondary.Collocation(np.array([0,5000]),np.array([1500.,1500.]),1,level)
  ck('Gauss_row_sums/'+str(level),np.max(abs(obj.A.sum(axis=1)-obj.ci))<1e-12);ck('Gauss_weights/'+str(level),abs(sum(obj.w)-1)<1e-14)
  T,der=obj.transfer(math.tau*201/1500,0,5);exact=np.array([[1,5*primary.S],[0,1]],complex);ck('constant_zero_potential_transfer/'+str(level),np.linalg.norm(T-exact)<1e-12)
  coeff=np.arange(1,obj.order+1);integrals=np.polynomial.legendre.legval(1,obj.anti)-obj.offset;ck('Gauss_dense_endpoint/'+str(level),np.max(abs(integrals-obj.w))<1e-12)
 shooter=primary.Shoot([0,20],[1500.,1500.],10,2e-10,2e-12,2.);endv,log,_=shooter.side(primary.S,0,20,[0,1]);target=np.array([20*primary.S,1.]);ck('DOP_zero_potential_transfer',np.linalg.norm(endv*np.exp(log)-target)<1e-11)
 for p in OUT.glob('*.py'):
  code=ast.parse(p.read_text(encoding='utf-8-sig'));ck('AST/'+p.name,True)
  if p.name in ['func2_dop853.py','func2_gauss.py','func2_main.py','func2_cold.py']:
   ck('no_old_eigen_or_external_calls/'+p.name,not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ['eig','eigs','eigsh','eigh','eigh_tridiagonal','Popen','run','fft','ifft'] for n in ast.walk(code)))
 ck('D_disk_budget',shutil.disk_usage('D:/').free>3*1024**3)
 historic=read(ROOT/'results/R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION/PREC3_COLD_PRUEFER.json');initial_roots={str(r['mode']):r['tight'] for r in historic}
 cfg={'stage':'T2-PREC3-FUNC2','parent_SHA':PARENT,'task_sha256':initial['task_sha256'],'inputs':inputs,'vectors':vectors,'env_path':prior['env_path'],'geometry_path':prior['geometry_path'],'probe_path':prior['probe_path'],'resources':prior['resources'],'modes':[1,481,961],'frequency_hz':201,'loss_scales':[0,.25,.5,1.],'initial_roots':initial_roots,'initial_root_source':'archived independent PREC3_COLD_PRUEFER tight; fixed mode identity; shared initial ranges only; each new method solves own roots','root_real_neighborhood_radius':3e-7,'root_imag_domain':[-1e-5,1e-12],'newton_derivative_steps':{'DOP853':2e-7,'GAUSS':1e-7},'max_newton_step':2e-6,'max_newton_iterations':10,'root_change_tolerance':1e-13,'root_mismatch_tolerance':2e-8,'primary':{'method':'segmented DOP853 two-sided shooting','segment_boundaries_m':'original251 knots20m','match_m':{'1':1300,'481':2500,'961':2500},'loose':{'rtol':2e-10,'atol':2e-12,'max_step_m':2},'tight':{'rtol':2e-12,'atol':2e-14,'max_step_m':1},'normalization_Gauss_nodes_per20m':32},'independent':{'method':'Gauss-Legendre implicit collocation transfer BVP with independent root search','subinterval_m':5,'stages':[8,10],'formal_orders':[16,20],'continuous_polynomial_evaluation':True,'normalization_Gauss_nodes_per5m':16},'k_budget':1.6666667e-9,'shape_budget':5e-4,'boundary_budget':1e-8,'boundary_definition':'max(|phi(0)|,|phi_prime(H)|/S)/max|phi|, S=omega/1500','match_budget':2e-8,'match_definition':'norm of scaled [phi,phi_prime/S] mismatch divided by left state norm','norm_budget':1e-6,'normalization':'integral phi(z)^2 dz = 1; not modulus squared; no renormalization of old vectors','comparison_grid_m':.25,'integral_cross_check_m':.5,'historical_old_function_threshold':.002,'depth_bands':'25 contiguous200m bands','weak_phi_amplitude':1e-6,'weak_local_L2_norm':1e-8,'weak_coupling_norm':1e-12,'coupling_budget':5e-4,'coupling_definition':'new independent paths and own accuracy settings: uncalibrated Q relative and projected shape both <=5e-4 for nonweak; weak values preserved, no relative certification','cold':{'method':'fresh source parse, fsum trapezoid, 8th-order sampled second derivative excluding1.1m around SSP knots; independent fixed-k RK45 mode481','rtol':3e-12,'atol':3e-14,'max_step_m':.5,'ODE_residual_budget':1e-6,'matching_budget':1e-7,'coupling_budget':5e-4,'no_final_root_shared_between_primary_methods':True},'budget':{'wall_s':3600,'main_wall_s':2900,'cold_wall_s':700,'tree_RSS_bytes':4*1024**3,'new_D_bytes':2*1024**3,'Git_bytes':64*1024**2,'main_calls':1,'cold_calls':1,'retries':0,'BLAS_threads':1,'external_acoustic_solvers':0,'old_FE_FD_Pruefer_fullmode_MC_waveform_DFT_FIM_localization':0},'exit_rules':{'A_failure':'FUNC2_ANALYTIC_REFERENCE_GATE_FAIL;no E-STD','B_failure':'FUNC2_HIGH_ORDER_CONTINUOUS_REFERENCE_NOT_ESTABLISHED;stop C/D','independent_unavailable':'FUNC2_INDEPENDENT_REFERENCE_NOT_READY','B_pass_old_error_above_0p002':'FUNC2_REFERENCE_SUPPORTED_BUT_OLD_MODAL_FUNCTION_NOT_ADMITTED','local_agreement':'FUNC2_LOCAL_MODAL_COUPLING_SUPPORTED_SUBSET_ONLY','science_failure':'retain all completed data, audit once, stop no retry','final_route':'explicit GO or NO_GO under existing physical data contract; no FUNC3 automatic'},'planned':{'analytic_rows':24,'estd_continuation_roots':48,'reference_modes':3,'old_functions':10,'old_depth_bands':250,'old_coupling_vectors':200,'reference_local_comparisons':180,'missing_FD961':2},'historical_results':'UNCHANGED','R4_percent':0,'frozen_code':[{'path':p.name,'sha256':sha(p),'bytes':p.stat().st_size} for p in OUT.glob('*.py')]}
 save(OUT/'FUNC2_DESIGN_FREEZE.json',cfg);save(OUT/'FUNC2_PREFLIGHT.json',{'PASS':True,'checks':checks,'scientific_executions':0,'sources':len(inputs),'vectors':len(vectors),'Python':sys.version,'NumPy':np.__version__,'unrelated_untracked_state':initial['status'],'writing_style':'No suitable writing references found in the returned results; follow supplied research prose','R4_percent':0})
 with (OUT/'FUNC2_SOURCE_IDENTITIES.csv').open('w',encoding='utf-8',newline='') as f:
  fields=['path','sha256','bytes','kind'];w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows([{**r,'kind':'PARENT_REPOSITORY_BLOB'} for r in inputs]);w.writerows([{'path':v['path'],'sha256':v['sha256'],'bytes':v['bytes'],'kind':'ORIGINAL_ARCHIVED_VECTOR'} for v in vectors])
 print('PREFLIGHT_PASS',len(checks),'sources',len(inputs),'vectors',len(vectors),'science executions0')
if __name__=='__main__':main()
