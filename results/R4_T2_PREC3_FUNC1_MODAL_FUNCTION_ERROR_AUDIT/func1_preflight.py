"""FUNC1 complete admission and freeze preparation, no historical error analysis."""
from pathlib import Path
import json,hashlib,subprocess,sys,csv,ast,shutil,datetime
import numpy as np
ROOT=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC1/worktree');LOCAL=ROOT.parent;IMPL=LOCAL/'implementation';OUT=ROOT/'results/R4_T2_PREC3_FUNC1_MODAL_FUNCTION_ERROR_AUDIT';PARENT='6b4bd489e7fd292adb64e25b72632948f20b2c4c'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
checks=[]
def ck(n,c):
 checks.append({'check':n,'PASS':bool(c)})
 if not c:
  write(LOCAL/'FUNC1_PREFLIGHT_STOP.json',{'PASS':False,'failed':n,'checks':checks,'main_calls':0,'cold_calls':0});raise RuntimeError('PREFLIGHT_STOP '+n)
ck('HEAD_parent',subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip()==PARENT)
ck('remote_parent',subprocess.check_output(['git','-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main'],cwd=ROOT).decode().split()[0]==PARENT)
ck('clean_isolated_tree',not subprocess.check_output(['git','status','--porcelain=v1'],cwd=ROOT).strip())
contexts=read(LOCAL/'CONTEXT_READ_INDEX.json');extras=['results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/environments/n2_K2_f201.env','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/ARRAY_GEOMETRY_TABLE.csv','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/T2_N2_SENTINELS_20261013_v1.csv','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/SCIENCE_OUTPUT_FREEZE.json','results/R4_T2_N2_INDEPENDENT_FORWARD_NUMERIC_BRACKET/SOLVER_PLAN.json']
paths=[r['path'] for r in contexts]+extras;query=''.join(PARENT+':'+p+'\n' for p in paths).encode();result=subprocess.run(['git','cat-file','--batch'],cwd=ROOT,input=query,capture_output=True,check=True).stdout;offset=0;inputs=[]
for p in paths:
 end=result.index(b'\n',offset);n=int(result[offset:end].split()[-1]);b=result[end+1:end+1+n];offset=end+1+n+1;ck('parent_blob/'+p,(ROOT/p).read_bytes()==b)
 original=next((r for r in contexts if r['path']==p),None)
 if original:ck('read_original/'+p,hashlib.sha256(b).hexdigest()==original['sha256'])
 inputs.append({'path':p,'sha256':hashlib.sha256(b).hexdigest(),'bytes':n})
old=ROOT/'results/R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION';index=read(old/'PREC3_LARGE_DATA_INDEX.json');vec=[];identities=[]
for method,modes in [('FE',[1,481,961]),('FD',[1,481])]:
 for m in modes:
  for N in [80000,160000]:
   name=f'{"gal" if method=="FE" else "fd"}_f201_m{m}_N{N}.npz';matches=[r for r in index['files'] if Path(r['path']).name==name];ck('index_unique/'+name,len(matches)==1);r=matches[0];p=Path(r['path']);ck('vector_present/'+name,p.is_file());ck('vector_sha/'+name,sha(p)==r['sha256'])
   with np.load(p,allow_pickle=False) as ar:
    ph=ar['phi'];h=ar['h'];ck('NPZ_keys/'+name,set(ar.files)=={'phi','h'});ck('phi_dtype_length/'+name,ph.dtype==np.dtype('complex128') and ph.shape==(N,));ck('h_dtype_scalar/'+name,h.dtype==np.dtype('float64') and h.shape==());ck('h_exact/'+name,float(h)==5000/N);ck('vector_finite/'+name,np.isfinite(ph).all() and np.isfinite(float(h)))
    v=np.r_[0j,ph];norm=float(h)/3*np.sum(v[:-1]**2+v[:-1]*v[1:]+v[1:]**2) if method=='FE' else float(h)*(sum(ph[:-1]**2)+.5*ph[-1]**2);ck('original_bilinear_norm/'+name,abs(norm-1)<=1e-8)
   vr={'method':method,'mode':m,'N':N,**r};vec.append(vr);identities.append({**vr,'phi_dtype':'complex128','phi_length':N,'h_m':float(h),'bilinear_norm_real':float(norm.real),'bilinear_norm_imag':float(norm.imag),'bilinear_norm_error':float(abs(norm-1)),'normalization_contract':'CONSISTENT_P1_MASS' if method=='FE' else 'FD_TRAPEZOID_PHYSICAL_VECTOR','PASS':True})
ck('exact10vector_set',len(vec)==10 and not any(r['method']=='FD' and r['mode']==961 for r in vec))
for p in IMPL.glob('*.py'):
 tree=ast.parse(p.read_text(encoding='utf-8'));ck('code_AST/'+p.name,True)
 if p.name in ['func1_main.py','func1_cold.py']:
  forbidden={'eig','eigh','eigsh','eigs','eigh_tridiagonal','solve_ivp','brentq','fft','ifft'}
  ck('no_new_solver_calls/'+p.name,not any(isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in forbidden for n in ast.walk(tree)))
sys.path.insert(0,str(IMPL));import func1_main as main;import func1_cold as cold
z=np.array([0.,1.,2.]);g=np.array([1+2j,2-1j,3+1j]);f=2j*g;w=main.weights(z);mm=main.metrics(f,g,w);cc=cold.metric(f,g,z)
ck('synthetic_sign_preserves_complex_phase',abs(mm['e_sign']-np.sqrt(5))<1e-14);ck('synthetic_complex_projection',mm['e_complex']<1e-14 and abs(mm['alpha_imag']-2)<1e-14);ck('synthetic_independent_QR',abs(mm['e_sign']-cc['e_sign'])<1e-14 and cc['e_complex']<1e-14)
ck('synthetic_minus_sign',main.metrics(-g,g,w)['e_sign']==0);ck('zero_reference_not_zero_error',main.metrics(g,np.zeros(3,complex),w)['e_sign'] is None)
grid=np.array([0.,.5,1.]);v=np.array([0j,1+1j,0j]);ck('synthetic_interpolation',np.allclose(cold.interpolate(v,.5,grid),v))
ck('D_reserve',shutil.disk_usage('D:/').free>1024**3)
import psutil
ck('budget_monitor_available',psutil.Process().memory_info().rss<2*1024**3)
ck('no_science_run_marker',not (OUT/'EXECUTION_STARTED.json').exists())
config={'stage':'T2-PREC3-FUNC1','parent_SHA':PARENT,'inputs':inputs,'vectors':vec,'env_path':extras[0],'geometry_path':extras[1],'probe_path':extras[2],'resources':['A0_0','A_P4M','B1_V4','B2_DUP'],'paired_modes':[1,481],'FE_only_modes':[961],'original_N':[80000,160000],'comparison_intervals_m':['SSP20_ORIGINAL_UNWEIGHTED','SSP20_TRAPEZOID',1.,.5,.25],'historical_error_absolute_tolerance':1e-8,'dense_stability_absolute_tolerance':2e-5,'original_PREC3_function_threshold':.002,'evaluation_bias_abs_floor':2e-5,'evaluation_bias_relative':.10,'evaluation_bias_scope':'DESCRIPTIVE; both10percent and2e-5 required; not a replacement Gate or statistical coverage claim','observed_second_order_interval':[1.8,2.2],'order_scope':'paired FE-FD discrepancy only; one own-method pair cannot determine own order','weak_phi_amplitude':1e-6,'weak_local_L2_norm':1e-8,'weak_coupling_norm':1e-12,'weak_scope':'design flags only; retain actual values and ratios; no denominator clipping; exact-zero denominators undefined','cold_abs_tolerance':1e-10,'cold_relative_tolerance':1e-9,'cold_derived_ratio_tolerance':1e-7,'projection_scope':'diagnostic only; never rescale archived functions or raw coupling; no channel-specific calibration; per-vector scalar only to report invariant shape','local_depth_grid':'180..220m inclusive,0.25m;25global200m bands','coupling_count':80,'coupling_channels':800,'mode961_policy':'two FE only; old cold values HISTORICAL_COLD_RESULT_ONLY; no FD reconstruction or eigen solve','budget':{'total_wall_s':1800,'main_wall_s':1200,'cold_wall_s':600,'RSS_bytes':2*1024**3,'new_D_bytes':2*1024**3,'Git_bytes':64*1024**2,'scientific_processes':1,'BLAS_threads':1,'main_calls':1,'cold_calls':1,'retries':0,'all_new_solvers_propagation_MC_waveform_DFT_FIM_localization':0},'historical_results':'UNCHANGED','R4_percent':0,'frozen_code':[{'path':p.name,'sha256':sha(p),'bytes':p.stat().st_size} for p in IMPL.glob('*.py')]}
OUT.mkdir(exist_ok=True);(OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8')
for p in IMPL.glob('*.py'):shutil.copyfile(p,OUT/p.name)
shutil.copyfile(LOCAL/'RESEARCH_LEAD_TASK.txt',OUT/'RESEARCH_LEAD_TASK.txt');shutil.copyfile(LOCAL/'WORKTREE_MATERIALIZATION.json',OUT/'WORKTREE_MATERIALIZATION.json')
write(OUT/'FUNC1_DESIGN_FREEZE.json',config)
with (OUT/'FUNC1_VECTOR_IDENTITIES.csv').open('w',encoding='utf-8',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in identities for k in r)),lineterminator='\n');wr.writeheader();wr.writerows(identities)
write(OUT/'FUNC1_PREFLIGHT.json',{'PASS':True,'checks':checks,'parent_SHA':PARENT,'science_executions':0,'versions':{'Python':sys.version,'NumPy':np.__version__},'source_index_is_parent_blob':True,'R4_percent':0})
print('PREFLIGHT_PASS',len(checks),'Git inputs',len(inputs),'vectors',len(vec),'scientific executions0')
if __name__=='__main__':pass
