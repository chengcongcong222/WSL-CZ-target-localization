from pathlib import Path
import json,hashlib,subprocess,sys,ast,shutil,struct,datetime
import numpy as np
from scipy.linalg import eigh
ROOT=Path('D:/ProjectStorage/WSL-CZ/T2_PREC3/worktree');LOCAL=ROOT.parent;OUT=ROOT/'results/R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION'
PARENT='d5c0da07c17b12c6e169713a23caf4f2bcfd317c'
def sha(b):return hashlib.sha256(b).hexdigest()
def save(p,v):p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
records=json.loads((LOCAL/'HISTORICAL_READ_INDEX.json').read_text(encoding='utf-8'));checks=[];normalized=[];inputs=[]
def ck(n,c):
 checks.append({'check':n,'PASS':bool(c)})
 if not c:raise RuntimeError('PRECHECK_STOP '+n)
ck('fixed_parent',subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip()==PARENT)
# Every historical context and input is bound to the parent Git blob, not a recreated source.
queries=''.join(PARENT+':'+r['path']+'\n' for r in records).encode()
proc=subprocess.run(['git','cat-file','--batch'],input=queries,capture_output=True,cwd=ROOT,check=True);data=proc.stdout;offset=0
for i,r in enumerate(records):
 end=data.index(b'\n',offset);header=data[offset:end].split();size=int(header[-1]);blob=data[end+1:end+1+size];offset=end+1+size+1;p=ROOT/r['path'];disk=p.read_bytes()
 ck('history_parent_content/'+r['path'],disk==blob or disk.replace(b'\r\n',b'\n')==blob.replace(b'\r\n',b'\n'))
 if disk!=blob:
  # Initial checkout representation only, before scientific freeze; canonical originals untouched.
  p.write_bytes(blob);normalized.append(r['path'])
 inputs.append({'path':r['path'],'sha256':sha(blob),'bytes':size})
 if i%100==0:print('PREFLIGHT_HISTORY',i,flush=True)
 for extra in ['t2_prec2_readonly.py','t2_prec2_cold_audit.py']:pass
base=OUT.parent/'R4_T2_PREC2_ESTD_EIGENVALUE_PRECISION_FEASIBILITY';cfg2=json.loads((base/'T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json').read_text(encoding='utf-8'))
for e in cfg2['archive_sets']:
 for k,rel in e['paths'].items():ck('archive18/'+rel,sha((ROOT/rel).read_bytes())==e['sha256'][k])
files=list(OUT.glob('*.py'))
for p in files:ast.parse(p.read_text(encoding='utf-8'));ck('AST/'+p.name,True)
sys.path.insert(0,str(OUT));import prec3_main as m
# Pure algebra prechecks only: no target SSP, registered roots, loss or physical field executions.
D=np.array([2.,3.,4.]);E=np.array([.1,.2]);MD=np.array([1.,1.2,1.3]);ME=np.array([.02,.03]);A=np.diag(D)+np.diag(E,1)+np.diag(E,-1);M=np.diag(MD)+np.diag(ME,1)+np.diag(ME,-1);roots=eigh(A,M,eigvals_only=True)
for x in [0.,2.1,2.8,5.]:ck('synthetic_inertia/'+str(x),m.greater(D,E,MD,ME,x)==sum(roots>x))
v=np.array([1.,-.2,.3]);ck('synthetic_matvec',np.allclose(m.tri_mv(D,E,v),A@v))
ck('synthetic_Richardson',abs(m.extrap([3+.2*x+.1*x*x+.03*x**3 for x in [1.,.25,.0625,.015625]])[0]-3)<1e-13)
# New task files must be the only dirt after checkout representation is normalized.
status=subprocess.check_output(['git','status','--porcelain=v1','--untracked-files=all'],cwd=ROOT).decode('utf-8')
ck('isolated_worktree_no_historical_changes',not subprocess.check_output(['git','diff','--binary'],cwd=ROOT).strip())
ck('only_new_stage_untracked',all(l.startswith('?? results/R4_T2_PREC3_COMPLEX_WAVENUMBER_AND_FIELD_VALIDATION/') for l in status.splitlines()))
ck('D_reserve',shutil.disk_usage('D:/').free>1024**3)
source=json.loads((LOCAL/'OFFICIAL_SOURCE_INDEX.json').read_text(encoding='utf-8'))
for rr in source:ck('official_pinned/'+rr['url'],sha(Path(rr['path']).read_bytes())==rr['sha256'])
bins=[]
for n in ['kraken.exe','krakenc.exe']:
 rel='tools/acoustics_toolbox/atWin10/at/bin/'+n;p=ROOT/rel;blob=subprocess.check_output(['git','show',PARENT+':'+rel],cwd=ROOT);ck('binary_parent/'+n,p.read_bytes()==blob);bins.append({'path':rel,'sha256':sha(blob),'bytes':len(blob),'build_correspondence':'NOT_ESTABLISHED'});inputs.append(bins[-1])
save(OUT/'PREC3_PREFLIGHT.json',{'PASS':True,'checks':checks,'parent_SHA':PARENT,'old_checkout_line_endings_only':normalized,'new_physical_model_tests':0,'versions':{'python':sys.version,'numpy':np.__version__},'binary_records':bins,'R4_percent':0})
save(OUT/'PREC3_INPUT_IDENTITIES.json',inputs)
prec2rows=list(m.rows(base/'PREC2_LOSSLESS_SUBSET_FDM.csv'))
fixed=[{'freq_hz':int(r['freq_hz']),'mode':int(r['mode'])} for r in prec2rows];ck('fixed15',len(fixed)==15)
design={'stage':'T2-PREC3','parent_SHA':PARENT,'authorization':'Research lead explicit approval of this new stage; no historical STOP revoked','inputs':inputs,'external_sources':source,'A_grids':[20000,40000,80000,160000],'A_method':'P1 consistent mass Galerkin4-point quadrature; generalized Sturm index44bisections; shift-invert eigsh; lambda Richardson h2/h4/h6','A_fixed_modes':fixed,'A_real_tolerance':1e-4/60000,'A_residual':'mass-dual Ritz residual scaled by2k, a discrete bound only; not continuous/full-field proof','A_convergence':'R4 coarse/fine and R6 vs independent PREC2;3 archived continuous roots plus fresh cold3','B_grids':[80000,160000],'loss_scales':[1.,.5,.25],'B_method':'independent first perturbation integral and directly solved complex Galerkin; paired complex-lambda shifts relative to same-grid lossless lambda extrapolated h2, added to verified lossless lambdaR6; not a splice with MOD imaginary parts','B_real_imag_tolerance':1e-4/60000,'B_source_operator':'public pinned CRCI at251nodes: c+i alpha c2/omega; interpolate complex c linearly; V=omega2/c_complex2; outgoing exp(-ikR) with negative Im(k)','B_loss_scaling':'each halving reduces error at least25percent unless within1e-11 m-1 registered numerical floor; floor-dominated test does not prove quadratic asymptotics','B_independent_FD_modes':[1,481,961],'modal_function_comparison_tolerance':.002,'C_scope':'only if A+B pass;60 fixed sentinel/resource/frequency diagnostics; both stored modal grids; C1 only15 certified k replacements; C2 only3 verified201Hz functions; no arbitrary gain fitting','C_geometry':'N2 P01..P05 all4resources 201235283; P04 retained','C3':'NOT_EVALUATED_FULL_REFERENCE_NOT_AVAILABLE_AND_ALL_MODE_CALCULATION_NOT_AUTHORIZED','C_gate':'no complete reference means no full-field certification regardless of partial changes;0.002 remains old spatial threshold, never rescored as old PASS','cold_tolerances':{'root':1e-4/60000,'scalar_field_QR':1e-10,'reconstructed_R6':1e-12},'budget':{'whole_wall_s':3600,'main_wall_s':3000,'cold_wall_s':600,'tree_RSS_bytes':4*1024**3,'local_bytes':2*1024**3,'git_result_bytes':64*1024**2,'scientific_processes':1,'BLAS_threads':1,'main_calls':1,'cold_calls':1,'retries':0,'new_KRAKEN_KRAKENC_FIELD_FIELDS_SCOOTER_MC_waveform_DFT_FIM_localization':0},'frozen_code':[{'path':p.name,'sha256':sha(p.read_bytes()),'bytes':p.stat().st_size} for p in files],'historical_results':'UNCHANGED','R4_percent':0}
save(OUT/'PREC3_DESIGN_FREEZE.json',design);save(LOCAL/'WORKTREE_INITIALISATION.json',{'checkout_newline_normalized_to_parent_blob_only':normalized,'canonical_untouched':True})
print('PRECHECK_PASS',len(checks),'inputfiles',len(inputs),'frozen_code',len(files))
