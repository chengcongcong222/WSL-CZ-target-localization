"""Preflight: historical identities, mathematical controls and runtime readiness.
No registered-state information experiment or Monte Carlo is executed here.
"""
from pathlib import Path
import json,csv,hashlib,subprocess,ast,sys,math,shutil,inspect
import numpy as np
import scipy,psutil,matplotlib
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T3_W1');PARENT='7d36d15ecbc9a3435c65e114585328a3bbbf94cc'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def git(*a):return subprocess.check_output(['git','-c','safe.directory='+ROOT.as_posix(),*a],cwd=ROOT).decode().strip()
checks=[]
def ck(name,value):
    checks.append({'check':name,'PASS':bool(value)})
    if not value:
        save(LOCAL/'PREFLIGHT_STOP.json',{'failed':name,'checks':checks,'main_calls':0,'cold_calls':0,'STOP':True});raise RuntimeError('PREFLIGHT_STOP '+name)

def main():
    ck('parent_HEAD',git('rev-parse','HEAD')==PARENT);ck('remote_parent',git('-c','http.version=HTTP/1.1','ls-remote','origin','refs/heads/main').split()[0]==PARENT);ck('no_tracked_or_staged_changes',not git('diff','--name-only') and not git('diff','--cached','--name-only'));ck('no_science_marker',not (OUT/'EXECUTION_STARTED.json').exists())
    names=['r4_a1_new_dynamic.py','r4_a1_new_estimator.py','r4_speed_diagnostic.py','r4_speed_bias_estimator.py','r4_speed_information.py','r4_speed_diagnostic_audit.py','r4_speed_information_audit.py','r4_aux_gate2a.py','r4_aux_gate2b.py','aux_d1_generate.py','aux_d1_model.py','aux_d1_evaluate.py','aux_d1_audit.py']
    folders={'R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC':['A1_NEW_TRUTH_PANEL.csv','A1_NEW_DESIGN_FREEZE.json','A1_NEW_DECISION.json','A1_NEW_CASE_SUMMARY.csv','GPT_SYNC.md'],
    'R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC':['GPT_SYNC.md','SPEED_DIAGNOSTIC_DECISION.json'],
    'R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY':['GPT_SYNC.md','SPEED_INFO_DECISION.json'],
    'R4_AUX_D1_PAIRED_HORIZONTAL_INCREMENT':['GPT_SYNC.md','DESIGN_FREEZE.json','DECISION.json','PUBLIC_OBSERVATION_MANIFEST.json','PRIVATE_GENERATOR_MANIFEST.json','CALIBRATION_THRESHOLDS.json'],
    'R4_AUX_GATE2A_NONIDEAL_ERROR_BUDGET':['AUX_GATE2A_DESIGN_FREEZE.json','AUX_GATE2A_DECISION.json','GPT_SYNC.md'],
    'R4_AUX_GATE2B_JOINT_ERROR_BUDGET':['AUX_GATE2B_DESIGN_FREEZE.json','AUX_GATE2B_DECISION.json','GPT_SYNC.md','AUX_CLIENT_REQUIREMENTS_UPDATED.md'],
    'R4_AUX_CLIENT_CONFIRMATION':['AUX_CLIENT_CONFIRMATION_MATRIX.csv','AUX_ARCHITECTURE_DECISION_AFTER_CLIENT.md','GPT_SYNC.md'],
    'R4_T2_PREC3_FUNC2_CONTINUOUS_MODAL_REFERENCE':['GPT_SYNC_FUNC2.md','T2_NUMERICAL_ROUTE_GO_NO_GO.md','FUNC2_FINAL_ADJUDICATION.json','STOP.json'],
    'R4_E1_G0_FREQUENCY_INFORMATION':['GPT_SYNC.md','E1_G0_DECISION.json','E1_G0_DESIGN_FREEZE.json']}
    # E1 source folder is discovered by exact filename only if that historical alias is absent.
    for folder,files in folders.items():
        for name in files:
            p=ROOT/'results'/folder/name
            ck('history_file_present/'+str(p.relative_to(ROOT)),p.is_file());names.append(p.relative_to(ROOT).as_posix())
    names.append('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv');query=''.join(PARENT+':'+p+'\n' for p in names).encode('utf-8');batch=subprocess.run(['git','-c','safe.directory='+ROOT.as_posix(),'cat-file','--batch'],cwd=ROOT,input=query,capture_output=True,check=True).stdout;pos=0;inputs=[];historical=[]
    for name in names:
        end=batch.index(b'\n',pos);head=batch[pos:end].split();ck('parent_blob_exists/'+name,len(head)==3 and head[1]==b'blob');n=int(head[2]);blob=batch[end+1:end+1+n];pos=end+1+n+1;physical=(ROOT/name).read_bytes()
        if name=='results/R4_MASTER/R4_EVIDENCE_LEDGER.csv':
            ck('ledger_Git_checkout_content',physical.replace(b'\r\n',b'\n')==blob);ck('ledger_original_physical_bytes',physical==(LOCAL/'LEDGER_BEFORE.bin').read_bytes());blob=physical;n=len(blob)
        else:ck('parent_content_identity_including_Git_EOL/'+name,physical==blob or physical.replace(b'\r\n',b'\n')==blob)
        text=physical.decode('utf-8-sig')
        if name.endswith('.py'):ast.parse(text)
        if name.endswith('.json'):
            content=json.loads(text);historical.append({'source':name,'top_level_fields':list(content),'historical_bindings_count':len(content.get('historical_bindings',{})),'scientific_fields':{k:v for k,v in content.items() if k not in ['historical_bindings','bindings','inputs','files'] and isinstance(v,(str,int,float,bool,type(None)))}})
        if name.endswith('.csv'):list(csv.DictReader(text.splitlines()))
        inputs.append({'path':name,'sha256':hashlib.sha256(physical).hexdigest(),'bytes':len(physical),'Git_blob_SHA256':hashlib.sha256(blob).hexdigest(),'checkout_EOL_only':physical!=blob})
    from t3_w1_model import Model,initialize
    import t3_w1_generate as generate
    import t3_w1_information as information
    from t3_w1_statistics import rank_quantile,evaluate
    qualities={'Q1':{'bearing_sigma_deg':.05,'nav_sigma_m':25.,'common_bound_deg':.05,'half_diff_bound_deg':.025},'Q2':{'bearing_sigma_deg':.075,'nav_sigma_m':37.5,'common_bound_deg':.075,'half_diff_bound_deg':.0375}}
    geometries={'B0':{'forward_km':0.,'lateral_km':5.},'B1':{'forward_km':0.,'lateral_km':7.},'B2':{'forward_km':2.,'lateral_km':10.},'B3':{'forward_km':3.,'lateral_km':15.}}
    t,p=generate.platform();ck('epochs_time_turn',len(t)==121 and t[0]==0 and t[-1]==1200 and np.all(np.diff(t)==10) and np.linalg.norm(p[60]-[1200,0])<1e-12)
    sys.path.insert(0,str(ROOT))
    import r4_a1_new_dynamic as old
    tt,pp=old.trajectories();ck('historical_platform_exact',np.array_equal(tt,t) and np.array_equal(pp,p))
    # Different, small deterministic geometry: algebra controls, no registeredH-panel and no noise draws.
    indices=np.array([0,20,60,90,120]);t=t[indices];P=np.stack([p[indices],p[indices]+[900,4000]],axis=1);y=np.zeros((5,2));q=np.array([.4,.6,.35,.6,.6,.4])
    for quality in qualities.values():
        m=Model(t,P,y,quality);f,J=m.evaluate(q);numeric=np.column_stack([(m.fun(q+np.eye(6)[j]*1e-6)-m.fun(q-np.eye(6)[j]*1e-6))/(2e-6) for j in range(6)]);err=float(np.max(abs(J-numeric)));ck('augmented_NLL_full_Jacobian_'+str(quality['bearing_sigma_deg']),err<2e-5)
        h=m.low+m.span*q;r,th,v,ps,c,d=h;th,ps=np.deg2rad([th,ps]);target=P[0,0]+1000*r*np.array([np.cos(th),np.sin(th)])+t[:,None]*v*np.array([np.cos(ps),np.sin(ps)]);D=target[:,None,:]-P;mu=np.arctan2(D[:,:,1],D[:,:,0])+np.deg2rad(c+d*np.array([-1,1]));rho2=np.sum(D*D,axis=-1);base=np.deg2rad(quality['bearing_sigma_deg'])**2;variance=base+quality['nav_sigma_m']**2/rho2;direct=.5*np.sum(mu*mu/variance+np.log(variance/base));ck('NLL_normalization_identity_'+str(quality['bearing_sigma_deg']),abs(.5*(f@f)-direct)<1e-7)
    state=np.array([52000.,700.,1.8,-.1,.0001,-.00003]);err1,err2=information.finite_difference_check(state,t,P,qualities['Q1']);ck('independent_Cartesian_mean_Jacobian',err1<1e-7);ck('independent_Cartesian_variance_Jacobian',err2<1e-6)
    ck('failure_INF_quantile',math.isinf(rank_quantile([1,2,math.inf],.95)));ck('nearest_rank_not_linear',rank_quantile([1,2,3,100],.5)==2)
    fake=[{'config':'Q1_B0','quality':'Q1','geometry':'B0','case_id':'UNIT','replicate':'0','side':'-1','draw_index':'0','accepted':'False','r_hat_km':'nan','theta_hat_deg':'nan','v_hat_mps':'nan','psi_hat_deg':'nan'}];tr=[{'case_id':'UNIT','r0_km':'52','theta0_deg':'0','v_mps':'2','psi_deg':'0'}];ck('failure_all_errors_INF',all(math.isinf(v) for k,v in evaluate(fake,tr)[0].items() if k in ['range_rel','bearing_deg','speed_abs_mps','speed_rel','heading_deg']))
    ck('runtime_monitor_psutil',psutil.Process().memory_info().rss>0 and psutil.__file__.startswith('D:'));ck('runtime_dependencies_import',scipy.__version__ and matplotlib.__version__)
    for prog in OUT.glob('*.py'):
        code=ast.parse(prog.read_text(encoding='utf-8-sig'));ck('AST/'+prog.name,True)
        if prog.name=='t3_w1_model.py':ck('estimator_no_truth_generator_access',not any(isinstance(n,(ast.Import,ast.ImportFrom)) and any(v in (getattr(n,'module','') or ' '.join(a.name for a in getattr(n,'names',[]))) for v in ['generate','information','main','statistics']) for n in ast.walk(code)))
    panelpath='results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv';panel=list(csv.DictReader((ROOT/panelpath).read_text(encoding='utf-8-sig').splitlines()));ck('registered12exact_cases',[r['case_id'] for r in panel]==['H%02d'%i for i in range(1,13)])
    values=np.array([[float(r[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']] for r in panel]);ck('all_states_in_inherited_physical_domain',np.all(values>=[45,-5,1,-15]) and np.all(values<=[60,5,3,15]));ck('D_disk_reserve',shutil.disk_usage('D:/').free>3*1024**3)
    cfg={'stage':'R4-T3-W1','parent_SHA':PARENT,'authorization':'researchleadexplicitindependentstage;historicalSTOPpreserved','inputs':inputs,'truth_panel':panelpath,'geometries':geometries,'qualities':qualities,'seed':202610151701,'seed_engine':'PCG64(SeedSequence([seed,case_index]));uniforms(500,2)thenstandard_normals(500,121,6);samearraysall8configurations','generation_order':'H01..H12;within Q1/Q2 thenB0..B3;allpublicblockssealedbeforefitting','estimation_order':'Q1/Q2;B0/B1/B2/B3;H01..H12;rep0..499','trials_per_config':6000,'trials_per_case':500,'mirrors_per_case_each':250,'base_innovations':6000,'method_configurations':48000,'platform':'2m/s;600s15degturn;fixedglobalauxoffsetmaintained;121epochs10s;noadditionaldeploymentanglejitter','coordinate':'globalmissionx,y m;atan2 CCWrad;directedbearing;relativeestimateoriginreportedMAINt0;truecaseorigintrueMAINt0','main_pair_identity':'times+reportedMAIN+observedMAINbytesidenticalacrossB0..B3forquality/case;standardizedinnovationssameacrossQ1/Q2','navigation':'independentnode/epoch/axisGaussian;reported=true+nav_sigma*normal;first-orderGaussianangularlikelihoodvar=sigma_bearing^2+sigma_nav^2/rho^2;NLLincludeslogvariance','bias':'bMAIN=c-d,bAUX=c+d;constantperrecord;independentUniformbounds;jointboundedflatpriorMLE;noextragaussianprior;biasinitial0','uniform_prior_Fisher':'harduniformboundisnonregular;noinversevariancetermadded;data-only6parameterFisherprofilesfreebias;knownbias4parametercontrolisidealuppercontrol','state_domain':{'r_km':[45,60],'theta_deg':[-5,5],'speed_mps':[1,3],'heading_deg':[-15,15]},'estimator':{'method':'oneobservedrayintersection+5IRLSinitializer;singleboundedTRFjoint6parameterfit','fallback':'physicaldomaincenter52.5km,0deg,2m/s,0deg;bias0;clippedstrictinsideunitbox','loss':'linear;exactfirstorderGaussianNAVequivalentNLLviaresidualpluslogvariance','max_nfev':160,'ftol':1e-10,'xtol':1e-10,'gtol':1e-10,'bounds':'unitcube;biasboundsdependQ','accept':'solver_success&finite;noaccuracythresholdorresidualtuning','boundary_flag_unit':1e-6,'failures':'INFforallerrors;retainreason;neverredraw,retry,orselecttruthbettercandidate','global_candidates':'NOT_CERTIFIED;singlelocalterminalnotcontinuoussupport'} ,'information_scale':information.SCALE.tolist(),'information':'fullGaussianFmean+Fvariance;nominaltruthonlyforinformation;noMCtruthinitialization;Schurviaorthogonalnuisanceprojection;zeroeigenvarianceforbidden','rank_relative_design_tolerance':1e-10,'effective_velocity_eigen_relative_tolerance':1e-12,'Jacobian_crosscheck':{'scaled_mean_abs':1e-7,'relative_variance_abs':1e-6,'step_scaled':1e-5},'quantiles':'nearest_rankceil(q*n);median/P90/P95/P99;conditionalandINF-unconditional;percase,mirrors,pooledgeomequal,worstcase','scientific_decision':'directionandmagnitudereportingwithout5or10percentperformancegate;conditions,weakdirections,costandstandardgeometrynoveltylimitrequired','information_prerequisite':'Jacobian/formula/identityvalidityonly;nocaseexcludedforweakinformation','cold':'noimportsfrommainscience;H01all4000recordsfromsavedinnovations;all48000metricsandlikelihood;allcase/globalquantilesandpaircounts;all192nominalFishergeometries;noadditionaloptimization','budget':{'science_wall_s':7200,'main_wall_s':5400,'cold_wall_s':1800,'tree_RSS_bytes':4*1024**3,'D_data_bytes':2*1024**3,'Git_results_bytes':64*1024**2,'main_calls':1,'cold_calls':1,'retries':0,'BLAS_threads':1},'forbidden_calls':{'KRAKEN':0,'KRAKENC':0,'FIELD':0,'FIELDS':0,'SCOOTER':0,'T2_recompute':0,'Doppler_feature':0,'depth':0,'waveform':0,'parameter_scan_after_freeze':0},'depth':'NOT_EVALUATED','R4_percent':0,'hardware_navigation_bearing':'UNKNOWN','synchronization_association':'ASSUMED_NOT_VALIDATED','frozen_code':[{'path':p.name,'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(OUT.glob('*.py'))]}
    save(OUT/'T3_W1_DESIGN_FREEZE.json',cfg);save(OUT/'HISTORICAL_CONTRACT_REVIEW.json',historical)
    save(OUT/'PREFLIGHT.json',{'PASS':True,'checks':checks,'source_count':len(inputs),'formal_main_calls':0,'formal_cold_calls':0,'versions':{'Python':sys.version,'NumPy':np.__version__,'SciPy':scipy.__version__,'psutil':psutil.__version__,'matplotlib':matplotlib.__version__},'initial_state':read(LOCAL/'INITIAL_STATE.json'),'R4_percent':0})
    with (OUT/'SOURCE_IDENTITIES.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['path','sha256','bytes','Git_blob_SHA256','checkout_EOL_only'],lineterminator='\n');w.writeheader();w.writerows(inputs)
    print('PREFLIGHT_PASS',len(checks),'historical sources',len(inputs),'formal main/cold0')
if __name__=='__main__':main()
