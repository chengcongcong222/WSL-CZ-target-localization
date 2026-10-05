"""Independent reconstruction of saved states, L2 objective and physical FIM inverse."""
from pathlib import Path
from collections import defaultdict
import argparse,csv,json,hashlib,math,subprocess,ast,xml.etree.ElementTree as ET
import numpy as np
OUT=Path('results/R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY')
OLD=Path('results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC')
DIAG=Path('results/R4_A1_NEW_SPEED_BOTTLENECK_DIAGNOSTIC')
SCALE=np.array([50000.,50000.,2.,2.])
METRICS=['range_error','bearing_error_deg','speed_error','heading_error_deg']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,data):
 with Path(p).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def dump(p,d):Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def flag(x):return x is True or x=='True'
def wrap(x):return (x+np.pi)%(2*np.pi)-np.pi
def nq(vals,q):return float(np.sort(np.asarray(vals,float))[math.ceil(len(vals)*q)-1])
def initial(t,n,b):
 u=np.stack([np.cos(b),np.sin(b)],axis=-1);delta=n[:,1]-n[:,0];det=u[:,0,0]*u[:,1,1]-u[:,0,1]*u[:,1,0]
 l1=(delta[:,0]*u[:,1,1]-delta[:,1]*u[:,1,0])/det;l2=(delta[:,0]*u[:,0,1]-delta[:,1]*u[:,0,0])/det
 xy=n[:,0]+l1[:,None]*u[:,0];valid=(abs(det)>1e-12)&(l1>0)&(l2>0)&np.isfinite(xy).all(axis=1);tt=t[valid];points=xy[valid]
 T=np.column_stack([np.ones(len(tt)),tt/1200]);coef=np.linalg.lstsq(T,points,rcond=None)[0]
 for _ in range(5):
  norm=np.linalg.norm(points-T@coef,axis=1);scale=max(10.,float(np.median(norm)));sw=np.sqrt(np.minimum(1.,1.5*scale/np.maximum(norm,1e-15)));coef=np.linalg.lstsq(T*sw[:,None],points*sw[:,None],rcond=None)[0]
 return np.r_[coef[0],coef[1]/1200]
def Hphysical(t,n,s,sigma):
 d=(s[:2]+t[:,None]*s[2:])[:,None,:]-n;rho=np.sum(d*d,axis=-1);gx=-d[:,:,1]/rho;gy=d[:,:,0]/rho
 return np.stack([gx,gy,gx*t[:,None],gy*t[:,None]],axis=-1).reshape(-1,4)/sigma

def audit():
 checks=[]
 def check(name,ok):
  if not ok:raise AssertionError(name)
  checks.append(dict(check=name,status='PASS'))
 def close(name,a,b,rtol=1e-6,atol=1e-8):check(name,bool(np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=rtol,atol=atol,equal_nan=False)))
 p=load(OUT/'SPEED_INFO_AUDIT_POLICY.json')
 for n,h in {**p['bindings'],**p['historical_bindings']}.items():check('binding:'+n,sha(n)==h)
 src=rows(DIAG/'SPEED_VARIANT_RUN_RESULTS.csv');source={(x['anchor'],int(x['draw_index'])):x for x in src if x['variant']=='D3'}
 check('E0_exact_D3_rows',rows(OUT/'E0_D3_IDENTITY.csv')==[x for x in src if x['variant']=='D3'])
 run=rows(OUT/'MLE_RUN_RESULTS.csv');cr=rows(OUT/'CRLB_GEOMETRY_RUN_RESULTS.csv');paired=rows(OUT/'MLE_PAIRED_BASIN_DIAGNOSTICS.csv');matrix=np.load(OUT/'CRLB_MATRIX_RECORDS.npz');matrix={k:matrix[k] for k in matrix.files}
 check('record_counts',len(run)==36000 and len(cr)==60000 and len(paired)==12000);check('matrix_shapes',matrix['F'].shape==matrix['C'].shape==(60000,4,4))
 ri={(x['anchor'],int(x['draw_index']),x['variant']):x for x in run};ci={(x['anchor'],int(x['draw_index']),x['window']):x for x in cr};pi={(x['anchor'],int(x['draw_index'])):x for x in paired}
 check('unique_keys',len(ri)==len(run) and len(ci)==len(cr) and len(pi)==len(paired))
 panel=rows(OLD/'A1_NEW_TRUTH_PANEL.csv');anchors=load(OLD/'A1_NEW_DESIGN_FREEZE.json')['anchors'];cold=defaultdict(list);coldcr=defaultdict(list)
 t=np.arange(121)*10.;pt=np.maximum(t-600,0);main=np.column_stack([2*np.minimum(t,600)+2*pt*np.cos(np.pi/12),2*pt*np.sin(np.pi/12)])
 windows=[('PREFIX_300',t<=300),('STRAIGHT_600',t<=600),('PREFIX_900',t<=900),('FULL_1200',t<=1200),('POST_TURN',t>600)]
 for a in ['A','B']:
  draws=np.load(OLD/f'A1_NEW_DRAWS_{a}.npz');uniforms=draws['uniforms'];normals=draws['normals'];anchor=anchors[a];sigma=math.radians(anchor['sigma_deg'])
  for ix in range(6000):
   truth=panel[ix//500];side=-1 if ix%2==0 else 1;u=uniforms[ix];z=normals[ix];beta=math.radians((2*u[2]-1)*anchor['beta_deg']);offset=anchor['baseline_m']*np.array([-math.sin(beta),side*math.cos(beta)]);nodes=np.stack([main,main+offset],axis=1)
   r=float(truth['r0_km'])*1000;theta=math.radians(float(truth['theta0_deg']));v=float(truth['v_mps']);psi=math.radians(float(truth['psi_deg']));true=np.array([r*math.cos(theta),r*math.sin(theta),v*math.cos(psi),v*math.sin(psi)])
   delta=(true[:2]+t[:,None]*true[2:])[:,None,:]-nodes;ideal=np.arctan2(delta[:,:,1],delta[:,:,0]);common=(2*u[0]-1)*anchor['common_deg'];diff=(2*u[1]-1)*anchor['half_diff_deg'];bias=np.deg2rad([common-diff,common+diff]);angles=(ideal+bias+sigma*z[:,:2])-bias
   init=initial(t,nodes,angles);name=a+str(ix)+':';states={}
   for variant in ['E0','E1','E2']:
    x=ri[a,ix,variant];check(name+variant+'metadata',x['case_id']==truth['case_id'] and int(x['realization'])==ix%500 and int(x['mirror'])==side)
    if variant=='E0':check(name+'E0_identity',all(x[k]==source[a,ix][k] for k in ['failed','failure_reason','x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps']+METRICS+['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg']))
    if flag(x['failed']):
     check(name+variant+'failure_INF',all(math.isinf(float(x[k])) for k in METRICS+['absolute_speed_error_mps','signed_speed_error_mps']));cold[a,truth['case_id'],variant].append(dict(failed=True,**{k:math.inf for k in METRICS+['absolute_speed_error_mps','signed_speed_error_mps']}));continue
    s=np.array([float(x[k]) for k in ['x_hat_m','y_hat_m','vx_hat_mps','vy_hat_mps']]);states[variant]=s;rh=np.linalg.norm(s[:2]);vh=np.linalg.norm(s[2:]);th=math.atan2(s[1],s[0]);ph=math.atan2(s[3],s[2])
    vals=dict(range_error=abs(rh-r)/r,bearing_error_deg=abs(math.degrees(wrap(th-theta))),speed_error=abs(vh-v)/v,heading_error_deg=abs(math.degrees(wrap(ph-psi))),absolute_speed_error_mps=abs(vh-v),signed_speed_error_mps=vh-v)
    for k,value in vals.items():close(name+variant+k,x[k],value)
    cold[a,truth['case_id'],variant].append(dict(failed=False,**vals))
    if variant!='E0':
     check(name+variant+'scope',x['scope']==('OBSERVATION_ONLY_DEVELOPMENT' if variant=='E1' else 'ORACLE_INITIALIZATION_DIAGNOSTIC_ONLY'));close(name+variant+'initializer',[x[k] for k in ['initial_x_m','initial_y_m','initial_vx_mps','initial_vy_mps']],init if variant=='E1' else true,rtol=1e-8,atol=1e-6)
     dd=(s[:2]+t[:,None]*s[2:])[:,None,:]-nodes;rr=(wrap(np.arctan2(dd[:,:,1],dd[:,:,0])-angles)/sigma).ravel();close(name+variant+'cost',x['cost'],.5*(rr@rr));close(name+variant+'RMS',x['residual_RMS_deg'],np.sqrt(np.mean(rr*rr))*anchor['sigma_deg'])
     close(name+variant+'optimality',x['optimality'],np.max(abs((Hphysical(t,nodes,s,sigma)*SCALE).T@rr)),rtol=1e-4,atol=1e-6);check(name+variant+'cap_success',flag(x['solver_success']) and int(x['nfev'])<=100)
   if 'E1' in states and 'E2' in states:
    x,y=ri[a,ix,'E1'],ri[a,ix,'E2'];pr=pi[a,ix];check(name+'paired_success',flag(pr['both_success']));close(name+'paired_speed',pr['absolute_signed_speed_difference_mps'],abs(float(x['v_hat_mps'])-float(y['v_hat_mps'])));close(name+'paired_cost',pr['absolute_cost_difference'],abs(float(x['cost'])-float(y['cost'])));close(name+'paired_state',pr['max_scaled_state_difference'],max(abs(states['E1']-states['E2'])/SCALE))
   for wi,(window,mask) in enumerate(windows):
    x=ci[a,ix,window];mi=int(x['matrix_index']);check(name+window+'index_epochs',mi==((0 if a=='A' else 6000)+ix)*5+wi and int(x['n_epochs'])==int(mask.sum()) and int(x['n_measurements'])==2*int(mask.sum()))
    H=Hphysical(t[mask],nodes[mask],true,sigma);J=H*SCALE;F=H.T@H;sv=np.linalg.svd(J,compute_uv=False);rank=int(np.sum(sv/sv[0]>1e-10));fsv=np.linalg.svd(F,compute_uv=False)
    check(name+window+'rank',int(x['rank_scaled_H'])==rank and flag(x['full_rank'])==(rank==4));close(name+window+'F',matrix['F'][mi],F,rtol=1e-8,atol=1e-12)
    close(name+window+'Jsv',[x[f'scaled_H_singular_{j+1}'] for j in range(4)],sv);close(name+window+'Fsv',[x[f'physical_F_singular_{j+1}'] for j in range(4)],fsv);close(name+window+'Jcondition',x['condition_scaled_H'],sv[0]/sv[-1]);close(name+window+'Fcondition',x['condition_physical_F'],fsv[0]/fsv[-1]);close(name+window+'ratio',x['smallest_largest_scaled_H'],sv[-1]/sv[0])
    if rank==4:C=np.linalg.inv(F);close(name+window+'inverseF',matrix['C'][mi],C,rtol=1e-6,atol=1e-6)
    else:
     check(name+window+'rank_not_bound',x['covariance_kind']=='PSEUDOINVERSE_NOT_FINITE_BOUND');C=matrix['C'][mi]
    gg=[np.array([true[0]/r,true[1]/r,0,0]),np.array([-true[1]/r**2,true[0]/r**2,0,0]),np.r_[0.,0.,true[2:]/v],np.array([0,0,-true[3]/v**2,true[2]/v**2])]
    variances=[max(0.,float(g@C@g)) for g in gg] if rank==4 else [math.inf]*4;sr,st,svv,sp=np.sqrt(variances);schur=F[2:,2:]-F[2:,:2]@np.linalg.solve(F[:2,:2],F[:2,2:])
    vals=dict(sigma_speed_mps=svv,sigma_speed_rel=svv/v,speed_rel_1645=1.645*svv/v,speed_rel_196=1.96*svv/v,heading_deg_196=1.96*math.degrees(sp),range_rel_196=1.96*sr/r,bearing_deg_196=1.96*math.degrees(st),speed_variance_mps2=variances[2],velocity_schur_trace=float(np.trace(schur)))
    for k,value in vals.items():close(name+window+k,x[k],value)
    coldcr[a,truth['case_id'],window].append(dict(**vals,full_rank=rank==4))
   if (ix+1)%500==0:print(f'Cold efficiency audit {a} {ix+1}/6000',flush=True)
 summary=rows(OUT/'MLE_CASE_SUMMARY.csv');eff=rows(OUT/'SPEED_EFFICIENCY_SUMMARY.csv');check('summary_counts',len(summary)==len(eff)==72);ei={(x['anchor'],x['case_id'],x['variant']):x for x in eff}
 for x in summary:
  key=x['anchor'],x['case_id'],x['variant'];g=cold[key];name=str(key);nf=sum(y['failed'] for y in g);check(name+'counts',len(g)==int(x['n_runs'])==500 and int(x['n_failed'])==nf);close(name+'failure_rate',x['failure_rate'],nf/500)
  for k in METRICS+['absolute_speed_error_mps']:
   for tag,q in [('median',.5),('P90',.9),('P95',.95),('P99',.99)]:close(name+k+tag,x[k+'_'+tag],nq([y[k] for y in g],q))
  close(name+'speed_alias',x['speed_rel_P95'],x['speed_error_P95']);close(name+'absolute_alias',x['speed_abs_P95_mps'],x['absolute_speed_error_mps_P95'])
  check(name+'gates',flag(x['PROJECT_speed_pass'])==(float(x['speed_error_P95'])<=.1) and flag(x['STRONG_speed_pass'])==(float(x['speed_error_P95'])<=.05) and flag(x['PROJECT_other_metrics_pass'])==all(float(x[m+'_P95'])<=v for m,v in [('range_error',.1),('bearing_error_deg',1),('heading_error_deg',5)]) and flag(x['route_failure_rate_pass'])==(nf/500<=.01))
  signed=np.array([y['signed_speed_error_mps'] for y in g if not y['failed']]);ee=ei[key];check(name+'nvalid',int(x['n_valid'])==int(ee['n_valid'])==len(signed))
  if len(signed)>1:
   stats=dict(mean_speed_bias_mps=float(signed.mean()),median_speed_bias_mps=float(np.median(signed)),speed_std_mps=float(signed.std(ddof=1)),speed_RMSE_mps=float(np.sqrt(np.mean(signed*signed))))
   for k,value in stats.items():close(name+k,x[k],value);close(name+'eff:'+k,ee[k],value)
  else:check(name+'no_finite_stats',x['speed_std_mps']=='NA');stats=dict(speed_std_mps=0)
  variance=np.mean([y['speed_variance_mps2'] for y in coldcr[key[0],key[1],'FULL_1200']]);close(name+'CRLBmean',ee['mean_conditional_local_CRLB_variance_mps2'],variance);eligible=nf==0 and math.isfinite(variance) and stats['speed_std_mps']>0;check(name+'eligible',flag(ee['efficiency_eligible'])==eligible)
  if eligible:close(name+'eta',ee['eta_v'],variance/stats['speed_std_mps']**2);close(name+'bias_std',ee['bias_to_std'],abs(stats['mean_speed_bias_mps'])/stats['speed_std_mps'])
  else:check(name+'eta_NA',ee['eta_v']=='NA')
 crcase=rows(OUT/'CRLB_CASE_SUMMARY.csv');prefix=rows(OUT/'CRLB_TIME_PREFIX.csv');segments=rows(OUT/'CRLB_SEGMENT_SUMMARY.csv');check('CRLBsummary_counts',len(crcase)==24 and len(prefix)==96 and len(segments)==72)
 for x in crcase+prefix+segments:
  g=coldcr[x['anchor'],x['case_id'],x['window']];name=x['anchor']+x['case_id']+x['window'];close(name+'rank_rate',x['full_rank_rate'],sum(y['full_rank'] for y in g)/500);check(name+'500_geometries',int(x['n_saved_geometries'])==500)
  for k in ['sigma_speed_mps','sigma_speed_rel','speed_rel_1645','speed_rel_196','heading_deg_196','range_rel_196','bearing_deg_196','speed_variance_mps2','velocity_schur_trace']:
   vals=[y[k] for y in g]
   for tag,value in [('min',min(vals)),('median',nq(vals,.5)),('P95',nq(vals,.95)),('max',max(vals))]:close(name+k+tag,x[k+'_'+tag],value)
  check(name+'CRLBgates',flag(x['PROJECT_CRLB_equivalent_pass'])==(float(x['speed_rel_196_max'])<=.1) and flag(x['STRONG_CRLB_equivalent_pass'])==(float(x['speed_rel_196_max'])<=.05))
 turns=rows(OUT/'CRLB_INFORMATION_ACCUMULATION.csv');check('12000_information_accumulation',len(turns)==12000)
 for x in turns:
  a=x['anchor'];ix=int(x['draw_index']);s=ci[a,ix,'STRAIGHT_600'];f=ci[a,ix,'FULL_1200'];post=ci[a,ix,'POST_TURN'];F=matrix['F'][int(f['matrix_index'])];Fs=matrix['F'][int(s['matrix_index'])];Fp=matrix['F'][int(post['matrix_index'])]
  check(a+str(ix)+'Fadd',np.allclose(F,Fs+Fp,rtol=1e-10,atol=1e-12));close(a+str(ix)+'Fadd_metric',x['F_additivity_max_abs'],np.max(abs(F-Fs-Fp)));close(a+str(ix)+'variance_gain',x['speed_variance_reduction'],1-float(f['speed_variance_mps2'])/float(s['speed_variance_mps2']));close(a+str(ix)+'trace_gain',x['effective_velocity_information_trace_gain'],float(f['velocity_schur_trace'])/float(s['velocity_schur_trace']))
 d=load(OUT/'SPEED_INFO_DECISION.json');A={v:[x for x in summary if x['anchor']=='A' and x['variant']==v] for v in ['E0','E1','E2']};full=[x for x in crcase if x['anchor']=='A'];passes={v:sum(flag(x['PROJECT_speed_pass']) for x in g) for v,g in A.items()};gap=max(abs(float(x['speed_rel_P95'])-float(next(y for y in A['E2'] if y['case_id']==x['case_id'])['speed_rel_P95'])) for x in A['E1']);above10=sum(float(x['speed_rel_196_max'])>.1 for x in full);above5=sum(float(x['speed_rel_196_max'])>.05 for x in full);rankok=all(float(x['full_rank_rate'])==1 for x in full);rescue=passes['E1']==12 and all(flag(x['PROJECT_other_metrics_pass']) and flag(x['route_failure_rate_pass']) for x in A['E1'])
 verdict='GAUSSIAN_MLE_SPEED_ROUTE_WORTH_FRESH_VALIDATION' if rescue else 'CRLB_NUMERICAL_RANK_NOT_CLOSED' if not rankok else 'CURRENT_1200S_GEOMETRY_INFORMATION_LIMIT_CONFIRMED_BY_CRLB' if passes['E1']<12 and gap<=.02 and above10>0 else 'ESTIMATOR_OR_NONLINEAR_EFFICIENCY_GAP_REMAINS' if passes['E1']<12 and above10==0 else 'OPTIMIZER_OR_INITIALIZER_INEFFICIENCY_PRESENT' if gap>.02 else 'INFORMATION_EFFICIENCY_AUDIT_INCONCLUSIVE'
 check('route_decision',d['route_decision']==verdict and d['PROJECT_speed_cases']==passes and d['CRLB_cases_above_10']==above10 and d['CRLB_cases_above_5']==above5);close('E1_E2_gap',d['E1_E2_max_case_P95_gap'],gap)
 check('zero_credit_no_new_MC',all(d[k]==0 for k in ['new_Monte_Carlo_realizations','new_truth_cases','new_bearing_noise_draws','R4_A1_percent','R4_percent']) and not d['automatic_next_stage'] and d['depth']=='NOT_OPENED')
 start=load(OUT/'SPEED_INFO_EXECUTION_START.json');check('pushed_before_CRLB',start['remote_verified_before_resolves_and_CRLB'] and start['policy_SHA']==d['policy_SHA'] and start['policy_sha256']==sha(OUT/'SPEED_INFO_AUDIT_POLICY.json'));check('policy_parent',subprocess.check_output(['git','rev-parse',d['policy_SHA']+'^'],text=True).strip()==p['parent_SHA'])
 xml=ET.parse(OUT/'SPEED_INFO_TESTS.xml').getroot();check('5tests_PASS',len(xml.findall('.//testcase'))==5 and not xml.findall('.//failure') and not xml.findall('.//error'))
 tree=ast.parse(Path('r4_speed_information.py').read_text(encoding='utf-8'));check('no_RNG_calls',not any(isinstance(x,ast.Call) and isinstance(x.func,ast.Attribute) and x.func.attr in ['default_rng','Generator','PCG64','standard_normal','random'] for x in ast.walk(tree)));e1=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='estimate_gaussian');check('E1_no_truth_API',len(e1.args.args)==4 and all('truth' not in x.arg for x in e1.args.args))
 for n in ['results/R4_MASTER/R4_PLAN.md','results/R4_MASTER/R4_EVIDENCE_LEDGER.csv']:
  before=subprocess.check_output(['git','show',p['parent_SHA']+':'+n]).decode('utf-8').replace('\r\n','\n');check('append_only:'+n,Path(n).read_text(encoding='utf-8').startswith(before))
 return checks

def present():
 d=load(OUT/'SPEED_INFO_DECISION.json');summary=rows(OUT/'MLE_CASE_SUMMARY.csv');cr=rows(OUT/'CRLB_CASE_SUMMARY.csv');prefix=rows(OUT/'CRLB_TIME_PREFIX.csv');eff=rows(OUT/'SPEED_EFFICIENCY_SUMMARY.csv');turns=rows(OUT/'CRLB_INFORMATION_ACCUMULATION.csv')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figs=OUT/'figures';figs.mkdir(exist_ok=True);plt.rcParams['svg.hashsalt']='SPEED_INFO'
 def save(fig,name):
  fig.tight_layout();fig.savefig(figs/(name+'.png'),dpi=150,bbox_inches='tight');fig.savefig(figs/(name+'.svg'),bbox_inches='tight',metadata={'Date':None});plt.close(fig)
  p=figs/(name+'.svg');p.write_text('\n'.join(x.rstrip() for x in p.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 for number in [1,2]:
  fig,axes=plt.subplots(1,2,figsize=(13,5))
  for ax,a in zip(axes,['A','B']):
   for v in (['E0','E1','E2'] if number==1 else ['E1','E2']):
    g=[x for x in summary if x['anchor']==a and x['variant']==v];ax.plot(range(1,13),[100*float(x['speed_rel_P95']) for x in g],'o-',label={'E0':'D3 soft-L1 (identity)','E1':'Gaussian L2 empirical','E2':'Oracle-init L2 empirical'}[v])
   if number==2:
    g=[x for x in cr if x['anchor']==a];ax.plot(range(1,13),[100*float(x['speed_rel_196_max']) for x in g],'s--',label='Local CRLB 1.96 sigma (saved max)')
   ax.axhline(10,color='black',linestyle='--',label='PROJECT 10%');ax.set(title='Anchor '+a+' development'+('; CRLB is not empirical P95' if number==2 else ''),xlabel='Frozen case',ylabel='Relative speed (%)');ax.legend();ax.grid(alpha=.2)
  save(fig,'FIG1_SOFT_L1_VS_GAUSSIAN_MLE' if number==1 else 'FIG2_MLE_VS_LOCAL_CRLB')
 fig,axes=plt.subplots(1,2,figsize=(13,5))
 for ax,a in zip(axes,['A','B']):
  for case in sorted({x['case_id'] for x in cr}):
   g=[next(x for x in prefix if x['anchor']==a and x['case_id']==case and x['window']==w) for w in ['PREFIX_300','STRAIGHT_600','PREFIX_900','FULL_1200']];ax.plot([300,600,900,1200],[100*float(x['speed_rel_196_max']) for x in g],alpha=.65,label=case)
  ax.axhline(10,color='black',linestyle='--');ax.set(title='Anchor '+a+' same saved trajectory prefixes',xlabel='Prefix duration (s)',ylabel='Local 1.96 sigma speed equivalent (%)');ax.grid(alpha=.2)
 save(fig,'FIG3_SPEED_CRLB_PREFIX')
 fig,axes=plt.subplots(1,2,figsize=(12,5))
 for ax,a in zip(axes,['A','B']):
  g=[x for x in cr if x['anchor']==a];fields=['range_rel_196_max','bearing_deg_196_max','speed_rel_196_max','heading_deg_196_max'];gates=[.1,1,.1,5];ax.bar(['Range','Bearing','Speed','Heading'],[max(float(x[k]) for x in g)/gate for k,gate in zip(fields,gates)]);ax.axhline(1,color='black',linestyle='--');ax.set(title='Anchor '+a+' local equivalent / PROJECT target',ylabel='Ratio; no empirical performance claim');ax.grid(axis='y',alpha=.2)
 save(fig,'FIG4_FOUR_METRIC_CRLB')
 lines=['# 速度信息效率审计','',d['route_decision'],'','DEVELOPMENT / INFORMATION EFFICIENCY AUDIT。复用 12000 个保存场景；E0 为 D3 逐 run 恒等，E1 只将 loss 改为 linear，E2 为 ORACLE_INITIALIZATION_DIAGNOSTIC_ONLY。新增 MC、truth case、bearing noise draw 均为 0；R4-A1=0%，R4=0%，depth NOT_OPENED。','',
 '| Anchor | Variant | Worst speed P95 % | PROJECT speed | Worst range P95 % | Bearing P95 deg | Heading P95 deg | Max failure rate |','|---|---|---:|---:|---:|---:|---:|---:|']
 for a in ['A','B']:
  for v in ['E0','E1','E2']:
   g=[x for x in summary if x['anchor']==a and x['variant']==v];lines.append(f"| {a} | {v} | {100*max(float(x['speed_rel_P95']) for x in g):.6f} | {sum(flag(x['PROJECT_speed_pass']) for x in g)}/12 | {100*max(float(x['range_error_P95']) for x in g):.6f} | {max(float(x['bearing_error_deg_P95']) for x in g):.6f} | {max(float(x['heading_error_deg_P95']) for x in g):.6f} | {max(float(x['failure_rate']) for x in g):.6f} |")
 lines+=['','## 模型与统计口径','',
 '四状态模型使用真实节点位置及已移除生成 bias 的同一批随机 bearing noise，共 242 measurements。E1 保留原 ray-intersection + IRLS initializer、状态缩放与解析 Jacobian，TRF、max_nfev=100、ftol/xtol/gtol=1e-10。没有裁点、truth warm start、truth 邻域边界、multistart、额外先验或事后调参。E2 唯一区别是真实 Cartesian 初值，不能形成 estimator claim。','',
 'F=HᵀH/sigma_theta²；原状态缩放 measurement Jacobian 的 s/smax>1e-10 为数值秩门限。已知 sigma，不用 fitted residual variance rescaling。physical F/C、奇异值及条件数按 matrix_index 保存在 CRLB_MATRIX_RECORDS.npz 与 CRLB_GEOMETRY_RUN_RESULTS.csv。满秩时缩放 SVD 恢复的 C 为 F 的逆，非满秩保留 pseudoinverse，但不确定度必须报告无穷，不能给有限下界。','',
 'CRLB 覆盖每个保存布放角及双镜像。每 case 报告 min/median/P95/max，route 使用 max。速度梯度 [0,0,vx/v,vy/v]，报告 1σ、1.645σ、1.96σ；主要比较为 1.96σ。所有等效值标为 LOCAL_LINEAR_GAUSSIAN_EQUIVALENT / NOT_EMPIRICAL_P95 / NOT_GLOBAL_GUARANTEE。r/theta/psi 使用 delta method，只作交叉诊断，不重定义 Gate。','',
 '效率比 eta_v 为平均条件局部 CRLB variance / empirical signed speed variance；sample std 使用 ddof=1，偏差、median、RMSE 同时报告。仅全部 run 有效且方差有限时定义；它是局部无偏方差参考，不是总误差 bound。失败保留 INF 并进入分位数/Gate；有限样本的偏差统计显式标为描述性统计。','',
 f"E1/E2 最大逐 case P95 差为 {100*d['E1_E2_max_case_P95_gap']:.9f} 个百分点，冻结门限为 2 个百分点。逐 run 状态、速度与 cost 差另存 paired basin diagnostics；P95 接近不等于证明全局 MLE 覆盖。",'']
 for a in ['A','B']:
  g=[x for x in cr if x['anchor']==a];lines.append(f"{a}: 最坏局部 1.96σ 速度等效值={100*max(float(x['speed_rel_196_max']) for x in g):.6f}%；超过 10% 的 case={sum(float(x['speed_rel_196_max'])>.1 for x in g)}/12，超过 5%={sum(float(x['speed_rel_196_max'])>.05 for x in g)}/12。")
  g=[x for x in eff if x['anchor']==a and x['variant']=='E1']
  for k in ['mean_speed_bias_mps','median_speed_bias_mps','speed_std_mps','speed_RMSE_mps','eta_v']:
   vals=[float(x[k]) for x in g if x[k]!='NA'];lines.append(f"{a} E1 {k}: case min..max={min(vals):.9f}..{max(vals):.9f}。" if vals else f"{a} E1 {k}: NA。")
 lines+=['','## 信息累积与转向窗口','',
 'T=300/600/900/1200 s 只截取同一 trajectory 做 analytic FIM，不运行 prefix optimizer，也不外推 1800/2400 s。STRAIGHT 包含 t<=600 的 61 epochs；POST_TURN 为 t>600 的 60 epochs；FULL 为 121 epochs。边界不重复，F_full=F_straight+F_post，所有窗口沿用同一 x0 和全局 t。','',
 '速度有效信息使用 Schur complement Fvv-Fvp Fpp^{-1} Fpv 消去初始位置。full 相对 straight 的增益同时包含新增时长、观测数和现有改变的几何，不是 isolated 15° turn 因果效应；没有额外直航 counterfactual，不能单独量化转向效应。','']
 for a in ['A','B']:
  for w in ['PREFIX_300','STRAIGHT_600','PREFIX_900','FULL_1200']:
   g=[x for x in prefix if x['anchor']==a and x['window']==w];lines.append(f"{a} {w}: 最坏局部 1.96σ 速度等效值={100*max(float(x['speed_rel_196_max']) for x in g):.6f}%。")
  g=[x for x in turns if x['anchor']==a];gain=[float(x['speed_variance_reduction']) for x in g];trace=[float(x['effective_velocity_information_trace_gain']) for x in g];lines.append(f"{a} full/straight: 局部速度方差下降 min/median/max={min(gain):.9f}/{nq(gain,.5):.9f}/{max(gain):.9f}；Schur trace 倍数={min(trace):.9f}/{nq(trace,.5):.9f}/{max(trace):.9f}。")
 lines+=['','## 路线与作用域','',
 'Outcome II 只在 E1 非全部通过、全部 primary case E1/E2 P95 差<=2pp、至少一个 primary case 的局部 1.96σ 速度等效值>10% 时成立。准确结论是 local Gaussian CRLB equivalent itself above target；不是所有算法的经验 P95 下界，也不是被动双节点的普遍物理不可辨识。B 不能替代 A。','',
 f"本轮路线判定：{d['route_decision']}；仅建议下一阶段 {d['next_recommended']}。OBSERVATION_DESIGN 对应 R4_A1_NEW_SPEED_OBSERVATION_DESIGN_GATE，需要另行授权。本轮禁止改变 duration/AUX 相对运动/baseline/turn/RMS/truth，不做 Kalman/smoother/Doppler/range-rate/TDOA/depth/acoustics；push 后 STOP。",'',
 '既有 range/bearing STRONG、heading PROJECT 结论保留。本轮为 PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT，不给 application credit。','',chr(96)*3+'json',json.dumps(d,indent=2),chr(96)*3]
 text='\n'.join(lines)+'\n'
 for name in ['SPEED_INFO_REPORT.md','GPT_SYNC.md']:(OUT/name).write_text(text,encoding='utf-8')
 dump('results/R4_MASTER/R4_SPEED_INFORMATION_PROGRESS.json',dict(stage=d['stage'],project_stage='PROJECT_APPLICATION_PRE_RESEARCH',range='STRONG_ESTABLISHED',bearing='STRONG_ESTABLISHED',heading='PROJECT_ESTABLISHED',speed=d['route_decision'],R4_A1_percent=0,R4_percent=0,depth='NOT_OPENED',recommended_next=d['next_recommended'],automatic_next_stage=False,audit_status='PENDING_RESEARCH_LEAD_INDEPENDENT_AUDIT'))
 with Path('results/R4_MASTER/R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## Speed information efficiency execution\n\n'+d['route_decision']+'; DEVELOPMENT ONLY, pending independent audit. [Report](../R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY/SPEED_INFO_REPORT.md). [Active status](R4_SPEED_INFORMATION_PROGRESS.json). New MC=0; R4-A1/R4=0%; depth NOT_OPENED. Next is a recommendation only; push then STOP.\n')
 with Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:csv.writer(f).writerow(['A1-NEW-SPEED-INFO','../R4_A1_NEW_SPEED_INFORMATION_EFFICIENCY/SPEED_INFO_DECISION.json',d['route_decision'],'Saved scenes; D3 model; L2 vs soft-L1; local CRLB; no new MC','PENDING_RESEARCH_LEAD_AUDIT;NO_SCIENTIFIC_CREDIT',0])

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--finalize',action='store_true');ap.add_argument('--verify-seal',action='store_true');args=ap.parse_args()
 if args.finalize:
  assert not (OUT/'VALIDATION.json').exists(),'One finalize only'
  checks=audit();present();write(OUT/'SPEED_INFO_INDEPENDENT_CHECKS.csv',checks)
  dump(OUT/'VALIDATION.json',dict(status='PASS',checks=len(checks),failed=0,tests_passed=5,saved_realizations=12000,MLE_records=36000,FIM_records=60000,new_Monte_Carlo_realizations=0,new_truth_cases=0,new_bearing_noise_draws=0,new_acoustic_propagation=0,new_depth_score=0))
  paths=[x for x in OUT.rglob('*') if x.is_file()]+[Path('results/R4_MASTER/R4_PLAN.md'),Path('results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'),Path('results/R4_MASTER/R4_SPEED_INFORMATION_PROGRESS.json')]
  dump(OUT/'SPEED_INFO_RESULTS_FREEZE.json',dict(bindings={x.as_posix():sha(x) for x in sorted(paths)}));print(json.dumps(dict(finalized='PASS',checks=len(checks))),flush=True)
 elif args.verify_seal:
  s=load(OUT/'SPEED_INFO_RESULTS_FREEZE.json')
  for n,h in s['bindings'].items():assert sha(n)==h,n
  print(json.dumps(dict(sealed='PASS',bindings=len(s['bindings']))),flush=True)
 else:ap.error('--finalize or --verify-seal required')
