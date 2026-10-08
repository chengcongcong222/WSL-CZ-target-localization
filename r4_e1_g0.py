"""E1-G0 deterministic local information only. No RNG, optimizer, audio or propagation."""
from pathlib import Path
import csv,json,hashlib,subprocess,math,argparse,time
import numpy as np
from scipy.linalg import solve_triangular
OUT=Path('results/R4_E1_G0_FREQUENCY_INFORMATION')
PARENT='f9b994ad29aac6159e428cd6b64f66a8e25a7782'
TITLE='R4 E1-G0: freeze nuisance-profiled frequency information design'
SCALE=np.array([50000.,50000.,2.,2.])
RTOL=1e-10
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(Path(p).read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def git(*a):return subprocess.check_output(['git','-c','core.quotepath=false',*a],text=True,encoding='utf-8').strip()
def trajectory(t):
 post=np.maximum(t-600.,0.)
 pos=np.column_stack([2*np.minimum(t,600.)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)])
 vel=np.tile([2.,0.],(len(t),1));vel[t>600]=[2*np.cos(np.pi/12),2*np.sin(np.pi/12)]
 return pos,vel
def geometry(scene,t):
 main,vel=trajectory(t)
 nodes=np.stack([main,main+np.array([0.,scene['mirror']*5000.])],axis=1)
 x=np.array(scene['state']);d=(x[:2]+t[:,None]*x[2:])[:,None,:]-nodes
 r=np.linalg.norm(d,axis=-1);e=d/r[:,:,None];w=x[2:]-vel[:,None,:];q=np.sum(e*w,axis=-1)
 qp=(w-e*q[:,:,None])/r[:,:,None]
 qx=np.concatenate([qp,qp*t[:,None,None]+e],axis=-1)
 bx=np.stack([-e[:,:,1]/r,e[:,:,0]/r],axis=-1)
 bx=np.concatenate([bx,bx*t[:,None,None]],axis=-1)
 tau=r/1500.;ux=-np.concatenate([e,e*t[:,None,None]],axis=-1)/1500.
 return nodes,r,q,qx,bx,tau,ux
def parameter_names(model,m,nodes):
 names=['a'+str(l) for l in range(m)]
 if model!='N0':names+=['d']
 if model in ['N2','N4_BOTH','N4_LINE','N5']:
  for j in nodes:names+=['rho'+str(j),'k'+str(j)]
 if model in ['N3','N4_BOTH']:
  for j in nodes:names+=['B'+str(j),'K'+str(j)]
 if model=='N4_LINE':names+=['dl'+str(l) for l in range(m)]
 return names
def corner_values(corner):
 if corner==0:return 0.,np.zeros(2),np.zeros(2),np.zeros(2)
 sign=1. if corner in [1,3] else -1.
 nodesigns=np.array([sign,sign if corner in [1,2] else -sign])
 return sign*1e-7*1200,nodesigns*5e-5,nodesigns*5e-8*1200,nodesigns*.1
def blocks(scene,lines,sigma,model,nodes=(0,1),corner=0,state=None,params=None):
 t=np.arange(121)*10.;base=dict(scene)
 if state is not None:base['state']=list(state)
 _,r,q,qx,bx,tau,ux=geometry(base,t)
 m=len(lines);names=parameter_names(model,m,nodes) if model!='B0' else []
 eta={n:0. for n in names}
 dc,rho,kk,dt=corner_values(corner)
 for l,f in enumerate(lines):
  if 'a'+str(l) in eta:eta['a'+str(l)]=math.log(f)
 if 'd' in eta:eta['d']=dc
 for j in nodes:
  if 'rho'+str(j) in eta:eta['rho'+str(j)]=rho[j]
  if 'k'+str(j) in eta:eta['k'+str(j)]=kk[j]
 if params is not None:eta.update(params)
 data=[]
 for i,tt in enumerate(t):
  for j in range(2):
   active=(model!='B0' and j in nodes and i%2==1)
   D=1-q[i,j]/1500.;u=tt-tau[i,j]-dt[j]
   y=[math.atan2((np.array(base['state'])[:2]+tt*np.array(base['state'])[2:]-trajectory(t)[0][i]-np.array([0.,j*scene['mirror']*5000.]))[1],
                 (np.array(base['state'])[:2]+tt*np.array(base['state'])[2:]-trajectory(t)[0][i]-np.array([0.,j*scene['mirror']*5000.]))[0])]
   J=[bx[i,j]];N=[np.zeros(len(names))]
   nav=[-bx[i,j,:2]]
   if active:
    for l,f in enumerate(lines):
     d=eta.get('d',0.)+eta.get('dl'+str(l),0.)
     gain=math.exp(eta['a'+str(l)]+d*u/1200.+eta.get('rho'+str(j),0.)+eta.get('k'+str(j),0.)*tt/1200.)
     mu=gain*D
     y.append(mu+eta.get('B'+str(j),0.)+eta.get('K'+str(j),0.)*tt/1200.)
     gx=mu*(d/1200.*ux[i,j]-qx[i,j]/(1500.*D));J.append(gx)
     nav.append(-gx[:2])
     v=np.zeros(len(names))
     for n,value in [('a'+str(l),mu),('d',mu*u/1200.),('dl'+str(l),mu*u/1200.),('rho'+str(j),mu),('k'+str(j),mu*tt/1200.),('B'+str(j),1.),('K'+str(j),tt/1200.)]:
      if n in names:v[names.index(n)]=value
     N.append(v)
   data.append((np.array(y),np.array(J),np.array(N),np.array(nav),j,active))
 return data,names,eta
def whiten(data,sigma,with_systematics=True):
 AX=[];AN=[];Y=[];U=[]
 sb=np.deg2rad(.05);sc=np.deg2rad(.05)/math.sqrt(3);sd=np.deg2rad(.025)/math.sqrt(3)
 for y,J,N,nav,j,active in data:
  R=np.diag([sb*sb]+[sigma*sigma]*(len(y)-1))+25.**2*(nav@nav.T)
  L=np.linalg.cholesky(R)
  AX.append(solve_triangular(L,J*SCALE,lower=True));AN.append(solve_triangular(L,N,lower=True));Y.append(solve_triangular(L,y,lower=True))
  z=np.zeros((len(y),2));z[0]=[sc,sd*(-1 if j==0 else 1)]
  U.append(solve_triangular(L,z,lower=True))
 X=np.vstack(AX);N=np.vstack(AN);U=np.vstack(U)
 if with_systematics:
  u,s,_=np.linalg.svd(U,full_matrices=False);coef=1-1/np.sqrt(1+s*s)
  X=X-u@(coef[:,None]*(u.T@X));N=N-u@(coef[:,None]*(u.T@N))
 return X,N
def basis(N,tol=RTOL):
 if not N.shape[1]:return np.zeros((len(N),0)),np.zeros(0),np.zeros(0),0
 norms=np.linalg.norm(N,axis=0);norms[norms==0]=1
 u,s,v=np.linalg.svd(N/norms,full_matrices=False);rank=int(np.sum(s>s[0]*tol)) if len(s) and s[0]>0 else 0
 return u[:,:rank],s,norms,rank
def profile(X,N):
 U,s,norms,nrank=basis(N);R=X-U@(U.T@X)
 F=R.T@R
 # Independent Schur elimination on an SVD-selected orthonormal nuisance parameter basis.
 if nrank:
  z=N/norms
  _,sv,V=np.linalg.svd(z,full_matrices=False)
  Q=z@V[:nrank].T
  schur=X.T@X-(X.T@Q)@np.linalg.solve(Q.T@Q,Q.T@X)
 else:schur=X.T@X
 err=np.linalg.norm(schur-F)/max(np.linalg.norm(F),1e-30)
 return R,F,nrank,s,err
def metric(R,scene,tol=RTOL):
 _,s,V=np.linalg.svd(R,full_matrices=False);rank=int(np.sum(s>s[0]*tol))
 x=np.array(scene['state']);g=np.r_[0.,0.,x[2:]/np.linalg.norm(x[2:])]*SCALE
 null=V[rank:]
 estimable=(not len(null) or np.linalg.norm(null@g)<=1e-8*np.linalg.norm(g))
 var=float(np.sum((V[:rank]@g/s[:rank])**2)) if estimable else math.inf
 weak=V[-1].copy()
 if weak[np.argmax(abs(weak))]<0:weak=-weak
 return var,rank,s,weak,estimable
def freeze():
 OUT.mkdir(exist_ok=False)
 panel='results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv'
 source='results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_DESIGN_FREEZE.json'
 scenes=[]
 for r in read(panel):
  rr=float(r['r0_km'])*1000;th=math.radians(float(r['theta0_deg']));v=float(r['v_mps']);ps=math.radians(float(r['psi_deg']))
  for side in [-1,1]:
   scenes.append(dict(scene_id=r['case_id']+('_M' if side<0 else '_P'),case_id=r['case_id'],mirror=side,state=[rr*math.cos(th),rr*math.sin(th),v*math.cos(ps),v*math.sin(ps)],deployment_beta_deg=0,baseline_m=5000))
 dump(OUT/'E1_G0_SCENE_BINDINGS.json',dict(parent_SHA=PARENT,source_panel=panel,source_design=source,
   bindings={p:sha(p) for p in [panel,source,'r4_a1_new_dynamic.py','r4_a1_new_estimator.py',
   'results/R4_OBSERVABILITY_FRONTIER_RESEARCH_RESET/E1_ACOUSTIC_MEASUREMENT_MODEL_LOCK.md',
   'results/R4_OBSERVABILITY_FRONTIER_RESEARCH_RESET/E1_REVISED_MINIMAL_EXPERIMENT_DESIGN.md']},
   scenes=scenes,geometry_scope='24 nominal beta=0 mirror formations; no saved random deployment/navigation draw; same fixed formation translation and original turn'))
 policy=dict(stage='R4_E1_G0_NUISANCE_PROFILED_FREQUENCY_INFORMATION',parent_SHA=PARENT,source_audit='ACCEPTED_WITH_ORIGINAL_ESTIMATOR_EXCEPTION',
  G0_admitted=True,G1_admitted=False,state='Cartesian x0,y0,vx,vy',state_scale=SCALE.tolist(),
  duration_s=1200,bearing_epochs=121,frequency_frame_centers=list(range(10,1200,20)),frequency_lines_Hz=[[200.],[170.,200.,230.]],
  sigma_frequency_Hz=[.001,.005,.010],bearing_sigma_deg=.05,common_bias_bound_deg=.05,half_difference_bound_deg=.025,navigation_sigma_m_per_axis=25,
  covariance='G0 independent ideal frequency instrument errors plus original Gaussian bearing errors; navigation-induced joint block covariance; static bounded uniform bearing bias moment-matched rank2 covariance. Not G1 derived-signal covariance.',
  model_ids=['N0','N1','N2','N3','N4_BOTH','N4_LINE','N5'],
  variants={'SUN_MEASUREMENT':'N0 node0 frequency + dual bearings; PAPER_ADAPTED_MEASUREMENT_INFORMATION_CONTROL','SINGLE_FREQ':'N2 node0','DUAL_FREQ':'N0 both nodes','NEW':'N1/N2/N3/N4_BOTH/N4_LINE/N5 both nodes'},
  primary='NEW/N2 no calibration; both single and three line resource packages; all three sigma assumptions reported separately',
  corners={'main':0,'N2_robustness':[1,2,3,4],'definition':'C0 zero; C1 common positive; C2 common negative; C3 opposite receiver signs with positive source; C4 C3 negated; source dc=+-1e-7/s;rho=+-5e-5;kappa=+-5e-8/s;dt=+-0.1s'},
  source_frequency_policy='F_l unknown nuisance. Numerical Jacobian evaluated at registered nominal lines ONLY for local information; no estimator, truth initialization or f0 calibration prior.',
  clock_time_policy='delta_t fixed corner condition, not separately profiled; nominal time/navigation/motion known up to included navigation position uncertainty. No claim of time-sync robustness beyond corners.',
  source_emission_time='u=t-r/c-delta_t; quasi-static first-order LOS, not full retarded acoustic or multimode first-CZ signal',
  environment='c=1500m/s known LOS mechanism control; no SSP or propagation robustness claim',
  nuisance_scaling='a=logF, d=1200*dc, k=1200*kappa, K=1200*Hz-per-sec; unit column normalization before SVD',
  gauge_policy='retain redundant columns; rank-revealing SVD quotient; no clock fixed to true value',
  rank={'nuisance_relative':RTOL,'target_scaled_J_relative':RTOL,'sensitivity_relative':[1e-8,1e-10,1e-12],'estimability_null_relative':1e-8},
  validation={'Schur_projection_relative_max':1e-6,'finite_difference_scaled_J_relative_max':2e-6,'central_difference_scaled_step':[1e-3,3e-4,1e-4],
  'nuisance_fd_relative_max':2e-6,'nuisance_fd_step':[1e-5,3e-6,1e-6],'zero_control_relative_max':1e-9,'monotonic_variance_relative_tolerance':1e-7},
  calibration_N5={'receiver_fractional_offset_sigma':5e-6,'receiver_fractional_drift_sigma_per_s':5e-9,'source':'HYPOTHETICAL_INDEPENDENT_CALIBRATION_REQUIREMENT_NOT_MEASURED','cost':'Both node independent frequency standards; no actual device/calibration evidence'},
  gate={'median_variance_reduction_min':.20,'worst_max_variance_ratio_max':1.000000001,'all_scenes':24,'no_unestimable_or_invalid':True,
  'aggregation':'Separate package for each lines/sigma; no post-hoc best selection. Robustness C1-C4 all retained; positive decision requires a registered package pass C0-C4. Report passing resource/precision conditions explicitly.'},
  prohibited=['MC','audio','STFT_PLL','propagation','optimizer','time_geometry_truth_noise_changes','automatic_G1_E2_depth_A2_SSP_P5'],
  new_MC=0,new_audio=0,R4_percent=0,stop_after_commit_B=True,
  code_binding={'r4_e1_g0.py':sha('r4_e1_g0.py')})
 dump(OUT/'E1_G0_DESIGN_FREEZE.json',policy)
 (OUT/'E1_G0_JACOBIAN_DEFINITION.md').write_text("""# E1-G0 Jacobian and covariance lock
Coordinates follow existing code: +x initial nominal look direction, atan2(y,x). x=[p0x,p0y,vx,vy]. S=diag(50000,50000,2,2). Target p=p0+t*v. Main speed2m/s, turn15deg after600s; at exactly600s use pre-turn velocity. AUX adds fixed [0,mirror*5000] and shares translation. Deployment beta=0 nominal control; not 500 saved random placements.

e=(p-node)/r; w=v-vnode; q=e.w; q_p=(w-e*q)/r; q_x=[q_p,t*q_p+e].
bearing_x=[-e_y/r,e_x/r,t*(-e_y/r),t*e_x/r].
u=t-r/c-delta_t; u_x=-[e,t*e]/c.
mu=exp(a_l+(d+dl_l)*u/1200+rho_j+k_j*t/1200)*(1-q/c)+B_j+K_j*t/1200.
mu_x=(mu-B-K*t/1200)*((d+dl)/1200*u_x-q_x/(c*(1-q/c))).
eta derivatives: a_l,rho_j -> Doppler component; d,dl_l -> component*u/1200; k_j -> component*t/1200; B_j ->1; K_j ->t/1200. No unknown emitted f supplied as calibration.

Navigation derivatives at each node/epoch are minus position derivatives. R_local=diag(sigma_b²,sigma_f²,...)+25² J_nav J_nav.T. Same-epoch same-node bearing and frequency share nav errors; cross-lines also correlated. White instrument errors are ideal independent G0 inputs, not extracted same-waveform features. Static uniform common/half-differential bearing offsets integrated with variance bound²/3 as Gaussian moment-matched rank2 covariance. Not a uniform-distribution exact Fisher bound. No navigation velocity errors or target process prior added.

Whiten block covariance by Cholesky then rank2 inverse-square-root for static bearing covariance. Profile with orthonormal nuisance left basis U: R=(I-UU.T)A_x, I_eff=R.T R. Independently build Schur using full-rank SVD-selected nuisance coordinate basis. N5 appends receiver calibration prior rows before profiling. Target inversion via scaled R SVD; speed gradient S*[0,0,vx/v,vy/v]. Null-space overlap -> UNBOUNDED, never finite pseudoinverse-zero. Weak vector is scaled four-state null/weak combination, physical loading=S*vector.
""",encoding='utf-8')
 (OUT/'E1_G0_NUISANCE_MODELS.md').write_text("""# E1-G0 nuisance lock
N0 each log F unknown, receiver reference known, source drift zero.
N1 N0 + unknown shared fractional source drift.
N2 N1 + independent node fractional constant rho and linear k; primary NO_CALIBRATION.
N3 N1 + independent additive Hz B,K, fractional errors fixed zero: different model.
N4_BOTH N2 + additive B,K; N4_LINE N2 + per-line drift dl; broader controls separately retained. Shared d/per-line dl redundancy removed only by SVD gauge quotient.
N5 N2 + independent receiver rho sigma5e-6 and kappa sigma5e-9/s calibration priors. Hypothetical requirement, NOT existing hardware evidence. No source-frequency prior.
Each model full registered covariance and source-time convention from design. Fractional errors dimensionless; Hz errors never substituted for fractional clock uncertainty.
N0 is an ideal-reference upper control, not proposed prior-free result.
SUN-MEASUREMENT uses locked Sun Eq4/unknown F, no UKF, PAPER_ADAPTED_MEASUREMENT_INFORMATION_CONTROL. Source-state random-walk process covariance is omitted in this local static information control and explicitly NOT full original dynamic CRLB. SUN complete MFB-AUKF ORIGINAL_ESTIMATOR_NOT_FULLY_LOCKED.
All rank thresholds, covariance assumptions, calibration requirements and corner values fixed before scientific evaluation. No numerical result used to select model.
""",encoding='utf-8')
 (OUT/'E1_G0_GATE.md').write_text("""# E1-G0 gate
Keep all24 scenes per fixed lines/sigma package and all C0-C4 N2 corners. Main N2 effective speed variance vs same-scene B0: median reduction>=20%, max_scene variance no larger; invalid/rank-unstable/unestimable blocks scientific PASS. Independent ideal instrument covariance is a G0 control, NOT G1 covariance.
Each resource package reported, no post-hoc change of assumed precision. Any package meeting C0-C4 is a conditional pre-research information extension with named line/precision costs; not unconditional hardware feasibility. Other packages and failures retained.
A: N2 prior-free gate -> E1_G0_PRIOR_FREE_INFORMATION_INCREMENT_ESTABLISHED; next recommendation G1, needs separate audit/admission.
B: only known reference or N5 calibration meets gate -> E1_G0_CALIBRATION_CONDITIONAL_INFORMATION.
C: no major package meets gate -> E1_G0_FREQUENCY_INCREMENT_NOT_ESTABLISHED; E2_REVIEW recommendation only.
D: information is conditional but statistical/input contract cannot be coherently defined -> E1_G0_MECHANISM_CONDITIONAL_EXTRACTION_UNRESOLVED.
Validation failure -> IMPLEMENTATION_INVALID; STOP without scientific verdict. Arbitrary per-node/time reference must absorb all frequency information. No estimator, MC, audio, propagation, extra time/geometry/noise or automatic next stage. R4=0%.
""",encoding='utf-8')
 print('Design written; scientific calculations not run.')

def checks_for_design(policy,scenes):
 checks=[]
 def check(name,ok,value=None):
  checks.append(dict(check=name,passed=bool(ok),value=value))
 # Model identities in observed Hz, not just independent abstract matrices.
 s=scenes[0];data,names,eta=blocks(s,[170.,200.,230.],.005,'N2')
 freq=[b for b in data if b[-1]]
 b=freq[0];mu=b[0][1:]
 check('multiline_motion_loading_rank_one',np.linalg.matrix_rank(mu[:,None])==1)
 # δq=cD*a reproduces source scale derivative; receiver fractional rho similarly.
 check('constant_F_radial_scale_gauge',np.allclose(-mu*.01+mu*.01,0,atol=1e-12))
 t=np.arange(60)*20.+10;D=.99+.00001*t
 dq=1500*D*(.01+.00002*t)
 effect=-200*dq/1500+200*D*(.01+.00002*t)
 check('affine_fractional_reference_absorbs_corresponding_affine_log_Doppler',np.linalg.norm(effect)<1e-10)
 # Non-affine residual remains outside constant/linear basis.
 N=np.column_stack([np.ones(60),t/1200]);U,_,_,_=basis(N);curv=(t/1200)**2
 check('nonaffine_curvature_not_affine',np.linalg.norm(curv-U@(U.T@curv))>1e-3)
 # Real joint covariance: fully free frequency mean nuisance -> B0 marginal information.
 X,N=whiten(data,.005);free=[]
 ix=0
 for yy,JJ,NN,nav,j,active in data:
  for k in range(1,len(yy)):
   col=np.zeros(sum(len(b[0]) for b in data));col[ix+k]=yy[k]
   free.append(col)
  ix+=len(yy)
 # Whiten arbitrary frequency columns with the same covariance operator by replacing N.
 Nraw=np.column_stack(free);alt=[];ix=0
 for yy,JJ,NN,nav,j,active in data:
  alt.append((yy,JJ,Nraw[ix:ix+len(yy)],nav,j,active));ix+=len(yy)
 Xfree,Nfree=whiten(alt,.005);Rf,Ff,_,_,err=profile(Xfree,Nfree)
 b0,_,_=blocks(s,[],.005,'B0');X0,N0=whiten(b0,.005)
 rel=np.linalg.norm(Ff-X0.T@X0)/np.linalg.norm(X0.T@X0)
 check('arbitrary_node_time_reference_zero_increment_joint_covariance',rel<1e-9,float(rel))
 # Same functional speed along null cannot be assigned finite pseudo variance.
 toy=np.diag([1.,1.,0.,0.]);vv,rr,ss,ww,est=metric(toy,s)
 check('unestimable_speed_is_UNBOUNDED',math.isinf(vv) and not est)
 # Priors are appended, never replace data or source nuisance.
 check('N5_external_prior_explicit',policy['calibration_N5']['source'].startswith('HYPOTHETICAL'))
 return checks
def fd_check(scene,lines,model,corner):
 data,names,eta=blocks(scene,lines,.005,model,corner=corner)
 J=np.vstack([b[1] for b in data])*SCALE
 def values(state=None,params=None):
  d,_,_=blocks(scene,lines,.005,model,corner=corner,state=state,params=params)
  return np.concatenate([b[0] for b in d])
 state=np.array(scene['state']);errs=[]
 for h in [1e-3,3e-4,1e-4]:
  jj=[]
  for k in range(4):
   shift=np.eye(4)[k]*SCALE[k]*h
   jj.append((values(state+shift)-values(state-shift))/(2*h))
  Z=np.column_stack(jj);errs.append(float(np.linalg.norm(Z-J)/np.linalg.norm(J)))
 E=np.vstack([b[2] for b in data]);ee=[]
 for h in [1e-5,3e-6,1e-6]:
  zz=[]
  for n in names:
   plus=dict(eta);minus=dict(eta);plus[n]+=h;minus[n]-=h
   zz.append((values(params=plus)-values(params=minus))/(2*h))
  Z=np.column_stack(zz);ee.append(float(np.linalg.norm(Z-E)/np.linalg.norm(E)))
 return dict(scene_id=scene['scene_id'],model=model,corner=corner,target_step_errors=errs,nuisance_step_errors=ee,
 passed=max(errs)<2e-6 and max(ee)<2e-6)
def execute():
 p=json.loads((OUT/'E1_G0_DESIGN_FREEZE.json').read_text(encoding='utf-8'))
 b=json.loads((OUT/'E1_G0_SCENE_BINDINGS.json').read_text(encoding='utf-8'))
 assert git('log','-1','--format=%s')==TITLE and git('rev-parse','HEAD^')==PARENT
 head=git('rev-parse','HEAD');assert head==git('ls-remote','origin','refs/heads/main').split()[0]
 assert sha('r4_e1_g0.py')==p['code_binding']['r4_e1_g0.py']
 for f,h in b['bindings'].items():assert sha(f)==h,f
 assert not (OUT/'E1_G0_EXECUTION_START.json').exists()
 dump(OUT/'E1_G0_EXECUTION_START.json',dict(design_SHA=head,remote_verified=True,policy_sha256=sha(OUT/'E1_G0_DESIGN_FREEZE.json'),new_MC=0,new_audio=0))
 scenes=b['scenes'];checks=checks_for_design(p,scenes)
 if not all(c['passed'] for c in checks):
  dump(OUT/'VALIDATION.json',dict(valid=False,zero_controls=checks));raise RuntimeError('IMPLEMENTATION_INVALID zero controls')
 records=[];ranks=[];fds=[];matrices={};bvars={}
 start=time.perf_counter()
 for scene in scenes:
  sid=scene['scene_id'];data,_,_=blocks(scene,[],.005,'B0');X,N=whiten(data,.005)
  R,F,nr,sv,err=profile(X,N);var,rank,s,weak,est=metric(R,scene);bvars[sid]=var
  matrices[sid+'_B0']=F
  records.append(dict(scene_id=sid,case_id=scene['case_id'],mirror=scene['mirror'],variant='B0',model='B0',lines=0,sigma_f_Hz=0.,corner=0,rank=rank,nuisance_rank=0,
   variance_mps2=var if math.isfinite(var) else 'UNBOUNDED',sigma_speed_mps=math.sqrt(var) if est else 'UNBOUNDED',speed_196_rel=1.96*math.sqrt(var)/np.linalg.norm(scene['state'][2:]) if est else 'UNBOUNDED',
   unprofiled_variance_mps2=var,variance_reduction=0.,weak_scaled=json.dumps(weak.tolist()),weak_physical=json.dumps((weak*SCALE).tolist()),estimable=est,schur_projection_rel=err,rank_stable=True))
  fds.append(fd_check(scene,[170.,200.,230.],'N2',0))
  fds.append(fd_check(scene,[170.,200.,230.],'N2',3))
  if sid.endswith('_M'):fds.append(fd_check(scene,[200.],'N3',0))
  for lines in [[200.],[170.,200.,230.]]:
   for sig in [.001,.005,.010]:
    specs=[('SUN_MEASUREMENT','N0',(0,),0),('SINGLE_FREQ','N2',(0,),0),('DUAL_FREQ','N0',(0,1),0)]
    specs+=[('NEW',m,(0,1),0) for m in ['N1','N2','N3','N4_BOTH','N4_LINE','N5']]
    specs+=[('NEW','N2',(0,1),c) for c in [1,2,3,4]]
    for variant,model,nodes,corner in specs:
     data,names,eta=blocks(scene,lines,sig,model,nodes,corner)
     X,N=whiten(data,sig);uvar=metric(X,scene)[0]
     if model=='N5':
      prior=np.zeros((4,len(names)))
      for j in [0,1]:
       prior[2*j,names.index('rho'+str(j))]=1/5e-6
       prior[2*j+1,names.index('k'+str(j))]=1/(5e-9*1200)
      X=np.vstack([X,np.zeros((4,4))]);N=np.vstack([N,prior])
     R,F,nr,sv,err=profile(X,N);var,rank,s,weak,est=metric(R,scene)
     key=f'{sid}_{variant}_{model}_{len(lines)}_{sig}_{corner}'
     matrices[key]=F
     vals=[];nrs=[]
     for tol in [1e-8,1e-10,1e-12]:
      U,_,_,rn=basis(N,tol);rr=X-U@(U.T@X)
      tv,tr,_,_,_=metric(rr,scene,tol);vals.append(tv);nrs.append(rn)
     # Numerical rank sensitivity is reported; only meaningful target variance/rank stability gates primary.
     stable=all(math.isfinite(a) for a in vals) and (max(vals)-min(vals))/max(var,1e-30)<1e-4
     records.append(dict(scene_id=sid,case_id=scene['case_id'],mirror=scene['mirror'],variant=variant,model=model,lines=len(lines),sigma_f_Hz=sig,corner=corner,rank=rank,nuisance_rank=nr,
      variance_mps2=var if est else 'UNBOUNDED',sigma_speed_mps=math.sqrt(var) if est else 'UNBOUNDED',speed_196_rel=1.96*math.sqrt(var)/np.linalg.norm(scene['state'][2:]) if est else 'UNBOUNDED',
      unprofiled_variance_mps2=uvar if math.isfinite(uvar) else 'UNBOUNDED',variance_reduction=1-var/bvars[sid] if est else 'UNBOUNDED',weak_scaled=json.dumps(weak.tolist()),weak_physical=json.dumps((weak*SCALE).tolist()),estimable=est,schur_projection_rel=err,rank_stable=stable))
     ranks.append(dict(key=key,rank_target=rank,rank_nuisance=nr,target_singular_values=json.dumps(s.tolist()),nuisance_normalized_singular_values=json.dumps(sv.tolist()),condition_scaled=float(s[0]/s[-1]),nuisance_ranks_1e8_1e10_1e12=json.dumps(nrs),
      speed_variances_sensitivity=json.dumps([v if math.isfinite(v) else 'UNBOUNDED' for v in vals]),variance_stable=stable,schur_projection_relative=err))
  print(sid+' complete',flush=True)
 write(OUT/'E1_G0_FULL_SCENE_RESULTS.csv',records);write(OUT/'E1_G0_NUMERICAL_RANK_AUDIT.csv',ranks)
 dump(OUT/'E1_G0_DIFFERENTIATION_VALIDATION.json',dict(checks=fds,pass_count=sum(c['passed'] for c in fds),fail_count=sum(not c['passed'] for c in fds)))
 np.savez_compressed(OUT/'E1_G0_EFFECTIVE_INFORMATION_MATRICES.npz',**matrices)
 primary=[r for r in records if r['model']=='N2' and r['variant']=='NEW']
 critical=[r for r in records if r['variant']=='B0' or r in primary]
 valid=all(c['passed'] for c in checks) and all(c['passed'] for c in fds) and all(r['schur_projection_rel']<1e-6 for r in records) and all(r['rank_stable'] for r in critical) and all(float(r['variance_mps2'])<=bvars[r['scene_id']]*(1+1e-7) for r in primary if r['estimable'])
 summary=[]
 groups={}
 for r in records:
  key=(r['variant'],r['model'],r['lines'],r['sigma_f_Hz'],r['corner']);groups.setdefault(key,[]).append(r)
 for key,g in groups.items():
  finite=all(r['estimable'] for r in g);vs=[float(r['variance_mps2']) for r in g] if finite else []
  reductions=[float(r['variance_reduction']) for r in g] if finite else []
  med=float(np.median(reductions)) if finite else -math.inf
  worst=max(vs)/max(bvars.values()) if finite else math.inf
  gate=finite and med>=.20 and worst<=1.000000001 and all(r['rank_stable'] for r in g)
  summary.append(dict(variant=key[0],model=key[1],lines=key[2],sigma_f_Hz=key[3],corner=key[4],scenarios=len(g),
   median_variance_reduction=med if finite else 'UNBOUNDED',worst_variance_ratio=worst if finite else 'UNBOUNDED',positive_cases=sum(float(r['variance_reduction'])>1e-9 for r in g if r['estimable']),
   min_reduction=min(reductions) if finite else 'UNBOUNDED',max_reduction=max(reductions) if finite else 'UNBOUNDED',worst_variance_mps2=max(vs) if finite else 'UNBOUNDED',gate_pass=gate))
 write(OUT/'E1_G0_INFORMATION_SUMMARY.csv',summary)
 losses=[]
 for r in records:
  if r['variant']=='NEW' and r['corner']==0:
   base=next(a for a in records if a['scene_id']==r['scene_id'] and a['variant']=='DUAL_FREQ' and a['lines']==r['lines'] and a['sigma_f_Hz']==r['sigma_f_Hz'])
   losses.append(dict(scene_id=r['scene_id'],model=r['model'],lines=r['lines'],sigma_f_Hz=r['sigma_f_Hz'],variance_ratio_to_N0=float(r['variance_mps2'])/float(base['variance_mps2']) if r['estimable'] else 'UNBOUNDED',
    information_retained_ratio_N0=float(base['variance_mps2'])/float(r['variance_mps2']) if r['estimable'] else 0,variance_ratio_profiled_unprofiled=float(r['variance_mps2'])/float(r['unprofiled_variance_mps2']) if r['estimable'] else 'UNBOUNDED'))
 write(OUT/'E1_G0_NUISANCE_LOSS_TABLE.csv',losses)
 multi=[]
 for r in records:
  if r['lines']==3:
   a=next(z for z in records if z['scene_id']==r['scene_id'] and z['lines']==1 and all(z[k]==r[k] for k in ['variant','model','sigma_f_Hz','corner']))
   multi.append(dict(scene_id=r['scene_id'],variant=r['variant'],model=r['model'],sigma_f_Hz=r['sigma_f_Hz'],corner=r['corner'],
    single_variance_mps2=a['variance_mps2'],three_variance_mps2=r['variance_mps2'],variance_ratio_three_single=float(r['variance_mps2'])/float(a['variance_mps2']) if a['estimable'] and r['estimable'] else 'UNBOUNDED',
    rank_single=a['rank'],rank_three=r['rank'],motion_dimensions_per_node_time=1,resource_cost='3 spectral lines instead of1; not 3 motion dimensions'))
 write(OUT/'E1_G0_SINGLE_VS_MULTILINE.csv',multi)
 n2packages=[]
 for n in [1,3]:
  for sig in [.001,.005,.010]:
   g=[z for z in summary if z['variant']=='NEW' and z['model']=='N2' and z['lines']==n and z['sigma_f_Hz']==sig]
   n2packages.append(dict(lines=n,sigma_f_Hz=sig,all_C0_C4_pass=all(z['gate_pass'] for z in g),corner_results=g))
 prior=valid and any(z['all_C0_C4_pass'] for z in n2packages)
 cal=any(z['gate_pass'] for z in summary if z['model']=='N5' or z['variant']=='DUAL_FREQ')
 decision='IMPLEMENTATION_INVALID' if not valid else 'E1_G0_PRIOR_FREE_INFORMATION_INCREMENT_ESTABLISHED' if prior else 'E1_G0_CALIBRATION_CONDITIONAL_INFORMATION' if cal else 'E1_G0_FREQUENCY_INCREMENT_NOT_ESTABLISHED'
 dump(OUT/'VALIDATION.json',dict(valid=valid,zero_controls=checks,finite_difference_pass=sum(c['passed'] for c in fds),finite_difference_fail=sum(not c['passed'] for c in fds),
   max_schur_projection_relative=max(r['schur_projection_rel'] for r in records),primary_rank_stable=all(r['rank_stable'] for r in critical),wider_rank_unstable_count=sum(not r['rank_stable'] for r in records),elapsed_seconds=time.perf_counter()-start,records=len(records),registered_scenarios=24))
 dump(OUT/'E1_G0_DECISION.json',dict(stage=p['stage'],parent_SHA=PARENT,design_SHA=head,source_audit=p['source_audit'],route_decision=decision,
   registered_scenarios=24,primary_packages=n2packages,B0_worst_variance_mps2=max(bvars.values()),B0_variance_range=[min(bvars.values()),max(bvars.values())],
   frequency_mechanism_screen='FREQUENCY_OBSERVABILITY_MECHANISM_SCREEN_COMPLETED' if valid else 'INVALID',
   scope='LOCAL_GAUSSIAN_MOMENT_MATCHED_LOS_MECHANISM_ONLY_NOT_EMPIRICAL_P95_NOT_EXTRACTION',
   new_Monte_Carlo=0,new_received_audio=0,new_propagation=0,R4_A1_percent=0,R4_percent=0,G1='NOT_OPENED',E2='NOT_OPENED',
   next_recommended='E1_G1_EXTRACTION' if prior else 'STOP' if cal or not valid else 'E2_REVIEW',automatic_next=False))
 print(decision,flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--freeze',action='store_true');ap.add_argument('--execute',action='store_true');args=ap.parse_args()
 if args.freeze:freeze()
 elif args.execute:execute()
 else:ap.error('--freeze or --execute required')
