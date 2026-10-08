"""Independent reconstruction and reporting; reads frozen run artifacts, never calls forward."""
from pathlib import Path
import csv,json,hashlib,subprocess,sys
import numpy as np
from scipy import linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path('results/R4_E2_G0_HLA_DIFFERENCE_INFORMATION')
SCALE=np.array([50000.,50000.,2.,2.,200.])
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(p,raw=False):
 b=Path(p).read_bytes();return hashlib.sha256(b if raw else b.replace(b'\r\n',b'\n')).hexdigest()
def git(*a):return subprocess.check_output(['git',*a],text=True,encoding='utf-8').strip()
def check(n,ok,detail=None):checks.append(dict(check=n,PASS=bool(ok),detail=detail))
def independent_whitening(raw,pairs,selected):
 ix={float(f):i for i,f in enumerate(raw)}
 pairs=[p for p in pairs if selected=='JOINT' or p['delta_hz']==selected]
 bp=np.zeros((len(pairs),len(raw)));bm=bp.copy()
 for i,p in enumerate(pairs):
  bp[i,ix[p['low_hz']]]=1;bp[i,ix[p['high_hz']]]=1
  bm[i,ix[p['low_hz']]]=-1;bm[i,ix[p['high_hz']]]=1
 def q(mat):
  e,v=np.linalg.eigh(mat.T@mat);return v[:,e>1e-10*max(e[-1],1)].T
 return q(bp),q(bm)
def helmert():
 h=np.zeros((7,8))
 for k in range(1,8):
  h[k-1,:k]=1/np.sqrt(k*(k+1));h[k-1,k]=-k/np.sqrt(k*(k+1))
 return h
def cold_design(j,qp,qm,c):
 # A distinct spatial basis and eigenvalue construction, followed by least squares nuisance removal.
 h=helmert()
 a=np.einsum('sm,tnmfp->tnsfp',h,j.real);b=np.einsum('sm,tnmfp->tnsfp',h,j.imag)
 a=np.einsum('qf,tnsfp->tnsqp',qp,a);b=np.einsum('qf,tnsfp->tnsqp',qm,b)
 if c=='C1':
  v=np.tile(qp@np.ones(j.shape[3]),j.shape[0])
  if np.linalg.norm(v)>1e-9:
   for node in range(j.shape[1]):
    for elem in range(7):
     aa=a[:,node,elem].reshape(-1,5)
     fit=linalg.lstsq(v[:,None],aa,lapack_driver='gelsy')[0]
     a[:,node,elem]=(aa-v[:,None]@fit).reshape(a.shape[0],a.shape[3],5)
 return np.r_[a.reshape(-1,5),b.reshape(-1,5)]
def effective(j,col,exact=False):
 if exact:j=j[:,:4]
 other=[i for i in range(j.shape[1]) if i!=col]
 fit=linalg.lstsq(j[:,other],j[:,col],cond=1e-10,lapack_driver='gelsd')[0]
 r=j[:,col]-j[:,other]@fit
 return float(r@r)

def cold_forward_reconstruction():
 import struct
 cache={}
 def mod(f,mesh):
  key=(f,mesh)
  if key in cache:return cache[key]
  path=OUT/'modal'/f'f{int(round(f*2)):04d}_s+0_n{mesh}.mod'
  buf=path.read_bytes();rec=4*struct.unpack('<i',buf[:4])[0]
  ndepth=struct.unpack('<i',buf[92:96])[0]
  z=np.frombuffer(buf[4*rec:4*rec+4*ndepth],dtype='<f4').astype(float)
  count=struct.unpack('<i',buf[5*rec:5*rec+4])[0]
  phi=np.empty((count,ndepth),complex)
  for n in range(count):
   phi[n]=np.frombuffer(buf[(7+n)*rec:(7+n)*rec+8*ndepth],dtype='<c8')
  k=np.frombuffer(buf[(7+count)*rec:(7+count)*rec+8*count],dtype='<c8').astype(complex)
  cache[key]=(z,phi,k);return cache[key]
 def p(f,mesh,depth,r):
  z,phi,k=mod(f,mesh);a=int(np.flatnonzero(z==depth)[0]);b=int(np.flatnonzero(z==200)[0])
  rr=np.asarray(r)[None,:]
  # Direct summation, separate scalar complex normalization via least squares.
  terms=(phi[:,a]*phi[:,b])[:,None]*np.sqrt(2*np.pi/(k[:,None]*rr))*np.exp(-1j*k[:,None]*rr-1j*np.pi/4)
  return np.sum(terms,axis=0)
 record=rows(OUT/'E2_FORWARD_CONVERGENCE.csv');diffs=[];cold=[]
 for row in record:
  f=float(row['frequency_hz']);z=float(row['source_depth_m'])
  if row['kind']=='mesh40001_vs20001':
   r=np.array([50000.,50007.,55000.,55007.,60000.,60007.])
   a=p(f,20001,z,r);b=p(f,40001,z,r)
  else:
   buf=(OUT/'modal'/f'field_f{int(f)}_z{int(z)}.shd').read_bytes()
   rec=4*struct.unpack('<i',buf[:4])[0]
   r=np.frombuffer(buf[9*rec:9*rec+12],dtype='<f4').astype(float)
   b=np.frombuffer(buf[10*rec:10*rec+24],dtype='<c8').astype(complex)
   a=p(f,20001,z,r)
  coeff=linalg.lstsq(a[:,None],b,lapack_driver='gelsd')[0][0]
  error=float(np.linalg.norm(b-coeff*a)/np.linalg.norm(b))
  diff=abs(error-float(row['error']));diffs.append(diff)
  cold.append(dict(kind=row['kind'],frequency_hz=f,source_depth_m=z,independent_error=error,record_error=float(row['error']),difference=diff,admission_PASS=error<.002))
 with (OUT/'E2_INDEPENDENT_FORWARD_RECONSTRUCTION.csv').open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(cold[0]));w.writeheader();w.writerows(cold)
 check('independent direct-sum/GELSD forward reconstruction',max(diffs)<1e-11,dict(checks=len(diffs),max_absolute_difference=max(diffs)))
 check('independent physical-gate verdicts',all((row['PASS']=='True')==r['admission_PASS'] for row,r in zip(record,cold)))
 mesh=[r for r in cold if r['kind']=='mesh40001_vs20001'];field=[r for r in cold if r['kind']=='independent_FIELD_aperture']
 status=dict(mesh_checks=len(mesh),mesh_FAIL=sum(not r['admission_PASS'] for r in mesh),FIELD_checks=len(field),FIELD_FAIL=sum(not r['admission_PASS'] for r in field),worst=max(cold,key=lambda r:r['independent_error']),maximum_record_difference=max(diffs))
 dump(OUT/'E2_INDEPENDENT_PHYSICS_SUMMARY.json',status)

def audit():
 global checks
 checks=[]
 design=read(OUT/'E2_G0_DESIGN_FREEZE.json');dec=read(OUT/'E2_G0_DECISION.json')
 started=read(OUT/'EXECUTION_STARTED.json');freq=read(OUT/'E2_REQUIRED_FREQUENCIES.json')
 scenes=read(OUT/'E2_SCENE_FREEZE.json')['scenes'];ids={s['scene_id'] for s in scenes}
 design_sha=started['design_SHA']
 check('24 frozen scenes',len(ids)==24)
 check('design commit descends from accepted parent',git('rev-parse',design_sha+'^')==design['parent_SHA'])
 check('two-stage title order',git('log','-1','--pretty=%s',design_sha)==design['commit_A_title'])
 for p,h in design['bindings'].items():
  raw=p.endswith(('.pdf','.exe','.mod')) or p.endswith('/.gitattributes')
  check('frozen binding '+p,sha(p,raw)==h)
 for p,h in design['source_raw_sha256'].items():check('source original bytes '+p,sha(p,True)==h)
 plane=read(OUT/'PLANE_CONTROL_RECORDS.json')
 for r in plane:check('null '+r['control'],r['PASS'],r['error'])
 check('no MC/audio',dec['new_Monte_Carlo']==0 and dec['new_received_audio']==0)
 check('R4 remains zero',dec['R4_percent']==0)
 check('science code identical to design',git('show',design_sha+':r4_e2_g0.py').replace('\r\n','\n')==Path('r4_e2_g0.py').read_text(encoding='utf-8').replace('\r\n','\n').strip())
 for filename in ('E2_FULL_SCENE_INFORMATION.csv','E2_NUISANCE_PROFILE_RESULTS.csv','E2_CANDIDATE_SEPARATION.csv','E2_CALIBRATION_SENSITIVITY.csv'):
  rr=rows(OUT/filename);check('all scenes '+filename,{r['scene_id'] for r in rr}==ids,len(rr))
 if (OUT/'E2_MODAL_GENERATION.csv').exists():
  rr=rows(OUT/'E2_MODAL_GENERATION.csv')
  check('modal hard call cap',len(rr)<=design['budgets']['modal_solver_calls_max'],len(rr))
  check('all solver timeout records within cap',all(float(r['elapsed_s'])<90.5 for r in rr))
  for r in rr:
   stem=f"f{int(round(float(r['frequency_hz'])*2)):04d}_s{int(r['ssp_offset_mps']):+d}_n{int(r['mesh'])}"
   check('modal hash '+stem,sha(OUT/'modal'/(stem+'.mod'),True)==r['mod_sha256'])
 if (OUT/'PHYSICS_ADMISSION.json').exists():
  admission=read(OUT/'PHYSICS_ADMISSION.json')
  rr=rows(OUT/'E2_FORWARD_CONVERGENCE.csv')
  check('physical admission aggregation',admission['PASS']==all(r['PASS']=='True' for r in rr))
  cold_forward_reconstruction()
 if (OUT/'INFORMATION_SUMMARY.json').exists():
  cov=read(OUT/'COVARIANCE_CONTROL.json');check('shared covariance reconstruction',cov['PASS'],cov)
  local=rows(OUT/'E2_LOCAL_NUMERICAL_CHECKS.csv');summary=read(OUT/'INFORMATION_SUMMARY.json')
  check('local admission aggregation',summary['checks_pass']==all(r['PASS']=='True' for r in local))
  ir=rows(OUT/'E2_FULL_SCENE_INFORMATION.csv')
  matrices={sel:independent_whitening(freq['raw_frequencies_hz'],freq['frequency_pairs'],sel) for sel in (2,5,10,'JOINT')}
  errors=[];count=0
  for sc in scenes:
   data=np.load(OUT/f"INPUT_{sc['scene_id']}.npz");j=data['relative_pressure_J']
   for nodes in (1,2):
    for sel in (2,5,10,'JOINT'):
     qp,qm=matrices[sel]
     for c in ('C0','C1'):
      z=cold_design(j[:,:nodes],qp,qm,c)
      ref=[r for r in ir if r['scene_id']==sc['scene_id'] and r['representation']==('P1' if nodes==1 else 'P2') and r['condition']==c and r['delta_hz']==str(sel) and r['sigma_relative']=='0.01'][0]
      calculated=effective(z,0)*2/.01**2/SCALE[0]**2
      expected=float(ref['range_information_depth_profiled'])
      err=abs(calculated-expected)/max(expected,1e-30)
      errors.append(err);count+=1
  check('independent eigen/Helmert/lstsq range information reconstruction',max(errors)<1e-7,dict(checks=count,max_relative_error=max(errors)))
  for rep in ('P1','P2'):
   for c in ('C0','C1','C2'):
    for sel in ('2','5','10','JOINT'):
     for sc in scenes:
      rr=[r for r in ir if r['scene_id']==sc['scene_id'] and r['representation']==rep and r['condition']==c and r['delta_hz']==sel]
      if len(rr)==2 and c!='C2':
       first=next(r for r in rr if r['sigma_relative']=='0.01');second=next(r for r in rr if r['sigma_relative']=='0.05')
       check('5percent scale '+sc['scene_id']+rep+c+sel,np.isclose(float(second['range_information_depth_profiled']),.04*float(first['range_information_depth_profiled']),rtol=1e-10))
 summary=dict(status='PASS' if all(c['PASS'] for c in checks) else 'FAIL',checks=len(checks),FAIL=sum(not c['PASS'] for c in checks),records=checks,independent_of_forward_generation=True,method='source hash checks + frozen records + independent eigen/Helmert/gelsd reconstruction when data exist; no new scientific configuration')
 dump(OUT/'VALIDATION.json',summary)
 print(json.dumps({k:summary[k] for k in ('status','checks','FAIL')},indent=2))
 if "--checks-only" not in sys.argv:report(design,dec,summary)
def report(design,dec,validation):
 full=rows(OUT/'E2_FULL_SCENE_INFORMATION.csv');is_evaluated='representation' in full[0]
 status={};details=''
 if is_evaluated:
  for rep in ('P1','P2'):
   primary=[r for r in full if r['representation']==rep and r['condition']=='C1' and r['delta_hz']=='JOINT' and r['sigma_relative']=='0.01']
   r=np.array([float(x['range_information_depth_profiled']) for x in primary]);z=np.array([float(x['depth_information']) for x in primary])
   status[rep]=dict(scenes=len(primary),positive_range_scenes=int(np.sum(r>1e-14)),positive_depth_scenes=int(np.sum(z>1e-10)),range_information_median_m_minus2=float(np.median(r)),range_information_min_m_minus2=float(np.min(r)),depth_information_median_m_minus2=float(np.median(z)))
  cand=rows(OUT/'E2_CANDIDATE_SEPARATION.csv')
  evaluated=[r for r in cand if r['status']=='FINITE_DIAGNOSTIC_NOT_GLOBAL_EXCLUSION']
  sep={}
  for sigma,name in ((.01,'noise_rms_1pct'),(.05,'noise_rms_5pct')):
   ratio=np.array([float(r['C1_CA_rms'])/float(r[name]) for r in evaluated])
   sep[str(sigma)]=dict(compared_candidates=len(ratio),below_one_marginal_noise_rms=int(np.sum(ratio<=1)),median_signal_noise_rms_ratio=float(np.median(ratio)),minimum_ratio=float(np.min(ratio)))
  status['candidate_noise_scale']=sep
  details=f"""
Single HLA C1 joint: {status['P1']}.
Dual HLA C1 joint: {status['P2']}.
Candidate separation versus marginal CA noise RMS: {sep}. Shared-pair covariance was used in FIM; these finite Frobenius ratios are marginal feature distances, not independent-pair likelihood significance or global exclusion.

These pressure/noise/array calculations assume exact navigation and three registered physical snapshots. The angular-state-profiled radial information is a propagation mechanism diagnostic, not an end-to-end bearing-fused estimator or a navigation-robust gain. Rank deficiency implies unbounded uncertainty; positive tiny residual information alone does not establish engineering usefulness. Environment controls are deterministic mismatch distances, not a robustness certificate. Unknown-calibration candidate adjustment is a feasible fit, not a certified nuisance minimum.
"""
 else:
  status={'P1':'NOT_EVALUATED','P2':'NOT_EVALUATED','candidate_noise_scale':'NOT_EVALUATED'}
  details='All24 scenes are retained in each required table with NOT_EVALUATED. Blank values are not zeros. No source/calibration/depth/candidate or physical-information claim follows from a passed plane-wave null.'
 dump(OUT/'AUDITED_SUMMARY.json',status)
 reporttext=f"""# E2-G0 HLA difference-frequency observability

**{dec['decision']}**

Parent:{design['parent_SHA']}. Design:{read(OUT/'EXECUTION_STARTED.json')['design_SHA']}.
Reason:{dec['reason']}

Plane-wave null:{'PASS' if dec['plane_null_PASS'] else 'FAIL'}. Complex-field admission error:{dec['physics_error']}.\nIndependent physical reconstruction:{read(OUT/'E2_INDEPENDENT_PHYSICS_SUMMARY.json') if (OUT/'E2_INDEPENDENT_PHYSICS_SUMMARY.json').exists() else 'NOT_EVALUATED'}.\nShared-pair covariance/FIM admission was NOT_REACHED when physical admission failed. The source/calibration/depth and noise-scale mechanisms are not assessed by the raw-field error checks.
Local execution audit:{validation['status']}, {validation['checks']} checks, {validation['FAIL']} FAIL. This is an execution/integrity audit; research-lead independent scientific acceptance remains pending.
Elapsed:{dec['elapsed_s']:.3f}s. New Monte Carlo0; new received audio0; R4-A1/R4=0%.

{details}

Numerical/physical admission must be PASS before any computed local information is scientifically admitted. If this decision is PHYSICS_OR_COVARIANCE_INCOMPLETE, any retained computed rows are diagnostic only and establish no range/depth capability.

Park2024 and Yin2026 full publisher originals recovered, formula/geometry scope checked and archived. Park supplies directly relevant negative HLA/CZ evidence; Yin's few-element configuration is a980m VLA, not a14m HLA. Source/method comparison:[E2_HLA_PRIMARY_PAPER_COMPARISON.md](E2_HLA_PRIMARY_PAPER_COMPARISON.md). FDSL/POMAP/CMAP strict original comparator remains pending; the evaluated construction is custom per-pair normalized spatial AP, not their complete estimator.

Array assumptions: two newly proposed coherent8-element HLAs,2m spacing14m aperture, global-x orientation, receiver200m and5km offset. Source support150–250Hz; three Δ2/5/10Hz,0.5Hz centers; both members stay inside band. Three snapshot epochs0/600/1200s sample the unchanged original trajectory. This finite mechanism screen does not characterize continuous-time/band coverage or synchronization/association/hardware performance.

Pre-run startup correction: the initial guard rejected a raw-CRLF .gitattributes hash because it used canonical-LF comparison. No scientific execution had started. Exact raw bytes match the freeze; the launcher verifies them under the correct raw-byte convention. Scientific code/model/provider/noise/candidates remain identical to the pushed design. [Record](PRE_EXECUTION_INFRASTRUCTURE_CORRECTION.json). This infrastructure exception is disclosed for lead audit.

No global uniqueness, branch reduction, acoustic recovery, empirical P95, or engineering route pass is declared. STOP after the execution commit/push. E1-G1 remains closed. E2 extraction, depth development, A2, SSP and P5 remain unopened. Hardware UNKNOWN is not the stopping reason.
"""
 (OUT/'E2_G0_REPORT.md').write_text(reporttext,encoding='utf-8')
 (OUT/'GPT_SYNC.md').write_text(reporttext,encoding='utf-8')
 figs=OUT/'figures';figs.mkdir(exist_ok=True)
 plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
 def save(name):
  plt.tight_layout();plt.savefig(figs/(name+'.png'),dpi=180);plt.savefig(figs/(name+'.svg'));plt.close()
  q=figs/(name+'.svg');q.write_text('\n'.join(x.rstrip() for x in q.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 if is_evaluated:
  names=[s['scene_id'] for s in read(OUT/'E2_SCENE_FREEZE.json')['scenes']]
  plt.figure(figsize=(12,4))
  for rep in ('P1','P2'):
   vals=[float(next(r for r in full if r['scene_id']==sid and r['representation']==rep and r['condition']=='C1' and r['delta_hz']=='JOINT' and r['sigma_relative']=='0.01')['range_information_depth_profiled']) for sid in names]
   plt.semilogy(names,vals,'o-',label=rep+' C1 depth-profiled')
  plt.xticks(rotation=65);plt.ylabel('Local radial information (m^-2), sigma1%');plt.legend();plt.title('Exact-navigation physical screen; three snapshots; not estimator performance');save('single_dual_HLA_information')
  plt.figure(figsize=(8,4));cc=rows(OUT/'E2_CALIBRATION_SENSITIVITY.csv')
  plt.bar([r['scene_id'] for r in cc],[float(r['trace_loss_fraction'])*100 for r in cc]);plt.xticks(rotation=65);plt.ylabel('Scaled local FIM trace lost (%)');plt.title('Unknown fixed gain amplitudes; gain phase cancels in AP');save('unknown_calibration_loss')
  selected=[r for r in evaluated if r['scene_id']=='H01_M' and r['scope']=='ABSOLUTE_FULL_REGISTERED_GRID']
  grid=np.array([[float(next(r for r in selected if float(r['range_m'])==rv and float(r['depth_m'])==z)['C1_CA_rms']) for rv in np.arange(50000.,60000.1,1000)] for z in (180,190,200,210,220)])
  plt.figure(figsize=(7,4));plt.imshow(grid,origin='lower',aspect='auto',extent=[49.5,60.5,175,225]);plt.colorbar(label='Feasible-gain-adjusted CA distance');plt.xlabel('Initial range (km)');plt.ylabel('Source depth (m)');plt.title('H01_M preselected complete finite grid; no global exclusion');save('range_depth_compensation')
 else:
  for name,label in [('single_dual_HLA_information','Single / dual HLA information'),('unknown_calibration_loss','Unknown calibration information loss'),('range_depth_compensation','Range-depth candidate separation')]:
   plt.figure(figsize=(8,3));plt.axis('off');plt.text(.5,.7,label,ha='center',fontsize=15);plt.text(.5,.5,'NOT EVALUATED: physical admission did not close',ha='center');plt.text(.5,.25,dec['reason'][:110],ha='center',wrap=True);save(name)

 if (OUT/'E2_FORWARD_CONVERGENCE.csv').exists():
  rr=rows(OUT/'E2_FORWARD_CONVERGENCE.csv')
  mesh=[r for r in rr if r['kind']=='mesh40001_vs20001']
  freq=sorted({float(r['frequency_hz']) for r in mesh})
  maximum=[max(float(r['error']) for r in mesh if float(r['frequency_hz'])==f) for f in freq]
  median=[np.median([float(r['error']) for r in mesh if float(r['frequency_hz'])==f]) for f in freq]
  plt.figure(figsize=(9,4));plt.plot(freq,np.array(maximum)*100,label='max across7 source depths');plt.plot(freq,np.array(median)*100,label='median across7 depths')
  plt.axhline(.2,color='r',ls='--',label='frozen admission0.2%')
  plt.axhline(1,color='k',ls=':',label='registered1% pressure error')
  plt.xlabel('Raw frequency (Hz)');plt.ylabel('Complex relative residual (%)');plt.title('20001 vs40001 mesh; one scalar across all6 registered ranges');plt.legend();save('complex_forward_convergence')
 master=Path('results/R4_MASTER')
 with (master/'R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## E2-G0 physical screen execution\n\n'+dec['decision']+'. [Report](../R4_E2_G0_HLA_DIFFERENCE_INFORMATION/E2_G0_REPORT.md). One bounded deterministic attempt; null and physical admission records retained. No MC/audio; R4=0%. STOP pending lead independent audit; no automatic extraction/depth/A2/SSP/P5.\n')
 with (master/'R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8') as f:f.write('E2-G0-SCREEN,../R4_E2_G0_HLA_DIFFERENCE_INFORMATION/E2_G0_DECISION.json,'+dec['decision']+',Frozen physics screen; all24 retained; no MC/audio,PENDING_RESEARCH_LEAD_AUDIT;STOP;NO_CREDIT,0\n')
 dump(master/'R4_E2_G0_PROGRESS.json',dict(status=dec['decision'],execution='COMPLETE_ONE_BOUNDED_ATTEMPT',lead_audit='PENDING',R4_A1_percent=0,R4_overall_percent=0,new_MC=0,new_audio=0,next='STOP',E2_extraction='NOT_OPENED'))
if __name__=='__main__':audit()
