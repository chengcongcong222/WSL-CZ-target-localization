"""Independent dense covariance audit and presentation; no added scientific design."""
from pathlib import Path
import json,csv,hashlib,subprocess,math
import numpy as np
from scipy.linalg import solve_triangular
O=Path('results/R4_E1_G0_FREQUENCY_INFORMATION')
def load(n):return json.loads((O/n).read_text(encoding='utf-8'))
def rows(n):
 with (O/n).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def dump(n,x):(O/n).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def write(n,g):
 with (O/n).open('w',encoding='utf-8',newline='') as f:w=csv.DictWriter(f,fieldnames=list(g[0]));w.writeheader();w.writerows(g)
def yes(x):return x=='True'
S=np.array([50000.,50000.,2.,2.])
def cold(scene,m):
 # Assemble dense covariance independently of block-plus-rank2 whitening.
 x=np.array(scene['state']);f=np.array([200.] if m==1 else [170.,200.,230.]);neta=m+5
 X=[];N=[];nav=[];sys=[]
 for i in range(121):
  t=10.*i;post=max(t-600,0.)
  node=np.array([2*min(t,600)+2*post*math.cos(math.pi/12),2*post*math.sin(math.pi/12)])
  vn=np.array([2.,0.] if t<=600 else [2*math.cos(math.pi/12),2*math.sin(math.pi/12)])
  for j in range(2):
   pj=node+np.array([0.,j*scene['mirror']*5000])
   d=x[:2]+t*x[2:]-pj;r=np.linalg.norm(d);e=d/r;w=x[2:]-vn;q=float(e@w)
   gp=(np.eye(2)-np.outer(e,e))@w/r
   qb=np.r_[gp,e+t*gp]
   bear=np.array([-e[1],e[0]])/r;bx=np.r_[bear,t*bear]
   XX=[bx*S];NN=[np.zeros(neta)];VV=[-bear];sy=[np.array([math.radians(.05),(-1 if j==0 else 1)*math.radians(.025)])/math.sqrt(3)]
   if i%2:
    u=t-r/1500
    for l,ff in enumerate(f):
     mu=ff*(1-q/1500);fj=-ff/1500*qb
     nn=np.zeros(neta);nn[l]=mu;nn[m]=mu*u/1200;nn[m+1+2*j]=mu;nn[m+2+2*j]=mu*t/1200
     XX.append(fj*S);NN.append(nn);VV.append(-fj[:2]);sy.append(np.zeros(2))
   X.extend(XX);N.extend(NN);nav.append(np.array(VV));sys.extend(sy)
 X=np.array(X);N=np.array(N);R=np.zeros((len(X),len(X)));ix=0
 for V in nav:
  k=len(V);R[ix:ix+k,ix:ix+k]=np.diag([math.radians(.05)**2]+[.001**2]*(k-1))+25**2*V@V.T;ix+=k
 U=np.array(sys);R+=U@U.T
 L=np.linalg.cholesky(R);A=solve_triangular(L,X,lower=True);E=solve_triangular(L,N,lower=True)
 norms=np.linalg.norm(E,axis=0);u,s,v=np.linalg.svd(E/norms,full_matrices=False);rank=sum(s>s[0]*1e-10)
 Q=A-u[:,:rank]@(u[:,:rank].T@A);F=Q.T@Q
 u,sv,V=np.linalg.svd(Q,full_matrices=False)
 g=np.r_[0.,0.,x[2:]/np.linalg.norm(x[2:])]*S
 var=float(np.sum((V@g/sv)**2))
 return F,var
def main():
 full=rows('E1_G0_FULL_SCENE_RESULTS.csv');summary=rows('E1_G0_INFORMATION_SUMMARY.csv');scenes=load('E1_G0_SCENE_BINDINGS.json')['scenes']
 checks=[]
 def check(n,ok,value=None):checks.append(dict(check=n,passed=bool(ok),value=value))
 check('1896_complete_records',len(full)==1896)
 check('24_original_scenes',len(scenes)==24)
 p=load('E1_G0_DESIGN_FREEZE.json');b=load('E1_G0_SCENE_BINDINGS.json')
 for path,sha in {**p['code_binding'],**b['bindings']}.items():
  check('frozen_binding:'+path,hashlib.sha256(Path(path).read_bytes().replace(b'\r\n',b'\n')).hexdigest()==sha)
 matrices=np.load(O/'E1_G0_EFFECTIVE_INFORMATION_MATRICES.npz')
 dense=[]
 for scene in scenes:
  for m in [1,3]:
   F,var=cold(scene,m);key=scene['scene_id']+'_NEW_N2_'+str(m)+'_0.001_0'
   r=next(r for r in full if r['scene_id']==scene['scene_id'] and r['variant']=='NEW' and r['model']=='N2' and int(r['lines'])==m and float(r['sigma_f_Hz'])==.001 and r['corner']=='0')
   ferr=float(np.linalg.norm(F-matrices[key])/np.linalg.norm(F));verr=abs(var-float(r['variance_mps2']))/var
   check('dense_F:'+key,ferr<2e-6,ferr);check('dense_variance:'+key,verr<2e-6,verr)
   dense.append(dict(scene_id=scene['scene_id'],lines=m,dense_variance=var,relative_F_error=ferr,relative_variance_error=verr))
 write('E1_G0_INDEPENDENT_DENSE_AUDIT.csv',dense)
 bv={r['scene_id']:float(r['variance_mps2']) for r in full if r['variant']=='B0'}
 for ss in summary:
  g=[r for r in full if all(r[k]==ss[k] for k in ['variant','model','lines','sigma_f_Hz','corner'])]
  check('summary_count:'+str(tuple(ss[k] for k in ['variant','model','lines','sigma_f_Hz','corner'])),len(g)==24)
  vs=np.array([float(r['variance_mps2']) for r in g]);delta=np.array([1-float(r['variance_mps2'])/bv[r['scene_id']] for r in g])
  check('summary_median',abs(float(np.median(delta))-float(ss['median_variance_reduction']))<1e-12)
  check('summary_worst',abs(max(vs)/max(bv.values())-float(ss['worst_variance_ratio']))<1e-12)
  gate=float(np.median(delta))>=.2 and max(vs)/max(bv.values())<=1.000000001 and all(yes(r['rank_stable']) for r in g)
  check('summary_gate',gate==yes(ss['gate_pass']))
 d=load('E1_G0_DECISION.json');v=load('VALIDATION.json')
 check('VALIDATION_PASS',v['valid'])
 check('NO_PACKAGE_MEETS_20_PERCENT',not any(yes(z['gate_pass']) for z in summary))
 check('correct_negative_scientific_decision',d['route_decision']=='E1_G0_FREQUENCY_INCREMENT_NOT_ESTABLISHED')
 check('zero_MC_audio_R4',all(d[k]==0 for k in ['new_Monte_Carlo','new_received_audio','new_propagation','R4_percent','R4_A1_percent']))
 fd=load('E1_G0_DIFFERENTIATION_VALIDATION.json')
 dump('E1_G0_INDEPENDENT_AUDIT.json',dict(checks=checks,pass_count=sum(z['passed'] for z in checks),fail_count=sum(not z['passed'] for z in checks),
   max_dense_F_relative=max(z['relative_F_error'] for z in dense),max_dense_variance_relative=max(z['relative_variance_error'] for z in dense),
   max_target_FD_error=max(max(z['target_step_errors']) for z in fd['checks']),max_nuisance_FD_error=max(max(z['nuisance_step_errors']) for z in fd['checks'])))
 assert all(z['passed'] for z in checks),[z for z in checks if not z['passed']]
 # Standalone publication-style figures; figures report all sigma/line packages.
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 figs=O/'figures';figs.mkdir(exist_ok=True)
 plt.rcParams['svg.hashsalt']='E1-G0'
 def save(fig,n):
  fig.tight_layout();fig.savefig(figs/(n+'.png'),dpi=180);fig.savefig(figs/(n+'.svg'),metadata={'Date':None});plt.close(fig); svg=figs/(n+'.svg');svg.write_text('\n'.join(z.rstrip() for z in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
 ids=[s['scene_id'] for s in scenes]
 fig,axes=plt.subplots(2,3,figsize=(17,9),sharex=True)
 for row,m in enumerate([1,3]):
  for col,sig in enumerate([.001,.005,.010]):
   ax=axes[row,col]
   for model in ['B0','N0','N1','N2','N3','N4_BOTH','N4_LINE','N5']:
    vals=[]
    for sid in ids:
     g=[r for r in full if r['scene_id']==sid and r['model']==model and (model=='B0' or (int(r['lines'])==m and float(r['sigma_f_Hz'])==sig and r['corner']=='0' and r['variant'] in ['DUAL_FREQ','NEW']))]
     vals.append(1/float(g[0]['variance_mps2']))
    ax.plot(range(24),vals,label=model,alpha=.8,linewidth=1.1)
   ax.set(title=f'{m} lines; ideal sigma_f={sig:g} Hz',ylabel='Effective speed information (s²/m²)',xlabel='H01-12, paired mirrors')
   ax.grid(alpha=.25)
 axes[0,0].legend(ncol=2,fontsize=8)
 save(fig,'FIG1_ALL_SCENES_EFFECTIVE_SPEED_INFORMATION')
 fig,axes=plt.subplots(1,3,figsize=(17,5))
 for ax,sig in zip(axes,[.001,.005,.010]):
  models=['N0','N1','N2','N3','N4_BOTH','N4_LINE','N5']
  for m in [1,3]:
   vals=[100*float(next(r for r in summary if r['model']==model and int(r['lines'])==m and float(r['sigma_f_Hz'])==sig and r['corner']=='0' and r['variant'] in ['DUAL_FREQ','NEW'])['median_variance_reduction']) for model in models]
   ax.plot(models,vals,'o-',label=f'{m} lines')
  ax.set(title=f'Ideal sigma_f={sig:g} Hz',ylabel='Median speed variance reduction vs B0 (%)')
  ax.tick_params(axis='x',rotation=45);ax.legend();ax.grid(alpha=.25)
 fig.suptitle('Nuisance losses; every package below the frozen 20% route threshold',fontsize=12)
 save(fig,'FIG2_SOURCE_REFERENCE_NUISANCE_LOSS')
 lines=['# E1-G0 未知量剖面后的频率信息验证','',d['route_decision'],'',
 '已复用H01–H12与双镜像24场景，1896条全量记录；无MC、录音、频率实际提取、传播或优化器运行。没有因Sun完整估计器未锁而停止观测机制研究。','',
 '24场景为既有5km固定全局横向布放的 nominal beta=0 控制；不代表500个保存随机布放的完整统计覆盖。原方位随机RMS0.05°、共模±0.05°、差模半差±0.025°、导航25m/轴均保留为局部Gaussian/二阶矩条件。','',
 '| 模型 | 单线：0.001 / 0.005 / 0.010 Hz，中位改善% | 三线：同档中位改善% |','|---|---|---|']
 for model in ['N0','N1','N2','N3','N4_BOTH','N4_LINE','N5']:
  vals=[]
  for m in [1,3]:
   vals.append(' / '.join(f"{100*float(next(r for r in summary if r['model']==model and int(r['lines'])==m and float(r['sigma_f_Hz'])==sig and r['corner']=='0' and r['variant'] in ['DUAL_FREQ','NEW'])['median_variance_reduction']):.6f}" for sig in [.001,.005,.010]))
  lines.append('| '+model+' | '+' | '.join(vals)+' |')
 n2=next(r for r in summary if r['model']=='N2' and r['lines']=='3' and float(r['sigma_f_Hz'])==.001 and r['corner']=='0' and r['variant']=='NEW')
 lines+=['',f"最精细登记误差档三线N2的中位方差下降={100*float(n2['median_variance_reduction']):.6f}%，最坏方差比={float(n2['worst_variance_ratio']):.9f}，有正增量24/24。改善为正但远小于20%，不能写成完全零信息。全部6个单/三线-精度资源包及N2五角点均未达到主Gate；N0/N5也未达到，故不是校准条件路线通过。",'',
 f"B0有效速度方差范围={d['B0_variance_range'][0]:.9f}–{d['B0_variance_range'][1]:.9f} (m/s)²。这是局部Gaussian二阶矩机制尺度，非经验P95、非最终10% Gate。",'',
 '## 信息来源与模型边界','',
 'SUN-MEASUREMENT为Sun Eq4未知F的单节点频率加当前双bearing控制，PAPER_ADAPTED_MEASUREMENT_INFORMATION_CONTROL；不加原论文频率随机游走过程先验，不报告MFB-AUKF RMSE或复现成功率。SINGLE-FREQ为单节点频率且采用N2未知量；DUAL-FREQ=N0双节点频率。NEW逐步增加源/参考未知量。共享F跨节点不是known emitted f0；局部truth状态/名义频率仅用于计算信息，不是估计器输入，本轮无估计器。','',
 '主N2未知logF、共享源漂移、独立节点fractional常数/线性参考；N3独立additive Hz控制；N4_BOTH同时允许两类参考，N4_LINE新增线特有漂移；N5两节点fractional频偏sigma5e-6、漂移sigma5e-9/s的独立外部校准假设。N5没有实际设备/独立标定证据，不给源F先验。','',
 '原方位偏差的bounded uniform仅以bound²/3二阶矩进入Gaussian信息；不是精确uniform likelihood CRLB。导航误差在同节点/时刻的bearing、各线frequency间造成相关，已联合纳入Σ；理想frequency白误差与bearing仪器白误差独立仅为G0观测假设，不代表同一录音衍生统计量已独立。','',
 'c=1500已知、单有效LOS传播，源时间u=t-r/c-delta_t为准静态近似；当前没有first-CZ模态传播。时间偏差是登记角点条件，不是被独立剖面的同步未知量，不能据此宣称同步或目标关联已验证。弱方向、奇异值、rank阈值敏感性、所有误差档和失败包完整保存。','',
 '三线可以提高精度，实际局部四状态秩仍为4，不能说三线增加三个独立运动维度。N4_LINE与共享模型的数值接近源于同源公共Doppler结构与已有线常数/共享时序约束；只限本模型，不推成任意漂移鲁棒。','',
 '## 数值证据','',
 f"7项零信息/结构控制通过；60项解析-有限差分交叉检查通过，最大target相对误差={max(max(z['target_step_errors']) for z in fd['checks']):.3e}，最大nuisance相对误差={max(max(z['nuisance_step_errors']) for z in fd['checks']):.3e}。",
 f"投影与Schur最大相对差={v['max_schur_projection_relative']:.3e}；所有登记rank敏感性1e-8/1e-10/1e-12的有效speed方差稳定。逐时任意节点reference控制相对B0矩阵差={v['zero_controls'][4]['value']:.3e}，不可估speed零方向输出UNBOUNDED。",
 f"独立稠密Σ重建24场景×单/三线共48项N2控制：F最大相对差={max(z['relative_F_error'] for z in dense):.3e}、方差最大相对差={max(z['relative_variance_error'] for z in dense):.3e}；与主程序的block+rank2白化算法独立。汇总与Gate逐组重建，0 FAIL。",
 '', '## 路线判定与停止','',
 '当前表示与1200s几何中，剖面未知源/参考后的非仿射频率信息太弱，未达到研发选路门槛。不继续用换滤波器、给true F、改频率精度、延长时间或调整布放救结果。该结论不是普遍物理不可辨识，也不是所有频率/传播方法都无效。','',
 'FREQUENCY_OBSERVABILITY_MECHANISM_SCREEN_COMPLETED；R4-A1=0%、R4=0%。建议下一阶段E2_REVIEW，须研究负责人独立审计与另行授权。E1-G1、E2、depth、A2、SSP、P5均未开放。','',
 '![全部场景有效信息](figures/FIG1_ALL_SCENES_EFFECTIVE_SPEED_INFORMATION.png)','',
 '![未知量逐步增加的信息损失](figures/FIG2_SOURCE_REFERENCE_NUISANCE_LOSS.png)','',
 '[全量表](E1_G0_FULL_SCENE_RESULTS.csv) · [汇总](E1_G0_INFORMATION_SUMMARY.csv) · [秩审计](E1_G0_NUMERICAL_RANK_AUDIT.csv) · [独立审计](E1_G0_INDEPENDENT_AUDIT.json) · [判定](E1_G0_DECISION.json)','']
 report='\n'.join(lines)
 (O/'E1_G0_REPORT.md').write_text(report,encoding='utf-8');(O/'GPT_SYNC.md').write_text(report,encoding='utf-8')
 master=Path('results/R4_MASTER')
 dumpdata=dict(stage=d['stage'],design_SHA=d['design_SHA'],parent_SHA=d['parent_SHA'],status=d['route_decision'],mechanism_screen='FREQUENCY_OBSERVABILITY_MECHANISM_SCREEN_COMPLETED',
  independent_audit='PENDING_RESEARCH_LEAD_AUDIT',R4_A1_percent=0,R4_percent=0,G1='NOT_OPENED',E2='NOT_OPENED',recommended_next='E2_REVIEW',automatic_next=False,new_MC=0,new_audio=0)
 (master/'R4_E1_G0_PROGRESS.json').write_text(json.dumps(dumpdata,indent=2)+'\n',encoding='utf-8')
 with (master/'R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\n\n## E1-G0 frequency mechanism screen\n\nSource measurement lock ACCEPTED_WITH_ORIGINAL_ESTIMATOR_EXCEPTION at f9b994a; G0 admitted, G1 not admitted. Design46242899 frozen/pushed before execution. [Report](../R4_E1_G0_FREQUENCY_INFORMATION/E1_G0_REPORT.md). 24 nominal mirror scenes,1896 deterministic information records; no MC/audio/propagation. E1_G0_FREQUENCY_INCREMENT_NOT_ESTABLISHED: best registered 3-line1mHz N2 C0 median variance reduction0.125681%, below20%; N0/N5 also below gate. All corners/packages retained. Numerical controls,finite differences,dense covariance reconstruction pass. Mechanism screen completed, no P95 claim; R4-A1/R4=0%. Recommended E2_REVIEW only; G1/E2/depth/A2/SSP/P5 unopened. Push execution commit then STOP for independent audit.\n')
 with (master/'R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:
  w=csv.writer(f);w.writerow(['E1-MEASUREMENT-LOCK-ACCEPTED','../R4_OBSERVABILITY_FRONTIER_RESEARCH_RESET/GPT_SYNC_E1_LOCK.md','E1_PRIMARY_MEASUREMENT_LOCK_ACCEPTED_WITH_ORIGINAL_ESTIMATOR_EXCEPTION','Original Sun Eq4 accepted; full estimator exception; G0 admitted; G1 not admitted','RESEARCH_LEAD_ACCEPTED_AT_f9b994a',0])
  w.writerow(['E1-G0-SCREEN','../R4_E1_G0_FREQUENCY_INFORMATION/E1_G0_DECISION.json',d['route_decision'],'24 nominal scenes; covariance/rank controls; no MC/audio; weak positive increment below20%','PENDING_RESEARCH_LEAD_AUDIT;NO_JOINT_GATE_CREDIT',0])
 files=sorted(p for p in O.rglob('*') if p.is_file() and p.name!='E1_G0_RESULT_MANIFEST.json')
 dump('E1_G0_RESULT_MANIFEST.json',dict(design_SHA=d['design_SHA'],files={p.as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files},hash_scope='raw worktree bytes; frozen code/source bindings separately normalize CRLF',new_MC=0,new_audio=0))
 print(json.dumps(dict(independent_checks=len(checks),failed=sum(not z['passed'] for z in checks),max_dense_variance_error=max(z['relative_variance_error'] for z in dense)),indent=2))
if __name__=='__main__':main()
