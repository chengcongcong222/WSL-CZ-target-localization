"""One independent cold audit. Does not import main/generation/model/statistics code."""
from pathlib import Path
import json,csv,hashlib,math,time,collections
import numpy as np
from independent_math import rebuild_panel,residual_and_jac,initial,fisher_design
O=Path(__file__).resolve().parent;D=O.parent;ROOT=O.parents[2];DATA=Path('D:/ProjectStorage/WSL-CZ/CZ_D2/RESUME_01_DATA')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def csvread(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  while b:=f.read(1048576):h.update(b)
 return h.hexdigest()
def ah(*arrays):
 h=hashlib.sha256()
 for a in arrays:
  a=np.ascontiguousarray(a);h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes())
 return h.hexdigest()
def emit(p,d):
 with Path(p).open('w',encoding='utf-8') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def num(x):return None if x in ('NOT_EVALUATED','',None) else float(x)
def eq(a,b,rtol=1e-10,atol=1e-10):
 if a is None or b is None:return a is b
 return bool(np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=False))
counts=collections.Counter();fail=[];fail_count=0

def ck(group,ok,detail=''):
 global fail_count
 counts[group]+=1
 if not ok:
  fail_count+=1
  if len(fail)<1000:fail.append(dict(group=group,detail=str(detail)))
def rankq(v,q):
 if len(v)==0:return None
 return float(sorted(v)[math.ceil(q*len(v))-1])
def idx(q,b,c,r,p):return ((q*2+b)*6000+c*500+r)*2+p
M=['range_abs_km','range_rel','bearing_circular_deg','speed_abs_mps','speed_rel','heading_circular_deg']

def main():
 tic=time.monotonic();cfg=read(D/'DESIGN_FREEZE.json');binding=read(O/'EFFECTIVE_INPUT_BINDING.json');reg=csvread(Path(binding['registry']['path']))
 ck('binding',sha(Path(binding['registry']['path']))==binding['registry']['sha256'] and len(reg)==762)
 for x in read(O/'PRODUCTION_CODE_FREEZE.json')['files']:ck('code_identity',sha(Path(x['path']))==x['sha256'],x['path'])
 saved=read(DATA/'private/PANEL_AND_REJECTIONS.json');reconstruction=rebuild_panel(cfg,reg);ck('panel_rejection_sequence',saved['states']==reconstruction['states'] and saved['attempts']==reconstruction['attempts']);states=np.array(saved['states'])
 ck('panel_seal',sha(DATA/'private/PANEL_AND_REJECTIONS.json')==read(O/'PANEL_SEAL.json')['sha256'])
 t=np.linspace(0,1200,121);post=np.maximum(0,t-600);pm=np.c_[2*np.minimum(t,600)+2*post*np.cos(np.pi/12),2*post*np.sin(np.pi/12)];sides=np.where(np.arange(500)<250,-1,1)
 manifest=csvread(O/'PUBLIC_RECORD_IDENTITIES.csv');lookup={(x['quality'],x['geometry'],x['case_id'],int(x['replicate'])):x for x in manifest};ck('public_record_count',len(lookup)==len(manifest)==24000)
 for x in read(O/'PUBLIC_ARCHIVE_SEAL.json')['files']:ck('public_file_SHA',sha(Path(x['path']))==x['sha256'],x['path'])
 priv={x['case']:x for x in read(O/'PRIVATE_ARCHIVE_INDEX.json')};ck('innovation_count',len(priv)==12)
 raw=np.load(DATA/'terminals.npy',mmap_mode='r');ck('terminal_full_plan_shape',raw.shape==(48000,));errors=np.full((48000,6),np.inf)
 with np.load(O/'TERMINALS_ARCHIVE.npz',allow_pickle=False) as pack:ck('terminal_pack_identity',np.array_equal(pack['terminals'],raw))
 maxima=dict(cost_difference=0.,initial_difference=0.,trajectory_difference_m=0.,bearing_difference_rad=0.)
 for ci in range(12):
  ip=Path(priv[ci]['path']);ck('private_file_SHA',sha(ip)==priv[ci]['sha256'])
  with np.load(ip,allow_pickle=False) as data:uu=data['uniforms'];nn=data['normals'];ss=data['side']
  ck('innovation_schema',uu.shape==(500,2) and nn.shape==(500,121,6) and np.array_equal(ss,sides) and uu.dtype==nn.dtype==np.dtype('float64') and np.isfinite(nn).all() and np.all((uu>=0)&(uu<1)))
  ih=ah(uu,nn,ss);ck('innovation_array_SHA',ih==priv[ci]['array_sha256'])
  truth=states[ci];target=np.array([truth[0]*1000*np.cos(np.deg2rad(truth[1])),truth[0]*1000*np.sin(np.deg2rad(truth[1]))])+t[:,None]*truth[2]*np.array([np.cos(np.deg2rad(truth[3])),np.sin(np.deg2rad(truth[3]))])
  for qi,Q in enumerate(['Q1','Q2']):
   q=cfg['quality'][Q];mainhashes=None
   for bi,B in enumerate(['B0','B3']):
    path=DATA/'public'/f'{Q}_{B}_D{ci+1:02d}.npz'
    with np.load(path,allow_pickle=False) as data:tt=data['times'];PP=data['reported'];yy=data['bearing']
    ck('public_schema',np.array_equal(tt,t) and PP.shape==(500,121,2,2) and yy.shape==(500,121,2) and PP.dtype==yy.dtype==np.dtype('float64'))
    mainhs=[]
    for rep in range(500):
     key=(Q,B,f'D{ci+1:02d}',rep);row=lookup[key];mainh=ah(t,PP[rep,:,0],yy[rep,:,0]);recordh=ah(t,PP[rep],yy[rep]);mainhs.append(mainh)
     ck('record_array_identity',row['record_sha256']==recordh and row['main_sha256']==mainh and row['innovation_sha256']==ih and int(row['side'])==int(sides[rep]),key)
     if ci==0:
      true=np.stack([pm,pm+np.array([cfg['geometry'][B]['forward_km']*1000,sides[rep]*cfg['geometry'][B]['lateral_km']*1000])],axis=1);wantP=true+q['nav_sigma_m']*nn[rep,:,2:].reshape(121,2,2);delta=target[:,None,:]-true;c=(2*uu[rep,0]-1)*q['common_bound_deg'];d=(2*uu[rep,1]-1)*q['half_diff_bound_deg'];wanty=np.arctan2(delta[...,1],delta[...,0])+np.deg2rad(c+d*np.array([-1,1]))+np.deg2rad(q['bearing_sigma_deg'])*nn[rep,:,:2];wanty=(wanty+np.pi)%(2*np.pi)-np.pi
      dp=float(abs(wantP-PP[rep]).max());dy=float(abs(np.arctan2(np.sin(wanty-yy[rep]),np.cos(wanty-yy[rep]))).max());maxima['trajectory_difference_m']=max(maxima['trajectory_difference_m'],dp);maxima['bearing_difference_rad']=max(maxima['bearing_difference_rad'],dy);ck('D01_full500_reconstruction',dp<=1e-8 and dy<=1e-12,key)
     initial_h,fb,ct=initial(t,PP[rep],yy[rep])
     for pi,P in enumerate(['P0','P1']):
      i=idx(qi,bi,ci,rep,pi);z=raw[i]
      if not z['executed']:continue
      qual=dict(q);qual['common_bound_deg']*=1+pi;qual['half_diff_bound_deg']*=1+pi
      h=z['h'];rr,J=residual_and_jac(h,t,PP[rep],yy[rep],qual);cost=.5*float(np.dot(rr,rr));maxima['cost_difference']=max(maxima['cost_difference'],abs(cost-z['cost']));ck('terminal_cost',eq(cost,z['cost'],1e-9,1e-9),i)
      ok=z['solver_status']>0 and np.isfinite(h).all() and np.isfinite(rr).all();ck('acceptance',bool(z['accepted'])==ok,i);ck('nfev_limit',0<z['nfev']<=160,i)
      md=float(abs(initial_h-z['initial']).max());maxima['initial_difference']=max(maxima['initial_difference'],md);ck('public_initial_identity',eq(initial_h,z['initial']) and fb==bool(z['fallback']) and ct==z['intersections'],i)
      low=np.array([45,-5,1,-15,-qual['common_bound_deg'],-qual['half_diff_bound_deg']]);span=np.array([15,10,2,30,2*qual['common_bound_deg'],2*qual['half_diff_bound_deg']]);unit=(h-low)/span;near=(unit<=1e-6)|(unit>=1-1e-6);ck('physical_bounds',np.all(unit>=-1e-12) and np.all(unit<=1+1e-12),i);ck('bound_flags',np.array_equal(near,z['near'].astype(bool)),i)
      original=np.array([q['common_bound_deg'],q['half_diff_bound_deg']]);active=original*(1+pi);ck('physical_bias_distances',eq(active-abs(h[4:]),z['active_bias_distance']) and eq(original-abs(h[4:]),z['original_bias_distance']) and np.array_equal(abs(h[4:])>original,z['outside_original'].astype(bool)),i)
      sv=np.linalg.svd(J*span,compute_uv=False);ck('fit_Jacobian_diagnostics',int(np.sum(sv/sv[0]>1e-10))==z['rank'] and eq(sv[-1]/sv[0],z['ratio'],1e-8,1e-12),i);ck('residual_rms',eq(np.sqrt(np.mean(rr[:242]**2)),z['rms']),i)
      if ok:
       dr=abs(h[0]-truth[0]);dv=abs(h[2]-truth[2]);angles=[abs((h[j]-truth[j]+180)%360-180) for j in (1,3)];errors[i]=[dr,dr/truth[0],angles[0],dv,dv/truth[2],angles[1]]
    if bi==0:mainhashes=mainhs
    else:ck('B0_B3_main_bytes',mainhashes==mainhs,(Q,ci))
  print('AUDIT_TERMINALS_CASE',ci+1,flush=True)
 with np.load(O/'EVALUATION_ARRAYS.npz') as e:ck('all_error_reconstruction',eq(errors,e['errors']) and np.array_equal(raw['executed'],e['executed']) and np.array_equal(raw['accepted'],e['accepted']))
 statuses=csvread(O/'CONFIGURATION_EXECUTION_STATUS.csv');ck('all48000_status_rows',len(statuses)==48000 and {int(x['configuration_id']) for x in statuses}==set(range(48000)))
 for row in statuses:
  i=int(row['configuration_id']);expected='NOT_EVALUATED' if not raw[i]['executed'] else 'ACCEPTED' if raw[i]['accepted'] else 'SOLVER_FAILED';ck('execution_status_identity',row['status']==expected)
 access=read(O/'ESTIMATOR_FILE_ACCESS.json');ck('private_read_isolation',access['denied']==[] and all('private' not in x.lower() and 'panel' not in x.lower() for x in access['paths'] if x.startswith('READ:')))
 actualpublic={x[5:] for x in access['paths'] if x.startswith('READ:') and x.endswith('.npz')};ck('all_public_files_actually_read',actualpublic=={str(DATA/'public'/f'{Q}_{B}_D{c+1:02d}.npz') for Q in ['Q1','Q2'] for B in ['B0','B3'] for c in range(12)})
 def group(config,case,side):
  Q,B,P=config.split('_');q=['Q1','Q2'].index(Q);b=['B0','B3'].index(B);p=['P0','P1'].index(P);cs=range(12) if case=='POOLED' else [int(case[1:])-1];reps=range(500) if side=='BOTH' else range(250) if side=='-1' else range(250,500);ids=np.array([idx(q,b,c,r,p) for c in cs for r in reps]);return ids,ids[raw['executed'][ids]==1]
 table=csvread(O/'STATISTICS.csv');ck('statistics_row_count',len(table)==3744)
 for row in table:
  planned,done=group(row['config'],row['case_id'],row['side']);sel=done if row['mode']=='unconditional' else done[raw['accepted'][done]==1];mi=M.index(row['metric']);v=errors[sel,mi]
  ck('stat_denominator',int(row['planned'])==len(planned) and int(row['n'])==len(sel) and (row['complete_panel']=='True')==(len(done)==len(planned)))
  for label,q in [('P50',.5),('P90',.9),('P95',.95),('P99',.99)]:ck('all_quantiles',eq(rankq(v,q),num(row[label])),(row['config'],row['case_id'],row['side'],row['metric'],row['mode'],label))
  ck('maximum',eq(float(max(v)) if len(v) else None,num(row['maximum'])))
 diag=csvread(O/'OUTPUT_AND_BOUNDARY_DIAGNOSTICS.csv');ck('diagnostic_count',len(diag)==312)
 for row in diag:
  allid,done=group(row['config'],row['case_id'],row['side']);valid=done[raw['accepted'][done]==1];n=len(allid);ck('diagnostic_denominator',int(row['planned'])==n and int(row['executed'])==len(done) and int(row['accepted'])==len(valid) and int(row['failed'])==len(done)-len(valid) and int(row['not_evaluated'])==n-len(done));ck('rates',eq(float(row['output_rate']),len(valid)/n) and eq(float(row['failed_rate']),(len(done)-len(valid))/n))
  for name,arr in [('state_at_bound',np.any(raw['near'][done,:4],axis=1)),('bias_at_bound',np.any(raw['near'][done,4:],axis=1)),('common_at_bound',raw['near'][done,4]),('diff_at_bound',raw['near'][done,5]),('outside_original_bias',np.any(raw['outside_original'][done],axis=1)),('initial_fallback',raw['fallback'][done]),('iteration_limit',raw['solver_status'][done]==0),('weak_fit_rank',raw['rank'][done]<6)]:ck('diagnostic_count_fields',int(row[name+'_count'])==int(np.sum(arr)))
  for mi,limit in [(1,.10),(2,.5),(4,.20),(5,5)]:bad=int(np.sum(errors[valid,mi]>limit));ck('tail_counts',int(row[M[mi]+'_large_success_count'])==bad and eq(float(row[M[mi]+'_large_or_failed_fraction']),(bad+len(done)-len(valid))/n))
 for row in csvread(O/'AGGREGATES.csv'):
  mi=M.index(row['metric']);sc=[];complete=True
  for c in range(12):
   ids,done=group(row['config'],f'D{c+1:02d}','BOTH');complete &= len(ids)==len(done);sc.append(errors[done,mi])
  if complete:
   med=[rankq(v,.5) for v in sc]
   for lab,q in [('P50',.5),('P90',.9),('P95',.95),('P99',.99)]:ck('worst_scene_quantile',eq(max(rankq(v,q) for v in sc),num(row['worst_scene_'+lab])))
   ck('equal_weight',eq(sum(med)/12,num(row['mean_of_12_scene_P50'])) and eq(rankq(med,.5),num(row['nearest_rank_median_of_12_scene_P50'])))
  else:ck('incomplete_global_NOT_EVALUATED',all(row[k]=='NOT_EVALUATED' for k in row if k.startswith('worst_scene_')))
  ids,done=group(row['config'],'POOLED','BOTH');ck('pooled',eq(rankq(errors[done,mi],.5),num(row['pooled_P50'])) and eq(rankq(errors[done,mi],.95),num(row['pooled_P95'])))
 pairs=csvread(O/'PAIRED_COMPARISONS.csv');ck('paired_rows',len(pairs)==1872)
 for row in pairs:
  aa,ad=group(row['baseline'],row['case_id'],row['side']);bb,bd=group(row['candidate'],row['case_id'],row['side']);ok=(raw['executed'][aa]==1)&(raw['executed'][bb]==1);mi=M.index(row['metric']);x=errors[aa[ok],mi];y=errors[bb[ok],mi];ff=np.isfinite(x)&np.isfinite(y)
  want=dict(planned=len(aa),paired_executed=len(x),unpaired_not_evaluated=len(aa)-len(x),better=int(np.sum(y<x)),worse=int(np.sum(y>x)),equal=int(np.sum(y==x)),both_failed=int(np.sum(~np.isfinite(x)&~np.isfinite(y))),baseline_failed_candidate_valid=int(np.sum(~np.isfinite(x)&np.isfinite(y))),candidate_failed_baseline_valid=int(np.sum(np.isfinite(x)&~np.isfinite(y))),both_finite=int(ff.sum()),ratio_n=int(np.sum(ff&(x>0))))
  ck('paired_integer_counts',all(int(row[k])==v for k,v in want.items()));ck('paired_difference',eq(rankq(y[ff]-x[ff],.5),num(row['finite_difference_P50'])) and eq(rankq(y[ff]-x[ff],.95),num(row['finite_difference_P95'])));rs=ff&(x>0);ck('paired_ratio',eq(rankq(y[rs]/x[rs],.5),num(row['finite_ratio_P50'])))
 # Independent nominal information reconstruction, never Uniform covariance as prior.
 if (O/'INFORMATION.csv').exists():
  it=csvread(O/'INFORMATION.csv');pack=np.load(O/'INFORMATION_MATRICES.npz');ck('information_rows',len(it)==192);scale=np.array(cfg['information']['scale'])
  for row in it:
   h=states[int(row['case_id'][1:])-1];state=np.r_[h[0]*1000*np.array([np.cos(np.deg2rad(h[1])),np.sin(np.deg2rad(h[1]))]),h[2]*np.array([np.cos(np.deg2rad(h[3])),np.sin(np.deg2rad(h[3]))]),0.,0.];geom=cfg['geometry'][row['geometry']];P=np.stack([pm,pm+[geom['forward_km']*1000,int(row['side'])*geom['lateral_km']*1000]],axis=1);q=cfg['quality'][row['quality']];A,H,V,var=fisher_design(state,t,P,q,scale);n=4 if row['bias_known']=='True' else 6;AA=A[:,:n];F=AA.T@AA;k=int(row['matrix_index']);ck('Fisher_matrix',eq(F,pack['Fscaled'][k,:n,:n],1e-8,1e-10),k)
   sv=np.linalg.svd(AA,compute_uv=False);ck('Fisher_rank',int(np.sum(sv/sv[0]>1e-10))==int(row['rank']),k);nuis=[0,1] if n==4 else [0,1,4,5];u,s,vh=np.linalg.svd(AA[:,nuis],full_matrices=False);rk=int(np.sum(s/s[0]>1e-10));Z=AA[:,2:4]-u[:,:rk]@(u[:,:rk].T@AA[:,2:4]);ck('profiled_velocity_matrix',eq(Z.T@Z,pack['velocity_effective'][k],1e-8,1e-10),k)
   # Finite differences from independent mean/variance, not primary information function.
   dm=0.;dv=0.
   for axis in range(6):
    step=scale[axis]*1e-5;plus=state.copy();minus=state.copy();plus[axis]+=step;minus[axis]-=step
    def mv(x):
     d=(x[:2]+t[:,None]*x[2:4])[:,None,:]-P;ss=np.sum(d*d,axis=-1);return (np.arctan2(d[...,1],d[...,0])+x[4]+x[5]*np.array([-1,1])).ravel(),((q['bearing_sigma_deg']*np.pi/180)**2+q['nav_sigma_m']**2/ss).ravel()
    ap,vp=mv(plus);am,vm=mv(minus);da=np.arctan2(np.sin(ap-am),np.cos(ap-am))/(2*step);dvar=(vp-vm)/(2*step);dm=max(dm,float(abs((da-H[:,axis])*scale[axis]).max()));dv=max(dv,float(abs((dvar-V[:,axis])*scale[axis]).max()/var.max()))
   ck('independent_Fisher_FD',dm<=1e-7 and dv<=1e-6,k)
 else:ck('information_missing_scope',np.sum(raw['executed'])<48000,'NOT_EVALUATED')
 if (O/'MAIN_OUTPUT_SEAL.json').exists():
  for x in read(O/'MAIN_OUTPUT_SEAL.json')['files']:ck('main_output_seal',sha(Path(x['path']))==x['sha256'],x['path'])
 emit(O/'COLD_AUDIT.json',dict(PASS=fail_count==0,checks=sum(counts.values()),PASS_count=sum(counts.values())-fail_count,FAIL_count=fail_count,groups=dict(counts),failures=fail,maximum_differences=maxima,executed=int(raw['executed'].sum()),not_evaluated=48000-int(raw['executed'].sum()),scope='INDEPENDENT_ARCHIVE_RECONSTRUCTION_NO_REFIT',elapsed_s=time.monotonic()-tic,new_optimization_calls=0,main_functions_imported=False))
 print('COLD_AUDIT',fail_count==0,'checks',sum(counts.values()),'fails',fail_count,flush=True)
 if fail_count:raise SystemExit(2)
if __name__=='__main__':main()
