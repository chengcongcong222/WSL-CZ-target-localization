"""Independent full cold reconstruction from original complex archives: QR/eigen and quadratic scores."""
from t2_n1_common import ROOT,OUT,LOCAL,OLD,ARCHIVE,RES,CATS,FREQ,PROBES,PAIRS,TOL,sha,read,rows,dump,table,utc,verify_identities
import numpy as np,time

def unit_independent(h):
 n=np.sqrt(np.einsum('cfm,cfm->cf',h.conj(),h).real)
 if np.any(n<=0) or not np.isfinite(h).all():raise ValueError('cold invalid input')
 return h/n[:,:,None],n

def qr_spectral_error(H,L):
 pair=np.stack([L,H],axis=-1);Q,R=np.linalg.qr(pair,mode='reduced');l=R[...,0];h=R[...,1];A=h[...,None]*h[...,None,:].conj()-l[...,None]*l[...,None,:].conj();ev=np.linalg.eigvalsh(A)
 return np.max(np.abs(ev),axis=-1),np.abs(R[...,1,1])

def quadratic_full_scores(u,probes):
 # observations are supplied explicitly as fixed normalized high-mesh arrays.
 w=probes;allscore=np.empty((len(w),len(u)))
 for begin in range(0,len(u),256):
  v=u[begin:begin+256];projector=v[...,None]*v[...,None,:].conj();q=np.einsum('pfm,cfmn,pfn->pcf',w.conj(),projector,w,optimize=True);allscore[:,begin:begin+len(v)]=q.real.mean(axis=-1)
 return allscore

def independent_order(s,ids=None):
 ids=np.arange(len(s)) if ids is None else ids
 return np.array(sorted(ids.tolist(),key=lambda i:(-float(s[i]),int(i))),int)

def main():
 if (OUT/'COLD_AUDIT_STARTED.json').exists():raise RuntimeError('cold audit already started; no repeat')
 started=time.monotonic();dump(OUT/'COLD_AUDIT_STARTED.json',{'utc':utc(),'original_archives_read_again':True,'main_score_or_normalization_function_imported':False});identities=verify_identities();design=read(OUT/'N1_DESIGN_FREEZE.json');checks=[];rank_notes=[];pair_cache={};maxe=0.;maxE=0.;maxscore=0.;maxbound=-float('inf');vectorentries=0;scoreentries=0;countabove=0;groupfails=0;probes_count=0
 def add(name,e,limit=TOL,detail='',elements=1):checks.append({'check':name,'error':float(e),'limit':limit,'PASS':bool(np.isfinite(e) and e<=limit),'detail':detail,'elements':int(elements)})
 for path,h in design['code_and_input_SHA256'].items():add('design_identity',0 if sha(ROOT/path)==h else 1,0,path)
 seal=read(OUT/'ANALYSIS_FREEZE.json')
 for path,h in seal['data_files'].items():add('analysis_array_identity',0 if sha(path)==h else 1,0,path)
 for name,h in seal['result_files'].items():add('analysis_result_identity',0 if sha(OUT/name)==h else 1,0,name)
 distribution=rows(OUT/'FIELD_SHAPE_DISTRIBUTION.csv');boundsummary=rows(OUT/'SCORE_PERTURBATION_BOUNDS_BY_RESOURCE.csv');pr=rows(OUT/'TEMPLATE_PROBE_RANKING.csv');normrows=rows(OUT/'FIELD_NORM_DECILE_STRESS.csv');strata=rows(OUT/'STRATIFIED_FAILURE_COUNTS.csv');joint=rows(OUT/'JOINT_GEOMETRY_FAILURE_COUNTS.csv');alltop=rows(OUT/'TOP_ERROR_CANDIDATES.csv');cold_max=0.
 identity=read(OUT/'INPUT_IDENTITY_AND_SCOPE.json');paths={r['filename']:r['path'] for r in identity['arrays']};near_tied_argmax_differences=0
 for cat in CATS:
  states=np.load(OLD/f'CATALOGUE_{cat}.npy',mmap_mode='r');n=len(states);ids=np.arange(n)
  for res in RES:
   high=np.load(paths[f'{res}_{cat}_n160001.npy'],mmap_mode='r');low=np.load(paths[f'{res}_{cat}_n80001.npy'],mmap_mode='r');H,nh=unit_independent(high);L,nl=unit_independent(low);e,er=qr_spectral_error(H,L);old=np.load(paths[f'SHAPE_ERROR_{res}_{cat}.npy'],mmap_mode='r');main=np.load(LOCAL/f'FIELDS_{res}_{cat}.npz');E=e.mean(axis=1);de=float(np.max(abs(e-old)));dmain=float(np.max(abs(e-main['e'])));dE=float(np.max(abs(E-main['E'])));maxe=max(maxe,de,dmain);maxE=max(maxE,dE);vectorentries+=e.size
   add('all_archived_e_by_QR_Hermitian_norm',de,TOL,res+'/'+cat,e.size);add('main_e_by_independent_QR',dmain,TOL,res+'/'+cat,e.size);add('QR_residual_vs_spectral_norm',float(np.max(abs(e-er))),TOL,res+'/'+cat,e.size);add('complete_E_from_eigen_norm',dE,TOL,res+'/'+cat,n)
   for fi,f in enumerate(FREQ):
    group=next(r for r in distribution if r['resource']==res and r['catalogue']==cat and int(r['frequency_hz'])==f);above=int(np.sum(e[:,fi]>.002));countabove+=above;groupfails+=int(e[:,fi].max()>.002);cold_max=max(cold_max,float(e[:,fi].max()));add('54_group_denominator',abs(int(group['N'])-n),0);add('54_above_threshold_count',abs(int(group['above_historical_0p002'])-above),0);add('54_error_median',abs(float(np.median(e[:,fi]))-float(group['median'])),TOL)
    order=np.array(sorted(ids.tolist(),key=lambda i:(float(nh[i,fi]),int(i))));labels=np.empty(n,int);labels[order]=(np.arange(n)*10//n)
    for d in range(10):
     r=next(r for r in normrows if r['resource']==res and r['catalogue']==cat and int(r['frequency_hz'])==f and int(r['norm_decile_low_to_high'])==d+1);mask=labels==d;add('norm_decile_N_and_above',abs(int(r['N'])-int(mask.sum()))+abs(int(r['above_historical_0p002'])-int(np.sum(e[mask,fi]>.002))),0)
    def reconstruct_mask(row):
     axis=row['axis'];value=row['stratum']
     if axis=='range':return (states[:,0]<50) if value=='[45,50)' else ((states[:,0]>=50)&(states[:,0]<55)) if value=='[50,55)' else states[:,0]>=55
     if axis=='angle':return (states[:,1]==0) if value=='theta=0' else ((states[:,1]>0)&(states[:,1]<5)) if value=='0<theta<5' else states[:,1]==5
     return states[:,2]==int(value)
    selected=[r for r in strata if r['resource']==res and r['catalogue']==cat and int(r['frequency_hz'])==f]
    for r in selected:
     mask=reconstruct_mask(r);add('full_marginal_stratum_count',abs(int(r['N'])-int(mask.sum()))+abs(int(r['above_historical_0p002'])-int(np.sum(e[mask,fi]>.002))),0)
    selected=[r for r in joint if r['resource']==res and r['catalogue']==cat and int(r['frequency_hz'])==f]
    for r in selected:
     mask=reconstruct_mask({'axis':'range','stratum':r['range_bin']})&reconstruct_mask({'axis':'angle','stratum':r['angle_bin']})&(states[:,2]==int(r['source_depth_m']));add('all_joint_stratum_count',abs(int(r['N'])-int(mask.sum()))+abs(int(r['above_historical_0p002'])-int(np.sum(e[mask,fi]>.002))),0)
    top=sorted([r for r in alltop if r['resource']==res and r['catalogue']==cat and int(r['frequency_hz'])==f],key=lambda r:int(r['rank']))
    for r in top:add('worst5_candidate_error',abs(e[int(r['candidate_id']),fi]-float(r['shape_error'])),TOL)
   pids=[]
   for pid,r,t,z in PROBES:
    match=np.flatnonzero(np.all(states==[r,t,z],axis=1));assert len(match)==1;pids.append(int(match[0]))
   w=H[pids].copy();SH=quadratic_full_scores(H,w);SL=quadratic_full_scores(L,w);stored=np.load(LOCAL/f'PROBE_SCORES_{res}_{cat}.npy');dsh=float(np.max(abs(SH-stored[:,0])));dsl=float(np.max(abs(SL-stored[:,1])));maxscore=max(maxscore,dsh,dsl);scoreentries+=SH.size+SL.size;add('all_high_quadratic_scores',dsh,TOL,res+'/'+cat,SH.size);add('all_low_quadratic_scores',dsl,TOL,res+'/'+cat,SL.size)
   excess=float(np.max(abs(SH-SL)-E[None,:]));maxbound=max(maxbound,excess);add('all_fixed_w_score_bounds',max(0.,excess),TOL,res+'/'+cat,SH.size)
   for pi,(pid,r,t,z) in enumerate(PROBES):
    row=next(v for v in pr if v['resource']==res and v['catalogue']==cat and v['probe']==pid);p=pids[pi];wrong=np.flatnonzero(abs(states[:,0]-r)>=1.);wh=int(independent_order(SH[pi],wrong)[0]);wl=int(independent_order(SL[pi],wrong)[0]);GH=float(SH[pi,p]-SH[pi,wh]);GL=float(SL[pi,p]-SL[pi,wl]);B=float(E[p]+E[wrong].max());add('probe_reference_and_wrong_max_scores',max(abs(SH[pi,p]-float(row['probe_score_H'])),abs(SL[pi,p]-float(row['probe_score_L'])),abs(SH[pi,wh]-float(row['strongest_wrong_score_H'])),abs(SL[pi,wl]-float(row['strongest_wrong_score_L']))),TOL,res+'/'+cat+'/'+pid,4);add('independent_gap_values',max(abs(GH-float(row['G_H'])),abs(GL-float(row['G_L']))),TOL,res+'/'+cat+'/'+pid,2);add('gap_B_independent',abs(B-float(row['B_gap'])),TOL);add('independent_gap_bound',max(0.,abs(GH-GL)-B),TOL)
    for mesh,s,coldid,key in [(160001,SH[pi],wh,'strongest_wrong_ID_H'),(80001,SL[pi],wl,'strongest_wrong_ID_L')]:
     mainid=int(row[key]);near=coldid!=mainid;delta=float(s[coldid]-s[mainid]);add('main_reported_wrong_argmax_cold_score',max(0.,delta),TOL,res+'/'+cat+'/'+pid+'/'+str(mesh));near_tied_argmax_differences+=int(near);rank_notes.append({'resource':res,'catalogue':cat,'probe':pid,'mesh':mesh,'main_wrong_id':mainid,'cold_quadratic_wrong_id':coldid,'ID_matches':not near,'cold_score_margin':delta,'within_declared_roundtrip_tolerance':abs(delta)<=TOL,'main_rank_not_overwritten':True})
    pair_cache[res,cat,pid]={'G_H':GH,'G_L':GL,'B':B};probes_count+=1
   print('COLD',res,cat,'all e/E and5 complete score vectors',flush=True)
 pairrows=rows(OUT/'PAIR_ADVANTAGE_NUMERIC_INTERVALS.csv');cert_count={};flip_count={};classification_notes=[]
 for row in pairrows:
  label=row['pair'];cat=row['catalogue'];pid=row['probe'];t=pair_cache[row['treated'],cat,pid];c=pair_cache[row['control'],cat,pid];aH=t['G_H']-c['G_H'];aL=t['G_L']-c['G_L'];B=t['B']+c['B'];add('20_pair_advantages_and_bounds',max(abs(aH-float(row['Adv_H'])),abs(aL-float(row['Adv_L'])),abs(B-float(row['B_advantage_sum']))),TOL,label+cat+pid,3);add('20_independent_pair_inequalities',max(0.,abs(aH-aL)-B),TOL);robust=abs(aH)>B;reported=row['sign_robust_by_strict_bound']=='True';add('20_strict_bound_certificates',0 if robust==reported else 1,0,label+cat+pid);cert_count[label]=cert_count.get(label,0)+int(robust);flip=bool(np.sign(aH)!=np.sign(aL));reportedflip=row['actual_advantage_sign_flip']=='True'
  if flip!=reportedflip:
   tiny=max(abs(aH),abs(aL),abs(float(row['Adv_H'])),abs(float(row['Adv_L'])))<=TOL;add('independent_advantage_sign_agreement_or_roundoff',0 if tiny else 1,0);classification_notes.append({'pair':label,'catalogue':cat,'probe':pid,'main_flip':reportedflip,'cold_flip':flip,'all_near_roundoff':tiny})
  flip_count[label]=flip_count.get(label,0)+int(flip)
 historical=read(OUT/'HISTORICAL_COUNT_REPRODUCTION.json');add('historical_total_vectors',abs(vectorentries-433620),0);add('historical8012_above',abs(countabove-8012),0);add('historical54_groupfails',abs(groupfails-54),0);add('historical_global_max',abs(cold_max-historical['shape_error_max']),TOL)
 decision=read(OUT/'DECISION.json')
 for label in ['A','B']:add('decision_strict_certificate_count',abs(cert_count[label]-decision['pair_analysis'][label]['strict_probe_certificates']),0)
 # Spectral norm identities and all bounds are independent; scores never call frozen main score code.
 table(OUT/'COLD_AUDIT.csv',checks);table(OUT/'COLD_RANK_ID_ROUNDOFF_DIAGNOSTICS.csv',rank_notes)
 if classification_notes:table(OUT/'COLD_CLASSIFICATION_ROUNDOFF_NOTES.csv',classification_notes)
 result={'PASS':all(r['PASS'] for r in checks),'check_rows':len(checks),'FAIL':sum(not r['PASS'] for r in checks),'source_identities':identities,'full_shape_error_entries':vectorentries,'full_E_candidate_entries':vectorentries//3,'full_probe_score_entries_both_models':scoreentries,'probe_groups':probes_count,'pair_advantage_groups':len(pairrows),'max_shape_error_absolute_discrepancy':maxe,'max_E_absolute_discrepancy':maxE,'max_score_absolute_discrepancy':maxscore,'max_signed_fixed_w_bound_excess':maxbound,'historical_above8012_reproduced':countabove,'historical_failed_max_groups_reproduced':groupfails,'independent_wrong_argmax_ID_differences':near_tied_argmax_differences,'independent_pair_certificate_counts':cert_count,'independent_pair_advantage_sign_flips':flip_count,'near_roundoff_classification_notes':len(classification_notes),'original_inputs_read_again':True,'main_score_function_used':False,'source_mutations_or_new_physical_inputs':0,'wall_s':time.monotonic()-started,'tolerances':{'e_score_inequalities':TOL},'scope':'two_archived_meshes_not_physical_accuracy_or_receiver_performance','R4_percent':0}
 dump(OUT/'VALIDATION.json',result);print('COLD COMPLETE',result,flush=True)
 if not result['PASS']:raise RuntimeError('INDEPENDENT_AUDIT_FAILURE_STOP')
if __name__=='__main__':main()