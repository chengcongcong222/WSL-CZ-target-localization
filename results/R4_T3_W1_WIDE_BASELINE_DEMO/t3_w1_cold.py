"""Independent T3-W1 cold audit. No generator, fitter, information or statistics imports.
Rebuilds observations, Cartesian Fisher, scalar fit likelihood, all quantiles and pairing.
No optimizer is run during cold verification.
"""
from pathlib import Path
import json,csv,hashlib,math,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent;LOCAL=Path('D:/ProjectStorage/WSL-CZ/T3_W1')
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def yes(v):return str(v).lower()=='true'
def quant(values,q):
    data=sorted(float(v) for v in values);return data[int(math.ceil(q*len(data)))-1] if data else None

def main():
    cfg=read(OUT/'T3_W1_DESIGN_FREEZE.json');checks=[];count_raw=0;public_rebuilt=0;metrics={};raw_groups={};public_hashes={};cost_errors=[];max_obs_error=0.;info_errors=[]
    def ck(name,error,tol,detail=''):
        checks.append({'check':name,'error':float(error),'tolerance':float(tol),'PASS':bool(np.isfinite(error) and error<=tol),'detail':str(detail)})
    def eq(name,a,b,tol=1e-9,detail=''):
        if a is None or b is None or b=='':ck(name,0 if a is None and (b is None or b=='') else 1,0,detail);return
        a,b=float(a),float(b)
        if not math.isfinite(a) or not math.isfinite(b):ck(name,0 if a==b else 1,0,detail)
        else:ck(name,abs(a-b),tol,detail)
    try:
        for r in cfg['inputs']:ck('source_identity',0 if sha(ROOT/r['path'])==r['sha256'] else 1,0,r['path'])
        for r in cfg['frozen_code']:ck('frozen_code_identity',0 if sha(OUT/r['path'])==r['sha256'] else 1,0,r['path'])
        if any(not r['PASS'] for r in checks):raise RuntimeError('COLD_IDENTITY_STOP')
        panel=rows(ROOT/cfg['truth_panel']);truth={r['case_id']:r for r in panel};public=rows(OUT/'PUBLIC_BLOCK_IDENTITIES.csv');private=rows(OUT/'PRIVATE_INNOVATION_IDENTITIES.csv');raw=rows(OUT/'ESTIMATOR_OUTPUTS.csv');error_rows=rows(OUT/'RUN_ERRORS.csv');rawmap={(r['config'],r['case_id'],int(r['replicate'])):r for r in raw};emap={(r['config'],r['case_id'],int(r['replicate'])):r for r in error_rows};names=['range_rel','bearing_deg','speed_abs_mps','speed_rel','heading_deg'];qs=[('median',.5),('P90',.9),('P95',.95),('P99',.99)]
        t=np.arange(0.,1201.,10.);mainp=np.empty((121,2));mainp[:,0]=2*np.minimum(t,600)+2*np.maximum(t-600,0)*math.cos(math.pi/12);mainp[:,1]=2*np.maximum(t-600,0)*math.sin(math.pi/12)
        for r in private:ck('private_block_identity',0 if sha(r['path'])==r['sha256'] else 1,0,r['case_id'])
        for block in public:
            ck('public_block_identity',0 if sha(block['path'])==block['sha256'] else 1,0,block['path']);quality=cfg['qualities'][block['quality']];scheme=cfg['geometries'][block['geometry']];config,case=block['config'],block['case_id'];tr=truth[case];htrue=np.array([float(tr[x]) for x in ['r0_km','theta0_deg','v_mps','psi_deg']]);g=int(case[1:])-1
            with np.load(block['path'],allow_pickle=False) as z:
                P=z['reported_positions'];y=z['bearings'];sides=z['sides'];tt=z['times'];ck('public_shapes',0 if P.shape==(500,121,2,2) and y.shape==(500,121,2) and np.isfinite(P).all() and np.isfinite(y).all() and np.array_equal(tt,t) else 1,0,config+'/'+case);ck('balanced_mirrors',0 if np.sum(sides==1)==250 and np.sum(sides==-1)==250 else 1,0,config+'/'+case)
                mh=hashlib.sha256(t.tobytes()+np.ascontiguousarray(P[:,:,0]).tobytes()+np.ascontiguousarray(y[:,:,0]).tobytes()).hexdigest();ck('main_hash_manifest',0 if mh==block['MAIN_byte_hash'] else 1,0,config+'/'+case)
                key=block['quality'],case
                if key in public_hashes:ck('paired_MAIN_exact_bytes',0 if mh==public_hashes[key] else 1,0,config+'/'+case)
                else:public_hashes[key]=mh
                # Full H01 input reconstruction: all500repetitions, all8resource-quality configurations.
                if case=='H01':
                    record=next(r for r in private if r['case_id']==case)
                    with np.load(record['path']) as pr:
                        innovations=pr['normals'];uniform=pr['uniforms'];rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([cfg['seed'],0])));ru=rng.random((500,2));rn=rng.standard_normal((500,121,6));ck('saved_uniform_seed_identity',float(np.max(abs(ru-uniform))),0,config);ck('saved_normal_seed_identity',float(np.max(abs(rn-innovations))),0,config)
                    r,angle,v,heading=htrue;angle,heading=np.deg2rad([angle,heading]);position=np.array([1000*r*math.cos(angle),1000*r*math.sin(angle)]);velocity=np.array([v*math.cos(heading),v*math.sin(heading)]);target=position+t[:,None]*velocity
                    true_nodes=np.repeat(mainp[None,:,None,:],500,axis=0);true_nodes=np.repeat(true_nodes,2,axis=2);true_nodes[:,:,1,0]+=scheme['forward_km']*1000;true_nodes[:,:,1,1]+=sides[:,None]*scheme['lateral_km']*1000
                    pp=true_nodes+quality['nav_sigma_m']*innovations[:,:,2:].reshape(500,121,2,2);c=(2*uniform[:,0]-1)*quality['common_bound_deg'];d=(2*uniform[:,1]-1)*quality['half_diff_bound_deg'];bias=np.deg2rad(np.column_stack([c-d,c+d]));delta=target[None,:,None,:]-true_nodes;yy=np.arctan2(delta[:,:,:,1],delta[:,:,:,0])+bias[:,None,:]+np.deg2rad(quality['bearing_sigma_deg'])*innovations[:,:,:2];yy=(yy+np.pi)%(2*np.pi)-np.pi
                    pe=float(np.max(abs(P-pp)));ae=float(np.max(abs(np.arctan2(np.sin(y-yy),np.cos(y-yy)))));ck('H01_all_public_positions_rebuilt',pe,1e-9,config);ck('H01_all_public_bearings_rebuilt',ae,1e-12,config);max_obs_error=max(max_obs_error,ae);public_rebuilt+=500
                for rep in range(500):
                    r=rawmap[config,case,rep];accepted=yes(r['accepted']);h=np.array([float(r[k]) for k in ['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg','common_hat_deg','half_diff_hat_deg']]);diff=h[:4]-htrue
                    if accepted:
                        angle_errors=np.abs(np.arctan2(np.sin(np.deg2rad(diff[[1,3]])),np.cos(np.deg2rad(diff[[1,3]]))))*180/math.pi
                        errors=[abs(diff[0])/htrue[0],angle_errors[0],abs(diff[2]),abs(diff[2])/htrue[2],angle_errors[1]]
                        a,p=math.radians(h[1]),math.radians(h[3]);target=P[rep,0,0]+1000*h[0]*np.array([math.cos(a),math.sin(a)])+t[:,None]*h[2]*np.array([math.cos(p),math.sin(p)]);dx=target[:,None,:]-P[rep];rho2=np.sum(dx*dx,axis=-1);mu=np.arctan2(dx[:,:,1],dx[:,:,0])+np.deg2rad(h[4]+h[5]*np.array([-1.,1.]));res=np.arctan2(np.sin(mu-y[rep]),np.cos(mu-y[rep]));base=math.radians(quality['bearing_sigma_deg'])**2;var=base+quality['nav_sigma_m']**2/rho2;nll=.5*math.fsum((res*res/var+np.log(var/base)).ravel());cost_errors.append(abs(nll-float(r['cost'])))
                    else:errors=[math.inf]*5
                    key=config,case,rep;metrics[key]=errors;raw_groups.setdefault((config,case),[]).append(r);rr=emap[key]
                    if any((math.isfinite(e) and abs(e-float(rr[m]))>1e-9) or (not math.isfinite(e) and float(rr[m])!=math.inf) for m,e in zip(names,errors)):raise RuntimeError('RAW_METRIC_OR_FAILURE_INF_MISMATCH '+str(key))
                    if int(r['side'])!=int(sides[rep]) or int(r['draw_index'])!=g*500+rep:raise RuntimeError('PAIR_OR_DRAW_ID_MISMATCH')
                    count_raw+=1
        ck('all_metric_records_rebuilt',abs(count_raw-48000),0);ck('max_likelihood_cost_rebuild_error',max(cost_errors,default=0),1e-6)
        case_summary=rows(OUT/'CASE_STATISTICS.csv')
        for summary in case_summary:
            config,case,side=summary['config'],summary['case_id'],summary['side_group'];rr=raw_groups[config,case];rr=[r for r in rr if side=='BOTH' or r['side']==side];good=[r for r in rr if yes(r['accepted'])];eq('case_n',len(rr),summary['n'],0);eq('case_accepted',len(good),summary['n_accepted'],0);eq('case_failure_rate',1-len(good)/len(rr),summary['failure_rate'],1e-14)
            for field,source in [('state_boundary_rate','state_at_bound'),('bias_boundary_rate','bias_at_bound'),('common_boundary_rate','common_at_bound'),('diff_boundary_rate','diff_at_bound')]:eq('boundary_rate/'+field,sum(yes(r[source]) for r in rr)/len(rr),summary[field],1e-14)
            for tag,sel in [('unconditional',rr),('accepted_only',good)]:
                for j,m in enumerate(names):
                    for quantile,q in qs:eq('case_quantile/'+m+'/'+tag+'/'+quantile,quant([metrics[r['config'],r['case_id'],int(r['replicate'])][j] for r in sel],q),summary[f'{m}_{tag}_{quantile}'],1e-9,config+'/'+case+'/'+side)
        globalrows=rows(OUT/'GLOBAL_STATISTICS.csv')
        for summary in globalrows:
            rr=[r for r in raw if r['config']==summary['config']];good=[r for r in rr if yes(r['accepted'])];eq('global_n',len(rr),summary['n'],0);eq('global_acceptance',len(good),summary['accepted'],0);eq('global_failure_rate',1-len(good)/len(rr),summary['failure_rate'],1e-14)
            for tag,sel in [('unconditional',rr),('accepted_only',good)]:
                for j,m in enumerate(names):
                    for qt,q in qs:eq('global_quantile/'+m+'/'+tag+'/'+qt,quant([metrics[r['config'],r['case_id'],int(r['replicate'])][j] for r in sel],q),summary[f'{m}_{tag}_{qt}'],1e-9,summary['config'])
            case_rows=[r for r in case_summary if r['config']==summary['config'] and r['side_group']=='BOTH']
            for m in names:
                for qt,q in qs:eq('worst_case_quantile',max(float(r[f'{m}_unconditional_{qt}']) for r in case_rows),summary[f'{m}_worst_case_unconditional_{qt}'],1e-12,summary['config']+'/'+m+'/'+qt)
        for r in rows(OUT/'PAIRED_CHANGES.csv'):
            j=names.index(r['metric']);a=np.array([metrics[r['quality']+'_B0',r['case_id'],i][j] for i in range(500)]);b=np.array([metrics[r['quality']+'_'+r['geometry'],r['case_id'],i][j] for i in range(500)]);finite=np.isfinite(a)&np.isfinite(b)
            for field,val in [('better',np.sum(b<a)),('worse',np.sum(b>a)),('equal',np.sum(b==a)),('both_finite',np.sum(finite)),('baseline_failed_candidate_valid',np.sum(~np.isfinite(a)&np.isfinite(b))),('candidate_failed_baseline_valid',np.sum(np.isfinite(a)&~np.isfinite(b))),('both_failed',np.sum(~np.isfinite(a)&~np.isfinite(b)))]:eq('paired_'+field,val,r[field],0,r['quality']+'/'+r['geometry']+'/'+r['case_id']+'/'+r['metric'])
            meda=quant(a,.5);medb=quant(b,.5);eq('paired_median_ratio',medb/meda if meda>0 and np.isfinite(meda) else None,r['median_error_ratio'],1e-9)
        # Fisher rebuilt using independent scalar row construction for every nominal geometry.
        info=rows(OUT/'LOCAL_INFORMATION.csv');matrix=np.load(OUT/'INFORMATION_MATRICES.npz');scale=np.array(cfg['information_scale']);Fmax=0.;Emax=0.;FDmax=0.
        for row in [r for r in info if r['control']=='BIAS_FREE_NO_PRIOR_LOCAL_DATA_INFORMATION']:
            quality=cfg['qualities'][row['quality']];scheme=cfg['geometries'][row['geometry']];tr=truth[row['case_id']];r,th,v,ps=[float(tr[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']];th,ps=map(math.radians,[th,ps]);state=np.array([1000*r*math.cos(th),1000*r*math.sin(th),v*math.cos(ps),v*math.sin(ps),0,0]);side=int(row['side']);P=np.stack([mainp,mainp+np.array([1000*scheme['forward_km'],1000*side*scheme['lateral_km']])],axis=1);H=[];D=[];variances=[]
            for epoch,tj in enumerate(t):
                target=state[:2]+tj*state[2:4]
                for node in [0,1]:
                    dx,dy=target-P[epoch,node];rho=dx*dx+dy*dy;meanrow=np.array([-dy/rho,dx/rho,-tj*dy/rho,tj*dx/rho,1,-1 if node==0 else 1]);deriv=-2*quality['nav_sigma_m']**2/rho**2*np.array([dx,dy,tj*dx,tj*dy,0,0]);vv=math.radians(quality['bearing_sigma_deg'])**2+quality['nav_sigma_m']**2/rho;H.append(meanrow);D.append(deriv);variances.append(vv)
            H=np.array(H);D=np.array(D);var=np.array(variances);A=np.concatenate([H/np.sqrt(var[:,None]),D/(np.sqrt(2)*var[:,None])])*scale;F=A.T@A;index=int(row['matrix_index']);targetF=matrix['F_data'][index]*scale[:,None]*scale[None,:];Fmax=max(Fmax,float(np.linalg.norm(F-targetF)/np.linalg.norm(F)))
            N=A[:,[0,1,4,5]];V=A[:,[2,3]];projection=np.linalg.lstsq(N,V,rcond=1e-10)[0];res=V-N@projection;E=res.T@res;Emax=max(Emax,float(np.linalg.norm(E-matrix['effective_velocity_scaled'][index])/np.linalg.norm(E)))
            # Independent finite difference of atan2 means in all6 scaled directions.
            for j,sc in enumerate(scale):
                step=sc*1e-5;a=state.copy();b=state.copy();a[j]+=step;b[j]-=step
                def bearing(s):
                    delta=(s[:2]+t[:,None]*s[2:4])[:,None,:]-P;return (np.arctan2(delta[:,:,1],delta[:,:,0])+s[4]+s[5]*np.array([-1,1])).ravel()
                da=bearing(a)-bearing(b);numer=np.arctan2(np.sin(da),np.cos(da))/(2*step);FDmax=max(FDmax,float(np.max(abs((numer-H[:,j])*sc))))
            for known in [False,True]:
                aa=A[:,:4] if known else A;_,ss,Vh=np.linalg.svd(aa,full_matrices=False);rank=int(np.sum(ss/ss[0]>1e-10));rr=next(vv for vv in info if vv['quality']==row['quality'] and vv['geometry']==row['geometry'] and vv['case_id']==row['case_id'] and vv['side']==row['side'] and vv['control']==('BIAS_KNOWN_INFORMATION_UPPER_CONTROL' if known else 'BIAS_FREE_NO_PRIOR_LOCAL_DATA_INFORMATION'));eq('Fisher_rank',rank,rr['rank'],0)
                if rank==aa.shape[1]:
                    sc=scale[:4] if known else scale;C=((Vh.T/ss**2)@Vh)*sc[:,None]*sc[None,:];
                    radius=math.hypot(*state[:2]);tail=[] if known else [0,0];gradients=[np.r_[state[0]/radius,state[1]/radius,0,0,tail],np.r_[-state[1]/radius**2,state[0]/radius**2,0,0,tail],np.r_[0,0,state[2]/v,state[3]/v,tail],np.r_[0,0,-state[3]/v**2,state[2]/v**2,tail]]
                    for gr,field,unit in zip(gradients,['sigma_range_m','sigma_bearing_deg','sigma_speed_mps','sigma_heading_deg'],[1,180/math.pi,1,180/math.pi]):
                        sigma=math.sqrt(max(0,float(gr@C@gr)))*unit;eq('profiled_uncertainty/'+field,sigma,rr[field],1e-6*max(1,sigma))
                else:eq('rankdef_speed_is_infinite',math.inf,rr['sigma_speed_mps'],0)
        ck('all_Fisher_matrices_relative_error',Fmax,1e-8);ck('all_velocity_profile_matrices_relative_error',Emax,1e-5);ck('independent_analytic_FD_Jacobian',FDmax,1e-7);matrix.close()
        # Synthetic failure bookkeeping control: rejected rows must affect unconditional quantiles.
        eq('failure_INF_nearest_rank_control',quant([1.,2.,math.inf],.95),math.inf,0);eq('circular_boundary_control',abs(math.degrees(math.atan2(math.sin(math.radians(359.8)),math.cos(math.radians(359.8))))),.2,1e-10)
        result={'PASS':all(r['PASS'] for r in checks),'checks':len(checks),'FAIL':sum(not r['PASS'] for r in checks),'raw_configurations_recomputed':count_raw,'public_H01_full_records_rebuilt':public_rebuilt,'all_public_block_identities':len(public),'quantiles':'allcase/bothside/pooled/worstcaseconditionalandunconditional;nearest-rank','paired_comparisons':'allB1/B2/B3vsB0qualitiescasesmetrics','Fisher_geometries_recomputed':192,'new_optimizations':0,'R4_percent':0,'STOP':True}
    except Exception as exc:
        result={'PASS':False,'checks':len(checks),'FAIL':sum(not r['PASS'] for r in checks)+1,'exception':str(exc),'raw_configurations_recomputed':count_raw,'R4_percent':0,'STOP':True};(OUT/'COLD_EXCEPTION.txt').write_text(traceback.format_exc(),encoding='utf-8')
    with (OUT/'COLD_CHECKS.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['check','error','tolerance','PASS','detail'],lineterminator='\n');writer.writeheader();writer.writerows(checks)
    (OUT/'INDEPENDENT_COLD_AUDIT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
