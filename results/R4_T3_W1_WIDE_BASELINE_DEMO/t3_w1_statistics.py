"""Registered error and paired-statistic rules. Postprocessing does not refit."""
from pathlib import Path
import csv,math,json
import numpy as np
METRICS=['range_rel','bearing_deg','speed_abs_mps','speed_rel','heading_deg']
QUANTILES=[('median',.5),('P90',.9),('P95',.95),('P99',.99)]

def csvread(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def csvwrite(p,rr,columns=None):
    fields=columns or list(dict.fromkeys(k for r in rr for k in r))
    with Path(p).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');w.writeheader();w.writerows(rr)
def rank_quantile(values,q):
    a=np.sort(np.asarray(values,float));return float(a[math.ceil(len(a)*q)-1]) if len(a) else None

def evaluate(records,panel):
    truths={r['case_id']:r for r in panel};data=[]
    for rr in records:
        tr=truths[rr['case_id']];ok=str(rr['accepted']).lower()=='true';h=np.array([float(rr[k]) for k in ['r_hat_km','theta_hat_deg','v_hat_mps','psi_hat_deg']]);truth=np.array([float(tr[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']])
        delta=h-truth;angular=np.abs((delta[[1,3]]+180)%360-180)
        error=[abs(delta[0])/truth[0],angular[0],abs(delta[2]),abs(delta[2])/truth[2],angular[1]] if ok else [math.inf]*5
        data.append({k:rr[k] for k in ['config','quality','geometry','case_id','replicate','side','draw_index','accepted']}|dict(zip(METRICS,error)))
    return data

def summarize(errors,raw):
    group={};diagnostics={};summary=[]
    for r in errors:
        for side in ['BOTH',str(r['side'])]:group.setdefault((r['config'],r['case_id'],side),[]).append(r)
    for r in raw:
        for side in ['BOTH',str(r['side'])]:diagnostics.setdefault((r['config'],r['case_id'],side),[]).append(r)
    for key,rr in sorted(group.items()):
        config,case,side=key;rg=diagnostics[key];valid=[r for r in rr if str(r['accepted']).lower()=='true'];n=len(rr)
        info={'config':config,'quality':rr[0]['quality'],'geometry':rr[0]['geometry'],'case_id':case,'side_group':side,'n':n,'n_accepted':len(valid),'n_failed':n-len(valid),'failure_rate':1-len(valid)/n,'state_boundary_rate':sum(str(r['state_at_bound']).lower()=='true' for r in rg)/n,'bias_boundary_rate':sum(str(r['bias_at_bound']).lower()=='true' for r in rg)/n,'common_boundary_rate':sum(str(r['common_at_bound']).lower()=='true' for r in rg)/n,'diff_boundary_rate':sum(str(r['diff_at_bound']).lower()=='true' for r in rg)/n,'weak_fit_rank_rate':sum(int(r['rank_augmented_likelihood_jac'])<6 for r in rg)/n}
        for mode,sel in [('unconditional',rr),('accepted_only',valid)]:
            for m in METRICS:
                for name,q in QUANTILES:info[f'{m}_{mode}_{name}']=rank_quantile([float(r[m]) for r in sel],q)
        summary.append(info)
    global_rows=[]
    for config in sorted({r['config'] for r in errors}):
        rr=[r for r in errors if r['config']==config];cases=[r for r in summary if r['config']==config and r['side_group']=='BOTH'];valid=[r for r in rr if str(r['accepted']).lower()=='true'];r={'config':config,'quality':rr[0]['quality'],'geometry':rr[0]['geometry'],'n':len(rr),'accepted':len(valid),'failure_rate':1-len(valid)/len(rr),'case_weighting':'equal500percase;pooledsameweight;P95stillfinitepaneldescriptive'}
        for mode,sel in [('unconditional',rr),('accepted_only',valid)]:
            for m in METRICS:
                for name,q in QUANTILES:r[f'{m}_{mode}_{name}']=rank_quantile([float(v[m]) for v in sel],q)
        for m in METRICS:
            for name,q in QUANTILES:r[f'{m}_worst_case_unconditional_{name}']=max(float(c[f'{m}_unconditional_{name}']) for c in cases)
        global_rows.append(r)
    paired=[];lookup={(r['config'],r['case_id'],int(r['replicate'])):r for r in errors}
    for quality in ['Q1','Q2']:
        for geometry in ['B1','B2','B3']:
            for case in sorted({r['case_id'] for r in errors}):
                for m in METRICS:
                    a=np.array([float(lookup[quality+'_B0',case,i][m]) for i in range(500)]);b=np.array([float(lookup[quality+'_'+geometry,case,i][m]) for i in range(500)]);both=np.isfinite(a)&np.isfinite(b)
                    paired.append({'quality':quality,'geometry':geometry,'case_id':case,'metric':m,'paired_n':500,'better':int(np.sum(b<a)),'worse':int(np.sum(b>a)),'equal':int(np.sum(b==a)),'both_finite':int(np.sum(both)),'baseline_failed_candidate_valid':int(np.sum(~np.isfinite(a)&np.isfinite(b))),'candidate_failed_baseline_valid':int(np.sum(np.isfinite(a)&~np.isfinite(b))),'both_failed':int(np.sum(~np.isfinite(a)&~np.isfinite(b))),'median_error_ratio':rank_quantile(b,.5)/rank_quantile(a,.5) if rank_quantile(a,.5)>0 and np.isfinite(rank_quantile(a,.5)) else None})
    return summary,global_rows,paired

def decision(summary,globals,information):
    rows=[]
    for quality in ['Q1','Q2']:
        for case in sorted({r['case_id'] for r in summary}):
            b=next(r for r in summary if r['config']==quality+'_B0' and r['case_id']==case and r['side_group']=='BOTH');c=next(r for r in summary if r['config']==quality+'_B3' and r['case_id']==case and r['side_group']=='BOTH')
            rows.append({'quality':quality,'case_id':case,'median_speed_improved':c['speed_rel_unconditional_median']<b['speed_rel_unconditional_median'],'P95_speed_improved':c['speed_rel_unconditional_P95']<b['speed_rel_unconditional_P95'],'P95_range_improved':c['range_rel_unconditional_P95']<b['range_rel_unconditional_P95'],'failure_rate_increase':c['failure_rate']-b['failure_rate']})
    return {'stage':'R4-T3-W1','scope':'CONDITIONAL_BEARING_NAVIGATION_FEATURE_LEVEL_DEMO','case_directions':rows,'summary':'WIDE_BASELINE_RESOURCE_INFORMATION_DEMO_ONLY','measurement_chain':'NOT_VALIDATED;directedbearing,sync,association,hardwareUNKNOWN','depth':'NOT_EVALUATED','physical_received_signal_evidence':False,'novelty':'baselineincreasealoneisstandardtriangulation;weakdirectionandbias/resourceconditionsaretheresearchquestion','R4_percent':0,'historical_Gates':'UNCHANGED','next_stage_authorized':False,'STOP':True}
