"""Independent saved-row audit, figures and B1 closeout. No forward modules imported."""
from pathlib import Path
from collections import defaultdict
import ast
import csv
import hashlib
import json
import math
import subprocess
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def js(name):return json.loads((OUT/name).read_text(encoding='utf-8-sig'))
def read_path(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def read(name):return read_path(OUT/name)
def write(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def istrue(x):return x is True or x=='True'

def audit():
    policy=js('B1B_DIAGNOSTIC_POLICY.json');decision=js('B1B_DECISION.json');tol=policy['numerical_tolerance_db'];checks=[]
    def check(name,ok):
        checks.append({'check':name,'pass':bool(ok)})
        if not ok:raise RuntimeError('Saved-data audit failed: '+name)
    original=read_path(ROOT/policy['depth_score_path']);source={}
    for line,row in enumerate(original,start=2):
        key=(row['case_id'],float(row['delta_r_km']),float(row['delta_v_mps']),int(row['z_label_m']))
        check('source_unique:'+str(key),key not in source);source[key]=(float(row['J_z_db']),line,row['J_z_db'])
    cases={c['case_id']:c for c in policy['cases']};labels=policy['depth_labels_m']
    main=read('B1B_PROFILED_DEPTH_CURVES.csv');secondary=read('B1B_MIN_NO_CENTER_DEPTH_CURVES.csv');rows=main+secondary
    check('curve_row_counts',len(main)==6048 and len(secondary)==2016)
    expected={(cid,R,V,a,z) for cid in cases for R in policy['range_rectangle_limits_km'] for V in policy['speed_rectangle_limits_mps'] for a in policy['aggregations'] for z in labels}
    curve_lookup={(r['case_id'],float(r['R_km']),float(r['V_mps']),r['aggregation'],int(r['z_label_m'])):r for r in rows}
    check('all_unique_curve_keys',len(curve_lookup)==8064 and set(curve_lookup)==expected)
    reconstructed={};groups=defaultdict(list)
    for key,row in curve_lookup.items():
        cid,R,V,a,z=key
        eligible=[(value[0],r,v,value[1],value[2]) for (c,r,v,depth),value in source.items() if c==cid and depth==z and abs(r)<=R and abs(v)<=V and (a!='MIN_NO_CENTER' or (r,v)!=(0.,0.))]
        eligible.sort()
        index=0 if a in ('MIN','MIN_NO_CENTER') else (len(eligible)+3)//4-1 if a=='Q25' else (len(eligible)+1)//2-1
        value,r,v,line,text=eligible[index]
        check('aggregate_from_saved_row:'+str(key),float(row['J_profiled_db'])==value and int(row['source_csv_line'])==line and row['source_score_text']==text and int(row['n_rv_nodes'])==len(eligible) and int(row['selected_order_statistic_ordinal'])==index+1)
        check('argmin_semantics:'+str(key),(float(row['argmin_dr_km'])==r and float(row['argmin_dv_mps'])==v) if a in ('MIN','MIN_NO_CENTER') else row['argmin_dr_km']==row['argmin_dv_mps']=='')
        reconstructed[key]=(value,r,v,line,text,eligible);groups[key[:4]].append((z,value))
    provenance=read('B1B_SCORE_SOURCE_PROVENANCE.csv');check('all_provenance_rows',len(provenance)==8064)
    seen=set()
    for row in provenance:
        key=row['case_id'],float(row['R_km']),float(row['V_mps']),row['aggregation'],int(row['z_label_m']);value,r,v,line,text,_=reconstructed[key]
        check('trace:'+str(key),key not in seen and int(row['source_csv_line'])==line and row['source_score_text']==text and float(row['J_profiled_db'])==value and float(row['source_delta_r_km'])==r and float(row['source_delta_v_mps'])==v and row['source_csv_sha256']==policy['bindings'][policy['depth_score_path']] and row['source_csv_path']==policy['depth_score_path'] and istrue(row['profile_value_is_existing_source_value']))
        seen.add(key)
    metric_rows=read('B1B_PROFILE_METRICS.csv');metric_lookup={(m['case_id'],float(m['R_km']),float(m['V_mps']),m['aggregation']):m for m in metric_rows}
    check('all_metric_keys',len(metric_rows)==len(metric_lookup)==384 and set(metric_lookup)==set(groups))
    def independent_metrics(z,j,truth):
        i=z.index(truth);b=min(range(21),key=lambda k:j[k]);low=j[b];span=max(j)-low
        width=None;comp=[]
        if span>tol:
            indices={k for k in range(21) if (j[k]-low)/span<=.5};lo=hi=b
            while lo-1 in indices:lo-=1
            while hi+1 in indices:hi+=1
            comp=z[lo:hi+1];width=z[hi]-z[lo]
        rank=1+sum(x<j[i] for x in j)
        return {'z_hat_grid_m':z[b],'true_depth_rank':rank,'tolerance_aware_true_depth_rank':1+sum(x<j[i]-tol for x in j),'J_true_db':j[i],'J_min_db':low,
                'min_tie_count':sum(x<=low+tol for x in j),'DeltaJ_second_db':min(j[k] for k in range(21) if k!=i)-j[i],
                'DeltaJ_neighbor_db':min(j[i-1],j[i+1])-j[i],'K_z_db':j[i-1]-2*j[i]+j[i+1],
                'strict_local_minima_count':sum(j[k]<j[k-1] and j[k]<j[k+1] for k in range(1,20)),
                'boundary_local_minima_count':int(j[0]<j[1])+int(j[20]<j[19]),'dynamic_range_db':span,'W50_z_m':width,'W50_component_labels':comp,
                'absolute_depth_error_m':abs(z[b]-truth),'practical_two_bin_pass':abs(z[b]-truth)<=10 and rank<=3}
    rebuilt={}
    for key,group in groups.items():
        group.sort();z=[x[0] for x in group];j=[x[1] for x in group];check('21_labels:'+str(key),z==labels)
        m=independent_metrics(z,j,cases[key[0]]['z_true_m']);rebuilt[key]=m;stored=metric_lookup[key]
        check('metric_scope:'+str(key),stored['metric_scope']=='PROFILED_RV_NUISANCE_DEPTH_METRICS')
        for name,value in m.items():
            valid=json.loads(stored[name])==value if isinstance(value,list) else istrue(stored[name])==value if isinstance(value,bool) else stored[name]=='' if value is None else math.isclose(float(stored[name]),value,rel_tol=1e-12,abs_tol=1e-12)
            check('metric:'+str(key)+'/'+name,valid)
    # Independent compensation maps and discrete node-jump diagnostics.
    maps=read('B1B_RV_COMPENSATION_MAP.csv')+read('B1B_MIN_NO_CENTER_RV_COMPENSATION_MAP.csv');check('compensation_map_rows',len(maps)==4032)
    maplookup={(r['case_id'],float(r['R_km']),float(r['V_mps']),r['aggregation'],int(r['z_label_m'])):r for r in maps}
    check('unique_map_keys',len(maplookup)==4032 and set(maplookup)=={key for key in reconstructed if key[3] in ('MIN','MIN_NO_CENTER')})
    for key,row in maplookup.items():
        cid,R,V,a,z=key;value,r,v,line,text,eligible=reconstructed[key];h0=source[cid,0.,0.,z][0];truth=cases[cid]['z_true_m'];i=labels.index(z)
        previous=None if i==0 else reconstructed[cid,R,V,a,labels[i-1]][1:3]
        ties={(rr,vv) for j,rr,vv,_,_ in eligible if j<=value+tol};savedties={(r['delta_r_km'],r['delta_v_mps']) for r in json.loads(row['numerically_equivalent_min_nodes'])}
        check('map:'+str(key),float(row['J_profiled_db'])==value and float(row['argmin_dr_km'])==r and float(row['argmin_dv_mps'])==v and float(row['H0_same_depth_J_db'])==h0 and math.isclose(float(row['improvement_vs_H0_at_same_depth_db']),h0-value,abs_tol=1e-12) and int(row['n_numerically_equivalent_min_nodes'])==len(ties) and savedties==ties)
        checks_bool={'is_wrong_depth':z!=truth,'selected_noncenter_node':(r,v)!=(0.,0.),'wrong_depth_compensation_supported':z!=truth and (r,v)!=(0.,0.) and h0-value>tol,
                     'argmin_node_changed_from_previous_depth':previous is not None and previous!=(r,v),'argmin_range_changed_from_previous_depth':previous is not None and previous[0]!=r,'argmin_speed_changed_from_previous_depth':previous is not None and previous[1]!=v}
        for name,value in checks_bool.items():check('map_flag:'+str(key)+'/'+name,istrue(row[name])==value)
    summaries=read('B1B_COMPENSATION_SUMMARY.csv');check('compensation_summaries',len(summaries)==96)
    for row in summaries:
        cid,R,V=row['case_id'],float(row['R_km']),float(row['V_mps']);truth=cases[cid]['z_true_m'];h0=independent_metrics(labels,[source[cid,0.,0.,z][0] for z in labels],truth);m=rebuilt[cid,R,V,'MIN']
        wrong=min((z for z in labels if z!=truth),key=lambda z:(reconstructed[cid,R,V,'MIN',z][0],z));j,r,v,*_=reconstructed[cid,R,V,'MIN',wrong]
        scope=[maplookup[cid,R,V,'MIN',z] for z in labels];tn=reconstructed[cid,R,V,'MIN',truth][1:3]
        values={'H0_z_hat_grid_m':h0['z_hat_grid_m'],'H0_DeltaJ_second_db':h0['DeltaJ_second_db'],'MIN_z_hat_grid_m':m['z_hat_grid_m'],'MIN_true_depth_rank':m['true_depth_rank'],'MIN_DeltaJ_second_db':m['DeltaJ_second_db'],'margin_reduction_vs_H0_db':h0['DeltaJ_second_db']-m['DeltaJ_second_db'],'best_wrong_depth_label_m':wrong,'best_wrong_depth_min_J_db':j,'best_wrong_depth_dr_km':r,'best_wrong_depth_dv_mps':v,'true_depth_argmin_is_center':tn==(0.,0.),'wrong_depth_compensation_label_count':sum(istrue(s['wrong_depth_compensation_supported']) for s in scope),'adjacent_depth_node_jumps':sum(istrue(s['argmin_node_changed_from_previous_depth']) for s in scope),'adjacent_depth_range_jumps':sum(istrue(s['argmin_range_changed_from_previous_depth']) for s in scope),'adjacent_depth_speed_jumps':sum(istrue(s['argmin_speed_changed_from_previous_depth']) for s in scope),'n_depths_equivalent_to_global_min_within_numerical_tolerance':m['min_tie_count']}
        for name,value in values.items():check('comp_summary:'+str((cid,R,V))+'/'+name,istrue(row[name])==value if isinstance(value,bool) else math.isclose(float(row[name]),value,rel_tol=1e-12,abs_tol=1e-12))
    def summarize(group):
        return {'n_cases':len(group),'exact_depth':sum(m['absolute_depth_error_m']==0 for m in group),'within_10m':sum(m['absolute_depth_error_m']<=10 for m in group),'rank_le3':sum(m['true_depth_rank']<=3 for m in group),'practical_pass':sum(m['practical_two_bin_pass'] for m in group),'worst_error_m':max(m['absolute_depth_error_m'] for m in group),'worst_true_rank':max(m['true_depth_rank'] for m in group),'minimum_second_margin_db':min(m['DeltaJ_second_db'] for m in group)}
    degradation=read('B1B_RECTANGLE_DEGRADATION_SUMMARY.csv');check('all_degradation_rows',len(degradation)==64 and {(float(r['R_km']),float(r['V_mps']),r['aggregation']) for r in degradation}=={(R,V,a) for R in policy['range_rectangle_limits_km'] for V in policy['speed_rectangle_limits_mps'] for a in policy['aggregations']})
    for row in degradation:
        R,V,a=float(row['R_km']),float(row['V_mps']),row['aggregation']
        for stage in ('ALL','S7','S8'):
            group=[m for (cid,r,v,agg),m in rebuilt.items() if r==R and v==V and agg==a and (stage=='ALL' or cases[cid]['historical_stage']==stage)]
            prefix='' if stage=='ALL' else stage+'_'
            for name,value in summarize(group).items():check('degradation:'+str((R,V,a))+'/'+prefix+name,math.isclose(float(row[prefix+name]),value,abs_tol=1e-12))
    target={a:[m for (cid,R,V,agg),m in rebuilt.items() if R==.25 and V==.05 and agg==a] for a in policy['aggregations']}
    for a,group in target.items():
        for name,value in summarize(group).items():check('target:'+a+'/'+name,math.isclose(decision['target_summaries'][a][name],value,abs_tol=1e-12))
    d1=all(m['true_depth_rank']<=3 for m in target['MIN']);d2=all(m['absolute_depth_error_m']<=10 for m in target['MIN'])
    qhits=sum(m['absolute_depth_error_m']<=10 for m in target['Q25']);s7hits=sum(m['absolute_depth_error_m']<=10 for (cid,R,V,a),m in rebuilt.items() if R==.25 and V==.05 and a=='Q25' and cases[cid]['historical_stage']=='S7');d3=qhits>=5 and s7hits>=2
    label='JOINT_RVZ_PROFILE_ROUTE_WORTH_FRESH_DESIGN' if d1 and d2 and d3 else 'JOINT_RVZ_PROFILE_ROUTE_NOT_SUPPORTED_BY_EXISTING_EVIDENCE'
    check('route_decision',decision['route_decision']==label and decision['D1_MIN_all_six_rank_le3']==d1 and decision['D2_MIN_all_six_within10m']==d2 and decision['D3_Q25_at_least_five_within10m_and_S7_at_least_two']==d3)
    check('no_scientific_credit',decision['B1_scientific_pass'] is False and decision['R4_B1_progress_percent']==decision['R4_overall_progress_percent']==0)
    for name in ('new_propagation_evaluations','new_depth_scores','new_rv_nodes','new_depth_labels','new_cases'):check('no_new_data:'+name,decision[name]==0)
    check('closed_engineering_scope',decision['current_B1_engineering_status']==('NOT_COMPLETED_PLUGIN_NOT_ROBUST_JOINT_ROUTE_DIAGNOSTIC_ONLY' if d1 and d2 and d3 else 'CURRENT_CONDITIONAL_DEPTH_ENGINEERING_ROUTE_CLOSED'))
    for path,digest in policy['bindings'].items():check('policy_binding:'+path,sha(ROOT/path)==digest)
    for path,digest in policy['protected_historical_files'].items():check('historical:'+path,sha(ROOT/path)==digest)
    for sourcefile in policy['sources']:
        tree=ast.parse((ROOT/sourcefile).read_text(encoding='utf-8-sig'))
        imported=[n.module if isinstance(n,ast.ImportFrom) else a.name for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom)) for a in (n.names if isinstance(n,ast.Import) else [None])]
        check('no_forward_imports:'+sourcefile,not any(name and (name.startswith('r4_b1_depth') or name.startswith('r4_b1_independent') or name.startswith('r3_') or name.startswith('scipy')) for name in imported))
    for name in ('R4_PLAN.md','R4_EVIDENCE_LEDGER.csv'):
        path='results/R4_MASTER/'+name;old=subprocess.check_output(['git','show',policy['parent_sha']+':'+path],cwd=ROOT).replace(b'\r\n',b'\n')
        check('append_only:'+name,(ROOT/path).read_bytes().replace(b'\r\n',b'\n').startswith(old))
    check('no_analysis_invalid_marker',not (OUT/'B1B_ANALYSIS_INVALID.json').exists())
    return checks

def finalize():
    if (OUT/'B1B_RESULTS_FREEZE.json').exists():raise RuntimeError('Already sealed')
    audit();d=js('B1B_DIAGNOSTIC_POLICY.json');result=js('B1B_DECISION.json');validation=js('B1B_ANALYSIS_VALIDATION.json')
    xml=Path('D:/ProjectStorage/WSL-CZ/b1b-20261005/tests.xml');suites=list(ET.parse(xml).getroot().iter('testsuite'));tests=sum(int(s.attrib['tests']) for s in suites);assert tests>0 and sum(int(s.attrib['failures'])+int(s.attrib['errors']) for s in suites)==0
    (OUT/'B1B_TESTS.xml').write_bytes(xml.read_bytes());write('B1B_TESTS.json',{'passed':tests,'failed':0})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    original=read_path(ROOT/d['depth_score_path']);curves=read('B1B_PROFILED_DEPTH_CURVES.csv')+read('B1B_MIN_NO_CENTER_DEPTH_CURVES.csv');metrics=read('B1B_PROFILE_METRICS.csv');figures=OUT/'figures';figures.mkdir(exist_ok=True)
    def save(fig,name):
        fig.tight_layout();fig.savefig(figures/(name+'.svg'));fig.savefig(figures/(name+'.png'),dpi=150);plt.close(fig)
        path=figures/(name+'.svg');path.write_text('\n'.join(x.rstrip() for x in path.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    fig,axes=plt.subplots(3,2,figsize=(12,12),sharex=True)
    for ax,case in zip(axes.T.ravel(),d['cases']):
        cid=case['case_id'];h0=sorted([r for r in original if r['case_id']==cid and float(r['delta_r_km'])==float(r['delta_v_mps'])==0],key=lambda r:int(r['z_label_m']))
        ax.plot([int(r['z_label_m']) for r in h0],[float(r['J_z_db']) for r in h0],label='H0',color='black',linewidth=2)
        for a in d['aggregations']:
            rows=sorted([r for r in curves if r['case_id']==cid and float(r['R_km'])==.25 and float(r['V_mps'])==.05 and r['aggregation']==a],key=lambda r:int(r['z_label_m']))
            ax.plot([int(r['z_label_m']) for r in rows],[float(r['J_profiled_db']) for r in rows],marker='.',label=a)
        ax.axvline(case['z_true_m'],color='gray',linestyle=':');ax.set_title(cid);ax.set_ylabel('Saved-lattice aggregate J (dB)');ax.grid(alpha=.2);ax.legend(fontsize=7,ncol=2)
    for ax in axes[-1]:ax.set_xlabel('Source-depth label (m)')
    fig.suptitle('H0 vs target rectangle saved-lattice profiles; route diagnostic only');save(fig,'FIG_1_H0_VS_TARGET_PROFILES')
    maps=read('B1B_RV_COMPENSATION_MAP.csv')
    for field,name,title in [('argmin_dr_km','FIG_2_DEPTH_ARGMIN_RANGE','MIN argmin range offset (km)'),('argmin_dv_mps','FIG_3_DEPTH_ARGMIN_SPEED','MIN argmin speed offset (m/s)')]:
        values=[]
        for case in d['cases']:
            rows=sorted([r for r in maps if r['case_id']==case['case_id'] and float(r['R_km'])==.25 and float(r['V_mps'])==.05],key=lambda r:int(r['z_label_m']))
            values.append([float(r[field]) for r in rows])
        fig,ax=plt.subplots(figsize=(14,4));im=ax.imshow(values,interpolation='nearest',aspect='auto',cmap='coolwarm');ax.set_xticks(range(21),[str(z) for z in d['depth_labels_m']],rotation=45);ax.set_yticks(range(6),[c['case_id'] for c in d['cases']]);ax.set_xlabel('Depth label (m)');ax.set_title(title+'; target rectangle, tested depth labels only');fig.colorbar(im,ax=ax);save(fig,name)
    fig,axes=plt.subplots(4,2,figsize=(15,13));rects=[(R,V) for R in d['range_rectangle_limits_km'] for V in d['speed_rectangle_limits_mps']]
    for i,a in enumerate(d['aggregations']):
        for j,field in enumerate(['absolute_depth_error_m','true_depth_rank']):
            values=[[float(next(m[field] for m in metrics if m['case_id']==c['case_id'] and float(m['R_km'])==R and float(m['V_mps'])==V and m['aggregation']==a)) for R,V in rects] for c in d['cases']]
            ax=axes[i,j];im=ax.imshow(values,interpolation='nearest',aspect='auto');ax.set_yticks(range(6),[c['case_id'] for c in d['cases']],fontsize=7);ax.set_xticks(range(16),[f'{R:g}/{V:g}' for R,V in rects],rotation=90,fontsize=7);ax.set_title(a+' '+('depth error (m)' if j==0 else 'true-depth rank'));fig.colorbar(im,ax=ax)
    fig.suptitle('All 16 frozen rectangles (R km / V m/s); discrete summaries, no interpolation');save(fig,'FIG_4_RECTANGLE_DEGRADATION')
    table='| Case | Aggregation | H0 z-hat | H0 margin (dB) | Profiled z-hat | True rank | Profiled margin (dB) |\n|---|---|---:|---:|---:|---:|---:|\n'
    for row in read('B1B_H0_VS_PROFILED_TARGET.csv'):
        table+=f"| {row['case_id']} | {row['aggregation']} | {row['H0_z_hat_grid_m']} | {float(row['H0_DeltaJ_second_db']):.6g} | {row['profiled_z_hat_grid_m']} | {row['profiled_true_depth_rank']} | {float(row['profiled_DeltaJ_second_db']):.6g} |\n"
    comp=result['target_compensation_summary'];smallest=min(r['MIN_DeltaJ_second_db'] for r in comp);largest=max(r['MIN_DeltaJ_second_db'] for r in comp)
    compensation_text=f"Target wrong-depth best joint margins span {smallest:.6g}-{largest:.6g} dB. MIN true-depth argmin uses center in {sum(r['true_depth_argmin_is_center'] for r in comp)}/6 configurations. At wrong depths, noncenter-node improvement over same-depth H0 occurs in {sum(r['wrong_depth_compensation_label_count'] for r in comp)}/120 case/depth entries. Adjacent-depth argmin jumps per case: {[r['adjacent_depth_node_jumps'] for r in comp]}; range jumps {[r['adjacent_depth_range_jumps'] for r in comp]}, speed jumps {[r['adjacent_depth_speed_jumps'] for r in comp]}. Depth labels numerically equivalent to global minimum per case: {[r['n_depths_equivalent_to_global_min_within_numerical_tolerance'] for r in comp]}."
    report=f"""# B1B saved joint r-v-z diagnostic

**{result['route_decision']}**; ROUTE_DIAGNOSTIC only, PENDING_RESEARCH_LEAD_AUDIT. R4-B1 and overall remain0%. Parent {d['parent_sha']} was independently accepted negative; policy {result['policy_commit']} was pushed before this analysis. New propagation evaluations=0, new depth scores=0, new cases/nodes/depth labels=0. Every aggregate value is exactly a selected existing B1A score with original CSV line/value and support-node provenance.

All original six configurations and16 rectangles are retained. S7 triple-frequency is primary interpretation; S8 four-line is the historical companion, not independent statistical confirmation. MIN is minimum across the saved rectangle; Q25 is the empirical inverse-CDF nearest-rank order statistic ceil(n/4), not interpolated percentile or a continuous optimizer. MEDIAN uses ceil(n/2), with odd support-node counts. Ties select ascending score/dr/dv. MIN_NO_CENTER removes only (0,0) and remains secondary. Q25/MEDIAN horizontal provenance is the selected order-statistic carrier, not an argmin; their argmin columns stay blank. This rule was frozen before execution; no quantile shopping or new rectangle.

Target rectangle R0.25 km/V0.05 m/s: {json.dumps(result['target_summaries'])}. D1={result['D1_MIN_all_six_rank_le3']}, D2={result['D2_MIN_all_six_within10m']}, D3={result['D3_Q25_at_least_five_within10m_and_S7_at_least_two']}. D3 requires Q25 within10m in at least5/6 and at least2/3 S7 cases; MEDIAN/leave-center cannot replace it. D4 all16 rectangles reported without demanding monotonicity. These route-admission diagnostics are not a B1 scientific completion or confirmation Gate.

{table}

{compensation_text}

The maps demonstrate finite-node range/speed compensation where noncenter nodes reduce the same wrong-depth cost, and discrete depth-dependent branch changes. They do not establish a smooth continuous depth-range/speed tradeoff or causal decomposition; no curve is fitted. Numerical-equivalence uses the inherited {d['numerical_tolerance_db']:.12g} dB tolerance, not an invented statistical interval. Positive wrong-depth margins report sampled separation only, never continuous uniqueness. Center is in every primary rectangle and is the exact generation state: preserving its near-zero true score is expected and cannot establish an observation-derived engineering route. Q25, MEDIAN and leave-center diagnose dependence on this exact node; none represents a calibrated uncertainty likelihood or confidence interval.

The previous plug-in target had13/150 practical passes, and all16 rectangles failed universally; the smallest had7/54. B1A preserved369 nonmonotonic witnesses. H0 unique true minima6/6 and positive margins/curvature remain accepted. B1A range/speed conditioning is not robust; prior theta+/-0.5 and psi+/-1 single-factor tests are comparatively stable at those scales, without a joint4D claim. New fine conditioning fractions, B2 depth interpolation, SSP/receiver-depth errors, optimizer or propagation are not run.

{validation['checks_passed']} analysis integrity checks pass; {tests} unit tests pass; separate saved-row audit reconstructs every aggregate selection, metric, compensation map,16-rectangle summary and route decision. Historical files/inputs and immutable source hashes are protected. Three main aggregations have6048 curve rows; secondary no-center has2016; all384 profiles and8064 source mappings are complete. Metadata and file hashing are the only other input reads; analysis imports no forward/depth-score module.

Current engineering status: {result['current_B1_engineering_status']}. A positive route-admission result would only recommend a separately instructed B1R_LOCAL_JOINT_RVZ_SUPPORT_VALIDATION fresh design. A negative result closes the current conditional-depth engineering route and preserves EXACT_HORIZONTAL_DEPTH_MECHANISM as oracle/mechanism evidence. B1R/B2/A2/A3/A4/B3/B4/C/P5 are NOT_OPENED, and CZ certified acquisition remains CLOSED. Stop after analysis commit/push for independent review.

See [curves](B1B_PROFILED_DEPTH_CURVES.csv), [no-center](B1B_MIN_NO_CENTER_DEPTH_CURVES.csv), [source provenance](B1B_SCORE_SOURCE_PROVENANCE.csv), [compensation map](B1B_RV_COMPENSATION_MAP.csv), [all rectangles](B1B_RECTANGLE_DEGRADATION_SUMMARY.csv), [B1 closeout](R4_B1_CLOSEOUT.md).
"""
    (OUT/'B1B_REPORT.md').write_text(report,encoding='utf-8');(OUT/'GPT_SYNC.md').write_text(report,encoding='utf-8')
    closeout=f"""# R4 B1 closeout

B1 engineering NOT_COMPLETED; management progress0/10%, overall R4=0%. Latest route decision {result['route_decision']} and resulting {result['current_B1_engineering_status']} are pending independent research-lead audit.

1. Exact horizontal: matched on-grid depth mechanism strongly identifiable, H0 unique global truth minima6/6, minimum second margin0.794481 dB, positive curvature6/6; W50=0-5 m is descriptive, not a confidence interval. Scope is three source-depth scenarios and two frequency configurations at a single horizontal reference, not six independent scenes.
2. Plug-in r/v: not robust in the tested lattice; no universally passing tested rectangular conditioning envelope. Frozen target0.25 km/0.05 m/s gives13/150 practical passes; smallest0.125/0.025 gives7/54; all16 candidate rectangles fail. No finer step search or lowered target is authorized.
3. Theta/psi: previously tested +/-0.5-degree theta and +/-1-degree psi give12/12 exact each; comparatively stable at frozen single-factor scales, not arbitrary-angle or joint4D robustness.
4. Range/speed: strong depth coupling, worst full-lattice70 m and369 preserved nonmonotonic witnesses. These are assumed upstream conditioned horizontal states, not sensor-error requirements or attainable estimates.
5. B2 off-grid continuous depth remains NOT_OPENED; SSP/B3 and B4 do not automatically advance.

B1B uses only saved B1A data with the policy-first MIN/Q25/MEDIAN and secondary leave-center diagnostics across the same16 rectangles. Target summaries: {json.dumps(result['target_summaries'])}. {compensation_text} Center-included MIN can preserve the exact-generation depth valley while depending on oracle support; it does not solve horizontal acquisition or establish joint estimator robustness. Q25 and leave-center outcomes remain route diagnostics, not confirmation.

{result['route_decision']}: {'a future B1R_LOCAL_JOINT_RVZ_SUPPORT_VALIDATION fresh/frozen design may be considered under a new explicit instruction; current B1 is still incomplete' if result['route_decision']=='JOINT_RVZ_PROFILE_ROUTE_WORTH_FRESH_DESIGN' else 'close CURRENT_CONDITIONAL_DEPTH_ENGINEERING_ROUTE; retain EXACT_HORIZONTAL_DEPTH_MECHANISM as oracle/mechanism evidence only'}. No current B1 PASS/10% claim, no off-grid z/interpolation, no new propagation/depth score/node/case. B1R/B2/A2/A3/A4/B3/B4/C/P5 are NOT_OPENED. CZ cold-start certified global search stays closed; A1 NOT_COMPLETED. No physics impossibility or globally unique continuous joint inverse claim. All historical evidence remains unchanged; only master plan/ledger appended at zero credited progress. Stop after commit/push and await lead audit.
"""
    (OUT/'R4_B1_CLOSEOUT.md').write_text(closeout,encoding='utf-8')
    master=ROOT/'results/R4_MASTER'
    with (master/'R4_PLAN.md').open('a',encoding='utf-8',newline='') as f:
        f.write('\n\n## B1B saved joint r-v-z diagnostic and B1 closeout\n\n'+f"{result['route_decision']}; {result['current_B1_engineering_status']}; PENDING_RESEARCH_LEAD_AUDIT. Policy{result['policy_commit']} pushed before saved-score analysis. Target MIN/Q25/MEDIAN/no-center summaries={json.dumps(result['target_summaries'])}. New propagation/depth scores/nodes/cases=0; all16 rectangles and six configurations retained. No B1 completion credit, B1=0/10%, R4=0%.\n\n"+'[B1 closeout](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/R4_B1_CLOSEOUT.md), [diagnostic report](../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/B1B_REPORT.md). Exact-horizontal mechanism retained; plug-in range/speed not robust, angles stable only at tested scales. B1R/B2/A2/A3/A4/B3/B4/C/P5 NOT_OPENED. A1 NOT_COMPLETED, CZ certified global acquisition CLOSED. Stop after commit/push. R4_PROGRESS.json unchanged.\n')
    with (master/'R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:
        csv.writer(f,lineterminator='\n').writerow(['B1B','../R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT/B1B_DECISION.json',result['route_decision'],'Saved lattice joint-profile route diagnostic and B1 closeout; no new propagation/score','PENDING_RESEARCH_LEAD_AUDIT',0])
    checks=audit()
    with (OUT/'B1B_INDEPENDENT_SAVED_AUDIT.csv').open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['check','pass'],lineterminator='\n');w.writeheader();w.writerows(checks)
    write('B1B_INDEPENDENT_SAVED_AUDIT.json',{'checks_passed':len(checks),'failed':0,'new_propagation':0,'new_depth_scores':0,'independent_source_order_statistic_metric_map_route_audit':True})
    bindings={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='B1B_RESULTS_FREEZE.json'}
    for name in ('R4_PLAN.md','R4_EVIDENCE_LEDGER.csv'):bindings[(master/name).relative_to(ROOT).as_posix()]=sha(master/name)
    write('B1B_RESULTS_FREEZE.json',{'policy_commit':result['policy_commit'],'policy_file_sha256':sha(OUT/'B1B_DIAGNOSTIC_POLICY.json'),'route_decision':result['route_decision'],'immutable_after_commit':True,'bindings':bindings})
    print(json.dumps({'route_decision':result['route_decision'],'saved_audit_checks':len(checks),'tests':tests,'sealed_bindings':len(bindings),'new_propagation':0,'new_depth_scores':0}))

if __name__=='__main__':
    if '--verify' in __import__('sys').argv:
        seal=js('B1B_RESULTS_FREEZE.json')
        for path,digest in seal['bindings'].items():assert sha(ROOT/path)==digest,path
        checks=audit();print(json.dumps({'sealed':'PASS','checks':len(checks),'failed':0,'new_propagation':0,'new_depth_scores':0}))
    else:finalize()
