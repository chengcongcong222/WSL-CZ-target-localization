"""Policy-gated saved B1A reanalysis. Imports no scientific forward/scoring modules."""
from pathlib import Path
import ast
import csv
import hashlib
import json
import subprocess
import time
import r4_b1b_diagnostic as stats
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_B1B_JOINT_RVZ_PROFILE_DIAGNOSTIC_AND_CLOSEOUT'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def js(name): return json.loads((OUT/name).read_text(encoding='utf-8-sig'))
def write(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as file: return list(csv.DictReader(file))
def csv_write(name,rows):
    if not rows: raise ValueError('Missing output data: '+name)
    with (OUT/name).open('w',encoding='utf-8',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),lineterminator='\n');writer.writeheader()
        for row in rows: writer.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in row.items()})

def verify_policy(push_required=True):
    d=js('B1B_DIAGNOSTIC_POLICY.json')
    for path,digest in d['bindings'].items():
        if sha(ROOT/path)!=digest: raise RuntimeError('Frozen source/input changed: '+path)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    if push_required:
        remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],cwd=ROOT,text=True).split()[0]
        if head!=remote: raise RuntimeError('Policy HEAD not pushed')
        if subprocess.check_output(['git','log','-1','--format=%s'],cwd=ROOT,text=True).strip()!='R4 B1B: freeze joint rvz diagnostic policy': raise RuntimeError('Wrong policy HEAD')
        if subprocess.check_output(['git','rev-parse','HEAD^'],cwd=ROOT,text=True).strip()!=d['parent_sha']: raise RuntimeError('Wrong parent')
        path=OUT/'B1B_DIAGNOSTIC_POLICY.json'; blob=subprocess.check_output(['git','show',head+':'+path.relative_to(ROOT).as_posix()],cwd=ROOT)
        if blob.replace(b'\r\n',b'\n')!=path.read_bytes().replace(b'\r\n',b'\n'): raise RuntimeError('Uncommitted policy')
    return d,head

def analyze():
    d,policy_commit=verify_policy()
    if (OUT/'B1B_ANALYSIS_START.json').exists(): raise RuntimeError('Saved-data analysis already started')
    write('B1B_ANALYSIS_START.json',{'policy_commit':policy_commit,'policy_file_sha256':sha(OUT/'B1B_DIAGNOSTIC_POLICY.json'),'remote_main_verified':True,'new_propagation_evaluations':0,'new_depth_scores':0})
    began=time.perf_counter(); checks=[]
    def check(name,ok,detail=''):
        checks.append({'check':name,'pass':bool(ok),'detail':detail})
        if not ok: raise RuntimeError('Analysis invalid: '+name+' '+str(detail))
    cases={c['case_id']:c for c in d['cases']}; original=read(ROOT/d['depth_score_path']); point_rows=read(ROOT/d['point_metrics_path'])
    check('all_original_score_rows',len(original)==10206,len(original)); check('all_original_point_rows',len(point_rows)==486,len(point_rows))
    records=[]; keys=set()
    for line,row in enumerate(original,start=2):
        r={'case_id':row['case_id'],'delta_r_km':float(row['delta_r_km']),'delta_v_mps':float(row['delta_v_mps']),'z_label_m':int(row['z_label_m']),
           'J_z_db':float(row['J_z_db']),'source_csv_line':line,'source_score_text':row['J_z_db']}
        key=r['case_id'],r['delta_r_km'],r['delta_v_mps'],r['z_label_m']
        check('unique_source:'+str(key),key not in keys); keys.add(key); records.append(r)
    expected={(cid,r,v,z) for cid in cases for r in d['range_offsets_km'] for v in d['speed_offsets_mps'] for z in d['depth_labels_m']}
    check('complete_saved_lattice',keys==expected)
    center={}
    for r in records:
        if (r['delta_r_km'],r['delta_v_mps'])==(0.,0.):center[r['case_id'],r['z_label_m']]=r
    curves=[]; no_center_curves=[]; provenance=[]; maps=[]; no_center_maps=[]; metrics=[]; compensation=[]
    tol=d['numerical_tolerance_db']; labels=d['depth_labels_m']
    for case in d['cases']:
        cid=case['case_id']; truth=case['z_true_m']
        for R in d['range_rectangle_limits_km']:
            for V in d['speed_rectangle_limits_mps']:
                support=[r for r in records if r['case_id']==cid and abs(r['delta_r_km'])<=R and abs(r['delta_v_mps'])<=V]
                by_depth={z:[r for r in support if r['z_label_m']==z] for z in labels}
                selections={}
                for aggregation in d['aggregations']:
                    rows=[]; previous=None
                    for z in labels:
                        selected,n,ordinal=stats.select_record(by_depth[z],aggregation)
                        source_center=center[cid,z]
                        base={'case_id':cid,'frequency_configuration':case['historical_stage'],'R_km':R,'V_mps':V,'aggregation':aggregation,'z_label_m':z,
                              'J_profiled_db':selected['J_z_db'],'n_rv_nodes':n,'argmin_dr_km':selected['delta_r_km'] if aggregation in ('MIN','MIN_NO_CENTER') else '',
                              'argmin_dv_mps':selected['delta_v_mps'] if aggregation in ('MIN','MIN_NO_CENTER') else '',
                              'source_csv_line':selected['source_csv_line'],'source_score_text':selected['source_score_text'],'selected_order_statistic_ordinal':ordinal}
                        rows.append(base)
                        (no_center_curves if aggregation=='MIN_NO_CENTER' else curves).append(base)
                        provenance.append({**base,'source_csv_path':d['depth_score_path'],'source_csv_sha256':d['bindings'][d['depth_score_path']],
                                           'source_delta_r_km':selected['delta_r_km'],'source_delta_v_mps':selected['delta_v_mps'],'profile_value_is_existing_source_value':True})
                        if aggregation in ('MIN','MIN_NO_CENTER'):
                            eligible=[r for r in by_depth[z] if aggregation!='MIN_NO_CENTER' or (r['delta_r_km'],r['delta_v_mps'])!=(0.,0.)]
                            tied=[{'delta_r_km':r['delta_r_km'],'delta_v_mps':r['delta_v_mps']} for r in eligible if r['J_z_db']<=selected['J_z_db']+tol]
                            node=selected['delta_r_km'],selected['delta_v_mps']; noncenter=node!=(0.,0.)
                            item={**base,'z_true_m':truth,'is_wrong_depth':z!=truth,'H0_same_depth_J_db':source_center['J_z_db'],
                                  'improvement_vs_H0_at_same_depth_db':source_center['J_z_db']-selected['J_z_db'],'selected_noncenter_node':noncenter,
                                  'wrong_depth_compensation_supported':z!=truth and noncenter and source_center['J_z_db']-selected['J_z_db']>tol,
                                  'n_numerically_equivalent_min_nodes':len(tied),'numerically_equivalent_min_nodes':tied,
                                  'argmin_node_changed_from_previous_depth':previous is not None and previous!=node,
                                  'argmin_range_changed_from_previous_depth':previous is not None and previous[0]!=node[0],
                                  'argmin_speed_changed_from_previous_depth':previous is not None and previous[1]!=node[1]}
                            (no_center_maps if aggregation=='MIN_NO_CENTER' else maps).append(item);previous=node
                    selections[aggregation]=rows
                    m=stats.curve_metrics(labels,[r['J_profiled_db'] for r in rows],truth,tol)
                    metrics.append({'case_id':cid,'frequency_configuration':case['historical_stage'],'R_km':R,'V_mps':V,'aggregation':aggregation,'metric_scope':'PROFILED_RV_NUISANCE_DEPTH_METRICS','z_true_m':truth,'n_rv_nodes':rows[0]['n_rv_nodes'],**m})
                h0=stats.curve_metrics(labels,[center[cid,z]['J_z_db'] for z in labels],truth,tol)
                minrows=selections['MIN']; bestwrong=min((r for r in minrows if r['z_label_m']!=truth),key=lambda r:(r['J_profiled_db'],r['z_label_m']))
                scope=[r for r in maps if r['case_id']==cid and r['R_km']==R and r['V_mps']==V]
                minmetrics=next(m for m in reversed(metrics) if m['case_id']==cid and m['R_km']==R and m['V_mps']==V and m['aggregation']=='MIN')
                truthrow=next(r for r in minrows if r['z_label_m']==truth)
                compensation.append({'case_id':cid,'frequency_configuration':case['historical_stage'],'R_km':R,'V_mps':V,
                                     'H0_z_hat_grid_m':h0['z_hat_grid_m'],'H0_DeltaJ_second_db':h0['DeltaJ_second_db'],
                                     'MIN_z_hat_grid_m':minmetrics['z_hat_grid_m'],'MIN_true_depth_rank':minmetrics['true_depth_rank'],
                                     'MIN_DeltaJ_second_db':minmetrics['DeltaJ_second_db'],'margin_reduction_vs_H0_db':h0['DeltaJ_second_db']-minmetrics['DeltaJ_second_db'],
                                     'best_wrong_depth_label_m':bestwrong['z_label_m'],'best_wrong_depth_min_J_db':bestwrong['J_profiled_db'],
                                     'best_wrong_depth_dr_km':bestwrong['argmin_dr_km'],'best_wrong_depth_dv_mps':bestwrong['argmin_dv_mps'],
                                     'true_depth_argmin_is_center':truthrow['argmin_dr_km']==truthrow['argmin_dv_mps']==0.,
                                     'wrong_depth_compensation_label_count':sum(r['wrong_depth_compensation_supported'] for r in scope),
                                     'adjacent_depth_node_jumps':sum(r['argmin_node_changed_from_previous_depth'] for r in scope),
                                     'adjacent_depth_range_jumps':sum(r['argmin_range_changed_from_previous_depth'] for r in scope),
                                     'adjacent_depth_speed_jumps':sum(r['argmin_speed_changed_from_previous_depth'] for r in scope),
                                     'n_depths_equivalent_to_global_min_within_numerical_tolerance':minmetrics['min_tie_count']})
    csv_write('B1B_PROFILED_DEPTH_CURVES.csv',curves);csv_write('B1B_MIN_NO_CENTER_DEPTH_CURVES.csv',no_center_curves)
    csv_write('B1B_SCORE_SOURCE_PROVENANCE.csv',provenance);csv_write('B1B_RV_COMPENSATION_MAP.csv',maps);csv_write('B1B_MIN_NO_CENTER_RV_COMPENSATION_MAP.csv',no_center_maps)
    csv_write('B1B_PROFILE_METRICS.csv',metrics);csv_write('B1B_COMPENSATION_SUMMARY.csv',compensation)
    target=[m for m in metrics if m['R_km']==.25 and m['V_mps']==.05]
    comparison=[]
    for case in d['cases']:
        cid=case['case_id'];h0=stats.curve_metrics(labels,[center[cid,z]['J_z_db'] for z in labels],case['z_true_m'],tol)
        for m in target:
            if m['case_id']==cid:comparison.append({'case_id':cid,'frequency_configuration':case['historical_stage'],'aggregation':m['aggregation'],
                                                    'H0_z_hat_grid_m':h0['z_hat_grid_m'],'H0_DeltaJ_second_db':h0['DeltaJ_second_db'],
                                                    'profiled_z_hat_grid_m':m['z_hat_grid_m'],'profiled_true_depth_rank':m['true_depth_rank'],
                                                    'profiled_DeltaJ_second_db':m['DeltaJ_second_db'],'profiled_absolute_depth_error_m':m['absolute_depth_error_m']})
    csv_write('B1B_H0_VS_PROFILED_TARGET.csv',comparison)
    degradation=[]
    for R in d['range_rectangle_limits_km']:
        for V in d['speed_rectangle_limits_mps']:
            for a in d['aggregations']:
                group=[m for m in metrics if m['R_km']==R and m['V_mps']==V and m['aggregation']==a]
                degradation.append({'R_km':R,'V_mps':V,'aggregation':a,**stats.summarize(group),**{'S7_'+k:v for k,v in stats.summarize([m for m in group if m['frequency_configuration']=='S7']).items()},**{'S8_'+k:v for k,v in stats.summarize([m for m in group if m['frequency_configuration']=='S8']).items()}})
    csv_write('B1B_RECTANGLE_DEGRADATION_SUMMARY.csv',degradation)
    check('main_curve_rows',len(curves)==6048,len(curves));check('no_center_curve_rows',len(no_center_curves)==2016,len(no_center_curves));check('profile_metrics',len(metrics)==384,len(metrics));check('complete_provenance',len(provenance)==8064,len(provenance))
    pointlookup={(r['case_id'],float(r['delta_r_km']),float(r['delta_v_mps'])):r for r in point_rows}
    for case in d['cases']:
        cid=case['case_id'];h0=stats.curve_metrics(labels,[center[cid,z]['J_z_db'] for z in labels],case['z_true_m'],tol);previous=pointlookup[cid,0.,0.]
        for field in ('z_hat_grid_m','true_depth_rank','DeltaJ_second_db','DeltaJ_neighbor_db','K_z_db'):
            check('accepted_H0:'+cid+'/'+field,abs(h0[field]-float(previous[field]))<=tol)
    for row in provenance:
        originalrow=original[row['source_csv_line']-2]
        check('source:'+str((row['case_id'],row['R_km'],row['V_mps'],row['aggregation'],row['z_label_m'])),row['source_score_text']==originalrow['J_z_db'] and row['J_profiled_db']==float(originalrow['J_z_db']))
    for path,digest in d['protected_historical_files'].items():check('protected:'+path,sha(ROOT/path)==digest)
    verify_policy();csv_write('B1B_ANALYSIS_INTEGRITY.csv',checks)
    route=stats.route_decision(target)
    summary={a:stats.summarize([m for m in target if m['aggregation']==a]) for a in d['aggregations']}
    targetcomp=[r for r in compensation if r['R_km']==.25 and r['V_mps']==.05]
    status='NOT_COMPLETED_PLUGIN_NOT_ROBUST_JOINT_ROUTE_DIAGNOSTIC_ONLY' if route['route_decision']=='JOINT_RVZ_PROFILE_ROUTE_WORTH_FRESH_DESIGN' else 'CURRENT_CONDITIONAL_DEPTH_ENGINEERING_ROUTE_CLOSED'
    write('B1B_DECISION.json',{'stage':d['stage'],'parent_sha':d['parent_sha'],'B1A_independent_audit':'ACCEPTED_NEGATIVE','policy_commit':policy_commit,'analysis_valid':True,
                              'diagnostic_scope':'DISCRETE_SAVED_LATTICE_PROFILE_DIAGNOSTIC; ROUTE_DIAGNOSTIC_ONLY; NOT_CONFIRMATION',
                              'primary_target_R_km':.25,'primary_target_V_mps':.05,'target_summaries':summary,**route,
                              'D4_full_rectangle_degradation_reported':len(degradation)==64,'target_compensation_summary':targetcomp,
                              'current_B1_engineering_status':status,'exact_horizontal_depth_mechanism':'SUPPORTED_ORACLE_MECHANISM_ONLY',
                              'B1_scientific_pass':False,'R4_B1_progress_percent':0,'R4_overall_progress_percent':0,'audit_status':'PENDING_RESEARCH_LEAD_AUDIT',
                              'new_propagation_evaluations':0,'new_depth_scores':0,'new_rv_nodes':0,'new_depth_labels':0,'new_cases':0,
                              'B1R':'NOT_OPENED_FRESH_DESIGN_REQUIRES_NEW_STAGE','B2':'NOT_OPENED','A2_A3_A4_B3_B4_C':'NOT_OPENED','P5':'NOT_OPENED','CZ_acquisition_route':'CLOSED_FOR_CERTIFIED_GLOBAL_SEARCH'})
    write('B1B_ANALYSIS_VALIDATION.json',{'analysis_valid':True,'checks_passed':len(checks),'failed':0,'main_curve_rows':6048,'no_center_curve_rows':2016,'metrics':384,'source_provenance_rows':8064,'protected_files':len(d['protected_historical_files']),'numerical_tolerance_db':tol,'new_propagation_evaluations':0,'new_depth_scores':0,'elapsed_seconds':time.perf_counter()-began,'policy_commit':policy_commit})
    print(json.dumps({'route_decision':route['route_decision'],'target_summaries':summary,'new_propagation':0,'new_depth_scores':0,'integrity_checks':len(checks)}))

if __name__=='__main__':
    try:analyze()
    except Exception as error:
        if (OUT/'B1B_ANALYSIS_START.json').exists():write('B1B_ANALYSIS_INVALID.json',{'decision':'ANALYSIS_INVALID','error':str(error),'new_propagation':0,'new_depth_scores':0,'R4_progress_percent':0})
        raise
