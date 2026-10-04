"""Execute only the pushed immutable B1 design, then save all profiles and audits."""
from pathlib import Path
import csv
import hashlib
import inspect
import json
import subprocess
import sys
import time
import numpy as np
import r4_b1_depth_profile as primary
import r4_b1_independent as independent
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_B1_CONDITIONAL_DEPTH_IDENTIFIABILITY'
FREEZE=OUT/'B1_DESIGN_FREEZE.json'

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read_json(path): return json.loads(path.read_text(encoding='utf-8-sig'))
def write_json(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csv_write(name,rows):
    assert rows
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with (OUT/name).open('w',newline='',encoding='utf-8') as file:
        writer=csv.DictWriter(file,fieldnames=fields,lineterminator='\n'); writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(list,dict)) else v for k,v in row.items()})

def verify_design(require_remote=True):
    design=read_json(FREEZE)
    for path,digest in design['bindings'].items():
        if sha(ROOT/path)!=digest: raise RuntimeError('Frozen design/input changed: '+path)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    if require_remote:
        remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],cwd=ROOT,text=True).split()[0]
        if head!=remote: raise RuntimeError('Design HEAD not pushed')
        message=subprocess.check_output(['git','log','-1','--format=%s'],cwd=ROOT,text=True).strip()
        if message!='R4 B1: freeze conditional depth identifiability design': raise RuntimeError('Not design-freeze commit')
        parent=subprocess.check_output(['git','rev-parse','HEAD^'],cwd=ROOT,text=True).strip()
        if parent!=design['parent_sha']: raise RuntimeError('Unexpected design parent')
        blob=subprocess.check_output(['git','show',head+':'+FREEZE.relative_to(ROOT).as_posix()],cwd=ROOT)
        if hashlib.sha256(blob.replace(b'\r\n',b'\n')).hexdigest()!=hashlib.sha256(FREEZE.read_bytes().replace(b'\r\n',b'\n')).hexdigest(): raise RuntimeError('Uncommitted design')
    return design,head

def metrics_rows(profiles,cases,tolerance):
    rows=[]
    for profile in profiles:
        reference=cases[profile['case_id']]
        result=primary.profile_metrics(list(range(150,251,5)),profile['scores'],reference['z_true_m'],tolerance)
        rows.append({**{k:v for k,v in profile.items() if k!='scores'},'z_true_m':reference['z_true_m'],**result})
    return rows

def execute():
    design,commit=verify_design()
    if (OUT/'B1_EXECUTION_START.json').exists(): raise RuntimeError('Single frozen execution already started; no rerun')
    write_json('B1_EXECUTION_START.json',{'design_sha':commit,'design_file_sha256':sha(FREEZE),'remote_main_verified':True,'new_B1_depth_evaluations_before_push':0})
    started=time.perf_counter(); models=primary.load_models(); other=independent.independent_models()
    cases={row['case_id']:row for row in design['cases']}; labels=design['depth_labels_m']
    observations={}; obs_rows=[]; checks=[]; errors=[]; repeated_errors=[]
    def check(name,condition,detail=''):
        checks.append({'check':name,'pass':bool(condition),'detail':detail})
        if not condition: raise RuntimeError('Execution invalid: '+name+' '+str(detail))
    for case in design['cases']:
        cid=case['case_id']; reference=case['horizontal_reference']; depth=case['z_true_m']
        raw=primary.raw_levels(reference,models,[depth]); obs=primary.centered(raw)
        independent_obs=independent.independent_center(independent.independent_levels(reference,other,[depth]))
        delta=max(float(np.max(np.abs(obs[f]-independent_obs[f]))) for f in case['frequencies_hz'])
        errors.append(delta); check(cid+'_independent_observation',delta<=design['reconstruction']['maximum_allowed_absolute_difference_db'],delta)
        observations[cid]={f:obs[f][0].copy() for f in case['frequencies_hz']}
        for f in case['frequencies_hz']:
            iz=int(np.argmin(np.abs(models[f]['depths']-depth)))
            for ti,t in enumerate(range(0,1201,10)):
                obs_rows.append({'case_id':cid,'frequency_hz':f,'time_s':t,'window_id':'W1' if ti<61 else 'W2',
                                 'raw_level_db':float(raw[f][0,ti]),'relative_level_db':float(obs[f][0,ti]),
                                 'effective_source_depth_m':float(models[f]['depths'][iz]),'provenance':'RECONSTRUCTED_ACCEPTED_R3_GENERATION_RULE_AFTER_DESIGN_PUSH'})
    csv_write('B1_OBSERVATIONS.csv',obs_rows)
    profiles=[]; profile_rows=[]; reconstruction_rows=[]
    for case in design['cases']:
        cid=case['case_id']; frequencies=case['frequencies_hz']; obs=observations[cid]
        for condition in design['horizontal_conditioning'][cid]:
            state=condition['horizontal_state']
            scores=primary.depth_profile(state,obs,models,labels,frequencies)
            repeat=primary.depth_profile(state,obs,primary.load_models(),labels,frequencies)
            cold=independent.independent_score(state,obs,other,labels,frequencies)
            rd=float(np.max(np.abs(scores-repeat))); cd=float(np.max(np.abs(scores-cold)))
            repeated_errors.append(rd); errors.append(cd)
            check(cid+'_'+condition['conditioning_id']+'_repeat',rd<=design['reconstruction']['maximum_allowed_absolute_difference_db'],rd)
            check(cid+'_'+condition['conditioning_id']+'_independent',cd<=design['reconstruction']['maximum_allowed_absolute_difference_db'],cd)
            base={'case_id':cid,'historical_stage':case['historical_stage'],'frequencies_hz':frequencies,**condition}
            profiles.append({**base,'scores':scores.tolist()})
            for iz,z in enumerate(labels):
                mapping={str(f):float(models[f]['depths'][np.argmin(np.abs(models[f]['depths']-z))]) for f in frequencies}
                check(cid+'_'+condition['conditioning_id']+'_z'+str(z)+'_finite',bool(np.isfinite(scores[iz])),scores[iz])
                profile_rows.append({**base,'z_label_m':z,'effective_source_depth_m':mapping,'effective_receiver_depth_m':{str(f):float(models[f]['depths'][np.argmin(np.abs(models[f]['depths']-200.))]) for f in frequencies},'J_z_db':float(scores[iz])})
                reconstruction_rows.append({'case_id':cid,'conditioning_id':condition['conditioning_id'],'z_label_m':z,
                                            'J_primary_db':float(scores[iz]),'J_repeat_db':float(repeat[iz]),'J_independent_db':float(cold[iz]),
                                            'absolute_repeat_difference_db':float(abs(scores[iz]-repeat[iz])),
                                            'absolute_independent_difference_db':float(abs(scores[iz]-cold[iz]))})
            print(json.dumps({'profile_complete':cid+'/'+condition['conditioning_id'],'rows':len(profile_rows),'repeat_difference_db':rd,'independent_difference_db':cd}),flush=True)
    csv_write('B1_DEPTH_PROFILES.csv',profile_rows); csv_write('B1_NUMERICAL_RECONSTRUCTION.csv',reconstruction_rows)
    maximum_j=max(max(p['scores']) for p in profiles)
    floor=design['reconstruction']['floating_floor_multiplier']*np.finfo(float).eps*max(1.,maximum_j)
    tolerance=max(float(floor),design['reconstruction']['residual_multiplier']*max(errors+repeated_errors))
    check('tolerance_below_execution_ceiling',tolerance<=design['reconstruction']['maximum_allowed_tolerance_db'],tolerance)
    metrics=metrics_rows(profiles,cases,tolerance); csv_write('B1_PROFILE_METRICS.csv',metrics)
    h0=[m for m in metrics if m['conditioning_type']=='H0']; h1=[m for m in metrics if m['conditioning_type']=='H1']
    h0_summary=[]; h1_summary=[]
    for m in h0:
        unique=m['z_hat_grid_m']==m['z_true_m'] and m['min_tie_count']==1
        separation=m['DeltaJ_second_db']>tolerance and m['DeltaJ_neighbor_db']>tolerance
        curvature=m['K_z_db'] is None or m['K_z_db']>tolerance
        h0_summary.append({**m,'G_B1_1_unique_truth_minimum':unique,'G_B1_2_positive_separation':separation,
                           'G_B1_3_positive_curvature':curvature,'G_B1_4_no_severe_global_tie':unique})
    for cid in cases:
        group=[m for m in h1 if m['case_id']==cid]
        h1_summary.append({'case_id':cid,'n_perturbations':len(group),'exact_depth_hits':sum(m['absolute_depth_error_m']==0 for m in group),
                           'within_5m':sum(m['absolute_depth_error_m']<=5 for m in group),
                           'within_10m':sum(m['absolute_depth_error_m']<=10 for m in group),
                           'true_rank_le3':sum(m['true_depth_rank']<=3 for m in group),
                           'joint_two_bin_rank_gate_count':sum(m['absolute_depth_error_m']<=10 and m['true_depth_rank']<=3 for m in group),
                           'worst_error_m':max(m['absolute_depth_error_m'] for m in group),
                           'H1_gate':all(m['absolute_depth_error_m']<=10 and m['true_depth_rank']<=3 for m in group)})
    csv_write('B1_H0_SUMMARY.csv',h0_summary); csv_write('B1_H1_SUMMARY.csv',h1_summary)
    check('complete_profiles',len(profiles)==54,len(profiles)); check('complete_depth_rows',len(profile_rows)==1134,len(profile_rows))
    check('truth_free_profile_signature',set(inspect.signature(primary.depth_profile).parameters)=={'horizontal','observed','models','depth_labels','frequencies'})
    # Rebuild the accepted final surviving psi=4 and psi=5 minima without selecting any new cases.
    with (ROOT/design['historical_final_survivors']).open(encoding='utf-8-sig',newline='') as file:
        old=list(csv.DictReader(file))
    historical_comparisons=[]
    for row in old:
        cid=next(c['case_id'] for c in design['cases'] if c['historical_stage']==row['stage'] and c['z_true_m']==float(row['z_true']))
        state=[float(row[a]) for a in ('r','theta','v','psi')]
        p=primary.depth_profile(state,observations[cid],models,labels,cases[cid]['frequencies_hz'])
        difference=abs(float(p.min())-float(row['J']))
        check(cid+'_historical_node_'+row['node_id'],difference<=design['reconstruction']['maximum_allowed_absolute_difference_db'] and labels[int(p.argmin())]==float(row['z_star']),difference)
        historical_comparisons.append({'case_id':cid,'node_id':row['node_id'],'historical_J_db':float(row['J']),'reconstructed_J_db':float(p.min()),'absolute_difference_db':difference,'historical_z_star_m':float(row['z_star']),'reconstructed_z_star_m':labels[int(p.argmin())]})
    csv_write('B1_HISTORICAL_SCORE_REPLAY.csv',historical_comparisons)
    for path,digest in design['protected_historical_files'].items(): check('protected:'+path,sha(ROOT/path)==digest)
    verify_design(require_remote=True)
    h0_pass=all(m['G_B1_1_unique_truth_minimum'] and m['G_B1_2_positive_separation'] and m['G_B1_3_positive_curvature'] and m['G_B1_4_no_severe_global_tie'] for m in h0_summary)
    h1_pass=all(m['H1_gate'] for m in h1_summary)
    decision='B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_ESTABLISHED' if h0_pass and h1_pass else 'B1_DEPTH_IDENTIFIABLE_ONLY_UNDER_TIGHT_HORIZONTAL_CONDITIONING' if h0_pass else 'B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_NOT_ESTABLISHED'
    csv_write('B1_VALIDATION_CHECKS.csv',checks)
    validation={'execution_valid':True,'checks_passed':len(checks),'checks_failed':0,'profiles':len(profiles),'depth_rows':len(profile_rows),
                'independent_score_rows':len(reconstruction_rows),'repeat_score_rows':len(reconstruction_rows),'historical_score_replays':len(historical_comparisons),
                'maximum_repeat_difference_db':max(repeated_errors),'maximum_independent_difference_including_observations_db':max(errors),
                'floating_floor_db':float(floor),'numerical_reconstruction_tolerance_db':tolerance,'tolerance_rule':design['reconstruction'],
                'historical_files_preserved':len(design['protected_historical_files']),'elapsed_seconds':time.perf_counter()-started,'design_commit':commit}
    write_json('B1_VALIDATION.json',validation)
    result={'stage':design['stage'],'parent_sha':design['parent_sha'],'design_freeze_sha':commit,'scientific_decision':decision,
            'CZ_acquisition_route':'CLOSED_FOR_CERTIFIED_GLOBAL_SEARCH','R4_A1':'BLOCKED / NOT_COMPLETED','R4_A1_progress_percent':0,
            'H0_case_count':len(h0),'H0_unique_true_minima':sum(m['G_B1_1_unique_truth_minimum'] for m in h0_summary),
            'H0_positive_second_best_margin':sum(m['DeltaJ_second_db']>tolerance for m in h0),
            'H0_positive_neighbor_margin':sum(m['DeltaJ_neighbor_db']>tolerance for m in h0),
            'H0_positive_curvature':sum(m['K_z_db'] is not None and m['K_z_db']>tolerance for m in h0),
            'H1_total':len(h1),'H1_exact_depth':sum(m['absolute_depth_error_m']==0 for m in h1),
            'H1_within_5m':sum(m['absolute_depth_error_m']<=5 for m in h1),'H1_within_10m':sum(m['absolute_depth_error_m']<=10 for m in h1),
            'H1_true_rank_le3':sum(m['true_depth_rank']<=3 for m in h1),'worst_H1_depth_error_m':max(m['absolute_depth_error_m'] for m in h1),
            'minimum_H0_DeltaJ_second_db':min(m['DeltaJ_second_db'] for m in h0),'minimum_H0_DeltaJ_neighbor_db':min(m['DeltaJ_neighbor_db'] for m in h0),
            'minimum_H1_DeltaJ_second_db':min(m['DeltaJ_second_db'] for m in h1),'minimum_H1_DeltaJ_neighbor_db':min(m['DeltaJ_neighbor_db'] for m in h1),
            'H0_interior_minima_counts':{m['case_id']:m['strict_interior_local_minima'] for m in h0},
            'H0_W50_m':{m['case_id']:m['W50_z_m'] for m in h0},'H1_W50_range_m':[min(m['W50_z_m'] for m in h1 if m['W50_z_m'] is not None),max(m['W50_z_m'] for m in h1 if m['W50_z_m'] is not None)],
            'numerical_tolerance_db':tolerance,'execution_valid':True,'primary_case_ids':list(cases),
            'historical_definition_correction':'Six configurations are S7 triple x3 depths plus S8 four-line x3 depths; three distinct synthetic source-depth scenarios, not six independent replicates. No new score selected.',
            'claim_scope':'ORACLE/EXACT-HORIZONTAL-CONDITIONED DEPTH MECHANISM; synthetic matched E0; on-grid labels; S7 triple and historical S8 four-line companion; no end-to-end performance claim',
            'proposed_R4_B1_progress_percent':10 if decision=='B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_ESTABLISHED' else 0,
            'proposed_R4_overall_progress_percent':10 if decision=='B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_ESTABLISHED' else 0,
            'confirmed_R4_progress_percent':0,'audit_status':'PENDING_RESEARCH_LEAD_AUDIT','B2':'NOT_OPENED','A2_A3_A4_B3_B4_C':'NOT_OPENED','P5':'NOT_OPENED'}
    write_json('B1_DECISION.json',result)
    print(json.dumps({'decision':decision,'H0_unique':result['H0_unique_true_minima'],'H1_within10':result['H1_within_10m'],'H1_total':len(h1),'tolerance_db':tolerance,'validation_checks':len(checks)}),flush=True)

if __name__=='__main__':
    try: execute()
    except Exception as error:
        if (OUT/'B1_EXECUTION_START.json').exists():
            write_json('B1_EXECUTION_INVALID.json',{'decision':'EXECUTION_INVALID','error':str(error),'scientific_decision_not_valid':True,'R4_progress_percent':0})
        raise
