"""One B1A execution, authorized only after the frozen design is pushed."""
from pathlib import Path
import csv
import hashlib
import inspect
import json
import subprocess
import time
import numpy as np
import r4_b1_depth_profile as primary
import r4_b1_independent as independent
import r4_b1a_boundary as boundary
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def js(name): return json.loads((OUT/name).read_text(encoding='utf-8-sig'))
def write(name,data): (OUT/name).write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as file: return list(csv.DictReader(file))
def csv_write(name,rows,empty_fields=None):
    fields=list(dict.fromkeys(k for row in rows for k in row)) if rows else empty_fields
    if not fields: raise ValueError('Missing CSV schema: '+name)
    with (OUT/name).open('w',encoding='utf-8',newline='') as file:
        writer=csv.DictWriter(file,fieldnames=fields,lineterminator='\n'); writer.writeheader()
        for row in rows: writer.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(list,dict)) else v for k,v in row.items()})

def verify_freeze(require_push=True):
    d=js('B1A_DESIGN_FREEZE.json')
    for path,digest in d['bindings'].items():
        if sha(ROOT/path)!=digest: raise RuntimeError('Frozen source/input changed: '+path)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    if require_push:
        remote=subprocess.check_output(['git','ls-remote','origin','refs/heads/main'],cwd=ROOT,text=True).split()[0]
        if head!=remote: raise RuntimeError('HEAD is not remote/main')
        title=subprocess.check_output(['git','log','-1','--format=%s'],cwd=ROOT,text=True).strip()
        if title!='R4 B1A: freeze range-speed conditioning boundary': raise RuntimeError('Not frozen design HEAD')
        parent=subprocess.check_output(['git','rev-parse','HEAD^'],cwd=ROOT,text=True).strip()
        if parent!=d['parent_sha']: raise RuntimeError('Wrong design parent')
        p=OUT/'B1A_DESIGN_FREEZE.json'
        blob=subprocess.check_output(['git','show',head+':'+p.relative_to(ROOT).as_posix()],cwd=ROOT)
        if blob.replace(b'\r\n',b'\n')!=p.read_bytes().replace(b'\r\n',b'\n'): raise RuntimeError('Uncommitted design')
    return d,head

def load_observations(path,case_frequencies):
    """Read only frozen observed levels; no source depth/reference information returned."""
    rows=read_csv(path); observed={}
    for cid,freqs in case_frequencies.items():
        observed[cid]={}
        for f in freqs:
            group=sorted((r for r in rows if r['case_id']==cid and int(r['frequency_hz'])==f),key=lambda r:int(r['time_s']))
            if [int(r['time_s']) for r in group]!=list(range(0,1201,10)): raise RuntimeError('Observation incomplete')
            observed[cid][f]=np.array([float(r['relative_level_db']) for r in group])
    return observed

def execute():
    d,design_commit=verify_freeze()
    if (OUT/'B1A_EXECUTION_START.json').exists(): raise RuntimeError('Single frozen execution already started')
    write('B1A_EXECUTION_START.json',{'design_commit':design_commit,'design_file_sha256':sha(OUT/'B1A_DESIGN_FREEZE.json'),'remote_main_verified_before_first_evaluation':True,'B1A_numerical_profiles_before_push':0})
    started=time.perf_counter(); cases={c['case_id']:c for c in d['cases']}
    obs=load_observations(ROOT/d['observation_path'],{c['case_id']:c['frequencies_hz'] for c in d['cases']})
    models=primary.load_models(); other=independent.independent_models(); labels=d['depth_labels_m']
    lattice=read_csv(OUT/'B1A_LATTICE.csv'); records=[]; score_rows=[]; reconstruction=[]; checks=[]; deltas=[]; repeats=[]
    def check(name,ok,detail=''):
        checks.append({'check':name,'pass':bool(ok),'detail':detail})
        if not ok: raise RuntimeError('Execution invalid: '+name+' '+str(detail))
    check('all_lattice_rows',len(lattice)==486,len(lattice))
    for n,item in enumerate(lattice):
        cid=item['case_id']; frequencies=cases[cid]['frequencies_hz']; horizontal=[float(item[k]) for k in ('conditioning_r_km','conditioning_theta_deg','conditioning_v_mps','conditioning_psi_deg')]
        scores=primary.depth_profile(horizontal,obs[cid],models,labels,frequencies)
        repeat=primary.depth_profile(horizontal,obs[cid],primary.load_models(),labels,frequencies)
        cold=independent.independent_score(horizontal,obs[cid],other,labels,frequencies)
        rd=float(np.max(np.abs(scores-repeat))); delta=float(np.max(np.abs(scores-cold)))
        repeats.append(rd); deltas.append(delta)
        check(cid+'/'+item['conditioning_id']+'_repeat',rd<=d['numerical_tolerance_rule']['maximum_allowed_absolute_difference_db'],rd)
        check(cid+'/'+item['conditioning_id']+'_independent',delta<=d['numerical_tolerance_rule']['maximum_allowed_absolute_difference_db'],delta)
        base={'case_id':cid,'frequency_configuration':cases[cid]['historical_stage'],'frequencies_hz':frequencies,
              'conditioning_id':item['conditioning_id'],'delta_r_km':float(item['delta_r_km']),'delta_v_mps':float(item['delta_v_mps']),
              'conditioning_r_km':horizontal[0],'conditioning_theta_deg':horizontal[1],'conditioning_v_mps':horizontal[2],'conditioning_psi_deg':horizontal[3]}
        records.append({**base,'scores':scores.tolist()})
        for iz,z in enumerate(labels):
            check(cid+'/'+item['conditioning_id']+'/'+str(z)+'_finite',bool(np.isfinite(scores[iz])),scores[iz])
            mapping={str(f):d['modal_depth_mapping'][str(f)]['source_label_mapping'][str(z)] for f in frequencies}
            score_rows.append({**base,'z_true_m':cases[cid]['z_true_m'],'z_label_m':z,'effective_source_depth_m':mapping,
                               'effective_receiver_depth_m':{str(f):d['modal_depth_mapping'][str(f)]['effective_receiver_depth_m'] for f in frequencies},'J_z_db':float(scores[iz])})
            reconstruction.append({'case_id':cid,'conditioning_id':item['conditioning_id'],'z_label_m':z,'J_primary_db':float(scores[iz]),'J_repeat_db':float(repeat[iz]),'J_independent_db':float(cold[iz]),'absolute_repeat_difference_db':float(abs(scores[iz]-repeat[iz])),'absolute_independent_difference_db':float(abs(scores[iz]-cold[iz]))})
        if (n+1)%27==0: print(json.dumps({'profiles_completed':n+1,'profiles_expected':486,'rows':len(score_rows),'maximum_independent_difference_db':max(deltas)}),flush=True)
    csv_write('B1A_DEPTH_PROFILES.csv',score_rows); csv_write('B1A_NUMERICAL_RECONSTRUCTION.csv',reconstruction)
    floor=d['numerical_tolerance_rule']['floating_floor_multiplier']*np.finfo(float).eps*max(1.,max(max(p['scores']) for p in records))
    tolerance=max(float(floor),d['numerical_tolerance_rule']['residual_multiplier']*max(deltas+repeats+[d['accepted_B1_reconstruction_residual_db']]))
    check('tolerance_execution_ceiling',tolerance<=d['numerical_tolerance_rule']['maximum_allowed_tolerance_db'],tolerance)
    points=[]
    for p in records:
        base={k:v for k,v in p.items() if k!='scores'}; truth=cases[p['case_id']]['z_true_m']
        m={**base,'z_true_m':truth,**primary.profile_metrics(labels,p['scores'],truth,tolerance)}
        m['local_minima_count']=m['strict_interior_local_minima']
        m['practical_two_bin_pass']=boundary.practical(m); m['strict_exact_pass']=boundary.strict(m,tolerance)
        points.append(m)
    csv_write('B1A_POINT_METRICS.csv',points)
    range_axis=[]; speed_axis=[]
    for axis,limits,result in [('range',d['range_rectangle_limits_km'],range_axis),('speed',d['speed_rectangle_limits_mps'],speed_axis)]:
        for level in limits:
            group=[p for p in points if p['delta_v_mps']==0 and abs(p['delta_r_km'])==level] if axis=='range' else [p for p in points if p['delta_r_km']==0 and abs(p['delta_v_mps'])==level]
            check(axis+'_axis_'+str(level)+'_count',len(group)==12,len(group))
            result.append({'axis':axis,'absolute_offset':level,**boundary.summarize(group)})
    csv_write('B1A_RANGE_AXIS_SUMMARY.csv',range_axis); csv_write('B1A_SPEED_AXIS_SUMMARY.csv',speed_axis)
    joint=[]
    for dr in d['range_offsets_km']:
        for dv in d['speed_offsets_mps']:
            group=[p for p in points if p['delta_r_km']==dr and p['delta_v_mps']==dv]
            check('joint_count_'+str((dr,dv)),len(group)==6,len(group))
            row={'delta_r_km':dr,'delta_v_mps':dv,**boundary.summarize(group)}
            for stage in ('S7','S8'):
                subset=[p for p in group if p['frequency_configuration']==stage]
                row.update({stage+'_'+k:v for k,v in boundary.summarize(subset).items()})
            joint.append(row)
    csv_write('B1A_JOINT_LATTICE_SUMMARY.csv',joint)
    rectangles=boundary.rectangles(points,d['range_rectangle_limits_km'],d['speed_rectangle_limits_mps']); csv_write('B1A_ROBUST_RECTANGLES.csv',rectangles)
    target=next(r for r in rectangles if r['is_primary_target']); check('target_coverage_count',target['n_points']==150,target['n_points'])
    configuration=[]
    for stage in ('S7','S8'):
        group=[p for p in points if p['frequency_configuration']==stage]
        target_group=[p for p in group if abs(p['delta_r_km'])<=.25 and abs(p['delta_v_mps'])<=.05]
        configuration.append({'frequency_configuration':stage,**boundary.summarize(group),**{'target_'+k:v for k,v in boundary.summarize(target_group).items()}})
    csv_write('B1A_FREQUENCY_CONFIGURATION_SUMMARY.csv',configuration)
    witnesses=boundary.nonmonotonic_witnesses(points)
    csv_write('B1A_NON_MONOTONIC_WITNESSES.csv',witnesses,['case_id','smaller_dr_km','smaller_dv_mps','smaller_error_m','smaller_true_rank','larger_dr_km','larger_dv_mps','larger_error_m','larger_true_rank','fixed_other_axis'])
    # Replay only the already-run origin and range/speed axis endpoint profiles from accepted B1.
    old=read_csv(ROOT/d['historical_profile_path']); lookup={(r['case_id'],r['conditioning_id'],int(r['z_label_m'])):float(r['J_z_db']) for r in old}
    historical=[]
    for p in records:
        dr,dv=p['delta_r_km'],p['delta_v_mps']; oldid=None
        if dr==dv==0: oldid='H0'
        elif dv==0 and abs(dr)==1: oldid='H1_r_km_'+('MINUS' if dr<0 else 'PLUS')
        elif dr==0 and abs(dv)==.2: oldid='H1_v_mps_'+('MINUS' if dv<0 else 'PLUS')
        if oldid is None: continue
        for z,value in zip(labels,p['scores']):
            accepted=lookup[p['case_id'],oldid,z]; difference=abs(value-accepted)
            check('B1_replay:'+p['case_id']+'/'+oldid+'/'+str(z),difference<=d['numerical_tolerance_rule']['maximum_allowed_absolute_difference_db'],difference)
            historical.append({'case_id':p['case_id'],'conditioning_id':p['conditioning_id'],'B1_conditioning_id':oldid,'z_label_m':z,'B1_J_db':accepted,'B1A_J_db':value,'absolute_difference_db':difference})
    csv_write('B1A_ACCEPTED_B1_PROFILE_REPLAY.csv',historical); check('B1_replay_count',len(historical)==630,len(historical))
    check('complete_profile_count',len(records)==486,len(records)); check('complete_score_rows',len(score_rows)==10206,len(score_rows))
    check('truth_free_profiler',set(inspect.signature(primary.depth_profile).parameters)=={'horizontal','observed','models','depth_labels','frequencies'})
    for path,digest in d['protected_historical_files'].items(): check('protected:'+path,sha(ROOT/path)==digest)
    verify_freeze()
    csv_write('B1A_VALIDATION_CHECKS.csv',checks)
    write('B1A_VALIDATION.json',{'execution_valid':True,'checks_passed':len(checks),'failed':0,'profiles':486,'depth_rows':10206,'independent_score_rows':10206,'cold_repeat_rows':10206,'accepted_B1_replay_rows':630,'accepted_B1_reconstruction_residual_db':d['accepted_B1_reconstruction_residual_db'],'maximum_repeat_difference_db':max(repeats),'maximum_independent_difference_db':max(deltas),'floating_floor_db':float(floor),'numerical_tolerance_db':tolerance,'numerical_tolerance_rule':d['numerical_tolerance_rule'],'historical_files_preserved':len(d['protected_historical_files']),'elapsed_seconds':time.perf_counter()-started,'design_commit':design_commit})
    decision=boundary.scientific_decision(target['all_points_pass'])
    write('B1A_DECISION.json',{'stage':d['stage'],'parent_sha':d['parent_sha'],'B1_audit_status':'ACCEPTED_CONDITIONAL','design_commit':design_commit,'execution_valid':True,'scientific_decision':decision,'primary_target_rectangle':target,'largest_tested_robust_rectangles':[r for r in rectangles if r['pareto_maximal_robust_rectangle']],'non_monotonic_response':bool(witnesses),'non_monotonic_status':'NON_MONOTONIC_DEPTH_CONDITIONING_RESPONSE' if witnesses else 'NOT_OBSERVED_IN_TESTED_LATTICE','non_monotonic_witness_count':len(witnesses),'fixed_other_axis_witness_count':sum(w['fixed_other_axis'] for w in witnesses),'frequency_configuration_comparison':configuration,'full_lattice':boundary.summarize(points),'range_axis_summary':range_axis,'speed_axis_summary':speed_axis,'profiles':486,'depth_rows':10206,'numerical_tolerance_db':tolerance,'proposed_R4_B1_progress_percent':10 if target['all_points_pass'] else 0,'proposed_R4_overall_progress_percent':10 if target['all_points_pass'] else 0,'confirmed_R4_progress_percent':0,'audit_status':'PENDING_RESEARCH_LEAD_AUDIT','claim_scope':'TESTED_LATTICE_BOUNDARY_ONLY; matched E0; six configurations, three source-depth scenarios; fixed theta=0/psi=5; upstream conditioned range/speed, not sensor specifications; no continuous interior/interpolation claim','B2':'NOT_OPENED','A2_A3_A4_B3_B4_C':'NOT_OPENED','P5':'NOT_OPENED','CZ_acquisition_route':'CLOSED_FOR_CERTIFIED_GLOBAL_SEARCH'})
    print(json.dumps({'scientific_decision':decision,'target_points':target['n_points'],'target_failures':target['n_failed'],'nonmonotonic_witnesses':len(witnesses),'validation_checks':len(checks)}),flush=True)

if __name__=='__main__':
    try: execute()
    except Exception as error:
        if (OUT/'B1A_EXECUTION_START.json').exists(): write('B1A_EXECUTION_INVALID.json',{'decision':'EXECUTION_INVALID','error':str(error),'scientific_decision_not_valid':True,'R4_progress_percent':0})
        raise
