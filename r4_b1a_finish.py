"""Independent reconstruction of B1A saved scores, metrics, summaries and Gate."""
from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import json
import math
import subprocess
import xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY'

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def js(name): return json.loads((OUT/name).read_text(encoding='utf-8-sig'))
def write(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def read(name):
    with (OUT/name).open(encoding='utf-8-sig',newline='') as file: return list(csv.DictReader(file))
def bool_field(value): return value if isinstance(value,bool) else value=='True'

def verify_saved():
    d=js('B1A_DESIGN_FREEZE.json'); decision=js('B1A_DECISION.json'); validation=js('B1A_VALIDATION.json')
    checks=[]
    def check(name,ok):
        checks.append({'check':name,'pass':bool(ok)})
        if not ok: raise RuntimeError('Saved reconstruction failed: '+name)
    lattice=read('B1A_LATTICE.csv'); expected={(r['case_id'],r['conditioning_id']):r for r in lattice}
    groups=defaultdict(list)
    for row in read('B1A_DEPTH_PROFILES.csv'): groups[row['case_id'],row['conditioning_id']].append(row)
    metric_rows=read('B1A_POINT_METRICS.csv')
    check('metric_row_count',len(metric_rows)==486)
    metrics={(row['case_id'],row['conditioning_id']):row for row in metric_rows}
    check('complete_lattice_keys',len(expected)==486 and set(groups)==set(metrics)==set(expected))
    cases={c['case_id']:c for c in d['cases']}; tol=validation['numerical_tolerance_db']; points=[]
    for key,rows in groups.items():
        rows.sort(key=lambda r:int(r['z_label_m'])); z=[int(r['z_label_m']) for r in rows]; j=[float(r['J_z_db']) for r in rows]
        check(str(key)+'_21_labels',z==list(range(150,251,5)))
        check(str(key)+'_finite',all(math.isfinite(x) for x in j))
        case=cases[key[0]]; true=case['z_true_m']; index=z.index(true); best=min(range(21),key=lambda i:j[i]); low=j[best]
        rank=1+sum(value<j[index] for value in j); second=min(j[i] for i in range(21) if i!=index)-j[index]
        neighbor=min(j[i] for i in (index-1,index+1))-j[index]; curvature=j[index-1]-2*j[index]+j[index+1]
        local=sum(j[i]<j[i-1] and j[i]<j[i+1] for i in range(1,20)); error=abs(z[best]-true)
        ties=sum(value<=low+tol for value in j); width=None; component=[]
        if max(j)-low>tol:
            mask=[(value-low)/(max(j)-low)<=.5 for value in j]; a=b=best
            while a>0 and mask[a-1]: a-=1
            while b<20 and mask[b+1]: b+=1
            component=z[a:b+1]; width=z[b]-z[a]
        p={'case_id':key[0],'frequency_configuration':case['historical_stage'],'delta_r_km':float(expected[key]['delta_r_km']),'delta_v_mps':float(expected[key]['delta_v_mps']),
           'z_hat_grid_m':z[best],'z_true_m':true,'true_depth_rank':rank,'absolute_depth_error_m':error,'DeltaJ_second_db':second,'DeltaJ_neighbor_db':neighbor,
           'K_z_db':curvature,'local_minima_count':local,'strict_interior_local_minima':local,'W50_z_m':width,
           'J_true_db':j[index],'J_min_db':low,'min_tie_count':ties,'dynamic_range_db':max(j)-low,
           'boundary_local_minima':int(j[0]<j[1])+int(j[-1]<j[-2]),
           'tolerance_aware_true_depth_rank':1+sum(value<j[index]-tol for value in j)}
        p['practical_two_bin_pass']=error<=10 and rank<=3
        p['strict_exact_pass']=z[best]==true and rank==1 and second>tol
        stored=metrics[key]
        for name,value in p.items():
            if isinstance(value,str): valid=stored[name]==value
            elif isinstance(value,bool): valid=bool_field(stored[name])==value
            elif value is None: valid=stored[name]==''
            else: valid=math.isclose(float(stored[name]),value,rel_tol=1e-12,abs_tol=1e-12)
            check(str(key)+'_metric_'+name,valid)
        check(str(key)+'_W50_component',json.loads(stored['W50_component_labels'])==component)
        for r in rows:
            check(str(key)+'_horizontal_'+r['z_label_m'],all(float(r[k])==float(expected[key][k]) for k in ['delta_r_km','delta_v_mps','conditioning_r_km','conditioning_theta_deg','conditioning_v_mps','conditioning_psi_deg']) and float(r['z_true_m'])==true and r['frequency_configuration']==case['historical_stage'])
            mapping=json.loads(r['effective_source_depth_m']); receiver=json.loads(r['effective_receiver_depth_m'])
            check(str(key)+'_depth_mapping_'+r['z_label_m'],mapping=={str(f):d['modal_depth_mapping'][str(f)]['source_label_mapping'][r['z_label_m']] for f in case['frequencies_hz']} and receiver=={str(f):d['modal_depth_mapping'][str(f)]['effective_receiver_depth_m'] for f in case['frequencies_hz']})
        points.append(p)
    def aggregate(group):
        return {'n_points':len(group),'n_exact':sum(p['absolute_depth_error_m']==0 for p in group),'n_within_5m':sum(p['absolute_depth_error_m']<=5 for p in group),'n_within_10m':sum(p['absolute_depth_error_m']<=10 for p in group),'n_rank_le3':sum(p['true_depth_rank']<=3 for p in group),'n_practical_pass':sum(p['practical_two_bin_pass'] for p in group),'n_strict_pass':sum(p['strict_exact_pass'] for p in group),'n_failed':sum(not p['practical_two_bin_pass'] for p in group),'all_points_pass':all(p['practical_two_bin_pass'] for p in group),'worst_depth_error_m':max(p['absolute_depth_error_m'] for p in group),'worst_true_depth_rank':max(p['true_depth_rank'] for p in group),'min_DeltaJ_second_db':min(p['DeltaJ_second_db'] for p in group),'min_DeltaJ_neighbor_db':min(p['DeltaJ_neighbor_db'] for p in group)}
    def compare(label,stored,rebuilt,prefix=''):
        for name,value in rebuilt.items():
            got=stored[prefix+name]
            check(label+prefix+name,bool_field(got)==value if isinstance(value,bool) else math.isclose(float(got),value,rel_tol=1e-12,abs_tol=1e-12))
    for name,axis in [('B1A_RANGE_AXIS_SUMMARY.csv','range'),('B1A_SPEED_AXIS_SUMMARY.csv','speed')]:
        rows=read(name); check(name+'_count',len(rows)==4)
        required=d['range_rectangle_limits_km'] if axis=='range' else d['speed_rectangle_limits_mps']
        check(name+'_keys',{float(r['absolute_offset']) for r in rows}==set(required))
        for row in rows:
            level=float(row['absolute_offset'])
            group=[p for p in points if p['delta_v_mps']==0 and abs(p['delta_r_km'])==level] if axis=='range' else [p for p in points if p['delta_r_km']==0 and abs(p['delta_v_mps'])==level]
            compare(name+str(level),row,aggregate(group)); check(name+str(level)+'_coverage',len(group)==12)
    joint=read('B1A_JOINT_LATTICE_SUMMARY.csv'); check('joint_summary_count',len(joint)==81)
    check('all_joint_keys',{(float(r['delta_r_km']),float(r['delta_v_mps'])) for r in joint}=={(r,v) for r in d['range_offsets_km'] for v in d['speed_offsets_mps']})
    for row in joint:
        dr,dv=float(row['delta_r_km']),float(row['delta_v_mps']); group=[p for p in points if p['delta_r_km']==dr and p['delta_v_mps']==dv]
        compare('joint'+str((dr,dv)),row,aggregate(group)); check('joint_coverage'+str((dr,dv)),len(group)==6)
        for stage in ('S7','S8'): compare('joint'+str((dr,dv)),row,aggregate([p for p in group if p['frequency_configuration']==stage]),stage+'_')
    rebuilt_rect=[]; rectangles=read('B1A_ROBUST_RECTANGLES.csv')
    check('all_16_rectangle_keys',len(rectangles)==16 and {(float(r['R_km']),float(r['V_mps'])) for r in rectangles}=={(r,v) for r in [.125,.25,.5,1] for v in [.025,.05,.1,.2]})
    for row in rectangles:
        r,v=float(row['R_km']),float(row['V_mps']); group=[p for p in points if abs(p['delta_r_km'])<=r and abs(p['delta_v_mps'])<=v]
        summary=aggregate(group); compare('rectangle'+str((r,v)),row,summary)
        rebuilt_rect.append({'R_km':r,'V_mps':v,**summary})
    frontier=[r for r in rebuilt_rect if r['all_points_pass'] and not any(q['all_points_pass'] and q['R_km']>=r['R_km'] and q['V_mps']>=r['V_mps'] and (q['R_km']>r['R_km'] or q['V_mps']>r['V_mps']) for q in rebuilt_rect)]
    for stored,row in zip(rectangles,rebuilt_rect):
        check('frontier:'+str((row['R_km'],row['V_mps'])),bool_field(stored['pareto_maximal_robust_rectangle'])==(row in frontier))
        check('target_marker:'+str((row['R_km'],row['V_mps'])),bool_field(stored['is_primary_target'])==(row['R_km']==.25 and row['V_mps']==.05))
    check('decision_frontier',{(r['R_km'],r['V_mps']) for r in decision['largest_tested_robust_rectangles']}=={(r['R_km'],r['V_mps']) for r in frontier})
    target=next(r for r in rebuilt_rect if r['R_km']==.25 and r['V_mps']==.05)
    compare('decision_target',decision['primary_target_rectangle'],aggregate([p for p in points if abs(p['delta_r_km'])<=.25 and abs(p['delta_v_mps'])<=.05]))
    check('target_count',target['n_points']==150)
    wanted='B1_CONDITIONAL_DEPTH_IDENTIFIABILITY_ESTABLISHED_WITH_RANGE_SPEED_ENVELOPE' if target['all_points_pass'] else 'B1_HORIZONTAL_CONDITIONING_ENVELOPE_BELOW_TARGET'
    check('scientific_decision',decision['scientific_decision']==wanted)
    check('proposed_progress',decision['proposed_R4_B1_progress_percent']==decision['proposed_R4_overall_progress_percent']==(10 if target['all_points_pass'] else 0))
    compare('all_points_decision',decision['full_lattice'],aggregate(points))
    config=read('B1A_FREQUENCY_CONFIGURATION_SUMMARY.csv'); check('config_count',len(config)==2)
    for row in config:
        group=[p for p in points if p['frequency_configuration']==row['frequency_configuration']]
        compare('config:'+row['frequency_configuration'],row,aggregate(group)); compare('config_target:'+row['frequency_configuration'],row,aggregate([p for p in group if abs(p['delta_r_km'])<=.25 and abs(p['delta_v_mps'])<=.05]),'target_')
    # Independent exhaustive check of non-monotonic point pairs, including axis-fixed subsets.
    witnesses=set(); lookup={(p['case_id'],p['delta_r_km'],p['delta_v_mps']):p for p in points}
    for small in points:
        if small['practical_two_bin_pass']: continue
        for large in points:
            if small['case_id']!=large['case_id'] or not large['practical_two_bin_pass']: continue
            a,b=small['delta_r_km'],small['delta_v_mps']; c,e=large['delta_r_km'],large['delta_v_mps']
            if a*c<0 or b*e<0: continue
            if abs(a)<=abs(c) and abs(b)<=abs(e) and (abs(a)<abs(c) or abs(b)<abs(e)): witnesses.add((small['case_id'],a,b,c,e))
    saved=read('B1A_NON_MONOTONIC_WITNESSES.csv'); stored={(r['case_id'],float(r['smaller_dr_km']),float(r['smaller_dv_mps']),float(r['larger_dr_km']),float(r['larger_dv_mps'])) for r in saved}
    check('all_nonmonotonic_witnesses',len(saved)==len(stored) and stored==witnesses)
    check('nonmonotonic_decision',decision['non_monotonic_response']==bool(witnesses) and decision['non_monotonic_witness_count']==len(witnesses))
    check('fixed_axis_witness_count',decision['fixed_other_axis_witness_count']==sum(a==c or b==e for _,a,b,c,e in witnesses))
    for row in saved:
        small=lookup[row['case_id'],float(row['smaller_dr_km']),float(row['smaller_dv_mps'])]; large=lookup[row['case_id'],float(row['larger_dr_km']),float(row['larger_dv_mps'])]
        check('witness_values:'+str((row['case_id'],row['smaller_dr_km'],row['smaller_dv_mps'],row['larger_dr_km'],row['larger_dv_mps'])),float(row['smaller_error_m'])==small['absolute_depth_error_m'] and int(row['smaller_true_rank'])==small['true_depth_rank'] and float(row['larger_error_m'])==large['absolute_depth_error_m'] and int(row['larger_true_rank'])==large['true_depth_rank'] and bool_field(row['fixed_other_axis'])==(small['delta_r_km']==large['delta_r_km'] or small['delta_v_mps']==large['delta_v_mps']))
    reconstruction=read('B1A_NUMERICAL_RECONSTRUCTION.csv'); check('all_reconstruction_rows',len(reconstruction)==10206)
    score_lookup={(key[0],key[1],int(r['z_label_m'])):float(r['J_z_db']) for key,rows in groups.items() for r in rows}
    check('unique_complete_reconstruction_keys',len({(r['case_id'],r['conditioning_id'],int(r['z_label_m'])) for r in reconstruction})==10206 and {(r['case_id'],r['conditioning_id'],int(r['z_label_m'])) for r in reconstruction}==set(score_lookup))
    errors=[]; repeats=[]
    for row in reconstruction:
        key=row['case_id'],row['conditioning_id'],int(row['z_label_m'])
        primary=float(row['J_primary_db']); repeat=float(row['J_repeat_db']); cold=float(row['J_independent_db']); rd=abs(primary-repeat); cd=abs(primary-cold)
        check('recomputed_score:'+str(key),primary==score_lookup[key] and max(rd,cd)<=d['numerical_tolerance_rule']['maximum_allowed_absolute_difference_db'])
        check('recomputed_residuals:'+str(key),rd==float(row['absolute_repeat_difference_db']) and cd==float(row['absolute_independent_difference_db']))
        errors.append(cd); repeats.append(rd)
    floor=d['numerical_tolerance_rule']['floating_floor_multiplier']*2.220446049250313e-16*max(1,max(score_lookup.values()))
    check('tolerance_rule',math.isclose(tol,max(floor,d['numerical_tolerance_rule']['residual_multiplier']*max(errors+repeats+[d['accepted_B1_reconstruction_residual_db']])),rel_tol=1e-14))
    check('maximum_residual',validation['maximum_independent_difference_db']==max(errors) and validation['maximum_repeat_difference_db']==max(repeats))
    for path,digest in d['bindings'].items(): check('frozen:'+path,sha(ROOT/path)==digest)
    for path,digest in d['protected_historical_files'].items(): check('protected:'+path,sha(ROOT/path)==digest)
    for name in ('R4_PLAN.md','R4_EVIDENCE_LEDGER.csv'):
        path='results/R4_MASTER/'+name; parent=subprocess.check_output(['git','show',d['parent_sha']+':'+path],cwd=ROOT).replace(b'\r\n',b'\n')
        check('append_only:'+name,(ROOT/path).read_bytes().replace(b'\r\n',b'\n').startswith(parent))
    check('no_invalid_execution_marker',not (OUT/'B1A_EXECUTION_INVALID.json').exists())
    return checks

def finalize():
    if (OUT/'B1A_RESULTS_FREEZE.json').exists(): raise RuntimeError('Result already sealed')
    verify_saved(); d=js('B1A_DESIGN_FREEZE.json'); result=js('B1A_DECISION.json'); v=js('B1A_VALIDATION.json')
    xml=Path('D:/ProjectStorage/WSL-CZ/b1a-20261005/tests.xml'); tree=ET.parse(xml).getroot(); suites=list(tree.iter('testsuite'))
    tests=sum(int(s.attrib['tests']) for s in suites); errors=sum(int(s.attrib['failures'])+int(s.attrib['errors']) for s in suites)
    assert tests>0 and errors==0
    (OUT/'B1A_TESTS.xml').write_bytes(xml.read_bytes()); write('B1A_TESTS.json',{'passed':tests,'failed':0})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    figures=OUT/'figures'; figures.mkdir(exist_ok=True)
    points=read('B1A_POINT_METRICS.csv'); joint=read('B1A_JOINT_LATTICE_SUMMARY.csv')
    def save(fig,name):
        fig.tight_layout(); fig.savefig(figures/(name+'.svg')); fig.savefig(figures/(name+'.png'),dpi=150); plt.close(fig)
        path=figures/(name+'.svg'); path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n',encoding='utf-8')
    for axis,name in [('range','FIG_A_RANGE_ONLY_DEPTH_ERROR'),('speed','FIG_B_SPEED_ONLY_DEPTH_ERROR')]:
        fig,ax=plt.subplots(figsize=(9,5))
        for case in d['cases']:
            rows=[p for p in points if p['case_id']==case['case_id'] and float(p['delta_v_mps'] if axis=='range' else p['delta_r_km'])==0]
            field='delta_r_km' if axis=='range' else 'delta_v_mps'; rows.sort(key=lambda r:float(r[field]))
            ax.plot([float(r[field]) for r in rows],[float(r['absolute_depth_error_m']) for r in rows],marker='o',linestyle='none',label=case['case_id'])
        ax.axhline(10,color='black',linestyle='--',label='10 m criterion'); ax.set_xlabel('Range conditioning offset (km)' if axis=='range' else 'Speed conditioning offset (m/s)'); ax.set_ylabel('Absolute depth-label error (m)'); ax.set_title('Tested axis nodes only; no interpolation'); ax.grid(alpha=.25); ax.legend(fontsize=8,ncol=2); save(fig,name)
    dr=d['range_offsets_km']; dv=d['speed_offsets_mps']; lookup={(float(p['delta_r_km']),float(p['delta_v_mps'])):p for p in joint}
    for field,name,title in [('worst_depth_error_m','FIG_C_JOINT_WORST_DEPTH_ERROR','Max depth error across six configurations (m)'),('worst_true_depth_rank','FIG_D_JOINT_WORST_TRUE_DEPTH_RANK','Worst true-depth rank across six configurations')]:
        values=np.array([[float(lookup[r,s][field]) for s in dv] for r in dr])
        fig,ax=plt.subplots(figsize=(10,8)); im=ax.imshow(values,origin='lower',interpolation='nearest',aspect='auto')
        ax.set_xticks(range(9),[f'{x:g}' for x in dv],rotation=45); ax.set_yticks(range(9),[f'{x:g}' for x in dr]); ax.set_xlabel('Speed conditioning offset (m/s)'); ax.set_ylabel('Range conditioning offset (km)'); ax.set_title(title+'\nEach block is one tested node; lattice-index spacing, no interpolation')
        for i in range(9):
            for j in range(9): ax.text(j,i,format(values[i,j],'.0f'),ha='center',va='center',fontsize=8,color='white' if values[i,j]<values.max()/2 else 'black')
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((1.5,1.5),5,5,fill=False,edgecolor='red',linewidth=2,label='Frozen target: |dr|<=0.25, |dv|<=0.05'))
        ax.legend(loc='upper center',bbox_to_anchor=(.5,-.18),fontsize=8); fig.colorbar(im,ax=ax); save(fig,name)
    table='| Axis | Absolute offset | Exact | <=5 m | <=10 m | Rank<=3 | Practical pass | Worst error (m) | Min second margin (dB) |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n'
    for row in result['range_axis_summary']+result['speed_axis_summary']:
        table+=f"| {row['axis']} | {row['absolute_offset']} | {row['n_exact']}/12 | {row['n_within_5m']}/12 | {row['n_within_10m']}/12 | {row['n_rank_le3']}/12 | {row['n_practical_pass']}/12 | {row['worst_depth_error_m']} | {row['min_DeltaJ_second_db']:.6g} |\n"
    target=result['primary_target_rectangle']; frontier=[{'R_km':r['R_km'],'V_mps':r['V_mps']} for r in result['largest_tested_robust_rectangles']]
    comparison='| Configuration | All-lattice practical pass | Target practical pass | Target failures | Target worst depth error (m) | Target worst true rank |\n|---|---:|---:|---:|---:|---:|\n'
    for c in result['frequency_configuration_comparison']:
        comparison+=f"| {c['frequency_configuration']} | {c['n_practical_pass']}/{c['n_points']} | {c['target_n_practical_pass']}/{c['target_n_points']} | {c['target_n_failed']} | {c['target_worst_depth_error_m']} | {c['target_worst_true_depth_rank']} |\n"
    report=f"""# B1A range-speed conditioning boundary

**{result['scientific_decision']}**. Execution valid, PENDING_RESEARCH_LEAD_AUDIT. Proposed B1/overall progress={result['proposed_R4_overall_progress_percent']}%; confirmed R4 remains0%. B2 stays NOT_OPENED.

Accepted parent B1 {d['parent_sha']}: exact-horizontal depth mechanism ESTABLISHED, conditioning-tolerance envelope NOT_ESTABLISHED. The B1A design was pushed first at {result['design_commit']}; full fixed dyadic Cartesian lattice was then executed once, without changing cases, score, frequency/window weighting, depth labels, numerical rules or Gate. Existing B1 observations were read byte-for-byte, not regenerated. Historical B1/R3 and the closed CZ acquisition branch remain unchanged.

Six configurations are S7 triple (201/235/283 Hz) and S8 four-line (adds338 Hz), each at180/200/220 m: three source-depth scenarios, not six independent statistical samples. Theta=0 degree, psi=5 degree fixed. Range offsets0,+/-0.125,+/-0.25,+/-0.5,+/-1 km and speed0,+/-0.025,+/-0.05,+/-0.1,+/-0.2 m/s derive from old1 km/0.2 m/s steps by1/8,1/4,1/2,1. All486 profiles and10206 depth rows are retained; no failed point or quadrant removed.

Primary criterion: depth-label error<=10 m AND exact saved-value true-depth rank<=3. Secondary strict criterion: exact true label, rank1 and positive second-best margin greater than numerical tolerance. Strict performance is diagnostic and never substitutes for the primary criterion.

{table}

Frozen primary target |dr|<=0.25 km, |dv|<=0.05 m/s includes150 profile points (25 signed lattice nodes x6 configurations). PASS requires all150, not an average. Actual practical pass={target['n_practical_pass']}/150, failures={target['n_failed']}, worst error={target['worst_depth_error_m']} m, worst rank={target['worst_true_depth_rank']}. All16 candidate rectangles, including axes and interior signed nodes, are explicitly tabulated. Largest nondominated TESTED_ROBUST_RECTANGLE(s): {json.dumps(frontier)}. A smaller passing rectangle cannot change the target decision or award10%.

{comparison}

NON_MONOTONIC_DEPTH_CONDITIONING_RESPONSE observed={result['non_monotonic_response']}; exhaustive same-case/same-closed-orthant coordinatewise smaller-failure/larger-success witnesses={result['non_monotonic_witness_count']}, including{result['fixed_other_axis_witness_count']} with the other axis fixed. Every witness is saved; sign-specific recovery islands are not smoothed into a monotonic boundary. This is a tested lattice boundary only: no spline, fitted boundary, between-node continuous safety guarantee or precise continuous tolerance claim. Four scientific figures are driven only by saved complete data. Heatmap blocks represent individual nodes on lattice indices; nonuniform physical dyadic spacing is explicitly labeled, not interpreted as uniform continuous cells.

B1 angle/heading evidence remains only the previously tested single-factor +/-0.5-degree theta and +/-1-degree psi controls (12/12 exact each). No angle sweep is repeated or expanded; these controls are not multiplied into a joint four-parameter robustness claim.

Numerical validation: repeat every full profile after fresh model reload, independently recompute all21 scores via the accepted independent struct/scalar-geometry/compensated-mode-sum/math.fsum path. Maximum cold repeat difference={v['maximum_repeat_difference_db']:.12g} dB; independent difference={v['maximum_independent_difference_db']:.12g} dB. Frozen tol=max(128*eps*max(1,max_abs_J),10*max(new reconstruction residual,accepted B1 observation/score residual)) gives{v['numerical_tolerance_db']:.12g} dB, not fitted to margins. {v['checks_passed']} runtime checks pass, including{v['historical_files_preserved']} historical file identities;630 score rows at old B1 origin and range/speed endpoints match accepted B1. {tests} unit tests pass. Separate stdlib-only saved-score/metric/rectangle/frontier/nonmonotonic/tolerance reconstruction checks provide a second audit without physical reruns.

This identifies an upstream horizontal-state conditioning requirement in a synthetic matched E0 mechanism experiment, not platform sensor accuracy, an attainable horizontal estimate, cold-start acquisition, end-to-end observed depth accuracy, environmental robustness or final project depth specification. A passing tested rectangle would remain finite node evidence, not certification of all continuous points. S7 and S8 are both required by Gate; no four-frequency-only rescue. Failure below the frozen target is preserved without a finer lattice or downgraded target. Source-depth labels150:5:250 map to nearest stored modal samples; effective source/receiver depths are saved, with no off-grid z or interpolation.

CZ acquisition route CLOSED_FOR_CERTIFIED_GLOBAL_SEARCH; A1 remains BLOCKED/NOT_COMPLETED. No C4/F4, changed enclosure/threshold/grid/budget, 21-case CZ development, SSP/bottom/receiver-depth mismatch, bearing/platform/noise sweep or Monte Carlo. B2 and A2/A3/A4/B3/B4/C/P5 NOT_OPENED. Stop after execution commit/push and wait research-lead independent audit.

See [profiles](B1A_DEPTH_PROFILES.csv), [points](B1A_POINT_METRICS.csv), [rectangles](B1A_ROBUST_RECTANGLES.csv), [decision](B1A_DECISION.json), [validation](B1A_VALIDATION.json), [range plot](figures/FIG_A_RANGE_ONLY_DEPTH_ERROR.svg), [speed plot](figures/FIG_B_SPEED_ONLY_DEPTH_ERROR.svg), [depth-error map](figures/FIG_C_JOINT_WORST_DEPTH_ERROR.svg), [rank map](figures/FIG_D_JOINT_WORST_TRUE_DEPTH_RANK.svg).
"""
    (OUT/'B1A_REPORT.md').write_text(report,encoding='utf-8'); (OUT/'GPT_SYNC.md').write_text(report,encoding='utf-8')
    master=ROOT/'results/R4_MASTER'
    with (master/'R4_PLAN.md').open('a',encoding='utf-8',newline='') as f:
        f.write('\n\n## B1A range-speed conditioning execution checkpoint\n\n'+f"{result['scientific_decision']}; PENDING_RESEARCH_LEAD_AUDIT. Design {result['design_commit']} pushed before all486 new profiles;10206 score rows retained. Frozen target +/-0.25 km, +/-0.05 m/s: {target['n_practical_pass']}/150 pass, {target['n_failed']} failures. Proposed B1/overall progress={result['proposed_R4_overall_progress_percent']}%; confirmed0%. Largest tested robust rectangle(s)={json.dumps(frontier)}; this cannot substitute for the primary target. Nonmonotonic witnesses retained. [B1A report](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_REPORT.md), [decision](../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_DECISION.json). Finite lattice, fixed theta/psi, matched environment only; no continuous safety/engineering sensor/end-to-end claim.\n\n"+'A1 NOT_COMPLETED, CZ certified-global-search route CLOSED. B2 and A2/A3/A4/B3/B4/C/P5 NOT_OPENED. No post-result lattice additions, score/tau/grid changes, SSP mismatch or numerical rerun. Stop after commit/push for independent audit. R4_PROGRESS.json remains historical and unchanged.\n')
    with (master/'R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8',newline='') as f:
        csv.writer(f,lineterminator='\n').writerow(['B1A','../R4_B1A_RANGE_SPEED_CONDITIONING_BOUNDARY/B1A_DECISION.json',result['scientific_decision'],'486 profiles; frozen range-speed target; all six frequency/depth configurations','PENDING_RESEARCH_LEAD_AUDIT',0])
    checks=verify_saved()
    with (OUT/'B1A_SAVED_EVIDENCE_RECONSTRUCTION.csv').open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['check','pass'],lineterminator='\n');writer.writeheader();writer.writerows(checks)
    write('B1A_SAVED_EVIDENCE_RECONSTRUCTION.json',{'checks_passed':len(checks),'failed':0,'physical_profiles_rerun':0,'independent_stdlib_saved_data_gate_audit':True})
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='B1A_RESULTS_FREEZE.json')
    bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files}
    for name in ('R4_PLAN.md','R4_EVIDENCE_LEDGER.csv'): bindings[(master/name).relative_to(ROOT).as_posix()]=sha(master/name)
    write('B1A_RESULTS_FREEZE.json',{'design_commit':result['design_commit'],'design_file_sha256':sha(OUT/'B1A_DESIGN_FREEZE.json'),'scientific_decision':result['scientific_decision'],'immutable_after_commit':True,'bindings':bindings})
    print(json.dumps({'scientific_decision':result['scientific_decision'],'tests_passed':tests,'saved_evidence_checks':len(checks),'sealed_bindings':len(bindings)}))

if __name__=='__main__':
    if '--verify' in __import__('sys').argv:
        seal=js('B1A_RESULTS_FREEZE.json')
        for path,digest in seal['bindings'].items(): assert sha(ROOT/path)==digest,path
        checks=verify_saved(); print(json.dumps({'sealed_verification':'PASS','saved_reconstruction_checks':len(checks),'failed':0}))
    else: finalize()
