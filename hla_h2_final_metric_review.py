"""Independent scalar output/error/quantile/decision rebuild, no optimization."""
import csv,math,json,time
from collections import defaultdict
import hla_h2_core as c
AXES=list(c.AXES)
ERR=list(c.ERROR_NAMES)
LABELS=['range_relative','bearing_deg','speed_relative','heading_deg']
def read(name):
    with (c.OUT/name).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
def q(v,p):
    a=sorted(v);return a[max(0,math.ceil(len(a)*p)-1)] if a else math.nan
def equal(a,b):
    return (math.isnan(a) and math.isnan(b)) or a==b or (math.isfinite(a) and math.isfinite(b) and abs(a-b)<=1e-12*max(1,abs(a),abs(b)))
def audit():
    truth=[]
    with (c.ROOT/'results/R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/OFFGRID_TRUTH_EVALUATION_ONLY.csv').open(newline='',encoding='utf-8') as f:
        truth=[[float(r[k]) for k in AXES] for r in csv.DictReader(f)]
    estimates={int(r['config_id']):r for r in read('ESTIMATES.csv')}
    rows=read('CONFIGURATION_RESULTS.csv');metrics=read('METRICS_BY_GEOMETRY.csv');summaries=read('SUMMARY_BY_CONDITION.csv')
    checks=0;fails=0;groups=defaultdict(list)
    for row in rows:
        cid=int(row['config_id']);est=estimates[cid];h=truth[int(row['geometry'])]
        if row['status']=='ACCEPTED':
            sh=[float(est[k]) for k in AXES]
            expected=[abs(sh[0]-h[0])/h[0],abs(math.degrees(math.atan2(math.sin(math.radians(sh[1]-h[1])),math.cos(math.radians(sh[1]-h[1]))))),abs(sh[2]-h[2])/h[2],abs(math.degrees(math.atan2(math.sin(math.radians(sh[3]-h[3])),math.cos(math.radians(sh[3]-h[3])))))]
        else:expected=[math.inf]*4
        for name,v in zip(ERR,expected):fails+=not equal(v,float(row[name]));checks+=1
        groups[(int(row['geometry']),row['method'],float(row['sigma']))].append(row)
    rebuilt={}
    for metric in metrics:
        key=(int(metric['geometry']),metric['method'],float(metric['sigma']))
        group=groups[key];accepted=[r for r in group if r['status']=='ACCEPTED']
        rebuilt[key]=dict(output_rate=len(accepted)/len(group),accepted=len(accepted))
        for keyname,value in [('planned',len(group)),('evaluated',len(group)),('accepted',len(accepted)),('abstain',len(group)-len(accepted)),('output_rate',len(accepted)/len(group))]:
            fails+=not equal(value,float(metric[keyname]));checks+=1
        for label,error in zip(LABELS,ERR):
            for prefix,subset in [('unconditional',group),('conditional',accepted)]:
                for quant,p in [('median',.5),('P90',.9),('P95',.95),('max',1.)]:
                    value=q([float(r[error]) for r in subset],p);field=f'{label}_{prefix}_{quant}'
                    rebuilt[key][field]=value;fails+=not equal(value,float(metric[field]));checks+=1
    def signal(method,sigma):
        count=0
        for g in range(8):
            b=rebuilt[g,'B',0];u=rebuilt[g,method,sigma]
            good=False
            for label in ['range_relative','speed_relative','heading_deg']:
                bv=b[label+'_unconditional_median'];uv=u[label+'_unconditional_median']
                if math.isfinite(bv) and bv>0:good|=math.isfinite(uv) and 1-uv/bv>=.30
                elif math.isinf(bv):good|=math.isfinite(u[label+'_unconditional_P90']) and u[label+'_unconditional_P90']<b[label+'_unconditional_P90']
            count+=u['output_rate']>=.9 and good
        return int(count)
    for summary in summaries:
        subset=[rebuilt[g,summary['method'],float(summary['sigma'])] for g in range(8)]
        for label in LABELS:
            for prefix in ['unconditional','conditional']:
                field=label+'_equal_scene_mean_'+prefix+'_median'
                value=sum(r[label+'_'+prefix+'_median'] for r in subset)/8
                fails+=not equal(value,float(summary[field]));checks+=1
                field=label+'_worst_scene_'+prefix+'_P95'
                value=max(r[label+'_'+prefix+'_P95'] for r in subset)
                fails+=not equal(value,float(summary[field]));checks+=1
    for row in read('PRECISION_REQUIREMENT_TABLE.csv'):
        index=LABELS.index(row['parameter']);reference=[.1,1,.1,5][index]
        value=rebuilt[int(row['geometry']),row['method'],float(row['sigma'])][row['parameter']+'_unconditional_P95']
        expected=value<=reference
        fails+=expected!=(row['measured_bin_attains_reference']=='True');checks+=1
    gain=signal('U',.05);signed=signal('S',.05);loss=0;velocity=0
    for g in range(8):
        u=rebuilt[g,'U',.05];b=rebuilt[g,'UB-matched',.05];base=rebuilt[g,'B',0]
        condition=u['output_rate']-b['output_rate']>=.10
        for label in ['range_relative','speed_relative','heading_deg']:
            old=u[label+'_unconditional_median'];new=b[label+'_unconditional_median']
            condition|=math.isfinite(old) and old>0 and new>=1.3*old
        loss+=condition
        bv=base['speed_relative_unconditional_median'];uv=u['speed_relative_unconditional_median']
        velocity+=math.isfinite(bv) and bv>0 and math.isfinite(uv) and 1-uv/bv>=.3
    expected=['HLA_H2_UNSIGNED_RADIAL_FEATURE_VALUE_SUPPORTED_CONDITIONALLY'] if gain>=5 else ['HLA_H2_BENEFIT_LIMITED_AT_TESTED_PRECISIONS']
    if gain<5 and signed>=5:expected.append('HLA_H2_SIGN_INFORMATION_REQUIRED')
    if loss>=5:expected.append('HLA_H2_RADIAL_REFERENCE_CONDITION_CRITICAL')
    if velocity>=5 and any(rebuilt[g,'U',.05]['range_relative_unconditional_P95']>.1 for g in range(8)):expected.append('HLA_H2_VELOCITY_GAIN_WITH_RANGE_LIMIT')
    if any(r['search_failure']=='True' for r in rows):expected.append('HLA_H2_CONTINUOUS_SEARCH_LIMITED')
    decision=json.loads((c.OUT/'DECISION.json').read_text())
    fails+=set(expected)!=set(decision['classifications']);checks+=1
    for name,v in [('primary_qualifying_geometries',gain),('signed_control_qualifying_geometries',signed),('radial_reference_loss_geometries',loss),('velocity_gain_geometries',int(velocity))]:
        fails+=v!=decision[name];checks+=1
    result=dict(PASS=fails==0,checks=checks,FAIL=int(fails),primary_qualifying_geometries=gain,signed_qualifying_geometries=signed,unknown_zero_loss_geometries=int(loss),scope='stdlib CSV + scalar source-panel error, order statistics, denominators, equal-scene summaries, tested reference bins, all final classifications; no fitting',estimated_features_actually_extracted=False)
    c.write_json(c.OUT/'INDEPENDENT_FINAL_METRIC_REVIEW.json',result)
    validation=json.loads((c.OUT/'VALIDATION.json').read_text());validation['final_scalar_metric_review']=result;validation['PASS']=validation['PASS'] and result['PASS']
    c.write_json(c.OUT/'VALIDATION.json',validation);print(json.dumps(result,indent=2))
    if fails:raise RuntimeError('Independent final metric review failed')
if __name__=='__main__':audit()
