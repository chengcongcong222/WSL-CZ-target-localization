"""Post-execution interpretation and manifest, no new research observations."""
import json,time
import numpy as np,pandas as pd
import hla_h1_core as c
from threadpoolctl import threadpool_limits
def run():
    f=pd.read_csv(c.OUT/'CONFIGURATION_RESULTS.csv');a=pd.read_csv(c.OUT/'PRIMARY_METRICS.csv')
    off=pd.read_csv(c.OUT/'OFFGRID_RESULTS.csv');truths=pd.read_csv(c.OUT/'OFFGRID_TRUTH_EVALUATION_ONLY.csv')
    logs=pd.read_csv(c.OUT/'OFFGRID_OPTIMIZER_LOG.csv')
    decision=json.loads((c.OUT/'DECISION.json').read_text());v=json.loads((c.OUT/'VALIDATION.json').read_text());full=json.loads((c.OUT/'FULL_FINITE_GRID_COLD_AUDIT.json').read_text())
    if not (v['PASS'] and full['PASS']):raise RuntimeError('Cold audit failed: cannot finalize a scientific gain conclusion')
    m=c.Model();data=np.load(c.OUT/'observations/OFFGRID.npz');diag=[]
    for g,t in truths.iterrows():
        s=t[list(c.AXES)].to_numpy(float);z=int((t.z_label-150)/5)
        truelevel=m.levels(s,True)[0,z];pred=c.geometry(s)[0][0]
        for rep in range(4):
            ci=g*4+rep;bc=np.square(c.wrap(pred-data['bearings'][ci])/np.radians(.1)).sum()
            for method in ('BEARING','M0','M1'):
                score=bc
                if method!='BEARING':
                    e=truelevel-data['levels'][ci]
                    for sl in (slice(0,61),slice(61,121)):
                        q=e[:,sl];q=q-q.mean(axis=-1,keepdims=True)
                        if method=='M1':q=q-q.mean(axis=-2,keepdims=True)
                        score+=np.square(q).sum()/.25**2
                est=off[(off.method==method)&(off.config==ci)].iloc[0]
                diag.append(dict(panel=g,replicate=rep,method=method,truth_joint_score=float(score),truth_joint_accepted=bool(score<=c.CUT[method]),search_score=float(est.score),search_accepted=bool(est.accepted_count>0),truth_used='POSTHOC_ONLY_NOT_ESTIMATOR_INPUT'))
    pd.DataFrame(diag).to_csv(c.OUT/'OFFGRID_TRUTH_SCORE_DIAGNOSTIC.csv',index=False)
    d=pd.DataFrame(diag);stats=[]
    for method,rows in off.groupby('method'):
        r=dict(method=method,configs=len(rows),accepted_configurations=int((rows.accepted_count>0).sum()),empty_search_envelope=int((rows.accepted_count==0).sum()),no_finite_outputs=int(rows.failed.sum()),optimizer_terminal_runs=len(logs[logs.method==method]),optimizer_successful_runs=int(logs[logs.method==method].success.sum()),truth_joint_accepted=int(d[d.method==method].truth_joint_accepted.sum()))
        for ax in c.AXES:
            r['median_error_'+ax]=float(rows['point_error_'+ax].median())
            r['max_error_'+ax]=float(rows['point_error_'+ax].max())
        r['median_relative_range']=float(rows.point_relative_r.median());r['max_relative_range']=float(rows.point_relative_r.max())
        r['median_relative_speed']=float(rows.point_relative_v.median());r['max_relative_speed']=float(rows.point_relative_v.max())
        stats.append(r)
    pd.DataFrame(stats).to_csv(c.OUT/'OFFGRID_SUMMARY.csv',index=False)
    ss=pd.DataFrame(stats).set_index('method')
    search_limited=bool(ss.loc['M1','accepted_configurations']<ss.loc['M1','truth_joint_accepted'])
    decision['finite_grid_classification']=decision['classification']
    decision['continuous_panel_classification']='HLA_H1_SEARCH_OR_NUMERICAL_LIMITED' if search_limited else 'FINITE_LOCAL_SEARCH_DIAGNOSTIC_ONLY'
    decision['continuous_success_not_established']=search_limited
    decision['continuous_M1_accepted_configurations']=int(ss.loc['M1','accepted_configurations'])
    decision['continuous_M1_truth_joint_accepted_posthoc']=int(ss.loc['M1','truth_joint_accepted'])
    decision['route_recommendation']='CONDITIONAL_MECHANISM_RESERVE; NOT_ESTABLISHED_CONTINUOUS_HORIZONTAL_ESTIMATOR'
    decision['full_finite_grid_direct_audit_PASS']=full['PASS']
    decision['independent_research_lead_audit']='PENDING'
    c.write_json(c.OUT/'DECISION.json',decision)
    paired=[]
    for g in range(6):
        archive=np.load(c.OUT/f'accepted/G{g+1:02}.npz')
        for sigma_index in range(3):
            for rep in range(32):
                ci=sigma_index*128+rep
                i0=archive[f'M1_{ci}_ids']
                for shift in (32,64):
                    i1=archive[f'M1_{ci+shift}_ids']
                    paired.append(dict(geometry=g+1,sigma_index=sigma_index,replicate=rep,source_shift=shift,same_set=bool(np.array_equal(i0,i1))))
    c.write_json(c.OUT/'EMPIRICAL_SOURCE_INVARIANCE.json',dict(checks=len(paired),all_common_source_sets_equal=all(x['same_set'] for x in paired),scope='Paired S0/S1/S2 finite-grid M1 sets; shared noise, not independent validation samples'))
    dimensions=[]
    for method,rows in f[(f.sigma==.25)&(f.source=='S2')].groupby('method'):
        for ax in c.AXES:
            dimensions.append(dict(method=method,parameter=ax,median_point_error=float(rows[f'point_error_{ax}'].median()),maximum_point_error=float(rows[f'point_error_{ax}'].max()),median_accepted_set_span=float(rows[f'span_{ax}'].median()) if rows[f'span_{ax}'].notna().any() else None,min_scene_median_span=float(a[(a.method==method)&(a.sigma==.25)&(a.source=='S2')][f'median_span_{ax}'].min()) if rows[f'span_{ax}'].notna().any() else None,max_scene_median_span=float(a[(a.method==method)&(a.sigma==.25)&(a.source=='S2')][f'median_span_{ax}'].max()) if rows[f'span_{ax}'].notna().any() else None))
    c.write_json(c.OUT/'FOUR_PARAMETER_PRIMARY_TABLE.json',dimensions)
    frozen=json.loads((c.OUT/'DESIGN_FREEZE.json').read_text())
    hash_pass=all(c.sha(c.ROOT/name)==digest for name,digest in frozen['sha256'].items())
    if not hash_pass:raise RuntimeError('Frozen bytes changed')
    v['full_finite_grid_audit']=full;v['frozen_input_code_hashes_unchanged']=hash_pass
    v['supplemental_structural_controls']=json.loads((c.OUT/'SUPPLEMENTAL_STRUCTURAL_CONTROLS.json').read_text())['PASS']
    v['integer_dtype_metric_limitation']='Preserved separately; all scored/generated pilot states are float64.'
    v['independent_research_lead_audit']='PENDING'
    c.write_json(c.OUT/'VALIDATION.json',v)
    manifest=json.loads((c.OUT/'EXECUTION_MANIFEST.json').read_text())
    c.write_json(c.OUT/'FORWARD_EVALUATION_ACCOUNTING.json',dict(new_KRAKEN_FIELD_BELLHOP_calls=0,new_modal_field_evaluations=True,full_grid_cache_horizontal_states=114576,full_grid_cache_shared_depth_states=114576*21,main_boundary_direct_horizontal_evaluations=int(f[f.method!='BEARING'].direct_boundary_states.sum()),main_best_direct_horizontal_evaluations=6912,offgrid_optimizer_terminal_direct_horizontal_evaluations=len(logs),offgrid_recorded_nfev=int(logs.nfev.sum()),offgrid_numerical_Jacobian_residual_calls_not_separately_instrumented=True,offgrid_residual_calls_upper_bound=5*int(logs.nfev.sum()),full_finite_grid_cold_direct_depth_states=full['direct_h_depth_states'],scope='Modal evaluations are new acoustic arithmetic, not new solver calls. nfev excludes numerical Jacobian calls; upper bound is accounting, not observed counter.'))
    md=['# 四参数结果与适用条件','',
        '主比较为0.25 dB、S1/S2。下表范围为各几何有限接受集跨度中位数的最小—最大值；空集没有跨度。宽度0只表示单个格点标签，不代表连续精度。','',
        '| 参数 | 方位基线集合跨度 | M0集合 | M1集合跨度 | M1离网格点误差：中位／最大 |',
        '|---|---:|---|---:|---:|']
    labels={'r_km':'距离 km','theta_deg':'方位 °','v_mps':'速度 m/s','psi_deg':'航向 °'}
    for ax in c.AXES:
        base=a[(a.method=='BEARING')&(a.sigma==.25)&(a.source=='S2')][f'median_span_{ax}']
        m1=a[(a.method=='M1')&(a.sigma==.25)&(a.source=='S2')][f'median_span_{ax}']
        m1off=ss.loc['M1']
        md.append(f"| {labels[ax]} | {base.min():.6g}—{base.max():.6g} | S1/S2全部空集 | {m1.min():.6g}—{m1.max():.6g} | {m1off['median_error_'+ax]:.6g}／{m1off['max_error_'+ax]:.6g} |")
    md+=['','S1和S2各自均6/6几何达到登记A/B描述性信号；M1真值保留率90.625%—100%，M0为0%。主矩阵的M1最优水平格点误差为0，但这依赖真值落在完整历史网格上。',
          '',f"离网格M1导出接受候选：{int(ss.loc['M1','accepted_configurations'])}/32；事后直接验证的真值(h,z)通过：{int(ss.loc['M1','truth_joint_accepted'])}/32。冻结有限搜索没有恢复这些可接受区域，不能解释为正确解不存在。所有误差保留，未开启第二次搜索。",
          '', 'S3在0.25 dB下M1保留率15.625%—50%。固定谱级可未知，公共时间变化可任意；非共同线变化仍构成失配。S0下M0/M1/ORACLE保留率与集合宽度见SOURCE_MODEL_TRADEOFF.csv；有限网格单标签不够量化格间信息损失。',
          '', '结论：公共源幅剖面化的条件机制可保留为备选研究线索；连续水平估计能力尚未建立。本轮不进入申请中的已验证精度主结论。名义匹配环境、可跟踪三频、精确导航、有符号0.1°方位及特征误差模型都是适用条件；没有真实接收提取、源谱稳定性、环境鲁棒性或连续覆盖认证。R4=0%。STOP。',
          '', '复现：冻结A的执行入口为hla_h1_execute.py；独立复核为hla_h1_full_grid_audit.py及hla_h1_audit.py，后者按保存深度区分原深度／重剖面记录。冻结报告原型保留，报告纠错说明见AUDIT_REPORTING_NOTE.json。']
    (c.OUT/'FOUR_PARAMETER_RESULTS.md').write_text('\n'.join(md)+'\n',encoding='utf-8')
    sync=['# HLA-H1同步','',f"parent：{frozen['parent']}",f"设计SHA：{manifest['design_sha']}",'执行SHA：包含本报告的Commit B（最终回传给出完整SHA）；remote/main核对后STOP。','',
          '有限网格：HLA_H1_SOURCE_PROFILED_HORIZONTAL_GAIN_SUPPORTED_CONDITIONALLY。',
          '离网格：HLA_H1_SEARCH_OR_NUMERICAL_LIMITED。只能保留为条件机制备选，不能声称连续估计已闭合。','',
          f"完整主矩阵9216=4×2304配置；192组基础创新配对复用。离网格96=3×32配置、3072起点；原深度与更优重剖面结果均保留。执行{manifest['elapsed_seconds']:.1f}s，峰值{manifest['peak_rss_bytes']/1024**3:.3f} GiB。新传播求解器0次；有新增存档模态声场评价。",'',
          '主比较0.25 dB，S1/S2分别6/6几何达到登记A/B信号。M1真值保留90.625%—100%；M0全部空集。M1有限网格最优点误差0，接受集为单标签；不解释成连续零误差。','']+md[4:10]+[
          '',f"M1离网格接受0/32，事后真值代价通过{int(ss.loc['M1','truth_joint_accepted'])}/32；正确区域存在但冻结搜索未命中。M0也接受0/32。离网格表中的误差为实际有限搜索最优点误差，接受包络为空，不能称作精度达标。",
          '', 'S3失配：M1保留率15.625%—50%。不对任意源谱鲁棒；S0信息代价逐场景完整保留。',
          '',f"冷复核：全有限网格6912个声学配置接受集一致，{full['direct_h_depth_states']:,}个直接(h,z)状态；误剔除0、误接受0。其余冷检查{v['checks']:,}项、0 FAIL。冻结输入/估计器哈希不变。审计报告修正了原深度记录的比较对象，未改变科学执行。整数数组误差助手限制另存，实际状态均float64。",
          '', '适用条件：名义环境匹配、原单HLA、三频可跟踪、精确导航、有符号0.1°方位、独立特征级Gaussian误差。深度收紧、实测信号提取、连续全域覆盖均未评价。R4=0%；不自动开展后续方法或增强架构。STOP。']
    (c.OUT/'GPT_SYNC.md').write_text('\n'.join(sync)+'\n',encoding='utf-8')
    # Append one row only; frozen master decisions are untouched.
    ledger=c.ROOT/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'
    stage='R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT'
    if stage not in ledger.read_text(encoding='utf-8'):
        import csv
        with ledger.open('a',encoding='utf-8',newline='') as stream:
            csv.writer(stream,lineterminator='\n').writerow([stage,'../R4_SINGLE_HLA_SOURCE_PROFILED_HORIZONTAL_PILOT/DECISION.json',decision['classification'],'Original single HLA conditional new-method pilot; finite grid positive; offgrid search not established','LOCAL_COLD_AUDIT_PASS_RESEARCH_LEAD_PENDING',0])
    size=sum(p.stat().st_size for p in c.OUT.rglob('*') if p.is_file())
    if size>512*1024**2:raise RuntimeError('Git delivery budget exceeded')
    c.write_json(c.OUT/'OUTPUT_MANIFEST.json',dict(sha256={str(p.relative_to(c.OUT)):c.sha(p) for p in c.OUT.rglob('*') if p.is_file() and p.name!='OUTPUT_MANIFEST.json'},new_code_sha256={p:c.sha(c.ROOT/p) for p in ['hla_h1_core.py','hla_h1_execute.py','hla_h1_preflight.py','hla_h1_report.py','hla_h1_audit.py','hla_h1_full_grid_audit.py','hla_h1_finalize.py']},bytes=size,R4_percent=0))
    print((c.OUT/'GPT_SYNC.md').read_text(encoding='utf-8'))
    print(ss.to_string())
if __name__=='__main__':
    with threadpool_limits(limits=4):run()
