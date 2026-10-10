"""Frozen metric definitions, terminal envelopes, figures, and research decision."""
import json,time
import numpy as np,pandas as pd
from scipy.stats import norm
from threadpoolctl import threadpool_limits
import hla_h2_core as c
PARAM_LABELS=('range_relative','bearing_deg','speed_relative','heading_deg')
REFERENCES=np.array([.10,1.,.10,5.])
def wilson(success,n):
    if n==0:return np.nan,np.nan
    p=success/n;z=norm.ppf(.975);den=1+z*z/n
    center=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return center-half,center+half
def statistics(results):
    rows=[]
    for (g,method,sigma),group in results.groupby(['geometry','method','sigma'],dropna=False):
        evaluated=group[group.status!='NOT_EVALUATED'];accepted=evaluated[evaluated.status=='ACCEPTED'];n=len(evaluated);k=len(accepted)
        low,high=wilson(k,n)
        r=dict(geometry=int(g),method=method,sigma=float(sigma),planned=len(group),evaluated=n,accepted=k,abstain=n-k,output_rate=k/n if n else np.nan,Wilson95_low=low,Wilson95_high=high,true_compatible=int(evaluated.true_compatible.sum()) if n else 0,search_failures=int(evaluated.search_failure.sum()) if n else 0,solver_exception_configurations=int((evaluated.solver_exceptions>0).sum()) if n else 0,nonconverged_terminal_runs=int(evaluated.solver_nonconverged.sum()) if n else 0)
        for name,label in zip(c.ERROR_NAMES,PARAM_LABELS):
            conditional=accepted[name].to_numpy(float)
            unconditional=evaluated[name].to_numpy(float)
            r[label+'_accepted_n']=k
            for prefix,v in [('conditional',conditional),('unconditional',unconditional)]:
                for quant,p in [('median',.5),('P90',.90),('P95',.95),('max',1.)]:r[f'{label}_{prefix}_{quant}']=c.nearest_rank(v,p)
        rows.append(r)
    return pd.DataFrame(rows)
def signal_rows(stats,method,sigma):
    output=[]
    for g in range(8):
        p=stats[(stats.geometry==g)&(stats.method==method)&(stats.sigma==sigma)]
        b=stats[(stats.geometry==g)&(stats.method=='B')]
        if not len(p) or not len(b):continue
        u=p.iloc[0];base=b.iloc[0]
        improvements={}
        qualifies=[]
        for label in ('range_relative','speed_relative','heading_deg'):
            a=u[label+'_unconditional_median'];v=base[label+'_unconditional_median']
            if np.isfinite(v) and v>0 and np.isfinite(a):
                reduction=1-a/v;passing=reduction>=.30
            elif np.isinf(v):
                reduction=np.nan;passing=bool(np.isfinite(u[label+'_unconditional_P90']) and u[label+'_unconditional_P90']<base[label+'_unconditional_P90'])
            else:reduction=np.nan;passing=False
            improvements[label+'_median_improvement']=reduction;improvements[label+'_qualifies']=passing;qualifies.append(passing)
        output.append(dict(geometry=g,method=method,sigma=sigma,complete=bool(u.evaluated==32 and base.evaluated==32),output_rate=u.output_rate,signal=bool(u.evaluated==32 and base.evaluated==32 and u.output_rate>=.9 and any(qualifies)),**improvements))
    return pd.DataFrame(output)
def terminal_envelopes(terminals,results,hidden):
    rows=[]
    for cid,group in terminals.groupby('config_id',sort=False):
        good=group[group['exception'].isna()];accepted=good[good.accepted_profile]
        result=results[results.config_id==cid].iloc[0]
        meta=dict(config_id=int(cid),geometry=int(result.geometry),replicate=int(result.replicate),method=result.method,sigma=float(result.sigma),scope='FINITE_TERMINAL_ENVELOPE',accepted_records=len(accepted))
        if not len(accepted):
            rows.append(dict(meta,clusters=0,distinct_expansion_branches=0));continue
        # Retain both original accepted b and complete profiled b, source records stay in raw file.
        candidates=[]
        for _,t in accepted.iterrows():
            h=t[list(c.AXES)].to_numpy(float)
            candidates.append((np.r_[h,t.profile_b],float(t.J_beta+t.J_profile),int(t.terminal_id)))
            if bool(t.accepted_original) and abs(t.b_original-t.profile_b)>0:
                candidates.append((np.r_[h,t.b_original],float(t.J_beta+t.J_folded),int(t.terminal_id)))
        candidates.sort(key=lambda v:(v[1],v[2]))
        reps=[]
        for p,cost,tid in candidates:
            if all(np.max(abs(p-r)/c.DEDUP)>1 for r in reps):reps.append(p)
        points=np.array([v[0] for v in candidates]);true=hidden['horizontal_truth'][int(result.geometry)]
        errors=c.error_vector(points[:,:4],true)
        r=dict(meta,clusters=len(reps),distinct_expansion_branches=accepted.branch_id.nunique())
        for j,ax in enumerate(c.AXES):
            r['minimum_'+ax]=float(points[:,j].min());r['maximum_'+ax]=float(points[:,j].max());r['span_'+ax]=float(np.ptp(points[:,j]))
            r['worst_terminal_'+c.ERROR_NAMES[j]]=float(errors[:,j].max())
        r['minimum_b']=float(points[:,4].min());r['maximum_b']=float(points[:,4].max());r['span_b']=float(np.ptp(points[:,4]))
        r['farthest_wrong_terminal_grid_scale']=float(np.linalg.norm((points[:,:4]-true)/c.DEDUP[:4],axis=1).max())
        rows.append(r)
    return pd.DataFrame(rows)
def format_value(x,scale=1):
    if pd.isna(x):return '未评价'
    if np.isinf(x):return 'INF'
    return f'{float(x)*scale:.4g}'
def report():
    validation=json.loads((c.OUT/'VALIDATION.json').read_text())
    manifest=json.loads((c.OUT/'EXECUTION_MANIFEST.json').read_text())
    results=pd.read_csv(c.OUT/'CONFIGURATION_RESULTS.csv');stats=statistics(results)
    stats.to_csv(c.OUT/'METRICS_BY_GEOMETRY.csv',index=False,float_format='%.17g')
    signals=[]
    for method in ('S','U'):
        for sigma in c.SIGMAS:signals.append(signal_rows(stats,method,sigma))
    signals.extend([signal_rows(stats,'U-mismatch',.05),signal_rows(stats,'UB-matched',.05)])
    signals=pd.concat(signals,ignore_index=True);signals.to_csv(c.OUT/'RESEARCH_SIGNAL_BY_GEOMETRY.csv',index=False)
    precision=[]
    for _,row in stats[stats.method.isin(['S','U'])].iterrows():
        for i,label in enumerate(PARAM_LABELS):
            precision.append(dict(geometry=int(row.geometry),method=row.method,sigma=float(row.sigma),parameter=label,reference=REFERENCES[i],unconditional_P95=row[label+'_unconditional_P95'],conditional_P95=row[label+'_conditional_P95'],accepted_n=int(row.accepted),planned=int(row.planned),output_rate=row.output_rate,measured_bin_attains_reference=bool(row.evaluated==32 and row[label+'_unconditional_P95']<=REFERENCES[i]),exploratory=True,no_between_bin_interpolation=True))
    pd.DataFrame(precision).to_csv(c.OUT/'PRECISION_REQUIREMENT_TABLE.csv',index=False,float_format='%.17g')
    stats[stats.method.isin(['S','U','U-mismatch','UB-matched'])].to_csv(c.OUT/'SIGN_AND_BIAS_BOUNDARY.csv',index=False,float_format='%.17g')
    hidden=np.load(c.OUT/'GENERATOR_TRUTH_AND_INNOVATIONS.npz')
    terminals=pd.read_csv(c.OUT/'ALL_TERMINALS.csv.gz')
    envelopes=terminal_envelopes(terminals,results,hidden);envelopes.to_csv(c.OUT/'TERMINAL_ENVELOPE_BY_CONFIGURATION.csv',index=False,float_format='%.17g')
    conditions=[]
    for (method,sigma),group in stats.groupby(['method','sigma']):
        r=dict(method=method,sigma=float(sigma),geometries=len(group),evaluated=int(group.evaluated.sum()),accepted=int(group.accepted.sum()),abstain=int(group.abstain.sum()),equal_scene_mean_output_rate=float(group.output_rate.mean()),minimum_scene_output_rate=float(group.output_rate.min()),worst_scene_abstain_rate=float(1-group.output_rate.min()))
        for label in PARAM_LABELS:
            r[label+'_equal_scene_mean_unconditional_median']=float(group[label+'_unconditional_median'].mean())
            r[label+'_worst_scene_unconditional_P95']=float(group[label+'_unconditional_P95'].max())
            r[label+'_equal_scene_mean_conditional_median']=float(group[label+'_conditional_median'].mean())
            r[label+'_worst_scene_conditional_P95']=float(group[label+'_conditional_P95'].max())
        conditions.append(r)
    summary=pd.DataFrame(conditions);summary.to_csv(c.OUT/'SUMMARY_BY_CONDITION.csv',index=False,float_format='%.17g')
    primary=signals[(signals.method=='U')&(signals.sigma==.05)]
    upper=signals[(signals.method=='S')&(signals.sigma==.05)]
    gain=int(primary.signal.sum())
    classifications=[]
    if gain>=5:classifications.append('HLA_H2_UNSIGNED_RADIAL_FEATURE_VALUE_SUPPORTED_CONDITIONALLY')
    else:classifications.append('HLA_H2_BENEFIT_LIMITED_AT_TESTED_PRECISIONS')
    if gain<5 and int(upper.signal.sum())>=5:classifications.append('HLA_H2_SIGN_INFORMATION_REQUIRED')
    losses=[]
    for g in range(8):
        u=stats[(stats.geometry==g)&(stats.method=='U')&(stats.sigma==.05)].iloc[0]
        b=stats[(stats.geometry==g)&(stats.method=='UB-matched')].iloc[0]
        loss=u.output_rate-b.output_rate>=.1
        for label in ('range_relative','speed_relative','heading_deg'):
            base=u[label+'_unconditional_median'];new=b[label+'_unconditional_median']
            if np.isfinite(base) and base>0:loss=loss or bool(new>=1.30*base)
        losses.append(bool(loss))
    if sum(losses)>=5:classifications.append('HLA_H2_RADIAL_REFERENCE_CONDITION_CRITICAL')
    velocity_count=int(primary.speed_relative_qualifies.sum())
    ustat=stats[(stats.method=='U')&(stats.sigma==.05)]
    if velocity_count>=5 and (ustat.range_relative_unconditional_P95>.1).any():classifications.append('HLA_H2_VELOCITY_GAIN_WITH_RANGE_LIMIT')
    if results.search_failure.fillna(False).any():classifications.append('HLA_H2_CONTINUOUS_SEARCH_LIMITED')
    if not validation['PASS']:classifications=['HLA_H2_IMPLEMENTATION_OR_INPUT_INVALID']
    if manifest['status']!='COMPLETE':classifications=['HLA_H2_PARTIAL_EXECUTION']
    decision=dict(stage='HLA-H2',classifications=classifications,primary='U / sigma_q .05 / zero bias',primary_qualifying_geometries=gain,signed_control_qualifying_geometries=int(upper.signal.sum()),radial_reference_loss_geometries=sum(losses),velocity_gain_geometries=velocity_count,INJECTED_FEATURE_ONLY=True,actual_receive_extraction='NOT_EVALUATED',COLD_START_ACOUSTIC_CAPABILITY='NOT_ESTABLISHED',continuous_global_certification=False,finite_envelope='FINITE_TERMINAL_ENVELOPE_ONLY',source_frequency_issue='E1_G0_FAIL_UNCHANGED',R4_percent=0,independent_research_lead_audit='PENDING',STOP=True)
    c.write_json(c.OUT/'DECISION.json',decision)
    table=['# HLA-H2四参数与特征精度需求','',
       '注入特征级试验，未从实际HLA信号提取径向量。以下为8个几何等权平均的无条件中位误差／最坏几何探索性P95；拒判记INF。各几何32次P95取31st顺序统计量。条件输出统计另在METRICS_BY_GEOMETRY.csv中保留。','',
       '| 条件 | 距离 % | 方位 ° | 速度 % | 航向 ° | 接受／配置 | 最低场景输出率 |',
       '|---|---:|---:|---:|---:|---:|---:|']
    for method,sigma in [('B',0),('U',.02),('U',.05),('U',.1),('U',.2),('S',.05),('U-mismatch',.05),('UB-matched',.05)]:
        r=summary[(summary.method==method)&(summary.sigma==sigma)].iloc[0]
        values=[]
        for label,scale in zip(PARAM_LABELS,(100,1,100,1)):
            values.append(format_value(r[label+'_equal_scene_mean_unconditional_median'],scale)+'／'+format_value(r[label+'_worst_scene_unconditional_P95'],scale))
        table.append('| '+method+f' {sigma:g}'+' | '+' | '.join(values)+f" | {int(r.accepted)}/{int(r.evaluated)} | {r.minimum_scene_output_rate:.2%} |")
    table+=['','剩余符号展开、连续终点聚类及全部参数包络见TERMINAL_ENVELOPE_BY_CONFIGURATION.csv。聚类尺度仅为重复结果去重尺度，不是物理可分辨尺度；所有原始终点保留。','',
       f"U主档达到登记描述性信号：{gain}/8几何；signed较强对照：{int(upper.signal.sum())}/8。共同零点未知明显损失：{sum(losses)}/8（预登记比较：输出率下降≥10pp或任一r/v/ψ无条件中位误差增加≥30%）。",
       '', '量测兼容概率属于理想连续集合：正确模型的β/q各0.975事件加union bound至少0.95；B单独0.975。有限128起点的导出候选不继承该覆盖保证。abs观测来自折叠高斯，不是幅值加白高斯；此残差与分支最小二乘也不是精确folded-normal MLE。','',
       '参考线是否达到逐参数逐档位记录，不在0.02/0.05/0.10/0.20之间插值临界精度。源相位、频率参考、相干窗、相速度、斜距映射及共同接收噪声均未认证。所有改善仅限原物理域、精确导航、带符号方位、CV目标与固定平台运动。','',
       'R4=0%；本轮STOP。建议仅根据本表决定是否投入实际特征提取；没有启动接收信号处理或增强架构。']
    (c.OUT/'FOUR_PARAMETER_RESULTS.md').write_text('\n'.join(table)+'\n',encoding='utf-8')
    # Figure coordinates mark INF explicitly above finite data, never silently drop it.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def draw(ax,x,y,**kw):
        y=np.asarray(y,float);valid=np.isfinite(y);top=max([float(v) for v in y[valid]]+[.01])
        ax.plot(np.asarray(x)[valid],y[valid],**kw)
        if (~valid).any():
            ax.scatter(np.asarray(x)[~valid],np.full((~valid).sum(),top*1.15),marker='^',color=kw.get('color','black'))
            for xx in np.asarray(x)[~valid]:ax.annotate('INF',(xx,top*1.15),fontsize=6)
    def footer(fig,methods,sigmas):
        lines=[]
        for method in methods:
            for sigma in sigmas:
                g=stats[(stats.method==method)&(stats.sigma==sigma)].sort_values('geometry')
                if len(g):lines.append(method+f' {sigma:g}: '+', '.join(f"G{int(r.geometry)+1:02} A/R/N={int(r.accepted)}/{int(r.abstain)}/{int(r.planned-r.evaluated)}" for _,r in g.iterrows()))
        fig.text(.02,.01,'A=accepted R=abstain N=not evaluated; all denominators 32 per scene\n'+'\n'.join(lines),fontsize=7,va='bottom',family='monospace')
        fig.tight_layout(rect=[0,.15 if len(lines)<4 else .20,1,.96])
    fig,axs=plt.subplots(2,2,figsize=(12,10))
    for pi,(ax,label,scale) in enumerate(zip(axs.flat,PARAM_LABELS,(100,1,100,1))):
        for g in range(8):
            group=stats[(stats.method=='U')&(stats.geometry==g)].sort_values('sigma')
            draw(ax,group.sigma,group[label+'_unconditional_P95']*scale,label=f'G{g+1:02}',color=plt.cm.tab10(g),marker='o')
        ax.axhline(REFERENCES[pi]*scale,color='black',linestyle='--',label='reference');ax.set_title(label+' exploratory unconditional P95');ax.set_xlabel('sigma_q m/s');ax.grid(alpha=.2)
    axs.flat[0].legend(ncol=3,fontsize=8);fig.suptitle('INJECTED FEATURE ONLY; all 8 scenes; INF = abstain included; no threshold interpolation');footer(fig,['U'],c.SIGMAS);fig.savefig(c.OUT/'PRECISION_VS_FOUR_PARAMETER_ERRORS.png',dpi=160);plt.close(fig)
    def footer(fig,methods,sigmas):
        lines=[]
        for method in methods:
            for sigma in sigmas:
                g=stats[(stats.method==method)&(stats.sigma==sigma)].sort_values('geometry')
                if len(g):lines.append(method+f' {sigma:g}: '+', '.join(f"G{int(r.geometry)+1:02} A/R/N={int(r.accepted)}/{int(r.abstain)}/{int(r.planned-r.evaluated)}" for _,r in g.iterrows()))
        fig.text(.02,.01,'A=accepted R=abstain N=not evaluated; all denominators 32 per scene\n'+'\n'.join(lines),fontsize=7,va='bottom',family='monospace')
        fig.tight_layout(rect=[0,.15 if len(lines)<4 else .20,1,.96])
    fig,axs=plt.subplots(2,2,figsize=(12,10))
    for ax,label,scale in zip(axs.flat,PARAM_LABELS,(100,1,100,1)):
        for method,col in [('S','#009e73'),('U','#0072b2')]:
            group=stats[(stats.method==method)&(stats.sigma==.05)].sort_values('geometry')
            draw(ax,group.geometry+1,group[label+'_unconditional_median']*scale,label=method,color=col,marker='o')
        ax.set_title(label+' unconditional median');ax.set_xlabel('all 8 development geometries');ax.grid(alpha=.2)
    axs.flat[0].legend();fig.suptitle('Known sign versus unsigned: 0.05 m/s; injected calibrated features only');footer(fig,['S','U'],[.05]);fig.savefig(c.OUT/'SIGNED_UNSIGNED_BOUNDARY.png',dpi=160);plt.close(fig)
    def footer(fig,methods,sigmas):
        lines=[]
        for method in methods:
            for sigma in sigmas:
                g=stats[(stats.method==method)&(stats.sigma==sigma)].sort_values('geometry')
                if len(g):lines.append(method+f' {sigma:g}: '+', '.join(f"G{int(r.geometry)+1:02} A/R/N={int(r.accepted)}/{int(r.abstain)}/{int(r.planned-r.evaluated)}" for _,r in g.iterrows()))
        fig.text(.02,.01,'A=accepted R=abstain N=not evaluated; all denominators 32 per scene\n'+'\n'.join(lines),fontsize=7,va='bottom',family='monospace')
        fig.tight_layout(rect=[0,.15 if len(lines)<4 else .20,1,.96])
    fig,axs=plt.subplots(2,2,figsize=(12,10))
    for ax,label,scale in zip(axs.flat,PARAM_LABELS,(100,1,100,1)):
        for method,col in [('U','#0072b2'),('U-mismatch','#cc6677'),('UB-matched','#e69f00')]:
            group=stats[(stats.method==method)&(stats.sigma==.05)].sort_values('geometry')
            draw(ax,group.geometry+1,group[label+'_unconditional_median']*scale,label=method,color=col,marker='o')
        ax.set_title(label+' unconditional median');ax.set_xlabel('all 8 development geometries');ax.grid(alpha=.2)
    axs.flat[0].legend();fig.suptitle('Zero-point boundary: clean U; +0.10 m/s bias ignored / unknown in [-0.20,+0.20]');footer(fig,['U','U-mismatch','UB-matched'],[.05]);fig.savefig(c.OUT/'RADIAL_ZERO_POINT_BOUNDARY.png',dpi=160);plt.close(fig)
    # Independent summary recomputation, including rank convention and denominators.
    audit_checks=0;audit_fail=0
    for _,row in stats.iterrows():
        group=results[(results.geometry==row.geometry)&(results.method==row.method)&(results.sigma==row.sigma)&(results.status!='NOT_EVALUATED')]
        for label,name in zip(PARAM_LABELS,c.ERROR_NAMES):
            for prefix,subset in [('unconditional',group),('conditional',group[group.status=='ACCEPTED'])]:
                x=sorted(float(v) for v in subset[name]);n=len(x)
                for quant,p in [('median',.5),('P90',.9),('P95',.95),('max',1.)]:
                    expected=x[max(0,int(np.ceil(p*n))-1)] if n else np.nan
                    actual=row[f'{label}_{prefix}_{quant}'];ok=bool((np.isnan(expected) and np.isnan(actual)) or expected==actual)
                    audit_fail+=not ok;audit_checks+=1
        k=int((group.status=='ACCEPTED').sum());n=len(group)
        expected=k/n if n else np.nan
        audit_fail+=not ((np.isnan(expected) and np.isnan(row.output_rate)) or expected==row.output_rate);audit_checks+=1
        if n:
            z=1.959963984540054;p=k/n;den=1+z*z/n
            bounds=((p+z*z/(2*n))-z*((p*(1-p)/n+z*z/(4*n*n))**.5))/den,((p+z*z/(2*n))+z*((p*(1-p)/n+z*z/(4*n*n))**.5))/den
            audit_fail+=not np.allclose(bounds,[row.Wilson95_low,row.Wilson95_high],rtol=0,atol=1e-14);audit_checks+=2
        audit_fail+=int(int(row.accepted)!=int((group.status=='ACCEPTED').sum()));audit_checks+=1
    # Second implementation of each registered scene-level research signal.
    for _,sig in signals.iterrows():
        u=stats[(stats.geometry==sig.geometry)&(stats.method==sig.method)&(stats.sigma==sig.sigma)].iloc[0]
        base=stats[(stats.geometry==sig.geometry)&(stats.method=='B')].iloc[0]
        good=False
        for label in ('range_relative','speed_relative','heading_deg'):
            bval=base[label+'_unconditional_median'];uval=u[label+'_unconditional_median']
            if np.isfinite(bval) and bval>0 and np.isfinite(uval):good=good or uval<=.70*bval
            elif np.isinf(bval):good=good or bool(np.isfinite(u[label+'_unconditional_P90']) and u[label+'_unconditional_P90']<base[label+'_unconditional_P90'])
        expected=bool(u.evaluated==32 and base.evaluated==32 and u.output_rate>=.9 and good)
        audit_fail+=expected!=bool(sig.signal);audit_checks+=1
    validation['summary_rebuild_checks']=audit_checks;validation['summary_rebuild_FAIL']=int(audit_fail);validation['PASS']=bool(validation['PASS'] and audit_fail==0)
    c.write_json(c.OUT/'VALIDATION.json',validation)
    if not validation['PASS']:
        decision['classifications']=['HLA_H2_IMPLEMENTATION_OR_INPUT_INVALID']
        c.write_json(c.OUT/'DECISION.json',decision)
    # Short research lead report, exact B SHA is supplied with the final return.
    main=summary[(summary.method=='U')&(summary.sigma==.05)].iloc[0]
    sync=['# HLA-H2同步','',f"parent：{c.PARENT}",'设计SHA：'+manifest['design_sha'],'执行SHA：包含本报告的Commit B；最终回传完整SHA及remote/main一致性。','',
       '原单HLA仅假设获得六个200 s净径向位移速率幅值；没有实际提取，没有声压、音频、复谱或传播求解。保留H1及E1原结论。导航精确、有符号0.1°方位、CV运动与独立特征噪声均为条件。','',
       f"完整2816配置、256组配对基础创新、{manifest['terminal_runs']:,}个优化终点。U/UB全部64符号展开，每支两起点；B17起点，S两起点。一次执行{manifest['elapsed_seconds']:.1f}s；所有拒判及不收敛记录均保留。",'',
       f"主U(0.05 m/s、零点已校准)接受{int(main.accepted)}/256，最低几何输出率{main.minimum_scene_output_rate:.2%}，最坏几何拒判率{main.worst_scene_abstain_rate:.2%}；{gain}/8几何达到登记信号。误差为各几何等权平均中位／最坏几何P95，拒判计INF：",'',
       '| 参数 | 无条件中位／P95 |','|---|---:|']
    for label,name,scale in zip(PARAM_LABELS,('距离 %','方位 °','速度 %','航向 °'),(100,1,100,1)):
        sync.append('| '+name+' | '+format_value(main[label+'_equal_scene_mean_unconditional_median'],scale)+'／'+format_value(main[label+'_worst_scene_unconditional_P95'],scale)+' |')
    sync+=['', '各已测档位U达到描述性信号的几何数：'+', '.join(f'{sigma:g}: {int(signals[(signals.method=="U")&(signals.sigma==sigma)].signal.sum())}/8' for sigma in c.SIGMAS)+'。完整条件输出P95、无条件INF和逐参数参考线见精度表；32次/几何的P95仅探索性，不在档位间插值。','',
       f"符号已知S主档信号{int(upper.signal.sum())}/8；未知偏置UB明显损失{sum(losses)}/8。U-mismatch/UB使用同批+0.10 m/s偏置观测；UB不读取真实偏置，只知道±0.20域。不能把S的符号条件借给U。",'',
       '判定：'+'；'.join(classifications)+'。','',
       f"冷复核每条终点、原/剖面代价、接受状态与全部分支；得分最大相对差{validation['max_relative_score_difference']:.3g}。梯度复核{validation['cold_Jacobian_checks']}条，汇总复算{audit_checks}项；FAIL={validation['terminal_score_FAIL']+int(audit_fail)}。外部研究负责人审计待定。",'',
       '全部包络只是有限终点包络，不是连续全域置信域。真实接收提取NOT_EVALUATED；未知发射频率、源跨时相位、200 s相干性、相速度/斜距映射和共享噪声仍有缺口。建议仅按已测价值档位决定是否投入提取试验，不宣布原HLA已达到精度。R4=0%；STOP。']
    (c.OUT/'GPT_SYNC.md').write_text('\n'.join(sync)+'\n',encoding='utf-8')
    c.write_json(c.OUT/'OUTPUT_MANIFEST.json',dict(sha256={str(p.relative_to(c.OUT)):c.sha(p) for p in c.OUT.rglob('*') if p.is_file() and p.name!='OUTPUT_MANIFEST.json'},new_code_sha256={p.name:c.sha(p) for p in c.ROOT.glob('hla_h2_*.py')},R4_percent=0))
    print((c.OUT/'GPT_SYNC.md').read_text(encoding='utf-8'))
if __name__=='__main__':
    with threadpool_limits(limits=4):report()
