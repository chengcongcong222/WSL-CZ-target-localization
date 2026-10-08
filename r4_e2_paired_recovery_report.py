"""Build handoff and figures from frozen once-only results. No new forward calculation."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import r4_e2_paired_recovery as p
def main():
 O=p.O;dec=p.read(O/'NUMERICAL_RECOVERABILITY_DECISION.json');val=p.read(O/'VALIDATION.json');assert val['FAIL']==0
 rows=p.base.rows
 ca=[r for r in rows(O/'FOUR_GRID_CA_CONVERGENCE.csv') if r['status']=='AVAILABLE'];fields=[r for r in rows(O/'THREE_GRID_FIELD_CONVERGENCE.csv') if r['status']=='AVAILABLE']
 cand=[r for r in rows(O/'CANDIDATE_NUMERICAL_STABILITY.csv') if r['status']=='AVAILABLE'];en=rows(O/'LOW_ENERGY_DIAGNOSTICS.csv');logs=rows(O/'PAIRED_MODAL_GENERATION.csv');md=rows(O/'MODAL_CONVERGENCE.csv')
 fine=[r for r in ca if int(r['coarse_mesh'])==80001]
 def vals(rs,key):return np.array([float(r[key]) for r in rs])
 def maximum(rs,key):return max((float(r[key]) for r in rs),default=None)
 def q(rs,key):return np.percentile(vals(rs,key),[0,50,95,100]).tolist() if rs else []
 bins=[]
 for low,high in [(0,5),(5,25),(25,100)]:
  rs=[r for r in fine if low<float(r['minimum_power_percentile'])<=high]
  bins.append(dict(power_percentile_low=low,power_percentile_high=high,rows=len(rs),CA_quantiles_min_median_p95_max=q(rs,'CA_Frobenius'),noise1_ratio_max=maximum(rs,'numerical_over_noise1')))
 summary=dict(field_comparisons={str(g)+'_'+str(h):dict(raw=q([r for r in fields if int(r['coarse_mesh'])==g],'raw_relative'),spatial=q([r for r in fields if int(r['coarse_mesh'])==g],'common_removed_relative')) for g,h in p.TRANS},CA_comparisons={str(g)+'_'+str(h):q([r for r in ca if int(r['coarse_mesh'])==g],'CA_Frobenius') for g,h in p.TRANS},CA_noise1_exceed=sum(float(r['numerical_over_noise1'])>=1 for r in fine),CA_noise5_exceed=sum(float(r['numerical_over_noise5'])>=1 for r in fine),candidate_bound_ratio_max=maximum(cand,'numerical_bound_over_min_separation'),candidate_separation_over_noise1_min=min((float(r['marginal_separation_over_1pct']) for r in cand),default=None),candidate_separation_over_noise5_min=min((float(r['marginal_separation_over_5pct']) for r in cand),default=None),array_norm_min=minimum_norm if (minimum_norm:=min((float(r['array_norm']) for r in en if 'array_norm' in r),default=None)) else None,normalization_gain_max=maximum(fine,'max_normalization_gain'),low_energy_bins=bins,elapsed_new_solver_s=sum(float(r['elapsed_s']) for r in logs if r['new_call']=='True'))
 p.dump(O/'RECOVERY_SUMMARY.json',summary)
 lowrows=[r for r in fine if r['low_signal_level_flag']=='True']
 # Every registered state remains; only descriptive bins and empirical quantiles.
 fig,ax=plt.subplots(1,2,figsize=(10,4))
 for g,h in p.TRANS:
  rr=[r for r in fields if int(r['coarse_mesh'])==g];cc=[r for r in ca if int(r['coarse_mesh'])==g]
  for ii,rs,key in [(0,rr,'common_removed_relative'),(1,cc,'CA_Frobenius')]:
   if rs:
    v=np.sort(vals(rs,key));ax[ii].plot(np.linspace(0,1,len(v)),np.maximum(v,1e-16),label=str(g)+'→'+str(h))
 for i in range(2):ax[i].set_yscale('log');ax[i].set_xlabel('Empirical fraction (correlated rows)');ax[i].legend()
 ax[0].set_ylabel('Common-removed P error');ax[1].set_ylabel('CA Frobenius error');fig.tight_layout();fig.savefig(O/'FIG_FIELD_CA_CONVERGENCE.png');plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4))
 if fine:ax.scatter(vals(fine,'minimum_power_percentile'),vals(fine,'numerical_over_noise1'),s=8,alpha=.25)
 ax.axhline(1,c='k',ls='--');ax.set_yscale('log');ax.set_xlabel('Minimum pair raw-frequency registered power percentile');ax.set_ylabel('Finest CA difference / 1% ideal noise');fig.tight_layout();fig.savefig(O/'FIG_LOW_ENERGY_SENSITIVITY.png');plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4))
 if cand:
  x=vals(cand,'separation_160001');y=vals(cand,'triangle_bound');ax.scatter(x,y,s=8,alpha=.3);lim=[min(x.min(),y.min()),max(x.max(),y.max())];ax.plot(lim,lim,'k--');ax.set_xscale('log');ax.set_yscale('log')
 ax.set_xlabel('Fixed local candidate CA separation');ax.set_ylabel('Truth + candidate CA grid difference');fig.tight_layout();fig.savefig(O/'FIG_CANDIDATE_SEPARATION.png');plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4))
 for g,h in p.TRANS:
  rs=[r for r in fields if int(r['coarse_mesh'])==g]
  if rs:
   v=np.sort(vals(rs,'raw_relative'));ax.plot(np.linspace(0,1,len(v)),np.maximum(v,1e-16),label=str(g)+'→'+str(h))
 ax.set_yscale('log');ax.set_xlabel('Empirical fraction (correlated rows)');ax.set_ylabel('Raw complex pressure relative difference');ax.legend();fig.tight_layout();fig.savefig(O/'FIG_RAW_FIELD_CONVERGENCE.png');plt.close(fig)
 dec['independent_audit']='PASS_LOCAL_COLD_RECONSTRUCTION';dec['validation_checks']=val['checks'];dec['research_lead_audit']='PENDING'
 p.dump(O/'NUMERICAL_RECOVERABILITY_DECISION.json',dec)
 text=f'''# E2成对差频高精度数值可恢复性报告

父提交{p.PARENT}。设计提交{dec['design_SHA']}。独立审计接受上一轮E2_CA_AND_RAW_FORWARD_UNSTABLE。本阶段只执行一次，结果提交后停止。

## 结论

{dec['decision']}。建议下一研究方向{dec['next']}，但新阶段执行NOT_AUTHORIZED。原E2准入FAIL_UNCHANGED；科学信息NOT_EVALUATED；R4=0%；新MC=0、新录音=0。即使有限频对稳定，也不外推569频对或完整距离/深度可辨识性。

## 冻结与执行

9组频对：150/152、150/155、150/160；198/200、200/205、195/205；244/246、244/249、240/250Hz。共13原始频率。20001/40001模态与4个已有80001频率复用。新增求解{dec['new_KRAKEN_calls']}/22次，提供器完整{dec['provider_complete']}，新增求解累计{summary['elapsed_new_solver_s']:.3f}s，总执行{dec['elapsed_s']:.3f}s。逐调用日志保留90s上限、状态和所有原始ENV/PRT/MOD/stdout。无重试、替换设置、补频或第五网格，160001不是真值。

24冻结场景、源深180/200/220m、0/600/1200s、2个8元2m间距HLA、5km横向基线未改。旧7源深与6距离额外回归记录在REGISTERED_SIX_RANGE_CONVERGENCE.csv。三次相邻网格的场误差、幅相和相位斜率在THREE_GRID_FIELD_CONVERGENCE.csv，模态数、波数、符号对齐函数差及身份不确定性在MODAL_CONVERGENCE.csv。全模式保留，无truth-guided模式筛选。有限13输出深度不能核验全水柱模态归一化积分。

三个相邻比较的最坏误差如下，保持原始复场与CA尺度分开：

| 比较 | 原始复场相对差 | 公共消除空间差 | CA Frobenius |
| --- | --- | --- | --- |
| 20001→40001 | 8.924909% | 6.630869% | 0.109270943 |
| 40001→80001 | 2.194641% | 1.614109% | 0.026711076 |
| 80001→160001 | 0.545666% | 0.401808% | 0.006659085 |

最终公共消除空间差最坏仍高于原0.2%门限。原Gate保持失败，Outcome A不替代它。全部阵列空间误差相邻细化下降；最终/前次比为0.073029–0.414221。全部CA最终/前次比0.090941–0.337144。误差下降不等于严格绝对误差界。

## 同质成对CA

每次CA来自相同网格的两个频率，不混用80001与40001。a=P高conj(P低)，u=a/||a||，CA=uuᴴ。节点独立归一化，无跨节点相干相位。原始谱样本共享时保持同一频率压力；均方描述汇总不把9个频对当作9份独立统计信息。

完整CA比较{dec['CA_comparison_rows']}行；完整状态频对{dec['complete_CA_state_pairs']}，两次相邻细化均不增长{dec['monotonically_decreasing_CA_state_pairs']}。冻结容差1e-12仅用于浮点消差，不是放宽科学Gate。最后80001→160001最坏Frobenius差{dec['finest_CA_worst']}；主子空间、阵元相对相位、单节点与非相干双节点全部记录。

最终CA差/1%理想噪声最坏{dec['finest_noise1_ratio_max']}，/5%最坏{dec['finest_noise5_ratio_max']}。超1%尺度{summary['CA_noise1_exceed']}项，超5%尺度{summary['CA_noise5_exceed']}项。proper复相对噪声采用原边际式σ sqrt(4(1−Σ|u|⁴))，此比较是有限数值可信度诊断，不是可提取性或定位性能证明。每个网格比较的最小/中位/P95/最坏误差见RECOVERY_SUMMARY.json；其中P95是误差经验分位数，不是项目定位P95。

## 整体低能量与空间形态

LOW_ENERGY_DIAGNOSTICS.csv保留所有阵列范数和功率。按相同原始频率、相同网格的432个登记状态计算中秩功率分位数及相对中位功率，不因功率低筛选样本。最小模型阵列范数{summary['array_norm_min']}；最终对照最大CA归一化增益{summary['normalization_gain_max']}。这是固定单位源模型幅值，未校准实际源级或接收声压。旧near_null_count=0不能排除整体多模态相消，本轮不再用它作无低能量状态证明。

最终CA记录中低信号描述标记{len(lowrows)}/{len(fine)}；低5%只是描述分组。分组的CA误差与理想噪声比见RECOVERY_SUMMARY.json与图。LOW_SIGNAL_LEVEL_SENSITIVITY表示低功率与误差共现；SPATIAL_MODAL_SHAPE_SENSITIVITY标签也不声称唯一根因。空间误差与整体场强可能同时敏感，未从相关分组建立因果。

相对1%/5%噪声意味着弱场处绝对噪声随信号降低，可能过于乐观；真实固定噪声底、源谱未知、阵元校准未知与观测提取都没有验证。本轮未选择SNR、删除失败或修改原准入。

## 固定局部候选

仅每场景旧8个局部扰动：±250m距离、±10m深度及4个组合。候选记录{dec['candidate_rows']}条，80001/160001排序变化{dec['candidate_rank_changes']}条，保守数值区间未排除零{dec['candidate_unresolved']}条。最坏数值界/最小分离{summary['candidate_bound_ratio_max']}；最小边际分离/1%理想尺度{summary['candidate_separation_over_noise1_min']}，/5%尺度{summary['candidate_separation_over_noise5_min']}。这些量使用节点与差频组的等权均方描述，不是标准化统计检验或独立信息累计。

有限局部分离不等于全域排除、完整距离/深度可辨识或特征可从真实信号提取。无全球候选搜索、优化或物理信息矩阵。

## 提供器、复核与停止

模态/身份/衰减登记异常{dec['mode_anomalies']}项。模式序号不明保留不确定性，无舍弃模式。前轮FIELD内部算术原因未进一步锁定，本轮不以此为修复认证。

独立复核{val['PASS']}/{val['checks']}通过，{val['FAIL']}失败。独立struct二进制解析、独立轨迹几何、全部存档状态/频率/网格冷压力求和、全部CA/噪声、完整提供器下全部192候选、功率分位数及Outcome重建均检查。冷压力最大相对差{val['cold_pressure_max_relative_difference']}，CA最大差{val['CA_max_difference']}，候选最大差{val['candidate_max_difference']}。完整性通过不改变原物理Gate失败。

停止。FULL_BAND_CERTIFICATION_REVIEW或H3均须另经研究负责人决定，不自动开展，也不开放A2/SSP/P5。所有新输出、日志与缓存位于D盘。
'''
 (O/'NUMERICAL_RECOVERABILITY_REPORT.md').write_text(text,encoding='utf-8')
 sync=dict(parent_SHA=p.PARENT,diagnostic_independent_audit='ACCEPTED_E2_CA_AND_RAW_FORWARD_UNSTABLE',paired_design_SHA=dec['design_SHA'],paired_execution_SHA='CONTAINING_RESULT_COMMIT',remote_main_SHA='VERIFY_AFTER_PUSH',**{k:v for k,v in dec.items() if k not in ['parent_SHA','design_SHA']})
 (O/'GPT_SYNC.md').write_text('# E2 paired-frequency handoff\n\n~~~json\n'+json.dumps(sync,ensure_ascii=False,indent=2)+'\n~~~\n\n结果提交后停止，下一阶段未经授权。\n',encoding='utf-8')
 master=Path('results/R4_MASTER')
 with (master/'R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\nE2 paired recovery completed: '+dec['decision']+'. Original admission FAIL_UNCHANGED; information NOT_EVALUATED; R4=0%. Suggested '+dec['next']+' only; new stage NOT_AUTHORIZED. STOP.\n')
 with (master/'R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8') as f:f.write('E2-PAIRED-RECOVERY,../R4_E2_PAIRED_FREQUENCY_RECOVERY/NUMERICAL_RECOVERABILITY_DECISION.json,'+dec['decision']+',9 frozen pairs; bounded provider; independent reconstruction,PENDING_RESEARCH_LEAD_AUDIT;STOP;NO_CREDIT,0\n')
 p.dump(master/'R4_E2_PAIRED_RECOVERY_PROGRESS.json',dict(stage='R4_E2_PAIRED_FREQUENCY_NUMERICAL_RECOVERABILITY',status=dec['decision'],integrity='PASS',original_E2_admission='FAIL_UNCHANGED',E2_information='NOT_EVALUATED',R4_percent=0,next_suggestion=dec['next'],next_stage_execution='NOT_AUTHORIZED'))
 manifest={}
 for path in sorted(O.rglob('*')):
  if not path.is_file() or path.name=='EXECUTION_ARTIFACT_MANIFEST.json':continue
  raw=path.suffix in ['.npz','.mod','.prt','.png'] or path.name.endswith('.stdout.txt')
  manifest[str(path)]=dict(kind='RAW' if raw else 'LF',sha256=p.sha(path,raw))
 for name in ['r4_e2_paired_recovery.py','r4_e2_paired_recovery_audit.py','r4_e2_paired_recovery_report.py']:manifest[name]=dict(kind='LF',sha256=p.sha(name))
 p.dump(O/'EXECUTION_ARTIFACT_MANIFEST.json',dict(design_SHA=dec['design_SHA'],artifacts=manifest,R4_percent=0))
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
