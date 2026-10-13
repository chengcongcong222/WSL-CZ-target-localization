"""Read-only reporting of the frozen N1 diagnostic; no performance interpretation."""
from t2_n1_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 validation=read(OUT/'VALIDATION.json');assert validation['PASS'];d=read(OUT/'DECISION.json');h=read(OUT/'HISTORICAL_COUNT_REPRODUCTION.json');bound=read(OUT/'SCORE_BOUND_CHECKS.json');pairs=rows(OUT/'PAIR_ADVANTAGE_NUMERIC_INTERVALS.csv');probes=rows(OUT/'TEMPLATE_PROBE_RANKING.csv');designSHA=read(OUT/'EXECUTION_STARTED.json')['design_SHA'];lines=['# T2-N1新增数值证据','',f'设计SHA：`{designSHA}`。执行SHA由外层回报，避免本文自引用。',f'任务SHA256：`{TASKSHA}`。','',f"存档身份：36个复场、18个旧误差数组及15项历史来源封存核对。逐频向量{h['frequency_vectors']}，超原0.002门槛{h['individual_above_0p002']}（{h['individual_above_0p002']/h['frequency_vectors']:.4%}），组内最大值失败{h['groupwise_max_failures']}/54；最大e={h['shape_error_max']:.12g}。",'','探针是160001档案中的固定单位模板，没有新声音、复谱或物理真值。H网格探针自己匹配得分约1为代数恒等，不能称定位成功。','', '|探针|A：AdvH|A：AdvL|A：B合计|A：严格符号证书|B：AdvH|B：AdvL|B：B合计|B：严格符号证书|','|---|---:|---:|---:|---|---:|---:|---:|---|']
 for pid,_,_,_ in PROBES:
  a=next(r for r in pairs if r['pair']=='A' and r['catalogue']=='primary' and r['probe']==pid);b=next(r for r in pairs if r['pair']=='B' and r['catalogue']=='primary' and r['probe']==pid);lines.append(f"|{pid}|{float(a['Adv_H']):.7g}|{float(a['Adv_L']):.7g}|{float(a['B_advantage_sum']):.7g}|{a['sign_robust_by_strict_bound']}|{float(b['Adv_H']):.7g}|{float(b['Adv_L']):.7g}|{float(b['B_advantage_sum']):.7g}|{b['sign_robust_by_strict_bound']}|")
 lines += ['', '上表仅主目录；细目录及全部90个资源探针、20个配对比较完整保存，未事后换主角或错误距离半径。Adv正值仅表示该固定模板对照中处理阵的远距区分裕量较大，不代表工程性能。','']
 for label in ['A','B']:
  r=d['pair_analysis'][label];lines += [f"{label}分析标签："+'；'.join(f'`{t}`' for t in r['analysis_tags'])+'。',f"严格上界符号证书{r['strict_probe_certificates']}/10，直接Adv符号翻转{r['actual_advantage_sign_flips']}/10；该对两资源最强错误候选ID变化{r['most_wrong_ID_changes_in_pair_resources']}/20。"]
 allchanges=sum(r['strongest_wrong_ID_changed']=='True' for r in probes);rchanges=sum(r['strongest_wrong_range_changed']=='True' for r in probes);lines += ['',f"全部资源最强错误候选ID变化{allchanges}/90，其中距离标签变化{rchanges}/90；不能把仅theta/z标签变化一概写成错误距离跳变。",f"固定w逐候选评分界检查{bound['candidate_probe_checks']}项，全部被E涵盖（数值容差1e-10）；G界90项、配对优势界20项均通过。E是三频e的平均，不除以sqrt(8)，也不是噪声或对物理真实场的误差界。",f"独立冷复核{validation['check_rows']}行检查，FAIL={validation['FAIL']}；重建全部{validation['full_shape_error_entries']}项e、{validation['full_E_candidate_entries']}项E及{validation['full_probe_score_entries_both_models']}项两网格评分。最大e差{validation['max_shape_error_absolute_discrepancy']:.3g}，评分差{validation['max_score_absolute_discrepancy']:.3g}。",f"独立公式的最强错误候选ID近同分差异：{validation['independent_wrong_argmax_ID_differences']}；原排名保留，差异另列，不默认为ID已认证。",'', '建议：'+('先另行设计带独立物理基准的传播精度认证，暂不投入新Monte Carlo。' if d['global_ranking_sensitive'] else '即使模板排序同号，也需另行认证传播精度、CHIGH和实际观测可用性，不能恢复旧G0。'),'','原G0：`FORWARD_NUMERIC_NOT_ADMITTED`不变；T2-A/B性能：`NOT_EVALUATED`；CHIGH：`UNCERTIFIED`。本轮没有接受集合、真值保留、接收SNR、定位误差或覆盖率。R4=0%；研究负责人独立审计待定；STOP。']
 (OUT/'GPT_SYNC.md').write_text('\n'.join(lines)+'\n',encoding='utf-8',newline='\n')
 nextdoc=['# 下一研究决策（未授权执行）','', '本轮只诊断两个存档表示之间的排序敏感性。', '']
 for label in ['A','B']:
  r=d['pair_analysis'][label];nextdoc.append(f"{label}：严格证书{r['strict_probe_certificates']}/10，直接优势符号翻转{r['actual_advantage_sign_flips']}/10，最强错误候选ID变化{r['most_wrong_ID_changes_in_pair_resources']}/20。")
 nextdoc+=['','优先另行设计独立传播精度认证：需要独立物理／解析基准、原始模态存储精度与相速度窗截断证据，不能把160001当真值，不能无上限加密网格。', 'CHIGH共同遗漏不在e/E/B的覆盖范围内；必要FIELD坐标／存储误差须单独定义，原FIELDGate不改写。', '即使若干固定探针获得严格符号证书，也仅适用于这两个网格和登记模板；不能准入真实接收、噪声、增益、动态定位或硬件选型。', '', '没有N2或旧G0恢复执行授权。本輪所有新求解、MC、源噪声、录音、复谱、优化器、FIM均0，R4=0%，STOP。']
 (OUT/'NEXT_RESEARCH_DECISION.md').write_text('\n'.join(nextdoc)+'\n',encoding='utf-8',newline='\n')
 (OUT/'figures').mkdir(exist_ok=True);fig,axs=plt.subplots(1,2,figsize=(12,4.6))
 for ax,label in zip(axs,['A','B']):
  for cat,dy,color in [('primary',-.13,'#487caf'),('fine',.13,'#b57f37')]:
   rr=[next(r for r in pairs if r['pair']==label and r['catalogue']==cat and r['probe']==pid) for pid,_,_,_ in PROBES];ax.errorbar([float(r['Adv_H']) for r in rr],np.arange(5)+dy,xerr=[float(r['B_advantage_sum']) for r in rr],fmt='o',capsize=3,color=color,label=cat)
  ax.axvline(0,color='black',ls='--');ax.set_yticks(range(5),[p[0] for p in PROBES]);ax.set_title(label+': archived template advantage +/- bound');ax.set_xlabel('AdvH and cross-mesh bound (dimensionless)');ax.legend()
 fig.suptitle('Fixed templates only; not reception performance');fig.tight_layout();fig.savefig(OUT/'figures/PAIR_ADVANTAGE_BOUNDS.png',dpi=160);plt.close(fig)
 print('REPORT complete; original G0 unchanged; STOP',flush=True)
if __name__=='__main__':main()