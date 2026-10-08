"""Audited report only; no scientific replay. Retained frozen FD values are explicitly unadmitted."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import r4_e2_nine_pair_pilot as p
def main():
 O=p.O;d=p.read(O/'PILOT_DECISION.json');v=p.read(O/'VALIDATION.json');info=p.base.rows(O/'LOCAL_INFORMATION_BY_SCENE.csv');st=p.base.rows(O/'NUMERICAL_INFORMATION_STABILITY.csv');cand=p.base.rows(O/'CANDIDATE_SEPARATION.csv');chain=p.base.rows(O/'ANALYTIC_CHAIN_RULE_AND_RANK_AUDIT.csv');cal=p.base.rows(O/'CALIBRATION_PROFILE.csv')
 assert d['decision']=='IMPLEMENTATION_INVALID' and v['FAIL']==60
 fine=[r for r in info if r['mesh']=='160001' and r['sigma']=='0.01' and r['bearing_profiled']=='True']
 ci={(r['scene_id'],r['noise_model'],r['calibration'],r['representation']):r for r in fine}
 scenes=[s['scene_id'] for s in p.read(O/'DESIGN_FREEZE.json')['scenes']]
 def quant(values):return np.percentile(np.array(values),[0,50,100]).tolist()
 ratios={}
 for rep,raw in [('P1','P0_SINGLE'),('P2','P0_DUAL')]:
  ratios[rep+'_over_'+raw]=quant([float(ci[s,'RELATIVE','C1',rep]['information_trace'])/float(ci[s,'RELATIVE','C1',raw]['information_trace']) for s in scenes])
 losses={}
 for key in ['information_trace','range_information_depth_profiled','depth_information_horizontal_profiled']:
  losses[key]=quant([float(ci[s,'RELATIVE','C1','P2'][key])/float(ci[s,'RELATIVE','C0','P2'][key]) for s in scenes])
 floor={}
 for key in ['range_information_depth_profiled','depth_information_horizontal_profiled']:
  floor[key]=quant([float(ci[s,'ABSOLUTE_FLOOR','C1','P2'][key])/float(ci[s,'RELATIVE','C1','P2'][key]) for s in scenes])
 strength={}
 for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
  rs=[r for r in fine if r['representation']=='P2' and r['calibration']=='C1' and r['noise_model']==noise]
  strength[noise]=dict(range_information_min=min(float(r['range_information_depth_profiled']) for r in rs),depth_information_min=min(float(r['depth_information_horizontal_profiled']) for r in rs),range_depth_unknown_over_known_quantiles=quant([float(r['range_information_depth_profiled'])/float(r['range_information_exact_depth']) for r in rs]),weak_direction_median_abs=np.median(np.abs([json.loads(r['scaled_weak_direction']) for r in rs]),axis=0).tolist())
 cs=[r for r in cand if r['mesh']=='160001']
 csum={}
 for noise in ['RELATIVE','ABSOLUTE_FLOOR']:
  rr=[r for r in cs if r['noise_model']==noise]
  csum[noise]=dict(min_CA_over_noise1=min(float(r['C1_CA_Frobenius_rms'])/float(r['CA_marginal_noise1_rms']) for r in rr),min_CA_over_noise5=min(float(r['C1_CA_Frobenius_rms'])/float(r['CA_marginal_noise5_rms']) for r in rr),numerical_positive=sum(r['numerical_interval_positive']=='True' for r in rr),max_numerical_bound_over_separation=max(float(r['numerical_triangle_bound'])/float(r['C1_CA_Frobenius_rms']) for r in rr),rank_changes=sum(r['rank']!=next(t['rank'] for t in cand if t['scene_id']==r['scene_id'] and t['candidate_id']==r['candidate_id'] and t['noise_model']==noise and t['mesh']=='80001') for r in rr))
 primary=[r for r in st if r['representation']=='P2' and r['noise_model']=='RELATIVE' and r['calibration']=='C1' and r['bearing_profiled']=='True']
 numerical={k:max(float(r[k]) for r in primary) for k in ['fd_relative','mesh_relative','range_information_relative','depth_information_relative']}
 summary=dict(status='RETAINED_FROZEN_FD_DIAGNOSTICS_NOT_SCIENTIFICALLY_ADMITTED',representation_trace_ratios=ratios,C1_over_C0=losses,C2_information=0,noise_floor_over_relative=floor,retained_dual_primary_information=strength,retained_candidate_diagnostics=csum,retained_primary_numerical_maxima=numerical,unresolved_stability_rows=sum(r['status']!='NUMERICALLY_STABLE' for r in st),total_stability_rows=len(st),single_chain_rank_mismatch_count=sum(r['representation']=='P1' and r['analytic_horizontal_chain_rank']!=r['frozen_fd_rank'] for r in chain),single_depth_chain_relative_max=max(float(r['depth_information_relative_difference']) for r in chain if r['representation']=='P1'),dual_range_chain_relative_max=max(float(r['range_information_relative_difference']) for r in chain if r['representation']=='P2'),dual_depth_chain_relative_max=max(float(r['depth_information_relative_difference']) for r in chain if r['representation']=='P2'))
 p.dump(O/'AUDITED_DIAGNOSTIC_SUMMARY.json',summary)
 fig,ax=plt.subplots(1,2,figsize=(10,4))
 for rep,color in [('P1','tab:orange'),('P2','tab:blue')]:
  rr=[r for r in chain if r['representation']==rep]
  ax[0].scatter([int(r['frozen_fd_rank']) for r in rr],[int(r['analytic_horizontal_chain_rank']) for r in rr],label=rep,c=color,alpha=.4)
  vv=sorted(float(r['depth_information_relative_difference']) for r in rr);ax[1].plot(np.linspace(0,1,len(vv)),vv,label=rep,c=color)
 ax[0].plot([3.5,5.5],[3.5,5.5],'k--');ax[0].set_xlabel('Frozen finite-difference rank');ax[0].set_ylabel('Analytic chain-rule rank');ax[0].legend()
 ax[1].axhline(.1,c='k',ls='--');ax[1].set_yscale('log');ax[1].set_xlabel('Registered scene/mesh fraction');ax[1].set_ylabel('Profiled-depth discrepancy');ax[1].legend();fig.tight_layout();fig.savefig(O/'FIG_PHYSICAL_FORMULA_FAILURE.png');plt.close(fig)
 fig,ax=plt.subplots(figsize=(7,4))
 for key in ['range_information_depth_profiled','depth_information_horizontal_profiled']:
  vv=[float(ci[s,'ABSOLUTE_FLOOR','C1','P2'][key])/float(ci[s,'RELATIVE','C1','P2'][key]) for s in scenes]
  ax.plot(range(24),vv,'o-',label=key.split('_information')[0])
 ax.axhline(1,c='k',ls='--');ax.set_yscale('log');ax.set_xlabel('Frozen scene index');ax.set_ylabel('Floor / relative local information');ax.set_title('Retained FD diagnostics only — pilot IMPLEMENTATION_INVALID');ax.legend();fig.tight_layout();fig.savefig(O/'FIG_FIXED_FLOOR_SENSITIVITY.png');plt.close(fig)
 (O/'PLANE_WAVE_AND_DATA_PROCESSING_CONTROLS.md').write_text(f'''# 公式、零信息与数据处理检查

冻结执行1785项平面波、固定增益吸收、自由响应C2、共享频率协方差、参考阵元对比坐标与CA解析图坐标、同未知量/同原始噪声下Fraw-Fca半正定检查均通过。原始实虚噪声按同13个频率传播，未按9个独立频率对累计。

这仍不足以放行。附加独立物理链式导数检查192项，60 FAIL：单HLA非方位剖面真实秩4，中心差分报告5；部分剖面深度信息偏差超过10%。初始1785 PASS不能覆盖这项失败。

解析水平导数使用∂P/∂r_i=Σ w_m h_m(r_i)(−ik_m−1/(2r_i))，再乘精确阵元距离对目标位置导数和初态链式几何导数。瞬时方位切向导数亦由解析梯度构造。所有已有模态保留，未新增场景、求解或频率。源深列沿用原冻结有限差分，仅水平链与切向消除采用独立解析检查。

在单节点3快照下，消除每快照瞬时方位后，水平部分只能通过3个中心距离变化进入；另有1个共享深度，所以秩上限4。中心差分截断残差形成虚假的第五方向。双节点可以具有5阶秩，但不能用双节点通过覆盖单节点公式失败。

标准重建检查{v['reconstruction_checks']}项通过，说明档案可复现；它不验证有限差分弱方向的物理秩。完整检查{v['checks']}项：{v['PASS']} PASS、{v['FAIL']} FAIL。最终IMPLEMENTATION_INVALID。

数据处理统计依据：[Pollard](https://arxiv.org/abs/1107.3797)。这里的信息只是一阶局部Gaussian图坐标模型；未声称非线性CA有限噪声精确似然。
''',encoding='utf-8')
 report=f'''# E2九组频对非方位信息pilot：独立审计不放行

父提交{p.PARENT}；设计提交{d['design_SHA']}。

## 最终判定

IMPLEMENTATION_INVALID。原E2-G0 FAIL_UNCHANGED；E2-G0科学信息NOT_EVALUATED；本阶段EXPLORATORY_PILOT_ONLY；R4=0%。不投入全频带认证，不执行H3。先审查最小公式修复，另行授权后才能复跑。

冻结自动规则原输出NINE_PAIR_NONBEARING_MECHANISM_PROMISING，已保存FROZEN_RULE_DECISION.json。独立物理公式检查否决，不能将自动A输出当作接受结论。全部原计算表、输入Jacobian和未知量剖面矩阵未修改。

## 执行与证据完整性

固定24场景、源深200m、3快照、双8元14mHLA、5km基线、9对/13频率、80001/160001已有模态。新增KRAKEN=0、MC=0、录音=0。只有一次注册科学执行，1800秒上限内完成；科学表存档时间跨度约{d['execution_elapsed_s_from_saved_artifact_timestamps']:.3f}s，属于存档时间估计而非精确计时器读数。

末尾JSON序列化因NumPy int32计数失败。仅由已保存表按相同冻结谓词恢复原自动判定；无科学重跑或冻结代码更改。OUTPUT_SERIALIZATION_RECOVERY.json保留该输出层异常，不将其与物理公式失败混为一项。

独立重建18290项全部通过：独立MOD struct解析、全部中心/半步冷压力、参考阵元图坐标和QR未知量消除、原Gaussian复场源谱投影、全部局部有效信息与固定增益候选。随后解析水平链式和秩审计额外192项、60失败；总18422 PASS/60 FAIL。

## 致命公式问题

单HLA剖面瞬时方位后的精确水平Jacobian最多3个独立距离方向，加共享深度最多4个状态方向。全部48个场景—网格检查真实秩4，冻结FD报告5。假第五方向由中心差分截断残差和过细谱秩阈值产生；不能解释为速度/方位信息。部分单阵剖面深度信息与解析链最大差{summary['single_depth_chain_relative_max']:.6%}，超过10%诊断口径。

双HLA主分析解析与FD秩均为5；其范围信息最大相对差{summary['dual_range_chain_relative_max']:.6%}、深度信息最大差{summary['dual_depth_chain_relative_max']:.6%}。这些事实保留为修复线索，不能替代整轮公式准入。全局Jacobian Frobenius小差、数据处理通过和网格收敛不能单独认证弱/零状态方向。

## 表示、未知量与信息来源

P0_SINGLE/P0_DUAL为同资源逐节点源谱未知的原场空间上限；P1主HLA、P2各节点归一化非相干联合。C0响应已知；C1逐节点阵元复增益跨13频率、3时刻固定；C2每阵元/频率/快照自由响应吸收控制。源谱采用旧保守逐节点/频率/快照自由系数，不使用跨节点相干益处。

使用L0=diag(Qs⊗I,Qs⊗I)，LCA=diag(Qs⊗B+,Qs⊗B−)，同原始协方差C=L diag(v)Lᵀ、Cholesky白化。固定增益联合消除；CA固定增益相位抵消，幅度在整个记录上投影。全处理为冻结局部Gaussian图坐标模型，不是精确CA有限噪声似然。1785项原始控制均通过，但秩失败使本pilot未获准。

未融合声学与方位FIM。BEARING_REFERENCE.csv保留原121时刻方位模型，仅作不同资源条件下的参照。额外方位切向投影是机制诊断，不是增加观测或改变物理未知量。

以下全部为保留的、尚未获准的FD诊断，见AUDITED_DIAGNOSTIC_SUMMARY.json：

- C1非方位P1/P0_SINGLE信息trace比最小/中位/最大：{ratios['P1_over_P0_SINGLE']}。
- C1非方位P2/P0_DUAL信息trace比：{ratios['P2_over_P0_DUAL']}。
- P2 C1/C0 trace保留比：{losses['information_trace']}；未知深度范围信息保留比：{losses['range_information_depth_profiled']}；深度信息保留比：{losses['depth_information_horizontal_profiled']}。
- C2目标信息0，相关不确定度INF，无伪逆零方差。
- P2 C1相对1%最小未知深度距离信息{strength['RELATIVE']['range_information_min']} m⁻²、剖面水平状态深度信息{strength['RELATIVE']['depth_information_min']} m⁻²；不能转成定位精度承诺。
- P2归一化弱方向绝对分量中位，顺序r0/theta0/vx/vy/z：{strength['RELATIVE']['weak_direction_median_abs']}。P1假第五方向禁止解释。
- P2主C1剖面方位24场景内，FD/半步Jacobian最坏{numerical['fd_relative']:.6%}、两网格Jacobian最坏{numerical['mesh_relative']:.6%}、范围信息差{numerical['range_information_relative']:.6%}、深度信息差{numerical['depth_information_relative']:.6%}。深度已接近10%上限。
- 全部1152个稳定性记录中43项标LOCAL_INFORMATION_NUMERICALLY_UNRESOLVED，全部保留。主双阵C1两个噪声模型各24场景通过原稳定性检查，不能外推为其他表示全部稳定。

## 候选与低能量

仅每场景旧8局部±250m/±10m及组合，768条网格/噪声记录。固定增益通过整个时间—频率的GLS对数图坐标联合拟合，不逐对调整。校正后的CA Frobenius只是可行增益调整，非全局最小值证书；图坐标尺度是联合一阶噪声度量，非统计显著性或全域分支排除。

最终网格相对噪声候选最小CA分离/1%尺度{csum['RELATIVE']['min_CA_over_noise1']}、/5%尺度{csum['RELATIVE']['min_CA_over_noise5']}，192项数值区间排除零；排序变化{csum['RELATIVE']['rank_changes']}项。固定底条件最小比{csum['ABSOLUTE_FLOOR']['min_CA_over_noise1']}与{csum['ABSOLUTE_FLOOR']['min_CA_over_noise5']}，192项区间排除零。这些有限候选诊断不弥补秩失败，也不构成正确分支保留证书。

绝对参考幅值1.0824612978073539e-4来自旧14976个阵元样本中位数，固定1%/5%比例。源级未知，不能称真实SNR。相对噪声可能过度有利，弱场未过滤；LOW_ENERGY_NOISE_SENSITIVITY.csv保留每场景原场范数、功率和两种噪声下信息。

逐场景固定底/相对噪声的距离信息比最小/中位/最大：{floor['range_information_depth_profiled']}；深度比：{floor['depth_information_horizontal_profiled']}。固定底下某些场景明显损失信息，不能用群体最小值变化声称固定底普遍更好。

## 最小修复方案（本轮不执行）

1. 水平状态与瞬时方位切向导数采用本轮独立核查的解析链式导数；源深导数保留精确输出的1m/0.5m检查。
2. 对P1非方位剖面强制验证物理秩上限4；弱方向使用结构零空间与有限差分误差界，避免用1e-10谱阈值把截断残差变成信息。
3. 重新审查全部P0/P1/P2、C0/C1/C2的弱方向与有效信息稳定性；保持9对、相同噪声和候选，不扩频或网格。
4. 修复须另经研究负责人授权，不在本轮自动重跑。

本结论是实现公式准入失败，不是否定E2/FDSL方法族，也不是物理不可辨识。完整全频带投入尚未由本pilot获得依据。提交推送后停止。
'''
 (O/'PILOT_REPORT.md').write_text(report,encoding='utf-8')
 (O/'GPT_SYNC.md').write_text('# 九频对pilot最终交接\n\n~~~json\n'+json.dumps(dict(d,execution_SHA='CONTAINING_RESULT_COMMIT',remote_SHA='VERIFY_AFTER_PUSH'),ensure_ascii=False,indent=2)+'\n~~~\n\n最终IMPLEMENTATION_INVALID，自动A判定已被独立物理公式检查否决。\n',encoding='utf-8')
 m=Path('results/R4_MASTER')
 with (m/'R4_PLAN.md').open('a',encoding='utf-8') as f:f.write('\nE2 nine-pair pilot final IMPLEMENTATION_INVALID: retained automatic A overridden by physical single-HLA rank/chain audit; 60 formula failures. Original E2-G0 FAIL_UNCHANGED/information NOT_EVALUATED; exploratory only; R4=0%. STOP for formula repair review; no full band/H3 execution.\n')
 with (m/'R4_EVIDENCE_LEDGER.csv').open('a',encoding='utf-8') as f:f.write('E2-NINE-PILOT,../R4_E2_NINE_PAIR_INFORMATION_PILOT/PILOT_DECISION.json,IMPLEMENTATION_INVALID,Automatic A rejected by analytic chain rank and depth checks; all tables retained,PENDING_RESEARCH_LEAD_AUDIT;STOP;NO_CREDIT,0\n')
 p.dump(m/'R4_E2_NINE_PAIR_PILOT_PROGRESS.json',dict(stage='R4_E2_NINE_PAIR_NONBEARING_INFORMATION_PILOT',status='IMPLEMENTATION_INVALID',original_admission='FAIL_UNCHANGED',E2_G0_information='NOT_EVALUATED',pilot='EXPLORATORY_PILOT_ONLY',R4_percent=0,next='STOP_FOR_FORMULA_REPAIR_REVIEW',automatic_next=False))
 manifest={}
 for q in sorted(O.rglob('*')):
  if not q.is_file() or q.name=='EXECUTION_ARTIFACT_MANIFEST.json':continue
  raw=q.suffix in ['.npz','.png'];manifest[str(q)]=dict(kind='RAW' if raw else 'LF',sha256=p.sha(q,raw))
 for q in ['r4_e2_nine_pair_pilot.py','r4_e2_nine_pair_pilot_audit.py','r4_e2_nine_pair_pilot_chain_audit.py','r4_e2_nine_pair_pilot_finalize.py','r4_e2_nine_pair_pilot_verdict.py','r4_e2_nine_pair_pilot_report.py']:manifest[q]=dict(kind='LF',sha256=p.sha(q))
 p.dump(O/'EXECUTION_ARTIFACT_MANIFEST.json',dict(artifacts=manifest,design_SHA=d['design_SHA'],R4_percent=0))
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
