"""Frozen report/label logic, evaluated only after independent cold audit."""
from common import *

def val(x):return None if x in ['NOT_EVALUATED','',None] else float(x)
def fmt(x,percent=False):
 if x is None:return 'NOT_EVALUATED'
 if math.isinf(x):return 'INF'
 return f'{x*(100 if percent else 1):.6g}'+('%' if percent else '')

def main():
 budget=js(O/'EXECUTION_BUDGET.json');audit=js(O/'COLD_AUDIT.json') if (O/'COLD_AUDIT.json').exists() else dict(PASS=False,scope='NOT_EVALUATED');completion=js(O/'MAIN_COMPLETION.json') if (O/'MAIN_COMPLETION.json').exists() else dict(executed=0,accepted=0,not_evaluated=48000)
 stats=csvread(O/'STATISTICS.csv') if (O/'STATISTICS.csv').exists() else [];agg=csvread(O/'AGGREGATES.csv') if (O/'AGGREGATES.csv').exists() else [];diag=csvread(O/'OUTPUT_AND_BOUNDARY_DIAGNOSTICS.csv') if (O/'OUTPUT_AND_BOUNDARY_DIAGNOSTICS.csv').exists() else []
 S={(r['config'],r['case_id'],r['side'],r['metric'],r['mode']):r for r in stats};A={(r['config'],r['metric']):r for r in agg};Dg={(r['config'],r['case_id'],r['side']):r for r in diag};integrity=bool(audit['PASS'] and not budget['failure']);complete=completion['executed']==48000
 labels=[]
 if complete and integrity:
  for Q in ['Q1','Q2']:
   for P in ['P0','P1']:
    b=Q+'_B0_'+P;c=Q+'_B3_'+P
    for m in METRICS:
     medpairs=[(val(S[b,f'D{i+1:02d}','BOTH',m,'unconditional']['P50']),val(S[c,f'D{i+1:02d}','BOTH',m,'unconditional']['P50'])) for i in range(12)]
     old=val(A[b,m]['worst_scene_P95']);new=val(A[c,m]['worst_scene_P95']);sideok=all(float(Dg[c,f'D{i+1:02d}',s]['output_rate'])>=float(Dg[b,f'D{i+1:02d}',s]['output_rate']) for i in range(12) for s in ['-1','1'])
     gain=all(y<=x for x,y in medpairs) and any(y<x for x,y in medpairs) and math.isfinite(new) and new<old and sideok
     tail=not sideok or new>old or any(not math.isfinite(val(S[c,f'D{i+1:02d}',s,m,'unconditional']['P95'])) for i in range(12) for s in ['BOTH','-1','1'])
     labels.append(dict(quality=Q,support=P,parameter=m,label='CONDITIONAL_HORIZONTAL_HOLDOUT_GAIN_SUPPORTED' if gain else 'UNIFORM_GAIN_CRITERION_NOT_MET',TAIL_OR_COVERAGE_NOT_ESTABLISHED=tail,median_nonworse_scenes=sum(y<=x for x,y in medpairs),median_better_scenes=sum(y<x for x,y in medpairs),worst_P95_B0=old,worst_P95_B3=new,both_side_output_nonworse=sideok))
  for Q in ['Q1','Q2']:
   for B in ['B0','B3']:
    for m in METRICS:
     reasons=[]
     for i in range(12):
      case=f'D{i+1:02d}';a=val(S[Q+'_'+B+'_P0',case,'BOTH',m,'unconditional']['P95']);b=val(S[Q+'_'+B+'_P1',case,'BOTH',m,'unconditional']['P95'])
      if math.isfinite(a) and (not math.isfinite(b) or (a>0 and b>=1.2*a)):reasons.append(case+':P95_INCREASE')
      if float(Dg[Q+'_'+B+'_P1',case,'BOTH']['output_rate'])<=float(Dg[Q+'_'+B+'_P0',case,'BOTH']['output_rate'])-.02:reasons.append(case+':OUTPUT_DROP')
     # Compare B3-B0 worst-scene P95 direction under each support, without INF subtraction.
     v=[val(A[Q+'_'+bb+'_'+pp,m]['worst_scene_P95']) for pp in ['P0','P1'] for bb in ['B0','B3']]
     if (v[1]<v[0] and v[3]>v[2]) or (v[1]>v[0] and v[3]<v[2]):reasons.append('B3_B0_WORST_P95_DIRECTION_FLIP')
     labels.append(dict(quality=Q,geometry=B,parameter=m,label='GAIN_PRIOR_SUPPORT_SENSITIVE' if reasons else 'REGISTERED_PRIOR_SENSITIVITY_FLAG_NOT_TRIGGERED',reasons=';'.join(reasons)))
 else:labels=[dict(label='IMPLEMENTATION_OR_NUMERIC_GATE_FAIL' if not integrity else 'TAIL_OR_COVERAGE_NOT_ESTABLISHED',executed=completion['executed'])]
 csvwrite(O/'RESEARCH_LABELS.csv',labels)
 decision=dict(stage='CZ-D2_COMMIT_B_RESUME_01',design_SHA=DESIGN,authorization_SHA=PARENT,registry_SHA256=V3HASH,integrity_pass=integrity,complete_panel=complete,completion=completion,scope='FROZEN_ESTIMATOR_BIAS_SUPPORT_SENSITIVITY;DUAL_NODE_SIGNED_BEARING_CONDITIONAL_FOUR_PARAMETER_DEMO',labels_file='RESEARCH_LABELS.csv',P1='ESTIMATOR_SUPPORT_DOUBLED;GENERATION_BIAS_REMAINS_P0;NOT_OUT_OF_SUPPORT_MISMATCH',new_panel='PREREGISTERED_NOT_SECRET_EXTERNAL_BLIND_PANEL',continuous_branch_coverage='NOT_ESTABLISHED',received_chain='NOT_VALIDATED',hardware='UNKNOWN',depth='NOT_EVALUATED',R4_percent=0,CZ_D3_authorized=False,STOP=True)
 dump(O/'CZ_D2_FINAL_DECISION.json',decision,True)
 lines=['# CZ-D2正式补强结果','',f'执行{completion["executed"]}/48000配置，接受{completion["accepted"]}；未执行{completion["not_evaluated"]}。独立冷审计PASS={audit["PASS"]}。', '', 'B0与B3均为双节点；P0为原偏置支持，P1仅扩大估计支持。以下为新12几何面板、冻结有符号方位/导航假设下的描述性结果。', '', '|配置|参数|池化无条件P50|12场景最坏无条件P95|','|---|---|---:|---:|']
 for conf in sorted({r['config'] for r in agg}):
  for m in ['range_rel','bearing_circular_deg','speed_rel','heading_circular_deg']:
   r=A[conf,m];pct=m in ['range_rel','speed_rel'];lines.append('|'+conf+'|'+m+'|'+fmt(val(r['pooled_P50']),pct)+'|'+fmt(val(r['worst_scene_P95']),pct)+'|')
 lines+=['','上述池化P50与“12个场景P50的算术均值”不同，后者及其nearest-rank中位数见AGGREGATES。逐场景、分侧、条件/无条件P50/P90/P95/P99和最大值见STATISTICS，完整误差与终点见NPZ归档。失败全部记INF；未执行项保留NOT_EVALUATED，不充当求解失败。','','|配置|接受/计划|偏置至少一维触界|P1超过原支持|','|---|---:|---:|---:|']
 for conf in sorted({r['config'] for r in diag}):
  x=Dg[conf,'POOLED','BOTH'];lines.append(f'|{conf}|{x["accepted"]}/{x["planned"]}|{x["bias_at_bound_count"]}|{x["outside_original_bias_count"]}|')
 lines+=['','P0/P1改变了偏置归一化尺度和数值求解轨迹，不能把全部差异解释为物理信息变化。原支持与当前支持的物理距离保存在每个终点；P1生成偏置仍位于P0范围。','','本轮不要求复现旧8.5659267%的最坏场景速度P95，不设5%或10%工程验收门。RESEARCH_LABELS仅按Commit A描述性规则判别；有利中位数不能抵消较差尾部、拒判或偏置支持依赖。','','信息：INFORMATION和INFORMATION_MATRICES分别保存96个几何的自由偏置数据Fisher及已知偏置上限，共192行。没有向Uniform支持添加Gaussian先验，也不将局部1.96σ解释成经验P95。局部秩/弱方向与单局部终点不能构成连续分支覆盖证书。','','资源：RESOURCE_GEOMETRY保留两节点全部轨迹的距离范围、域外点数和5/15.297km部署成本。45–60km仅为历史研究几何范围，不证明第一CZ实际可检测性。同步、关联、有符号测向、导航与系统误差分布仍是假设。','','身份与审计：原Commit A、两次来源STOP、REPAIR-01/02不变；671 PASS/1 FAIL与14项补充闭合仍分别保留；MFP政策及H2严格相等FAIL未改。新绑定明确读取762条V3；生产代码在新面板生成前封存。估计器实际文件访问记录在ESTIMATOR_FILE_ACCESS；独立冷审计没有导入主生成/拟合/信息/统计函数或重复优化。', '', f'预算：主/冷调用={budget["main_calls"]}/{budget["cold_calls"]}，重试0；总壁钟{budget["total_wall_s"]:.3f}s，采样进程树RSS峰值{budget["peak_tree_RSS_bytes_sampled"]}字节，D盘数据{budget["D_bytes"]}字节。0.2s轮询是采样监视，不宣称连续峰值保证。', '', '原始公开观测和隔离私有档案保存在D:/ProjectStorage/WSL-CZ/CZ_D2/RESUME_01_DATA；索引、全数组/文件SHA和完整终点/误差归档随Git交付。源代码、版本与执行锁均保留。','','本轮只交双节点四参数条件证据。单HLA、深度、五维联合性能没有升级；深度NOT_EVALUATED，R4=0%。完成推送后STOP，不自动启动CZ-D3。']
 (O/'CZ_D2_RESEARCH_REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 sync=['# CZ-D2 Commit B 同步报告','',f'设计SHA：`{DESIGN}`；恢复授权SHA：`{PARENT}`。','',f'正式配置：{completion["executed"]}/48000；接受：{completion["accepted"]}；未执行：{completion["not_evaluated"]}。',f'主分析/冷审计：{budget["main_calls"]}/{budget["cold_calls"]}次；重试0；冷审计PASS={audit["PASS"]}。','', '主结果（12场景最坏无条件速度P95）：','', '|条件|B0|B3|','|---|---:|---:|']
 for Q in ['Q1','Q2']:
  for P in ['P0','P1']:
   k=(Q+'_B0_'+P,'speed_rel');j=(Q+'_B3_'+P,'speed_rel');sync.append('|'+Q+'/'+P+'|'+(fmt(val(A[k]['worst_scene_P95']),True) if k in A else 'NOT_EVALUATED')+'|'+(fmt(val(A[j]['worst_scene_P95']),True) if j in A else 'NOT_EVALUATED')+'|')
 sync+=['','全部四参数、侧别、拒判、触界与配对恶化见CZ_D2_RESEARCH_REPORT.md及完整CSV。P1是冻结处理器的支持敏感性，不是无先验或真实偏置超界实验。新面板不是保密外部盲测。','', '原main、历史Gate和R4=0%保持；深度NOT_EVALUATED。最终SHA由本次Git提交封印，推送后核对remote research/cz-d2。STOP，等待独立审计，不启动CZ-D3。']
 (O/'GPT_SYNC_CZ_D2_B.md').write_text('\n'.join(sync)+'\n',encoding='utf-8')
 print('REPORT_WRITTEN',completion,flush=True)
if __name__=='__main__':main()
