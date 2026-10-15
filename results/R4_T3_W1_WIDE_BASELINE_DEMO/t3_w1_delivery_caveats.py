"""Text-only caveats update, based on already sealed statistics. No scientific rerun."""
from pathlib import Path
import hashlib,json
OUT=Path(__file__).resolve().parent
p=OUT/'GPT_SYNC_T3_W1.md';before=p.read_bytes();s=before.decode('utf-8-sig')
needle='每组最大场景偏置触界率和物理边界触碰见FOUR_PARAMETER_DEMO_TABLE.csv。触界干扰估计不是已完成校准；偏置和运动状态仍可能耦合。'
new='''B3在Q1和Q2的全部6000次输出中，均至少一个偏置参数触界。没有额外Gaussian先验并不意味着没有先验约束：Uniform有限支持和四个物理参数域都是条件。偏置触界不是校准成功，不能把较好的MC速度误差解释为数据单独恢复了偏置。未知偏置完全无独立支持时，最坏局部1.96σ速度等效误差仍为Q1约23.32%、Q2约33.93%；它与有界且可能有偏的MC输出属于不同统计对象。

Q2/H12的距离P95在B0与B3均为约0.882724%，是持平而非变差；其余11个场景改善。方位最坏P95由Q1的约0.1308°变为0.1413°、Q2的约0.1974°变为0.2119°，不能写成所有参数都改善。逐记录距离、速度与航向也均存在变差，完整计数在上方和PAIRED_CHANGES保留。

各物理参数和偏置触界率见FOUR_PARAMETER_DEMO_TABLE及逐场景表。B0最大单场景物理触界率Q1为57.6%、Q2为75.4%，B3对应0.6%/13.6%。本轮有界处理器与1200s条件下的资源增量不等于无边界或所有估计器的性能差距。'''
assert s.count(needle)==1;s=s.replace(needle,new)
s=s.replace('B3横向15 km并前置3 km','前置方向固定为公开全局+x轴，不按每个真实目标方位调整。B3横向15 km并前置3 km')
s=s.replace('增大基线的几何信息收益本身不是新算法。','新种子结果没有被要求复现负责人探索中的8.88%或0.40%：探索数据及处理器并非本轮冻结输入。本轮Q1的B3速度P95约8.57%、距离P95约0.70%，是当前有界偏置拟合和导航等效似然下的独立结果，不能把差别当成同处理器恒等回归。\n\n增大基线的几何信息收益本身不是新算法。')
p.write_text(s,encoding='utf-8');after=p.read_bytes()
q=OUT/'RESEARCH_ROUTE_AND_LIMITS.md';r=q.read_text(encoding='utf-8');r=r.replace('偏置触界率必须与误差表并读，触界的干扰状态不视为已校准。','B3的两个质量档中，全部输出均至少一个偏置参数触界；约束支持对结果有实质作用。偏置触界率必须与误差表并读，触界的干扰状态不视为已校准。');q.write_text(r,encoding='utf-8')
(OUT/'DELIVERY_TEXT_CAVEATS_PROVENANCE.json').write_text(json.dumps({'scope':'text-only clarification from sealed statistics; no refits or new information calculations','GPT_SYNC_before_sha256':hashlib.sha256(before).hexdigest(),'GPT_SYNC_after_sha256':hashlib.sha256(after).hexdigest(),'correction':'Q2/H12 rangeP95 equal, not worse; allB3biasbound; knownpriorandphysicalboundarylimitations','main_code_changed':False,'main_outputs_changed':False,'new_scientific_executions':0},indent=2)+'\n',encoding='utf-8')
# Refresh delivery file index only.
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'files':[{'path':p.relative_to(OUT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='OUTPUT_MANIFEST.json'],'D_raw_data_indices':['PUBLIC_BLOCK_IDENTITIES.csv','PRIVATE_INNOVATION_IDENTITIES.csv'],'scope':'allGitdeliveryfilesexceptmanifest;largeobservationsandinnovationsremainDindexed','R4_percent':0,'STOP':True}
(OUT/'OUTPUT_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print('Text caveats recorded; scientific data/code unchanged')
