"""Read-only report/figures from AUX-D1 sealed tables, no new inference."""
import json,math,pathlib
from aux_d1_common import *
PARAM=['r','theta','v','psi'];NAMES=['距离','方位','速度','航向'];UNITS=['%','°','%','°'];SCALE=[100,1,100,1]
def num(x,s=1):
 if x is None or x=='':return 'NA'
 if isinstance(x,str) and x=='INF':return 'INF'
 x=float(x);return 'INF' if math.isinf(x) else format(x*s,'.6g')
def write(n,t):(OUT/n).write_text(t.strip()+'\n',encoding='utf-8',newline='\n')
def report():
 verify();d=read('DECISION.json');v=read('VALIDATION.json') if (OUT/'VALIDATION.json').exists() else {'PASS':False,'FAIL':['audit incomplete']};summary=d['summary'];seal=read('INFERENCE_SEAL.json');clock=read('SCIENCE_CLOCK.json');table=['| 参数 | U0无条件中位汇总 | U5无条件中位汇总 | U0最坏几何无条件P95 | U5最坏几何无条件P95 | 改善≥30%且输出条件成立的几何 |','|---|---:|---:|---:|---:|---:|']
 for j,p in enumerate(PARAM):
  a,b=summary;cnt=len(d['primary']['signal_geometries_by_parameter'].get(p,[])) if p!='theta' else '不触发主判'
  table.append('|'+NAMES[j]+' '+UNITS[j]+'|'+num(a[p+'_unconditional_mean_geometry_median'],SCALE[j])+'|'+num(b[p+'_unconditional_mean_geometry_median'],SCALE[j])+'|'+num(a[p+'_unconditional_worst_geometry_P95'],SCALE[j])+'|'+num(b[p+'_unconditional_worst_geometry_P95'],SCALE[j])+'|'+str(cnt)+'/12|' if p!='theta' else '|'+NAMES[j]+' '+UNITS[j]+'|'+num(a[p+'_unconditional_mean_geometry_median'],SCALE[j])+'|'+num(b[p+'_unconditional_mean_geometry_median'],SCALE[j])+'|'+num(a[p+'_unconditional_worst_geometry_P95'],SCALE[j])+'|'+num(b[p+'_unconditional_worst_geometry_P95'],SCALE[j])+'|'+str(cnt)+'|')
 text='''# AUX-D1四参数配对结果

ABSOLUTE_BEARING_SIGN_CONDITION_ASSUMED。只比较单／双节点方位基线；新增节点收益包括额外观测、基线和登记处理器，不是纯基线因果。固定偏置没有被四参数soft_l1估计器消除。本轮为新的假设观测，不是恢复旧档案或海试数据。有限终点不是连续置信域，没有径向、深度或接收链验收。

'''+ '\n'.join(table)+'\n\n'
 text+='统计口径：每几何128次、两側各64次；上表中位汇总是十二几何无条件nearest-rank中位的等权平均，P95是最坏几何的无条件nearest-rank P95。拒判／求解失败记INF，未执行记NOT_EVALUATED。角度是绝对圆周误差，r/v是相对误差。每个场景和两侧的中位/P90/P95/P99/max，以及成功输出条件统计，分别保存在PAIRED_RESULTS_BY_CASE.csv、PAIRED_RESULTS_BY_SIDE.csv。P95仍为探索性，不证明真实总体成功率。\n\n'
 text+='输出与尾部：\n\n|资源|接受／已执行／计划|任一参数大误差输出|存在更好接受终点而未选|\n|---|---:|---:|---:|\n'
 for s in summary:text+=f"|{s['resource']}|{s['accepted']}/{s['evaluated']}/{s['planned']}|{s['large_any_output_count']}|{s['better_terminal_not_selected']}|\n"
 text+='\n大误差参考为距离/速度>10%、方位>1°、航向>5°，只是评价标签，未参与选解。逐参数大误差数量和配对改善／变差、双方仅一方有输出等完整分母见SUMMARY_BY_RESOURCE.csv与PAIRED_OUTCOME_COUNTS.csv。不能将接受次数称为精度达标次数。\n\n'
 text+='分类：'+ '; '.join(d['classifications'])+'。符合主信号的几何为'+str(d['primary']['qualifying_geometry_indices'])+'；逐参数为'+str(d['primary']['signal_geometries_by_parameter'])+'。U5速度P95>10%的几何为'+str(d['U5_speed_P95_above10percent_geometries'])+'。中位收益不能覆盖拒判、大误差或尾部恶化，速度限制也不是全部算法的物理上限。\n\n'
 text+='共享起点子集仅从0及2–127号的同批终点重建，分类为'+d['shared_start_subset']['primary_classification']+'，逐参数信号为'+str(d['shared_start_subset']['signal_geometries_by_parameter'])+'。资源特有1号的选点改变、输出改变及误差保存在SHARED_START_SUBSET_DIAGNOSTIC.csv；没有重跑或重校准，未替换主判。贪心簇尺度只作去重；窄终点簇和单簇不能解释成置信区间或全部物理解的分支数。\n\n'
 th=read('CALIBRATION_THRESHOLDS.json')['thresholds'];text+=f"独立校准：512组非测试目标与新误差，U0/U5各512个真物理状态处Q，不调用优化器。第500顺序值为U0={th['U0']:.17g}、U5={th['U5']:.17g}；门槛在测试生成前冻结，是包含假设偏置／导航影响的经验筛选，不是每个几何97.5%覆盖或输出准确率保证。\n\n"
 text+='条件：1200 s、121历元、主平台2 m/s且600 s后15°转向；5 km固定全局侧向辅助位移。随机方位0.05°，共模±0.05°、半差±0.025°、布放±10°均独立Uniform且每记录固定；导航每节点/历元/轴独立Gaussian25 m。同步、目标关联、有符号测向是假设，未覆盖阵姿态、真实INS相关性、海流和目标机动。两侧是布放侧，不能称原HLA左右舷已解决。\n\n'
 text+='计算与复核：'+str(seal['complete_pairs'])+'/1536对，'+str(seal['complete_configurations'])+'/3072配置，'+str(seal['terminals'])+'/393216终点；主耗时'+num(seal['main_elapsed_s'])+' s，实际nfev='+str(seal['actual_nfev'])+'，峰值RSS='+str(seal['peak_RSS_bytes'])+' bytes。相同起点与nfev上限不代表CPU成本相同，U5每次有双倍角样本。冷复核PASS='+str(v['PASS'])+'，最大Q差='+num(v.get('maximum_Q_absolute_difference'))+'；近同分仅记排序敏感，不换原ID。冷复核没有运行优化器；其通过不认证全域覆盖或真实硬件。\n\n'
 text+='历史Anchor A的1.30%/0.085°/28.94%/2.56°仅旁栏：本轮增加有界域、多起点和新拒判校准，不能声称复现旧6000次或计算其性能提升。原单阵全部声学方法并非U0基线，旧R3/H1/H2/RX/SP1/CLOSE1状态均未改。\n\n'
 text+='研究决策：'+('值得把辅助节点水平增量作为条件预研节点，后续另行设计接收资源、符号与同步／导航系统条件验证，保留速度尾部限制。' if len(d['primary']['qualifying_geometry_indices'])>=8 and d['execution_status']=='COMPLETE' else '登记处理器下的增量或执行完整性不足，保留全部边界；下一工作需独立审查，不自动修补搜索或更换条件。')+'本轮新增传播0，R4=0%；提交推送核对后STOP，不自动AUX-D2、垂向、径向或真实接收生成。'
 write('FOUR_PARAMETER_PAIRED_RESULTS.md',text)
 sync='''# AUX-D1同步报告

Parent：'''+PARENT+'；设计SHA：'+d['design_SHA']+'；执行SHA：CONTAINING_COMMIT，最终回报核对remote/main。任务SHA256：'+TASKSHA+'。\n\n'+text.split('\n\n')[1]+'\n\n'+'\n'.join(table)+'\n\n'
 sync+='主试验'+str(seal['complete_pairs'])+'/1536组成对观测，'+str(seal['complete_configurations'])+'/3072资源配置，'+str(seal['terminals'])+'个终点。每个几何128次，双侧各64；几何为已见开发状态，数据为新假设观测。U0和U5主节点输入逐字节相同，不称“单HLA全部声学方法对比”。\n\n'
 sync+='输出U0 '+str(summary[0]['accepted'])+'/'+str(summary[0]['evaluated'])+'、U5 '+str(summary[1]['accepted'])+'/'+str(summary[1]['evaluated'])+'；任一项大误差输出分别'+str(summary[0]['large_any_output_count'])+'、'+str(summary[1]['large_any_output_count'])+'。大误差参考r/v10%、θ1°、ψ5°，接受不等于达标，拒判及失败计INF。成功输出条件分位与两侧结果另列，未执行不是0。\n\n'
 sync+='判定：'+ '; '.join(d['classifications'])+'。主信号几何'+str(len(d['primary']['qualifying_geometry_indices']))+'/12；距离、速度、航向信号数分别'+str([len(d['primary']['signal_geometries_by_parameter'][p]) for p in ['r','v','psi']])+'。速度P95仍超10%的几何数'+str(len(d['U5_speed_P95_above10percent_geometries']))+'。中位收益和尾部并列，不要求四项同时10%，也不推断物理上限。\n\n'
 sync+='独立512组校准只评价真物理状态Q，无优化。U0/U5第500/512门槛为'+num(th['U0'])+'/'+num(th['U5'])+'，在测试生成前锁定；它是完整假设误差模型的经验门槛，不是逐几何97.5%保证。真值不进入估计器初始化或选点。\n\n'
 sync+='公共127起点＋1资源起点；共享子集分类'+d['shared_start_subset']['primary_classification']+'。所有公开方位、报告轨迹、时间、原始哈希、起点和float64终点已保存，私有真值／创新独立归档。有限接受终点和去重簇不是连续置信区间；评价发现的更好终点不会替换主解。\n\n'
 sync+='条件保持1200 s、121样本、600 s后15°转向、5 km固定全局侧向位移；随机方位0.05°、共模±0.05°、半差±0.025°、布放±10°、导航25 m/轴；有符号方位、正确同步关联均是假设，soft_l1不消除固定偏差。新增节点收益包含额外观测、几何与处理器，不能称纯基线因果。\n\n'
 sync+='主运行耗时'+num(seal['main_elapsed_s'])+' s；冷复核PASS='+str(v['PASS'])+'，最大Q差'+num(v.get('maximum_Q_absolute_difference'))+'，无新增优化；原始选点不变。旧双节点数值仅历史旁栏，旧科学结论与PENDING字段不升级。'+('建议继续独立设计接收资源与系统条件验证。' if len(d['primary']['qualifying_geometry_indices'])>=8 else '有限增量和限制用于下一路线审查。')+'新传播、径向、垂向、声压／音频／复谱、FIM全部0。D盘工作，R4=0%，推送后STOP。'
 write('GPT_SYNC.md',sync)
 write('AUX_D1_CLIENT_UPDATE.md','''# 辅助测向节点的条件增量预研

本轮用完全相同的主节点观测，比较只用主节点方位与加入约5 km辅助测向节点。结果仅限登记的有符号方位、导航与固定偏置模型，尚不是实际接收提取或设备验收。

'''+ '\n'.join(table)+'\n\n上表为各几何中位误差等权汇总与最坏几何探索性P95，拒判计INF。辅助节点的观测数量、几何和处理器共同产生条件收益；不与原单HLA全部声学方法作优劣比较。\n\n'+ ('已有足够的条件水平增量信号，值得下一步审查实际测向、导航同步和目标关联要求，同时保留速度尾部限制。' if len(d['primary']['qualifying_geometry_indices'])>=8 else '当前结果保留条件与限制，下一步需独立路线审查。')+'本轮没有启用垂向通道、径向特征、声学参考或改变原观测时长，不预先承诺四参数或五参数同精度。')
 # Only three descriptive figures from saved tables; no additional model analysis.
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 case=csvread('PAIRED_RESULTS_BY_CASE.csv');pair=csvread('PAIRED_OUTCOME_COUNTS.csv');side=csvread('PAIRED_RESULTS_BY_SIDE.csv');figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
 fig,axs=plt.subplots(4,2,figsize=(14,14),constrained_layout=True)
 for j,p in enumerate(PARAM):
  for k,metric in enumerate(['median','P95']):
   ax=axs[j,k]
   for R,color in [('U0','tab:blue'),('U5','tab:orange')]:
    rr=sorted([r for r in case if r['variant']=='PRIMARY' and r['resource']==R],key=lambda r:int(r['geometry_index']));ys=[float(r[p+'_unconditional_'+metric])*SCALE[j] if r[p+'_unconditional_'+metric] else math.nan for r in rr];finite=[(g,y) for g,y in enumerate(ys) if math.isfinite(y)];ax.plot([x for x,y in finite],[y for x,y in finite],'o-',color=color,label=R);bad=[g for g,y in enumerate(ys) if not math.isfinite(y)];ax.text(.01,.95 if R=='U0' else .86,R+' INF/NE: '+str(bad),transform=ax.transAxes,va='top',color=color,fontsize=8)
   ax.set_title(p+' unconditional '+metric+' ('+['percent','deg','percent','deg'][j]+')');ax.set_xlabel('Geometry index');ax.grid(alpha=.3);ax.legend()
 fig.savefig(figdir/'FIG1_FOUR_PARAMETERS.png',dpi=140);fig.savefig(figdir/'FIG1_FOUR_PARAMETERS.svg');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(14,5),constrained_layout=True);r=next(x for x in pair if x['geometry']=='ALL');labels=['Both','U0 only','U5 only','Neither'];axs[0].bar(labels,[int(r[x]) for x in ['both_accepted','U0_only','U5_only','neither']]);axs[0].set_title('Paired output outcomes (evaluated pairs)');counts=[[s[p+'_large_output_count'] for p in PARAM]+[s['large_any_output_count']] for s in summary]
 for k,R in enumerate(['U0','U5']):axs[1].bar([j+(.0 if k==0 else .35) for j in range(5)],counts[k],width=.35,label=R)
 axs[1].set_xticks([j+.175 for j in range(5)],PARAM+['any']);axs[1].set_title('Accepted outputs exceeding evaluation references');axs[1].legend();fig.savefig(figdir/'FIG2_OUTPUT_AND_TAIL.png',dpi=140);fig.savefig(figdir/'FIG2_OUTPUT_AND_TAIL.svg');plt.close(fig)
 fig,axs=plt.subplots(1,2,figsize=(14,5),constrained_layout=True)
 for R in ['U0','U5']:
  rates=[sum(int(r['accepted']) for r in side if r['resource']==R and int(r['side'])==s)/768 for s in [-1,1]];axs[0].plot([-1,1],rates,'o-',label=R)
  ee=[int(r['cluster_count']) for r in csvread('FINITE_TERMINAL_ENVELOPES.csv') if r['resource']==R and r['cluster_count']!='NOT_EVALUATED'];axs[1].hist(ee,bins=range(0,max(ee,default=1)+2),alpha=.5,label=R)
 axs[0].set_title('Output rate by auxiliary deployment side');axs[0].set_xticks([-1,1]);axs[0].legend();axs[1].set_title('Finite terminal cluster count (not physical branch count)');axs[1].legend();fig.savefig(figdir/'FIG3_SIDE_AND_FINITE_CLUSTERS.png',dpi=140);fig.savefig(figdir/'FIG3_SIDE_AND_FINITE_CLUSTERS.svg');plt.close(fig)
 payload={'task_SHA256':TASKSHA,'parent_SHA':PARENT,'design_SHA':d['design_SHA'],'execution_SHA':'CONTAINING_COMMIT','files':[filemeta(p) for p in OUT.rglob('*') if p.is_file() and p.name!='OUTPUT_MANIFEST.json'],'new_scientific_run':False}
 dump('OUTPUT_MANIFEST.json',payload);print('REPORTS_WRITTEN',flush=True)
if __name__=='__main__':report()