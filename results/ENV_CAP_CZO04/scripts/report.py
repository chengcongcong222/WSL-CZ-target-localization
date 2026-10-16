"""Reporting only: reads sealed CSV/NPZ; no model calls, no science rerun."""
from pathlib import Path
import csv,json,hashlib,subprocess,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path(__file__).resolve().parents[1];R=O.parents[1];W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP_CZO04')
def rows(name):
 with (O/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def rd(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def js(p,x):Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def md(p,x):Path(p).write_text(x,encoding='utf-8')
def table(p,x):
 with Path(p).open('w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(x[0]));w.writeheader();w.writerows(x)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
S=rows('FINITE_CANDIDATE_SUPPORT.csv');L=rows('LOCAL_EFFECTIVE_INFORMATION.csv');E=rows('SOURCE_DRIFT_INFORMATION_LOSS.csv');B=rows('BEARING_FINITE_GRID_BASELINE.csv');A=rd(O/'CZO04_INDEPENDENT_COLD_AUDIT.json');F=rd(O/'DESIGN_FREEZE.json')
idx={(q['scene'],q['truth_profile'],q['rho'],q['condition']):q for q in L};summary=[];stable=[]
for rho in ['0','0.9']:
 for sigma in ['0.1','0.25','0.5']:
  for c in ['C0','C1','C2']:
   aa=[q for q in S if q['rho']==rho and q['sigma_db']==sigma and q['condition']==c];nn=[q for q in aa if int(q['I2_retained'])>0];sh=[q for q in nn if float(q['range_width_km'])<float(q['bearing_range_width_km'])];good=[q for q in sh if int(q['disagreement'])==0 and idx[(q['scene'],q['truth_profile'],rho,c)]['numerically_resolved']=='True'];stable.extend([dict(q) for q in good])
   summary.append(dict(rho=rho,sigma_db=sigma,condition=c,cases=len(aa),nonempty=len(nn),empty=len(aa)-len(nn),range_shrink=len(sh),identical_I1_I2_candidate_masks=sum(int(q['disagreement'])==0 for q in aa),local_resolved_and_mask_identical_shrink=len(good),nonempty_width_median_km=float(np.median([float(q['range_width_km']) for q in nn])) if nn else '',continuous_domain='UNRESOLVED'))
table(O/'CONDITION_AND_NOISE_SUMMARY.csv',summary)
prim={a['condition']:a for a in summary if a['rho']=='0' and a['sigma_db']=='0.25'}
primary_stable=[a for a in stable if a['rho']=='0' and a['sigma_db']=='0.25' and a['condition']=='C1']
table(O/'NUMERICALLY_RESOLVED_FINITE_SHRINK_WITNESSES.csv',stable)
both=[a for a in E if a['rho']=='0' and a['both_numerically_resolved']=='True' and a['retained_fraction']!=''];fraction=[float(a['retained_fraction']) for a in both]
empty=[a['scene'] for a in B if a['bearing_retained']=='0'];nonregular=sorted({a['scene'] for a in L if a['bearing_regular']=='False'})
# Post-seal finite-grid evaluation controls only: no new forward values or refitting.
grid=np.round(np.array([45.,-5.,1.,-15.])+np.indices((16,21,11,31)).reshape(4,-1).T*np.array([1.,.5,.2,1.]),10);retention=[]
for b in B:
 if b['truth_on_grid']!='True':continue
 z=np.load(W/(b['scene']+'_SCORES.npz'));ii=z['bearing_ids'];loc=np.flatnonzero(np.all(np.isclose(grid[ii],z['state'],rtol=0,atol=1e-10),axis=1));assert len(loc)==1
 for p in range(9):
  for c in ['C0','C1','C2']:
   retention.append(dict(scene=b['scene'],truth_profile=p,condition=c,true_grid_candidate_present=True,residual_squared_db=float(z[f'I2_rho0_{c}'][p,loc[0]]),noiseless_regression_PASS=bool(z[f'I2_rho0_{c}'][p,loc[0]]<1e-8),scope='NOISELESS_MATCHED_GRID_CONTROL_NOT_STATISTICAL_COVERAGE'))
table(O/'NOISELESS_GRID_POINT_RETENTION_CONTROL.csv',retention)
assert all(a['noiseless_regression_PASS'] for a in retention)
mainusage=rd(O/'ANALYSIS_PROCESS_USAGE.json');coldusage=rd(O/'COLD_PROCESS_USAGE.json')
decision=dict(stage='ENV_CAP_CZO04',decision='LIMITED_CONDITIONAL_LEVEL_INFORMATION_GENERAL_RANGE_ANCHOR_NOT_ESTABLISHED',parent_SHA=F['parent_SHA'],design_SHA=subprocess.check_output(['git','-c','safe.directory='+str(R),'rev-parse','HEAD'],cwd=R).decode().strip(),cold_audit_PASS=A['PASS'],labels=dict(CZ_LEVEL_INFORMATION_PRESENT_CONDITIONALLY=dict(supported=len(primary_stable)>0,scope='one on-grid historical geometry G02 and its mirror; JA-generated source-depth250m; all nine estimator labels profiled; only two archived incoherent representations and finite distance labels'),SOURCE_DRIFT_ERODES_INFORMATION=dict(supported=True,scope='same feature, nested nuisance model; 29 local numerically resolved paired upper-bound entries, not independent experiments',median_C1_to_C0_information_increment_ratio=float(np.median(fraction)),resolved_pairs=len(both)),CURRENT_CZ_LEVEL_ROUTE_NOT_SUPPORTED=dict(global_physical_negative_claim=False,current_general_reliable_anchor_established=False,scope='limited conditional witness preserved; broad performance and continuous coverage not established'),OBSERVABLE_EXTRACTION_NOT_ESTABLISHED=True),primary=dict(case_count=324,C0=prim['C0'],C1=prim['C1'],C2=prim['C2']),primary_stable_witnesses=[dict(scene=a['scene'],truth_profile=int(a['truth_profile']),range_labels_km=json.loads(a['range_components']),retained_states=int(a['I2_retained'])) for a in primary_stable],finite_grid_empty_bearing_scenes=empty,unsigned_bearing_nonregular_scenes=nonregular,continuous_depth_and_environment='NOT_EVALUATED',true_grid_regression=dict(cases=108,per_condition_true_point_retained=108,statistical_coverage='NOT_EVALUATED'),new_external_solvers=0,new_MC=0,new_receive_records=0,science_wall_s=mainusage['elapsed_s']+coldusage['elapsed_s'],science_peak_RSS_bytes=max(mainusage['tree_peak_RSS_bytes'],coldusage['tree_peak_RSS_bytes']),R4='0%',main_merged=False,status='STOP',next_research='No automatic continuation. Before any receiver-performance study, independently define/extract the single-frequency incoherent-energy observable and establish source+gain stability and environment support. Do not resume Snell, T2, coherent providers, or candidate densification.')
js(O/'CZO04_FINAL_DECISION.json',decision)
# Tables for a direct scientific reading.
lines=['|源级/增益条件|均值案例|非空有限集|距离标签跨度缩小|I1/I2完全相同|局部数值可解释且集合一致的收紧|非空跨度中位 km|','|---|---:|---:|---:|---:|---:|---:|']
for c in ['C0','C1','C2']:
 a=prim[c];lines.append(f"|{c}|{a['cases']}|{a['nonempty']}|{a['range_shrink']}|{a['identical_I1_I2_candidate_masks']}|{a['local_resolved_and_mask_identical_shrink']}|{a['nonempty_width_median_km']:.1f}|")
main_table='\n'.join(lines)
md(O/'GPT_SYNC_CZO04.md',f'''# CZO04：第一会聚区强度信息与资源条件

**结论：存在很窄的条件性信息线索，尚未建立一般性的可靠距离锚。** 本轮已完成一次主研究、一次独立冷审计，提交后STOP。R4=0%，main不合并。

父SHA：`{F['parent_SHA']}`。设计SHA：`{decision['design_SHA']}`。最终SHA为本报告所在Commit B，远端`research/env-cap-czo04`核对记录保存在D盘，避免文件自引用。

## 最重要的研究结果

200 Hz、实际接收深度200 m的完整121时刻能量级模板，在全部九个环境/源深标签共同剖面后比较完整历史114576点。36个固定水平几何/镜像×9个生成标签=324个**观测均值案例**，不是324条独立航迹或随机试验；没有生成接收录音、传播场或MC。

主比较σ=.25 dB、独立特征误差假设下：

{main_table}

声级变化不普遍收紧距离。C1下52个案例的有限距离跨度缩小，只有**G02及其镜像、JA源深250 m生成条件**同时满足本轮局部数值诊断、两存档表示的候选一致性与收紧。对应各保留11个水平点，距离标签仅47 km；这是1 km步长有限目录的结果，**绝不是连续零宽度、零误差或工程测距精度**。估计侧仍枚举全部九标签，不知道生成标签。

C1局部主条目29/324满足数值解释条件。29个C0/C1均可解释的配对条目中，C1保留的有效距离信息增量占C0的中位比为{np.median(fraction):.3g}；其余不稳定结果保留，不能作强物理结论。线性漂移会同时消去大部分缓慢传播变化；这项代价不能以稳定源的高信息掩盖。

**C2每历元自由源级下，全部条件声级增量严格归零、候选严格等于方位基线。** 常数源/接收增益可消，但任意时间变化不可能由单频强度自行区分。

## 边界与资源

- 四个离网格场景`{', '.join(empty)}`的方位有限网格即为空，全部标记未决；没有用真值加点、选侧或补回分支。离网格空集不解释为物理无解。108个在网格内的生成标签控制每种条件均保留真实网格点，只是匹配均值回归，不是噪声覆盖率。
- 三个气候标签仍是一族WOA资料；深度仅150/200/250 m有限标签。局部Fisher为固定标签上限，连续深度/环境未知后的信息未评价。阵轴处无符号方位非正则另列。
- 只使用I类SHD非相干能量。**静态非相干平均能量不等于真实200 Hz线谱功率**；实际路径交叉项能否被时窗/统计平均消去没有接收数据支持。点接收模板也不是8元HLA波束产品。动态映射、AGC、波束增益、源级曲率、噪声协方差和源线稳定性尚需观测证据。
- 资源表给出最近≥1 km错误有限候选的白化间隙、漂移拟合量和噪声筛选尺度，不是最低硬件要求或真实SNR。现有1200 s登记轨迹径向行程约29.6–988.2 m，并非都能穿越充分的CZ形状。

冷审计**{A['checks']:,}项通过、0失败**；从原始SHD、显式几何与独立最小二乘重建全部病例，分数最大差{A['max_raw_cost_difference']:.3g} dB²，有效信息最大相对差{A['max_effective_info_relative_difference']:.3g}。科学累计{decision['science_wall_s']:.2f} s，采样峰值RSS约{decision['science_peak_RSS_bytes']/1024**2:.1f} MiB；新求解/MC/录音均0，所有新增数据在D盘。

建议：保留有限条件机制证据，**不据此直接投入大规模接收精度试验，不追加同类参数扫描或数值修补**。若另行研究，先取得源级+增益稳定性和单频能量观测定义的证据，并确定能独立约束这些干扰的观测。NUM03及T2/FUNC2保持收口，Yang/MMAR、相干提供器不自动启动。

完整证据：[科研决策](CZO04_RESEARCH_DECISION.md)、[误差条件汇总](CONDITION_AND_NOISE_SUMMARY.csv)、[逐案例候选](FINITE_CANDIDATE_SUPPORT.csv)、[局部信息](LOCAL_EFFECTIVE_INFORMATION.csv)、[资源尺度](RESOURCE_SEPARATION_SCALES.csv)、[冷审计](CZO04_INDEPENDENT_COLD_AUDIT.json)。
''')
md(O/'CZO04_RESEARCH_DECISION.md',f'''# CZO04：科研去留与申请证据

## 1. 真正解决的问题

本轮将未知源级处理与传播信息损失放入同一可复核模型：常数未知量可消去；加入线性漂移会丢失近似线性的传播变化；任意逐时源级完全吸收单频传播曲线。两个含义要分开：模型中存在弱/强传播信息，不等于接收系统已经取得该信息；信息为正，也不等于完整候选域已经收紧。

零控制及数据处理不等式通过，独立冷重建{A['checks']}项0失败。改变的是当前新研究证据，旧R2可观测链缺口、CZ envelope覆盖预算失败和T2声场精度失败均不追溯改判。

## 2. 完整有限目录得到什么

{main_table}

该表是σ=.25 dB、ρ=0下的同一组确定性均值，不能称保留率/成功率。主C1有52次有限范围缩小，其中50次未同时获得局部解释与两表示一致，保留为数值/表示限制。两次更完整证据来自同一个历史几何及镜像、同一个生成环境/源深标签，不算两个独立海域。主条件的非空距离跨度中位15 km与原域一样宽，显示普通单频强度模板并非普遍距离锚。

保持ρ=.9相关误差时，C0/C1范围缩小的均值案例分别只有28/18，说明独立121个样本的假设会显著影响可分离尺度。此处所有σ/.9相关参数为设计假设，没有由UUV数据估计。

源深和环境标签在候选端全部枚举，无真值锁定；未知连续深度、SSP、底质、接收阵深和波束变化均未包含，不得升级为一般可辨识性或实际覆盖声明。原±25 m平均和9点距离平滑没有用于本轮主表示。

## 3. 漂移与数值边界

全部C0/C1配对局部上限中，C1/C0新增有效距离信息比值中位约0.391；但大部分条目的Jacobian/有效信息跨网格尚不稳定，不能把这个总体中位当可信物理标定。两者都满足登记数值诊断的29项中，该比值中位{np.median(fraction):.8g}，范围{min(fraction):.3g}–{max(fraction):.3g}。这是有限案例内部的条件上限诊断，不是独立随机样本估计。

局部模型在125 m档案节点之间用分段线性dB插值，导数只是该表示的斜率。跨网格和两步长检查有解释力，但不是实际传播导数精度认证。C1模板网格差RMS中位约0.00552 dB，而Jacobian相对差中位约17.3%；模板级差小不自动保证导数稳定。跨网格最坏模板投影RMS约0.312 dB，最坏Jacobian相对差约40.5，完整失败保留。

无符号夹角在阵轴处非正则，G01两镜像的相关条目不准入局部Fisher解释。H01–H12中的4个离网格场景有限目录方位为空；它们的连续兼容域仍未求解。网格内108个生成标签控制均保留真点，只验证匹配模型的均值回归。

## 4. 为什么不能直接称真实接收强度信息

原数据来自静态I类射线/束非相干能量。它可以作为条件性期望能量模板，不能充当未知源CW在实际接收时窗的功率真值；真实多径交叉项、时间/频率平均条件、线谱关联和噪声偏差未认证。没有相干波数/模态，也没有原始波形，不能量化实际处理误差。

当前目标运动通过静态能量曲线沿距离采样处理，是准静态假设。源深仅3个离散标签，环境也仅同族3个气候剖面；这轮不新增深度精度，不替代Yang/MMAR。总能量中折射/边界族贡献未分离，几何证据仍有效，能量提取需另论。

## 5. 资源与研究去留

数学必要条件：单频每时刻源级/接收增益必须受独立约束；完全自由时任何传播曲线皆可解释。常数/线性假设是明确先验，不是估计器的鲁棒性贡献。

假设性上限：RESOURCE_SEPARATION_SCALES.csv给出登记标签与有限候选的最小间隙、半间隙扰动尺度和噪声筛选尺度；只有真点在目录、模型正确且扰动满足相应范数条件时，半间隙才有有限目录最近模板区分意义。对离网格情况只作描述。没有连续真点、输入误差界和可提取性，不能反推设备最低σ、真实SNR或所需最短观测时间。

已有条件配置：G02及镜像在一项JA/250 m生成条件留下可继续讨论的距离标签区分。其余场景和相关噪声压力结果一起保留，不能以这一个条件替代全体。部署层面不新增阵元和辅助节点，但需要单频功率归一化元数据、源级/增益稳定性、源深与环境范围资料、同步导航及可用径向行程；这些都不能无成本默认获得。

正式建议：不把本轮作为一般距离锚或工程精度Demo，不启动新的接收性能MC或传播求解。若未来另行授权，应先检验可跟踪线谱的实际功率统计/增益稳定性，或者引入能独立约束源级与传播的观测；不继续提高Snell精度、扩大现有目录或开发旧T2链。

申请可表述为：第一会聚区非相干强度模板的有限目录距离区分在局部条件成立，源/增益时间自由度显著控制信息保留；实际可提取性和连续候选覆盖是明确的后续研究问题。不能表述为原HLA已完成距离定位、真实鲁棒性或五维精度。

完成后STOP，全部历史Gate和R4=0%保持。
''')
md(O/'ZERO_CONTROLS_AND_INFORMATION_DERIVATION.md',r'''# 单频未知量剖面与零控制

令y=P(h,z,e)+Xβ+η，Σ=σ²R，W=chol(R)^(-1)。常数源及接收增益的和共用同一列，不重复引入可分离参数。C0 X=1；C1 X=[1,(t-600)/600]；C2 X=I。唯一主表示始终是121维原生单深度dB时间曲线。

令Q为WX列空间的正交基，Π=I-QQᵀ。剖面残差为ΠW(y-P)，信息JᵀWᵀΠWJ/σ²。因为0≤Π≤I，剖面化不增加信息；白化与投影的次序不能交换。未用独立噪声公式处理相关中心化误差，也没有把中心化后的121个量作为独立观测重复累计。

C2时WX满维，Π=0，所以任何候选传播变化都能被源级/增益吸收。C2的候选条件和局部增量必须严格等于仅方位结果。平坦P(r)的距离导数为0；常数和线性序列分别被C0/C1消去。合成控制与独立最小二乘重建均验证这些关系。

局部距离有效信息用最小二乘剖面另外三水平方向：||j_r-Proj(col(J_other))j_r||²。秩和最弱方向由预登记尺度下的SVD报告。缺秩方向不赋零方差，不输出未经认证的CRLB/P95。无符号方位在0夹角处不正则，相关条目不作为Fisher正结论。

局部信息逐一固定离散传播标签，因此只构成条件上限。有限候选则对所有环境/深度标签做最小残差；两种对象不能混称连续深度/环境剖面。连续源深与环境误差信息NOT_EVALUATED。

候选每次从完整历史114576点出发，先使用同一方位筛选，再为全部相容点计算所有9个声级标签。任意声级约束均与方位约束取交集，因此方位已不相容点不可能在本规则下被声级补救；这不是优化器丢失分支。不同左右角符号都在网格中，未向估计器提供真实侧向。

方位/声级分别使用χ²(.975,121)与χ²(.975,120/119)尺度。若今后有严格登记Gaussian随机观测且真点/真标签在候选中，固定真点检验可用联合界解释；本轮只有均值分离，没有实施这种随机覆盖验证。多标签并集、有限目录及连续点间区域均不能被称已认证95%定位区间。

I1/I2表示差异会影响接受标签，结果逐项保留两者、并集/交集与分歧数；两表示相同只能说明已测表示间一致，不能排除共同的传播偏差。完整连续状态域仍UNRESOLVED。
''')
md(O/'OBSERVABLE_AND_EXTRACTION_GAPS.md','''# 能量模板到接收观测的缺口

|层次|本轮证据|未建立的条件|
|---|---|---|
|首折返几何|继承NUM03，不重算|真实海域、多环境适用范围|
|总能量|实际200 m接收行、200 Hz I类SHD|总能量各路径族贡献、实海况统计|
|动态时序|静态曲线沿登记目标/平台距离采样|实际多径相干/非平稳与接收窗映射|
|线谱功率|假设可跟踪的特征级模板|原始波形、非相干平均算子、噪声底、失锁、关联|
|接收增益|与源级合并为常数/线性/自由项|AGC和波束权重、阵形/阵深、通道稳定性|
|源信号|未知常数或线性变化的条件模型|真实UUV功率波动、方向性及曲率统计|
|环境/深度|九个标签并集|连续深度、SSP/底质/海深失配与独立环境资料|
|定位集合|完整历史有限点枚举，未决区域保留|连续外包络、真值覆盖和工程精度|

官方手册明确区分I非相干和C相干TL：[BELLHOP官方说明](https://oalib-acoustics.org/website_resources/AcousticsToolbox/manual/node61.html)。稳定单频源本身并不使多径自动非相干。本轮没有验证实际时间平均会消掉交叉项，不能把存档I类能量直接宣布为可从目标录音取得的声级特征。

若每时刻的源级和增益之和完全自由，新增单频强度不提供信息，这是数学零控制；需要受独立观测支持的稳定性限制，或不同类型的额外观测，才能解除混淆。多频/空间/时延或外部协同可作为另行研究候选，不在本轮执行，也没有在这里证明其资源下界。

CZO04不修改旧H1、RX1Z、SP1及两航向名义候选成果：历史名义距离50 km、方位0°、速度2 m/s、航向4°/5°双候选保留；该三频历史结果与本轮200 Hz新气候模板不可直接比较。
''')
# Input and new scientific data identities, both Git and D.
rawindex=[]
for a in F['execution_inputs']:rawindex.append(dict(kind=a.get('source','frozen_input'),path=a['path'],bytes=Path(a['path']).stat().st_size,sha256=a['sha256']))
for p in sorted(W.glob('*_SCORES.npz')):rawindex.append(dict(kind='CZO04_ARCHIVED_FINITE_SCORES',path=str(p),bytes=p.stat().st_size,sha256=sha(p)))
table(O/'RAW_DATA_INDEX.csv',rawindex)
js(O/'STOP_RECORD.json',dict(stage='ENV_CAP_CZO04',status='STOP',analysis_launches=1,cold_launches=1,forbidden_automatic_continuations=['receiver performance MC','new propagation solvers','NUM03 Snell repair','T2 FUNC2','Yang MMAR','coherent provider','finer or enlarged candidate grid'],history_and_main='UNCHANGED',R4='0%'))
# Standalone research figures, all from saved tables only.
plt.rcParams.update({'font.size':9,'savefig.dpi':160})
fig,ax=plt.subplots(figsize=(8,4));xx=np.arange(3)
for i,rho in enumerate(['0','0.9']):
 a=[x for x in summary if x['rho']==rho and x['sigma_db']=='0.25'];ax.bar(xx+(i-.5)*.35,[x['range_shrink'] for x in a],.35,label='rho='+rho)
ax.set_xticks(xx,['C0 constant','C1 constant+linear','C2 free']);ax.set_ylabel('Finite mean cases with smaller range-label span / 324');ax.set_title('Conditional shrinkage: sigma = 0.25 dB; not coverage');ax.legend();fig.tight_layout()
fig.savefig(O/'FIG_SOURCE_CONDITION_AND_CORRELATION.png');fig.savefig(O/'FIG_SOURCE_CONDITION_AND_CORRELATION.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4));v=np.array(fraction);ax.semilogy(np.arange(len(v))+1,np.maximum(v,1e-12),'o');ax.axhline(np.median(v),color='orange',label='median');ax.set_ylim(1e-10,2);ax.set_xlabel('29 numerically resolved paired entries; not independent samples');ax.set_ylabel('C1 / C0 incremental effective range information');ax.set_title('Linear source/gain drift absorbs propagation trends');ax.legend();fig.tight_layout();fig.savefig(O/'FIG_DRIFT_INFORMATION_LOSS.png');fig.savefig(O/'FIG_DRIFT_INFORMATION_LOSS.pdf');plt.close(fig)
fig,axs=plt.subplots(3,3,figsize=(11,8),sharex=True);cur=rows('SINGLE_DEPTH_INPUT_LEVELS.csv')
for p,ax in enumerate(axs.flat):
 a=[q for q in cur if int(q['profile'])==p and 45000<=float(q['range_m'])<=60000]
 ax.plot([float(q['range_m'])/1000 for q in a],[float(q['I2_level_db']) for q in a],label='I2',lw=1)
 ax.plot([float(q['range_m'])/1000 for q in a],[float(q['I1_level_db']) for q in a],label='I1',lw=.8,ls='--');ax.set_title(['JW','JA','ANNUAL'][p//3]+f' source { [150,200,250][p%3] } m');ax.set_xlabel('Range (km)');ax.set_ylabel('Energy proxy (dB)')
axs[0,0].legend();fig.suptitle('Native 200 m receiver, 200 Hz: no depth-band or range smoothing');fig.tight_layout();fig.savefig(O/'FIG_ALL_NINE_NATIVE_POINT_CURVES.png');fig.savefig(O/'FIG_ALL_NINE_NATIVE_POINT_CURVES.pdf');plt.close(fig)
mat=np.array([[float(idx[(b['scene'],str(p),'0','C1')]['relative_J_grid_difference']) for p in range(9)] for b in B]);fig,ax=plt.subplots(figsize=(8,9));im=ax.imshow(np.log10(np.maximum(mat,1e-6)),aspect='auto',vmin=-3,vmax=2,cmap='magma');ax.set_yticks(range(36),[q['scene'] for q in B],fontsize=7);ax.set_xticks(range(9),[f'{e}/{z}' for e in ['JW','JA','ANNUAL'] for z in [150,200,250]],rotation=45);fig.colorbar(im,ax=ax,label='log10 relative I1/I2 projected Jacobian difference');ax.set_title('C1 local representation sensitivity, rho=0');fig.tight_layout();fig.savefig(O/'FIG_LOCAL_JACOBIAN_NUMERICAL_LIMITS.png');fig.savefig(O/'FIG_LOCAL_JACOBIAN_NUMERICAL_LIMITS.pdf');plt.close(fig)
print(json.dumps(dict(primary=prim,paired_resolved=len(both),drift_ratio_median=float(np.median(fraction)),decision=decision['decision']),ensure_ascii=False))
