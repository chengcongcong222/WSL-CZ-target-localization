"""Assemble the REPAIR-02 audit delivery; no research computations."""
from pathlib import Path
import json,hashlib,csv,collections,datetime
O=Path(__file__).resolve().parent;D=O.parent;R=D.parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,s):
 with (O/n).open('x',encoding='utf-8') as f:f.write(s.rstrip()+'\n')
def j(n,x):save(n,json.dumps(x,ensure_ascii=False,indent=2))
final=json.loads((O/'SOURCE_COMPLETENESS_FINAL.json').read_text());assert final['SOURCE_SEMANTIC_COMPLETENESS']=='PASS_FIXED_SNAPSHOT_AND_BOUND_ARCHIVES'
with (O/'HISTORICAL_STATE_EXCLUSIONS_V3.csv').open(encoding='utf-8') as f:counts=collections.Counter(x['kind'] for x in csv.DictReader(f))
save('MFP_APPLICABILITY_ADJUDICATION.md','''# MFP部分真值适用性补充裁决

依据本轮研究负责人的REPAIR-02授权保存，原Commit A不变。

`PARTIAL_RZ_GENERATING_TRUTH_RECORDED_NOT_4D_EXCLUSION`

原MFP三组实际生成真值为(47.6 km,187 m)、(50 km,200 m)、(52.3 km,213 m)。S0/S1/S2重复使用，原始表9行、不同物理状态3组。完整行号、来源SHA256、代码入口见`MFP_PARTIAL_HISTORY.csv`。

`pressure_spectrum(mode,freqs,r_m,z_s,z_r)`未定义水平theta、v、psi。本轮不为缺失维度填0，不新增一维距离排除，也不把MFP的声压作为CZ-D2公开观测。部分历史曾被使用的事实永久保留；不能宣传新面板每个距离值都首次出现。

完整四维真值仍采用原四维邻域及theta/psi联合镜像；静态r/theta采用原二维投影及方位镜像。候选网格、数值参考、噪声、无物理坐标代数控制与设备校准分开记账。判断依据是实际生成用途及可定义的状态身份，不是文件名含truth或记录数量。

已知径向孔径/端射测深参考另有具体处理：相对径向速度不能冒充全局目标速度，掠射角不能冒充水平方位。代码明确的+x端射位置投影受原theta=0粗网格保护；浅海5.01–10 km阵列角度参考在原45–60 km域外。两者均保留实际生成用途，不能将其一概称为候选网格。见来源补充表。

本裁决只解决来源适用性；不构成正式实验执行授权。
''')
save('SOURCE_SEMANTIC_CLOSEOUT_REPORT.md','''# CZ-D2历史来源语义收口报告

本轮将固定历史快照及已绑定私有档案中的可能生成真值来源完成分类。结论为`PASS_FIXED_SNAPSHOT_AND_BOUND_ARCHIVES`，不声称未提供的档案或仓库中每个单一坐标数值均已被证明未使用。

## 范围与证据链

固定历史快照为a580499475cc4148f1195e6256c4813c6e225726，REPAIR-02父提交为65f9cec61c05e43b3a628f91b078bfdcc4b5318b。沿用REPAIR-01的14,582份文件发现目录、CSV字段、JSON键、数组头和私有档案入口。此次新增结构化JSON复核找到1,899个带状态字段的对象；它们包括源真值、域定义、候选和重复标签，不能当作1,899个独立真值。

首版语义表分类3,791份候选文件，通过21个来源家族的生成代码、实际输出和复用关系绑定。独立审计另扫描全部499份历史Python源码，得到61个与随机状态或径向采样相关的语句。11份脚本未在首版逐文件表中显式映射，因此审计如实保留一项FAIL；随后逐项补充代码用途和范围依据，最终分类表3,802份，无未裁定来源。并不是重新运行3,802个科研实验。

11份历史脚本不能由当前Python3.10 AST解析，均保留在`HISTORICAL_PYTHON_PARSE_STATUS.csv`；所有源码仍纳入词法扫描和依赖审阅。历史文件没有被修复或导入执行。源码静态审阅不是一般程序行为的形式化证明，本结论限定于已封存实验和可追溯输出。

REPAIR-01中4,818份私有档案身份核验沿用已接受的记录；本轮没有再次哈希全部4,818份文件。此次直接核查了3,791份候选来源的实际SHA及补充11份源码身份，并重新读取AUX512校准状态、SP1私有标签和H2归档状态等必要数据。不得将继承证据写成第二次全量私有数据复算。

## 真实生成状态与候选区分

| 来源 | 裁定与依据 |
| --- | --- |
| 原216条、AUX校准512、T2校准16 | V2的744条全部保留；原来源/行号和空缺维度不变 |
| P2测向误差与航向边界试验 | 实际生成了观测，显式补7个完整状态；航向扫描并非全部是候选网格 |
| P3困难候选对 | ref生成观测，alt只被评分；ref状态由V3覆盖；省略坐标的汇总行按pair_id回连冻结表 |
| P4、RC23闭环、B1 | 固定名义真值或原随机试验复用；B1变化的r/v条件网格不当作新生成真值 |
| AUX静态Gate0/1/2 | 代码定义target=(R,0)；补11个50–60 km静态位置，不填v/psi |
| H1/H2/RX/SP1、E1/E2/H3、T3 | 通过原面板、场景绑定及镜像关系追溯复用；不因新噪声或新配置再当作未见状态 |
| MFP | 3组r,z实际生成真值独立登记，不触发4D或静态r/theta排除 |
| 径向孔径测深参考 | 保留真实生成用途，但不编造全局速度或航向；明确端射投影和浅海参考的原粗网格/域外覆盖 |
| CAP06 | MAIN_TRACK回连16状态和镜像；无距离的定向校准、非匀速应力记录的零占位及P1代数控制分别记账 |
| T2 N1/N2及精度哨兵 | 存档声场/模板的确定性数值参考，不是接收性能面板；校准16条另行保留 |
| 软件预检、Jacobian控制 | 用于公式及代码一致性，不是正式性能或拒判阈值校准人口；不能将随机残差或导数探针登记成目标 |

V3为762条：7个P2完整状态和11个AUX静态位置新增，均已受原粗网格保护。显式登记用于来源追溯，不改变阈值、状态域、原面板SHA算法或排除逻辑。注册行数不是独立物理状态数量，原跨来源重复项仍完整保留。

## 独立复核与异常保留

完整来源审计记录为671项PASS、1项FAIL。唯一FAIL是上述11份脚本的逐文件语义映射缺口，未发生数值身份或邻域覆盖失败。补充闭合14项检查全部PASS，形成最终来源完整性判定；原FAIL不覆盖，也不对外合并成“全部首次通过”。

另有两次只读审计实现中止，均保留代码、启动锁和原因：第一次未处理P3汇总表的省略字段；第二次文本修订误改了truth_checks初始化并重复创建输出目录。最终审计版本v3与补充闭合程序独立留档。此前一次相对路径启动因runtime脚本改变cwd未打开程序，也已记录。均未生成CZ-D2面板或执行科研。

H2严格浮点相等FAIL和首次列/行存储读取错误继续保留。正确读取Fortran布局后最大差7.105427357601002e-15、最多2 ULP，在原冻结邻域内。此处结论是排除覆盖不受影响，不是逐位相等通过。

## 不变边界与下一步

原Commit A、原STOP、REPAIR-01、CZ-D1和历史科学Gate未修改。B0/B3、Q1/Q2、P0/P1、种子、面板算法、估计器、统计口径、科学预算全部未改。本轮没有物化新面板、随机创新、公开观测或估计输出，没有运行优化、Fisher、MC或正式冷审计。

`FORMAL_EXECUTION_AUTHORIZED=false`。48,000配置全部NOT_EVALUATED。来源收口只允许交由负责人独立审阅；仍需单独授权才能恢复Commit B。CZ-D3未启动，R4=0%，提交后STOP。
''')
save('GPT_SYNC_CZ_D2_REPAIR_02.md','''# CZ-D2-HISTORY-REPAIR-02 同步报告

本轮完成部分真值适用性裁定和固定历史范围内的来源收口。正式精度试验未恢复。

| 项目 | 结果 |
| --- | --- |
| SOURCE_SEMANTIC_COMPLETENESS | PASS_FIXED_SNAPSHOT_AND_BOUND_ARCHIVES |
| REGISTRY_IDENTITY_AND_COVERAGE | PASS；有效V3共762条 |
| MFP_PARTIAL_STATE_POLICY | PARTIAL_RZ_GENERATING_TRUTH_RECORDED_NOT_4D_EXCLUSION |
| DESIGN_CONTRACT_UNCHANGED | true |
| FORMAL_EXECUTION_AUTHORIZED | false |

V3完整继承原216条和V2的744条，显式补入7条P2四维真值、11条AUX静态位置。这18条均原已受粗网格保护，没有改变排除阈值或新面板生成规则。

MFP的9条输出对应3组r,z真值，独立表永久登记，不虚构theta/v/psi，不建立仅距离排除。H2逐位相等FAIL及读取错误保留；正确差异最多2 ULP、在原邻域内，排除覆盖不受影响。

最终语义表覆盖3,802份候选来源。独立来源审计671项PASS、1项FAIL：独立扫描发现11份脚本缺显式映射。逐项审阅后14项补充闭合检查PASS，原FAIL完整保留。两次只读审计程序中止及修正记录也保留，没有将其包装为一次无异常通过。4,818份私有档案哈希采用REPAIR-01已接受证据，本轮未再次全量哈希。

结论限固定a580499历史快照与已绑定可访问档案，不证明每个单独参数值从未出现。没有新增精度结果或新的物理结论。

审阅入口：`SOURCE_SEMANTIC_CLOSEOUT_REPORT.md`、`SOURCE_COMPLETENESS_FINAL.json`、`SOURCE_SEMANTIC_ADJUDICATION_FINAL.csv`、`V2_V3_ADDITIONS.csv`、`MFP_PARTIAL_HISTORY.csv`、`AUDIT_FINAL/INDEPENDENT_SOURCE_COMPLETENESS_AUDIT.json`、`INDEPENDENT_DISCOVERY_GAP_CLOSURE.json`。

父提交65f9cec61c05e43b3a628f91b078bfdcc4b5318b；本轮在research/cz-d2提交，提交SHA以Git封印为准。main保持a9a80191ee4d0f0cecd9db296198606b63a68080。原设计和历史STOP未修改。

新面板/创新/观测/优化/Fisher/MC/正式冷审计均为0；48,000配置全部NOT_EVALUATED，R4=0%。提交推送核对后STOP，等待负责人独立审计与单独恢复授权，不启动CZ-D3。
''')
j('REPAIR_02_INPUT_BINDING.json',dict(parent='65f9cec61c05e43b3a628f91b078bfdcc4b5318b',historical_snapshot='a580499475cc4148f1195e6256c4813c6e225726',original_A='dc91982e67921eeb7a6aac4e12c482a980d5ad89',effective_registry=dict(path='HISTORICAL_STATE_EXCLUSIONS_V3.csv',rows=762,sha256=sha(O/'HISTORICAL_STATE_EXCLUSIONS_V3.csv'),kinds=dict(counts)),supersedes_only_registry_input_after_future_authorization=True,design_manifest_sha256=sha(D/'DESIGN_FILE_MANIFEST.json'),old_registry_v2_sha256=sha(D/'HISTORY_REPAIR_01/HISTORICAL_STATE_EXCLUSIONS_V2.csv'),source_table_sha256=sha(O/'SOURCE_SEMANTIC_ADJUDICATION_FINAL.csv'),original_exclusion_thresholds_unchanged=True,FORMAL_EXECUTION_AUTHORIZED=False))
j('REPAIR_02_STOP_RECORD.json',dict(status='STOP_AWAIT_RESEARCH_LEAD_AUDIT',SOURCE_SEMANTIC_COMPLETENESS=final['SOURCE_SEMANTIC_COMPLETENESS'],FORMAL_EXECUTION_AUTHORIZED=False,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),formal_configurations=48000,evaluated=0,NOT_EVALUATED=48000,new_panel=0,new_innovations=0,new_observations=0,optimizer_calls=0,Fisher_calls=0,MC_calls=0,formal_cold_audit_calls=0,science_threshold_changes=0,CZ_D3='NOT_STARTED',R4='0%',previous_STOP_records_preserved=True))
print('Reports written. Kinds:',dict(counts))
