"""Assemble source-coverage and stop reports from saved read-only results."""
from pathlib import Path
import csv,json,hashlib,collections
R=Path(__file__).resolve().parents[3];O=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text(encoding='utf-8-sig'))
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig',newline='')))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def j(n,x):
 with (O/n).open('x',encoding='utf-8') as f:json.dump(x,f,ensure_ascii=False,indent=2);f.write('\n')
def md(n,s):
 with (O/n).open('x',encoding='utf-8') as f:f.write(s)
audit=read(O/'INDEPENDENT_REPREFLIGHT_RESULT.json');build=read(O/'REPAIR_BUILD_CHECKS.json');diag=read(O/'H2_REPRESENTATION_DIAGNOSTIC_V2.json');inventory=rows(O/'SNAPSHOT_INVENTORY.csv');cov={x['path']:x for x in rows(O/'SOURCE_COVERAGE_ADJUDICATION.csv')};freeze=read(O/'REPAIR_CODE_FREEZE.json')
assert all(sha(R/x['path'])==x['sha256'] for x in freeze['files'])
coverage=[]
for x in inventory:
 c=cov.get(x['path']);coverage.append(dict(path=x['path'],sha256=x['sha256'],discovery=x['discovery'],semantic_status=c['decision'] if c else 'DISCOVERY_ONLY_NOT_COMPLETENESS_CERTIFIED',purpose=c['purpose'] if c else '',reason=c['reason'] if c else 'Inventory retained. Do not infer exclusion, inclusion or complete source certification from metadata scan alone. Semantic adjudication stopped at unresolved range-only truth mapping.'))
with (O/'FULL_DISCOVERY_COVERAGE_STATUS.csv').open('x',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(coverage[0]),lineterminator='\n');w.writeheader();w.writerows(coverage)
final=dict(stage='CZ-D2-HISTORY-REPAIR-01',status='CZ_D2_HISTORY_REPAIR_PARTIAL_SOURCE_MAPPING_UNRESOLVED',parent_STOP_SHA='3c7291bd6d4683796deccb352d1e8c105db63d87',original_design_SHA='dc91982e67921eeb7a6aac4e12c482a980d5ad89',historical_snapshot='a580499475cc4148f1195e6256c4813c6e225726',V2_rows=744,original_preserved_rows=216,added_AUX_calibration=512,added_T2_static_calibration=16,V2_admitted_for_panel_generation=False,source_completeness='UNRESOLVED',independent_identity_and_rule_check=audit['status'],independent_check_PASS=audit['PASS_count'],independent_check_FAIL=audit['FAIL_count'],build_exact_equality_failed_check='H2_eight_existing_truths; retained. Separate corrected diagnostic shows <=2ULP and all32 coordinate differences inside original exclusion widths.',supplemental_diagnostic_reader_error='Preserved and withdrawn via DIAGNOSTIC_READER_ERRATUM.json; did not affect V2 or independent re-preflight.',full_snapshot_files=len(inventory),external_private_files_checked=build['external_files'],source_discovery_is_not_complete_semantic_adjudication=True,new_panel=0,new_noise=0,new_public_observations=0,formal_optimizer=0,Fisher=0,MC=0,formal_main_calls=0,formal_cold_calls=0,planned_configurations=48000,configuration_status='ALL_NOT_EVALUATED',main_changed=False,historical_Gates_changed=False,R4_percent=0,CZ_D3='NOT_OPENED',STOP=True,next_action='Research-lead adjudication of radial-only historical truth applicability; finish remaining semantic source review under authorized rule before any future panel materialization.')
j('REPAIR_FINAL_DECISION.json',final)
j('REPAIR_STOP_RECORD.json',dict(reason=final['status'],all_new_science_prohibited=True,old_locks_untouched=True,new_read_only_locks=['DISCOVERY_STARTED.json','REPAIR_STARTED.json','INDEPENDENT_REPREFLIGHT_STARTED.json'],execution_of_existing_lock_reset=False,formal_execution_resume_requires_new_authorization=True))
md('SOURCE_DISCOVERY_AND_COVERAGE.md',f'''# 历史来源发现与覆盖裁定

本轮以CZ-D1冻结快照 `{final['historical_snapshot']}` 为唯一历史目录边界。扫描了14,582份文件，记录逐文件Git对象与SHA256；读取40,140个数组头、4,495份CSV表头、JSON来源引用及生成代码文本。解析错误为0。这是来源发现，不是14,582份文件全部完成了科学语义审查。

## 已核实的修订

| 来源 | V2处理 | 证据 |
|---|---|---|
| 原12来源216条 | 原列值逐条继承 | 原始CSV/JSON行、单位换算、原设计哈希 |
| AUX-D1独立阈值校准 | 增加512条四维真值 | `h_truth[0:512,:]`，原阈值记录、1024条评分 |
| T2-G0静态校准 | 增加16条二维位置 | `CALIBRATION_GEOMETRIES.csv`；噪声参考冻结及18份CAL声场档案 |

V2共744条来源记录，不能解释为744组独立物理场景。它是尚未获完整性准入的修订候选表。原216条、原Commit A、原STOP及全部原锁保持不变。

T2-G0校准位置曾实际用于参考功率和前向校准计算；后续接收性能没有准入，不会使这些已使用位置重新成为未见。只纳入r、theta，v和psi留空，不补造静态目标的运动参数。

## 私有归档与复用关系

本轮读取并核验4,818份外部私有或校准档案：CAP06为4,416份，SP1为384份独立创新档案，T2为18份CAL声场。全部SHA身份符合原索引。CAP06主航迹的存储状态与原16状态及联合镜像一致；压力控制中的零状态是占位值，不能作为零距离目标登记。其定向压力波形不是同一匀速四维航迹。CAP06的4个校准入射角没有距离真值，P1为无物理坐标的数学玩具，均不能补造四维状态。

SP1实际校准标签对应原16静态位置，测试标签对应H1的8个离网格位置。AUX-D1测试私有1536行对应原12状态重复128次。H2私有8状态与H1 CSV存在最多2 ULP的读取表示差，见专门诊断；原排除邻域已覆盖，不能将其当作8个新开发几何。

## 未解决的来源映射：必须STOP

`results/R3_RC3_method_recovery/MFP_UPPER_BOUND_summary.csv`保存真实生成并评价过的三组径深状态：47.6 km/187 m、50 km/200 m、52.3 km/213 m。`r3_0_method_recovery.py`的`run_mfp_upper_bound`确实用这些状态调用`pressure_spectrum(mode,freqs,r_m,z_s,z_r)`生成观测。该模型没有theta、v或psi输入。

因此它们不是搜索网格，也不能称为未用于实验的参考状态。但原冻结合同只定义完整四维邻域和(r,theta)静态二维投影，没有定义仅r,z的历史真值如何排除新四维状态。将theta任意置0会虚构身份；将整条距离环排除则会增加一维排除规则。本轮未授权变更规则，故将此来源列为UNRESOLVED，不擅自采用任何一种解释。

在该歧义处停止语义裁定，所有未完成的来源审查保留在`FULL_DISCOVERY_COVERAGE_STATUS.csv`，状态为`DISCOVERY_ONLY_NOT_COMPLETENESS_CERTIFIED`。这些行不是已排除来源，也不是已通过来源。JSON发现的1,734个唯一外部引用同样不能冒称全部重新读取；实际核验范围以4,818行私有档案清单为准。

## 搜索候选与真值分开

MFP的`MFP_UPPER_BOUND_rz.csv`是距离—深度评分面，生成真值位于summary；T2的CATALOGUE是副本搜索目录，校准位置另表登记。不得把它们整表纳入已见真值。数值参考探针、求导扰动及候选起点需要按用途识别，不能仅凭列数为4就称为四维实验真值。

## 审计范围

独立程序未导入修订程序或科研生成器，使用标准库NPY读取、CSV/JSON原行复算、独立圆周距离和网格表达式，完成786项身份与规则检查，0失败；包括744×744个历史状态间邻域判别的两实现交叉核查。该PASS只证明已核实修订的机械正确性，不证明完整来源覆盖。整体完整性仍为UNRESOLVED。
''')
md('REPAIR_CHANGE_AND_EXCEPTION_RECORD.md',f'''# 修订差异、失败保留与诊断勘误

唯一拟改变的正式实验输入是历史排除表：216→744条，新增AUX512和T2静态16。当前只是候选V2，`DESIGN_FREEZE.json`及其中registry路径没有改写。生成算法、seed、阈值、估计器、统计和预算逐字节保持。

修订构建共8项聚合检查，其中7项通过、1项严格相等检查失败：H2 NPZ与Python CSV解析的32个数不是逐位一致。最大绝对差为{diag['max_absolute_difference']:.17g}，最大2 ULP。原H2生成器使用pandas.read_csv；来源面板和NPZ的历史SHA全部一致，32个差值全部处于原排除邻域内。保留FAIL，不将其重写为逐位一致。`SOURCE_COVERAGE_ADJUDICATION.csv`中H2的“equals”应解释为来源行与既有排除覆盖对应，严格数值相等不成立；本勘误优先。

补充诊断第一版忽略NPY的`fortran_order=True`，误将列优先数组解释为行优先，产生错误的大差值。保留第一版程序和输出，`DIAGNOSTIC_READER_ERRATUM.json`明确撤回该数值解释。第二版先验证固定整数布局控制，再按头信息解码，得到上述2 ULP结果。修订主程序使用numpy.load正确处理布局，V2未受影响；独立预检、正式实验均没有重跑。

没有通过修改科学门槛恢复准入。当前STOP的实质来源是径向真值到冻结排除规则的映射未决，以及该处停止后的剩余来源语义审查未完成。
''')
md('GPT_SYNC_CZ_D2_REPAIR_01.md',f'''# CZ-D2-HISTORY-REPAIR-01 同步报告

**结论：受控修订已封存；来源完整性UNRESOLVED，正式仿真继续STOP。**

父STOP提交：`{final['parent_STOP_SHA']}`。原设计：`{final['original_design_SHA']}`。本报告的结果提交SHA以包含本文件的Git提交为准，最终回传核对远端分支，不改main。

- V2为744条：原216条完整保留，新增AUX-D1校准512条，以及本次新发现的T2-G0静态校准16条。原表、原设计与原STOP不变。
- 冻结快照14,582份文件完成元数据及来源发现；4,818份相关外部私有/校准档案SHA核验通过。不能把发现扫描称为全量语义完整性通过。
- 独立只读预检786项通过、0失败，覆盖已登记来源身份、V2继承、新增真值、镜像/邻域/粗网格和原合同未变。它不是正式科学冷审计。
- 来源完整性未闭合：早期MFP确实使用了3组径深真值，但没有定义方位和运动状态。当前四维/静态二维排除合同没有规定仅距离真值的处理方式。未补造theta，未新增一维距离排除，也未将该来源擅自忽略。其后待裁定来源仍明确留为未完成。
- 保留H2严格浮点相等检查的1项FAIL；独立补充诊断显示最多2 ULP，已被原邻域覆盖，历史SHA一致。补充诊断曾有列优先读取错误，已保留原输出并单列勘误，未影响V2或正式实验。

**需要负责人裁定的事项**：仅r,z的历史生成数据，是否属于本轮四维“已见状态”的适用来源；若纳入，需要明确独立的一维投影排除政策。任何新规则都需另行授权。该裁定完成后仍需补齐剩余来源语义审查，不能直接恢复48,000配置。

新目标面板、随机创新、公开观测、主优化、Fisher、Monte Carlo、正式冷审计均为0。48,000配置全部NOT_EVALUATED。CZ-D3未启动；历史Gate、main、R4=0%保持。

入口：`SOURCE_DISCOVERY_AND_COVERAGE.md`、`HISTORICAL_STATE_EXCLUSIONS_V2.csv`、`SOURCE_COMPLETENESS_GAPS.json`、`INDEPENDENT_REPREFLIGHT_RESULT.json`、`REPAIR_CHANGE_AND_EXCEPTION_RECORD.md`、`REPAIR_FINAL_DECISION.json`。提交并核对科研分支后STOP。
''')
print(json.dumps(dict(V2_rows=744,decision=final['status'],files=len(list(O.iterdir())))))
