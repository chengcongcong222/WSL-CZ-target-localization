"""AUX-D1 design packaging only, without noisy panel performance execution."""
import shutil,sys,platform,subprocess,json,pathlib,ast
import numpy as np,scipy
from scipy.stats import qmc
from aux_d1_common import *
from aux_d1_model import LOW,HIGH
PKG=pathlib.Path('D:/Downloads/AUX_D1_TASK_PACKAGE_20261012_v1')
def prepare():
 OUT.mkdir(exist_ok=True);LOCAL.mkdir(exist_ok=True)
 if (OUT/'PREPARED.json').exists():raise RuntimeError('Package already prepared')
 assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==PARENT
 package=json.loads((PKG/'AUX_D1_PACKAGE_MANIFEST_20261012_v1.json').read_text(encoding='utf-8-sig'))
 for name,meta in package['files'].items():
  assert sha(PKG/name)==meta['sha256'] and (PKG/name).stat().st_size==meta['bytes'];shutil.copyfile(PKG/name,OUT/name)
 shutil.copyfile(PKG/'AUX_D1_PACKAGE_MANIFEST_20261012_v1.json',OUT/'AUX_D1_PACKAGE_MANIFEST_20261012_v1.json')
 sources=[]
 for item in read('AUX_D1_SOURCE_ENTRYPOINTS_20261012_v1.json')['sources']:
  p=item['path'];raw=(ROOT/p).read_bytes();blob=subprocess.check_output(['git','show',PARENT+':'+p]);bid=subprocess.check_output(['git','rev-parse',PARENT+':'+p],text=True).strip()
  if item['git_blob_sha_when_available']:assert bid==item['git_blob_sha_when_available']
  sources.append({'source_id':item['source_id'],'path':p,'read_at_commit':PARENT,'git_blob_id':bid,'git_blob_byte_sha256':hashlib.sha256(blob).hexdigest(),'working_raw_byte_sha256':hashlib.sha256(raw).hexdigest(),'purpose':item['purpose'],'full_text_read':True})
 dump('SOURCE_MANIFEST.json',{'sources':sources,'task_SHA256':TASKSHA})
 panel=list(csv.DictReader((ROOT/'results/R4_A1_NEW_AUGMENTED_OFFGRID_DYNAMIC/A1_NEW_TRUTH_PANEL.csv').read_text(encoding='utf-8-sig').splitlines()));h=np.array([[float(r[k]) for k in ['r0_km','theta0_deg','v_mps','psi_deg']] for r in panel]);assert len(h)==12 and np.all((h>=LOW)&(h<=HIGH));dump('PANEL_DOMAIN_PREFLIGHT.json',{'all12_original_physical_states_in_domain':True,'cases':[r['case_id'] for r in panel],'truth_file_used_only_for_domain_check_and_generator_evaluation':True})
 points=qmc.Sobol(d=4,scramble=True,seed=2026101204).random_base2(7);np.save(OUT/'SOBOL_128_RAW.npy',points);np.save(OUT/'SHARED_FIXED_126.npy',LOW+(HIGH-LOW)*points[:126])
 plan=[{'configuration':2*(rep*12+g)+r,'paired_record_index':rep*12+g,'geometry_index':g,'replicate':rep,'resource':['U0','U5'][r],'starts':128,'side':-1 if rep%2==0 else 1} for rep in range(128) for g in range(12) for r in range(2)];csvwrite('CONFIGURATION_PLAN.csv',plan)
 (OUT/'.gitattributes').write_text('* -text\n',encoding='utf-8');dump('PREPARED.json',{'utc':utc(),'task_id':package['task_id'],'task_SHA256':TASKSHA,'source_count':len(sources),'no_calibration_or_test_observations_generated':True})
 docs={
 'SOURCE_AND_SCOPE.md':'''# AUX-D1来源与作用域

ABSOLUTE_BEARING_SIGN_CONDITION_ASSUMED。本轮比较U0单节点方位与U5双节点方位，同主节点输入逐字节相同。不是原HLA全部声学方法对比，不是控制观测数量后的纯基线因果。

源文件按9060ca4固定读取，Git blob与工作区原始字节分开记录。旧原始方位／报告轨迹缺失不重搜、不从真值声称恢复。新数据是新独立仿真，12状态为已见开发几何，512校准状态为一次独立Uniform抽样。原三档、R3/H1/H2/RX/SP1和原PENDING审计字段不改。

SciPy least_squares官方文档说明其为局部稳健最小二乘，固定损失并不消除整段偏差：https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html 。本机冻结SciPy1.15.3，显式选项见完整任务与代码，不临时升级。Sobol完整128点再取126点，不宣称截取后仍有完整平衡性质：https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.random_base2.html 。

估计器以报告MAIN_t0为参照，评价真值以真实MAIN_t0为参照，保持导航参照影响。未知深度不在此水平模型内，不称已知或已估计。同步、关联、阵姿态与被动有符号方位为条件，不是硬件实测。没有传播、径向、源频相位、垂向、FIM或真实接收计算。
''',
 'DRAFT_TO_EXECUTION_CHANGES.md':'''# 从CLOSE1草案到本轮执行

本轮新增授权只在新目录生效。历史草案PROPOSED_NOT_FROZEN继续保持，不回填旧文档。

新增新种子完整观测档案；512独立校准记录只评价真物理状态Q，第500顺序值分别校准U0/U5，在测试生成前锁定；四状态soft_l1不估计真实偏差；共同有界域与128起点明确为新处理器；127共同与1资源起点、现有终点共享子集、严格float64 Q/start_id选择、经验拒判和大误差标签全部事前登记。

细化实现约定：按历元优先／节点其次展开残差；stable terms以NumPy float64 sum聚合，Q直接算，不用solver_cost替代选点。边界触碰阈值为归一化u距边界≤1e-8，只标记不拒判。128起点映射后clip至[1e-12,1-1e-12]；重复不增补。报告与真值隔离由独立生成／推断／评价／审计进程及封存哈希执行。科学耗时累计为各科学进程实际墙钟耗时总和，编辑/Git不占科学预算；单进程运行时含其全部等待与I/O。

PAIRED_TRADEOFF_PRESENT的描述性触发记录总面板的逐参数大误差计数、任一项大误差或拒判是否增加，同时全量保留逐几何变化。共享子集若主分类不同才追加INITIALIZATION_DEPENDENCE_OBSERVED；参数信号差异即使未改变主分类仍完整报告。大误差参考仅评价，不影响主选点。没有事后调门槛。
''',
 'REPRODUCTION.md':'''# AUX-D1登记顺序

依赖冻结现有Python3.10.11、NumPy1.26.4、SciPy1.15.3；所有输出在D盘。核心单计算进程，BLAS线程固定1（上限4）。不并行诊断，不安装或升级科学依赖。

一次执行：先aux_d1_preflight.py与设计提交A；随后aux_d1_generate.py calibration冻结阈值，aux_d1_generate.py test生成封存测试；aux_d1_execute.py单次推断；封存后aux_d1_evaluate.py、aux_d1_audit.py，再aux_d1_report.py只读报告。各生成／执行阶段有START守卫，存在时禁止重做；评价或报告读取错误仅允许另记纯恢复，不重新计算优化或改变选点。

主预算18000 s，科学累计21600 s；每起点100个残差调用，393216终点上限；RSS8GiB，新增本地2GiB，提交512MiB。科学预算已耗时在SCIENCE_CLOCK.json，各对checkpoint，单次残差守卫时间硬限制。RSS/存储每配对检查；达到限制当前已执行部分全部保存，未执行NOT_EVALUATED。所有终点使用无损float64 NPY，其字段schema随结果提供；CSV为17位有效数字，校准／选点使用原float64。

公开方位与报告轨迹在public_inputs，私有生成扰动在evaluation_only。估计模块无生成器／私有真值导入，U0输入只有主轴。校准器可用校准真值但优化0；测试评价真值必须在INFERENCE_SEAL后读取。冷复核不运行优化器，仅重建输入、Cartesian得分、接受／选择／去重与已有统计。连续全域证书、真实接收能力及R4验收不由本轮产生。B提交推送并核对后STOP。
'''}
 for n,t in docs.items():(OUT/n).write_text(t.strip()+'\n',encoding='utf-8',newline='\n')
 print('PREPARED',len(plan),'configurations;sourcehashesverified',flush=True)
def freeze():
 pref=read('PREFLIGHT.json');assert pref['PASS'];files=[p for p in ROOT.glob('aux_d1_*.py')]+[OUT/n for n in ['AUX_D1_TASK_SPEC_20261012_v1.json','AUX_D1_CODEX_TASK_20261012_v1.md','SOBOL_128_RAW.npy','SHARED_FIXED_126.npy','CONFIGURATION_PLAN.csv','SOURCE_MANIFEST.json','PREFLIGHT.json','STATISTICS_PREFLIGHT.json']];hashes={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files}
 dump('DESIGN_FREEZE.json',{'stage':'AUX-D1','task_id':'AUX_D1_PAIRED_BEARING_RESOURCE_INCREMENT_V1','task_SHA256':TASKSHA,'parent_SHA':PARENT,'design_commit':'CONTAINING_COMMIT','utc':utc(),'code_and_input_sha256':hashes,'runtime':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'BLAS_threads':1},'science_model':'Fulltask MD andTASK_SPEC;bounded4parameterobservation-onlyprocessor','no_test_or_calibration_execution_before_freeze':True,'execution_authorized':True,'next_stage_after_B_authorized':False,'preflight_summary':{'fits':pref['fits'],'gradient_checks':pref['gradient_checks']},'budgets':read('AUX_D1_TASK_SPEC_20261012_v1.json')['budgets'],'R4_percent':0})
 print('DESIGN_FROZEN',len(hashes),'identities',flush=True)
if __name__=='__main__':prepare() if sys.argv[1]=='prepare' else freeze()