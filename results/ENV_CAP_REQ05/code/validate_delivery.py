import pathlib,json,csv,hashlib,re,subprocess,time
R=pathlib.Path(r'D:\XiaomiMiMoProjects\WSL系统汇聚区目标定位');O=R/'results/ENV_CAP_REQ05';W=pathlib.Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_REQ05')
for src,dst in [('deliver.py','build_delivery.py'),('freeze.py','development_freeze.py')]: (O/'code'/dst).write_bytes((W/src).read_bytes())
p=O/'RESOURCE_AND_RESEARCH_DECISION.md';s=p.read_text(encoding='utf-8');s += '\n未知量投影、比值消元和增加方位基线本身均不作为本项目创新成果。可进一步研究的内容是：第一CZ下哪些可提取观测在源/标定/环境未知量同时存在时仍解除状态退化，以及该结论需要怎样的资源与不确定性约束。本轮只建立论证框架，尚未完成这一新方法的性能证明。\n';p.write_text(s,encoding='utf-8')
(O/'REPRODUCIBILITY.md').write_text('''# 复核入口

首先阅读GPT_SYNC_REQ05.md和RESOURCE_AND_RESEARCH_DECISION.md。历史证据索引为EVIDENCE_INDEX.csv（检出SHA256及父提交blob SHA256分列，CRLF差异单列），方法原文使用封存PDF；不重新下载数据。

Commit A冻结代码与18个无量纲矩阵。analyze.py使用NumPy正交投影；cold_audit.py不导入主程序，使用Fraction精确消元/逆矩阵及独立标量差分。ALGEBRA_RESULTS.json和EXACT_NULLSPACES.json提供全部矩阵和零空间，INDEPENDENT_COLD_AUDIT.json逐条记录124项。

run_once.py使用D盘排他启动标记；正式主分析和冷审计已各执行一次，不得在当前阶段删除标记或重试。开发脚本仅用于追溯冻结前构建过程，不是继续科研入口。build_delivery.py只格式化已保存结果/编写报告，不调用科学分析。后续外部审核如需独立重算，应使用另行批准的审计副本及独立输出位置，禁止覆盖本轮正式结果。

运行环境：Windows，Python 3.10，NumPy 1.26.4，psutil 5.9.8。PYTHONUTF8=1、PYTHONDONTWRITEBYTECODE=1；环境与临时目录来自D:/ProjectStorage/WSL-CZ/runtime-environment.ps1；依赖只读路径D:/ProjectStorage/WSL-CZ/T2_PREC3_FUNC2/runtime_dependencies。使用该依赖目录不运行FUNC2。全部计算为小型确定性代数控制，无传播求解、MC或实际观测生成。

冻结前两次文档/身份工具错误已记录：CRLF/LF检出差异与Windows GBK默认解码。均在正式Commit A前处理；未修改任何历史输入和科学阈值。预检与正式研究分别记账。

结果文件不自引用其所在提交SHA；最终执行SHA以Git中包含报告的Commit B及远端核对为准。DOC_DELIVERY_VALIDATION是交付文件/论断一致性检查，不是第二次科学冷审计。
''',encoding='utf-8')
checks=[]
def ck(n,b,detail=''):checks.append({'name':n,'PASS':bool(b),'detail':str(detail)})
freeze=json.loads((O/'DESIGN_FREEZE.json').read_text());cold=json.loads((O/'INDEPENDENT_COLD_AUDIT.json').read_text())
for q in freeze['code_and_spec']:ck('frozen:'+q['path'],hashlib.sha256((O/q['path']).read_bytes()).hexdigest()==q['sha256'])
def rows(n):return list(csv.DictReader((O/n).open(encoding='utf-8',newline='')))
e=rows('EVIDENCE_INDEX.csv');ids={r['id'] for r in e}
for q in e:
 ck('source_current:'+q['id'],hashlib.sha256((R/q['path']).read_bytes()).hexdigest()==q['sha256'])
 blob=subprocess.check_output(['git','-c','safe.directory='+str(R).replace('\\','/'),'show',freeze['parent_SHA']+':'+q['path']],cwd=R)
 ck('source_parent:'+q['id'],hashlib.sha256(blob).hexdigest()==q['parent_blob_sha256'])
obs=rows('OBSERVATION_EVIDENCE_LAYERS.csv');ck('six_observations',set(r['id'] for r in obs)=={'O1','O2','O3','O4','O5','O6'})
for r in obs:ck('evidence_ids:'+r['id'],set(r['evidence'].split(';'))<=ids)
m=rows('FIVE_STATE_INFORMATION_MATRIX.csv');ck('all_30_state_pairs',len(m)==30 and len(set((r['observation'],r['state']) for r in m))==30)
req=rows('PARTNER_DATA_REQUEST.csv');ck('19_requests_no_threshold',len(req)==19 and all(r['numeric_acceptance_threshold']=='NOT_SET' for r in req))
recs=json.loads((O/'NEXT_EXPERIMENT_RECOMMENDATIONS.json').read_text());ck('two_conditional_candidates',len(recs)==2 and all(r['status']=='PROPOSED_NOT_FROZEN_NUMERICAL_EXECUTION_NOT_AUTHORIZED' for r in recs))
for r in recs:ck('candidate_prerequisites:'+r['id'],set(r['required'])<=set(x['id'] for x in req))
ck('no_new_depth_accuracy',next(r['status'] for r in m if r['observation']=='O6' and r['state']=='z')=='STRUCTURAL_ZERO_IN_HORIZONTAL_BEARING_MODEL')
ck('cold_audit_pass',cold['PASS'] and cold['checks']==124 and cold['failed']==0)
order=json.loads((O/'EXECUTION_ORDER.json').read_text());ck('freeze_then_analysis_then_audit',time.mktime(time.strptime(freeze['freeze_UTC'],'%Y-%m-%dT%H:%M:%SZ'))-8*3600 < order['analysis_started_unix'] < order['audit_started_unix'])
# Above ordering also verified using exact Git commit timestamp, independent of local timezone.
ct=int(subprocess.check_output(['git','-c','safe.directory='+str(R).replace('\\','/'),'show','-s','--format=%ct',order['design_SHA']],cwd=R))
ck('git_design_precedes_science',ct<=order['analysis_started_unix'])
for p in O.glob('*.md'):
 for target in re.findall(r'\]\(([^)]+)\)',p.read_text(encoding='utf-8')):
  if not target.startswith(('http','mailto:')):ck('link:'+p.name+':'+target,(p.parent/target.split('#')[0]).exists())
size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());wsize=sum(p.stat().st_size for p in W.rglob('*') if p.is_file())
ck('storage_budget',size<16*1024**2 and size+wsize<128*1024**2,size+wsize)
changed=subprocess.check_output(['git','-c','safe.directory='+str(R).replace('\\','/'),'diff','--name-only',freeze['parent_SHA']],cwd=R).decode('utf-8').splitlines()
ck('only_req05_tracked_changes',all(p.startswith('results/ENV_CAP_REQ05/') for p in changed))
main=subprocess.check_output(['git','-c','safe.directory='+str(R).replace('\\','/'),'rev-parse','main'],cwd=R).decode().strip()
ck('main_unchanged',main=='a9a80191ee4d0f0cecd9db296198606b63a68080',main)
out={'scope':'DOCUMENT_DELIVERY_AND_IDENTITY_CHECK_NOT_SCIENCE_RERUN','PASS':all(c['PASS'] for c in checks),'checks':checks,'failed':sum(not c['PASS'] for c in checks),'new_repo_bytes_before_manifest':size,'new_D_total_bytes_before_manifest':size+wsize}
(O/'DOC_DELIVERY_VALIDATION.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps({'PASS':out['PASS'],'checks':len(checks),'failed':out['failed'],'bytes':size+wsize}))
if not out['PASS']:print([c for c in checks if not c['PASS']]);raise SystemExit(1)
# Manifest excludes itself; Git commit protects the manifest.
fs=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='FILE_HASHES_SHA256.csv':fs.append(dict(path=p.relative_to(O).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
with (O/'FILE_HASHES_SHA256.csv').open('w',encoding='utf-8',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['path','bytes','sha256']);w.writeheader();w.writerows(fs)
