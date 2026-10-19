import pathlib,json,csv,hashlib,importlib.util,subprocess,sys,time
R=pathlib.Path(r'D:\XiaomiMiMoProjects\WSL系统汇聚区目标定位');O=R/'results/ENV_CAP_REQ05'
def put(n,s):(O/n).write_text(s.strip()+'\n',encoding='utf-8')
s=(O/'EVIDENCE_README.md').read_text(encoding='utf-8').replace('表中各文件与父提交Git字节逐一核对','各来源仅允许Git检出导致的CRLF/LF差异，分别保存实际检出字节SHA256和父提交blob SHA256；内容逐一核对')
put('EVIDENCE_README.md',s)
def mod(name):
 spec=importlib.util.spec_from_file_location(name,O/'code'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
main=mod('analyze');cold=mod('cold_audit');spec=json.loads((O/'ALGEBRA_SPEC.json').read_text());r=main.evaluate(spec);checks=[]
for c,z in zip(spec['controls'],r['controls']):
 a,F,rank,n=cold.exact(c);checks.append(dict(id=c['id'],expected=c['expected_rank'],exact_rank=rank,numerical_rank=z['rank'],PASS=rank==z['rank']==c['expected_rank']))
pre={'scope':'PRE_FREEZE_SYNTHETIC_CODE_TEST_ONLY','checks':checks,'PASS':all(c['PASS'] for c in checks),'formal_runs':0,'historical_science_rerun':False,'development_events':[{'event':'input check rejected CRLF checkout against LF Git blob','resolution':'record both SHA256 identities; require exact equality modulo newline encoding only; no source modification'},{'event':'Chinese document reader used Windows default GBK','resolution':'PYTHONUTF8=1, UTF-8 output; no scientific code or thresholds changed'}]}
put('DEVELOPMENT_PRECHECK.json',json.dumps(pre,indent=2))
assert pre['PASS']
# Preserve preparation utilities and development provenance, no production imports.
for f in ['prepare.py','prepare_review.py']:(O/'code'/('development_'+f)).write_bytes((pathlib.Path(r'D:\ProjectStorage\WSL-CZ\ENV_CAP_REQ05')/f).read_bytes())
freeze={'stage':'ENV_CAP_REQ05','parent_SHA':'d658f660d6bd6c9b3a585f0e5b8812c49ee8b0c4','branch':'research/env-cap-req05','freeze_UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'formal_analysis_runs':1,'cold_audit_runs':1,'budget':{'wall_seconds_total':600,'process_tree_GiB':4,'D_data_MiB':128,'git_MiB':16,'acoustic_solver_calls':0,'MC':0,'new_recordings':0},'code_and_spec':[],'paper_online_status':'DOI access failed; archived primary PDFs and prior audited identity used','state_order':['r','theta','z','v','psi'],'performance_gate':'NOT_DEFINED_NOT_A_PERFORMANCE_STAGE','next_candidate_limit':2}
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='DESIGN_FREEZE.json':freeze['code_and_spec'].append({'path':p.relative_to(O).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
put('DESIGN_FREEZE.json',json.dumps(freeze,indent=2))
print(json.dumps({'precheck_PASS':pre['PASS'],'controls':len(checks),'frozen_files':len(freeze['code_and_spec'])}))
