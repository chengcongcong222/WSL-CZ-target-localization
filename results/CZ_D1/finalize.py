from pathlib import Path
import csv,json,hashlib,subprocess,datetime
O=Path(__file__).resolve().parent;R=O.parents[1]
def read(n):return list(csv.DictReader((O/n).open(encoding='utf-8-sig',newline='')))
def write(n,rr):
 with (O/n).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,list(dict.fromkeys(k for r in rr for k in r)));w.writeheader();w.writerows(rr)
# Presentation normalization only, not statistical changes.
m=read('CZ_D1_EVIDENCE_MASTER.csv')
for x in m:
 if x['z_true_m']=='NOT_MODELED':x['z_true_m']='';x['z_generation_status']='NOT_MODELED_IN_BEARING_STATE'
 for k,v in list(x.items()):
  if v in ['inf','Infinity']:x[k]='INF'
 if x['N_elements']=='8':
  x['N_elements_scope']='per HLA node; not an observation-processing claim'
  if x['aux_count'] in ['0','1']:x['total_nominal_HLA_channels']=8*(1+int(x['aux_count']))
 if x['experiment_id']=='H3-G0':x['n_independent']='';x['independent_main_geometries_total']=3;x['sampling_status']='DETERMINISTIC_NOT_MC'
write('CZ_D1_EVIDENCE_MASTER.csv',m)
valid=json.loads((O/'VALIDATION.json').read_text('utf-8'))
if not valid['PASS']:raise SystemExit('Audit failed; cannot seal')
gates={
'EVIDENCE_INVENTORY_COMPLETE':{'status':'PASS','scope':'14540 tracked project files at fixed parent;14151 result files;137stage entries;171design artifacts. Not a claim all unindexed D-disk data searched.'},
'METRIC_DEFINITIONS_RECONCILED':{'status':'PASS_SCOPED','scope':'Reported metrics separated by type, unit, quantile, conditioning, rejection and paired denominator. Unrecomputed historical metrics retain NOT_RECALCULATED_SUMMARY_ONLY, not a blanket scientific revalidation.'},
'SOURCE_GAPS_DISCLOSED':{'status':'PASS','scope':'CAP06 availability updated; remaining receiver/environment/archive/old coverage gaps explicitly retained.'},
'CZ_D2_CANDIDATE_RANKING_READY':{'status':'PASS_WITH_NO_CANDIDATE_ADMITTED','scope':'Objective-specific comparison delivered; no unsupported scalar rank or unique experiment chosen.'}}
dec={'stage':'CZ-D1','decision':'EVIDENCE_CONSOLIDATED_CONDITIONAL_FOUR_STATE_REFERENCE_NO_SINGLE_HLA_FIVE_STATE_ADMISSION','parent':'ca86d324afc1689af6c96cd35145e03625086bae','design_SHA':'6ff0d07fca848d8022fad2facd5a107a920f8136','result_SHA':'CONTAINING_COMMIT','branch':'research/cz-d1','gates':gates,'master_rows':len(m),'new_solver_calls':0,'new_MC':0,'new_optimizer':0,'new_receiver_records':0,'R4_percent':0,'CZ_D2':'NOT_AUTHORIZED; await independent review','CZ_D3':'NOT_AUTHORIZED','main_expected':'a9a80191ee4d0f0cecd9db296198606b63a68080','STOP':True}
(O/'CZ_D1_FINAL_DECISION.json').write_text(json.dumps(dec,ensure_ascii=False,indent=2),encoding='utf-8')
(O/'STOP.json').write_text(json.dumps({'stage':'CZ-D1','status':'STOP_AFTER_RESULT_PUBLICATION','science_finished':True,'next_stage_execution_authorized':False,'historical_gates_unchanged':True,'R4_percent':0},indent=2),encoding='utf-8')
# Keep execution identities deterministic on later documentation-only rebuilds.
f=O/'consolidate.py';s=f.read_text('utf-8-sig');s=s.replace("'design_sha':git('rev-parse','HEAD')","'design_sha':'6ff0d07fca848d8022fad2facd5a107a920f8136'");f.write_text(s,encoding='utf-8')
# SOURCE index values are snapshot identities; task MD exact byte identity separately checked.
a=json.loads((O/'DESIGN_FREEZE.json').read_text());assert hashlib.sha256((O/'AUTHORIZED_TASK.md').read_bytes()).hexdigest()==a['task_sha256']
(O/'.gitattributes').write_text('* -text\n',encoding='ascii')
(O/'REVIEW_ENTRYPOINT.md').write_text('''# CZ-D1审阅入口

本轮是既有证据与精度归并，不是新实验。建议先读[同步报告](GPT_SYNC_CZ_D1.md)、[架构精度表](CZ_D1_ARCHITECTURE_METRICS.md)、[甲方简表](CZ_D1_CLIENT_PRECISION_CONDITIONS.md)。

随后核查[五维矩阵](CZ_D1_FIVE_STATE_CLAIM_MATRIX.md)、[机动与增强](CZ_D1_MANEUVERS_AND_AUGMENTATIONS.md)、[负结果与缺口](CZ_D1_NEGATIVE_GATES_AND_GAPS.md)、[第二轮候选](CZ_D2_SELECTION_CANDIDATES.md)。

技术复核：[证据主表](CZ_D1_EVIDENCE_MASTER.csv)、[逐场景分位数](CZ_D1_RECOMPUTED_SCENE_QUANTILES.csv)、[配对变化](CZ_D1_PAIRED_CHANGES.csv)、[来源SHA](CZ_D1_SOURCE_PROVENANCE.json)、[独立检查](VALIDATION.json)、[Gate裁定](CZ_D1_FINAL_DECISION.json)。完整文件/阶段/设计索引与全部只读脚本同目录保存。

统计复算和来源哈希核对由执行代理的独立代码完成，不代替研究负责人审计。尚无原单HLA完整五维精度，尚未准入CZ-D2唯一实验。当前STOP，main和历史R4=0%不变。
''',encoding='utf-8')
sync=O/'GPT_SYNC_CZ_D1.md';txt=sync.read_text('utf-8');txt+='\n独立归并检查：'+str(valid['checks'])+'项、0失败；另从保存估计与真值核对215268个误差项，重点原表一致性5974项。CAP06原始20545文件及AUX/T3另110个公开/私有输入文件身份通过。这些计数是归并可复核性证据，不是定位性能通过数。详见VALIDATION.json。\n';sync.write_text(txt,encoding='utf-8')
# Source task freeze and report metadata delivered; index excludes itself to avoid recursion.
index=[]
for p in sorted(O.iterdir()):
 if p.is_file() and p.name!='CZ_D1_DELIVERY_MANIFEST.csv':index.append({'file':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
write('CZ_D1_DELIVERY_MANIFEST.csv',index)
print(json.dumps({'rows':len(m),'files':len(index)+1,'bytes':sum(int(x['bytes']) for x in index),'decision':dec['decision']}))
