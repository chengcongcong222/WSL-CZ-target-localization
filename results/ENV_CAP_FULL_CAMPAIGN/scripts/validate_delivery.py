"""Final delivery validation only. Does not run scientific analyses."""
from env_provider import *
import ast,struct,shutil,subprocess
import pypdfium2 as pdf
# Preserve all operational preparation code; no import/execution of it here.
infra=O/'scripts/infrastructure';infra.mkdir(exist_ok=True)
for p in W.glob('*.py'):
 if p.name!='e0_source_cold_audit.py':shutil.copy2(p,infra/p.name)
index=[]
for p in sorted(W.rglob('*')):
 if p.is_file():index.append({'path':str(p),'local_relative_path':p.relative_to(W).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p),'kind':'INPUT' if ('sources' in p.parts or p.suffix=='.env') else 'RAW_OR_LOG_OR_OPERATIONAL_CODE','in_git':False,'immutable_after_STOP':True})
table(O/'raw_local_archive_index.csv',index)
required='E0_BASELINE_AUDIT.md E0_FREEZE.json E1_ENV_SOURCE_INDEX.csv E1_CZ_ENVIRONMENTS.csv E1_CZ_GATE.json E1_CZ_PHYSICS_EVIDENCE.md E1_CZ_RANGE_SHAPE_DIAGNOSTICS.csv E1_CZ_NEGATIVE_EXAMPLES.csv E2_PAPER_METHOD_IDENTITY.md E2_YANG_ERRATUM_FORMULA_AUDIT.md E2_YANG_MMAR_COMPARISON.csv E2_MODE_DENSITY_AND_APERTURE.csv E2_METHOD_GATE.json E3_PERFORMANCE_BY_CASE.csv E3_OBSERVABLE_EXTRACTION.csv E3_NUISANCE_PROFILED_INFORMATION.csv E3_ENV_SENSITIVITY.csv E3_GATE.json E4_RESOURCE_INVERSION.csv E4_NEEDS_FROM_CLIENT.csv E4_PARETO_RESOURCE_FEASIBILITY.csv E4_GATE.json ENV_CAP_EVIDENCE_LEDGER.csv ENV_CAP_FINAL_DECISION.json ENV_CAP_FINAL_SCIENTIFIC_REPORT.md GPT_SYNC_ENV_CAP_FINAL.md APPLICATION_USABLE_CLAIMS.md APPLICATION_FORBIDDEN_CLAIMS.md REVIEW_ENTRYPOINT.md RUN_STATE.json SESSION_HANDOFF.md INDEPENDENT_AUDIT_SUMMARY.json STOP.json raw_local_archive_index.csv'.split()
checks=[]
def ck(name,ok,detail=''):checks.append({'check':name,'PASS':bool(ok),'detail':str(detail),'scope':'artifact_integrity_only_not_science_gate'})
for name in required:ck('required_'+name,(O/name).is_file() and (O/name).stat().st_size>0)
for p in (O/'scripts').rglob('*.py'):
 try:ast.parse(p.read_text(encoding='utf-8-sig'));good=True
 except SyntaxError:good=False
 ck('syntax_'+p.relative_to(O).as_posix(),good)
for unit,freeze in [('E2','E2_UNIT_FREEZE.json'),('E4','E4_UNIT_FREEZE.json')]:
 f=json.loads((O/freeze).read_text(encoding='utf-8-sig'));ck(unit+'_frozen_code_unchanged',sha(O/f'scripts/{"e2_controls.py" if unit=="E2" else "e4_resources.py"}')==f['code_sha256'])
for a,h in json.loads((O/'E0_CODE_FREEZE.json').read_text(encoding='utf-8-sig')).items():ck('E0_frozen_'+a,sha(O/a)==h)
for p in sorted((O/'figures').glob('*.pdf')):
 doc=pdf.PdfDocument(str(p));text=''.join(doc[i].get_textpage().get_text_range() for i in range(len(doc)));ck('PDF_'+p.name,len(doc)==1 and len(text)>50,{'pages':len(doc),'text_chars':len(text)});doc.close()
ck('seven_PDFs_and_PNGs',len(list((O/'figures').glob('*.pdf')))==7 and len(list((O/'figures').glob('*.png')))==7)
for a in csv.DictReader((O/'FIGURE_DATA_PROVENANCE.csv').open()):
 paths=a['source_paths'].split(';');hashes=a['source_SHA256'].split(';');ck('figure_data_SHA_'+a['figure'],all(sha(O/p if (O/p).exists() else R/p)==h for p,h in zip(paths,hashes)))
state=json.loads((O/'RUN_STATE.json').read_text(encoding='utf-8-sig'));stop=json.loads((O/'STOP.json').read_text(encoding='utf-8-sig'));ck('STOP_true',state['STOP'] and stop['STOP']);ck('all_phases_terminal',all(x in ['PASS','LIMITED','FAIL','NOT_EVALUATED'] for x in state['phases'].values()))
ck('old_ledger_unchanged',json.loads((O/'HISTORICAL_LEDGER_PRESERVATION.json').read_text(encoding='utf-8-sig'))['after_SHA256']==sha(R/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv'))
changes=subprocess.check_output(['git','-c',f'safe.directory={R.as_posix()}','diff','--name-only','-z','a9a80191ee4d0f0cecd9db296198606b63a68080']).decode('utf-8').strip('\0').split('\0');ck('no_historical_tracked_changes',all(p.startswith('results/ENV_CAP_FULL_CAMPAIGN/') for p in changes))
ck('scientific_failures_not_hidden',not json.loads((O/'E1_INDEPENDENT_COLD_AUDIT.json').read_text(encoding='utf-8-sig'))['PASS'] and not json.loads((O/'E2_METHOD_GATE.json').read_text(encoding='utf-8-sig'))['analytic_and_stored_control_PASS'] and not json.loads((O/'INDEPENDENT_AUDIT_SUMMARY.json').read_text(encoding='utf-8-sig'))['campaign_scientific_admission_PASS'])
ck('R4_zero',json.loads((O/'ENV_CAP_FINAL_DECISION.json').read_text(encoding='utf-8-sig'))['R4_percent']==0)
Dbytes=sum(a['bytes'] for a in index);Gitbytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());ck('storage_budgets_measured_scope',Dbytes+Gitbytes<=30*1024**3 and Gitbytes<=120*1024**2 and max(p.stat().st_size for p in O.rglob('*') if p.is_file())<=20*1024**2)
P=O/'task_package/ENV_CAP_FULL_CAMPAIGN_20261016_v1';manifest=json.loads((P/'FILE_MANIFEST_SHA256.json').read_text(encoding='utf-8-sig'));ck('all32_original_package_bytes_unchanged',all(sha(P/a['path'])==a['sha256'] and (P/a['path']).stat().st_size==a['bytes'] for a in manifest['files']))
result={'artifact_integrity_PASS':all(c['PASS'] for c in checks),'checks':len(checks),'FAIL':sum(not c['PASS'] for c in checks),'scientific_admission_PASS':False,'scope':'Required outputs, code syntax and frozen identity, sevenPDF/PNG pairs, figure sourceSHA, STOP and negative-gate preservation. No scientific reruns.','rows':checks}
js(O/'DELIVERY_VALIDATION.json',result)
if not result['artifact_integrity_PASS']:print(json.dumps([c for c in checks if not c['PASS']],indent=2));raise SystemExit(2)
# SHA CSV excludes itself to avoid recursive hashing; every other delivered file is included.
records=[]
for p in sorted(O.rglob('*')):
 if p.is_file() and p.name!='FILE_HASHES_SHA256.csv':records.append({'path':p.relative_to(O).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)})
table(O/'FILE_HASHES_SHA256.csv',records)
print(json.dumps({'artifact_integrity_PASS':True,'checks':len(checks),'science_PASS':False,'indexed_delivered_files':len(records),'indexed_raw_files':len(index),'archive_GiB':Dbytes/1024**3,'repo_MiB':Gitbytes/1024**2},indent=2))

