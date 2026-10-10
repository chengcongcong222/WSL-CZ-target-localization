"""Documentation-only check: archived literal statistics and bytes, no science imports."""
import csv,hashlib,json,math,pathlib,re,statistics,subprocess,ast
OUT=pathlib.Path(__file__).resolve().parent;ROOT=OUT.parent.parent
checks=[]
def check(name,ok,detail=''):checks.append({'check':name,'PASS':bool(ok),'detail':detail})
def sha(b):return hashlib.sha256(b).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def rows(p):return list(csv.DictReader(p.read_text(encoding='utf-8-sig').splitlines()))
def git(*args):return subprocess.check_output(['git',*args],cwd=ROOT)
m=load(OUT/'SOURCE_MANIFEST.json');parent=m['parent_SHA'];sid={s['source_id']:s for s in m['sources']}
for s in m['sources']:
 b=git('show',parent+':'+s['path']);raw=(ROOT/s['path']).read_bytes()
 check('source:'+s['source_id']+':git_identity',sha(b)==s['git_blob_byte_sha256'] and git('rev-parse',parent+':'+s['path']).decode().strip()==s['git_blob_id'])
 if s['source_id']!='E22':check('source:'+s['source_id']+':working_unchanged',sha(raw)==s['working_raw_byte_sha256'])
 else:check('ledger_prefix',len(raw)>=m['ledger_prefix_bytes'] and sha(raw[:m['ledger_prefix_bytes']])==m['ledger_prefix_sha256'])
 if 'archived_status_fields' in s:
  x=json.loads(b.decode('utf-8-sig'));check('status_identity:'+s['source_id'],all(x.get(k)==v for k,v in s['archived_status_fields'].items()))
r3=load(ROOT/sid['E01']['path']);check('R3_canonical',r3['final_nominal_candidate_set']=={'r':50,'theta':0,'v':2,'psi':[4,5],'units':{'r':'km','theta':'deg','v':'m/s','psi':'deg'}})
rr=rows(ROOT/sid['E02']['path']);check('R3_all_survivor_labels',len(rr)==12 and all(float(r['r'])==50 and float(r['theta'])==0 and float(r['v'])==2 and float(r['psi']) in [4,5] for r in rr))
ev=rows(OUT/'EVIDENCE_MATRIX.csv');h2=rows(ROOT/sid['E39']['path']);sp1=rows(ROOT/sid['E11']['path']);aux=rows(ROOT/sid['E30']['path']);vert=rows(ROOT/sid['E36']['path'])
for r in ev:
 check('matrix_source:'+r['stage']+':'+r['statistic_name'],r['source_id'] in sid and all(r[k] for k in r))
 expected=None
 if r['source_id']=='E39':
  parts=dict(x.split('=',1) for x in r['exact_field'].split(';') if '=' in x);sel=[x for x in h2 if x['method']==parts['method'] and float(x['sigma'])==float(parts['sigma'])];assert len(sel)==8
  if r['statistic_name']=='accepted_configurations':expected=sum(int(x['accepted']) for x in sel)
  else:
   field=r['exact_field'].split(';')[-1];vals=[float(x[field]) for x in sel];expected=(statistics.mean(vals) if field.endswith('median') else max(vals))*(100 if r['unit']=='%' else 1)
 if r['source_id']=='E11':
  parts=r['exact_field'].split(';');sel=dict(x.split('=',1) for x in parts[:-1]);rec=next(x for x in sp1 if all(x[k]==v for k,v in sel.items()));expected=float(rec[parts[-1]])
 if r['source_id']=='E30':
  anchor=r['exact_field'].split(';')[0].split('=')[1];rec=next(x for x in aux if x['anchor']==anchor);expected=float(rec[r['exact_field'].split(';')[1]])*(100 if r['unit']=='%' else 1)
 if r['source_id']=='E36':
  parts=r['exact_field'].split(';');sel=dict(x.split('=',1) for x in parts[:-1]);rec=next(x for x in vert if all(x[k]==v for k,v in sel.items()));expected=float(rec['B1_to_B2'])
 if expected is not None:
  actual=float(r['value']);check('matrix_value:'+r['stage']+':'+r['target_parameter']+':'+r['statistic_name'],actual==expected if math.isinf(expected) else math.isclose(actual,expected,rel_tol=2e-11,abs_tol=1e-12))
# Same-definition stored-summary consistency; no individual trajectory or likelihood rebuild.
case=rows(ROOT/sid['E31']['path']);print('CASE_COLUMNS',list(case[0]))
for a in aux:
 for field in ['range_error','bearing_error_deg','speed_error','heading_error_deg']:
  key=field+'_P95';c=[r for r in case if r['anchor']==a['anchor']];check('aux_saved_case_max:'+a['anchor']+field,math.isclose(max(float(r[key]) for r in c),float(a[field+'_worst_case_P95']),rel_tol=1e-12))
check('H3_B1_B2_bounds',math.isclose(min(float(r['B1_to_B2']) for r in vert),2.340748922039478) and math.isclose(max(float(r['B1_to_B2']) for r in vert),4.4743040030997))
h1=rows(ROOT/sid['E40']['path']);sel=[r for r in h1 if r['method']=='M1' and abs(float(r['sigma'])-.25)<1e-12 and r['source'] in ['S1','S2']]
check('H1_main_retention_bounds',len(sel)==12 and min(float(r['truth_retention']) for r in sel)==.90625 and max(float(r['truth_retention']) for r in sel)==1)
check('H1_M0_empty_ranking_distinct',all(float(r['empty_rate'])==1 and float(r['median_point_error_r_km'])==0 for r in h1 if r['method']=='M0' and abs(float(r['sigma'])-.25)<1e-12 and r['source'] in ['S1','S2']))
check('H1_common_source_pairing',load(ROOT/sid['E05']['path'])['all_common_source_sets_equal'])
check('E2_final_stop',load(ROOT/sid['E23']['path'])['decision']=='E2_FORMULA_REPAIR_FAILED' and load(ROOT/sid['E23']['path'])['development_regressions_executed']==0)
d=load(OUT/'CLOSEOUT_DECISION.json');check('authorizations_false',d['next_numeric_execution_authorized'] is False and d['vertical_numeric_execution_authorized'] is False and d['all_methods_exhausted'] is False and d['physical_no_go_claim'] is False)
check('science_zero_R4_stop',all(d[k]==0 for k in ['new_propagation_calls','new_observations_or_noise','new_optimizer_or_FIM_runs','new_MC','new_recordings','R4_percent','R4_delta']) and d['STOP'] is True)
req=['SINGLE_HLA_STAGE_CLOSEOUT.md','THREE_TIER_EVIDENCE_AND_CONDITIONS.md','CLIENT_DISCUSSION_BRIEF.md','NEXT_AUX_D1_DESIGN_DRAFT.md','VERTICAL_OPTION_DESIGN_NOTE.md','EVIDENCE_MATRIX.csv','METHOD_COVERAGE_DELTA.csv','SOURCE_MANIFEST.json','CONFLICTS_AND_GAPS.md','CLOSEOUT_DECISION.json','GPT_SYNC.md','NEXT_CHAT_HANDOFF.md']
check('required_deliverables',all((OUT/p).is_file() for p in req))
lengths={}
for name in req:
 p=OUT/name;t=p.read_text(encoding='utf-8');check('UTF8:'+name,'\ufffd' not in t and '\x00' not in t)
 if p.suffix=='.md':
  lengths[name]=len(re.findall('[\u4e00-\u9fff]',t))
  for target in re.findall(r'\]\(([^)]+)\)',t):
   if '://' not in target:check('local_link:'+name+':'+target,(p.parent/target.split('#')[0]).is_file())
   if '#e' in target:check('source_anchor:'+target,target.rsplit('#',1)[1].upper() in sid)
for name in ['NEXT_AUX_D1_DESIGN_DRAFT.md','VERTICAL_OPTION_DESIGN_NOTE.md']:
 t=(OUT/name).read_text(encoding='utf-8');check('proposal_only:'+name,all(s in t for s in ['PROPOSED_NOT_FROZEN','NUMERICAL_EXECUTION_NOT_AUTHORIZED','false']) and not re.search(r'(?im)^\s*(python|kraken|field|wsl)\s+',t))
check('GPT_sync_short',lengths['GPT_SYNC.md']<=900)
check('task_hash',sha((OUT/'HLA_CLOSE1_CODEX_TASK_20261010_v1.md').read_bytes())==d['task_SHA256']==m['task_SHA256'])
check('independent_checker_stdlib_only',not any(isinstance(n,(ast.Import,ast.ImportFrom)) and any(w in ast.unparse(n) for w in ['numpy','scipy','r4_','hla_h1','hla_h2']) for n in ast.walk(ast.parse(pathlib.Path(__file__).read_text(encoding='utf-8')))))
ledger=(ROOT/sid['E22']['path']).read_bytes();suffix=ledger[m['ledger_prefix_bytes']:].decode('utf-8');check('ledger_one_doc_row',len(suffix.strip().splitlines())==1 and all(s in suffix for s in ['DOCUMENTATION_AND_PLAN_COMPLETE','scientific_runs=0','R4_delta=0']))
changed=git('diff','--name-only',parent).decode('utf-8').splitlines();check('no_historical_science_change',all(p==sid['E22']['path'] or p.startswith('results/R4_SINGLE_HLA_CLOSEOUT_AND_THREE_TIER_PLAN/') for p in changed))
result={'stage':'HLA-CLOSE1','validation_scope':'DOCUMENTATION_ONLY; stored-statistics extraction, source bytes and prospective boundaries. No historical scientific rerun.','PASS':all(c['PASS'] for c in checks),'checks':checks,'source_records':len(sid),'evidence_rows':len(ev),'document_chinese_character_counts':lengths,'new_scientific_runs':0,'historical_full_audits_rerun':False,'style_skill':'write-like-me: retrieval tools unavailable; used provided Chinese task references; no fetched writing-profile claim','source_gaps_retained':d['source_gaps'],'parent_SHA':parent}
(OUT/'DOCUMENTATION_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps({'PASS':result['PASS'],'FAILS':[c for c in checks if not c['PASS']],'lengths':lengths},ensure_ascii=False));raise SystemExit(0 if result['PASS'] else 1)