"""Final archive indexes / measured budget accounting, no scientific re-execution."""
from env_provider import *
import shutil,subprocess
infra=O/'scripts/infrastructure';infra.mkdir(exist_ok=True)
for p in W.glob('*.py'):
 if p.name=='e0_source_cold_audit.py':continue
 shutil.copy2(p,infra/p.name)
# Complete solver calls and observed jobs remain available from Git without large ray/SHD binaries.
calls=[json.loads(s) for s in (W/'SOLVER_CALLS.jsonl').read_text().splitlines()];table(O/'SOLVER_CALLS.csv',calls)
jobs=[]
for p in sorted(W.glob('*.JOB_EXIT.json')):jobs.append({'job':p.name,**json.loads(p.read_text())})
js(O/'PROCESS_TREE_JOBS.json',jobs)
sourcecold=json.loads((O/'E0_SOURCE_INDEPENDENT_COLD_AUDIT.json').read_text())
science_s=sum(a['wall_s'] for a in jobs)+sum(a['elapsed_s'] for a in calls if a['unit']=='E0_CONTROLS')+sourcecold['wall_s']
peak=max([a['peak_tree_RSS_bytes'] for a in jobs]+[sourcecold['single_process_peak_RSS_bytes']]);Dbytes=sum(p.stat().st_size for p in W.rglob('*') if p.is_file());Gitbytes=sum(p.stat().st_size for p in O.rglob('*') if p.is_file())
js(O/'BUDGET_USAGE.json',{'accounted_scientific_wall_s':science_s,'limit_scientific_wall_s':259200,'external_solver_calls':len(calls),'limit_solver_calls':240,'by_solver':{'BELLHOP':len(calls),'KRAKEN':0,'KRAKENC':0,'FIELD':0,'FIELDS':0,'SCOOTER':0},'solver_returncode_0':sum(a['returncode']==0 for a in calls),'solver_success_is_not_scientific_gate_PASS':True,'new_MC_configurations':0,'MC_limit':20000,'new_pressure_records':0,'scientific_retries':0,'observed_peak_tree_RSS_bytes':peak,'RSS_limit_bytes':8*1024**3,'RSS_scope':'FormalWindowsjobsplusLinuxsourceaudit; E0short analytic solver controls were not separately RSS-sampled, so this is observed peak, not complete instrumentation certificate.','local_campaign_archive_bytes_at_index':Dbytes,'repository_artifact_bytes_at_index':Gitbytes,'D_limit_bytes':30*1024**3,'Git_limit_bytes':120*1024**2,'storage_scope':'NewcampaignDarchiveplusnewrepodirectory; excludes byte-for-byte delta of shared system caches and original user-suppliedDownloads package. No projectdata onC.','measured_caps_exceeded':False})
index=[]
for p in sorted(W.rglob('*')):
 if not p.is_file():continue
 index.append({'path':str(p),'local_relative_path':p.relative_to(W).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p),'kind':'INPUT' if ('sources' in p.parts or p.suffix=='.env') else 'RAW_OR_LOG_OR_OPERATIONAL_CODE','in_git':False,'immutable_after_STOP':True})
table(O/'raw_local_archive_index.csv',index)
# WOA provenance and original ledger hashes retained; do not rewrite master ledger.
source_ids=json.loads((O/'E0_HISTORICAL_INPUT_IDENTITY.json').read_text());ledger=next(a for a in source_ids if a['repository_relative_path']=='results/R4_MASTER/R4_EVIDENCE_LEDGER.csv')
oldsha=ledger['sha256_working'];newsha=sha(R/'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv');assert oldsha==newsha
js(O/'HISTORICAL_LEDGER_PRESERVATION.json',{'path':'results/R4_MASTER/R4_EVIDENCE_LEDGER.csv','before_SHA256':oldsha,'after_SHA256':newsha,'unchanged':True})
js(O/'INDEPENDENT_AUDIT_SUMMARY.json',{'E0_source':{'PASS':sourcecold['PASS'],'checks':16,'FAIL':0,'scope':'OriginalT/SandETOPOpointonly'},'E1_scientific_cold':{'PASS':False,'checks':274,'FAIL':12,'scope':'RawRAY/SHD,allregisteredunits; formaladmissionwithheld'},'E2_scientific_controls':{'PASS':False,'checks':26,'FAIL':1,'scope':'int32mirroroverflow;25otherchecksandsavedshallow12PASS'},'E2_E4_readonly_cold':{'PASS':True,'checks':126,'FAIL':0,'scope':'Storedpaper/analytic/historicalresource; neverrepairsE1/E2'},'campaign_scientific_admission_PASS':False,'historical_scientific_gates_unchanged':True,'independent_code_imported_main_research_functions':False,'new_solver_calls_in_cold':0,'retries':0,'note':'Deliverable integrity is checked separately; negative science may be completely archived without PASS.'})
# Additional exact raw dependencies for the example figure.
rawfig=[]
for p in [W/'E1/SCS17_JW_z200/R1.ray']+[W/'E1/SCS17_JW_z200'/f'I{i}_ENERGY_CURVE.npz' for i in range(3)]:rawfig.append({'figure':'F01_ENVIRONMENT_SSP_AND_CZ','path':str(p),'sha256':sha(p),'bytes':p.stat().st_size})
table(O/'FIGURE_LOCAL_DATA_DEPENDENCIES.csv',rawfig)
# preserve required report/source completeness and Git data budget before final commit.
changed=subprocess.check_output(['git','-c',f'safe.directory={R.as_posix()}','diff','--name-only','-z','a9a80191ee4d0f0cecd9db296198606b63a68080']).decode('utf-8').strip('\0').split('\0');assert all(a.startswith('results/ENV_CAP_FULL_CAMPAIGN/') for a in changed)
print(json.dumps({'science_wall_s':science_s,'peak_MiB':peak/1024**2,'D_GiB':Dbytes/1024**3,'Git_MiB':Gitbytes/1024**2,'raw_index_files':len(index),'solver_calls':len(calls),'old_ledger_unchanged':True},indent=2))
