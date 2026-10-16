from pathlib import Path
import json,hashlib,csv,subprocess,sys
R=Path.cwd(); O=R/'results/ENV_CAP_FULL_CAMPAIGN'; W=Path('D:/ProjectStorage/WSL-CZ/ENV_CAP'); P=O/'task_package/ENV_CAP_FULL_CAMPAIGN_20261016_v1'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def js(p,v): p.write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
records=[]
for row in csv.DictReader((P/'10_REFERENCE_SOURCE_INDEX.csv').open(encoding='utf-8-sig')):
 p=R/row['repository_relative_path']; b=subprocess.check_output(['git','-c',f'safe.directory={R.as_posix()}','show',f'a9a80191ee4d0f0cecd9db296198606b63a68080:{row["repository_relative_path"]}'])
 records.append({**row,'sha256_working':sha(p),'sha256_git_blob':hashlib.sha256(b).hexdigest(),'bytes':p.stat().st_size,'normalized_match':p.read_bytes().replace(b'\r\n',b'\n')==b.replace(b'\r\n',b'\n')})
 # Complete source text snapshot for researcher review only, historic content unchanged.
 if p.suffix in ['.md','.csv','.py','.json']:
  (W/'sources'/('reference_'+str(len(records)).zfill(2)+'.txt')).write_bytes(p.read_bytes())
assert all(r['normalized_match'] for r in records)
js(O/'E0_HISTORICAL_INPUT_IDENTITY.json',records)
js(O/'E0_ACQUISITION_SELECTION_FREEZE.json',{'independent_site':{'lat':25.,'lon':140.,'period':'annual','product':'WOA23 decav91C0 0.25 deg','reason':'fixed northwestern Pacific comparison; not selected on acoustic outputs','sal_temp_indices':{'lat':[459,460],'lon':[1279,1280]},'mirror':'https://thredds.rda.ucar.edu/thredds/dodsC/files/d285000/woa23_netcdf/','bathymetry_source':'ETOPO2022 NOAA CoastWatch ERDDAP at registered site'},'rule':'Attempt fixed site once per T/S/bathy request; no acoustic screening or site substitution if access/mask fails. WOA source-grid 17N118E monthly January over seasonal winter/annual are preregistered composite subcases of ONE provenance family. 18N114E diagnostic reconstruction cannot be counted as an independently sourced full-depth observed SSP.'})
js(O/'E0_PRECHECK.json',{'package_content_verified':32,'package_content_total':32,'package_tests_pass':3,'master_sha256':sha(P/'01_CODEX_MASTER_PROMPT.md'),'zip_sha256':'163f7f03819ddadc0f46c68256156f653447982364172b293e5cc7a46b0b8b2e','HEAD_and_remote_main':'a9a80191ee4d0f0cecd9db296198606b63a68080','tracked_changes_before':False,'historical_sources':28,'source_normalized_content_match':True,'pre_science_infrastructure_issue':'First WSL cp command rejected by shell quoting before execution; direct cp completed; no science executed.','C_writes_authorized':False})
(O/'.gitattributes').write_text('* -text\n',encoding='ascii',newline='\n')
print('28 source identities matched; acquisition selection frozen')

