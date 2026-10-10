"""Final byte/storage closure; no scientific calculation or changes."""
from hla_sp1_core import *
import subprocess

def main():
 checks=[]
 def check(name,ok):checks.append(dict(check=name,PASS=bool(ok)))
 design=read(OUT/'DESIGN_FREEZE.json')
 for n,s in design['code_sha256'].items():check('code:'+n,sha(ROOT/n)==s)
 for n,s in design['input_sha256'].items():check('input:'+n,sha(OUT/n)==s)
 for r in rows(OUT/'INPUT_PROVENANCE.csv'):check('historical:'+r['path'],sha(ROOT/r['path'])==r['working_sha256'])
 for r in read(OUT/'PROVIDER_PROVENANCE.json')['calls']:
  stem=f"sp1_f{int(r['frequency_hz'])}_n{r['mesh']}";check('MOD:'+stem,sha(LOCAL/'modal'/(stem+'.mod'))==r['mod_sha256']);check('ENV:'+stem,sha(LOCAL/'modal'/(stem+'.env'))==r['env_sha256'])
 for kind in ['calibration','test']:
  manifest=read(OUT/(kind.upper()+'_RECORD_MANIFEST.json'));check('public_index:'+kind,sha(manifest['public_index'])==manifest['public_index_sha256']);check('private_manifest:'+kind,sha(manifest['private_manifest'])==manifest['private_manifest_sha256']);p=read(LOCAL/kind/'PRIVATE_MANIFEST.json')['records'];unique={r['private_innovations']:r['private_sha256'] for r in p}
  for path,s in unique.items():check('innovations:'+path,sha(path)==s)
  df=read(OUT/(kind.upper()+'_DFT_FREEZE.json'));check('DFT:'+kind,sha(df['file'])==df['sha256'])
 infer=read(OUT/'INFERENCE_FREEZE.json')
 for n,s in infer['score_files'].items():check('scores:'+n,sha(LOCAL/n)==s)
 for n,s in infer['full_support_files'].items():check('support:'+n,sha(OUT/n)==s)
 check('threshold',sha(OUT/'CALIBRATION_THRESHOLDS.csv')==read(OUT/'THRESHOLD_FREEZE.json')['threshold_sha256'])
 localbytes=sum(p.stat().st_size for p in LOCAL.rglob('*') if p.is_file());repobytes=sum(p.stat().st_size for p in OUT.iterdir() if p.is_file())+sum(p.stat().st_size for p in ROOT.glob('hla_sp1_*.py'));check('local_storage_cap',localbytes<=16*1024**3);check('repo_added_cap',repobytes<=256*1024**2);budget=guard();check('wall_RSS_cap',budget[0]<=14400 and budget[1]<=8*1024**3);v=read(OUT/'VALIDATION.json');check('independent_validation',v['PASS']);dump(OUT/'FINAL_BYTE_AND_BUDGET_CLOSURE.json',{'time_utc':stamp(),'checks':checks,'all_PASS':all(r['PASS'] for r in checks),'local_bytes':localbytes,'repo_added_bytes':repobytes,'wall_s':budget[0],'scope':'immutability and operational budget, not scientific PASS'});print('FINAL CLOSURE',len(checks),'FAIL',sum(not r['PASS'] for r in checks),'local bytes',localbytes,'repo bytes',repobytes)
if __name__=='__main__':main()
