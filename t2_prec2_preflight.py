#!/usr/bin/env python3
"""Pre-Commit-A gate: parent Git source bytes + immutable package SHA manifest."""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,os,datetime,shutil

def sha(blob):return hashlib.sha256(blob).hexdigest()
def cmd(cwd,*args):
 p=subprocess.run(args,cwd=cwd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,check=False)
 if p.returncode:raise RuntimeError(f'command failed: {args} rc={p.returncode}, stderr={p.stderr[:500]!r}')
 return p.stdout

def run(repo:Path,package:Path,out:Path):
    errors=[];checks=[];data=None
    def check(label,condition,actual='',expected=''):
        checks.append({'item':label,'PASS':bool(condition),'actual':str(actual),'expected':str(expected)})
        if not condition:errors.append(label)
    try:
        config=json.loads((package/'T2_PREC2_EXPERIMENT_SPEC_20261014_v1.json').read_text(encoding='utf-8'))
        manifest=json.loads((package/'T2_PREC2_FILE_MANIFEST_20261014_v1.json').read_text(encoding='utf-8'))
        for entry in manifest['files']:
            f=package/entry['filename'];exists=f.is_file()
            check('package_file_present/'+entry['filename'],exists)
            if exists:
                b=f.read_bytes();check('package_sha256/'+entry['filename'],sha(b)==entry['sha256'],sha(b),entry['sha256'])
                check('package_bytes/'+entry['filename'],len(b)==entry['bytes'],len(b),entry['bytes'])
        # One parent, no hidden mutations; inspect before producing any preflight output.
        actual=cmd(repo,'git','rev-parse','HEAD').decode().strip()
        check('HEAD_equals_frozen_parent',actual==config['parent_sha'],actual,config['parent_sha'])
        dirty=cmd(repo,'git','status','--porcelain=v1','--untracked-files=all').decode('utf-8',errors='replace')
        check('clean_independent_worktree',not dirty,dirty[:180],'no tracked/untracked files')
        for src in config['archive_sets']:
            for kind in ('env','mod','prt'):
                rel=src['paths'][kind];file=repo/rel
                present=file.is_file();check('source_exists/'+rel,present)
                if not present:continue
                local=file.read_bytes();expected=src['sha256'][kind]
                check('local_sha256/'+rel,sha(local)==expected,sha(local),expected)
                gitblob=cmd(repo,'git','show',config['parent_sha']+':'+rel)
                check('parent_sha256/'+rel,sha(gitblob)==expected,sha(gitblob),expected)
                check('worktree_equals_parent/'+rel,local==gitblob)
        try:
            import numpy,scipy
            data={'python':sys.version,'numpy':numpy.__version__,'scipy':scipy.__version__,
                  'platform':sys.platform,'visible_blas_threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')}}
            check('scipy_numpy_available',True)
            check('BLAS_thread_cap',all(os.environ.get(k) in (None,'1') for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')),data['visible_blas_threads'],'1 or unset')
        except Exception as e:check('scipy_numpy_available',False,str(e))
        check('reserve_storage_GiB',shutil.disk_usage(repo).free>=1*1024**3,shutil.disk_usage(repo).free,'>=1GiB')
    except Exception as exc:
        errors.append('PREFLIGHT_EXCEPTION:'+repr(exc))
    result={'stage':'T2-PREC2-PREFLIGHT','time_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'admitted':not errors,'errors':errors,'checks':checks,'environment':data,
            'next_step':'Commit A only if admitted; never auto execute','R4_percent':0}
    out.mkdir(parents=True,exist_ok=True)
    (out/'PREC2_PREFLIGHT.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PREC2 PREFLIGHT', 'PASS' if not errors else 'STOP',len(checks),'checks',len(errors),'errors')
    if errors:raise SystemExit(1)

if __name__=='__main__':
 pa=argparse.ArgumentParser();pa.add_argument('--repo',type=Path,required=True);pa.add_argument('--package',type=Path,required=True);pa.add_argument('--out',type=Path,required=True)
 a=pa.parse_args();run(a.repo.resolve(),a.package.resolve(),a.out.resolve())
