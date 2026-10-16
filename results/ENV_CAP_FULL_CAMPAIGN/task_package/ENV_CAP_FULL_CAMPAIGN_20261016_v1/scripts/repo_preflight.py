#!/usr/bin/env python3
"""Readonly environment + Git-parent + evidence-path guard for Codex handoff.

Checks Git objects at the pinned parent without changing the user's checkout.
Does not claim to audit scientific outputs or external environmental datasets.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import subprocess
import sys

DEFAULT_PARENT='a9a80191ee4d0f0cecd9db296198606b63a68080'

def run(repo, *args):
    p=subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True,encoding='utf-8',errors='replace')
    return p.returncode,p.stdout.strip(),p.stderr.strip()


def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--repo',required=True)
    pa.add_argument('--expected-parent',default=DEFAULT_PARENT)
    pa.add_argument('--source-index',default=None)
    p=pa.parse_args()
    repo=Path(p.repo).resolve()
    source=Path(p.source_index).resolve() if p.source_index else Path(__file__).resolve().parents[1]/'10_REFERENCE_SOURCE_INDEX.csv'
    checks=[];failures=[]
    rc, head,err=run(repo,'rev-parse','HEAD')
    checks.append({'check':'git_HEAD_available','PASS':rc==0,'value':head,'detail':err})
    rc0,real,err0=run(repo,'cat-file','-t',p.expected_parent)
    checks.append({'check':'frozen_parent_is_commit','PASS':rc0==0 and real=='commit','value':real,'detail':err0})
    rcr,origin,er=run(repo,'rev-parse','refs/remotes/origin/main')
    if rcr==0:
        checks.append({'check':'local_origin_main_matches_frozen_parent','PASS':origin==p.expected_parent,'value':origin,'detail':'Local tracking ref; remote freshness not certified.'})
    else:
        checks.append({'check':'origin_main_local_tracking_available','PASS':False,'value':'NOT_AVAILABLE','detail':er})
    if len(p.expected_parent)!=40:
        failures.append({'error':'MALFORMED_SHA1','value':p.expected_parent})
    try:
        with source.open(encoding='utf-8-sig',newline='') as f:
            rows=list(csv.DictReader(f))
    except OSError as e:
        print(json.dumps({'PASS':False,'error':str(e)},indent=2));return 2
    for item in rows:
        path=item['repository_relative_path']
        exists=(repo/path).is_file()
        rcb,_,errb=run(repo,'cat-file','-e',f'{p.expected_parent}:{path}')
        check={'check':'git_source_at_parent','path':path,'mandatory':item['must_exist']=='YES','PASS':rcb==0,'checkout_present':exists}
        checks.append(check)
        if not exists and item['must_exist']=='YES':
            failures.append({'error':'MANDATORY_SOURCE_NOT_CHECKED_OUT','path':path})
    mandatory_fail=[x for x in checks if not x['PASS'] and (x.get('mandatory',True))]
    failures.extend({'error':'CHECK_FAILED','check':x['check'],'path':x.get('path')} for x in mandatory_fail)
    result={'PASS':not failures,'expected_parent':p.expected_parent,'HEAD':head,'source_rows':len(rows),'checks':checks,'failures':failures,'warning':'Read-only Git file presence, not scientific or source SHA authentication. Verify authoritative remote via git fetch/ls-remote; report local EOL differences separately.'}
    print(json.dumps(result,ensure_ascii=False,indent=2))
    return 0 if result['PASS'] else 2

if __name__=='__main__':
    sys.exit(main())
