#!/usr/bin/env python3
"""Read-only Git preflight for Codex execution checkout. No checkout/reset/clean commands."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

PARENT="dddad1f66e6d61668ecd5af38cad3e3af9fc9245"
BRANCH="research/env-cap-cap06"

def git(*args,cwd):
    p=subprocess.run(["git",*args],cwd=cwd,capture_output=True,text=True,check=False)
    if p.returncode:raise RuntimeError(f"git {' '.join(args)} error: {p.stderr.strip()}")
    return p.stdout.strip()

def main(argv=None):
    a=argparse.ArgumentParser()
    a.add_argument("--repo",required=True)
    a.add_argument("--pack",default=str(Path(__file__).resolve().parent.parent))
    p=a.parse_args(argv)
    repo=Path(p.repo)
    top=git("rev-parse","--show-toplevel",cwd=repo)
    sha=git("rev-parse","HEAD",cwd=repo)
    branch=git("branch","--show-current",cwd=repo)
    git("merge-base","--is-ancestor",PARENT,sha,cwd=repo)
    files=[]
    for row in csv.DictReader((Path(p.pack)/"11_HISTORICAL_EVIDENCE_INDEX.csv").open(encoding="utf-8",newline="")):
        missing=not (Path(top)/row["git_path"]).is_file()
        files.append({"id":row["id"],"missing":missing,"path":row["git_path"]})
    status=git("status","--porcelain",cwd=repo)
    d={"head":sha,"branch":branch,"expected_branch":BRANCH,"parent_ancestor":True,"missing_source_paths": [x for x in files if x["missing"]],"dirty_worktree":bool(status),"PASS":branch==BRANCH and not status and not any(x["missing"] for x in files)}
    print(json.dumps(d,ensure_ascii=False,indent=2));return 0 if d["PASS"] else 2

if __name__=="__main__":sys.exit(main())
