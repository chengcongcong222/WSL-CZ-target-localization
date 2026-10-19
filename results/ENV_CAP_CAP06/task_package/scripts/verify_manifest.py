#!/usr/bin/env python3
"""Check SHA256/length of every CAP06 task-pack content file listed in FILE_MANIFEST_SHA256.json."""
import hashlib
import json
from pathlib import Path
import sys


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()


def verify(folder):
    folder=Path(folder)
    d=json.loads((folder/"FILE_MANIFEST_SHA256.json").read_text(encoding="utf-8"))
    fails=[]
    paths={row["path"] for row in d["files"]}
    actual={p.relative_to(folder).as_posix() for p in folder.rglob("*") if p.is_file() and p.name!="FILE_MANIFEST_SHA256.json" and "__pycache__" not in p.parts}
    for row in d["files"]:
        p=folder/row["path"]
        if not p.is_file() or p.stat().st_size!=row["bytes"] or sha256(p)!=row["sha256"]:fails.append(row["path"])
    if paths!=actual:fails.append(f"file-set mismatch extra={sorted(actual-paths)} missing={sorted(paths-actual)}")
    print(json.dumps({"PASS":not fails,"count":len(paths),"failures":fails},ensure_ascii=False,indent=2))
    return not fails

if __name__=="__main__":sys.exit(0 if verify(Path(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parent.parent)) else 2)
