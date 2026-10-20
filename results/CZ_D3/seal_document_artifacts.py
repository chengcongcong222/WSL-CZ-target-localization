"""Seal staged publication artifacts; record checkout bytes separately from Git bytes."""
from pathlib import Path
import hashlib,json,subprocess
out=Path(__file__).resolve().parent
root=out.parents[1]
entries=[]
for p in sorted(out.rglob('*')):
    if p.is_file() and p.name!='DELIVERY_MANIFEST.json':
        rel=p.relative_to(root).as_posix()
        data=subprocess.check_output(['git','show',':'+rel],cwd=root)
        working=p.read_bytes()
        entries.append({'path':p.relative_to(out).as_posix(),'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data),'working_tree_sha256':hashlib.sha256(working).hexdigest(),'working_tree_bytes':len(working)})
(out/'DELIVERY_MANIFEST.json').write_text(json.dumps({'scope':'sha256 and bytes refer to staged Git blobs; working tree fields separately preserve line-ending differences; commit seals this manifest','files':entries,'total_git_bytes_without_manifest':sum(r['bytes'] for r in entries)},indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(len(entries),'Git artifact hashes sealed')