#!/usr/bin/env python3
"""Byte-for-byte integrity check for the original ENV-CAP task package.

No remote calls, no file repairs; prints a deterministic summary, exits nonzero on mismatch.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

MANIFEST = 'FILE_MANIFEST_SHA256.json'

def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def verify(root):
    root = Path(root).resolve()
    man = json.loads((root / MANIFEST).read_text(encoding='utf-8'))
    expected = {rec['path']: rec for rec in man['files']}
    current = {p.relative_to(root).as_posix(): p for p in root.rglob('*') if p.is_file() and p.name != MANIFEST and '__pycache__' not in p.parts}
    problems = []
    for name, rec in sorted(expected.items()):
        src = current.get(name)
        if src is None:
            problems.append(('MISSING', name))
        elif src.stat().st_size != rec['bytes']:
            problems.append(('SIZE_MISMATCH', name))
        elif digest(src) != rec['sha256']:
            problems.append(('SHA256_MISMATCH', name))
    for name in sorted(set(current) - set(expected)):
        if name.endswith('.pyc'):
            continue
        problems.append(('EXTRA_FILE', name))
    result = {'PASS': not problems, 'registered': len(expected), 'verified': len(expected)-sum(p[0] != 'EXTRA_FILE' for p in problems), 'problems': problems}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', default='.')
    args = parser.parse_args()
    try:
        result = verify(args.root)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        result = {'PASS': False, 'error': str(exc)}
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2))
    return 0 if result['PASS'] else 2

if __name__ == '__main__':
    sys.exit(main())
