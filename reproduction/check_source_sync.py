#!/usr/bin/env python3
"""Check canonical research R/Fortran against a package source directory."""
from pathlib import Path
import argparse,hashlib,json
p=argparse.ArgumentParser();p.add_argument('package',type=Path);a=p.parse_args()
root=Path(__file__).resolve().parents[1]
pairs=[(f,a.package/'src'/f.name) for f in (root/'src/fortran/u11').glob('*.f90')]
pairs += [(f,a.package/'R'/f.name) for f in (root/'R/u11').glob('*.R')]
for source,target in pairs:
    if not target.exists() or source.read_bytes()!=target.read_bytes():raise SystemExit(f'MISMATCH: {target}')
for record in json.loads((root/'src/fortran/u11/PROVENANCE.json').read_text()):
    assert hashlib.sha256((root/record['path']).read_bytes()).hexdigest()==record['sha256']
print(f'PASS: {len(pairs)} package/research files identical; archived U11 hashes unchanged')
