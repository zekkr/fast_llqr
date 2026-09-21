from pathlib import Path
import json,hashlib
root=Path(__file__).resolve().parent
for name in ('FROZEN_SHA256.json','ARCHIVE_SHA256.json'):
    records=json.loads((root/name).read_text())
    for r in records:
        if hashlib.sha256((root/r['path']).read_bytes()).hexdigest()!=r['sha256']:
            raise SystemExit('Hash mismatch: '+r['path'])
    print('PASS',name,len(records),'files')
for r in json.loads((root.parent/'src/fortran/u11/PROVENANCE.json').read_text()):
    assert hashlib.sha256((root.parent/r['path']).read_bytes()).hexdigest()==r['sha256']
print('PASS archived U11 core fingerprints')
