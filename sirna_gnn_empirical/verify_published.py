"""Read-only SHA-256 verification of the scientific GitHub contribution."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
manifest=json.loads((ROOT/'docs/sirna_gnn_empirical/publication/v1/files.json').read_text())
for name,expected in manifest.items():
    p=ROOT/name
    if not p.is_file():raise SystemExit(f'Missing published artifact: {name}')
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    if p.stat().st_size!=expected['bytes'] or h.hexdigest()!=expected['sha256']:
        raise SystemExit(f'Published artifact integrity mismatch: {name}')
print(f'Verified {len(manifest)} published scientific artifacts; no computation or training rerun.')
