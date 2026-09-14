"""Read-only SHA-256 verification of the current scientific publication."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
manifest={}
# Later manifests override only explicitly revised documentation/assets.
for version in ['v1','v3']:
    manifest.update(json.loads((ROOT/f'docs/sirna_gnn_empirical/publication/{version}/files.json').read_text()))
for name,expected in manifest.items():
    p=ROOT/name
    if not p.is_file():raise SystemExit(f'Missing published artifact: {name}')
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    if p.stat().st_size!=expected['bytes'] or h.hexdigest()!=expected['sha256']:
        raise SystemExit(f'Published artifact integrity mismatch: {name}')
print(f'Verified {len(manifest)} current published scientific artifacts; no training or inference.')
