"""Administrative copying, JSON inspection and hashing only; imports no project code."""
from pathlib import Path
import hashlib
import json
import re
import shutil
import time
from urllib.parse import unquote

ROOT = Path('/home/shadi/iclr2027')
OUT = ROOT / 'docs/interaction_recoverability/contribution_review/v1'
started_wall, started_cpu = time.perf_counter(), time.process_time()

def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()

def links(p):
    text = p.read_text()
    text = re.sub(r'^```[^\n]*\n.*?^```\s*$', '', text, flags=re.M | re.S)
    return [unquote(m.group(1).split(' "')[0].strip('<>').split('#')[0])
            for m in re.finditer(r'!?\[[^\]\n]+\]\(([^\)\n]+)\)', text)]

seeds = set()
for d in ['warp_family/v1', 'warp_family_lower_bound/v2', 'external_review/v2',
          'practical_feasibility/v1', 'efficiency_audit/v1', 'learned_target_v1']:
    seeds.update((ROOT / 'docs/interaction_recoverability' / d).glob('*.md'))
for p in OUT.glob('*.md'):
    for link in links(p):
        if link.startswith('evidence/'):
            seeds.add(ROOT / link.removeprefix('evidence/'))
for directory in ['interaction_recoverability/learned_target', 'scmp', 'certmp', 'intervention']:
    seeds.update((ROOT / directory).rglob('*.py'))
run = ROOT / 'runs/interaction_recoverability/learned-target-v1-20260912-164000'
seeds.update(run.glob('predictions/*/predictions.json'))
seeds.update((run / 'predictions/d000').glob('*.json'))
seeds.update([run / 'evaluation/prediction_hashes.json', ROOT / 'Codex_Combined_Review_GNN_siRNA.md'])
seeds.update((ROOT / 'runs/revision').glob('p4_gpu-20260911-171200-784773-*.json'))

queue = list(seeds)
seen = set()
missing = []
while queue:
    p = queue.pop().resolve()
    if p in seen:
        continue
    if not p.is_relative_to(ROOT) or p.is_relative_to(OUT) or not p.is_file():
        missing.append(str(p))
        continue
    seen.add(p)
    if p.suffix == '.md':
        for target in links(p):
            if not target or re.match(r'[a-zA-Z][\w+.-]*:', target):
                continue
            queue.append(p.parent / target)
if missing:
    raise RuntimeError(f'Missing or invalid evidence dependency: {missing}')
records = []
for p in sorted(seen):
    rel = p.relative_to(ROOT)
    dest = OUT / 'evidence' / rel
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(p, dest)
    src_hash = digest(p)
    assert src_hash == digest(dest), str(rel)
    records.append({'original_path': str(rel), 'packet_path': 'evidence/' + str(rel),
                    'bytes': p.stat().st_size, 'sha256': src_hash})

cfg = json.loads((run / 'frozen_config.json').read_text())
commitments = []
for rel, expected in cfg['learner_code_sha256'].items():
    actual = digest(ROOT / rel)
    commitments.append({'path': rel, 'expected': expected, 'actual': actual, 'match': actual == expected})
for dataset, expected in json.loads((run / 'evaluation/prediction_hashes.json').read_text()).items():
    p = run / 'predictions' / dataset / 'predictions.json'
    actual = digest(p)
    commitments.append({'path': str(p.relative_to(ROOT)), 'expected': expected, 'actual': actual, 'match': actual == expected})

cycle = json.loads((ROOT / 'runs/revision/p4_gpu-20260911-171200-784773-cycle.json').read_text())
for item in cycle['artifacts'].values():
    p = ROOT / item['path']
    actual = digest(p)
    commitments.append({'path': item['path'], 'expected': item['sha256'], 'actual': actual,
                        'match': actual == item['sha256'], 'checkpoint_not_deserialized': p.suffix == '.pt'})
p = ROOT / cycle['config_path']
commitments.append({'path': cycle['config_path'], 'expected': cycle['config_sha256'],
                    'actual': digest(p), 'match': digest(p) == cycle['config_sha256']})
assert all(c['match'] for c in commitments), 'Historical commitment mismatch'

omitted = []
for name in ['siRNA_Design_Codex_Research_Brief.md', 'siRNA_Transfer_Source_Evidence.json',
             'data/davis2025/gkaf479_supplemental_files.zip']:
    omitted.append({'original_path': name, 'status': 'not found in inspected evidence',
                    'reason': 'biological evidence gap; no substitute fabricated'})
for p in sorted((ROOT / 'runs/revision').glob('p4_gpu-20260911-171200-784773-*.pt')):
    omitted.append({'original_path': str(p.relative_to(ROOT)), 'status': 'exists; omitted from static packet',
                    'reason': 'checkpoint inference not required or authorized', 'bytes': p.stat().st_size, 'sha256': digest(p)})
for rel in ['runs/interaction_recoverability/prep-20260912-082154/unknown-profiles-v1/trabs_v2.txt',
            'runs/interaction_recoverability/prep-20260912-082154/sources/dextri_v1.txt',
            'runs/interaction_recoverability/prep-20260912-082154/sources/strong_functionals.txt']:
    p = ROOT / rel
    omitted.append({'original_path': rel, 'status': 'primary text inspected; not redistributed',
                    'reason': 'versioned primary links in review', 'bytes': p.stat().st_size, 'sha256': digest(p)})
report = {'evidence': records, 'missing_active_evidence_links': [], 'omitted': omitted,
          'historical_commitments': commitments,
          'cost_scope': 'this helper: traversal, copying, hashing and JSON parsing; excludes prior inspection and authoring',
          'process_cpu_seconds_before_report_write': time.process_time() - started_cpu,
          'wall_seconds_before_report_write': time.perf_counter() - started_wall}
(OUT / 'administrative/evidence_inventory.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'evidence_files': len(records), 'evidence_bytes': sum(x['bytes'] for x in records),
                  'historical_hash_checks': len(commitments), 'all_match': all(c['match'] for c in commitments),
                  'process_cpu_seconds': time.process_time() - started_cpu,
                  'wall_seconds': time.perf_counter() - started_wall}))
