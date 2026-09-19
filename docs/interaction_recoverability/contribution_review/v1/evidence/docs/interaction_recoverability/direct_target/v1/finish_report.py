"""Author the final note from completed records; overhead already conservatively debited."""
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
audit=json.loads((HERE/'resource_audit_job_31.json').read_text())
ledger=ROOT/'runs/interaction_recoverability/learned-target-v1-20260912-164000/ledger.jsonl'
assert hashlib.sha256(ledger.read_bytes()).hexdigest()==audit['ledger_sha256']
assert not audit['changed_preserved_files']
note=HERE/'commands_and_resources.md'
text=note.read_text()
assert '**Final accounting through current-ledger job 31.**' not in text
text+='''

**Final accounting through current-ledger job 31.** Values below come from [resource_audit_job_31.json](resource_audit_job_31.json) and the completed [snapshot 31](../../../../runs/interaction_recoverability/learned-target-v1-20260912-164000/resource_snapshot-31.json).

| Quantity | This continuation / final status |
|---|---:|
| Measured wrapper plus waited-child CPU | 7.365251 seconds = 0.122754 core-minutes |
| Aggregate admitted job wall time | 7.391338 seconds = 0.123189 minutes |
| Conservative task charge, including 150 seconds overhead | 157.716733 seconds = 2.628612 minutes |
| Additional cap, independently applied to CPU and job wall | 600 seconds each |
| Remaining under conservative task charge | 442.283267 seconds |
| Cumulative current-ledger conservative charge | 1128.109330 seconds |
| Remaining under original 7200-second caps | 6071.890670 seconds |
| Protected original reporting reserve | 600 seconds, preserved |
| GPU / numerical workers / numerical threads / affinity cores | 0 / 1 / 1 / 1 |
| Existing files verified unchanged | 336 |

Jobs 29 (identity checks), 30 (static code/source read), and 31 (static final check) all completed with status zero, no timeout and no observed affinity violation. No fit or new numerical experiment followed job 29. The final ledger SHA-256 is `7fdd296abb73fa66456019c43685def0a02cac09966394f7c10356dadba908b7`. The historical preparation ledger hash and all frozen learner/prediction/configuration commitments remain unchanged. The current ledger was intentionally appended under this continuation authorization; previous resource snapshots and reports remain historical records.

The 150-second debit is conservative unmeasured overhead, not 150 measured CPU seconds. It includes registration, independent static-agent reads, final verification and this note's authoring. Elapsed research/session time is neither CPU time nor aggregate compute-job wall time. A sandbox failure prevented an attempted file-tool wording update; the original text was retained until this final authoring step corrected the two phrases. This was not a numerical failure and did not change any historical file. [finish_report.py](finish_report.py) writes only new continuation notes and their hash manifest from the completed snapshot; it does not append the ledger or run numerical work.

The new [bundle hashes](bundle_hashes.json), [preservation check](final_checks.json), [exact admitted commands](commands.sh), and [static discovery commands](discovery_commands.sh) provide the audit trail. No automatic continuation or additional training authorization is requested.
'''
note.write_text(text)
identity=HERE/'target_identity.md'
s=identity.read_text().replace('The truncated target equals','The truncated target `I12-R_K` equals').replace('Dropping `(1-chi_B)` from the integrand incurs at most','Multiplying the direct functional by chi_B incurs error at most')
identity.write_text(s)
commands=HERE/'commands.sh'
with commands.open('a') as f:
    f.write('PYTHONDONTWRITEBYTECODE=1 python -B docs/interaction_recoverability/direct_target/v1/finish_report.py\n')
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
        for p in HERE.rglob('*') if p.is_file() and p.name not in ['task.lock','bundle_hashes.json']}
module=ROOT/'interaction_recoverability/direct_target_v1.py'
hashes[str(module.relative_to(ROOT))]=hashlib.sha256(module.read_bytes()).hexdigest()
with (HERE/'bundle_hashes.json').open('x') as f:json.dump(hashes,f,indent=2)
print(json.dumps({'report_finalized':True,'hashed_new_files':len(hashes),'ledger_unchanged_since_final_verification':True}))
