"""Static preservation and report validation, with saved-result summaries only."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def main():
    path=HERE/'final_checks.json'
    if path.exists():raise FileExistsError('Preserve earlier validation')
    prior=json.loads((HERE/'preservation.json').read_text())
    changed=[r for r,h in prior.items() if not (ROOT/r).is_file() or hashlib.sha256((ROOT/r).read_bytes()).hexdigest()!=h]
    assert not changed,changed
    parsed=[]
    for p in sorted(HERE.glob('*.py')):
        ast.parse(p.read_text());parsed.append(p.name)
    required=['variance_analysis.md','numerical_stability.md','feasibility_decision.md','arithmetic_addendum.md']
    assert all((HERE/p).is_file() and (HERE/p).stat().st_size for p in required)
    missing=[]
    for p in HERE.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('https:','http:','#')) or link=='commands_and_resources.md':continue
            if not (p.parent/link.split('#')[0]).exists():missing.append([p.name,link])
    assert not missing,missing
    cases=['variance_base','variance_depth','variance_integration','variance_tail','variance_deep','variance_final','variance_accumulator']
    records={c:json.loads((HERE/(c+'.json')).read_text()) for c in cases}
    assert all(r['status']=='passed' for r in records.values())
    codehash=hashlib.sha256((HERE/'influence_variance.py').read_bytes()).hexdigest()
    assert all(r['code_sha256']==codehash for r in records.values())
    assert json.loads((HERE/'error_and_precision.json').read_text())['status']=='passed'
    summary={'scope':'Differences of saved diagnostic outputs, not error certificates.',
        'integration_sd_difference':abs(records['variance_deep']['asymptotic_sd_constant']-records['variance_final']['asymptotic_sd_constant']),
        'accumulator_variance_difference':abs(records['variance_accumulator']['complete_variance']-records['variance_final']['complete_variance']),
        'cases':[{k:r[k] for k in ('case','complete_variance','asymptotic_sd_constant','child_cpu_s','child_wall_s','child_peak_rss_kib')} for r in records.values()]}
    (HERE/'diagnostic_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    result={'status':'passed','preserved_file_count':len(prior),'changed_preserved_files':changed,
      'parsed_python':parsed,'missing_links':missing,'raw_result_code_hashes_match':True,
      'numerical_cases_retained':cases,'scope':'Static validation only. No influence computation rerun.'}
    path.write_text(json.dumps(result,indent=2)+'\n')
    manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir() if p.is_file() and p.name not in ('task.lock','bundle_hashes.json')}
    (HERE/'bundle_hashes.json').write_text(json.dumps({'scope':'Completed analysis/code/raw bundle before final accounting files.','sha256':manifest},indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
