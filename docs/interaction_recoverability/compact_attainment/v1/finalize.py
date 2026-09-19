"""Static bundle validation and preservation check; no numerical experiment."""
import ast
import hashlib
import json
from pathlib import Path
import re

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]


def main():
    out=HERE/'final_checks.json'
    if out.exists():
        raise FileExistsError('Keep the previous static validation record')
    preserved=json.loads((HERE/'preservation.json').read_text())
    changed=[rel for rel,sha in preserved.items()
             if not (ROOT/rel).is_file() or hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()!=sha]
    syntax=[]
    for p in sorted(HERE.glob('*.py')):
        ast.parse(p.read_text(),filename=str(p)); syntax.append(p.name)
    tree=ast.parse((HERE/'estimator.py').read_text())
    sample_entry=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='estimate_from_sources')
    assert isinstance(sample_entry.body[0],ast.Raise)
    required=['construction.md','theorem.md','algorithm.md','decision.md','quadrature_proof.md',
              'review_addendum.md','novelty.md','static_review.md','checks.json','practical_counts.json']
    assert all((HERE/p).is_file() and (HERE/p).stat().st_size>0 for p in required)
    missing=[]
    for p in HERE.glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('http:','https:','#')):
                continue
            link=link.split('#')[0]
            if link=='commands_and_resources.md':
                continue  # Authored from the completed final job and gate snapshot.
            if not (p.parent/link).exists():
                missing.append({'document':p.name,'link':link})
    assert not changed,changed
    assert not missing,missing
    assert json.loads((HERE/'checks.json').read_text())['status']=='passed'
    assert json.loads((HERE/'practical_counts.json').read_text())['status']=='passed'
    result={'status':'passed','preserved_file_count':len(preserved),'changed_preserved_files':changed,
            'parsed_python_files':syntax,'unconditional_sample_entry_refusal':True,
            'missing_document_links':missing,'scope':'Static verification only; no mathematical experiment rerun.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    bundle={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir()
            if p.is_file() and p.name not in ('task.lock','bundle_hashes.json')}
    (HERE/'bundle_hashes.json').write_text(json.dumps({'scope':'Completed mathematical/code/check bundle before final accounting note.',
            'sha256':bundle},indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
