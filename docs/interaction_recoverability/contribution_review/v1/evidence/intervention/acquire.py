"""Pin public data and repository inventories, without executing author code."""
from __future__ import annotations
import argparse
import datetime
import json
from pathlib import Path
import requests
from .accounting import elapsed, snapshot, write_bytes, write_json

ROOT = Path('data/intervention_sources')
REPOSITORIES = {
    'agile': ('bowang-lab/AGILE', 'aee76246c90634500fd73579ac7ab2f94f91cc04'),
    'lantern': ('AsalMehradfar/LANTERN', '11240f29ef92323649ae60d177b21df77e2d428b'),
    'transma': ('wklix/TransMA', '11de7ec979fa90df91f1efd4334358288fcce00c'),
    'lipobart': ('Sanofi-Public/LipoBART', None),
    'rnagym': ('MarksLab-DasLab/RNAGym', None),
}


def fetch(url):
    r = requests.get(url, timeout=40)
    r.raise_for_status()
    if len(r.content) > 15_000_000:
        raise ValueError('source exceeds 15 MB bounded audit limit')
    return r


def run():
    start = snapshot()
    manifest = {'retrieved_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'source_evidence': 'GNN_RNA_Data_Research_Evidence.json',
                'files': [], 'repositories': {}, 'training_steps': 0}
    evidence = json.loads(Path(manifest['source_evidence']).read_text())
    for item in evidence['retrievals']:
        r = fetch(item['url'])
        file = write_bytes(ROOT / (item['name'] + '.csv'), r.content)
        assert file['sha256'] == item['sha256'], f'input hash mismatch: {item["name"]}'
        manifest['files'].append({**item, **file, 'hash_matches_prior_inspection': True})
    for name, (repo, commit) in REPOSITORIES.items():
        if commit is None:
            commit = fetch(f'https://api.github.com/repos/{repo}/commits?per_page=1').json()[0]['sha']
        tree = fetch(f'https://api.github.com/repos/{repo}/git/trees/{commit}?recursive=1').json()
        if tree.get('truncated'):
            raise ValueError('incomplete source inventory')
        tree_record = write_bytes(ROOT / name / 'tree.json', json.dumps(tree, indent=2).encode())
        manifest['repositories'][name] = dict(repo=repo, commit=commit, tree=tree_record,
                                              files=[], checkpoint_loaded=False)
        for entry in tree['tree']:
            p = entry['path']
            if entry['type'] != 'blob':
                continue
            keep = (p.lower().endswith(('.py', '.yaml', '.yml')) or
                    Path(p).name.lower() in ('readme.md', 'license', 'license.md', 'license.txt'))
            if name in ('rnagym','lipobart'):
                keep = keep or (p.lower().endswith(('.csv', '.tsv')) and entry.get('size', 0) < 1_000_000)
            if not keep or entry.get('size', 0) > 400_000 or '/.ipynb_checkpoints/' in p:
                continue
            url = f'https://raw.githubusercontent.com/{repo}/{commit}/{p}'
            try:
                r = fetch(url)
                record = write_bytes(ROOT / name / p, r.content)
                manifest['repositories'][name]['files'].append(dict(url=url, **record))
            except requests.RequestException as exc:
                manifest['repositories'][name]['files'].append(dict(url=url, error=str(exc)))
        print(name, commit, len(manifest['repositories'][name]['files']), flush=True)
    manifest['resources'] = elapsed(start)
    write_json(ROOT / 'manifest-v1.json', manifest)
    print(json.dumps(manifest['resources']), flush=True)


if __name__ == '__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    run()
