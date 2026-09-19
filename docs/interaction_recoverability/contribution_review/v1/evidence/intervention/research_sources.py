"""Bounded retrieval of primary papers for a local, hash-pinned audit."""
from __future__ import annotations
import concurrent.futures
import datetime
import json
from pathlib import Path
import re
import subprocess
import requests
from .accounting import snapshot, elapsed, write_bytes, write_json

ARXIV = {
    'pgexplainer': '2011.04573', 'fastshap': '2107.07436',
    'stochastic_amortization': '2401.15866', 'graphshapiq': '2501.16944',
    'spex': '2502.13870', 'proxyspex': '2505.17495',
    'selective': '2405.19562', 'rex': '2507.12900',
    'degenerate_gnn': '2601.20815', 'gflowexplainer': '2303.02448',
    'cnp': '1807.01613', 'cgnp': '1812.05212', 'mpnp': '2009.13895',
    'dad': '2103.02438', 'tndp': '2411.02064',
    'lantern_preprint': '2507.03209', 'transma_preprint': '2407.05736',
}
OTHER = {
    'icexplainer': 'https://ojs.aaai.org/index.php/AAAI/article/download/39108/43070',
    'agile_published': 'https://www.nature.com/articles/s41467-024-50619-z.pdf',
    'lantern_published': 'https://www.nature.com/articles/s44488-026-00007-x.pdf',
    'squid': 'https://www.nature.com/articles/s42256-024-00851-5.pdf',
    'transma_published': 'https://academic.oup.com/bib/article-pdf/26/3/bbaf307/63544989/bbaf307.pdf',
    'lipobart': 'https://academic.oup.com/bioinformatics/article/40/6/btae342/7684951',
    'recent_uc_explainer': 'https://www.nature.com/articles/s44387-026-00135-w.pdf',
    'recent_np_costs': 'https://raw.githubusercontent.com/mlresearch/v341/main/assets/young26a/young26a.pdf',
}


def retrieve(item):
    name, url = item
    rec = dict(name=name, url=url)
    try:
        r = requests.get(url, timeout=40)
        rec.update(status_code=r.status_code, final_url=r.url)
        r.raise_for_status()
        if len(r.content) > 20_000_000:
            raise ValueError('bounded source limit exceeded')
        if not r.content.startswith(b'%PDF'):
            rec['status'] = 'not_pdf_or_access_blocked; not parsed as paper'
            return rec
        file = Path('data/intervention_sources/papers') / (name + '.pdf')
        rec.update(write_bytes(file, r.content))
        textfile = file.with_suffix('.txt')
        if not textfile.exists():
            subprocess.run(['pdftotext', '-layout', str(file), str(textfile)], check=True,
                           timeout=20, capture_output=True)
        text = textfile.read_text()
        rec['arxiv_version_in_text'] = re.findall(r'arXiv:\s*(\d{4}\.\d{4,5}v\d+)', text[:10000])
        rec['code_links'] = sorted(set(re.findall(r'https?://github.com/[^\s<>]+', text)))
        rec['status'] = 'retrieved_full_text'
    except Exception as exc:
        rec.update(status='unavailable', error=str(exc))
    return rec


def main():
    start = snapshot()
    urls = {**{k: 'https://arxiv.org/pdf/' + v for k,v in ARXIV.items()}, **OTHER}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(retrieve, urls.items()))
    write_json('data/intervention_sources/papers/manifest-v1.json',
               dict(retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    sources=records, resources=elapsed(start)))
    for r in records:
        print(r['name'], r['status'], r.get('arxiv_version_in_text'), r.get('code_links'), flush=True)


if __name__ == '__main__':
    main()
