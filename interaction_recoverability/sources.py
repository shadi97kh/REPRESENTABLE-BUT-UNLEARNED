"""Retrieve a fixed primary-source list serially; no dependencies are installed."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import requests

SOURCES={
 'dextri_v1':'https://arxiv.org/pdf/2608.19849v1',
 'pdae_v3':'https://arxiv.org/pdf/2504.18522v3',
 'bilinear_colt2023':'https://proceedings.mlr.press/v195/simchowitz23a/simchowitz23a.pdf',
 'extrapolation_jmlr2026':'https://www.jmlr.org/papers/volume27/24-0838/24-0838.pdf',
 'engression_chen_v1':'https://arxiv.org/pdf/2607.27723v1',
 'engression_huang_v1':'https://arxiv.org/pdf/2606.01002v1',
 'strong_functionals':'https://arxiv.org/pdf/2208.08291',
 'donoho_liu_II_technical_report':'https://digicoll.lib.berkeley.edu/record/85827/files/120.pdf',
}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',required=True)
    out=Path(parser.parse_args().out);out.mkdir(parents=True,exist_ok=False)
    manifest=[]
    for name,url in SOURCES.items():
        row=dict(name=name,url=url)
        try:
            r=requests.get(url,timeout=12);r.raise_for_status()
            if not r.content.startswith(b'%PDF'):raise ValueError('response is not a PDF')
            p=out/(name+'.pdf');p.write_bytes(r.content)
            row.update(path=str(p),bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
            subprocess.run(['pdftotext','-layout',str(p),str(p.with_suffix('.txt'))],check=True,timeout=12)
            row['status']='retrieved and text extracted'
        except Exception as e:row.update(status='unresolved retrieval',error=str(e))
        manifest.append(row);print(name,row['status'],flush=True)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__':main()
