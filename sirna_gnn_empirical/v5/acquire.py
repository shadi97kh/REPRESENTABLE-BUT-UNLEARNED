"""Read-only public primary-text acquisition. No fitting or external contact."""
from pathlib import Path
import requests, hashlib, json, time, concurrent.futures, subprocess
ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/(Path(__file__).with_name('active_run.txt').read_text().strip())
DEST=RUN/'literature'
URLS={
 'harsanyinet':'https://proceedings.mlr.press/v202/chen23s/chen23s.pdf',
 'ncm':'https://arxiv.org/pdf/2210.00035v1',
 'consistency':'https://arxiv.org/pdf/2405.15673v3',
 'purification':'https://proceedings.mlr.press/v108/lengerich20a/lengerich20a.pdf',
 'purification_supp':'https://proceedings.mlr.press/v108/lengerich20a/lengerich20a-supp.pdf',
 'gnavar':'https://arxiv.org/pdf/2606.08390v1',
 'guidelines':'https://iclr.cc/Conferences/2027/AuthorGuidelines',
 'ai_policy':'https://iclr.cc/Conferences/2027/AIPolicyForAuthors',
}
def fetch(item):
 key,url=item;t=time.monotonic();ext='.html' if key in ('guidelines','ai_policy') else '.pdf';p=DEST/(key+ext)
 try:
  r=requests.get(url,timeout=90);r.raise_for_status();p.write_bytes(r.content)
  if ext=='.pdf':
   assert r.content.startswith(b'%PDF');subprocess.run(['pdftotext','-layout',str(p),str(p.with_suffix('.txt'))],check=True)
  return dict(key=key,url=url,status='acquired',bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),wall_s=time.monotonic()-t)
 except Exception as e:return dict(key=key,url=url,status='failed',error=repr(e),wall_s=time.monotonic()-t)
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:rows=list(ex.map(fetch,URLS.items()))
 (DEST/'acquisition.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
