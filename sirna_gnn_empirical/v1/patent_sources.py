import argparse,hashlib,json,urllib.request,subprocess,concurrent.futures
from pathlib import Path
from bs4 import BeautifulSoup
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
records=[]
def fetch(patent):
 rec={'patent':patent}
 try:
  soup=BeautifulSoup((r/'raw'/f'{patent}.html').read_text(),'lxml');url=soup.find('meta',attrs={'name':'citation_pdf_url'})['content'];rec['url']=url
  with urllib.request.urlopen(url,timeout=45) as response:data=response.read(64*1024*1024+1)
  assert data.startswith(b'%PDF') and len(data)<=64*1024*1024
  pdf=r/'raw'/f'{patent}.pdf';pdf.write_bytes(data);text=r/'source_inspection'/f'{patent}_pdf.txt'
  subprocess.run(['pdftotext','-layout',str(pdf),str(text)],check=True)
  rec.update(status='acquired',pdf=str(pdf.relative_to(r)),text=str(text.relative_to(r)),sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
  lines=text.read_text().splitlines()
  print('PATENT',patent,'BYTES',len(data),'TEXT',text.stat().st_size)
  for i,line in enumerate(lines):
   if ('Table 4' in line or 'Table 7' in line or 'Table 14' in line or 'Table 16' in line or 'Table 17' in line) and len(line)<180:print(i,line)
 except Exception as e:rec.update(status='unavailable',error=repr(e))
 return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:records=list(pool.map(fetch,['WO2020132227A2','WO2022165172A1']))
(r/'patent_pdf_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
(r/'patent_sources.outputs.json').write_text(json.dumps(['patent_pdf_manifest.json']+[v[k] for v in records if v['status']=='acquired' for k in ['pdf','text']]))
print(json.dumps(records,indent=2))
