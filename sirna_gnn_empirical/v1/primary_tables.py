import argparse,json,urllib.request,hashlib,concurrent.futures
from pathlib import Path
from bs4 import BeautifulSoup
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
records=[];outputs=[]
for patent in ['US20220002724A1','US20240117349A1']:
 url=f'https://patents.google.com/patent/{patent}/en';rec={'patent':patent,'url':url,'role':'Related national publication used to recover machine-readable tables, never independent measurements.'}
 try:
  with urllib.request.urlopen(url,timeout=40) as response:data=response.read(32*1024*1024)
  path=r/'raw'/f'{patent}.html';path.write_bytes(data);outputs.append(str(path.relative_to(r)))
  soup=BeautifulSoup(data,'lxml');desc=soup.find('section',itemprop='description') or soup;text=desc.get_text(' ',strip=True)
  tp=r/'source_inspection'/f'{patent}.txt';tp.write_text(text);outputs.append(str(tp.relative_to(r)))
  tables=[]
  for i,t in enumerate(desc.find_all('table')):
   rows=[[c.get_text(' ',strip=True) for c in row.find_all(['td','th'],recursive=False)] for row in t.find_all('tr')]
   tables.append({'index':i,'id':t.get('id'),'rows':rows})
  jp=r/'source_inspection'/f'{patent}_tables.json';jp.write_text(json.dumps(tables,indent=2));outputs.append(str(jp.relative_to(r)))
  rec.update(status='acquired',sha256=hashlib.sha256(data).hexdigest(),tables=len(tables),text_bytes=len(text))
  print('PATENT',patent,rec)
  for t in tables:print(t['index'],t['id'],len(t['rows']),str(t['rows'][:3])[:220])
 except Exception as e:rec.update(status='unavailable',error=repr(e))
 records.append(rec)
(r/'primary_table_manifest.json').write_text(json.dumps(records,indent=2));outputs.append('primary_table_manifest.json')
(r/'primary_tables.outputs.json').write_text(json.dumps(outputs))
