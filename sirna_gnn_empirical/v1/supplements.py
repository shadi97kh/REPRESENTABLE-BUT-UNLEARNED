import argparse,json,hashlib,urllib.request,zipfile
from pathlib import Path
from bs4 import BeautifulSoup
import pandas as pd
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);root=ap.parse_args().run
soup=BeautifulSoup((root/'raw/davis2025.html').read_text(),'lxml')
urls=list(dict.fromkeys(a.get('href') for a in soup.find_all('a') if 'gkaf479_supplemental_files.zip' in (a.get('href','') or '')))
records=[];outputs=[]
for url in urls[::-1]:
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=35) as response:data=response.read(32*1024*1024)
  dest=root/'raw/davis2025_supplemental_files.zip';dest.write_bytes(data)
  with zipfile.ZipFile(dest) as z:
   assert z.testzip() is None
   for name in z.namelist():assert not name.startswith('/') and '..' not in Path(name).parts
   assert sum(i.file_size for i in z.infolist())<64*1024*1024
   z.extractall(root/'raw/davis2025_supplements')
   print('ZIP MEMBERS',z.namelist())
  records.append({'status':'acquired','url':url,'path':str(dest.relative_to(root)),'sha256':hashlib.sha256(data).hexdigest()});outputs.append(str(dest.relative_to(root)));break
 except Exception as e:records.append({'status':'unavailable','url':url,'error':repr(e)})
for p in (root/'raw/davis2025_supplements').rglob('*'):
 if p.is_file():outputs.append(str(p.relative_to(root)))
 if p.suffix.lower()=='.xlsx':
  print('FILE',p.name)
  book=pd.ExcelFile(p)
  for sheet in book.sheet_names:
   d=pd.read_excel(book,sheet_name=sheet,header=None);print('SHEET',sheet,d.shape);print(d.head(9).to_string(index=False,header=False)[:5000])
(root/'supplement_acquisition.json').write_text(json.dumps(records,indent=2)+'\n');outputs.append('supplement_acquisition.json')
(root/'supplements.outputs.json').write_text(json.dumps(outputs,indent=2))
print('STATUS',[(x['status'],x.get('error','')) for x in records])
