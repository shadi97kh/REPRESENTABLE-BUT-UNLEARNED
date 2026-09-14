import argparse,json,re
from pathlib import Path
from collections import Counter
import pandas as pd
from bs4 import BeautifulSoup
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);root=ap.parse_args().run
out=root/'source_inspection';out.mkdir(exist_ok=True)
cm=pd.read_csv(root/'raw/cmsirnadb_APP.tsv',sep='\t',keep_default_na=False)
print('CM',cm.shape)
for k in ['patent_ID','Cell_Type','Concentration','Time_of_administration']:
 print(k,cm[k].value_counts().to_dict())
print(cm.groupby(['patent_ID','Cell_Type','Concentration','Time_of_administration']).size().to_string())
for patent in cm.patent_ID.unique():
 soup=BeautifulSoup((root/'raw'/f'{patent}.html').read_text(),'lxml')
 desc=soup.find('section',itemprop='description') or soup
 text=desc.get_text(' ',strip=True)
 (out/f'{patent}.txt').write_text(text)
 tables=[]
 for i,t in enumerate(desc.find_all('table')):
  rows=[[c.get_text(' ',strip=True) for c in row.find_all(['td','th'],recursive=False)] for row in t.find_all('tr')]
  tables.append({'index':i,'id':t.get('id'),'rows':rows})
 (out/f'{patent}_tables.json').write_text(json.dumps(tables,indent=2))
 print('PATENT',patent,'text length',len(text),'tables',len(tables))
 for t in tables:print(t['index'],t['id'],len(t['rows']),str(t['rows'][:3])[:350])
 for term in ['transfection','Table 8','TABLE 8','Be(2)','BE(2)','Lipofectamine','LIPOFECTAMINE','Table 13','TABLE 13']:
  hits=[m.start() for m in re.finditer(re.escape(term),text)]
  if hits: print('TERM',term,[text[max(0,i-200):i+600] for i in hits[:2]])
ense=[]
for p in sorted((root/'raw/ensi/ENsiRNA-mod/dataset').glob('*.xlsx')):
 d=pd.read_excel(p).fillna('');print('ENSI',p.name,d.shape);ense.append(d)
 print('MODS',Counter(str(v).strip() for col in ['sense mod','anti mod'] for v in d[col]).most_common(12))
print('DAVIS HEAD',BeautifulSoup((root/'raw/davis2025.html').read_text(),'lxml').get_text(' ',strip=True)[:500])
soup=BeautifulSoup((root/'raw/davis2025.html').read_text(),'lxml')
print('DAVIS SUPPLEMENT',[(a.get_text(' ',strip=True),a.get('href')) for a in soup.find_all('a') if 'supp' in (a.get('href','') or '').lower()])
(root/'inspect_sources.outputs.json').write_text(json.dumps([str(p.relative_to(root)) for p in out.iterdir()]))
