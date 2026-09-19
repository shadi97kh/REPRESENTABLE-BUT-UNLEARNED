"""Primary bibliography evidence acquisition, separate from scientific fitting."""
from common import *
import re,requests,concurrent.futures,html
P=ROOT/'papers/interaction_recoverability_iclr2027/v6';out=RUN/'sources/citations';out.mkdir(exist_ok=True)
bib='\n'.join((P/n).read_text() for n in ['references.bib','references_new.bib']);entries=re.split(r'(?=@(?:article|misc|inproceedings)\{)',bib);records=[]
alternatives={'gneiting2007':'https://sites.stat.washington.edu/people/raftery/Research/PDF/Gneiting2007jasa.pdf','ridge1970':'https://www.stat.cmu.edu/technometrics/70-79/VOL-12-01/v1201055.pdf','xu2019':'https://arxiv.org/abs/1810.00826','zaheer2017':'https://arxiv.org/abs/1703.06114','numpy2020':'https://numpy.org/citing-numpy/','matplotlib2007':'https://matplotlib.org/stable/project/citing.html','bramsen2009':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2685080/fullTextXML','davis2025':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12205987/fullTextXML','ensi2025':'https://pubmed.ncbi.nlm.nih.gov/40194620/','meg2026':'https://pubmed.ncbi.nlm.nih.gov/42228641/'}
def get(entry):
 m=re.match(r'@\w+\{([^,]+),',entry)
 if not m:return None
 key=m.group(1);du=re.search(r'doi=\{([^}]+)\}',entry);uu=re.search(r'url=\{([^}]+)\}',entry);url=alternatives.get(key,uu.group(1) if uu else 'https://doi.org/'+du.group(1));record=out/(key+'.json')
 if record.exists():return json.loads(record.read_text())
 t=time.perf_counter();rec=dict(key=key,primary_url=url,bib_entry=entry)
 try:
  z=requests.get(url,timeout=40);p=out/(key+('.pdf' if z.content.startswith(b'%PDF') else '.html'));p.write_bytes(z.content);rec.update(status=z.status_code,final_url=z.url,path=str(p.relative_to(RUN)),sha256=sha(p),bytes=len(z.content))
  if z.status_code==200 and not z.content.startswith(b'%PDF'):
   from bs4 import BeautifulSoup
   soup=BeautifulSoup(z.content,'html.parser');ab=soup.find(class_='abstract') or soup.find('abstract');text=soup.get_text(' ',strip=True);(out/(key+'.txt')).write_text(text);rec['abstract_excerpt']=ab.get_text(' ',strip=True) if ab else ''
 except Exception as e:rec['error']=repr(e)
 rec['wall_s']=time.perf_counter()-t;write(record,rec);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
 for r in pool.map(get,entries):
  if r:records.append(r);print(r['key'],r.get('status'),flush=True)
write(out/'acquisition_manifest.json',records)
