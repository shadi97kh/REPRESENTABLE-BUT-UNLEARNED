"""Download pinned public measurements and primary-source material, never execute upstream code."""
import argparse,concurrent.futures,hashlib,json,time,urllib.request
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);args=ap.parse_args();root=args.run
repo=Path(__file__).resolve().parents[2]; evidence=json.loads((repo/'siRNA_Transfer_Source_Evidence.json').read_text())
items=[{'name':'ensi/'+x['path'],'url':x['url'],'expected_sha256':x['sha256'],'required':True} for x in evidence['data_inspections']['ensi_inspection.json']]
cm=evidence['data_inspections']['cmsirnadb_app_inspection.json'];items.append({'name':'cmsirnadb_APP.tsv','url':cm['url'],'expected_sha256':cm['sha256'],'required':True})
commit=evidence['data_inspections']['tanwenchong_ENsiRNA.json']['sha']
for path in ['ENsiRNA-mod/README.md','ENsiRNA-mod/data/mod_utils.py','ENsiRNA-mod/data/dataset.py','ENsiRNA-mod/models/mod_model.py','ENsiRNA-mod/train.py']:
 items.append({'name':'ensi/'+path,'url':f'https://raw.githubusercontent.com/tanwenchong/ENsiRNA/{commit}/{path}','required':False})
for patent in ['WO2020132227A2','WO2022165172A1','US20220380773A1']:
 items.append({'name':patent+'.html','url':f'https://patents.google.com/patent/{patent}/en','required':False})
items.extend([
 {'name':'cmsirnadb_paper.html','url':'https://link.springer.com/article/10.1186/s12859-025-06359-y','required':False},
 {'name':'davis2025.html','url':'https://academic.oup.com/nar/article/53/12/gkaf479/8171869','required':False},
 {'name':'ensi_tree.json','url':f'https://api.github.com/repos/tanwenchong/ENsiRNA/git/trees/{commit}?recursive=1','required':False}])
def fetch(item):
 start=time.monotonic();p=root/'raw'/item['name'];p.parent.mkdir(parents=True,exist_ok=True)
 try:
  if p.exists():data=p.read_bytes();item['reused']=True
  else:
   req=urllib.request.Request(item['url'],headers={'User-Agent':'Mozilla/5.0 (public scientific source retrieval)'})
   with urllib.request.urlopen(req,timeout=35) as response:
    data=response.read(32*1024*1024+1);assert len(data)<=32*1024*1024,'source exceeds per-file cap'
   p.write_bytes(data)
  actual=hashlib.sha256(data).hexdigest()
  if item.get('expected_sha256'):assert actual==item['expected_sha256'],'pinned checksum mismatch'
  item.update(status='acquired',sha256=actual,bytes=len(data),path=str(p.relative_to(root)))
 except Exception as e:item.update(status='unavailable',error=repr(e))
 item['wall_s']=time.monotonic()-start;return item
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(fetch,items))
(root/'acquisition_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
(root/'acquire.outputs.json').write_text(json.dumps(['acquisition_manifest.json']+[x['path'] for x in records if x['status']=='acquired'],indent=2)+'\n')
for r in records:print(r['name'],r['status'],r.get('bytes'),r.get('error',''))
assert all(x['status']=='acquired' for x in records if x['required']),'Required pinned data unavailable; preserve partial downloads'
