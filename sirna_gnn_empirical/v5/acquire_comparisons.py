from acquire import *
URLS.clear()
URLS.update({
 'optimal_recovery':'https://arxiv.org/pdf/2006.03706',
 'multiplicity':'https://proceedings.mlr.press/v119/marx20a/marx20a.pdf',
 'selective':'https://proceedings.mlr.press/v97/geifman19a/geifman19a.pdf',
 'gin':'https://arxiv.org/pdf/1810.00826',
 'gnnexplainer':'https://ai.stanford.edu/~jure/pubs/gnnexplainer-neurips19.pdf',
})
rows=[fetch(x) for x in URLS.items()]
(DEST/'comparison_acquisition.json').write_text(json.dumps(rows,indent=2))
repos=['csluchen/harsanyinet','CausalAILab/NCMCounterfactuals','Jiyuan-Tan/NeuralPartialID','microsoft/interpret']
records=[]
for repo in repos:
 try:
  meta=requests.get('https://api.github.com/repos/'+repo,timeout=30).json();branch=meta['default_branch']
  commit=requests.get(f'https://api.github.com/repos/{repo}/commits/{branch}',timeout=30).json()['sha']
  tree=requests.get(f'https://api.github.com/repos/{repo}/git/trees/{commit}?recursive=1',timeout=30).json()
  d=DEST/'official_code'/repo.replace('/','__');d.mkdir(parents=True,exist_ok=True)
  (d/'tree.json').write_text(json.dumps(tree));selected=[]
  for entry in tree['tree']:
   p=entry['path'];low=p.lower()
   if entry['type']!='blob' or entry.get('size',0)>100000:continue
   if not (low.endswith(('.py','.cpp','.h','.md'))):continue
   if repo=='microsoft/interpret' and not ('purif' in low or low=='readme.md'):continue
   if repo!='microsoft/interpret' and not (low.endswith('.py') or low=='readme.md'):continue
   r=requests.get(f'https://raw.githubusercontent.com/{repo}/{commit}/{p}',timeout=30);r.raise_for_status();dest=d/p;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(r.content)
   selected.append(dict(path=p,sha256=hashlib.sha256(r.content).hexdigest()))
  records.append(dict(repo=repo,commit=commit,files=selected,status='acquired source; not executed'))
 except Exception as e:records.append(dict(repo=repo,error=repr(e),status='incomplete'))
(DEST/'official_code_manifest.json').write_text(json.dumps(records,indent=2));print([(r.get('repo'),r.get('commit'),len(r.get('files',[]))) for r in records])
