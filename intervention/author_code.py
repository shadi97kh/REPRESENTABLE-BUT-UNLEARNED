"""Pin small, relevant author-code excerpts for inspection, never execution."""
from pathlib import Path
import concurrent.futures
import requests
from .accounting import snapshot,elapsed,write_bytes,write_json

REPOS={
    'pgexplainer':('flyingdoog/PGExplainer',['Explainer','PGExplainer']),
    'fastshap':('iancovert/fastshap',['fastshap.py','surrogate.py']),
    'stochastic':('chanwkimlab/amortized-attribution',['amortized','train.py']),
    'graphshapiq':('FFmgll/GraphSHAP-IQ',['graphshapiq.py','graph.py','graph_shap']),
    'spex':('basics-lab/spectral-explain',['spex.py','lasso.py','sparse','explainer.py']),
    'selective':('LucasMonteiroPaes/selective-explanations',['select','uncertainty','initial']),
    'cnp':('google-deepmind/neural-processes',['conditional_neural_process.ipynb']),
    'dad':('ae-foster/dad',['design','networks','layers']),
    'tndp':('huangdaolang/amortized-decision-aware-bed',['model','tndp','encoder']),
    'squid':('evanseitz/squid-nn',['surrogate.py','mutagenizer.py','impress.py']),
    'degenerate_gnn':('steveazzolin/gnn_deg_expl',['extension','sufficiency','faithfulness','metric']),
}


def retrieve(item):
    name,(repo,terms)=item
    out=dict(repo=repo,files=[])
    try:
        r=requests.get(f'https://api.github.com/repos/{repo}/commits?per_page=1',timeout=30)
        r.raise_for_status();commit=r.json()[0]['sha']
        out['commit']=commit;out['commit_date']=r.json()[0]['commit']['committer']['date']
        r=requests.get(f'https://api.github.com/repos/{repo}/git/trees/{commit}?recursive=1',timeout=30)
        r.raise_for_status();tree=r.json()
        out['tree']=write_bytes(f'data/intervention_sources/author_code/{name}/tree.json',r.content)
        entries=[x for x in tree['tree'] if x['type']=='blob' and x.get('size',0)<300000 and
                 (x['path'].lower().endswith(('.py','.ipynb')) or x['path'].lower()=='readme.md')]
        selected=[x for x in entries if any(t.lower() in x['path'].lower() for t in terms) or x['path'].lower()=='readme.md']
        if len(selected)<2: selected=entries[:4]
        for x in selected[:8]:
            url=f'https://raw.githubusercontent.com/{repo}/{commit}/{x["path"]}'
            rr=requests.get(url,timeout=30);rr.raise_for_status()
            out['files'].append(dict(url=url,**write_bytes(f'data/intervention_sources/author_code/{name}/{x["path"]}',rr.content)))
        out['status']='pinned_not_executed'
    except Exception as e: out['status']='unavailable';out['error']=str(e)
    return name,out


if __name__=='__main__':
    start=snapshot()
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records=dict(pool.map(retrieve,REPOS.items()))
    write_json('data/intervention_sources/author_code/manifest-v1.json',dict(repositories=records,resources=elapsed(start)))
    for n,r in records.items(): print(n,r.get('commit'),r['status'],[x['path'] for x in r['files']],flush=True)
