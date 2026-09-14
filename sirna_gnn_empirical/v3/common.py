"""Version 3 review-directed empirical campaign utilities; historical inputs read only."""
from pathlib import Path
import json,hashlib,os,time
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).resolve().parent
RUN=ROOT/(HERE/'active_run.txt').read_text().strip()
OLD=ROOT/'runs/sirna_gnn_empirical/v1-20260913T185027Z'
PREV=ROOT/'runs/sirna_gnn_empirical/v2-20260913T220847Z'
SEEDS=[1103,2207,3301,4409,5519,6637,7753,8867,9973,11027]
DEVSEEDS=[701,907]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(x,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x),allow_nan=False)+'\n');tmp.replace(p)
def readlines(p):return [json.loads(s) for s in Path(p).read_text().splitlines() if s.strip()]
def rows_all():return sum([readlines(OLD/f'{s}_observations.jsonl') for s in ['train','validation','test_internal','test_APP']],[])
def weights(g):
 g=np.asarray(g);u,n=np.unique(g,return_counts=True);d=dict(zip(u,n));return np.array([len(g)/(len(u)*d[x]) for x in g],dtype=np.float64)
def export(df,name):
 p=RUN/name;p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,index=False);return p
