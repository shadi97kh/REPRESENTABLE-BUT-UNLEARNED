"""Source-held-out prediction of primary measured B3 interactions; no new fit."""
from common import *
from engine import *
from metrics import metrics
out=RUN/'b3_primary';primary=readlines(out/'observations.jsonl');rectangles=json.loads((out/'rectangles.json').read_text());rows,raw=load();n=len(rows);praw=dict(np.load(out/'graphs.npz'));both=rows+primary;data={k:np.concatenate([raw[k],praw[k]]) for k in raw};test=np.arange(n,n+len(primary));results=[];endpoints=[];checks=[]
def km(seq):return {seq.replace('T','U')[i:i+13] for i in range(len(seq)-12)}
pset=set().union(*(km(r[st]['sequence']) for r in primary for st in ['guide','passenger']))
for fp in sorted((RUN/'fits/activity/outer0/final').rglob('fit.json')):
 rec=json.loads(fp.read_text());tr=np.array(rec['signature']['train']);va=np.array(rec['signature']['validation']);method=rec['method'];seed=rec['seed'];assert all(rows[i]['source_family']!='19282453' for i in np.r_[tr,va]);assert all(not(pset&km(rows[i][st]['sequence'])) for i in np.r_[tr,va] for st in ['guide','passenger']);dd,sup=transform(data,both,tr);assert sup==json.loads((fp.parent/'preprocessing.json').read_text())
 for name,h in rec['hashes'].items():assert sha(fp.parent/name)==h
 if rec['status']!='completed':values=np.full(len(primary),rec['training_equal_group_mean'])
 elif method in ['chemistry_tree','guide_ridge','pairwise_ridge']:
  with (fp.parent/'model.pkl').open('rb') as f:st=pickle.load(f)
  xx=classical_features(dd,method)[:,st['columns']]
  if method!='chemistry_tree':xx=(xx-st['mean'])/st['sd']
  values=st['model'].predict(xx[test])
 else:
  device=setup();st=torch.load(fp.parent/'best.pt',map_location=device,weights_only=False);m=network(method,sup).to(device);m.load_state_dict(st['model']);values=predict(m,dd,test,device,st['target_mean'],st['target_sd']);del m
 ep=dict(zip([x['record_id'] for x in primary],values));checks.append(dict(fit=rec['name'],source_overlap=0,sequence13mer_overlap=0,trained_updates=0,checkpoint_hash=rec['hashes'].get('best.pt',rec['hashes'].get('model.pkl'))))
 for row,value in zip(primary,values):endpoints.append(dict(record_id=row['record_id'],method=method,seed=seed,prediction=float(value),activity=row['activity'],measurement_sd=row['measurement_sd'],fit=rec['name']))
 for z in rectangles:
  ys={k:float(ep[rid]) for k,rid in z['record_ids'].items()};results.append(dict(rectangle_id=z['rectangle_id'],AS=z['AS'],SS=z['SS'],method=method,seed=seed,prediction=ys['11']-ys['10']-ys['01']+ys['00'],observed_interaction=z['observed_interaction'],**{f'prediction_{k}':v for k,v in ys.items()},**{f'observed_{k}':v for k,v in z['observed'].items()}))
export(pd.DataFrame(endpoints),'b3_primary/endpoint_predictions.csv');rr=pd.DataFrame(results);en=rr.groupby(['rectangle_id','AS','SS','method'],as_index=False).agg(prediction=('prediction','mean'),observed_interaction=('observed_interaction','first'));en['seed']='ensemble';allp=pd.concat([rr,en],ignore_index=True);mm=[]
for (method,seed),z in allp.groupby(['method','seed']):mm.append(dict(method=method,seed=seed,**metrics(z.observed_interaction,z.prediction)))
y=np.array([r['observed_interaction'] for r in rectangles]);mm.append(dict(method='zero_interaction',seed='deterministic',**metrics(y,np.zeros(len(y)))));export(allp,'b3_primary/interaction_predictions.csv');export(pd.DataFrame(mm),'tables/b3_metrics.csv');write(out/'source_holdout_checks.json',checks);write(out/'evaluation_summary.json',dict(rectangles=len(rectangles),conditions=len(primary),source_families=1,sequence_backgrounds=1,fitted_objects_reused=len(checks),new_fits=0,new_updates=0,uncertainty='Descriptive only: all rectangles share one background and a common measured reference. Endpoint covariance unresolved.'));print(pd.DataFrame(mm)[pd.DataFrame(mm).seed.isin(['ensemble','deterministic'])].to_string(index=False))
