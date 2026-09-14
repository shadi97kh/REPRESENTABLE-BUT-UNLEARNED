"""Group-separated measured bundle-effect supervision and non-graph controls."""
from common import *
from engine import *
rows,raw=load();protocol=json.loads((RUN/'protocol.json').read_text());pairs=readlines(OLD/'eligible_pairs.jsonl');lu={x['record_id']:i for i,x in enumerate(rows)};endpointids={p[k] for p in pairs for k in ['reference_id','changed_id']};predictions=[];endpoint_predictions=[];selections=[];members=[]
def idx(groups):return np.array([i for i,x in enumerate(rows) if x['record_id'] in endpointids and x['sequence_group'] in groups],int)
def eligible(ii):
 ids={rows[i]['record_id'] for i in ii};return [p for p in pairs if p['reference_id'] in ids]
def choose(dev,grid):
 scores={c['name']:float(np.mean([r['validation_score'] if r['status']=='completed' else float('inf') for r in dev if r['config']['name']==c['name']])) for c in grid};return min(grid,key=lambda c:(scores[c['name']],c['name'])),scores

def pairclassical(name,method,cfg,seed,tr,va,te):
 out=RUN/'fits'/name;out.mkdir(parents=True,exist_ok=True);sig=dict(code=sha(HERE/'pair.py'),engine=sha(HERE/'engine.py'),protocol=sha(RUN/'protocol.json'),inputs=sha(RUN/'preservation_before.json'),train=tr.tolist(),validation=va.tolist(),test=te.tolist(),method=method,config=cfg,seed=seed)
 if (out/'fit.json').exists():
  r=json.loads((out/'fit.json').read_text());assert r['signature']==sig
  for p,h in r['hashes'].items():assert sha(out/p)==h
  return r,dict(np.load(out/'predictions.npz'))
 t=time.perf_counter();cpu=time.process_time();data,sup=transform(raw,rows,tr);trainp=eligible(tr);valp=eligible(va);testp=eligible(te)
 def features(pp):
  a=np.array([lu[p['reference_id']] for p in pp],int);b=np.array([lu[p['changed_id']] for p in pp],int)
  return np.concatenate([(data['X'][b]-data['X'][a]).reshape(len(pp),64*93),data['X'][a,:32,:5].reshape(len(pp),160)],1)
 xx=features(trainp);columns=np.flatnonzero(np.ptp(xx,axis=0)>0);xx=xx[:,columns];mu=xx.mean(0);sd=xx.std(0);sd[sd<1e-6]=1;xx=(xx-mu)/sd;yy=np.array([p['observed_difference'] for p in trainp]);ww=weights([rows[lu[p['reference_id']]]['sequence_group'] for p in trainp]);m=Ridge(alpha=cfg['alpha'],solver='lsqr') if method=='pair_ridge' else ExtraTreesRegressor(n_estimators=300,min_samples_leaf=cfg['leaf'],max_features=.3,n_jobs=2,random_state=seed);m.fit(xx,yy,sample_weight=ww)
 pr={k:m.predict((features(pp)[:,columns]-mu)/sd) for k,pp in [('validation',valp),('test',testp)] if pp};vscore=float(np.average((pr['validation']-np.array([p['observed_difference'] for p in valp]))**2,weights=weights([rows[lu[p['reference_id']]]['sequence_group'] for p in valp])))
 with (out/'model.pkl').open('wb') as f:pickle.dump(dict(model=m,columns=columns,mean=mu,sd=sd,preprocessing=sup,training_pair_ids=[p['pair_id'] for p in trainp]),f)
 np.savez_compressed(out/'predictions.npz',**pr);write(out/'preprocessing.json',sup);r=dict(name=name,method=method,config=cfg,seed=seed,signature=sig,status='completed',failure=None,validation_score=vscore,training_n=len(trainp),validation_n=len(valp),test_n=len(testp),selected_epoch=0,selected_updates=0,executed_updates=0,epochs_executed=0,history=[],checkpoint='model.pkl',wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu,hashes={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='fit.json'});write(out/'fit.json',r);print(name,vscore,flush=True);return r,pr

def add(pop,method,seed,pp,z,fitname,status='completed'):
 for p,value in zip(pp,z):predictions.append(dict(population=pop,method=method,seed=seed,pair_id=p['pair_id'],reference_id=p['reference_id'],changed_id=p['changed_id'],background_id=p['background_id'],sequence_group=rows[lu[p['reference_id']]]['sequence_group'],observed_difference=p['observed_difference'],prediction=float(value),measurement_sd_reference=p['measurement_sd_reference'],measurement_sd_changed=p['measurement_sd_changed'],mean_difference_se=None,fit=fitname,status=status))
for pop in protocol['pair_populations']:
 test=idx(pop['test_groups']);final=pop['inner'][pop['final_validation_fold']];tr,va=idx(final['train_groups']),idx(final['validation_groups']);testpairs=eligible(test);trainpairs=eligible(tr)
 for role,ii in [('train',tr),('validation',va),('test',test)]:
  for i in ii:members.append(dict(population=pop['name'],role=role,record_id=rows[i]['record_id'],sequence_group=rows[i]['sequence_group']))
 dy=np.array([p['observed_difference'] for p in trainpairs]);ww=weights([rows[lu[p['reference_id']]]['sequence_group'] for p in trainpairs]);mean=float(np.average(dy,weights=ww));sign=np.where(abs(dy)<=.02,0,np.sign(dy));majority=min([-1,0,1],key=lambda k:(-ww[sign==k].sum(),k))
 for m,c in [('zero',0.),('training_mean',mean),('training_majority',majority)]:add(pop['name'],m,0,testpairs,np.repeat(c,len(testpairs)),'training-only control')
 for method in ['corrected_gnn','corrected_no_message']:
  dev=[]
  for j,inner in enumerate(pop['inner']):
   dtr,dva=idx(inner['train_groups']),idx(inner['validation_groups'])
   for cfg in protocol['neural_configurations']:
    for seed in DEVSEEDS:
     name=f'pair/{pop["name"]}/activity_dev/{method}/i{j}/{cfg["name"]}/s{seed}';r,p=fit(name,method,cfg,seed,rows,raw,dtr,dva,np.array([],int),groupkey='sequence_group',pairs=pairs,lam=0.,selection='activity');dev.append(r)
  cfg,scores=choose(dev,protocol['neural_configurations']);selections.append(dict(population=pop['name'],method=method,kind='activity_config',config=cfg,scores=scores));write(RUN/'pair_selection.json',selections)
  ld=[]
  for j,inner in enumerate(pop['inner']):
   dtr,dva=idx(inner['train_groups']),idx(inner['validation_groups'])
   for lam in protocol['pair_lambdas']:
    for seed in DEVSEEDS:
     name=f'pair/{pop["name"]}/pair_dev/{method}/i{j}/l{lam:g}/s{seed}';r,p=fit(name,method,cfg,seed,rows,raw,dtr,dva,np.array([],int),groupkey='sequence_group',pairs=pairs,lam=lam,selection='pair');ld.append(r)
  ls={lam:float(np.mean([r['validation_score'] if r['status']=='completed' else float('inf') for r in ld if r['lambda_pair']==lam])) for lam in protocol['pair_lambdas']};lam=min(ls,key=lambda k:(ls[k],k));selections.append(dict(population=pop['name'],method=method,kind='pair_lambda',lambda_pair=lam,scores=ls));write(RUN/'pair_selection.json',selections)
  tasks=[('activity_gnn',0.,'activity'),('pair_gnn',lam,'pair')] if method=='corrected_gnn' else [('pair_no_message',lam,'pair')]
  for label,ll,sel in tasks:
   for seed in SEEDS:
    name=f'pair/{pop["name"]}/final/{label}/s{seed}';r,p=fit(name,method,cfg,seed,rows,raw,tr,va,test,groupkey='sequence_group',pairs=pairs,lam=ll,selection=sel);ep=dict(zip([rows[i]['record_id'] for i in test],p['test']));add(pop['name'],label,seed,testpairs,[ep[p['changed_id']]-ep[p['reference_id']] for p in testpairs],name,r['status'])
    for i,z in zip(test,p['test']):endpoint_predictions.append(dict(population=pop['name'],method=label,seed=seed,record_id=rows[i]['record_id'],prediction=float(z),activity=rows[i]['activity'],sequence_group=rows[i]['sequence_group'],fit=name,status=r['status']))
  export(pd.DataFrame(predictions),'pair_predictions.csv');export(pd.DataFrame(endpoint_predictions),'pair_endpoint_predictions.csv')
 for method in ['pair_ridge','pair_tree']:
  grid=protocol['ridge_configurations'] if method=='pair_ridge' else protocol['tree_configurations'];dev=[]
  for j,inner in enumerate(pop['inner']):
   dtr,dva=idx(inner['train_groups']),idx(inner['validation_groups'])
   for cfg in grid:
    for seed in [0] if method=='pair_ridge' else DEVSEEDS:
     name=f'pair/{pop["name"]}/ridge_dev/i{j}/{cfg["name"]}' if method=='pair_ridge' else f'pair/{pop["name"]}/tree_dev/i{j}/{cfg["name"]}/s{seed}';r,p=pairclassical(name,method,cfg,seed,dtr,dva,np.array([],int));dev.append(r)
  cfg,scores=choose(dev,grid);selections.append(dict(population=pop['name'],method=method,config=cfg,scores=scores));write(RUN/'pair_selection.json',selections)
  for seed in [0] if method=='pair_ridge' else SEEDS:
   name=f'pair/{pop["name"]}/final/pair_ridge' if method=='pair_ridge' else f'pair/{pop["name"]}/final/pair_tree/s{seed}';r,p=pairclassical(name,method,cfg,seed,tr,va,test);add(pop['name'],method,seed,testpairs,p['test'],name)
 export(pd.DataFrame(predictions),'pair_predictions.csv');export(pd.DataFrame(endpoint_predictions),'pair_endpoint_predictions.csv');export(pd.DataFrame(members),'pair_membership.csv')
