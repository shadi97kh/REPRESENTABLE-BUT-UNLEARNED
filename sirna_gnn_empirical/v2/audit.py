from core import *
import subprocess,platform,sys
rows,data=load();device=setup();train=np.array([i for i,x in enumerate(rows) if x['split']=='train']);valid=np.array([i for i,x in enumerate(rows) if x['split']=='validation']);test=np.array([i for i,x in enumerate(rows) if x['split'].startswith('test')]);manifest=json.loads((OLD/'graph_manifest.json').read_text())
zero=(data['X'][train][data['mask'][train].astype(bool)]==0).all(0);cz=(data['C'][train]==0).all(0)
rel=[];features=[];context=[]
for s in SPLITS:
 ii=np.array([i for i,x in enumerate(rows) if x['split']==s]);X=data['X'][ii][data['mask'][ii].astype(bool)]
 for j,name in enumerate(manifest['edge_relations']):rel.append(dict(split=s,relation=j,name=name,edges=int(data['A'][ii,j].sum())))
 for j in range(93):features.append(dict(split=s,column=j,valid_nonzero_nodes=int(np.count_nonzero(X[:,j])),min=float(X[:,j].min()),max=float(X[:,j].max()),identically_zero_training=bool(zero[j])))
 for j,name in enumerate(manifest['context_features']):context.append(dict(split=s,column=j,name=name,min=float(data['C'][ii,j].min()),max=float(data['C'][ii,j].max()),identically_zero_training=bool(cz[j])))
pd.DataFrame(rel).to_csv(RUN/'relation_support.csv',index=False);pd.DataFrame(features).to_csv(RUN/'feature_support.csv',index=False);pd.DataFrame(context).to_csv(RUN/'context_support.csv',index=False)
assert np.all(data['A'][train][:,[0,1,3,4]]==0) and np.all(data['A'][valid][:,[0,1,3,4]]==0)
assert sum(z['edges'] for z in rel if z['split']=='train' and z['relation'] in [2,5])==215050
chem=[]
for name in ['2-O-hexadecyl','Glycol nucleic acid','5-vinylphosphonate']:
 j=39+manifest['modification_vocabulary'].index(name)
 for s in SPLITS:
  ii=[i for i,x in enumerate(rows) if x['split']==s];chem.append(dict(name=name,split=s,rows=int(np.any(data['X'][ii,:,j],axis=1).sum())))
pd.DataFrame(chem).to_csv(RUN/'chemistry_support.csv',index=False)
fits=json.loads((OLD/'final_fit_results.json').read_text());fits=fits if isinstance(fits,list) else fits['fits']
records=[];gradients=[];replay=[];oldpred=pd.read_csv(OLD/'heldout_predictions_by_seed.csv')
for kind in ['gnn','nongraph','token_cnn']:
 for f in fits:
  if f['kind']!=kind:continue
  st=torch.load(OLD/f['checkpoint'],map_location=device,weights_only=False);m=model(kind,False,None).to(device);m.load_state_dict(st['model']);m.eval()
  # Pure data loss; no AdamW decay, optimizer step or modification of archived state.
  if kind=='gnn':
   m.zero_grad(set_to_none=True);target=np.array([rows[i]['activity'] for i in train]);w=weights([rows[i]['study_group'] for i in train]);
   for j in range(0,len(train),128):
    ii=train[j:j+128];dd={k:torch.as_tensor(v[ii],device=device) for k,v in data.items()};p=m(**dd)*st['target_std']+st['target_mean'];loss=(torch.tensor(w[j:j+128],device=device)*(p-torch.tensor(target[j:j+128],device=device))**2).sum()/len(train);loss.backward()
   for layer,ww in enumerate(m.message_weights):
    for k in range(8):gradients.append(dict(seed=f['seed'],type='relation',layer=layer,column=k,l2=float(ww.grad[k].norm())))
   for k in range(93):gradients.append(dict(seed=f['seed'],type='node',layer=0,column=k,l2=float(m.embed.weight.grad[:,k].norm())))
   for k in range(8):gradients.append(dict(seed=f['seed'],type='context',layer=0,column=k,l2=float(m.readout[0].weight.grad[:,-8+k].norm())))
  for fallback in ['original','backbone','zero_support','combined']:
   dd=dict(data)
   if fallback in ['backbone','combined']:
    dd['A']=data['A'].copy();before=dd['A'].sum(1);dd['A'][:,2]+=dd['A'][:,0]+dd['A'][:,1];dd['A'][:,5]+=dd['A'][:,3]+dd['A'][:,4];dd['A'][:,[0,1,3,4]]=0;assert np.array_equal(before,dd['A'].sum(1))
   if fallback in ['zero_support','combined']:
    dd['X']=data['X'].copy();dd['X'][:,:,zero]=0;dd['C']=data['C'].copy();dd['C'][:,cz]=0
    assert np.array_equal(dd['X'][train],data['X'][train]) and np.array_equal(dd['C'][train],data['C'][train])
   pp=pred(m,dd,test,device,st['target_mean'],st['target_std'])
   for i,p in zip(test,pp):records.append(dict(record_id=rows[i]['record_id'],split=rows[i]['split'],model=kind,arm=fallback,seed=f['seed'],prediction=float(p)))
   if fallback=='original':
    actual=dict(zip([rows[i]['record_id'] for i in test],pp));ref=oldpred[(oldpred.model==kind)&(oldpred.seed==f['seed'])];err=max(abs(actual[x.record_id]-x.prediction) for x in ref.itertuples());assert err<1e-6;replay.append(dict(model=kind,seed=f['seed'],maximum_replay_error=err))
 pd.DataFrame(records).to_csv(RUN/'frozen_diagnostic_predictions.csv',index=False)
grad=pd.DataFrame(gradients);assert grad[(grad.type=='relation')&grad.column.isin([0,1,3,4])].l2.max()==0
assert grad[(grad.type=='node')&grad.column.isin(np.flatnonzero(zero))].l2.max()==0
assert grad[(grad.type=='context')&grad.column.isin(np.flatnonzero(cz))].l2.max()==0
grad.to_csv(RUN/'data_gradients.csv',index=False);write(RUN/'frozen_replay.json',replay)
effects=pd.read_csv(OLD/'chemical_effect_predictions.csv');g=effects[effects.model=='gnn'];obs=np.sign(g.observed_difference)* (abs(g.observed_difference)>.02);pp=np.sign(g.predicted_difference)*(abs(g.predicted_difference)>.02);active=obs!=0;calls=active&(pp!=0)
check=dict(observed_negative=int((obs<0).sum()),observed_ties=int((obs==0).sum()),observed_positive=int((obs>0).sum()),calls=int(calls.sum()),correct_calls=int((calls&(obs==pp)).sum()),predicted_min=float(g.predicted_difference.min()),predicted_max=float(g.predicted_difference.max()),mean_observed=float(g.observed_difference.mean()),mean_predicted=float(g.predicted_difference.mean()),unsupported_relation_scalars=8192)
assert [check[k] for k in ['observed_negative','observed_ties','observed_positive','calls','correct_calls']]==[111,19,26,94,18]
write(RUN/'support_checks.json',check)
write(RUN/'environment.json',dict(python=sys.version,platform=platform.platform(),torch=torch.__version__,numpy=np.__version__,pandas=pd.__version__,cuda=torch.version.cuda,device=device,affinity=sorted(os.sched_getaffinity(0)),threads=torch.get_num_threads(),gpu_memory_fraction=.30,packages=subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True).splitlines(),nvidia_smi=subprocess.check_output(['nvidia-smi'],text=True),cgroup=subprocess.check_output(['cat','/proc/self/cgroup'],text=True)))
print(check)
