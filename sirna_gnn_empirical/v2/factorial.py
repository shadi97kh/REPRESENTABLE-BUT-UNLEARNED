from core import *
rows,raw=load();protocol=json.loads((RUN/'protocol.json').read_text());tr=np.array([i for i,x in enumerate(rows) if x['split']=='train']);va=np.array([i for i,x in enumerate(rows) if x['split']=='validation']);te=np.array([i for i,x in enumerate(rows) if x['split'].startswith('test')]);results=[];predictions=[];selection=[]
for kind,R,C in [('gnn',False,False),('gnn',True,False),('gnn',False,True),('gnn',True,True),('chemistry_tree',True,True),('token_cnn',True,True),('nongraph',True,True),('gnn_no_chemistry',True,True),('context_ridge',True,True)]:
 arm=f'R{int(R)}C{int(C)}';configs=TREE_CONFIGS if kind=='chemistry_tree' else [dict(name='alpha1')] if kind=='context_ridge' else CONFIGS;dev=[]
 for cfg in configs:
  f,p=fit(f'factorial/development/{arm}',kind,R,C,cfg,701,rows,raw,tr,va,np.array([],dtype=int));results.append(f);dev.append(f)
 selected=min(dev,key=lambda f:(f['validation_score'],f['config']['name']))['config'];selection.append(dict(model=kind,arm=arm,config=selected))
 for seed in ([1103] if kind=='context_ridge' else SEEDS):
  f,p=fit(f'factorial/final/{arm}',kind,R,C,selected,seed,rows,raw,tr,va,te);results.append(f)
  for i,z in zip(te,p['test']):predictions.append(dict(record_id=rows[i]['record_id'],split=rows[i]['split'],model=kind,arm=arm,seed=seed,prediction=float(z)))
# Both distinct constants are trained on the original training labels only.
y=np.array([rows[i]['activity'] for i in tr]);ww=weights([rows[i]['study_group'] for i in tr])
for name,value in [('mean_unweighted',y.mean()),('mean_equal_study',np.average(y,weights=ww))]:
 for i in te:predictions.append(dict(record_id=rows[i]['record_id'],split=rows[i]['split'],model=name,arm='R1C1',seed=0,prediction=float(value)))
write(RUN/'factorial_selection.json',selection);write(RUN/'factorial_fits.json',results);pd.DataFrame(predictions).to_csv(RUN/'factorial_predictions.csv',index=False)
