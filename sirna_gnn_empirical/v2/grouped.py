from core import *
rows,raw=load();protocol=json.loads((RUN/'protocol.json').read_text());results=[];predictions=[];selection=[]
for fold in protocol['grouped_folds']:
 f=fold['fold'];idx=lambda role:np.array([i for i,x in enumerate(rows) if x['study_group'] in fold[role+'_groups']],int);tr,va,te=[idx(s) for s in ['train','validation','test']]
 for kind in ['gnn','chemistry_tree','token_cnn','nongraph']:
  dev=[]
  for cfg in TREE_CONFIGS if kind=='chemistry_tree' else CONFIGS:
   rec,p=fit(f'grouped/fold{f}/development',kind,True,True,cfg,701,rows,raw,tr,va,np.array([],int));results.append(rec);dev.append(rec)
  cfg=min(dev,key=lambda r:(r['validation_score'],r['config']['name']))['config'];selection.append(dict(fold=f,kind=kind,config=cfg))
  for seed in SEEDS:
   rec,p=fit(f'grouped/fold{f}/final',kind,True,True,cfg,seed,rows,raw,tr,va,te);results.append(rec)
   for i,z in zip(te,p['test']):predictions.append(dict(record_id=rows[i]['record_id'],outer_fold=f,model=kind,seed=seed,prediction=float(z)))
write(RUN/'grouped_fits.json',results);write(RUN/'grouped_selection.json',selection);pd.DataFrame(predictions).to_csv(RUN/'grouped_predictions.csv',index=False)
# Freeze final all-ENsiRNA models from inner development evidence only, before APP scoring.
final=[];finalpred=[];choices=[];tr=np.array([i for i,x in enumerate(rows) if x['dataset']=='ENsiRNA']);te=np.array([i for i,x in enumerate(rows) if x['split']=='test_APP']);va=np.array([],int)
for kind in ['gnn','chemistry_tree','token_cnn','nongraph']:
 dev=[f for f in results if f['kind']==kind and '/development' in f['tag']];names=sorted({f['config']['name'] for f in dev});name=min(names,key=lambda n:(np.mean([f['validation_score'] for f in dev if f['config']['name']==n]),n));chosen=[f for f in dev if f['config']['name']==name];cfg=chosen[0]['config'];epochs=max(1,int(np.floor(np.median([f['selected_epoch'] for f in chosen])+.5)));choices.append(dict(kind=kind,config=cfg,fixed_epochs=epochs,inner_scores=[f['validation_score'] for f in chosen]));write(RUN/'deployment_selection.json',choices)
 for seed in SEEDS:
  rec,p=fit('deployment',kind,True,True,cfg,seed,rows,raw,tr,va,te,fixed_epochs=epochs);final.append(rec)
  for i,z in zip(te,p['test']):finalpred.append(dict(record_id=rows[i]['record_id'],split=rows[i]['split'],model=kind,arm='all_ENsi_inner_selected',seed=seed,prediction=float(z)))
write(RUN/'deployment_fits.json',final);pd.DataFrame(finalpred).to_csv(RUN/'deployment_predictions.csv',index=False)
