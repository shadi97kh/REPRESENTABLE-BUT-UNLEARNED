from core import *
rows,raw=load();protocol=json.loads((RUN/'protocol.json').read_text());pairs=readlines(OLD/'eligible_pairs.jsonl');lookup={x['record_id']:i for i,x in enumerate(rows)};endpointids={p[k] for p in pairs for k in ['reference_id','changed_id']};results=[];prediction=[];membership=[];selection=[]
# Admission depends on full verified assay context, not shared source or guide alone.
outside=[x for x in rows if x['split']!='test_APP'];eligible=[x for x in outside if x['assay_id'] and x['dose_unit']=='nM' and x['time_h'] is not None and x['cell'] and x['delivery'] and x['source_table_verified']]
write(RUN/'b2_admission.json',dict(outside_APP_rows=len(outside),outside_APP_rows_with_complete_comparability=len(eligible),admitted_independent_training_pairs=0,status='No outside-APP training pair admitted; ENsiRNA assay identity, concentration units, time and delivery unresolved. Published shared-source rows cannot stand in for primary matched evidence.',pair_count=len(pairs),endpoint_count=len(endpointids),sequence_components=12,backgrounds=len({p['background_id'] for p in pairs}),shared_control='AD1955 is a named assay normalizer; exact shared-well identity/covariance not reported. Fold grouping covers known sequence dependence, not certified independent measurement noise.',B3_complete_rectangles=0))
assert not eligible
for fold in protocol['b2_folds']:
 f=fold['fold'];idx=lambda role:np.array([i for i,x in enumerate(rows) if x['record_id'] in endpointids and x['sequence_group'] in fold[role+'_groups']],int);tr,va,te=[idx(s) for s in ['train','validation','test']];trids={rows[i]['record_id'] for i in tr};teids={rows[i]['record_id'] for i in te};trpairs=[p for p in pairs if p['reference_id'] in trids];tepairs=[p for p in pairs if p['reference_id'] in teids]
 for role,ii in [('train',tr),('validation',va),('test',te)]:
  for i in ii:membership.append(dict(outer_fold=f,role=role,record_id=rows[i]['record_id'],sequence_group=rows[i]['sequence_group'],assay_id=rows[i]['assay_id']))
 def add(method,seed,values):
  for p,z in zip(tepairs,values):prediction.append(dict(outer_fold=f,model=method,seed=seed,pair_id=p['pair_id'],reference_id=p['reference_id'],changed_id=p['changed_id'],sequence_group=rows[lookup[p['reference_id']]]['sequence_group'],observed_difference=p['observed_difference'],predicted_difference=float(z),measurement_sd_reference=p['measurement_sd_reference'],measurement_sd_changed=p['measurement_sd_changed'],mean_difference_se=None))
 g=[rows[lookup[p['reference_id']]]['sequence_group'] for p in trpairs];ww=weights(g);dy=np.array([p['observed_difference'] for p in trpairs]);constant=float(np.average(dy,weights=ww));directions=np.sign(dy)*(abs(dy)>.02);mass={s:float(ww[directions==s].sum()) for s in [-1,0,1]};majority=min(mass,key=lambda s:(-mass[s],s))
 for name,value in [('zero',0.),('training_mean',constant),('training_majority',majority)]:add(name,0,np.repeat(value,len(tepairs)))
 # Direct pair ridge, intentionally low-capacity, with training-fold feature support and no outer endpoints.
 dd,sup=transform(raw,rows,tr,True,True)
 def feats(pp):
  ar=np.array([lookup[p['reference_id']] for p in pp]);bb=np.array([lookup[p['changed_id']] for p in pp]);return np.concatenate([(dd['X'][bb]-dd['X'][ar]).reshape(len(pp),-1),dd['X'][ar,:32,:5].reshape(len(pp),-1)],1)
 xx=feats(trpairs);columns=np.flatnonzero(np.ptp(xx,axis=0)>0);ridge=Ridge(alpha=1.);ridge.fit(xx[:,columns],dy,sample_weight=ww);out=RUN/'fits'/f'b2/fold{f}/pair_ridge';out.mkdir(parents=True,exist_ok=True)
 with (out/'model.pkl').open('wb') as ff:pickle.dump(dict(model=ridge,columns=columns,support=sup,training_pair_ids=[p['pair_id'] for p in trpairs]),ff)
 add('pair_ridge',0,ridge.predict(feats(tepairs)[:,columns]));del dd
 for kind,method,lambdas in [('gnn','gnn_activity',[0.]),('gnn','gnn_pair',[1.,10.]),('nongraph','nongraph_pair',[1.,10.])]:
  dev=[]
  for lam in lambdas:
   rec,p=fit(f'b2/fold{f}/development/{method}',kind,True,True,CONFIGS[0],701,rows,raw,tr,va,np.array([],int),groupkey='sequence_group',pairs=pairs,lam=lam);results.append(rec);dev.append(rec)
  lam=min(dev,key=lambda r:(r['validation_score'],r['lambda_pair']))['lambda_pair'];selection.append(dict(outer_fold=f,method=method,lambda_pair=lam,training_mean_effect=constant,training_majority_direction=majority,train_pairs=len(trpairs),test_pairs=len(tepairs)))
  for seed in SEEDS:
   rec,p=fit(f'b2/fold{f}/final/{method}',kind,True,True,CONFIGS[0],seed,rows,raw,tr,va,te,groupkey='sequence_group',pairs=pairs,lam=lam);results.append(rec);endpoint=dict(zip([rows[i]['record_id'] for i in te],p['test']));add(method,seed,[endpoint[p['changed_id']]-endpoint[p['reference_id']] for p in tepairs])
write(RUN/'b2_fits.json',results);write(RUN/'b2_selection.json',selection);pd.DataFrame(membership).to_csv(RUN/'b2_membership.csv',index=False);pd.DataFrame(prediction).to_csv(RUN/'b2_predictions.csv',index=False)
