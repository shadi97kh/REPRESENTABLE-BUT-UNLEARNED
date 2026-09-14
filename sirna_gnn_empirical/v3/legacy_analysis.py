"""Execute missing diagnostics from unchanged saved predictions before new fits."""
from common import *
from metrics import *
rows=rows_all();ext=readlines(PREV/'external_observations.jsonl');meta=pd.DataFrame(rows+ext);lookup=meta.set_index('record_id');protocol=json.loads((PREV/'protocol.json').read_text());result=[];predictions=[];pairs_out=[];seed_rows=[]
allp=pd.read_csv(PREV/'all_B1_predictions.csv');datasets=[]
for (phase,arm,split),p in allp.groupby(['phase','arm','split']):datasets.append((f'v2_{phase}_{arm}',split,p,None))
for fname,label,cohort in [('grouped_predictions.csv','v2_grouped','ENsiRNA_grouped'),('external_predictions.csv','v2_deployment','Davis_S7')]:
 p=pd.read_csv(PREV/fname);datasets.append((label,cohort,p,'grouped' if 'grouped' in fname else 'deployment'))
for label,cohort,p,role in datasets:
 ids=list(p.record_id.unique());m=lookup.loc[ids].reset_index();merged=p.merge(m,on='record_id',suffixes=('','_metadata'),validate='many_to_one')
 ensemble=merged.groupby(['record_id','model'],as_index=False).prediction.mean();ensemble=ensemble.merge(m,on='record_id',validate='many_to_one');ensemble['protocol']=label;ensemble['cohort']=cohort
 for weightname,key in [('rows',None),('equal_sequence','sequence_group'),('equal_study','study_group')]:
  if key and (key not in m or m[key].isna().any()):continue
  for model,z in ensemble.groupby('model'):
   d=metrics(z.activity,z.prediction,weights(z[key]) if key else None);result.append(dict(protocol=label,cohort=cohort,model=model,weighting=weightname,sequence_components=z.sequence_group.nunique(),source_families=z.source_family.nunique() if 'source_family' in z else None,study_linked_groups=z.study_group.nunique() if 'study_group' in z else None,**d))
 constants=[]
 for rid in ids:
  if role=='grouped':
   f=int(p.loc[p.record_id==rid,'outer_fold'].iloc[0]);groups=protocol['grouped_folds'][f]['train_groups'];train=meta[meta.study_group.isin(groups)]
  elif label.startswith('v2_deployment') or role=='deployment':train=meta[meta.dataset=='ENsiRNA']
  else:train=meta[meta.split=='train']
  constants.append((rid,float(train.activity.mean()),float(np.average(train.activity,weights=weights(train.study_group)))))
 for k,name in [(1,'training_row_mean'),(2,'training_equal_study_mean')]:
  z=m.copy();z['model']=name;z['prediction']=[a[k] for a in constants];z['protocol']=label;z['cohort']=cohort;ensemble=pd.concat([ensemble,z],ignore_index=True)
  for wn,gk in [('rows',None),('equal_sequence','sequence_group'),('equal_study','study_group')]:
   if gk and (gk not in z or z[gk].isna().any()):continue
   result.append(dict(protocol=label,cohort=cohort,model=name,weighting=wn,**metrics(z.activity,z.prediction,weights(z[gk]) if gk else None)))
 for wn,gk in [('rows',None),('equal_sequence','sequence_group'),('equal_study','study_group')]:
  if gk and (gk not in m or m[gk].isna().any()):continue
  ww=weights(m[gk]) if gk else np.ones(len(m));c=np.average(m.activity,weights=ww);result.append(dict(protocol=label,cohort=cohort,model='oracle_test_mean_diagnostic',weighting=wn,**metrics(m.activity,np.full(len(m),c),ww)))
 for (model,seed),z in merged.groupby(['model','seed']):seed_rows.append(dict(protocol=label,cohort=cohort,model=model,seed=seed,**metrics(z.activity,z.prediction)))
 predictions.append(ensemble)
 wide=ensemble.pivot(index='record_id',columns='model',values='prediction').reset_index().merge(m,on='record_id',validate='one_to_one')
 for a,b in [('gnn','chemistry_tree'),('gnn','nongraph'),('gnn','training_equal_study_mean')]:
  if a not in wide or b not in wide:continue
  d,z=bootstrap_difference(wide,a,b);d.update(protocol=label,cohort=cohort);pairs_out.append(d);z['protocol']=label;z['cohort']=cohort;z['a']=a;z['b']=b;export(z,f'legacy/losses/{label}_{cohort}_{a}_minus_{b}.csv')
allensembles=pd.concat(predictions,ignore_index=True);export(allensembles,'legacy/activity_predictions.csv');export(pd.DataFrame(result),'legacy/activity_metrics.csv');export(pd.DataFrame(seed_rows),'legacy/seed_metrics.csv');export(pd.DataFrame(pairs_out),'legacy/direct_comparisons.csv')
rank=[];overlap=[]
for (pro,co,model,assay),z in allensembles[allensembles.split=='test_APP'].groupby(['protocol','cohort','model','assay_id']):
 if len(z)>=5:rank.append(dict(protocol=pro,cohort=co,model=model,assay_id=assay,**ranking(z)))
app=meta[meta.split=='test_APP']; pools={a:set(z.sequence_group) for a,z in app.groupby('assay_id') if len(z)>=5}
for a in pools:
 for b in pools:
  if a<b:overlap.append(dict(pool_a=a,pool_b=b,shared_sequence_components=len(pools[a]&pools[b])))
export(pd.DataFrame(rank),'legacy/ranking_by_pool.csv');export(pd.DataFrame(overlap),'legacy/pool_overlap.csv')
b2=pd.read_csv(PREV/'b2_predictions.csv');bm=[];be=b2.groupby(['pair_id','model'],as_index=False).agg(predicted_difference=('predicted_difference','mean'));be=be.merge(b2.drop_duplicates('pair_id').drop(columns=['model','seed','predicted_difference']),on='pair_id',validate='many_to_one')
for model,z in be.groupby('model'):
 for wn,w in [('rows',None),('equal_sequence',weights(z.sequence_group))]:
  d=metrics(z.observed_difference,z.predicted_difference,w);obs=np.where(abs(z.observed_difference)<=.02,0,np.sign(z.observed_difference));pr=np.where(abs(z.predicted_difference)<=.02,0,np.sign(z.predicted_difference));bm.append(dict(model=model,weighting=wn,sign_accuracy=float(np.mean(obs==pr)),non_tie_correct=int(((obs==pr)&(obs!=0)).sum()),non_ties=int((obs!=0).sum()),negative_calls=int((pr<0).sum()),zero_calls=int((pr==0).sum()),positive_calls=int((pr>0).sum()),**d))
export(be,'legacy/b2_predictions.csv');export(pd.DataFrame(bm),'legacy/b2_metrics.csv')
checks=[]
for x in readlines(OLD/'eligible_pairs.jsonl'):
 a=lookup.loc[x['reference_id']];b=lookup.loc[x['changed_id']];assert a.sequence_group==b.sequence_group;assert abs(b.activity-a.activity-x['observed_difference'])<1e-9
 checks.append(dict(pair_id=x['pair_id'],reference_id=x['reference_id'],changed_id=x['changed_id'],sequence_group=a.sequence_group,background_id=x['background_id'],orientation_verified=True,replicate_n_reference=a.replicate_n,replicate_n_changed=b.replicate_n,measurement_sd_reference=x['measurement_sd_reference'],measurement_sd_changed=x['measurement_sd_changed']))
export(pd.DataFrame(checks),'legacy/pair_identity_audit.csv')
write(RUN/'legacy/summary.json',dict(activity_metric_rows=len(result),seed_rows=len(seed_rows),conditional_comparisons=len(pairs_out),ranking_rows=len(rank),measured_pairs=len(checks),review_corrections=['S7 token modestly exceeds oracle test mean descriptively; oracle mean is not deployable or a lower bound.','Ensemble and mean member MSE are distinct. Fixed-ensemble bootstrap is valid conditionally.','No source establishes that the support defect alone fully predicts failure.','All three listed corrected GNN seed MSE values exceed .1627; exact seed exports supplied.','Only one varying law does not invalidate the full five-law lower-bound subexperiment.']))
print(pd.DataFrame(result).query("cohort=='Davis_S7' and weighting=='rows'").to_string(index=False))
