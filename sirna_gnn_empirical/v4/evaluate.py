"""Evaluate completed focused fits and reuse the verified Bramsen-excluded outer0 models."""
from engine import *
from metrics import metrics,bootstrap_difference,ranking
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';out=RUN/'evaluation';out.mkdir(exist_ok=True)
assert (RUN/'fit_completion.json').exists(),'Do not evaluate incomplete final-seed blocks'
rows,raw=load();lookup={x['record_id']:i for i,x in enumerate(rows)};targets=pd.read_csv(RUN/'adjudication/primary_reanchored_targets.csv').set_index('record_id').primary_activity.to_dict();protocol=json.loads((RUN/'protocol.json').read_text());methods=protocol['methods'];df=pd.read_csv(RUN/'focused_seed_predictions.csv');reused=[];checks=[];cache={};device=setup()
test=np.array([i for i,x in enumerate(rows) if x.get('split')=='test_APP' or x['record_id'].startswith('davis:')],int)
for method in methods:
 for seed in SEEDS:
  fp=V3/f'fits/activity/outer0/final/{method}/s{seed}/fit.json';rec=json.loads(fp.read_text());tr=np.array(rec['signature']['train']);va=np.array(rec['signature']['validation']);assert all(rows[i]['source_family']!='19282453' for i in np.r_[tr,va]);assert not set(np.r_[tr,va])&set(test)
  for name,h in rec['hashes'].items():assert sha(fp.parent/name)==h
  k=tuple(tr)
  if k not in cache:cache[k]=transform(raw,rows,tr)
  data,sup=cache[k];assert sup==json.loads((fp.parent/'preprocessing.json').read_text())
  if rec['status']!='completed':values=np.full(len(test),rec['training_equal_group_mean'])
  elif method=='chemistry_tree':
   with (fp.parent/'model.pkl').open('rb') as f:st=pickle.load(f)
   xx=classical_features(data,method)[:,st['columns']];values=st['model'].predict(xx[test]);del xx,st
  else:
   st=torch.load(fp.parent/'best.pt',map_location=device,weights_only=False);m=network(method,sup).to(device);m.load_state_dict(st['model']);values=predict(m,data,test,device,st['target_mean'],st['target_sd']);del m,st
  for i,z in zip(test,values):
   x=rows[i];reused.append(dict(scenario='source_excluded',population='deployment_reused_outer0',cohort='APP' if x.get('split')=='test_APP' else 'Davis_S7',record_id=x['record_id'],source=x['source_family'],study_group=x['study_group'],sequence_group=x['sequence_group'],method=method,seed=seed,prediction=float(z),activity=x['activity'],status=rec['status'],fit=str(fp.parent.relative_to(ROOT)),training_row_mean=rec['training_row_mean'],training_equal_group_mean=rec['training_equal_group_mean']))
  saved=dict(np.load(fp.parent/'predictions.npz'))
  for i,z in zip(saved['test_indices'],saved['test']):
   x=rows[int(i)]
   if x['record_id'] not in targets:continue
   reused.append(dict(scenario='primary_reanchored',population='outer0_reused',cohort='ENsiRNA_grouped',record_id=x['record_id'],source=x['source_family'],study_group=x['study_group'],sequence_group=x['sequence_group'],method=method,seed=seed,prediction=float(z),activity=targets[x['record_id']],status=rec['status'],fit=str(fp.parent.relative_to(ROOT)),training_row_mean=rec['training_row_mean'],training_equal_group_mean=rec['training_equal_group_mean']))
  checks.append(dict(method=method,seed=seed,fit=str(fp.relative_to(ROOT)),fit_record_sha256=sha(fp),checkpoint_sha256=rec['hashes'][rec['checkpoint']],Bramsen_training_count=0,Bramsen_validation_count=0,original_selected_updates=rec['selected_updates'],new_updates=0,scope='Reuse original full-grid outer0 selection and retained validation; retrospective new-cohort inference and primary-label rescoring only'))
df=pd.concat([df,pd.DataFrame(reused)],ignore_index=True);assert not df.duplicated(['scenario','cohort','method','seed','record_id']).any()
# Constants depend on training partition, not held-out means; identical across seeds/methods for a population.
constants=[]
for (scenario,population),g in df.groupby(['scenario','population']):
 for name in ['training_row_mean','training_equal_group_mean']:
  assert g[name].max()-g[name].min()<1e-12
  q=g.drop_duplicates('record_id').copy();q['method']=name;q['seed']=0;q['prediction']=g[name].iloc[0];q['fit']='training-partition constant';constants.append(q)
df=pd.concat([df]+constants,ignore_index=True);export(df,'evaluation/seed_predictions.csv');write(out/'reused_fits.json',checks)
keys=['scenario','cohort','method','record_id'];en=df.groupby(keys,as_index=False).agg(prediction=('prediction','mean'),activity=('activity','first'),source=('source','first'),study_group=('study_group','first'),sequence_group=('sequence_group','first'),population=('population','first'),members=('seed','nunique'))
assert all(en[en.method==m].members.eq(10).all() for m in methods)
export(en,'evaluation/ensemble_predictions.csv');scores=[];seeds=[];sources=[]
for (scenario,cohort,method),g in en.groupby(['scenario','cohort','method']):
 for weighting in ['rows','equal_study_component']:
  w=None if weighting=='rows' else weights(g.study_group);scores.append(dict(scenario=scenario,cohort=cohort,method=method,weighting=weighting,sources=g.source.nunique(),sequence_components=g.sequence_group.nunique(),study_components=g.study_group.nunique(),**metrics(g.activity,g.prediction,w)))
 if cohort=='ENsiRNA_grouped':
  for source,q in g.groupby('source'):sources.append(dict(scenario=scenario,source=source,method=method,study_components=q.study_group.nunique(),sequence_components=q.sequence_group.nunique(),**metrics(q.activity,q.prediction)))
for (scenario,cohort,method,seed),g in df.groupby(['scenario','cohort','method','seed']):
 for weighting in ['rows','equal_study_component']:seeds.append(dict(scenario=scenario,cohort=cohort,method=method,seed=seed,weighting=weighting,**metrics(g.activity,g.prediction,None if weighting=='rows' else weights(g.study_group))))
# Oracle means are distinct diagnostic baselines; no seed variability assigned.
for (scenario,cohort),g in en[en.method=='corrected_gnn'].groupby(['scenario','cohort']):
 for weighting in ['rows','equal_study_component']:
  w=np.ones(len(g)) if weighting=='rows' else weights(g.study_group);mu=np.average(g.activity,weights=w);scores.append(dict(scenario=scenario,cohort=cohort,method='oracle_test_mean_diagnostic',weighting=weighting,sources=g.source.nunique(),sequence_components=g.sequence_group.nunique(),study_components=g.study_group.nunique(),**metrics(g.activity,np.full(len(g),mu),w)))
sc=pd.DataFrame(scores);sm=pd.DataFrame(seeds);export(sc,'evaluation/activity_metrics.csv');export(sm,'evaluation/seed_metrics.csv');export(pd.DataFrame(sources),'evaluation/per_source_metrics.csv');st=sm.groupby(['scenario','cohort','method','weighting'],as_index=False).agg(seed_count=('seed','nunique'),mean_seed_mse=('mse','mean'),sd_seed_mse=('mse','std'),min_seed_mse=('mse','min'),max_seed_mse=('mse','max'));export(st,'evaluation/seed_stability.csv')
contrasts=[];paired=[]
for (scenario,cohort),g in en.groupby(['scenario','cohort']):
 meta=g.drop_duplicates('record_id').set_index('record_id');wide=g.pivot(index='record_id',columns='method',values='prediction').join(meta[['activity','sequence_group','study_group']]).reset_index()
 for b in ['corrected_no_message','chemistry_tree','token_cnn','training_equal_group_mean']:
  for weighting in ['rows','equal_groups']:
   res,_=bootstrap_difference(wide,'corrected_gnn',b,weighting=weighting);contrasts.append(dict(scenario=scenario,cohort=cohort,**res))
 for b in ['corrected_no_message','chemistry_tree','token_cnn']:
  z=sm[(sm.scenario==scenario)&(sm.cohort==cohort)&(sm.weighting=='rows')];a=z[z.method=='corrected_gnn'].set_index('seed').mse;bb=z[z.method==b].set_index('seed').mse
  for seed in SEEDS:paired.append(dict(scenario=scenario,cohort=cohort,a='corrected_gnn',b=b,seed=seed,mse_difference=a[seed]-bb[seed]))
export(pd.DataFrame(contrasts),'evaluation/direct_comparisons.csv');export(pd.DataFrame(paired),'evaluation/seed_paired_differences.csv')
# Identical rows for released-vs-primary-label scoring; no claim that rescoring repairs training.
old=pd.read_csv(V3/'evaluated/activity_ensemble_predictions.csv');old=old[(old.cohort=='ENsiRNA_grouped')&old.method.isin(methods)].copy();newids=set(en[(en.scenario=='primary_reanchored')&(en.cohort=='ENsiRNA_grouped')].record_id);matched=[]
for method,g in old[old.record_id.isin(newids)].groupby('method'):
 for target in ['released','primary_reanchored']:
  yy=np.array([targets.get(r,y) if target=='primary_reanchored' else y for r,y in zip(g.record_id,g.activity)]);matched.append(dict(training='historical_released',evaluation_target=target,method=method,**metrics(yy,g.prediction)))
export(pd.DataFrame(matched),'evaluation/matched_subset_rescoring.csv')
# Candidate ranking for the sensitivity deployments, joined to the original exact assay IDs.
assays={x['record_id']:x.get('assay_id') for x in rows};rr=[]
for (scenario,method),g in en[en.cohort=='APP'].groupby(['scenario','method']):
 g=g.copy();g['assay_id']=g.record_id.map(assays)
 for pool,q in g.groupby('assay_id'):rr.append(dict(scenario=scenario,method=method,pool=pool,**ranking(q)))
export(pd.DataFrame(rr),'evaluation/ranking_pools.csv');write(out/'summary.json',dict(new_fits=json.loads((RUN/'fit_completion.json').read_text())['actual_fits'],reused_final_fits=len(checks),source_excluded_grouped_rows=1003,primary_reanchored_grouped_rows=2607,APP_rows=1839,S7_rows=118,B2_fits_reused=True,B3_fits_reused=True,seed_predictions=len(df),ensemble_predictions=len(en),interpretation='Retrospective grouped and deployment sensitivity; no fresh independent study, no source-causal attribution when selection/partition roles differ.'))
print(sc[(sc.weighting=='rows')&sc.method.isin(methods)][['scenario','cohort','method','mse','r_squared','prediction_sd']].to_string(index=False))
