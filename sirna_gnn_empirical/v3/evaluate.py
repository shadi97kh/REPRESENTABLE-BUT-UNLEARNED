"""Metrics, direct matched uncertainty, all-seed stability and exact-tie rankings."""
from common import *
from metrics import *
rows=rows_all()+readlines(PREV/'external_observations.jsonl');meta=pd.DataFrame(rows);pred=pd.read_csv(RUN/'activity_predictions.csv');pred=pred.merge(meta,on='record_id',validate='many_to_one');allmetrics=[];seedmetrics=[];ensembles=[];comparisons=[];seedcontrasts=[];losses=[]
for cohort,p in pred.groupby('cohort'):
 # Each original observation is tested exactly once per method/seed in outer grouped evaluation.
 assert not p.duplicated(['record_id','method','seed']).any()
 for (method,seed),z in p.groupby(['method','seed']):
  for wn,g in [('rows',None),('equal_study','study_group'),('equal_sequence','sequence_group')]:seedmetrics.append(dict(cohort=cohort,method=method,seed=seed,weighting=wn,**metrics(z.activity,z.prediction,weights(z[g]) if g else None)))
 e=p.groupby(['record_id','method'],as_index=False).prediction.mean().merge(meta,on='record_id',validate='many_to_one');e['cohort']=cohort
 c=p.drop_duplicates('record_id')
 for col,name in [('training_row_mean','training_row_mean'),('training_equal_group_mean','training_equal_group_mean')]:
  q=c[meta.columns.tolist()].copy();q['method']=name;q['prediction']=c[col].to_numpy();q['cohort']=cohort;e=pd.concat([e,q],ignore_index=True)
 for method,z in e.groupby('method'):
  for wn,g in [('rows',None),('equal_study','study_group'),('equal_sequence','sequence_group')]:allmetrics.append(dict(cohort=cohort,method=method,weighting=wn,sequence_components=z.sequence_group.nunique(),study_linked_groups=z.study_group.nunique(),source_families=z.source_family.nunique(),**metrics(z.activity,z.prediction,weights(z[g]) if g else None)))
 c=e.drop_duplicates('record_id')
 for wn,g in [('rows',None),('equal_study','study_group'),('equal_sequence','sequence_group')]:
  w=weights(c[g]) if g else np.ones(len(c));mu=np.average(c.activity,weights=w);allmetrics.append(dict(cohort=cohort,method='oracle_test_mean_diagnostic',weighting=wn,sequence_components=c.sequence_group.nunique(),study_linked_groups=c.study_group.nunique(),source_families=c.source_family.nunique(),**metrics(c.activity,np.repeat(mu,len(c)),w)))
 ensembles.append(e);wide=e.pivot(index='record_id',columns='method',values='prediction').reset_index().merge(meta,on='record_id',validate='one_to_one')
 pairs=[('original_gnn','original_no_message'),('corrected_gnn','corrected_no_message'),('corrected_gnn','chemistry_tree'),('original_gnn','chemistry_tree'),('corrected_gnn','original_gnn'),('corrected_no_message','original_no_message'),('corrected_gnn','token_cnn'),('corrected_gnn','training_equal_group_mean'),('corrected_gnn','training_row_mean'),('corrected_gnn','guide_ridge'),('corrected_gnn','pairwise_ridge'),('corrected_gnn','original_no_message')]
 for a,b in pairs:
  for wn,group in [('rows','sequence_group'),('equal_groups','study_group' if cohort=='ENsiRNA_grouped' else 'sequence_group')]:
   d,l=bootstrap_difference(wide,a,b,group,wn);d.update(cohort=cohort,resampled_unit=group);comparisons.append(d);l=l.merge(meta[['record_id']+[k for k in ['sequence_group','study_group','source_family','assay_id'] if k not in l]],on='record_id',validate='one_to_one');l['cohort']=cohort;l['a']=a;l['b']=b;l['weighting']=wn;losses.append(l)
  aa=p[p.method==a];bb=p[p.method==b]
  if b in ['training_row_mean','training_equal_group_mean']:
   bb=aa.copy();bb['prediction']=bb[b]
  deterministic_b=b in ['guide_ridge','pairwise_ridge']
  for seed in sorted(set(aa.seed) if deterministic_b else set(aa.seed)&set(bb.seed)):
   aa1=aa[aa.seed==seed];bb1=bb[bb.seed==(0 if deterministic_b else seed)];q=aa1[['record_id','activity','prediction','study_group','sequence_group']].merge(bb1[['record_id','prediction']],on='record_id',suffixes=('_a','_b'),validate='one_to_one')
   for wn,g in [('rows',None),('equal_study','study_group'),('equal_sequence','sequence_group')]:
    w=weights(q[g]) if g else np.ones(len(q));seedcontrasts.append(dict(cohort=cohort,a=a,b=b,seed=seed,weighting=wn,difference=float(np.average((q.prediction_a-q.activity)**2-(q.prediction_b-q.activity)**2,weights=w))))
 # Controlled difference-in-differences in squared loss, conditional on four realized ensembles.
 y=wide.activity;v=(wide.corrected_gnn-y)**2-(wide.corrected_no_message-y)**2-(wide.original_gnn-y)**2+(wide.original_no_message-y)**2
 # Bootstrap the signed loss contrast directly, never a difference of separately constructed CIs.
 u=pd.DataFrame(dict(group=wide.sequence_group,value=v));a=u.groupby('group').value.agg(['sum','count']);rng=np.random.default_rng(90413);vals=[]
 for _ in range(100):
  ix=rng.integers(len(a),size=(100,len(a)));vals.extend(a['sum'].to_numpy()[ix].sum(1)/a['count'].to_numpy()[ix].sum(1))
 comparisons.append(dict(cohort=cohort,a='support_by_message_interaction',b='zero',difference=float(v.mean()),lower=float(np.quantile(vals,.025)),upper=float(np.quantile(vals,.975)),groups=len(a),weighting='rows',resampled_unit='sequence_group',bootstrap_resamples=10000,analysis_seed=90413,uncertainty='conditional four-ensemble loss interaction'))
ens=pd.concat(ensembles,ignore_index=True);identity=['record_id','activity','sequence_group','study_group','source_family','assay_id','dataset','split'];export(ens[['cohort','method','prediction']+[c for c in identity if c in ens]],'evaluated/activity_ensemble_predictions.csv');export(pred[[c for c in ['cohort','population','method','seed','prediction','training_row_mean','training_equal_group_mean','fit','status']+identity if c in pred]],'evaluated/activity_seed_predictions.csv');export(pd.DataFrame(allmetrics),'tables/activity_metrics.csv');export(pd.DataFrame(seedmetrics),'tables/activity_seed_metrics.csv');export(pd.DataFrame(comparisons),'tables/direct_comparisons.csv');export(pd.DataFrame(seedcontrasts),'tables/seed_paired_differences.csv');export(pd.concat(losses,ignore_index=True),'evaluated/observation_loss_differences.csv')
sm=pd.DataFrame(seedmetrics);ss=sm.groupby(['cohort','method','weighting']).mse.agg(['count','mean','std','min','max']).reset_index().rename(columns={'mean':'mean_seed_mse','std':'sd_seed_mse','min':'min_seed_mse','max':'max_seed_mse'});export(ss,'tables/seed_stability.csv')
sc=pd.DataFrame(seedcontrasts);export(sc.groupby(['cohort','a','b','weighting']).difference.agg(['count','mean','std','min','max']).reset_index(),'tables/seed_comparison_stability.csv')
rank=[]
for (method,assay),z in ens[ens.cohort=='APP'].groupby(['method','assay_id']):
 if len(z)>=5:rank.append(dict(method=method,pool=assay,**ranking(z)))
rk=pd.DataFrame(rank);export(rk,'tables/ranking_pools.csv');export(rk.groupby('method').agg(pools=('pool','nunique'),top5_mean_percentile=('top5_mean_percentile','mean'),selected_mean_activity=('selected_mean_activity','mean')).reset_index(),'tables/ranking_summary.csv')
w=rk.pivot(index='pool',columns='method',values='top5_mean_percentile');w['tree_minus_corrected_gnn']=w.chemistry_tree-w.corrected_gnn;w['corrected_gnn_minus_no_message']=w.corrected_gnn-w.corrected_no_message;export(w.reset_index(),'tables/ranking_direct_pool_differences.csv')
# Preserve dependence metadata: no pool-independent significance calculation.
pool={a:set(z.sequence_group) for a,z in ens[(ens.cohort=='APP')&(ens.method=='corrected_gnn')].groupby('assay_id')};export(pd.DataFrame([dict(pool_a=a,pool_b=b,shared_sequence_components=len(pool[a]&pool[b])) for a in pool for b in pool if a<b]),'tables/ranking_pool_overlap.csv')
bp=pd.read_csv(RUN/'pair_predictions.csv');be=bp.groupby(['pair_id','method'],as_index=False).prediction.mean().merge(bp.drop_duplicates('pair_id').drop(columns=['method','seed','prediction','fit','status']),on='pair_id',validate='many_to_one');bmetrics=[];bseeds=[];bcomparisons=[]
def bmetric(z,method):
 obs=np.where(abs(z.observed_difference)<=.02,0,np.sign(z.observed_difference));calls=np.where(abs(z.prediction)<=.02,0,np.sign(z.prediction));out=[]
 for wn,g in [('rows',None),('equal_sequence','sequence_group'),('equal_background','background_id')]:
  w=weights(z[g]) if g else np.ones(len(z));d=metrics(z.observed_difference,z.prediction,w)
  if method=='training_majority':
   for k in ['mse','mae','r_squared','correlation','prediction_sd','mean_bias','prediction_mean']:d[k]=None
  out.append(dict(method=method,weighting=wn,**d,sign_accuracy=float(np.average(obs==calls,weights=w)),non_tie_correct=int(((obs==calls)&(obs!=0)).sum()),non_ties=int((obs!=0).sum()),negative_calls=int((calls<0).sum()),zero_calls=int((calls==0).sum()),positive_calls=int((calls>0).sum())))
 return out
for method,z in be.groupby('method'):bmetrics.extend(bmetric(z,method))
for (method,seed),z in bp.groupby(['method','seed']):bseeds.extend([dict(seed=seed,**x) for x in bmetric(z,method)])
bwide=be.pivot(index='pair_id',columns='method',values='prediction').reset_index().merge(be.drop_duplicates('pair_id')[['pair_id','sequence_group','observed_difference']],on='pair_id',validate='one_to_one').rename(columns={'pair_id':'record_id','observed_difference':'activity'})
for a,b in [('pair_gnn','training_mean'),('pair_gnn','activity_gnn'),('pair_gnn','pair_no_message'),('pair_gnn','pair_ridge'),('pair_gnn','pair_tree'),('pair_gnn','zero')]:
 for wn in ['rows','equal_groups']:
  d,z=bootstrap_difference(bwide,a,b,weighting=wn);bcomparisons.append(d);z=z.merge(be.drop_duplicates('pair_id')[['pair_id','population','background_id','reference_id','changed_id']].rename(columns={'pair_id':'record_id'}),on='record_id',validate='one_to_one');export(z,f'evaluated/B2_{a}_minus_{b}_{wn}.csv')
export(be,'evaluated/pair_ensemble_predictions.csv');export(pd.DataFrame(bmetrics),'tables/pair_metrics.csv');export(pd.DataFrame(bseeds),'tables/pair_seed_metrics.csv');export(pd.DataFrame(bcomparisons),'tables/pair_comparisons.csv')
ep=pd.read_csv(RUN/'pair_endpoint_predictions.csv');em=[]
for (method,seed),z in ep.groupby(['method','seed']):em.append(dict(method=method,seed=seed,**metrics(z.activity,z.prediction,weights(z.sequence_group))))
for method,z in ep.groupby('method'):
 z=z.groupby('record_id',as_index=False).agg(prediction=('prediction','mean'),activity=('activity','first'),sequence_group=('sequence_group','first'));em.append(dict(method=method,seed='ensemble',**metrics(z.activity,z.prediction,weights(z.sequence_group))))
export(pd.DataFrame(em),'tables/pair_endpoint_metrics.csv')
fitrows=[];curves=[];gradients=[]
for p in sorted((RUN/'fits').rglob('fit.json')):
 r=json.loads(p.read_text());row={k:r.get(k) for k in ['name','method','seed','status','failure','training_n','validation_n','test_n','parameters','selected_epoch','selected_updates','executed_updates','epochs_executed','validation_score','target_mean','target_sd','lambda_pair','selection','wall_s','cpu_s','gpu_region_wall_s']};row['configuration']=r['config']['name'];fitrows.append(row)
 for h in r.get('history',[]):curves.append(dict(fit=r['name'],method=r['method'],seed=r['seed'],**h))
 gp=p.parent/'gradient_audit.json'
 if gp.exists():
  for g in json.loads(gp.read_text()):gradients.append(dict(fit=r['name'],method=r['method'],seed=r['seed'],**g))
f=pd.DataFrame(fitrows);export(f,'tables/all_fit_summary.csv');export(pd.DataFrame(curves),'evaluated/all_learning_curves.csv');export(pd.DataFrame(gradients),'tables/relation_gradients_by_fit.csv');export(f[f.status!='completed'],'tables/failed_fits.csv')
expected=set(pd.read_csv(RUN/'experiment_matrix.csv').name);actual=set(f.name);assert actual==expected,(len(expected-actual),len(actual-expected));assert f[f.name.str.contains('/final/')&~f.method.isin(['guide_ridge','pairwise_ridge','pair_ridge'])].groupby(['name']).size().min()==1
write(RUN/'evaluation_summary.json',dict(actual_fits=len(f),failed_fits=int((f.status!='completed').sum()),executed_updates=int(f.executed_updates.sum()),selected_updates=int(f.selected_updates.sum()),final_seeds=SEEDS,core_matrix_completed=True,all_finite_unfavorable_seeds_retained=True,original_APP_B3_complete_rectangles=0,primary_B3=json.loads((RUN/'b3_primary/evaluation_summary.json').read_text()),published_baseline='See sources/reproduction_status.json; no exact reproduction claimed from missing dependencies',scientific_claim='Interpret actual tables; no positive conclusion is presumed.'))
print(pd.DataFrame(allmetrics).query("weighting=='rows'").to_string(index=False));print(pd.DataFrame(bmetrics).query("weighting=='rows'").to_string(index=False))
