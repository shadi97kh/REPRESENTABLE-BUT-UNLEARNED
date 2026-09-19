"""Exact-ID cross-campaign reconciliation under one fractional-tie ranking definition."""
from common import *
from metrics import ranking
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';out=RUN/'ranking';out.mkdir(exist_ok=True)
meta=pd.DataFrame([dict(record_id=x['record_id'],activity=x['activity'],assay_id=x['assay_id'],sequence_group=x['sequence_group']) for x in readlines(OLD/'test_APP_observations.jsonl')]);assert meta.record_id.is_unique and meta.assay_id.nunique()==17
versions=[]
z=pd.read_csv(OLD/'heldout_predictions.csv');z=z[z.record_id.isin(meta.record_id)].copy();z['campaign']='v1';z['phase']='original';z['arm']='original';z['member_count']=z['number_of_fits'];versions.append(z[['campaign','phase','model','arm','record_id','prediction','member_count']])
z=pd.read_csv(PREV/'all_B1_predictions.csv');z=z[z.record_id.isin(meta.record_id) & pd.to_numeric(z.seed,errors='coerce').notna()];z=z.groupby(['phase','model','arm','record_id'],dropna=False).agg(prediction=('prediction','mean'),member_count=('seed','nunique')).reset_index();z['campaign']='v2';versions.append(z[['campaign','phase','model','arm','record_id','prediction','member_count']])
z=pd.read_csv(V3/'evaluated/activity_ensemble_predictions.csv');z=z[z.cohort=='APP'].copy();z['campaign']='v3';z['phase']='group_selected';z['model']=z.method;z['arm']=z.method;z['member_count']=z.method.map(lambda m:1 if 'ridge' in m or 'mean' in m else 10);versions.append(z[['campaign','phase','model','arm','record_id','prediction','member_count']])
old1=pd.read_csv(OLD/'assay_metrics.csv');old2=pd.read_csv(PREV/'ranking_by_assay.csv');old2=old2[old2.seed.astype(str)=='ensemble'];old3=pd.read_csv(V3/'tables/ranking_pools.csv')
allz=pd.concat(versions,ignore_index=True).merge(meta,on='record_id',validate='many_to_one');assert not allz.duplicated(['campaign','phase','model','arm','record_id']).any();records=[];coverage=[]
for (v,phase,model,arm),data in allz.groupby(['campaign','phase','model','arm']):
 ids=set(data.record_id);coverage.append(dict(campaign=v,phase=phase,model=model,arm=arm,predicted_rows=len(ids),missing_expected_rows=len(set(meta.record_id)-ids),exact_same_1839_rows=ids==set(meta.record_id)))
 for pool,g in data.groupby('assay_id'):
  stats=ranking(g);prior=None
  if v=='v3':
   q=old3[(old3.method==model)&(old3.pool==pool)]
   if len(q)==1:prior=float(q.iloc[0].top5_mean_percentile);assert abs(prior-stats['top5_mean_percentile'])<1e-12
  elif v=='v2':
   q=old2[(old2.phase==phase)&(old2.model==model)&(old2.arm==arm)&(old2.assay_id==pool)]
   if len(q)==1:prior=float(q.iloc[0].top5_percentile)
  else:
   if 'assay_id' in old1 and 'model' in old1:
    q=old1[(old1.model==model)&(old1.assay_id==pool)]
    col=next((c for c in ['top5_percentile','top5_mean_percentile'] if c in q),None)
    if len(q)==1 and col:prior=float(q.iloc[0][col])
  records.append(dict(campaign=v,phase=phase,model=model,arm=arm,pool=pool,membership_sha256=hashlib.sha256('\n'.join(sorted(g.record_id)).encode()).hexdigest(),ensemble_members_min=int(g.member_count.min()),ensemble_members_max=int(g.member_count.max()),legacy_reported_top5=prior,common_minus_legacy=None if prior is None else stats['top5_mean_percentile']-prior,**stats))
records=pd.DataFrame(records);export(records,'ranking/all_pool_comparisons.csv');export(pd.DataFrame(coverage),'ranking/coverage.csv')
summary=records.groupby(['campaign','phase','model','arm'],as_index=False).agg(pools=('pool','nunique'),common_top5_percentile=('top5_mean_percentile','mean'),historical_top5_percentile=('legacy_reported_top5','mean'),selected_activity=('selected_mean_activity','mean'),minimum_members=('ensemble_members_min','min'),maximum_members=('ensemble_members_max','max'))
export(summary,'ranking/campaign_summary.csv')
# All reported pool memberships agree across model/campaign paths; pools still overlap in sequences.
assert records.groupby('pool').membership_sha256.nunique().max()==1
pools={k:set(g.sequence_group) for k,g in meta.groupby('assay_id')};overlap=[dict(pool_a=a,pool_b=b,shared_sequence_components=len(pools[a]&pools[b])) for i,a in enumerate(sorted(pools)) for b in sorted(pools)[i+1:]];export(pd.DataFrame(overlap),'ranking/pool_overlap.csv');export(meta,'ranking/exact_pool_membership.csv')
for name,p in [('v1_historical_assay_metrics.csv',OLD/'assay_metrics.csv'),('v2_historical_rankings.csv',PREV/'ranking_by_assay.csv'),('v3_historical_rankings.csv',V3/'tables/ranking_pools.csv')]:
 import shutil;shutil.copy2(p,out/name)
write(out/'summary.json',dict(pools=17,observations=1839,identical_record_and_pool_membership=True,response_direction='larger measured inhibition preferred; labels unchanged',common_definition='Average-rank outcome percentile; fractional selection mass among exact prediction ties at top-five boundary; constant expectation 0.5',pool_dependence='Related patent family and overlapping sequence components; no independent-pool significance test',protocol_changes=['v1 original limited training/validation schedule and three final neural seeds','v2 frozen diagnostic arms, factorial refits and separate all-ENsiRNA grouped-selected deployment; three seeds','v3 corrected support-by-message control, training-only feature/context masks, grouped selection, 500-epoch cap and ten final seeds'],causal_attribution='Several protocol and ensemble choices changed jointly. Ranking changes are descriptive and limit stability claims; they do not identify one cause or render the statistic meaningless.',new_fits=0))
print(summary[(summary.model.isin(['gnn','chemistry_tree','nongraph','token_cnn','corrected_gnn','corrected_no_message'])) & ~((summary.campaign=='v2')&(summary.phase=='frozen'))].to_string(index=False))
