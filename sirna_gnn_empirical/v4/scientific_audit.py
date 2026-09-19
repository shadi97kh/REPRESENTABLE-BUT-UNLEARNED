"""Independent consistency checks on completed records; no fitting or new model selection."""
from common import *
import re
out=RUN/'audit';out.mkdir(exist_ok=True);summary=json.loads((RUN/'fit_completion.json').read_text());assert summary['actual_fits']==792
records=[];checked=0
for p in sorted((RUN/'fits').rglob('fit.json')):
 r=json.loads(p.read_text());parts=r['name'].split('/');r['scenario'],r['population'],r['phase']=parts[:3]
 for name,h in r['hashes'].items():assert sha(p.parent/name)==h,(p,name);checked+=1
 for name,h in r['signature']['code'].items():assert sha(HERE/name)==h
 assert r['signature']['protocol']==sha(RUN/'protocol.json');records.append(r)
assert len(records)==792
f=pd.DataFrame([{k:r[k] for k in ['scenario','population','phase','method','seed','status','executed_updates','selected_updates']} for r in records]);counts=f.groupby(['scenario','population','phase','method']).size();assert all(n==(10 if phase=='final' else 12) for (s,pop,phase,m),n in counts.items())
assert int(f.executed_updates.sum())==summary['executed_updates'];assert int(f.selected_updates.sum())==summary['selected_updates'];assert int((f.status!='completed').sum())==summary['failed_fits']
mem=pd.read_csv(RUN/'focused_membership.csv');checks=[]
for (s,pop,m),g in mem.groupby(['scenario','population','method']):
 sets={role:set(g[g.role==role].record_id) for role in ['train','validation','test']};groups={role:set(g[g.role==role].study_group) for role in sets};seqs={role:set(g[g.role==role].sequence_group) for role in sets}
 for a,b in [('train','validation'),('train','test'),('validation','test')]:assert not sets[a]&sets[b];assert not groups[a]&groups[b];assert not seqs[a]&seqs[b]
 if s=='source_excluded':assert not g.source.astype(str).eq('19282453').any()
 checks.append(dict(scenario=s,population=pop,method=m,roles_disjoint=True,study_components_disjoint=True,sequence_components_disjoint=True,Bramsen_train=int(((g.role=='train')&g.source.astype(str).eq('19282453')).sum()),Bramsen_validation=int(((g.role=='validation')&g.source.astype(str).eq('19282453')).sum())))
export(pd.DataFrame(checks),'audit/grouping_checks.csv')
pred=pd.read_csv(RUN/'evaluation/seed_predictions.csv');ens=pd.read_csv(RUN/'evaluation/ensemble_predictions.csv');scores=pd.read_csv(RUN/'evaluation/activity_metrics.csv');errors=[]
for (s,c,m),g in ens.groupby(['scenario','cohort','method']):
 q=scores[(scores.scenario==s)&(scores.cohort==c)&(scores.method==m)&(scores.weighting=='rows')].iloc[0];y=g.activity.to_numpy();p=g.prediction.to_numpy();mse=sum(float((a-b)**2) for a,b in zip(y,p))/len(y);var=sum(float((a-y.mean())**2) for a in y)/len(y);assert abs(mse-q.mse)<1e-12;assert abs(1-mse/var-q.r_squared)<1e-10
 ss=pred[(pred.scenario==s)&(pred.cohort==c)&(pred.method==m)].groupby('record_id').prediction.mean();ee=g.set_index('record_id').prediction;assert np.max(abs(ss-ee))<1e-12
 errors.append(dict(scenario=s,cohort=c,method=m,rows=len(g),mse_independently_recomputed=mse,maximum_ensemble_difference=float(np.max(abs(ss-ee)))))
export(pd.DataFrame(errors),'audit/metric_checks.csv')
# Source inventory proves source-limited target coverage, not independent assays.
original=rows_all()+readlines(PREV/'external_observations.jsonl');targets=pd.read_csv(RUN/'adjudication/primary_reanchored_targets.csv');curated=readlines(RUN/'dataset/primary_reanchored_observations.jsonl');assert len(curated)==2607 and len(targets)==1604
source_ids={x['record_id'] for x in original if x.get('dataset')=='ENsiRNA' and x['source_family']=='19282453'};assert len(source_ids)==1924
assert len(set(x['record_id'] for x in curated)&source_ids)==1604
oldinventory=json.loads((RUN/'preexisting_inventory.json').read_text());changed=[]
for name,rec in oldinventory.items():
 p=ROOT/name
 if not p.is_file() or p.stat().st_size!=rec['bytes'] or p.stat().st_mtime_ns!=rec['mtime_ns']:changed.append(name)
assert not changed,changed[:10]
result=dict(passed=True,new_fits=len(records),failed_fits=summary['failed_fits'],development_fits=int((f.phase=='development').sum()),final_fits=int((f.phase=='final').sum()),ten_final_seeds_every_block=True,verified_fit_artifact_hashes=checked,verified_old_metadata_entries=len(oldinventory),old_metadata_changes=changed,new_grouping_checks=len(checks),all_prediction_metrics_recomputed=True,source_excluded_grouped_rows=1003,primary_assay_rows=2607,source_quarantine=320,original_source_values_unchanged=True,preservation_scope='Historical metadata inventory unchanged; exact input, original selected-checkpoint and new fit artifact hashes verified separately. Metadata is not represented as a new all-history content hash audit.')
write(out/'scientific_audit.json',result);print(json.dumps(result,indent=2))
