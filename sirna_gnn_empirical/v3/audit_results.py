"""Independent export/identity audits after the frozen campaign; no fitting."""
from common import *
import subprocess
from metrics import metrics,ranking
D=ROOT/'docs/sirna_gnn_empirical/v3';checks=[]
def check(name,condition,details=None):
 checks.append(dict(check=name,passed=bool(condition),details=details));assert condition,(name,details)
expected=pd.read_csv(RUN/'experiment_matrix.csv');fits=[];filecount=0
for p in sorted((RUN/'fits').rglob('fit.json')):
 r=json.loads(p.read_text());fits.append(r)
 for n,h in r['hashes'].items():check_path=p.parent/n;assert sha(check_path)==h,check_path;filecount+=1
 sig=r['signature'];sets=[set(sig[k]) for k in ['train','validation','test']]
 check('fit partition row exclusion: '+r['name'],not any(sets[i]&sets[j] for i,j in [(0,1),(0,2),(1,2)]))
check('Exact frozen fit matrix',set(expected.name)=={r['name'] for r in fits},dict(expected=len(expected),actual=len(fits)))
rows=rows_all()+readlines(PREV/'external_observations.jsonl');meta=pd.DataFrame(rows)
for r in fits:
 sig=r['signature'];key='sequence_group' if r['name'].startswith('pair/') else 'study_group'
 sets=[{rows[i][key] for i in sig[k]} for k in ['train','validation','test']]
 check('protected grouping: '+r['name'],not any(sets[i]&sets[j] for i,j in [(0,1),(0,2),(1,2)]))
 if r['executed_updates'] and r['status']=='completed':
  check('optimization trajectory: '+r['name'],r['epochs_executed']>=50 and 0<r['selected_updates']<=r['executed_updates'])
  h=r['history'];check('history update count: '+r['name'],h[-1]['updates']==r['executed_updates'] and len(h)==r['epochs_executed'])
final=[r for r in fits if '/final/' in r['name']];blocks={}
for r in final:
 # Pair/outer/method and activity/outer/method are explicit directory components.
 q=r['name'].split('/');key='/'.join(q[:q.index('final')+2]);blocks.setdefault(key,[]).append(r)
for key,rs in blocks.items():
 deterministic=rs[0]['method'] in ['guide_ridge','pairwise_ridge','pair_ridge'];check('Final seed block: '+key,sorted(r['seed'] for r in rs)==([0] if deterministic else sorted(SEEDS)))
# Every published rounded metric is recomputed from row predictions with the stated weights.
ens=pd.read_csv(RUN/'evaluated/activity_ensemble_predictions.csv');table=pd.read_csv(RUN/'tables/activity_metrics.csv');maxerror=0.
for q in table.to_dict('records'):
 z=ens[(ens.cohort==q['cohort'])&(ens.method==(q['method'] if q['method']!='oracle_test_mean_diagnostic' else 'corrected_gnn'))];g={'rows':None,'equal_study':'study_group','equal_sequence':'sequence_group'}[q['weighting']];w=weights(z[g]) if g else np.ones(len(z));y=z.activity.to_numpy(float);p=z.prediction.to_numpy(float)
 if q['method']=='oracle_test_mean_diagnostic':p=np.repeat(np.average(y,weights=w),len(y))
 ym=np.average(y,weights=w);pm=np.average(p,weights=w);vy=np.average((y-ym)**2,weights=w);vp=np.average((p-pm)**2,weights=w);mse=np.average((p-y)**2,weights=w);m=dict(mse=mse,mae=np.average(abs(p-y),weights=w),r_squared=1-mse/vy if vy>1e-20 else None,correlation=np.average((y-ym)*(p-pm),weights=w)/np.sqrt(vy*vp) if vy*vp>1e-20 else None,prediction_sd=np.sqrt(vp),label_sd=np.sqrt(vy),mean_bias=pm-ym)
 for k in ['mse','mae','r_squared','correlation','prediction_sd','label_sd','mean_bias']:
  assert (m[k] is None)==pd.isna(q[k]),(q,k,'undefined mismatch')
  if m[k] is not None:maxerror=max(maxerror,abs(m[k]-q[k]));assert abs(m[k]-q[k])<1e-10,(q,k,m[k])
check('All activity metrics reproduce from full-precision exports',maxerror<1e-10,dict(maximum_absolute_error=maxerror))
ranks=pd.read_csv(RUN/'tables/ranking_pools.csv')
for q in ranks.to_dict('records'):
 z=ens[(ens.cohort=='APP')&(ens.method==q['method'])&(ens.assay_id==q['pool'])];v=ranking(z);assert abs(v['top5_mean_percentile']-q['top5_mean_percentile'])<1e-12
check('Both training constants have exact random-tie ranking',np.allclose(ranks[ranks.method.isin(['training_row_mean','training_equal_group_mean'])].top5_mean_percentile,.5,atol=1e-12))
pp=pd.read_csv(RUN/'pair_predictions.csv');lu=meta.set_index('record_id');unique=pp.drop_duplicates('pair_id');maxpair=0
for q in unique.itertuples():
 a,b=lu.loc[q.reference_id],lu.loc[q.changed_id];maxpair=max(maxpair,abs((b.activity-a.activity)-q.observed_difference));assert a.sequence_group==b.sequence_group==q.sequence_group
check('Every measured B2 endpoint identity/orientation',maxpair<1e-9,dict(pairs=len(unique),maximum_absolute_error=maxpair))
# Predictor contrasts must match their saved endpoint values, not synthetic labels.
ep=pd.read_csv(RUN/'pair_endpoint_predictions.csv');lut=ep.set_index(['method','seed','record_id']).prediction
maximum_real_difference_error=0.; maximum_binary32_csv_error=0.
for q in pp[pp.method.isin(['activity_gnn','pair_gnn','pair_no_message'])].itertuples():
 b=float(lut.loc[(q.method,q.seed,q.changed_id)]);a=float(lut.loc[(q.method,q.seed,q.reference_id)])
 rounded=np.float32(np.float32(b)-np.float32(a))
 assert np.float32(q.prediction)==rounded,(q.method,q.seed,q.pair_id,'binary32 subtraction mismatch')
 error=abs((b-a)-q.prediction);bound=.5*float(np.spacing(abs(rounded)))+1e-15
 assert error<=bound,(q.method,q.seed,q.pair_id,error,bound)
 maximum_real_difference_error=max(maximum_real_difference_error,error);maximum_binary32_csv_error=max(maximum_binary32_csv_error,abs(float(rounded)-q.prediction))
check('B2 predicted differences reproduce binary32 endpoint subtraction',True,dict(maximum_float64_recomputation_error=maximum_real_difference_error,maximum_binary32_csv_error=maximum_binary32_csv_error,check='Exact binary32 result equality and half-ULP rounding bound plus 1e-15 CSV parse allowance. The initial 1e-10 real-arithmetic assertion failed and is retained; no prediction changed.'))
br=pd.read_csv(RUN/'b3_primary/interaction_predictions.csv');finite=br[br.seed.astype(str)!='ensemble'];delta=finite.prediction_11-finite.prediction_10-finite.prediction_01+finite.prediction_00
check('B3 four-forward identity',np.allclose(delta,finite.prediction,atol=1e-10))
check('B3 source and 13mer held out',all(x['source_overlap']==x['sequence13mer_overlap']==0 for x in json.loads((RUN/'b3_primary/source_holdout_checks.json').read_text())))
# Check provenance links and rounded cells against the exact generator exports.
for name in ['figure_provenance.json','table_provenance.json']:
 for item in json.loads((RUN/name).read_text()):
  for source in item.get('sources',item.get('source',[])):assert sha(RUN/source['path'])==source['sha256'],source
check('Every figure/table source hash matches',True)
pa=pd.read_csv(RUN/'tables/selected_parameter_activity_summary.csv');expected_neural=sum(r['executed_updates']>0 and r['status']=='completed' for r in final);check('Every completed final neural checkpoint has a realized parameter audit',len(pa)==expected_neural,dict(expected=expected_neural,actual=len(pa)))
check('Owned scalar counts partition into zero/nonzero activity gradients',bool((pa.parameters==pa.selected_activity_gradient_nonzero_scalars+pa.selected_activity_gradient_zero_scalars).all()))
check('Actual-objective and changed-scalar counts respect owned dimensions',bool(((pa.selected_training_objective_nonzero_scalars<=pa.parameters)&(pa.changed_from_initialized_scalars<=pa.parameters)).all()))
upper=pd.read_csv(RUN/'tables/structural_parameter_counts.csv');joined=pa.merge(upper[['fit','train_reachable_upper_bound']],on='fit',validate='one_to_one');check('Realized nonzero gradients respect structural reachability upper bounds',bool((joined.selected_training_objective_nonzero_scalars<=joined.train_reachable_upper_bound).all()))
write(RUN/'final_scientific_checks.json',dict(checks=checks,fit_artifacts_hashed=filecount,fit_count=len(fits),failures=sum(r['status']!='completed' for r in fits),final_fits=len(final),development_fits=len(fits)-len(final),neural_fits=sum(r['executed_updates']>0 for r in fits),selected_updates=sum(r['selected_updates'] for r in fits),executed_updates=sum(r['executed_updates'] for r in fits),scope='Independent export and saved-record checks; no new fitting or checkpoint inference.'))
(D/'scientific_checks.md').write_text('# Final scientific checks\n\nAll '+str(len(checks))+' recorded checks passed. Exact fit matrix, training-only grouping, all ten final seed blocks, optimizer histories, artifact hashes, full activity metrics, tied constants, B2 endpoint differences and B3 four-forward identities were checked. Figure/table source hashes match their full-precision exports. This verifies the recorded computation, not biological independence or transport.\n\nDetailed checks: `'+str((RUN/'final_scientific_checks.json').relative_to(ROOT))+'`.\n')
print('Passed',len(checks),'checks;',len(fits),'fits;',filecount,'hashed fit artifacts')
