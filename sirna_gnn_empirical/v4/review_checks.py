"""Independent finite-sample review checks on saved v3 predictions; no fitting."""
from common import *
from metrics import metrics
import collections
os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';out=RUN/'review_checks';out.mkdir(exist_ok=True)
df=pd.read_csv(V3/'evaluated/activity_ensemble_predictions.csv');grouped=df[df.cohort=='ENsiRNA_grouped'];results=[]
for subset,z in [('all_sources',grouped),('Bramsen_only',grouped[grouped.source_family.astype(str)=='19282453']),('other_sources',grouped[grouped.source_family.astype(str)!='19282453'])]+[(str(s),z) for s,z in grouped.groupby('source_family')]:
 for method,g in z.groupby('method'):
  for weighting in ['rows','equal_source','equal_study_component']:
   w=None if weighting=='rows' else weights(g.source_family if weighting=='equal_source' else g.study_group)
   results.append(dict(subset=subset,method=method,weighting=weighting,sources=g.source_family.nunique(),study_components=g.study_group.nunique(),sequence_components=g.sequence_group.nunique(),**metrics(g.activity,g.prediction,w)))
export(pd.DataFrame(results),'review_checks/source_scores.csv')
base=grouped[grouped.method=='corrected_gnn'];concentration=base.groupby('source_family').agg(n=('record_id','size'),sequence_components=('sequence_group','nunique'),study_components=('study_group','nunique')).reset_index();concentration['fraction']=concentration.n/len(base);export(concentration,'review_checks/source_concentration.csv')
cal=[]
for (cohort,method),g in df.groupby(['cohort','method']):
 y=g.activity.to_numpy(float);p=g.prediction.to_numpy(float);b=p.mean()-y.mean();vy=np.mean((y-y.mean())**2);vp=np.mean((p-p.mean())**2);cov=np.mean((p-p.mean())*(y-y.mean()));mse=np.mean((p-y)**2)
 assert abs(mse-(b*b+vp+vy-2*cov))<1e-12
 cal.append(dict(cohort=cohort,method=method,n=len(g),mse=mse,squared_mean_bias=b*b,prediction_variance=vp,label_variance=vy,covariance=cov,centered_error=vp+vy-2*cov,mean_bias=b,scope='exact observed-sample identity; not an explanation of the origin of bias'))
cal=pd.DataFrame(cal);export(cal,'review_checks/calibration_decomposition.csv')
app=cal[cal.cohort=='APP'].set_index('method');gap=app.loc['chemistry_tree','mse']-app.loc['corrected_gnn','mse'];biasgap=app.loc['chemistry_tree','squared_mean_bias']-app.loc['corrected_gnn','squared_mean_bias']
# Six prespecified contrasts, independent RNG; resample components, not rows/rectangles.
comparisons=[('APP','corrected_gnn','chemistry_tree'),('APP','corrected_gnn','corrected_no_message'),('APP','corrected_gnn','original_no_message'),('Davis_S7','corrected_gnn','original_gnn'),('Davis_S7','corrected_gnn','corrected_no_message'),('ENsiRNA_grouped','corrected_gnn','chemistry_tree')];boot=[]
for cohort,a,b in comparisons:
 z=df[df.cohort==cohort];meta=z.drop_duplicates('record_id').set_index('record_id');wide=z.pivot(index='record_id',columns='method',values='prediction');loss=(wide[a]-meta.activity)**2-(wide[b]-meta.activity)**2;tmp=pd.DataFrame(dict(loss=loss,group=meta.sequence_group));gs=tmp.groupby('group').loss.agg(['sum','count']);rng=np.random.default_rng(914221);vals=[]
 for _ in range(100):
  ix=rng.integers(len(gs),size=(100,len(gs)));vals.extend(gs['sum'].to_numpy()[ix].sum(1)/gs['count'].to_numpy()[ix].sum(1))
 boot.append(dict(cohort=cohort,a=a,b=b,difference=loss.mean(),lower=np.quantile(vals,.025),upper=np.quantile(vals,.975),groups=len(gs),rng_seed=914221,resamples=10000,scope='fixed-prediction component-conditional; separate training-seed variation'))
export(pd.DataFrame(boot),'review_checks/independent_conditional_comparisons.csv')
# Actual role trace for every old fit: includes inner selection, not just final training.
rows=rows_all()+readlines(PREV/'external_observations.jsonl');trace=[]
for fp in sorted((V3/'fits').rglob('fit.json')):
 rec=json.loads(fp.read_text());sig=rec['signature'];counts={role:sum(rows[i]['source_family']=='19282453' for i in sig.get(key,[])) for role,key in [('training','train'),('validation','validation'),('test','test')]};trace.append(dict(fit=rec['name'],method=rec['method'],seed=rec['seed'],**counts,used_for_model_selection='/dev/' in rec['name'] or '_dev/' in rec['name'],source_excluded_every_role=counts['training']==counts['validation']==0,fit_record_sha256=sha(fp)))
export(pd.DataFrame(trace),'review_checks/historical_training_source_roles.csv')
summary=dict(source_rows=int(concentration[concentration.source_family.astype(str)=='19282453'].n.iloc[0]),all_rows=len(base),source_fraction=1924/2927,APP_tree_minus_GNN_mse=gap,APP_tree_minus_GNN_squared_bias=biasgap,APP_bias_fraction_of_gap=biasgap/gap,legacy_fits_traced=len(trace),B2_Bramsen_training_validation_count=sum(x['training']+x['validation'] for x in trace if x['fit'].startswith('pair/')),independent_bootstrap_seed=914221,review_handoff='Named Review_Validity_Assessment.md/reanalyze_review.py not supplied at execution time; equivalent checks implement the explicit listed review claims. No claim to have read missing reviews.',new_fits=0)
write(out/'summary.json',summary);print(json.dumps(summary,indent=2));print(pd.DataFrame(results).query("subset == 'other_sources' and weighting == 'rows' and method in ['corrected_gnn','chemistry_tree']")[['method','mse','r_squared']].to_string(index=False))
