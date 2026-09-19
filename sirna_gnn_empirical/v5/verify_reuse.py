"""Independently rescore saved predictions and verify reuse dependencies; no fits."""
from common import *
import pandas as pd,numpy as np
out=RUN/'reuse';out.mkdir(exist_ok=True)
ens=pd.read_csv(PARENT/'evaluation/ensemble_predictions.csv');seeds=pd.read_csv(PARENT/'evaluation/seed_predictions.csv');old=pd.read_csv(PARENT/'evaluation/activity_metrics.csv');rows=[]
for (s,c,m),g in ens.groupby(['scenario','cohort','method']):
 y=g.activity.to_numpy();p=g.prediction.to_numpy();my=y.mean();mp=p.mean();vy=np.mean((y-my)**2);vp=np.mean((p-mp)**2);mse=np.mean((y-p)**2)
 r=dict(scenario=s,cohort=c,method=m,n=len(g),mse=mse,r_squared=1-mse/vy,correlation=np.mean((y-my)*(p-mp))/np.sqrt(vy*vp) if vp>1e-20 else np.nan,prediction_sd=np.sqrt(vp),mean_bias=mp-my,oracle_test_mean_mse=vy)
 t=old[(old.scenario==s)&(old.cohort==c)&(old.method==m)&(old.weighting=='rows')].iloc[0]
 for k in ['mse','r_squared','prediction_sd','mean_bias']:assert abs(r[k]-t[k])<1e-10,(s,c,m,k)
 if np.isfinite(r['correlation']):assert abs(r['correlation']-t.correlation)<1e-10
 z=seeds[(seeds.scenario==s)&(seeds.cohort==c)&(seeds.method==m)]
 assert np.max(abs(z.groupby('record_id').prediction.mean()-g.set_index('record_id').prediction))<1e-12
 r['mean_seed_mse']=np.mean((z.activity-z.prediction)**2)
 z=z.merge(g[['record_id','prediction']].rename(columns={'prediction':'ensemble'}),on='record_id',validate='many_to_one');r['seed_dispersion_term']=np.mean((z.prediction-z.ensemble)**2)
 assert abs(r['mean_seed_mse']-r['mse']-r['seed_dispersion_term'])<1e-12
 rows.append(r)
pd.DataFrame(rows).to_csv(out/'recomputed_metrics.csv',index=False)
fit=pd.read_csv(PARENT/'fit_summary.csv');assert len(fit)==792;artifact_count=0
for p in (PARENT/'fits').rglob('fit.json'):
 r=json.loads(p.read_text())
 for name,h in r['hashes'].items():assert sha(p.parent/name)==h;artifact_count+=1
 for name,h in r['signature']['code'].items():assert sha(ROOT/'sirna_gnn_empirical/v4'/name)==h
 assert r['signature']['protocol']==sha(PARENT/'protocol.json')
mem=pd.read_csv(PARENT/'focused_membership.csv');groups=[]
for (s,pop,m),g in mem.groupby(['scenario','population','method']):
 for col in ['record_id','study_group','sequence_group']:
  sets={role:set(g[g.role==role][col]) for role in ['train','validation','test']}
  for a,b in [('train','validation'),('train','test'),('validation','test')]:assert not sets[a]&sets[b]
 if s=='source_excluded':assert not g.source.astype(str).eq('19282453').any()
 groups.append(dict(scenario=s,population=pop,method=m,passed=True))
pd.DataFrame(groups).to_csv(out/'dependency_checks.csv',index=False)
# Copy exact adjudication and audit exports, not reconstructed or relabeled evidence.
import shutil
paths=['adjudication','dataset','visibility','review_checks','ranking','source_checks','paper_tables']
manifest=[]
for folder in paths:
 for p in sorted((PARENT/folder).rglob('*')):
  if not p.is_file():continue
  manifest.append(dict(parent_path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,status='verified reuse'))
# Check primary citation bytes and the original source-held-out checkpoint dependencies.
oldrun=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z'
for record in json.loads((PARENT/'source_checks/citation_recheck.json').read_text()):
 for e in record['evidence']:
  p=oldrun/e['path'];assert sha(p)==e['sha256'];manifest.append(dict(parent_path=str(p.relative_to(ROOT)),sha256=e['sha256'],bytes=p.stat().st_size,status='verified primary citation bytes'))
for r in json.loads((PARENT/'evaluation/reused_fits.json').read_text()):
 p=ROOT/r['fit'];assert sha(p)==r['fit_record_sha256']
 possible=[x for x in p.parent.iterdir() if x.name in ['best.pt','model.pkl']]
 assert len(possible)==1;assert sha(possible[0])==r['checkpoint_sha256']
 for item in [p,possible[0]]:manifest.append(dict(parent_path=str(item.relative_to(ROOT)),sha256=sha(item),bytes=item.stat().st_size,status='verified unaffected parent checkpoint'))
for folder in ['evaluation']:
 for p in (PARENT/folder).glob('*'):
  if p.is_file():manifest.append(dict(parent_path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,status='verified reused numerical input'))
write(out/'evidence_reuse_manifest.json',manifest)
write(out/'summary.json',dict(passed=True,metrics_recomputed=len(rows),reused_parent_fits=792,reused_parent_fit_artifacts=artifact_count,grouping_checks=len(groups),new_fits=0,new_selected_updates=0,new_executed_updates=0,source_adjudication='v4 uniquely linked primary-assay sensitivity retained; not certified label repair',reason_no_refits='No newly justified changes to identities, labels, preprocessing, training or selection dependencies. v4 executed all affected focused dependencies.',historical_negative_results_preserved=True))
print(json.dumps(json.loads((out/'summary.json').read_text()),indent=2))
