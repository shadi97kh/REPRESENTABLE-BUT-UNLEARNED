"""Prespecified retrospective evaluation of frozen predictions against measured outcomes."""
import argparse,collections,datetime,hashlib,itertools,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr
from common import read_jsonl,write_json
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
commit=json.loads((r/'prediction_commit.json').read_text());assert not (r/'evaluation_results.json').exists(),'Evaluation already frozen'
for name,h in commit['files'].items():assert hashlib.sha256((r/name).read_bytes()).hexdigest()==h,name
assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==commit['evaluation_code_sha256']
start=datetime.datetime.now(datetime.timezone.utc).isoformat();pred=pd.read_csv(r/'heldout_predictions.csv');seedpred=pd.read_csv(r/'heldout_predictions_by_seed.csv');obs=read_jsonl(r/'test_internal_observations.jsonl')+read_jsonl(r/'test_APP_observations.jsonl');byid={x['record_id']:x for x in obs};models=sorted(pred.model.unique());reference=json.loads((r/'frozen_selection.json').read_text())['development_selected_strongest_baseline'];rows=[];poolrows=[];grouprows=[];seedrows=[]
def rho(y,p):return float(spearmanr(y,p).statistic) if np.ptp(y)>0 and np.ptp(p)>0 and len(y)>2 else None
def order_accuracy(y,p):
 ii,jj=np.triu_indices(len(y),1);dy=y[ii]-y[jj];dp=p[ii]-p[jj];use=np.abs(dy)>1e-12
 return float(np.mean((np.sign(dy[use])*np.sign(dp[use])+1)/2)) if use.any() else None
for split in ['test_internal','test_APP']:
 ids=sorted(x['record_id'] for x in obs if x['split']==split);y=np.array([byid[i]['activity'] for i in ids]);studies=np.array([byid[i]['study_group'] for i in ids]);assays=np.array([byid[i]['assay_id'] or byid[i]['study_group'] for i in ids]);seq=np.array([byid[i]['sequence_group'] for i in ids]);groups=sorted(set(seq));losses={};group_losses={}
 for model in models:
  pp=pred[(pred.split==split)&(pred.model==model)].set_index('record_id').loc[ids].prediction.to_numpy();e=(pp-y)**2;losses[model]=e
  for group in sorted(set(studies)):
   q=studies==group;grouprows.append({'split':split,'model':model,'study_group':group,'n':int(q.sum()),'mse':float(e[q].mean()),'spearman':rho(y[q],pp[q])})
  prs=[]
  for assay in sorted(set(assays)):
   q=np.flatnonzero(assays==assay);verified=all(byid[ids[i]]['source_table_verified'] for i in q)
   pr={'split':split,'model':model,'assay_id':assay,'n':len(q),'primary_verified':verified,'mse':float(e[q].mean()),'spearman':rho(y[q],pp[q]),'pairwise_order_accuracy':order_accuracy(y[q],pp[q])}
   if verified and len(q)>=20:
    order=sorted(q,key=lambda i:(-pp[i],ids[i]));chosen=order[:5];u=(rankdata(y[q],method='average')-1)/(len(q)-1);local={j:i for i,j in enumerate(q)};pr.update(top5_percentile=float(np.mean([u[local[i]] for i in chosen])),selected_activity=float(y[chosen].mean()),selection_regret=float(np.sort(y[q])[-5:].mean()-y[chosen].mean()),selected_ids='|'.join(ids[i] for i in chosen));prs.append(pr)
   poolrows.append(pr)
  result={'split':split,'model':model,'n':len(y),'mse':float(e.mean()),'mae':float(np.abs(pp-y).mean()),'spearman':rho(y,pp),'equal_assay_mse':float(np.mean([e[assays==a].mean() for a in sorted(set(assays))])),'equal_assay_spearman':float(np.mean([rho(y[assays==a],pp[assays==a]) for a in sorted(set(assays)) if rho(y[assays==a],pp[assays==a]) is not None])),'ranking_pools':len(prs),'top5_percentile':float(np.mean([v['top5_percentile'] for v in prs])) if prs else None,'selected_activity':float(np.mean([v['selected_activity'] for v in prs])) if prs else None,'selection_regret':float(np.mean([v['selection_regret'] for v in prs])) if prs else None,'prediction_min':float(pp.min()),'prediction_max':float(pp.max()),'nonfinite_predictions':0}
  rows.append(result);group_losses[model]=np.array([[e[seq==g].sum(),sum(seq==g)] for g in groups])
  for seed in sorted(seedpred[seedpred.model==model].seed.unique()):
   sp=seedpred[(seedpred.split==split)&(seedpred.model==model)&(seedpred.seed==seed)].set_index('record_id').loc[ids].prediction.to_numpy();seedrows.append({'split':split,'model':model,'seed':int(seed),'mse':float(np.mean((sp-y)**2)),'spearman':rho(y,sp)})
 # Conditional paired sequence-cluster resampling. Same draws for all models.
 rng=np.random.default_rng(30260913);draws=rng.integers(0,len(groups),(2000,len(groups)));boot={m:group_losses[m][draws,0].sum(1)/group_losses[m][draws,1].sum(1) for m in models}
 for result in rows:
  if result['split']!=split:continue
  model=result['model'];result['mse_sequence_bootstrap_95']=np.quantile(boot[model],[.025,.975]).tolist();result['mse_difference_vs_reference']=float(losses[model].mean()-losses[reference].mean());result['paired_mse_difference_95']=np.quantile(boot[model]-boot[reference],[.025,.975]).tolist();result['bootstrap_sequence_components']=len(groups)
 # All paired model differences allow explicit comparison with strongest observed held-out control without retuning.
 pairout=[]
 for a,b in itertools.combinations(models,2):pairout.append({'model_a':a,'model_b':b,'mse_a_minus_b':float(losses[a].mean()-losses[b].mean()),'paired_95':np.quantile(boot[a]-boot[b],[.025,.975]).tolist()})
 write_json(r/f'{split}_paired_comparisons.json',pairout)
pairs=read_jsonl(r/'eligible_pairs.jsonl');effects=[];effectmetrics=[]
for model in models:
 pp=pred[pred.model==model].set_index('record_id').prediction.to_dict();mm=[]
 for pair in pairs:
  a,b=pair['reference_id'],pair['changed_id'];assert byid[a]['split']==byid[b]['split']=='test_APP';v=pp[b]-pp[a];z={**pair,'model':model,'predicted_difference':float(v),'sequence_group':byid[a]['sequence_group'],'difference_error':float(v-pair['observed_difference'])};effects.append(z);mm.append(z)
 if mm:
  y=np.array([x['observed_difference'] for x in mm]);p=np.array([x['predicted_difference'] for x in mm]);g=np.array([x['sequence_group'] for x in mm]);groups=sorted(set(g));e=(p-y)**2;active=np.abs(y)>.02;predsign=np.where(np.abs(p)<=.02,0,np.sign(p));sgn=(predsign==np.sign(y))
  weights=np.array([1/sum(g==a) for a in g]);weights/=weights.sum();gm=np.array([e[g==a].mean() for a in groups]);rng=np.random.default_rng(30260913);boot=gm[rng.integers(0,len(groups),(2000,len(groups)))].mean(1)
  coef=np.polyfit(p,y,1).tolist() if np.ptp(p)>1e-12 else None
  effectmetrics.append({'model':model,'pairs':len(mm),'exact_sequence_backgrounds':len(set(x['background_id'] for x in mm)),'sequence_components':len(groups),'pair_mse':float(e.mean()),'equal_sequence_component_mse':float(gm.mean()),'equal_component_mse_95':np.quantile(boot,[.025,.975]).tolist(),'mae':float(np.mean(np.abs(p-y))),'spearman':rho(y,p),'observed_non_ties':int(active.sum()),'sign_accuracy_on_observed_non_ties':float(sgn[active].mean()) if active.any() else None,'observed_ties':int((~active).sum()),'predicted_ties':int((np.abs(p)<=.02).sum()),'calibration_slope_intercept':coef,'measurement_uncertainty':'Reported endpoint SDs preserved; no fabricated SE or shared-control covariance.'})
# Pairwise effect comparisons preserve the same sequence grouping and all six assay/dose appearances.
effect_comparisons=[]
if effects:
 ef=pd.DataFrame(effects);base=ef[ef.model==reference].set_index('pair_id')
 for model in models:
  a=ef[ef.model==model].set_index('pair_id').loc[base.index];g=a.sequence_group.to_numpy();diff=a.difference_error.to_numpy()**2-base.difference_error.to_numpy()**2;groups=sorted(set(g));v=np.array([diff[g==k].mean() for k in groups]);rng=np.random.default_rng(30260913);boot=v[rng.integers(0,len(v),(2000,len(v)))].mean(1);effect_comparisons.append({'model':model,'reference':reference,'equal_component_mse_difference':float(v.mean()),'paired_95':np.quantile(boot,[.025,.975]).tolist(),'sequence_components':len(v)})
 pd.DataFrame(effects).drop(columns=['changes']).to_csv(r/'chemical_effect_predictions.csv',index=False)
for name,values in [('performance.csv',rows),('assay_metrics.csv',poolrows),('study_metrics.csv',grouprows),('seed_metrics.csv',seedrows),('chemical_effect_metrics.csv',effectmetrics)]:pd.DataFrame(values).to_csv(r/name,index=False)
write_json(r/'chemical_effect_paired_comparisons.json',effect_comparisons)
rectangles=read_jsonl(r/'eligible_rectangles.jsonl');assert len(rectangles)==0,'B3 scorer needs explicit implementation if eligible rectangle count changes'
result={'evaluation_started_utc':start,'prediction_commit_utc':commit['committed_utc'],'reference_selected_on_validation':reference,'performance':rows,'chemical_effects':effectmetrics,'chemical_effect_comparisons':effect_comparisons,'B3':{'eligible_rectangles':0,'status':'incomplete: no complete two-single-position measured factorial rectangle in admitted matched blocks'},'uncertainty_scope':'2000 paired sequence-component bootstrap draws conditional on these observed source families. APP has one related source-family group; intervals do not estimate new-study or new-laboratory variation. Two internal held-out study components and only three validation sequence components limit external validity. Ranking endpoints have per-pool evidence and descriptive aggregates; no independent-assay confidence interval is asserted.','ranking_scope':'Primary-table-defined target/cell/dose/time pools, all n>=20. Different doses and related studies are not independent biological cohorts.','B2_scope':'Retrospective differences of primary-table means for observed matched chemical patterns; no claimed causal mechanism or calibrated population effect.','wet_lab_validation_performed':False,'nonfinite_or_failed_prediction_rows':0,'all_negative_results_preserved':True,'strongest_test_baseline_posthoc':min((x for x in rows if x['split']=='test_APP' and x['model'] not in ['gnn','gnn_no_chemistry']),key=lambda x:x['mse'])['model']}
write_json(r/'evaluation_results.json',result);write_json(r/'evaluate.outputs.json',['evaluation_results.json','performance.csv','assay_metrics.csv','study_metrics.csv','seed_metrics.csv','chemical_effect_metrics.csv','chemical_effect_predictions.csv','chemical_effect_paired_comparisons.json','test_internal_paired_comparisons.json','test_APP_paired_comparisons.json']);print(json.dumps(result,indent=2))
