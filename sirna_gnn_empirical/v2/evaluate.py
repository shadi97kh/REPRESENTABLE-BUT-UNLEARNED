from core import *
from scipy.stats import spearmanr,rankdata
rows,_=load();meta=pd.DataFrame([{k:x.get(k) for k in ['record_id','activity','split','study_group','sequence_group','assay_id','source_family']} for x in rows]);pairs=readlines(OLD/'eligible_pairs.jsonl');lookup={x['record_id']:x for x in rows}
# These hashes commit completed predictions before any new metric table is computed.
write(RUN/'prediction_commit.json',dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),predictions={p.name:sha(p) for p in RUN.glob('*predictions.csv')},evaluation_code_sha256=sha(__file__)))
def ensemble(df,keys,target):
 a=df.copy();a['seed']=a.seed.astype(str);en=a.groupby(keys,as_index=False)[target].mean();en['seed']='ensemble';return pd.concat([a,en],ignore_index=True)
def rho(y,p):return float(spearmanr(y,p).statistic) if np.ptp(y)>0 and np.ptp(p)>1e-12 else np.nan
def effect_metric(df):
 y=df.observed_difference.to_numpy();p=df.predicted_difference.to_numpy();e=p-y;w=weights(df.sequence_group);o=np.sign(y)*(abs(y)>.02);d=np.sign(p)*(abs(p)>.02);non=o!=0;call=non&(d!=0);correct=call&(o==d)
 return dict(n=len(df),sequence_components=df.sequence_group.nunique(),pair_mse=np.mean(e**2),equal_component_mse=np.average(e**2,weights=w),pair_mae=np.mean(abs(e)),equal_component_mae=np.average(abs(e),weights=w),mean_observed=np.mean(y),mean_predicted=np.mean(p),predicted_min=min(p),predicted_max=max(p),predicted_mean_absolute=np.mean(abs(p)),observed_negative=sum(o<0),observed_tie=sum(o==0),observed_positive=sum(o>0),predicted_negative=sum(d<0),predicted_tie=sum(d==0),predicted_positive=sum(d>0),raw_positive=sum(p>0),observed_nonties=sum(non),calls=sum(call),correct_calls=sum(correct),coverage=sum(call)/sum(non) if sum(non) else np.nan,correctness_among_calls=sum(correct)/sum(call) if sum(call) else np.nan,correctness_all_nonties=sum(correct)/sum(non) if sum(non) else np.nan)
allmetrics=[];rankings=[];allpred=[];effects=[]
for phase,file in [('frozen','frozen_diagnostic_predictions.csv'),('factorial','factorial_predictions.csv'),('deployment','deployment_predictions.csv')]:
 df=pd.read_csv(RUN/file);df=ensemble(df,['record_id','split','model','arm'],'prediction');df['phase']=phase;allpred.append(df)
 for (kind,arm,seed),part in df.groupby(['model','arm','seed']):
  record=dict(phase=phase,model=kind,arm=arm,seed=seed);pp=part.merge(meta.drop(columns='split'),on='record_id',validate='many_to_one')
  for split,g in pp.groupby('split'):
   e=g.prediction-g.activity;prefix='APP' if split=='test_APP' else 'internal';record[prefix+'_mse']=np.mean(e**2);record[prefix+'_mae']=np.mean(abs(e));record[prefix+'_rho']=rho(g.activity,g.prediction)
   if split=='test_APP':
    rr=[]
    for aid,a in g.groupby('assay_id'):
     a=a.copy();a['percentile']=(rankdata(a.activity,method='average')-1)/(len(a)-1);sel=a.sort_values(['prediction','record_id'],ascending=[False,True]).head(5)
     z=dict(phase=phase,model=kind,arm=arm,seed=seed,assay_id=aid,n=len(a),spearman=rho(a.activity,a.prediction),top5_percentile=sel.percentile.mean(),selected_activity=sel.activity.mean(),oracle_top5_activity=a.nlargest(5,'activity').activity.mean(),pool_mean_activity=a.activity.mean());rankings.append(z);rr.append(z)
    record['APP_equal_assay_rho']=np.nanmean([z['spearman'] for z in rr]) if any(np.isfinite(z['spearman']) for z in rr) else np.nan;record['APP_top5_percentile']=np.mean([z['top5_percentile'] for z in rr]);record['APP_selected_activity']=np.mean([z['selected_activity'] for z in rr])
  pmap=dict(zip(part.record_id,part.prediction));ee=[]
  for p in pairs:
   if p['reference_id'] in pmap:
    ee.append(dict(pair_id=p['pair_id'],sequence_group=lookup[p['reference_id']]['sequence_group'],observed_difference=p['observed_difference'],predicted_difference=pmap[p['changed_id']]-pmap[p['reference_id']]))
  if ee:
   metrics=effect_metric(pd.DataFrame(ee));record['B2_pair_mse']=metrics['pair_mse'];record['B2_equal_component_mse']=metrics['equal_component_mse']
   for z in ee:effects.append(dict(phase=phase,model=kind,arm=arm,seed=seed,**z))
  allmetrics.append(record)
pd.concat(allpred).to_csv(RUN/'all_B1_predictions.csv',index=False);pd.DataFrame(allmetrics).to_csv(RUN/'B1_comparison.csv',index=False);pd.DataFrame(rankings).to_csv(RUN/'ranking_by_assay.csv',index=False);pd.DataFrame(effects).to_csv(RUN/'B2_transfer_predictions.csv',index=False)
# Add exact fit accounting to each per-seed/ensemble B1 row, retaining frozen fits as archived counts.
fits=[]
for file in ['factorial_fits.json','grouped_fits.json','deployment_fits.json','b2_fits.json']:
 fits+=json.loads((RUN/file).read_text())
fitdf=pd.DataFrame([{k:v for k,v in f.items() if k not in ['history']} for f in fits]);fitdf.to_csv(RUN/'all_fit_summary.csv',index=False)
summary=pd.read_csv(RUN/'B1_comparison.csv');oldfits=json.loads((OLD/'final_fit_results.json').read_text())
for i,row in summary.iterrows():
 if row.phase=='frozen':ff=[f for f in oldfits if f['kind']==row.model]
 else:ff=[f for f in fits if f['kind']==row.model and ((row.phase=='factorial' and f['tag']==f'factorial/final/{row.arm}') or (row.phase=='deployment' and f['tag']=='deployment'))]
 if str(row.seed)!='ensemble':ff=[f for f in ff if str(f['seed'])==str(row.seed)]
 summary.loc[i,'parameters_each']=ff[0].get('parameters',np.nan) if ff else (1 if row.model.startswith('mean_') else np.nan);summary.loc[i,'selected_epochs']='/'.join(str(f.get('selected_epoch',0)) for f in ff);summary.loc[i,'executed_updates']=sum(f.get('optimizer_updates',0) for f in ff);summary.loc[i,'selected_updates']=sum(f.get('selected_updates',f.get('selected_checkpoint_optimizer_updates',0)) for f in ff)
summary.to_csv(RUN/'B1_comparison.csv',index=False)
# Grouped outer predictions and exact component results; no outer label drives any fit.
group=ensemble(pd.read_csv(RUN/'grouped_predictions.csv'),['record_id','outer_fold','model'],'prediction').merge(meta,on='record_id');group.to_csv(RUN/'grouped_predictions_with_labels.csv',index=False);gm=[];detail=[]
for (kind,seed),g in group.groupby(['model','seed']):
 e=g.prediction-g.activity;comp=g.assign(squared_error=e**2,absolute_error=abs(e));gm.append(dict(model=kind,seed=seed,n=len(g),components=g.study_group.nunique(),pooled_mse=np.mean(e**2),equal_component_mse=comp.groupby('study_group').squared_error.mean().mean(),equal_fold_mse=comp.groupby('outer_fold').squared_error.mean().mean(),mae=np.mean(abs(e))))
 for (fold,study),a in comp.groupby(['outer_fold','study_group']):detail.append(dict(model=kind,seed=seed,outer_fold=fold,study_group=study,n=len(a),source_labels=';'.join(sorted(a.source_family.unique())),mse=a.squared_error.mean(),mae=a.absolute_error.mean()))
pd.DataFrame(gm).to_csv(RUN/'grouped_comparison.csv',index=False);pd.DataFrame(detail).to_csv(RUN/'grouped_component_results.csv',index=False)
b=ensemble(pd.read_csv(RUN/'b2_predictions.csv'),['outer_fold','model','pair_id','reference_id','changed_id','sequence_group','observed_difference','measurement_sd_reference','measurement_sd_changed'],'predicted_difference');b.to_csv(RUN/'B2_supervised_predictions.csv',index=False);bm=[];bd=[]
for (kind,seed),g in b.groupby(['model','seed']):
 z=effect_metric(g)
 if kind=='training_majority':
  for key in ['pair_mse','equal_component_mse','pair_mae','equal_component_mae','mean_predicted','predicted_min','predicted_max','predicted_mean_absolute']:z[key]=np.nan
 bm.append(dict(model=kind,seed=seed,**z))
 for fold,a in g.groupby('outer_fold'):bd.append(dict(model=kind,seed=seed,outer_fold=fold,**effect_metric(a)))
pd.DataFrame(bm).to_csv(RUN/'B2_supervised_comparison.csv',index=False);pd.DataFrame(bd).to_csv(RUN/'B2_fold_results.csv',index=False)
# Conditional paired component bootstrap, identical draws within each comparison. Descriptive post-design uncertainty.
boot=[];rng=np.random.default_rng(30260914)
def compare(df,idcol,groupcol,ycol,pcol,modelcol,reference,phase,weighting):
 table=df.pivot(index=idcol,columns=modelcol,values=pcol);md=df.drop_duplicates(idcol).set_index(idcol).loc[table.index];y=md[ycol];groups=md[groupcol];names=sorted(groups.unique());draws=rng.integers(0,len(names),size=(2000,len(names)))
 for name in table.columns:
  if name==reference:continue
  diff=(table[name]-y)**2-(table[reference]-y)**2;num=np.array([diff[groups==g].sum() for g in names]);den=np.array([(groups==g).sum() for g in names]);vv=(num[draws].sum(1)/den[draws].sum(1)) if weighting=='rows' else (num/den)[draws].mean(1);estimate=diff.mean() if weighting=='rows' else np.mean(num/den);lo,hi=np.quantile(vv,[.025,.975]);boot.append(dict(phase=phase,model=name,reference=reference,weighting=weighting,difference_mse=estimate,lower95=lo,upper95=hi,groups=len(names),bootstrap_replicates=2000))
df=pd.concat(allpred);ens=df[(df.seed=='ensemble')&(df['split']=='test_APP')].merge(meta.drop(columns='split'),on='record_id');ens['modelarm']=ens.model+'_'+ens.arm
for phase in ['factorial','deployment']:
 q=ens[ens.phase==phase];ref='chemistry_tree_R1C1' if phase=='factorial' else 'chemistry_tree_all_ENsi_inner_selected';compare(q,'record_id','sequence_group','activity','prediction','modelarm',ref,phase,'rows')
compare(group[group.seed=='ensemble'],'record_id','study_group','activity','prediction','model','chemistry_tree','grouped','components')
compare(b[(b.seed=='ensemble')&(b.model!='training_majority')],'pair_id','sequence_group','observed_difference','predicted_difference','model','training_mean','B2_supervised','components');pd.DataFrame(boot).to_csv(RUN/'paired_comparisons.csv',index=False)
print(summary[summary.seed=='ensemble'][['phase','model','arm','internal_mse','APP_mse','APP_top5_percentile','B2_pair_mse']].to_string(index=False));print(pd.DataFrame(gm).query("seed=='ensemble'").to_string(index=False));print(pd.DataFrame(bm).query("seed=='ensemble'").to_string(index=False))
