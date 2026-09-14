"""Figures and descriptive audits from immutable measured-data predictions; no fitting."""
import argparse,collections,hashlib,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import read_jsonl,write_json,digest
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
assert (r/'evaluate.done.json').exists();out=r/'figures';out.mkdir(exist_ok=True);data=out/'data';data.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':180})
names={'guide_ridge':'Guide ridge','guide_pairwise_ridge':'Pairwise ridge','chemistry_tree':'Chemistry tree','token_cnn':'Token CNN','nongraph':'No message passing','gnn':'Chemistry GNN','gnn_no_chemistry':'GNN without chemistry'}
order=list(names);colors={m:('#006d77' if m=='gnn' else '#df7e32' if m=='chemistry_tree' else '#777f8c') for m in order};manifest=[]
def save(fig,name,sources):
 fig.tight_layout()
 for ext in ['pdf','svg','png']:fig.savefig(out/f'{name}.{ext}',bbox_inches='tight')
 plt.close(fig);manifest.append({'figure':name,'sources':sources,'files':[f'figures/{name}.{ext}' for ext in ['pdf','svg','png']]})
def copydata(name):shutil.copy2(r/name,data/name)
for name in ['performance.csv','assay_metrics.csv','seed_metrics.csv','chemical_effect_metrics.csv','chemical_effect_predictions.csv','test_APP_paired_comparisons.json','test_internal_paired_comparisons.json']:copydata(name)
e=json.loads((r/'evaluation_results.json').read_text());perf=pd.DataFrame(e['performance']);assay=pd.read_csv(r/'assay_metrics.csv');seed=pd.read_csv(r/'seed_metrics.csv');effects=pd.read_csv(r/'chemical_effect_predictions.csv');pred=pd.read_csv(r/'heldout_predictions.csv');obs=read_jsonl(r/'test_internal_observations.jsonl')+read_jsonl(r/'test_APP_observations.jsonl');lookup={x['record_id']:x for x in obs}
fig,axs=plt.subplots(1,2,figsize=(7.1,2.8),sharey=True)
for ax,sp,title in zip(axs,['test_internal','test_APP'],['Internal source holdout (198)','APP source-family holdout (1,839)']):
 vals=perf[perf.split==sp].set_index('model').loc[order];x=vals.mse.to_numpy();ci=np.array(vals.mse_sequence_bootstrap_95.tolist());ax.errorbar(x,np.arange(7),xerr=np.maximum(0,np.stack([x-ci[:,0],ci[:,1]-x])),fmt='none',ecolor='#a3a3a3',capsize=2)
 for i,m in enumerate(order):ax.scatter(x[i],i,c=colors[m],s=30)
 ax.set_yticks(range(7),[names[m] for m in order]);ax.set_title(title);ax.set_xlabel('Measured activity MSE');ax.grid(axis='x',alpha=.15)
axs[0].invert_yaxis();save(fig,'main_performance',['performance.csv'])
fig,axs=plt.subplots(1,2,figsize=(7.1,2.65))
a=perf[perf.split=='test_APP'].set_index('model').loc[order]
for i,m in enumerate(order):axs[0].barh(i,a.loc[m,'top5_percentile'],color=colors[m]);axs[1].barh(i,a.loc[m,'selection_regret'],color=colors[m])
for ax in axs:ax.set_yticks(range(7),[names[m] for m in order]);ax.invert_yaxis()
axs[0].axvline(.5,ls=':',c='black',lw=1);axs[0].set_xlim(0,1);axs[0].set_xlabel('Top-five measured midrank percentile');axs[0].set_title('17 primary-verified assay pools')
axs[1].set_xlabel('Top-five activity regret');axs[1].set_title('Equal-pool descriptive averages');save(fig,'main_ranking',['performance.csv','assay_metrics.csv'])
fig,axs=plt.subplots(1,2,figsize=(7.1,2.65))
x=effects[effects.model=='gnn'];axs[0].scatter(x.observed_difference,x.predicted_difference,s=12,alpha=.5,color=colors['gnn']);lo=min(x.observed_difference.min(),-.1);hi=max(x.observed_difference.max(),.1);axs[0].plot([lo,hi],[lo,hi],ls='--',c='gray',lw=1);axs[0].axhline(0,c='gray',lw=.5);axs[0].set_xlabel('Measured chemistry difference');axs[0].set_ylabel('GNN predicted difference');axs[0].set_title('156 pairs / 26 backgrounds / 12 clusters')
comps=json.loads((r/'test_APP_paired_comparisons.json').read_text());v=[]
for m in order:
 if m=='gnn':continue
 q=next(z for z in comps if {z['model_a'],z['model_b']}=={m,'gnn'});sgn=1 if q['model_a']=='gnn' else -1;ci=np.sort(np.array(q['paired_95'])*sgn);v.append((m,q['mse_a_minus_b']*sgn,ci))
for i,(m,z,ci) in enumerate(v):axs[1].errorbar(z,i,xerr=[[max(0,z-ci[0])],[max(0,ci[1]-z)]],fmt='o',c=colors[m],capsize=2)
axs[1].axvline(0,c='gray',ls='--',lw=1);axs[1].set_yticks(range(len(v)),[names[z[0]] for z in v]);axs[1].invert_yaxis();axs[1].set_xlabel('APP MSE: GNN minus comparator');axs[1].set_title('Paired conditional 95% intervals');save(fig,'main_effects_ablations',['chemical_effect_predictions.csv','test_APP_paired_comparisons.json'])
# Seven additional full diagnostics from saved evidence.
fig,axs=plt.subplots(2,4,figsize=(10,5),sharex=True,sharey=True)
for ax,m in zip(axs.flat,order):
 z=pred[(pred.model==m)&(pred.split=='test_APP')];yy=np.array([lookup[i]['activity'] for i in z.record_id]);ax.scatter(yy,z.prediction,s=3,alpha=.25,c=colors[m]);ax.plot([-2,1.4],[-2,1.4],c='gray',ls='--',lw=.7);ax.set_title(names[m]);ax.set_xlabel('Measured inhibition');ax.set_ylabel('Prediction')
axs.flat[-1].axis('off');save(fig,'appendix_prediction_scatter',['heldout_predictions.csv','test_APP_observations.jsonl'])
fig,axs=plt.subplots(1,2,figsize=(9,3))
for ax,sp in zip(axs,['test_internal','test_APP']):
 for i,m in enumerate(order):
  z=seed[(seed.model==m)&(seed.split==sp)];ax.scatter(np.arange(len(z))*.08+i-.08,z.mse,c=colors[m],s=25)
 ax.set_xticks(range(7),[names[m] for m in order],rotation=35,ha='right');ax.set_ylabel('Per-seed MSE');ax.set_title(sp.replace('test_','')+'; seeds are not biological replicates')
save(fig,'appendix_seed_variation',['seed_metrics.csv'])
q=assay[(assay.split=='test_APP')].pivot(index='assay_id',columns='model',values='top5_percentile').loc[:,order];q.to_csv(data/'assay_top5_matrix.csv');fig,ax=plt.subplots(figsize=(9,6));im=ax.imshow(q.to_numpy(),vmin=0,vmax=1,cmap='viridis',aspect='auto');ax.set_xticks(range(7),[names[m] for m in order],rotation=30,ha='right');ax.set_yticks(range(len(q)),[s.replace('US20220002724A1','2020 family').replace('US20240117349A1','2022 family').replace('Primary ','').replace(' cell line','').replace('|24h','') for s in q.index],fontsize=7);fig.colorbar(im,ax=ax,label='Selected top-five measured percentile');save(fig,'appendix_assay_ranking',['assay_metrics.csv'])
fig,axs=plt.subplots(2,4,figsize=(10,5),sharex=True,sharey=True)
for ax,m in zip(axs.flat,order):
 z=effects[effects.model==m];ax.scatter(z.observed_difference,z.predicted_difference,s=8,alpha=.5,c=colors[m]);ax.plot([lo,hi],[lo,hi],c='gray',ls='--',lw=.7);ax.set_title(names[m]);ax.set_xlabel('Measured difference');ax.set_ylabel('Predicted difference')
axs.flat[-1].axis('off');save(fig,'appendix_all_effects',['chemical_effect_predictions.csv'])
fits=json.loads((r/'final_fit_results.json').read_text());fits=fits['fits'] if isinstance(fits,dict) else fits
hist=[];fig,axs=plt.subplots(1,4,figsize=(10,2.6))
for ax,m in zip(axs,['token_cnn','nongraph','gnn','gnn_no_chemistry']):
 for f in fits:
  if f['kind']!=m:continue
  ax.plot([h['epoch'] for h in f['history']],[h['validation_equal_study_mse'] for h in f['history']],label=str(f['seed']));hist.extend([{'model':m,'seed':f['seed'],**h} for h in f['history']])
 ax.set_title(names[m],fontsize=8);ax.set_xlabel('Executed epoch');ax.set_ylabel('Validation MSE');ax.legend(fontsize=6)
pd.DataFrame(hist).to_csv(data/'final_training_histories.csv',index=False);save(fig,'appendix_learning_curves',['final_fit_results.json'])
graph=json.loads((r/'graph_manifest.json').read_text());support=[{'split':k,**v} for k,v in graph['statistics'].items()];pd.DataFrame(support).to_csv(data/'feature_support.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(9,3));labels=[x['split'] for x in support];axs[0].bar(labels,[x['observations'] for x in support],color='#777f8c');axs[0].set_ylabel('Admitted observations');axs[1].bar(labels,[x['rows_with_unseen_modification_names']/x['observations'] for x in support],color='#df7e32');axs[1].set_ylabel('Fraction with chemistry name absent in training');axs[1].set_ylim(0,1.05)
for ax in axs:ax.tick_params(axis='x',rotation=20)
save(fig,'appendix_data_support',['graph_manifest.json'])
# Diagnostic context shifts: no model transformation is changed after evaluation.
context=[]
for sp in ['train','validation','test_internal','test_APP']:
 z=np.load(r/f'{sp}_graphs.npz')['C'];context.extend({'split':sp,'feature':n,'minimum':float(z[:,j].min()),'maximum':float(z[:,j].max()),'mean':float(z[:,j].mean())} for j,n in enumerate(graph['context_features']))
pd.DataFrame(context).to_csv(data/'context_shift.csv',index=False)
fig,ax=plt.subplots(figsize=(8,3));xx=np.arange(8)
for k,sp in enumerate(['validation','test_internal','test_APP']):
 a=[z for z in context if z['split']==sp];ax.plot(xx,[max(abs(z['minimum']),abs(z['maximum'])) for z in a],marker='o',label=sp)
ax.set_yscale('symlog',linthresh=1);ax.set_xticks(xx,graph['context_features'],rotation=30,ha='right',fontsize=7);ax.set_ylabel('Maximum absolute standardized input');ax.legend();save(fig,'appendix_context_shift',['graph_manifest.json','*_graphs.npz'])
# Confirm B3 absence is not merely an implementation's single-edit restriction.
blocks=collections.defaultdict(set)
for x in obs:
 if x['split']!='test_APP':continue
 key=(x['assay_id'],x['guide']['sequence'],x['passenger']['sequence'],digest([x['guide']['terminal'],x['passenger']['terminal']]))
 blocks[key].add(digest([x['guide'],x['passenger']]))
checks={'matched_block_count':len(blocks),'chemistry_states_per_block_histogram':dict(collections.Counter(len(v) for v in blocks.values())),'maximum_chemistry_states_per_block':max(map(len,blocks.values())),'no_chemistry_effect_prediction_max_absolute':float(effects[effects.model=='gnn_no_chemistry'].predicted_difference.abs().max()),'interpretation':'No-chemistry B2 differences at float32 roundoff are not chemical sensitivity; its raw near-constant correlation/calibration diagnostics are numerically unstable and uninterpretable. Retain original values without treating them as evidence.','figure_design':'Chosen after test evaluation; metrics, split, fits and predictions unchanged.'}
write_json(r/'post_evaluation_descriptive_audit.json',checks)
write_json(out/'manifest.json',{'plots':manifest,'source_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in data.iterdir() if p.is_file()},'no_new_fits':True});files=[str(p.relative_to(r)) for p in out.rglob('*') if p.is_file()]+['post_evaluation_descriptive_audit.json'];write_json(r/'figures.outputs.json',files);print(json.dumps({'figures':len(manifest),'audit':checks}))
