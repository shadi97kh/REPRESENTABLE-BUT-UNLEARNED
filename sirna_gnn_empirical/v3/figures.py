"""Result figures with explicit machine-readable plotted-mark provenance."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import shutil
PAPER=ROOT/'papers/interaction_recoverability_iclr2027/v6';OUT=RUN/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.bbox':'tight','axes.titleweight':'bold','axes.labelsize':11})
colors=['#9bbad9','#a9d6b8','#d3b5d9','#e9b6a8','#e3ca8a','#9dd2d5','#bbc3cd','#c9ca9c']
methods=['original_gnn','corrected_gnn','original_no_message','corrected_no_message','chemistry_tree','token_cnn','guide_ridge','pairwise_ridge'];short=['R0 G','R1 G','R0 N','R1 N','Tree','CNN','Guide','Pairwise'];names=dict(zip(methods,short));names.update(activity_gnn='Activity GNN',pair_gnn='Pair GNN',pair_no_message='Pair no-msg',pair_tree='Pair tree',pair_ridge='Pair ridge',training_mean='Training mean',zero='Zero',training_row_mean='Training row mean',training_equal_group_mean='Training group mean',oracle_test_mean_diagnostic='Oracle test mean')
provenance=[]
def read(name):return pd.read_csv(RUN/name)
def save(fig,name,marks,sources,description):
 for ext in ['pdf','svg','png']:fig.savefig(OUT/f'{name}.{ext}',dpi=240)
 plt.close(fig);marks.to_csv(OUT/f'{name}_marks.csv',index=False)
 for ext in ['pdf','svg','png']:shutil.copy2(OUT/f'{name}.{ext}',PAPER/'figures'/f'{name}.{ext}')
 provenance.append(dict(figure=name,sources=[dict(path=s,sha256=sha(RUN/s)) for s in sources],marks=f'figures/{name}_marks.csv',transformation=description,code_sha256=sha(HERE/'figures.py')))
activity=read('tables/activity_metrics.csv');sm=read('tables/activity_seed_metrics.csv');ensemble=read('evaluated/activity_ensemble_predictions.csv');grad=read('tables/relation_gradients_by_fit.csv');rel=read('tables/final_relation_counts.csv')
# Figure 2: counts and exact data gradients; logarithm represents zero separately.
z=rel[(rel.fit.str.contains('activity/deployment/final/original_gnn/'))&(rel.seed==1103)&(rel.role.isin(['train','test']))];gg=grad[(grad.fit.str.contains('activity/deployment/final/'))&grad.method.isin(['original_gnn','corrected_gnn'])];ga=gg.groupby(['method','parameter','relation'],as_index=False).data_gradient_norm.mean();fig,ax=plt.subplots(1,2,figsize=(7.3,3.25));rnames=['PO+','PS+','?+','PO−','PS−','?−','Pair','Mismatch'];marks=[]
for j,role in enumerate(['train','test']):
 q=z[z.role==role].sort_values('relation');ax[0].bar(np.arange(8)+(j-.5)*.36,q.directed_edges,width=.36,color=colors[j],label='ENsiRNA train' if role=='train' else 'APP + S7');marks.extend([dict(panel='counts',role=role,**row) for row in q[['relation','directed_edges']].to_dict('records')])
 for x,y in zip(np.arange(8)+(j-.5)*.36,q.directed_edges):
  if y<=2:ax[0].text(x,max(1.3,y*1.6),str(int(y)),ha='center',fontsize=10)
ax[0].set_yscale('symlog',linthresh=1);ax[0].set_xticks(range(8),rnames,rotation=30);ax[0].set_ylabel('Directed edge count');ax[0].set_title('A. Exposure by relation');
for j,m in enumerate(['original_gnn','corrected_gnn']):
 for k,layer in enumerate(sorted(ga.parameter.unique())):
  q=ga[(ga.method==m)&(ga.parameter==layer)].sort_values('relation');ax[1].plot(q.relation,q.data_gradient_norm,marker='o' if k==0 else 's',linestyle='-' if k==0 else '--',color=colors[j],label=names[m]+f' layer {k+1}');marks.extend([dict(panel='gradient',method=m,layer=k+1,**row) for row in q[['relation','data_gradient_norm']].to_dict('records')])
ax[1].set_yscale('symlog',linthresh=1e-6);ax[1].set_xticks(range(8),rnames,rotation=30);ax[1].set_ylabel('Data-gradient norm');ax[1].set_title('B. Data gradients');ha,la=ax[0].get_legend_handles_labels();hb,lb=ax[1].get_legend_handles_labels();fig.legend(ha+hb,la+lb,loc='lower center',bbox_to_anchor=(.5,-.04),ncol=3,fontsize=9);fig.tight_layout(rect=(0,.12,1,1));save(fig,'main_support',pd.DataFrame(marks),['tables/final_relation_counts.csv','tables/relation_gradients_by_fit.csv'],'Counts for deployment seed1103, train and combined APP/S7 test. Gradient mean across all ten deployment seeds, by exact layer/relation, activity loss only.')
# Figure 3: seed risk and descriptive calibration on shifted cohorts.
fig,ax=plt.subplots(2,2,figsize=(7.3,4.8));marks=[]
for col,cohort in enumerate(['APP','Davis_S7']):
 for j,m in enumerate(methods[:6]):
  q=sm[(sm.cohort==cohort)&(sm.method==m)&(sm.weighting=='rows')];x=j+np.linspace(-.13,.13,len(q));ax[0,col].scatter(x,q.mse,c=colors[j],s=15,edgecolor='#555555',linewidth=.3);risk=activity[(activity.cohort==cohort)&(activity.method==m)&(activity.weighting=='rows')].iloc[0].mse;ax[0,col].plot(j,risk,'k_',ms=11);marks.extend([dict(panel='seed risk',**row) for row in q.to_dict('records')]);marks.append(dict(panel='ensemble risk',cohort=cohort,method=m,mse=risk))
 for m,style in [('training_equal_group_mean','--'),('oracle_test_mean_diagnostic',':')]:
  value=activity[(activity.cohort==cohort)&(activity.method==m)&(activity.weighting=='rows')].iloc[0].mse;ax[0,col].axhline(value,color='#555555',ls=style,lw=1,label=names[m]);marks.append(dict(panel='constant',cohort=cohort,method=m,mse=value))
 ax[0,col].set_xticks(range(6),short[:6],rotation=30,ha='right');ax[0,col].set_ylabel('MSE');ax[0,col].set_title(cohort.replace('_',' '));
 for j,m in enumerate(['corrected_gnn','corrected_no_message','chemistry_tree']):
  q=ensemble[(ensemble.cohort==cohort)&(ensemble.method==m)].copy();q['bin']=pd.qcut(q.prediction,5,duplicates='drop');v=q.groupby('bin',observed=True).agg(prediction=('prediction','mean'),activity=('activity','mean'),n=('record_id','size')).reset_index(drop=True);ax[1,col].plot(v.prediction,v.activity,marker='o',color=colors[[1,3,4][j]],label=names[m]);marks.extend([dict(panel='calibration',cohort=cohort,method=m,**row) for row in v.to_dict('records')])
 limits=ax[1,col].get_xlim();yr=ax[1,col].get_ylim();lo=min(limits[0],yr[0]);hi=max(limits[1],yr[1]);ax[1,col].plot([lo,hi],[lo,hi],color='#777777',ls=':',lw=1);ax[1,col].set_xlabel('Mean prediction (up to five bins)');ax[1,col].set_ylabel('Mean observed activity');ax[1,col].legend(fontsize=10)
fig.legend(*ax[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,1.02),ncol=2,fontsize=9);fig.tight_layout(rect=(0,0,1,.95));save(fig,'main_activity',pd.DataFrame(marks),['tables/activity_metrics.csv','tables/activity_seed_metrics.csv','evaluated/activity_ensemble_predictions.csv'],'Row MSE by seed; ensemble MSE black marks; distinct trained/oracle constants; five-quantile partition with duplicate cuts dropped for descriptive calibration, not fitted recalibration.')
# Figure 4: ranking and pair predictions, no false independent-pool interval.
rk=read('tables/ranking_pools.csv');bp=read('evaluated/pair_ensemble_predictions.csv');fig,ax=plt.subplots(1,2,figsize=(7.3,3.65));marks=[]
for j,m in enumerate(['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn']):
 q=rk[rk.method==m];ax[0].scatter(np.full(len(q),j)+np.linspace(-.16,.16,len(q)),q.top5_mean_percentile,color=colors[[1,3,4,5][j]],s=18,edgecolor='#555555',linewidth=.3);ax[0].plot(j,q.top5_mean_percentile.mean(),'k_',ms=15);marks.extend([dict(panel='ranking',**row) for row in q.to_dict('records')]);marks.append(dict(panel='ranking_mean',method=m,pools=len(q),top5_mean_percentile=float(q.top5_mean_percentile.mean())))
ax[0].axhline(.5,color='#777777',ls=':',label='Random tie expectation');ax[0].set_xticks(range(4),['GNN','No-msg','Tree','CNN']);ax[0].set_ylabel('Top-five mean percentile');ax[0].set_ylim(-.025,1.025);ax[0].legend(fontsize=9,loc='lower left');ax[0].set_title('A. APP candidate ranking')
for j,(fold,q) in enumerate(bp[bp.method=='pair_gnn'].groupby('population')):
 ax[1].scatter(q.observed_difference,q.prediction,color=colors[j],s=15,alpha=.7,label='Fold '+fold[-1]);marks.extend([dict(panel='B2',population=fold,pair_id=t.pair_id,method='pair_gnn',observed_difference=t.observed_difference,prediction=t.prediction) for t in q.itertuples()])
q=bp[bp.method=='training_mean'];ax[1].scatter(q.observed_difference,q.prediction,color='#777777',s=15,alpha=.6,label='Training-fold mean',marker='x');marks.extend([dict(panel='B2',population=t.population,pair_id=t.pair_id,method='training_mean',observed_difference=t.observed_difference,prediction=t.prediction) for t in q.itertuples()])
lo=min(bp.observed_difference.min(),bp[bp.method=='pair_gnn'].prediction.min());hi=max(bp.observed_difference.max(),bp[bp.method=='pair_gnn'].prediction.max());ax[1].plot([lo,hi],[lo,hi],color='#777777',ls=':',lw=1);ax[1].axhline(0,color='#aaaaaa',lw=.7);ax[1].axvspan(-.02,.02,color='#eeeeee');ax[1].set_xlabel('Measured changed − reference');ax[1].set_ylabel('Predicted difference');ax[1].set_title('B. Measured B2 differences');plotted=bp[bp.method.isin(['pair_gnn','training_mean'])].prediction;ypmin=min(float(plotted.min()),0.);ypmax=max(float(plotted.max()),0.);ypad=.08*max(ypmax-ypmin,.1);ax[1].set_ylim(ypmin-ypad,ypmax+ypad);ax[1].legend(fontsize=8,ncol=2,loc='upper right',framealpha=.85);fig.tight_layout();save(fig,'main_pair_ranking',pd.DataFrame(marks),['tables/ranking_pools.csv','evaluated/pair_ensemble_predictions.csv'],'Each pool score is exact tie-aware top five; black mark is unweighted pool mean; B2 fixed ensembles and training-fold means on all eligible pairs.')
# Complete calibration/discrimination panels for all cohorts and principal comparators.
fig,ax=plt.subplots(3,3,figsize=(7.3,7.2));marks=[]
for i,c in enumerate(['ENsiRNA_grouped','APP','Davis_S7']):
 for j,m in enumerate(['corrected_gnn','corrected_no_message','chemistry_tree']):
  q=ensemble[(ensemble.cohort==c)&(ensemble.method==m)];ax[i,j].scatter(q.activity,q.prediction,s=4,alpha=.25,color=colors[methods.index(m)]);lo=min(q.activity.min(),q.prediction.min());hi=max(q.activity.max(),q.prediction.max());ax[i,j].plot([lo,hi],[lo,hi],color='#555555',ls=':',lw=.8);ax[i,j].set_title(c.replace('ENsiRNA_grouped','Grouped')+' / '+names[m],fontsize=10);ax[i,j].set_xlabel('Observed');ax[i,j].set_ylabel('Predicted');marks.extend(q[['cohort','method','record_id','activity','prediction']].to_dict('records'))
fig.tight_layout();save(fig,'appendix_calibration',pd.DataFrame(marks),['evaluated/activity_ensemble_predictions.csv'],'Every unaltered observation prediction, no clipping or fitted intercept correction.')
# All seed distributions and ensemble-member decomposition.
fig,ax=plt.subplots(3,1,figsize=(7.3,6.5));marks=[]
for i,c in enumerate(['ENsiRNA_grouped','APP','Davis_S7']):
 for j,m in enumerate(methods):
  q=sm[(sm.cohort==c)&(sm.method==m)&(sm.weighting=='rows')];ax[i].scatter(np.full(len(q),j)+np.linspace(-.15,.15,len(q)),q.mse,c=colors[j],s=20,edgecolor='#555555',linewidth=.3);e=activity[(activity.cohort==c)&(activity.method==m)&(activity.weighting=='rows')].iloc[0];ax[i].plot(j,e.mse,'k_',ms=15);marks.extend([dict(mark='member',**v) for v in q.to_dict('records')]);marks.append(dict(mark='ensemble',cohort=c,method=m,mse=e.mse))
 ax[i].set_xticks(range(8),short);ax[i].set_title(c.replace('_',' '));ax[i].set_ylabel('MSE')
fig.tight_layout();save(fig,'appendix_seed_risk',pd.DataFrame(marks),['tables/activity_seed_metrics.csv','tables/activity_metrics.csv'],'All ten stochastic seed MSEs, one deterministic ridge risk, and black ensemble-risk markers.')
# Detailed selected versus executed optimization trajectories, final deployment and pair protocols.
curves=read('evaluated/all_learning_curves.csv');fits=read('tables/all_fit_summary.csv');fig,ax=plt.subplots(2,3,figsize=(7.3,4.9));marks=[]
for a,m in zip(ax.flat,[*methods[:4],'token_cnn','pair_gnn']):
 ff=fits[fits.name.str.contains('pair/b2outer0/final/pair_gnn/') if m=='pair_gnn' else fits.name.str.contains('activity/deployment/final/')&(fits.method==m)]
 for row in ff.itertuples():
  z=curves[curves.fit==row.name];a.plot(z.epoch,z.validation_score,alpha=.7,lw=.8,color=plt.get_cmap('tab10')(SEEDS.index(row.seed)),label=str(row.seed));a.plot(row.selected_epoch,row.validation_score,'o',ms=2,color='#555555');marks.extend([dict(mark='trajectory',**v) for v in z.to_dict('records')]);marks.append(dict(mark='selected_checkpoint',fit=row.name,method=m,seed=row.seed,epoch=row.selected_epoch,validation_score=row.validation_score))
 a.set_title(names[m]+(' / fold 0' if m=='pair_gnn' else ''));a.set_xlabel('Executed epoch');a.set_ylabel('Validation MSE')
fig.legend(*ax.flat[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,-.03),ncol=5,fontsize=9);fig.tight_layout(rect=(0,.06,1,1));save(fig,'appendix_learning_curves',pd.DataFrame(marks),['evaluated/all_learning_curves.csv','tables/all_fit_summary.csv'],'Every final deployment neural seed trajectory and ten outer0 pair-GNN trajectories; selected checkpoint dots. Tree training has no neural epochs.')
# Full B2 model scatter, preserving unfavorable predictions.
fig,ax=plt.subplots(2,3,figsize=(7.3,4.9));marks=[]
for a,m in zip(ax.flat,['activity_gnn','pair_gnn','pair_no_message','pair_tree','pair_ridge','training_mean']):
 q=bp[bp.method==m];a.scatter(q.observed_difference,q.prediction,s=9,alpha=.65,color=colors[{'activity_gnn':0,'pair_gnn':1,'pair_no_message':3,'pair_tree':4,'pair_ridge':7,'training_mean':6}[m]]);lo=min(q.observed_difference.min(),q.prediction.min());hi=max(q.observed_difference.max(),q.prediction.max());a.plot([lo,hi],[lo,hi],ls=':',color='#555555',lw=.7);a.axhline(0,lw=.5,color='#777777');a.set_title(names[m]);a.set_xlabel('Measured difference');a.set_ylabel('Prediction');marks.extend(q[['pair_id','method','observed_difference','prediction']].to_dict('records'))
fig.tight_layout();save(fig,'appendix_B2_models',pd.DataFrame(marks),['evaluated/pair_ensemble_predictions.csv'],'Complete out-of-fold B2 predictions including training-fold mean; shared endpoints/assays remain dependent.')
# Every APP pool and direct matched differences.
w=read('tables/ranking_direct_pool_differences.csv');fig,ax=plt.subplots(2,1,figsize=(7.3,5));x=np.arange(len(w));marks=[]
for j,m in enumerate(['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn']):ax[0].plot(x,w[m],marker='o',ms=3,lw=.8,color=colors[methods.index(m)],label=names[m])
ax[0].axhline(.5,color='#888888',ls=':');ax[0].legend(ncol=4,fontsize=10);ax[0].set_ylabel('Top-five percentile');ax[0].set_xticks(x,[f'P{i+1}' for i in x]);ax[1].bar(x-.18,w.tree_minus_corrected_gnn,width=.36,color=colors[4],label='Tree − R1 GNN');ax[1].bar(x+.18,w.corrected_gnn_minus_no_message,width=.36,color=colors[1],label='R1 GNN − R1 no-msg');ax[1].axhline(0,color='#888888',lw=.6);ax[1].set_xticks(x,[f'P{i+1}' for i in x]);ax[1].set_ylabel('Matched pool difference');ax[1].legend(fontsize=10);fig.tight_layout();save(fig,'appendix_ranking_pools',w,['tables/ranking_direct_pool_differences.csv'],'Sorted pool table order, all exact tie-aware scores and descriptive differences; no independent-pool test.')
# Measured B3 heatmap and fixed source-held-out predictions.
b3=read('b3_primary/interaction_predictions.csv');q=b3[(b3.seed.astype(str)=='ensemble')&(b3.method=='corrected_gnn')];fig=plt.figure(figsize=(7.3,4.1));grid=fig.add_gridspec(1,3,width_ratios=[1,1,.05]);ax=[fig.add_subplot(grid[0,0]),fig.add_subplot(grid[0,1])];cax=fig.add_subplot(grid[0,2]);limit=max(abs(q.observed_interaction).max(),abs(q.prediction).max());marks=[]
for a,c,title in [(ax[0],'observed_interaction','Measured'),(ax[1],'prediction','R1 GNN')]:
 mat=q.pivot(index='AS',columns='SS',values=c);im=a.imshow(mat,cmap='RdBu_r',vmin=-limit,vmax=limit,aspect='auto');a.set_xticks(range(len(mat.columns)),mat.columns,rotation=70,fontsize=9);a.set_yticks(range(len(mat.index)),mat.index,fontsize=9);a.set_title(title);a.set_xlabel('SS chemistry state');a.set_ylabel('AS chemistry state' if a is ax[0] else '');
 if a is ax[1]:a.set_yticklabels([])
fig.colorbar(im,cax=cax,label='Activity-scale interaction');fig.subplots_adjust(left=.15,right=.89,bottom=.25,top=.91,wspace=.15);save(fig,'appendix_B3_panel',q,['b3_primary/interaction_predictions.csv'],'All 140 measured rectangles and corresponding mean of fixed outer0 seed predictions; common symmetric color range; one background, no independence assumption.')
fig,ax=plt.subplots(2,3,figsize=(7.3,5));marks=[]
for a,m in zip(ax.flat,methods[:6]):
 q=b3[(b3.seed.astype(str)=='ensemble')&(b3.method==m)];a.scatter(q.observed_interaction,q.prediction,s=8,color=colors[methods.index(m)],alpha=.7);a.axhline(0,c='#888888',lw=.6);a.set_title(names[m]);a.set_xlabel('Measured interaction');a.set_ylabel('Predicted interaction');marks.extend(q.to_dict('records'))
fig.tight_layout();save(fig,'appendix_B3_predictions',pd.DataFrame(marks),['b3_primary/interaction_predictions.csv'],'All source-held-out B3 ensemble predictions, without independent-rectangle confidence intervals.')
# Primary/release disagreement is shown, not silently cleaned away.
d=read('b3/primary_reconciliation.csv');d=d[d.primary_activity.notna()];fig,ax=plt.subplots(figsize=(5.2,3.5));ax.scatter(d.primary_activity,d.archived_activity,s=8,alpha=.5,color=colors[2]);lo=min(d.primary_activity.min(),d.archived_activity.min());hi=max(d.primary_activity.max(),d.archived_activity.max());ax.plot([lo,hi],[lo,hi],ls=':',color='#555555');ax.set_xlabel('Primary activity (1 − relative eGFP)');ax.set_ylabel('Released ENsiRNA activity');ax.set_title('Identity-matched primary/release audit');fig.tight_layout();save(fig,'appendix_source_reconciliation',d,['b3/primary_reconciliation.csv'],'All uniquely matched primary/released pairs, including disagreements; identity mapping does not use model predictions.')
write(RUN/'figure_provenance.json',provenance);print('Generated figures:',len(provenance))
