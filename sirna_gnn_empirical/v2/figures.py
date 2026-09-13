from core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import rankdata
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
out=RUN/'figures';out.mkdir(exist_ok=True)
def save(fig,name,data=None):
 fig.tight_layout()
 for ext in ['pdf','svg','png']:fig.savefig(out/(name+'.'+ext),dpi=180,bbox_inches='tight')
 if data is not None:data.to_csv(out/(name+'.csv'),index=False)
 plt.close(fig)
name={'gnn':'GNN','chemistry_tree':'Chemistry tree','token_cnn':'Token CNN','nongraph':'No-message','gnn_no_chemistry':'No chemistry','gnn_activity':'Activity only','gnn_pair':'Pair GNN','nongraph_pair':'Pair no-message','training_mean':'Fold mean','training_majority':'Fold majority','pair_ridge':'Pair ridge','zero':'Zero'}
b1=pd.read_csv(RUN/'B1_comparison.csv');b1=b1[b1.seed=='ensemble'];ci=pd.read_csv(RUN/'paired_comparisons.csv');g=pd.read_csv(RUN/'grouped_comparison.csv');g=g[g.seed=='ensemble'];b=pd.read_csv(RUN/'B2_supervised_comparison.csv');b=b[b.seed=='ensemble']
fig,ax=plt.subplots(1,2,figsize=(7,2.65));q=ci[(ci.phase=='factorial')&ci.model.isin(['gnn_R1C1','nongraph_R1C1','token_cnn_R1C1','gnn_no_chemistry_R1C1'])].copy();xx=np.arange(len(q));ax[0].errorbar(q.difference_mse,xx,xerr=[q.difference_mse-q.lower95,q.upper95-q.difference_mse],fmt='o',color='#1c6d91');ax[0].axvline(0,color='gray',ls='--');ax[0].set_yticks(xx,[name[x.replace('_R1C1','')] for x in q.model]);ax[0].set_xlabel('APP MSE minus corrected tree');ax[0].set_title('Original-split corrected refits')
xx=np.arange(len(g));ax[1].bar(xx-.18,g.pooled_mse,.36,label='Pooled rows',color='#1c6d91');ax[1].bar(xx+.18,g.equal_component_mse,.36,label='Equal component',color='#da8d43');ax[1].set_xticks(xx,[name[x] for x in g.model],rotation=25,ha='right');ax[1].set_ylim(0,.09);ax[1].set_ylabel('ENsiRNA outer-fold MSE');ax[1].legend(fontsize=8);save(fig,'main_b1',q)
r=pd.read_csv(RUN/'ranking_by_assay.csv');r=r[(r.seed=='ensemble')&(r.phase=='deployment')];pools=sorted(r.assay_id.unique());fig,ax=plt.subplots(figsize=(7,2.5))
for kind in ['gnn','chemistry_tree','token_cnn','nongraph']:
 a=r[r.model==kind].set_index('assay_id').loc[pools];ax.plot(range(1,18),a.top5_percentile,marker='o',ms=3,label=name[kind])
ax.axhline(.5,c='gray',ls='--');ax.set_xticks(range(1,18));ax.set_xlabel('Primary APP assay pool (IDs in appendix)');ax.set_ylabel('Top-five measured percentile');ax.set_ylim(0,1);ax.legend(ncol=4,fontsize=8,loc='lower center');save(fig,'main_ranking',r)
bp=pd.read_csv(RUN/'B2_supervised_predictions.csv');bp=bp[bp.seed=='ensemble'];gg=bp[bp.model=='gnn_pair'];fig,ax=plt.subplots(1,2,figsize=(7,2.65));
for f,a in gg.groupby('outer_fold'):ax[0].scatter(a.observed_difference,a.predicted_difference,s=14,alpha=.7,label=f'Fold {f+1}')
ax[0].plot([-1.7,.6],[-1.7,.6],c='gray',ls='--');ax[0].axhline(0,c='gray',lw=.5);ax[0].set_xlabel('Measured bundle difference');ax[0].set_ylabel('Predicted difference');ax[0].legend(fontsize=7,ncol=2)
vals=[(gg.observed_difference<-.02).sum(),(abs(gg.observed_difference)<=.02).sum(),(gg.observed_difference>.02).sum()];preds=[(gg.predicted_difference<-.02).sum(),(abs(gg.predicted_difference)<=.02).sum(),(gg.predicted_difference>.02).sum()];xx=np.arange(3);ax[1].bar(xx-.18,vals,.36,label='Measured',color='#1c6d91');ax[1].bar(xx+.18,preds,.36,label='Pair GNN',color='#da8d43');ax[1].set_xticks(xx,['Negative','Within band','Positive']);ax[1].set_ylabel('Pairs');ax[1].legend(fontsize=8);save(fig,'main_b2',gg)
# Every source group, not just fold-level averages.
gd=pd.read_csv(RUN/'grouped_component_results.csv');gd=gd[gd.seed=='ensemble'];pv=gd.pivot(index='study_group',columns='model',values='mse');fig,ax=plt.subplots(figsize=(7,3));pv.rename(columns=name).plot.bar(ax=ax);ax.set_xticklabels([str(i+1) for i in range(len(pv))],rotation=0);ax.set_xlabel('Study-linked component (mapping in table)');ax.set_ylabel('Outer-test MSE');ax.legend(fontsize=8,ncol=4);save(fig,'appendix_group_results',gd)
# Coherent full-fold B2 scatter.
fig,axs=plt.subplots(2,2,figsize=(7,5))
for f,ax in enumerate(axs.flat):
 for kind in ['gnn_pair','nongraph_pair','training_mean']:
  a=bp[(bp.model==kind)&(bp.outer_fold==f)];ax.scatter(a.observed_difference,a.predicted_difference,s=10,alpha=.7,label=name[kind])
 ax.plot([-1.7,.6],[-1.7,.6],c='gray',ls='--');ax.set_title(f'Outer fold {f+1}');ax.set_xlabel('Measured difference');ax.set_ylabel('Predicted difference')
axs.flat[0].legend(fontsize=7);save(fig,'appendix_b2_folds',bp)
fig,axs=plt.subplots(2,3,figsize=(9,5))
for kind,ax in zip(['zero','training_mean','pair_ridge','gnn_activity','gnn_pair','nongraph_pair'],axs.flat):
 a=bp[bp.model==kind];ax.scatter(a.observed_difference,a.predicted_difference,s=9,c=a.outer_fold,cmap='tab10');ax.plot([-1.7,.6],[-1.7,.6],c='gray',ls='--');ax.set_title(name[kind]);ax.set_xlabel('Measured difference');ax.set_ylabel('Predicted difference')
save(fig,'appendix_b2_controls',bp)
rel=pd.read_csv(RUN/'relation_support.csv');pv=rel.pivot(index='name',columns='split',values='edges');fig,ax=plt.subplots(figsize=(7,3));im=ax.imshow(np.log10(pv+1),aspect='auto',cmap='Blues');ax.set_yticks(range(8),pv.index);ax.set_xticks(range(4),pv.columns);fig.colorbar(im,ax=ax,label='log10(edge count + 1)');save(fig,'appendix_relation_support',rel)
ctx=pd.read_csv(RUN/'context_support.csv');fig,ax=plt.subplots(figsize=(7,3));p=ctx.assign(magnitude=np.maximum(abs(ctx['min']),abs(ctx['max']))).pivot(index='name',columns='split',values='magnitude');p.plot.barh(ax=ax,logx=True);ax.set_xlabel('Maximum absolute standardized input (log scale)');save(fig,'appendix_context_support',ctx)
ff=pd.read_csv(RUN/'all_fit_summary.csv');f=ff[(ff.kind=='gnn')&ff.tag.str.contains('/final')];fig,ax=plt.subplots(figsize=(7,3));ax.scatter(range(len(f)),f.epochs_executed,label='Executed',marker='x');ax.scatter(range(len(f)),f.selected_epoch,label='Selected',s=15);ax.set_xlabel('Final GNN fit, all tasks/folds/seeds');ax.set_ylabel('Epoch');ax.legend();save(fig,'appendix_selection',f)
fig,ax=plt.subplots(figsize=(7,3));p=b1[(b1.phase=='factorial')&(b1.model=='gnn')];ax.bar(p.arm,p.APP_mse,color=['#bbb','#87adc0','#dcac7c','#1c6d91']);ax.set_ylabel('Retrospective APP MSE');save(fig,'appendix_factorial',p)
frozen=b1[(b1.phase=='frozen')&(b1.model=='gnn')];fig,ax=plt.subplots(figsize=(7,3));ax.bar(frozen.arm,frozen.APP_mse,color='#1c6d91');ax.set_ylabel('APP MSE, frozen modified inference');save(fig,'appendix_frozen',frozen)
# Deterministic ID ties do not establish ranking by a constant predictor. Report exact random tie expectation separately.
a=pd.read_csv(RUN/'all_B1_predictions.csv');a=a[a.seed=='ensemble'];rows,_=load();md=pd.DataFrame([{k:x.get(k) for k in ['record_id','activity','assay_id']} for x in rows]);a=a[a['split']=='test_APP'].merge(md,on='record_id');tie=[]
for keys,part in a.groupby(['phase','model','arm']):
 for aid,g in part.groupby('assay_id'):
  ranks=(rankdata(g.activity,method='average')-1)/(len(g)-1);gg=g.assign(percentile=ranks).sort_values('prediction',ascending=False);cut=gg.prediction.iloc[4];above=gg[gg.prediction>cut];tied=gg[gg.prediction==cut];left=5-len(above);tie.append(dict(phase=keys[0],model=keys[1],arm=keys[2],assay_id=aid,tie_count_at_boundary=len(tied),expected_top5_percentile=(above.percentile.sum()+left*tied.percentile.mean())/5,expected_selected_activity=(above.activity.sum()+left*tied.activity.mean())/5))
pd.DataFrame(tie).to_csv(RUN/'ranking_tie_expectations.csv',index=False)
if (RUN/'external_predictions_with_labels.csv').exists():
 e=pd.read_csv(RUN/'external_predictions_with_labels.csv');e=e[e.seed=='ensemble'];fig,axs=plt.subplots(1,4,figsize=(10,2.8))
 for kind,ax in zip(['gnn','chemistry_tree','token_cnn','nongraph'],axs):
  z=e[e.model==kind];ax.scatter(z.activity,z.prediction,s=9,alpha=.7);ax.plot([0,1],[0,1],c='gray',ls='--');ax.set_title(name[kind]);ax.set_xlabel('Measured activity');ax.set_ylabel('Predicted activity')
 save(fig,'appendix_external',e)
write(RUN/'figure_manifest.json',[dict(file=str(p.relative_to(RUN)),sha256=sha(p)) for p in sorted(out.iterdir())]);print('figures',len(list(out.glob('*.pdf'))))
