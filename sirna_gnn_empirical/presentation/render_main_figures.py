"""Presentation-only: redraw the committed plot inputs; no fitting or metric estimation."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent;D=P/'figure_data';O=P/'figures'
for name,h in json.loads((D/'source_hashes.json').read_text()).items():
 assert hashlib.sha256((D/name).read_bytes()).hexdigest()==h
plt.rcParams.update({'font.size':9.5,'axes.titlesize':10,'axes.labelsize':9.5,'xtick.labelsize':8.5,'ytick.labelsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
name={'gnn':'GNN','chemistry_tree':'Tree','token_cnn':'Token CNN','nongraph':'No-message','gnn_no_chemistry':'No chemistry'}
blue='#1c6d91';orange='#da8d43';colors={'gnn':blue,'chemistry_tree':orange,'token_cnn':'#4f8a5b','nongraph':'#9468a0'}
def save(fig,n):
 fig.savefig(O/(n+'.pdf'),bbox_inches='tight',pad_inches=.04)
 fig.savefig(O/(n+'.svg'),bbox_inches='tight',pad_inches=.04)
 plt.close(fig)
q=pd.read_csv(D/'main_b1.csv');g=pd.read_csv(D/'grouped_comparison.csv');g=g[g.seed=='ensemble'];order=['gnn','chemistry_tree','token_cnn','nongraph'];g=g.set_index('model').loc[order]
fig,ax=plt.subplots(1,2,figsize=(5.5,3.05),layout='constrained',gridspec_kw={'width_ratios':[1.12,1]})
y=np.arange(len(q));ax[0].errorbar(q.difference_mse,y,xerr=[q.difference_mse-q.lower95,q.upper95-q.difference_mse],fmt='o',capsize=3,color=blue)
ax[0].axvline(0,color='0.5',ls='--',lw=.8);ax[0].set_yticks(y,[name[x.replace('_R1C1','')] for x in q.model]);ax[0].set_xlabel('MSE minus corrected tree');ax[0].set_title('(a) APP refits');ax[0].set_ylim(-.6,len(q)-.4);ax[0].grid(axis='x',alpha=.17)
x=np.arange(len(g));ax[1].bar(x-.18,g.pooled_mse,.36,label='Pooled rows',color=blue);ax[1].bar(x+.18,g.equal_component_mse,.36,label='Equal component',color=orange);ax[1].set_xticks(x,[name[k] for k in order],rotation=38,ha='right');ax[1].set_ylim(0,.09);ax[1].set_ylabel('Outer-test MSE');ax[1].set_title('(b) Grouped ENsiRNA');ax[1].legend(fontsize=8,frameon=False,loc='upper right');ax[1].grid(axis='y',alpha=.17);ax[1].set_axisbelow(True);save(fig,'main_b1')
r=pd.read_csv(D/'main_ranking.csv');pools=sorted(r.assay_id.unique());assert len(pools)==17
fig,ax=plt.subplots(figsize=(5.5,2.9),layout='constrained')
for k in order:
 a=r[r.model==k].set_index('assay_id').loc[pools];ax.plot(range(1,18),a.top5_percentile,marker='o',ms=4,lw=1.4,label=name[k],color=colors[k])
ax.axhline(.5,c='0.5',ls='--',lw=.8);ax.set_xticks(range(1,18));ax.set_xlabel('APP assay pool (IDs and sizes in appendix)');ax.set_ylabel('Top-five measured percentile');ax.set_ylim(0,1);ax.legend(ncol=4,fontsize=8.5,loc='upper center',bbox_to_anchor=(.5,1.16),frameon=False,columnspacing=1);ax.grid(axis='y',alpha=.17);save(fig,'main_ranking')
gg=pd.read_csv(D/'main_b2.csv');scores=pd.read_csv(D/'B2_supervised_comparison.csv');row=scores[(scores.model=='gnn_pair')&(scores.seed=='ensemble')].iloc[0]
fig,ax=plt.subplots(1,2,figsize=(5.5,3.05),layout='constrained')
for fold,a in gg.groupby('outer_fold'):
 ax[0].scatter(a.observed_difference,a.predicted_difference,s=18,alpha=.7,label=f'Fold {fold+1}')
ax[0].plot([-1.7,.6],[-1.7,.6],c='0.5',ls='--',lw=.8);ax[0].axhline(0,c='0.5',lw=.5);ax[0].set_xlabel('Measured difference');ax[0].set_ylabel('Predicted difference');ax[0].set_title('(a) Matched pair effects');ax[0].legend(fontsize=8,ncol=2,loc='lower right',frameon=False);ax[0].grid(alpha=.15)
observed=[row[f'observed_{k}'] for k in ['negative','tie','positive']];predicted=[row[f'predicted_{k}'] for k in ['negative','tie','positive']]
x=np.arange(3);a=ax[1].bar(x-.18,observed,.36,label='Measured',color=blue);b=ax[1].bar(x+.18,predicted,.36,label='Pair GNN',color=orange);ax[1].bar_label(a,padding=2,fontsize=8);ax[1].bar_label(b,padding=2,fontsize=8);ax[1].set_xticks(x,['Negative','Within\nband','Positive']);ax[1].set_ylabel('Pairs');ax[1].set_ylim(0,190);ax[1].set_title('(b) Direction counts');ax[1].legend(fontsize=8,frameon=False,loc='upper right');ax[1].grid(axis='y',alpha=.15);ax[1].set_axisbelow(True);save(fig,'main_b2')
print('Redrew three figures from unchanged committed data; no metrics or models recomputed.')
