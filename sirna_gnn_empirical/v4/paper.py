"""Generate compact tables and final manuscript from verified saved evaluations; no fits."""
from common import *
import re, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';P=ROOT/'papers/interaction_recoverability_iclr2027/v9';S=P/'source';O=RUN/'paper_tables';O.mkdir(exist_ok=True)
assert (RUN/'evaluation/summary.json').exists() and (RUN/'visibility/conventions.json').exists()
methods=['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn'];names=dict(corrected_gnn='R1 GNN',corrected_no_message='R1 no-message',chemistry_tree='Chemistry tree',token_cnn='Token CNN',training_row_mean='Training row mean',training_equal_group_mean='Training group mean',oracle_test_mean_diagnostic='Oracle test mean')
scenarios=['source_excluded','primary_reanchored'];sn={'source_excluded':'Source-excluded','primary_reanchored':'Primary-assay'};cohorts=['ENsiRNA_grouped','APP','Davis_S7'];cn={'ENsiRNA_grouped':'Grouped','APP':'APP','Davis_S7':'S7'}
read=lambda p:pd.read_csv(p)
activity=read(RUN/'evaluation/activity_metrics.csv');direct=read(RUN/'evaluation/direct_comparisons.csv');seeds=read(RUN/'evaluation/seed_metrics.csv');counts=read(RUN/'visibility/input_counts.csv');oldact=read(V3/'tables/activity_metrics.csv');oldpair=read(V3/'tables/pair_metrics.csv');source=read(RUN/'review_checks/source_scores.csv');conc=read(RUN/'review_checks/source_concentration.csv');cal=read(RUN/'review_checks/calibration_decomposition.csv');fit=read(RUN/'fit_summary.csv');completion=json.loads((RUN/'fit_completion.json').read_text());assert len(fit)==792
mapping=[]
def fmt(x,d=4):
 return '--' if pd.isna(x) else f'{x:.{d}f}'
def tex(s):
 return str(s).replace('&',r'\&').replace('_',r'\_').replace('%',r'\%').replace('#',r'\#')
def table(key,headers,rows,caption,label,align=None):
 pd.DataFrame(rows,columns=headers).to_csv(O/(key+'.csv'),index=False)
 col=align or 'l'+'r'*(len(headers)-1)
 result='\\begin{table}[ht]\n\\centering\\small\n\\begin{tabular}{@{}'+col+'@{}}\\toprule\n'+' & '.join(headers)+r'\\\midrule'+'\n'
 result+='\n'.join(' & '.join(map(str,row))+r'\\' for row in rows)
 result+='\n\\bottomrule\\end{tabular}\n\\caption{'+caption+'}\\label{'+label+'}\n\\end{table}\n'
 mapping.append(dict(object=label,table_csv='paper_tables/'+key+'.csv',scope=caption));return result
repl={}
rr=[]
for _,q in conc.iterrows():
 sid=str(q.source_family); ss=source[(source.subset=='source:'+sid)&(source.weighting=='rows')]
 # source subset keys are verified explicitly below
 if ss.empty:ss=source[(source.subset==sid)&(source.weighting=='rows')]
 assert len(ss),sid
 v=ss.set_index('method');label='Patent source (594)' if q.n==594 else sid
 rr.append([label,int(q.n),int(q.sequence_components),int(q.study_components),fmt(v.loc['corrected_gnn','r_squared'],3),fmt(v.loc['chemistry_tree','r_squared'],3)])
repl['SOURCE_TABLE']=table('source_counts_scores',['Released source','Rows','Seq.','Study','$R^2_G$','$R^2_T$'],rr,'All seventeen source labels and their saved out-of-fold row-weighted scores. Seq./Study are protected component counts within each source; they cannot be summed as independent studies. The patent row denotes US20120088815A1 / EP2415869A1. These are scoring subsets of released-trained models, not source-excluded refits. Bramsen (19282453) remains subject to the adjacent primary-source audit.','tab:data')
app_rows=readlines(OLD/'test_APP_observations.jsonl');s7_rows=readlines(PREV/'external_observations.jsonl');b3_rows=readlines(V3/'b3_primary/observations.jsonl');pairs=readlines(OLD/'eligible_pairs.jsonl')
repl['COHORT_TABLE']=table('evaluation_units',['Evidence','Unit','Count','Seq. comp.'],[['APP B1','Activity rows',len(app_rows),len({x['sequence_group'] for x in app_rows})],['Davis S7','Activity rows',len(s7_rows),len({x['sequence_group'] for x in s7_rows})],['APP B2','Measured pairs',len(pairs),12],['Primary B3','Measured states',len(b3_rows),1]],'Additional measured evidence, separate from the released grouped benchmark. APP has one related patent-family cohort and seventeen overlapping ranking pools. B2 has 26 exact duplex backgrounds; B3 yields 140 dependent four-condition contrasts. Row, pair and state counts are different units, not biological replication counts.','tab:cohorts','llrr')
repl['ADJUDICATION_TABLE']=table('adjudication',['Admitted Bramsen category','Rows'],[['Unique reported identity / primary measurement',1604],['No reported-identity match',201],['Multiple released rows for primary identity',76],['No unique primary measurement',43],['Total admitted / raw source rows','1,924 / 1,972']], 'Mutually exclusive admitted-source adjudication. The primary-assay sensitivity retains 1,604 linked records and quarantines the other 320; no model accuracy or outcome agreement selects a mapping. Original numeric fields are preserved.','tab:adjudication','lr')
ff=fit.assign(scenario=fit.name.str.split('/').str[0],phase=fit.name.str.split('/').str[2]);rr=[]
for s in scenarios:
 for m in methods:
  q=ff[(ff.scenario==s)&(ff.method==m)];rr.append([sn[s] if m==methods[0] else '',names[m],int((q.phase=='development').sum()),int((q.phase=='final').sum()),f'{q.executed_updates.sum():,}',f'{q.selected_updates.sum():,}'])
repl['FIT_TABLE']=table('new_fits',['Sensitivity','Method','Dev.','Final','Exec.','Selected'],rr,'Actual new fits and neural optimizer updates. Trees are fitted but have zero neural updates. Forty unaffected final checkpoints are reused separately. Every final stochastic block has ten seeds. The original 2,312 fits are preserved and are not counted again.','tab:optimization','llrrrr')
rr=[]
for s in scenarios:
 for m in methods+['training_row_mean','oracle_test_mean_diagnostic']:
  q=activity[(activity.scenario==s)&(activity.method==m)&(activity.weighting=='rows')].set_index('cohort');rr.append([sn[s] if m==methods[0] else '',names[m]]+[fmt(q.loc[c,k],4 if k=='mse' else 3) for c in cohorts for k in ['mse','r_squared']])
repl['FOCUSED_TABLE']=table('focused_activity',['Sensitivity','Method',r'\shortstack{Grouped\\MSE}','$R^2$',r'\shortstack{APP\\MSE}','$R^2$',r'\shortstack{S7\\MSE}','$R^2$'],rr,'Row-weighted ensemble prediction on matched rows within each sensitivity. Grouped source-excluded coverage is 1,003; primary-assay coverage is 2,607. APP/S7 retain 1,839/118 observations and unchanged measured targets. Source-excluded deployment reuses verified outer0 checkpoints and their original selection grid. Constants use training rows; oracle means use evaluation labels only as diagnostics.','tab:activity','llrrrrrr')
rr=[]
for c in ['APP','Davis_S7']:
 for m in methods:
  q=oldact[(oldact.cohort==c)&(oldact.method==m)&(oldact.weighting=='rows')].iloc[0];rr.append([cn[c] if m==methods[0] else '',names[m],fmt(q.mse),fmt(q.r_squared,3),fmt(q.correlation,3),fmt(q.prediction_sd),fmt(q.mean_bias)])
repl['CALIBRATION_TABLE']=table('historical_calibration',['Cohort','Original v3 model','MSE','$R^2$','$r$','Pred. SD','Bias'],rr,'Unchanged released-trained activity results on original response scales. Label SD is 0.4028 on APP and 0.3000 on S7. Small prediction spread and mean calibration are distinct from discrimination. No test-fitted adjustment is used.','tab:calibration','llrrrrr')
ranks=read(RUN/'ranking/campaign_summary.csv');rr=[]
for ca,phase,arm in [('v1','original','original'),('v2','factorial','R1C1'),('v2','deployment','all_ENsi_inner_selected'),('v3','group_selected',None)]:
 q=ranks[(ranks.campaign==ca)&(ranks.phase==phase)];q=q if arm is None else q[q.arm==arm]; mm=methods if ca=='v3' else ['gnn','nongraph','chemistry_tree','token_cnn'];rr.append([ca+' '+phase.replace('_',' ')]+[fmt(q[q.model==m].common_top5_percentile.iloc[0],3) for m in mm])
rnew=read(RUN/'evaluation/ranking_pools.csv')
for s in scenarios:rr.append([sn[s]]+[fmt(rnew[(rnew.scenario==s)&(rnew.method==m)].top5_mean_percentile.mean(),3) for m in methods])
repl['RANKING_TABLE']=table('ranking_history',['Protocol','GNN','No-message','Tree','CNN'],rr,'Common tie-aware APP mean top-five percentile across identical seventeen pools; larger is preferred. Historical v1/v2 ensembles have three actual seeds, v3 and sensitivities ten. Original reported values and every alternative historical arm remain in the evidence archive.','tab:ranking','lrrrr')
rr=[]
for m,label in [('activity_gnn','Activity GNN'),('pair_gnn','Pair GNN'),('pair_no_message','Pair no-message'),('pair_tree','Pair tree'),('pair_ridge','Pair ridge'),('zero','Zero'),('training_mean','Training mean')]:
 q=oldpair[(oldpair.method==m)&(oldpair.weighting=='rows')].iloc[0];eq=oldpair[(oldpair.method==m)&(oldpair.weighting=='equal_sequence')].iloc[0];rr.append([label,fmt(q.mse),fmt(eq.mse),fmt(q.correlation,3),fmt(q.prediction_sd,3),f'{int(q.non_tie_correct)}/{int(q.non_ties)}'])
repl['B2_TABLE']=table('B2_unchanged',['Method','Row MSE','Comp. MSE','$r$','Pred. SD','Sign'],rr,'Unchanged held-out measured B2 errors. Comp. gives each of twelve sequence components equal mass. Sign counts use 137 non-ties under the operational 0.02 band; the majority-direction control also achieves 111/137. Its regression error is undefined.','tab:pairranking','lrrrrr')
b3=read(RUN/'visibility/B3_stratified_metrics.csv');forced=read(RUN/'visibility/B3_rectangle_visibility.csv');rr=[]
for m in methods:
 q=b3[(b3.method==m)&(b3.stratum=='all')&(b3.predictor=='model')].iloc[0];ic=counts[(counts.case=='B3_original_source_heldout')&(counts.method==m)].iloc[0];f=forced[forced.method==m];rr.append([names[m],int(ic.exact_inputs),int(f.forced_exact_input_cancellation.sum()),int(f.forced_routed_real_arithmetic_cancellation.sum()),fmt(q.mse,6),fmt(q.prediction_sd,4)])
rr.append(['Zero interaction','--','--','--',fmt(b3[(b3.predictor=='zero')&(b3.stratum=='all')].mse.iloc[0],6),'0'])
repl['B3_TABLE']=table('B3_visibility',['Method','Inputs','Forced','Routed','MSE','Pred. SD'],rr,'Primary B3: distinct exact inputs among 165 states, forced-cancellation rectangles among 140, and separate algebraic routing cancellations. Model scores use unchanged source-held-out ensemble predictions. Counts concern finite observed inputs; routed equality is qualified by floating-point reassociation.','tab:B3','lrrrrr')
# Narrow, mechanically generated result prose is reviewed after execution.
prose=[]
for s in scenarios:
 q=activity[(activity.scenario==s)&(activity.cohort=='ENsiRNA_grouped')&(activity.weighting=='rows')].set_index('method');best=q.loc[methods,'mse'].idxmin();g=q.loc['corrected_gnn'];d=direct[(direct.scenario==s)&(direct.cohort=='ENsiRNA_grouped')&(direct.b=='corrected_no_message')&(direct.weighting=='rows')].iloc[0]
 prose.append(f"{sn[s]} grouped GNN MSE is {g.mse:.6f} ($R^2={g.r_squared:.3f}$); the lowest observed error among the four methods is {names[best]} at {q.loc[best,'mse']:.6f}. GNN minus no-message is {d.difference:+.6f} [{d.lower:+.6f},{d.upper:+.6f}] under conditional sequence resampling.")

for s in scenarios:
 d=direct[(direct.scenario==s)&(direct.cohort=='Davis_S7')&(direct.b=='corrected_no_message')&(direct.weighting=='rows')].iloc[0]
 prose.append(f"On S7, {sn[s].lower()} GNN minus no-message is {d.difference:+.6f} [{d.lower:+.6f},{d.upper:+.6f}]. Both sensitivities favor GNN in this fixed-ensemble contrast, while tree and CNN have lower observed MSE.")
repl['FOCUSED_TEXT']='\n\n'.join(prose)
repl['ABSTRACT_NEW']=f"The focused study completes {len(fit):,} new fits and reuses forty unaffected final checkpoints. The tree leads grouped prediction; GNN improves over no-message conditionally on S7, without beating the strongest observed baseline."
repl['CONCLUSION_NEW']='Both focused S7 comparisons favor GNN over matched no-message conditionally, while tree and CNN retain lower observed MSE. Grouped source-excluded prediction favors the tree. These target-specific findings do not establish a general GNN advantage or erase the historical source limitations.'
bc=counts[(counts.case=='B3_original_source_heldout')&(counts.method=='corrected_gnn')].iloc[0];fc=forced[forced.method=='corrected_gnn'];ss=counts[(counts.case=='S7_original_deployment')&(counts.method=='corrected_gnn')].iloc[0]
repl['VISIBILITY_TEXT']=f"The frozen GNN has {int(bc.exact_inputs)} exact B3 inputs among 165 states, with {int(fc.forced_exact_input_cancellation.sum())} literal forced cancellations and {int(fc.forced_routed_real_arithmetic_cancellation.sum())} algebraic forced-cancellation flags. Its original S7 deployment has {int(ss.exact_inputs)} exact inputs among 118 rows. Chemistry-name masks and complete input identity are audited separately; distinctness alone does not establish accurate chemistry prediction."
# Three caption-free numerical figures; all marks derive from named exports.
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
colors=['#65a89b','#9d88be','#d39773','#699fca'];fig,axes=plt.subplots(1,2,figsize=(7,2.35),layout='constrained')
for ax,s in zip(axes,scenarios):
 for j,m in enumerate(methods):
  z=seeds[(seeds.scenario==s)&(seeds.cohort=='ENsiRNA_grouped')&(seeds.method==m)&(seeds.weighting=='rows')];ax.scatter(np.full(len(z),j)+np.linspace(-.12,.12,len(z)),z.mse,s=10,c=colors[j],alpha=.65)
  v=activity[(activity.scenario==s)&(activity.cohort=='ENsiRNA_grouped')&(activity.method==m)&(activity.weighting=='rows')].mse.iloc[0];ax.scatter(j,v,marker='D',c='black',s=23,zorder=5)
 ref=activity[(activity.scenario==s)&(activity.cohort=='ENsiRNA_grouped')&(activity.method=='training_row_mean')&(activity.weighting=='rows')].mse.iloc[0];ax.axhline(ref,ls='--',c='#666',lw=.9,label='Training row mean');ax.set_xticks(range(4),['GNN','No-msg','Tree','CNN']);ax.set_title(sn[s]);ax.set_ylabel('Grouped MSE');ax.legend(fontsize=7,loc='best')
fig.savefig(S/'figures/main_robustness.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(7,2.25),layout='constrained');q=cal[cal.cohort=='APP'].set_index('method');a=q.loc['chemistry_tree'];b=q.loc['corrected_gnn'];vals=[a.squared_mean_bias-b.squared_mean_bias,a.prediction_variance-b.prediction_variance,-2*(a.covariance-b.covariance)];axes[0].barh(['Squared bias','Prediction variance','−2 covariance'],vals,color=[colors[2],colors[0],colors[1]]);axes[0].axvline(0,c='k',lw=.6);axes[0].set_xlabel('APP tree − GNN MSE components');axes[0].ticklabel_format(axis='x',style='sci',scilimits=(-2,2))
oldcomp=read(V3/'tables/direct_comparisons.csv');chosen=oldcomp[(oldcomp.a=='corrected_gnn')&oldcomp.b.isin(['corrected_no_message','chemistry_tree'])&(oldcomp.weighting=='rows')&oldcomp.cohort.isin(['APP','Davis_S7'])];labs=[]
for j,(_,r) in enumerate(chosen.iterrows()):axes[1].errorbar(r.difference,j,xerr=[[r.difference-r.lower],[r.upper-r.difference]],fmt='o',c=colors[0],capsize=3);labs.append(cn[r.cohort]+' / '+('no-msg' if r.b=='corrected_no_message' else 'tree'))
axes[1].set_yticks(range(len(labs)),labs);axes[1].axvline(0,c='k',lw=.7);axes[1].set_xlabel('Original GNN − control MSE');fig.savefig(S/'figures/main_calibration.pdf');plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(7,2.35),layout='constrained');q=counts[counts.case=='B3_original_source_heldout'].set_index('method');axes[0].bar(range(4),[q.loc[m,'exact_inputs'] for m in methods],color=colors);axes[0].axhline(165,c='#666',ls='--',lw=.8);axes[0].set_xticks(range(4),['GNN','No-msg','Tree','CNN']);axes[0].set_ylabel('Distinct B3 encoded inputs');axes[0].set_ylim(0,175);q=forced[forced.method=='corrected_gnn'];axes[1].scatter(q.observed_interaction,q.predicted_interaction,s=14,alpha=.6,color=colors[0]);lo=min(q.observed_interaction.min(),q.predicted_interaction.min());hi=max(q.observed_interaction.max(),q.predicted_interaction.max());axes[1].plot([lo,hi],[lo,hi],c='#888',ls='--',lw=.8);axes[1].axhline(0,c='#888',lw=.7);axes[1].set_xlabel('Measured four-term contrast');axes[1].set_ylabel('Fixed GNN contrast');fig.savefig(S/'figures/main_visibility.pdf');plt.close(fig)
# Preserve the originally inlined material as authoring inputs outside the minimal source.
backup=HERE/'authoring/initial_inlined';backup.mkdir(exist_ok=True)
for name in ['main.tex','appendix.tex']:
 if not (backup/name).exists():shutil.copy2(S/name,backup/name)
oldmain=(backup/'main.tex').read_text();pre=oldmain.split('\\maketitle')[0];body=(HERE/'authoring/main_body.tex').read_text()
for key,val in repl.items():assert '@@'+key+'@@' in body;body=body.replace('@@'+key+'@@',val)
assert '@@' not in body
statements=oldmain[oldmain.index('\\section*{AI use statement}'):];statements=statements.replace('A separate source archive contains every dependency needed to compile this manuscript.','The editable source archive has exactly four entries. Official style dependencies and the external build procedure are supplied separately; a clean archive build is verified.')
text=pre+body+statements;text=text.replace('Do Graphs Learn the Chemistry? Controlled Evidence from siRNA Prediction','Do Graphs Learn the Chemistry? Source Integrity and Calibration in siRNA Prediction').replace('Do Graphs Learn the Chemistry?\\\\Controlled Evidence from siRNA Prediction','Do Graphs Learn the Chemistry?\\\\Source Integrity and Calibration in siRNA Prediction');(S/'main.tex').write_text(text)
write(O/'main_object_evidence.json',mapping);write(O/'generation.json',dict(fit_count=len(fit),fit_completion=completion,inputs={str(p.relative_to(ROOT)):sha(p) for p in [RUN/'evaluation/activity_metrics.csv',RUN/'evaluation/direct_comparisons.csv',RUN/'visibility/input_counts.csv',RUN/'review_checks/source_scores.csv',RUN/'ranking/campaign_summary.csv',V3/'tables/pair_metrics.csv']},figures=['main_robustness.pdf','main_calibration.pdf','main_visibility.pdf']))
import subprocess,sys
subprocess.run([sys.executable,str(HERE/'appendix_author.py')],check=True)
print('Generated complete main manuscript, appendix, figures and numerical tables.')
