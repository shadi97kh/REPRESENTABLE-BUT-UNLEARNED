"""Inline complete retained methods and compact numerical summaries outside scientific fitting."""
from common import *
import re
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';P=ROOT/'papers/interaction_recoverability_iclr2027/v9';S=P/'source';backup=HERE/'authoring/initial_inlined';old=(backup/'appendix.tex').read_text();read=lambda p:pd.read_csv(p)
methods=['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn'];names={'corrected_gnn':'R1 GNN','corrected_no_message':'R1 no-msg','original_gnn':'R0 GNN','original_no_message':'R0 no-msg','chemistry_tree':'Tree','token_cnn':'CNN','guide_ridge':'Guide ridge','pairwise_ridge':'Pairwise ridge','training_row_mean':'Train row mean','training_equal_group_mean':'Train group mean','oracle_test_mean':'Oracle mean','oracle_test_mean_diagnostic':'Oracle mean'}
def esc(s):return str(s).replace('&',r'\&').replace('_',r'\_').replace('%',r'\%')
def fmt(v,n=4):return '--' if pd.isna(v) else f'{v:.{n}f}'
def tab(headers,rows,caption,label):
 cols=''.join('l' if j<2 else 'r' for j in range(len(headers)))
 h=' & '.join(headers)+r'\\\midrule'
 return '\n{\\small\n\\begin{longtable}{@{}'+cols+'@{}}\n\\toprule '+h+'\\endfirsthead\n\\toprule '+h+'\\endhead\n'+'\n'.join(' & '.join(map(str,r))+r'\\' for r in rows)+'\n\\bottomrule\n\\caption{'+caption+'}\\label{'+label+'}\\\\\n\\end{longtable}}\n'
def path(p):return r'\protect\path{'+p+'}'
# Replace obsolete main-object roadmap; definitions/proofs remain untouched below it.
start=old.index('\\section{Data, chemistry');a=old[:old.index('\\section{Appendix roadmap')]+r'''\section{Appendix roadmap and complete dependency map}\label{app:roadmap}
The source preserves complete finite-model definitions and proofs. Historical campaign v3 is explicitly distinguished from the focused v4 sensitivity. Supplement paths beginning \texttt{v3/} or \texttt{v4/} identify the corresponding run subtree in the indexed evidence archive. Raw records are external; no mathematical proof ends at a file pointer.
'''
maprows=[['Workflow (Fig.~\\ref{fig:workflow})','\\ref{app:architecture}, \\ref{app:training}','Forward pass, masks, losses, both derivatives'],['Source tables (Tables~\\ref{tab:data}, \\ref{tab:cohorts})','\\ref{app:data}, \\ref{app:robustness}','All sources, grouping, scoring roles'],['Adjudication (Table~\\ref{tab:adjudication})','\\ref{app:robustness}','Cells, identities, ambiguity, target rules'],['New fits (Table~\\ref{tab:optimization})','\\ref{app:training}, \\ref{app:robustness}','Frozen selection, reuse, updates'],['Prediction (Table~\\ref{tab:activity}; Fig.~\\ref{fig:robustness})','\\ref{app:metrics}, \\ref{app:focusedresults}','Matched rows, seeds, constants, uncertainty'],['Calibration (Table~\\ref{tab:calibration}; Fig.~\\ref{fig:calibration})','\\ref{app:metrics}, \\ref{app:robustness}','Exact decomposition, conditional contrasts'],['Ranking (Table~\\ref{tab:ranking})','\\ref{app:metrics}, \\ref{app:focusedresults}','Fractional ties, pools, protocol history'],['B2 (Table~\\ref{tab:pairranking})','\\ref{app:data}, \\ref{app:training}, \\ref{app:results}','Measured endpoints, losses, controls'],['B3 (Table~\\ref{tab:B3}; Fig.~\\ref{fig:visibility})','\\ref{app:b3}, \\ref{app:inputcollision}','States, cancellation proof, dependence']]
a+='\\begin{longtable}{p{.39\\textwidth}p{.16\\textwidth}p{.34\\textwidth}}\\toprule Main object & Appendix & Complete support\\\\\\midrule\n'+'\n'.join(' & '.join(r)+r'\\' for r in maprows)+'\n\\bottomrule\\caption{Main-object dependency map; detailed evidence identifiers appear in the final index.}\\\\\\end{longtable}\n'
a+=old[start:old.index('\\section{Complete current results')]
a=a.replace('The current experiment matrix enumerates','The historical v3 experiment matrix enumerates').replace('The exact selected configurations and every competing score appear in the result tables.','Every selected configuration and competing score is preserved in '+path('v3/tables/selection_scores.csv')+'.')
a=a.replace('AP remains unresolved and excluded.','AP was excluded from the original B3 admission; its primary sugar alias and duplicate release mappings are adjudicated in Appendix~\\ref{app:robustness}.')
# Avoid misrepresenting historical proof/data definitions as new-target guarantees.
a=a.replace('Primary versus released outcomes}', 'Original primary versus released comparison}')
a+='\n\\section{Preserved results and negative evidence}\\label{app:results}\nThe following compact tables report the frozen released-label campaign. They are not new fits. Cohorts, weighting and target definitions are retained. Focused source-excluded and primary-assay results follow in Appendix~\\ref{app:focusedresults}.\n'
q=read(V3/'tables/activity_metrics.csv');rr=[]
for (cohort,method),g in q[q.weighting=='rows'].groupby(['cohort','method']):
 r=g.iloc[0];rr.append([cohort.replace('ENsiRNA_grouped','Grouped').replace('Davis_S7','S7'),names.get(method,esc(method)),fmt(r.mse),fmt(r.r_squared,3),fmt(r.correlation,3),fmt(r.prediction_sd),fmt(r.mean_bias)])
a+=tab(['Cohort','Model','MSE','$R^2$','$r$','Pred. SD','Bias'],rr,'All original ensemble activity procedures with row weighting. Full alternate weighting and metrics: '+path('v3/tables/activity_metrics.csv')+'.','tab:oldactivity')
q=read(V3/'tables/pair_metrics.csv');rr=[]
for _,r in q.iterrows():rr.append([esc(r.method).replace('\\_',' '),r.weighting.replace('equal_sequence','equal seq.').replace('equal_background','equal background'),fmt(r.mse),fmt(r.correlation,3),fmt(r.prediction_sd),str(int(r.non_tie_correct))+'/'+str(int(r.non_ties))])
a+=tab(['B2 method','Weighting','MSE','$r$','Pred. SD','Signs'],rr,'All original B2 procedure/weighting combinations; zero and majority are distinct controls. Endpoint scores, within-assay comparisons and conditional contrasts remain in '+path('v3/tables/pair_endpoint_metrics.csv')+', '+path('v3/tables/B2_within_assay_variation.csv')+' and '+path('v3/tables/pair_comparisons.csv')+'.','tab:oldpairs')
# Keep every actual appendix diagnostic graphic, with its separate caption and provenance.
figures=re.findall(r'\\begin\{figure\}.*?\\end\{figure\}',old,re.S)
a+='\\subsection{Preserved diagnostic figures}\nThese plots use only the original v3 predictions and optimization records. They provide the evidence behind the historical comparisons; the new sensitivities do not overwrite them.\n'
for f in figures:
 f=f.replace('The complete seed distributions are retained','The original complete seed distributions are retained').replace('Small discrepancies can reflect rounding; larger ones require further curation.','Small discrepancies can reflect rounding; source adjudication and the separate primary-assay sensitivity are in Appendix~\\ref{app:robustness}.')
 a+=f+'\n'
a+='\\clearpage\n'
# Historical reproduction/resource section retained, with version scope corrected.
r=old[old.index('\\section{Reproduction, arithmetic'):].replace('\\label{appendix:end}','')
r=r.replace('\\section{Reproduction, arithmetic, resources and preservation}','\\section{Historical reproduction, arithmetic and resources}')
a+=r+'\n'+(HERE/'authoring/robustness_appendix.tex').read_text()
a+='\n\\section{Complete focused numerical summaries and evidence index}\\label{app:focusedresults}\nAll following tables derive from completed v4 fits or explicitly reused fixed predictions. Detailed rows are indexed separately; no unfavorable finite seed or failed attempt is removed.\n'
# New activity: full cohort/scenario, methods and constants, both weighting schemes.
q=read(RUN/'evaluation/activity_metrics.csv');rr=[]
for (scenario,weighting),g in q.groupby(['scenario','weighting']):
 a+='\\subsection{'+('Source-excluded' if scenario=='source_excluded' else 'Primary-assay')+' activity: '+esc(weighting).replace('\\_',' ')+'}\n';rr=[]
 for _,r in g.iterrows():rr.append([r.cohort.replace('ENsiRNA_grouped','Grouped').replace('Davis_S7','S7'),names[r.method],fmt(r.mse),fmt(r.r_squared,3),fmt(r.correlation,3),fmt(r.prediction_sd),fmt(r.mean_bias)])
 a+=tab(['Cohort','Method','MSE','$R^2$','$r$','Pred. SD','Bias'],rr,'Complete procedure metrics for this target and weighting. Full precision, MAE, label spread and means: '+path('v4/evaluation/activity_metrics.csv')+'.','tab:new:'+scenario+':'+weighting)
q=read(RUN/'evaluation/seed_stability.csv');q=q[(q.weighting=='rows')&q.method.isin(methods)];rr=[]
for _,r in q.iterrows():rr.append([('Excluded' if r.scenario=='source_excluded' else 'Primary')+'/'+r.cohort.replace('ENsiRNA_grouped','group').replace('Davis_S7','S7'),names[r.method],int(r.seed_count),fmt(r.min_seed_mse),fmt(r.mean_seed_mse),fmt(r.max_seed_mse),fmt(r.sd_seed_mse)])
a+=tab(['Target/cohort','Method','Seeds','Min','Mean','Max','SD'],rr,'All final-seed MSE summaries. Mean member risk is not ensemble risk; the exact decomposition is proved in Appendix~\\ref{app:metrics}. Full seed and fold records: '+path('v4/evaluation/seed_metrics.csv')+' and '+path('v4/fit_summary.csv')+'.','tab:newseeds')
q=read(RUN/'evaluation/direct_comparisons.csv');q=q[q.weighting=='rows'];rr=[]
for _,r in q.iterrows():rr.append([('Excluded' if r.scenario=='source_excluded' else 'Primary')+'/'+r.cohort.replace('ENsiRNA_grouped','group').replace('Davis_S7','S7'),names[r.b],fmt(r.difference,5),fmt(r.lower,5),fmt(r.upper,5),int(r.groups)])
a+=tab(['Target/cohort','GNN minus','Difference','Lower','Upper','Groups'],rr,'Fixed-ensemble row-weighted sequence-component comparisons, 10,000 draws. No multiplicity adjustment or independent-study coverage is asserted. Equal-group variants and separate seed contrasts: '+path('v4/evaluation/direct_comparisons.csv')+' and '+path('v4/evaluation/seed_paired_differences.csv')+'.','tab:newcontrasts')
q=read(RUN/'evaluation/matched_subset_rescoring.csv');rr=[]
for _,r in q.iterrows():rr.append([r.evaluation_target.replace('primary_reanchored','Primary').replace('released','Released'),names[r.method],int(r.n),fmt(r.mse),fmt(r.r_squared,3),fmt(r.correlation,3)])
a+=tab(['Evaluation target','Released-trained model','Rows','MSE','$R^2$','$r$'],rr,'Identical 2,607 rows rescored with unchanged released-trained predictions. This is distinct from the primary-assay refits above.','tab:matchedtargets')
q=read(RUN/'visibility/input_counts.csv');rr=[]
short={'B3_original_source_heldout':'B3 / held-out','S7_original_deployment':'S7 / original','S7_source_excluded_reused':'S7 / excluded','S7_primary_reanchored':'S7 / primary'}
for _,r in q.iterrows():rr.append([short[r.case],names[r.method],int(r.rows),int(r.exact_inputs),int(r.routed_real_arithmetic_inputs),int(r.rows_in_collision_classes),fmt(r.finite_sample_collision_MSE_floor,5)])
a+=tab(['Frozen preprocessor','Method','Rows','Inputs','Routed','Collided','Floor'],rr,'Every audited method and deployment. Collided counts rows in non-singleton exact-input classes. The floor is the finite empirical within-class variance, not a population bound. All ten preprocessors per method are checked.','tab:inputcounts')
q=read(RUN/'visibility/B3_leave_one_margin.csv');rr=[]
for (m,axis),g in q.groupby(['method','omitted_axis']):rr.append([names[m],axis,len(g),fmt(g.difference.min(),6),fmt(g.difference.median(),6),fmt(g.difference.max(),6)])
a+=tab(['Method','Omitted margin','Cases','Min $\\Delta$','Median','Max $\\Delta$'],rr,'Dependent leave-one-margin model-minus-zero MSE. All fourteen AS and ten SS omissions are included per method. These are descriptive sensitivities, not independent fitted folds.','tab:b3leave')
q=read(RUN/'visibility/B3_shared_reference_sensitivity.csv');rr=[]
for _,r in q.iterrows():rr.append([names[r.method],fmt(r.reference_shift,5),fmt(r.model_mse,6),fmt(r.zero_mse,6),fmt(r.model_minus_zero,6)])
a+=tab(['Method','Reference shift','Model MSE','Zero MSE','Difference'],rr,'Common-reference perturbations by one reported endpoint SD in either direction, with the unperturbed comparison. These shifts are not confidence limits for a biological mean.','tab:b3reference')
q=read(RUN/'ranking/campaign_summary.csv');rr=[]
for _,r in q.iterrows():rr.append([esc(r.campaign+'/'+r.phase).replace('\\_',' '),esc(r.model).replace('\\_',' '),esc(r.arm).replace('\\_',' '),int(r.minimum_members),fmt(r.historical_top5_percentile),fmt(r.common_top5_percentile)])
# use fixed widths for long protocol strings
rt=tab(['Campaign','Method','Arm','Seeds','Original','Common'],rr,'Every historical protocol and arm under a common exact-tie definition. Constant and chemistry-invariant historical scores can change when ordering-dependent tie handling is replaced; the original values are preserved. All 17 per-pool values and overlaps: '+path('v4/ranking/all_pool_comparisons.csv')+' and '+path('v4/ranking/pool_overlap.csv')+'.','tab:allranking').replace('@{}llrrrr@{}','@{}p{.14\\textwidth}p{.19\\textwidth}p{.19\\textwidth}rrr@{}');a+=rt
# All-source scores, every method including constants, remain indexed. Main table covers every source.
a+='\\subsection{Evidence index and complete row-level records}\nThe separately indexed scientific archive uses paths relative to its root. Every listed table retains its full precision, column schema and SHA256 in the archive manifest. Figures are generated from the same named records.\n'
index=[('v4/adjudication/row_reconciliation.csv','All 1,972 raw source rows (1,924 admitted): cells, chemistry, cardinality, original/primary values, categories'),('v4/adjudication/direct_cell_review.csv','Every 299 large discrepancy and 32 deterministic controls'),('v4/dataset/primary_reanchored_observations.jsonl','2,607 versioned observations, original values and graph indices'),('v4/review_checks/source_scores.csv','All-source, Bramsen, other-source, every individual source and constant; three weightings'),('v4/evaluation/per_source_metrics.csv','Every source after focused refitting; exact matched coverage'),('v4/evaluation/seed_predictions.csv','All new/reused final members and training-only constants'),('v4/evaluation/ensemble_predictions.csv','Exact member means and cohort identifiers'),('v4/fit_summary.csv','Every new development/final fit, seed, failure and optimizer count'),('v4/fits/','Best/last checkpoints, optimizer states, histories, selected transforms and memberships'),('v4/visibility/state_hashes.csv','Canonical state identities for every method and audited population'),('v4/visibility/chemical_feature_visibility.csv','Each observed chemistry name and frozen feature support'),('v4/visibility/B3_rectangle_visibility.csv','All method/rectangle identities and forced-cancellation flags'),('v4/visibility/B3_stratified_metrics.csv','Visible/forced strata versus zero-interaction controls'),('v4/ranking/all_pool_comparisons.csv','Every historical pool/arm, original and reconciled score'),('v4/evaluation/ranking_pools.csv','Both focused deployments and every matched APP pool'),('v3/tables/all_fit_summary.csv','All original 2,312 fits, preserved without new fitting'),('v3/tables/relation_gradients_by_fit.csv','Complete original arm/layer/seed gradient records'),('v3/tables/selected_parameter_activity.csv','Complete original parameter-level data-gradient audit'),('v3/tables/selection_scores.csv','Every original grouped candidate and selected configuration'),('v3/tables/per_assay_metrics.csv','All original per-assay efficacy comparisons'),('v3/b3_primary/','Original primary panel, rectangles, fixed predictions and source-exclusion checks')]
a+='\\begin{longtable}{p{.45\\textwidth}p{.45\\textwidth}}\\toprule Archive path & Preserved evidence\\\\\\midrule\n'+'\n'.join(path(p)+' & '+esc(desc)+r'\\' for p,desc in index)+'\n\\bottomrule\\caption{Indexed detailed evidence; no repetitive per-fit dump is printed in the manuscript.}\\\\\\end{longtable}\n'
a+='\\subsection{Current execution and source layout}\n'+r'''The new resumable entry point is \path{sirna_gnn_empirical/v4/run_all.sh}. Executed stages and their commands, hashes, wall time, child CPU time and maximum reported child RSS are in the evidence archive. The recorded fit matrix contains 792 actual fits; 40 previously fitted final models are reused. Checkpoint verification and inference are not additional fits. Original campaign costs remain separate. Administrative inspection, authoring, browsing, rendering, hashing and compilation have additional nonzero costs and are not GPU training measurements.

The source directory contains exactly \path{main.tex}, \path{appendix.tex}, \path{references.bib} and \path{figures/}. All main tables are inline, and the only manuscript inclusion is the complete appendix. Unmodified official ICLR style/bibliography files and the external build script are supplied separately. Build products, audits, source data and scientific code are outside this source directory. The delivered ZIP is independently compiled with those external dependencies. The original manuscripts, datasets, checkpoints, ledgers and unrelated proofs remain preserved.
'''
a+='\n\\label{appendix:end}\n'
f=json.loads((RUN/'fit_completion.json').read_text());stage=next(x for x in readlines(RUN/'commands.jsonl') if x['stage']=='fit')
resource_rows=[['New development / final fits','432 / 360'],['New failed fits',f['failed_fits']],['Executed neural updates',f['executed_updates']],['Selected-checkpoint updates',f['selected_updates']],['Summed fit CPU seconds',fmt(f['fit_cpu_s'],3)],['Summed fit wall seconds',fmt(f['fit_wall_s'],3)],['GPU fitting-region elapsed seconds',fmt(f['gpu_region_wall_s'],3)],['Fit-stage child CPU seconds',fmt(stage['child_cpu_s'],3)],['Fit-stage wall seconds',fmt(stage['wall_s'],3)],['Maximum fit-stage child RSS (KiB)',stage['max_rss_kib']]]
a=a.replace(r'\label{appendix:end}','')
a+=tab(['Focused campaign quantity','Recorded value'],resource_rows,'New fitting cost only; overlapping fit/stage/GPU regions are not added. GPU elapsed is not kernel time or energy. Evaluation, visibility, source checks and administrative work have additional nonzero costs, recorded in separate command receipts. Original costs above remain unchanged.','tab:newresources')
a+=r'\label{appendix:end}'
(S/'appendix.tex').write_text(a)

# Remove only unreferenced copied artwork in the fresh source; historical assets are intact.
refs=set(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{figures/([^}]+)\}',(S/'main.tex').read_text()+a))
for p in (S/'figures').iterdir():
 if p.name not in refs:p.unlink()
print('Wrote complete appendix;',len(refs),'referenced figure assets')
