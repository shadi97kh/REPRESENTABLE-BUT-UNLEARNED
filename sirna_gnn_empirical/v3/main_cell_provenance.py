"""Closed, machine-readable row/field map for every main-table data cell."""
from common import *
P=ROOT/'papers/interaction_recoverability_iclr2027/v6';maps={Path(q['output']).stem:q for q in json.loads((RUN/'table_provenance.json').read_text())};cells=[]
def cell(table,row,column,path,selector,fields,operation='identity',digits=None):
 df=pd.read_csv(RUN/path);z=df
 for k,v in selector.items():z=z[z[k].astype(str)==str(v)]
 assert len(z)==1,(table,selector,len(z));q=z.iloc[0];values={f:None if pd.isna(q[f]) else q[f].item() if isinstance(q[f],np.generic) else q[f] for f in fields};cells.append(dict(table=table,row=row+1,column=column+1,display=maps[table]['rows'][row][column],source=path,source_sha256=sha(RUN/path),selector=selector,fields=values,operation=operation,decimal_places=digits))
methods=['original_gnn','corrected_gnn','original_no_message','corrected_no_message','chemistry_tree','token_cnn','guide_ridge','pairwise_ridge','training_row_mean','training_equal_group_mean','oracle_test_mean_diagnostic']
for i,m in enumerate(methods):
 for j,c in enumerate(['ENsiRNA_grouped','APP','Davis_S7']):
  for k,f in enumerate(['mse','r_squared']):cell('main_activity',i,1+2*j+k,'tables/activity_metrics.csv',dict(cohort=c,method=m,weighting='rows'),[f],digits=4 if f=='mse' else 3)
for i,c in enumerate(['ENsiRNA_grouped','APP','Davis_S7']):
 for j,m in enumerate(['corrected_no_message','chemistry_tree']):
  row=2*i+j;sel=dict(cohort=c,a='corrected_gnn',b=m,weighting='rows')
  cell('main_comparisons',row,2,'tables/direct_comparisons.csv',sel,['difference','lower','upper'],'point estimate and formatted interval on two lines',4)
  cell('main_comparisons',row,3,'tables/seed_comparison_stability.csv',sel,['mean','std'],'mean and parenthesized sample SD on two lines',4)
for i,m in enumerate(['activity_gnn','pair_gnn','pair_no_message','pair_tree','pair_ridge','zero','training_mean','training_majority']):
 for col,fields,w,digits,op in [(1,['mse'],'rows',4,'identity'),(2,['mse'],'equal_sequence',4,'identity'),(3,['correlation'],'rows',3,'identity'),(4,['prediction_sd'],'rows',3,'identity'),(5,['non_tie_correct','non_ties'],'rows',None,'formatted correct/non-tie counts')]:cell('main_pair_ranking',i,col,'tables/pair_metrics.csv',dict(method=m,weighting=w),fields,op,digits)
for i,m in enumerate(['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn','zero_interaction']):
 for j,f in enumerate(['mse','correlation','prediction_sd','mean_bias']):cell('main_B3',i,j+1,'tables/b3_metrics.csv',dict(method=m,seed='deterministic' if m=='zero_interaction' else 'ensemble'),[f],digits=4 if f=='mse' else 3)
for i,c in enumerate(['ENsiRNA_grouped','APP','Davis_S7']):
 for j,f in enumerate(['n','sequence_components','study_linked_groups','source_families']):cell('main_data',i,j+1,'tables/activity_metrics.csv',dict(cohort=c,method='corrected_gnn',weighting='rows'),[f],digits=0)
b3=json.loads((RUN/'b3_primary/evaluation_summary.json').read_text())
for j,f in enumerate(['conditions','sequence_backgrounds','source_families','source_families']):cells.append(dict(table='main_data',row=4,column=j+2,display=maps['main_data']['rows'][3][j+1],source='b3_primary/evaluation_summary.json',source_sha256=sha(RUN/'b3_primary/evaluation_summary.json'),json_field=f,value=b3[f],operation='single primary background/source; study-link count is one for this entire panel'))
for i,c in enumerate(['ENsiRNA_grouped','APP','Davis_S7']):
 for j,m in enumerate(['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn']):
  for k,f in enumerate(['correlation','prediction_sd','label_sd','mean_bias']):cell('main_calibration',4*i+j,k+2,'tables/activity_metrics.csv',dict(cohort=c,method=m,weighting='rows'),[f],digits=3 if f=='correlation' else 4)

opt=pd.read_csv(RUN/'tables/optimization_by_task.csv')
for i,q in opt.iterrows():
 for j,f in enumerate(['development_fits','final_fits','executed_updates','selected_updates']):cell('main_optimization',i,j+2,'tables/optimization_by_task.csv',dict(task=q.task,method=q.method),[f],digits=0)

# Non-numeric identifiers and headings have a closed literal mapping, not an unexplained data source.
covered={(q['table'],q['row'],q['column']) for q in cells}
for name in ['main_data','main_activity','main_comparisons','main_pair_ranking','main_B3','main_optimization','main_calibration']:
 for i,row in enumerate(maps[name]['rows']):
  for j,value in enumerate(row):
   if (name,i+1,j+1) not in covered:cells.append(dict(table=name,row=i+1,column=j+1,display=value,operation='literal method/cohort identifier; full identifier dictionaries in paper.py',generator_sha256=sha(HERE/'paper.py')))
write(RUN/'main_cell_provenance.json',dict(cells=cells,all_main_data_cells=len(cells),table_generator_sha256=sha(HERE/'paper.py'),mapping_generator_sha256=sha(Path(__file__)),metric_definitions='metrics.py and manuscript Appendix A5; weighting and selectors explicit per cell',rounding='Python fixed decimal formatting in paper.py; full-precision source CSV is authoritative; undefined shown as --'))
print('Mapped',len(cells),'main-table cells')
