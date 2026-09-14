from common import *
from engine import *
rows,raw=load();pr=json.loads((RUN/'protocol.json').read_text());results=[];preds=[];selections=[];membership=[]
def idx(groups):return np.array([i for i,x in enumerate(rows) if x.get('study_group') in groups],int)
for pop in pr['activity_populations']:
 test=idx(pop['test_groups']) if pop['name']!='deployment' else np.array([i for i,x in enumerate(rows) if x.get('split')=='test_APP' or x['record_id'].startswith('davis:')],int)
 for method in pr['methods']+['chemistry_tree','guide_ridge','pairwise_ridge']:
  grid=pr['neural_configurations'] if method in pr['methods'] else pr['tree_configurations'] if method=='chemistry_tree' else pr['ridge_configurations'];dev=[]
  for j,fold in enumerate(pop['inner']):
   tr,va=idx(fold['train_groups']),idx(fold['validation_groups'])
   for cfg in grid:
    for seed in DEVSEEDS if method in pr['methods']+['chemistry_tree'] else [0]:
     name=f'activity/{pop["name"]}/dev/{method}/i{j}/{cfg["name"]}/s{seed}';r,p=fit(name,method,cfg,seed,rows,raw,tr,va,np.array([],int));dev.append(r)
  scores={c['name']:np.mean([r['validation_score'] if r['status']=='completed' else float('inf') for r in dev if r['config']['name']==c['name']]) for c in grid};chosen=min(grid,key=lambda c:(scores[c['name']],c['name']));sel=dict(population=pop['name'],method=method,config=chosen,inner_scores=scores);selections.append(sel);write(RUN/'activity_selection.json',selections)
  fold=pop['inner'][pop['final_validation_fold']];tr,va=idx(fold['train_groups']),idx(fold['validation_groups'])
  for role,ii in [('train',tr),('validation',va),('test',test)]:
   for i in ii:membership.append(dict(population=pop['name'],method=method,role=role,record_id=rows[i]['record_id'],study_group=rows[i].get('study_group'),sequence_group=rows[i]['sequence_group']))
  for seed in SEEDS if method in pr['methods']+['chemistry_tree'] else [0]:
   name=f'activity/{pop["name"]}/final/{method}/s{seed}';r,p=fit(name,method,chosen,seed,rows,raw,tr,va,test);results.append(r)
   for i,z in zip(test,p['test']):preds.append(dict(record_id=rows[i]['record_id'],population=pop['name'],cohort='Davis_S7' if rows[i]['record_id'].startswith('davis:') else 'APP' if rows[i].get('split')=='test_APP' else 'ENsiRNA_grouped',method=method,seed=seed,prediction=float(z),status=r['status'],fit=name,training_row_mean=r['training_row_mean'],training_equal_group_mean=r['training_equal_group_mean']))
  export(pd.DataFrame(preds),'activity_predictions.csv');export(pd.DataFrame(membership),'activity_membership.csv')
write(RUN/'activity_final_fits.json',results)
