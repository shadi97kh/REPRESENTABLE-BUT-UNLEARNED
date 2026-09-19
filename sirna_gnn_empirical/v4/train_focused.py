"""Only affected grouped dependencies; frozen 792-fit maximum, no test selection."""
from engine import *
import copy
protocol=json.loads((RUN/'protocol.json').read_text());assert sha(RUN/'protocol.json')==json.loads((RUN/'protocol_freeze.json').read_text())['protocol_sha256']
original_rows,raw=load();curated=pd.read_csv(RUN/'adjudication/primary_reanchored_targets.csv');targets=dict(zip(curated.record_id,curated.primary_activity));records=[];predictions=[];selections=[];membership=[]
for pop in protocol['populations']:
 scenario=pop['scenario'];rows=copy.deepcopy(original_rows)
 if scenario=='primary_reanchored':
  for x in rows:
   if x['record_id'] in targets:x['activity']=float(targets[x['record_id']])
 eligible={i for i,x in enumerate(rows) if x['dataset']=='ENsiRNA' and (x['source_family']!='19282453' or scenario=='primary_reanchored' and x['record_id'] in targets)}
 def idx(groups):return np.array([i for i in sorted(eligible) if rows[i]['study_group'] in groups],int)
 test=idx(pop['test_groups']) if pop['name']!='deployment' else np.array([i for i,x in enumerate(rows) if x.get('split')=='test_APP' or x['record_id'].startswith('davis:')],int)
 for method in protocol['methods']:
  grid=protocol['tree_configurations'] if method=='chemistry_tree' else protocol['neural_configurations'];dev=[]
  for j,fold in enumerate(pop['inner']):
   tr,va=idx(fold['train_groups']),idx(fold['validation_groups']);assert len(tr)>0 and len(va)>0
   assert not(set(tr)&set(va) or set(tr)&set(test) or set(va)&set(test))
   if scenario=='source_excluded':assert all(rows[i]['source_family']!='19282453' for i in np.r_[tr,va,test])
   for cfg in grid:
    for seed in DEVSEEDS:
     name=f'{scenario}/{pop["name"]}/development/{method}/i{j}/{cfg["name"]}/s{seed}';rec,_=fit(name,method,cfg,seed,rows,raw,tr,va,np.array([],int));dev.append(rec);records.append(rec)
  scores={c['name']:float(np.mean([r['validation_score'] if r['status']=='completed' else float('inf') for r in dev if r['config']['name']==c['name']])) for c in grid};chosen=min(grid,key=lambda c:(scores[c['name']],c['name']));selections.append(dict(scenario=scenario,population=pop['name'],method=method,config=chosen,scores={k:v if np.isfinite(v) else None for k,v in scores.items()}));write(RUN/'selection.json',selections)
  f=pop['inner'][pop['final_validation_fold']];tr,va=idx(f['train_groups']),idx(f['validation_groups'])
  for role,ii in [('train',tr),('validation',va),('test',test)]:
   for i in ii:membership.append(dict(scenario=scenario,population=pop['name'],method=method,role=role,record_id=rows[i]['record_id'],source=rows[i]['source_family'],study_group=rows[i]['study_group'],sequence_group=rows[i]['sequence_group']))
  for seed in SEEDS:
   name=f'{scenario}/{pop["name"]}/final/{method}/s{seed}';rec,pr=fit(name,method,chosen,seed,rows,raw,tr,va,test);records.append(rec)
   for i,z in zip(test,pr['test']):
    x=rows[i];predictions.append(dict(scenario=scenario,population=pop['name'],cohort='APP' if x.get('split')=='test_APP' else 'Davis_S7' if x['record_id'].startswith('davis:') else 'ENsiRNA_grouped',record_id=x['record_id'],source=x['source_family'],study_group=x['study_group'],sequence_group=x['sequence_group'],method=method,seed=seed,prediction=float(z),activity=x['activity'],status=rec['status'],fit=name,training_row_mean=rec['training_row_mean'],training_equal_group_mean=rec['training_equal_group_mean']))
   export(pd.DataFrame(predictions),'focused_seed_predictions.csv')
  export(pd.DataFrame(membership),'focused_membership.csv')
  export(pd.DataFrame([{k:v for k,v in r.items() if k not in ['history','signature','hashes']} for r in records]),'fit_summary.csv')
assert len(records)<=protocol['maximum_new_fits'];write(RUN/'fit_completion.json',dict(actual_fits=len(records),failed_fits=sum(r['status']!='completed' for r in records),executed_updates=sum(r['executed_updates'] for r in records),selected_updates=sum(r['selected_updates'] for r in records),fit_cpu_s=sum(r['cpu_s'] for r in records),fit_wall_s=sum(r['wall_s'] for r in records),gpu_region_wall_s=sum(r['gpu_region_wall_s'] for r in records),fitting_scope='432 development and 360 final fits; four methods and ten final seeds; no B2 refits, no B3 label training'))
