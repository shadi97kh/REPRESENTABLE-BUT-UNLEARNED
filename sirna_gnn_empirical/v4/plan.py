from common import *
import collections,datetime
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';original=json.loads((V3/'protocol.json').read_text());rows=rows_all();cur=pd.read_csv(RUN/'adjudication/primary_reanchored_targets.csv');eligible=set(cur.record_id);br='59dde2a1f1d4bf70';methods=['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn']
def inner(groups,counts):
 b=[[] for _ in range(min(3,len(groups)))];n=[0]*len(b)
 for g in sorted(groups,key=lambda g:(-counts[g],g)):
  j=min(range(len(n)),key=lambda j:(n[j],j));b[j].append(g);n[j]+=counts[g]
 return [dict(train_groups=sorted(set(groups)-set(v)),validation_groups=v) for v in b]
tasks=[]
for scenario in ['source_excluded','primary_reanchored']:
 admitted=[x for x in rows if x['dataset']=='ENsiRNA' and (x['source_family']!='19282453' or scenario=='primary_reanchored' and x['record_id'] in eligible)];counts=collections.Counter(x['study_group'] for x in admitted)
 for p in original['activity_populations']:
  if p['name']=='outer0':continue # saved models exclude Bramsen; source scores/primary labels can be reevaluated without refitting
  if scenario=='source_excluded' and p['name']=='deployment':continue # full original outer0 design is an existing source-excluded deployment fit
  groups=set(counts)-set(p['test_groups']);ii=inner(groups,counts);fv=min(range(len(ii)),key=lambda j:(sum(counts[g] for g in ii[j]['validation_groups']),j));tasks.append(dict(scenario=scenario,name=p['name'],test_groups=p['test_groups'],inner=ii,final_validation_fold=fv,counts=dict(counts)))
configs=[dict(name='lr0.0001_wd0.0001',lr=.0001,weight_decay=.0001),dict(name='lr0.001_wd0.0001',lr=.001,weight_decay=.0001)];trees=[dict(name='leaf2',leaf=2),dict(name='leaf8',leaf=8)]
protocol=dict(version='robustness-v4',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),methods=methods,neural_configurations=configs,tree_configurations=trees,seeds=SEEDS,development_seeds=DEVSEEDS,populations=tasks,adjudication_sha256=sha(RUN/'adjudication/summary.json'),curated_targets_sha256=sha(RUN/'adjudication/primary_reanchored_targets.csv'),maximum_new_fits=len(tasks)*(4*3*2*2+4*10),development_fits=len(tasks)*48,final_fits=len(tasks)*40,optimizer=dict(epochs_max=500,epochs_min=50,patience=50,improvement_threshold=1e-8,batch_size=128,gradient_clip=5,dropout=.1,training_only_preprocessing=True),reuse=dict(outer0_final=40,outer0_development_selection='All original source-excluded outer0 selection and final validation retained for deployment/B3/curated source scoring. This inherits a broader six-neural-configuration historical grid; no new outcome-driven selection.',B2='All pair fits use APP endpoints only. Bramsen has no role; reuse unchanged.',B3='Existing outer0 source-held-out checkpoints; no retraining on panel labels.'),grouping='Original protected outer component assignments retained; inner groups balanced deterministically using eligible row counts only. No component split.',target='primary_reanchored is a named primary-assay sensitivity target on unique reported chemical identities; remaining 320 Bramsen rows quarantined. Not a certified replacement of the release or full molecule/assay provenance.',evaluation='Retrospective sensitivity throughout. APP/S7/B3 and prior released benchmarks were inspected historically; no independent confirmation claim.',resource_limits=dict(cpu_threads=2,gpu_memory_fraction=.30,gpu_device=0,parallel_fitting_workers=1,paid_resources=False),resource_estimate='792-fit upper bound, two candidate settings, two development seeds, three inner folds, ten final seeds. Approximate 1-3 hours at prior campaign throughput; actual cost recorded, no optional positive-result search.',no_push=True)
assert len(tasks)==9 and protocol['maximum_new_fits']==792
write(RUN/'protocol.json',protocol)
matrix=[]
for task in tasks:
 for method in methods:
  for j in range(3):
   for c in trees if method=='chemistry_tree' else configs:
    for seed in DEVSEEDS:matrix.append(dict(scenario=task['scenario'],population=task['name'],method=method,phase='development',inner=j,configuration=c['name'],seed=seed))
  for seed in SEEDS:matrix.append(dict(scenario=task['scenario'],population=task['name'],method=method,phase='final',inner=task['final_validation_fold'],configuration='inner_selected',seed=seed))
export(pd.DataFrame(matrix),'experiment_matrix.csv');write(RUN/'protocol_freeze.json',dict(protocol_sha256=sha(RUN/'protocol.json'),matrix_sha256=sha(RUN/'experiment_matrix.csv'),new_fits_before_freeze=0));print('Frozen',len(matrix),'new fits maximum;',protocol['development_fits'],'development;',protocol['final_fits'],'final; 40 unaffected final fits reused.')
