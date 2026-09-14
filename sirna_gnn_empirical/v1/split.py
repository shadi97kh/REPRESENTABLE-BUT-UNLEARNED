import argparse,json,collections,random,hashlib,datetime
from pathlib import Path
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run;rows=read_jsonl(r/'verified_observations.jsonl')
study=DSU(len(rows));seqgroup=DSU(len(rows));sources={};kmers={}
for i,x in enumerate(rows):
 if x['source_family'] in sources:study.join(i,sources[x['source_family']])
 else:sources[x['source_family']]=i
 for strand in ['guide','passenger']:
  seq=x[strand]['sequence'].replace('T','U')
  for k in set(seq[j:j+13] for j in range(len(seq)-12)):
   if k in kmers:study.join(i,kmers[k]);seqgroup.join(i,kmers[k])
   else:kmers[k]=i
comps=collections.defaultdict(list);seqs=collections.defaultdict(list)
for i,x in enumerate(rows):comps[study.find(i)].append(i);seqs[seqgroup.find(i)].append(i)
cids={k:digest(sorted(rows[i]['record_id'] for i in v))[:16] for k,v in comps.items()};sids={k:digest(sorted(set(rows[i]['guide']['sequence'] for i in v)))[:16] for k,v in seqs.items()}
external=[k for k,v in comps.items() if any(rows[i]['dataset']=='CMsiRNAdb_APP' for i in v)];assert len(external)==1 and all(rows[i]['dataset']=='CMsiRNAdb_APP' for i in comps[external[0]])
internal=sorted((k for k in comps if k not in external),key=lambda k:(-len(comps[k]),cids[k]));train=internal[:2];rest=sorted(internal[2:],key=lambda k:cids[k]);random.Random(20260913).shuffle(rest);dev=rest[:2];test=rest[2:4];train+=rest[4:]
assignment={k:'train' for k in train}|{k:'validation' for k in dev}|{k:'test_internal' for k in test}|{k:'test_APP' for k in external}
manifest=[]
for i,x in enumerate(rows):
 x.update(split=assignment[study.find(i)],study_group=cids[study.find(i)],sequence_group=sids[seqgroup.find(i)])
 manifest.append({k:x[k] for k in ['record_id','split','study_group','sequence_group','source_family','dataset']})
for key in ['source_family','study_group','sequence_group']:
 seen=collections.defaultdict(set)
 for x in rows:seen[x[key]].add(x['split'])
 assert all(len(v)==1 for v in seen.values()),key
save_jsonl(r/'split_manifest.jsonl',manifest)
for split in ['train','validation','test_internal','test_APP']:
 group=[x for x in rows if x['split']==split];save_jsonl(r/f'{split}_observations.jsonl',group)
 # Prediction takes only label-free inputs. Outcomes remain in a separate evaluation file.
 save_jsonl(r/f'{split}_inputs.jsonl',[{k:v for k,v in x.items() if k in ['record_id','split','study_group','sequence_group','guide','passenger','pairing','pairing_fraction','dose_reported','dose_unit','target','cell','time_h','delivery']} for x in group])
summary={split:{'observations':sum(x['split']==split for x in rows),'study_components':len({x['study_group'] for x in rows if x['split']==split}),'sequence_components':len({x['sequence_group'] for x in rows if x['split']==split}),'sources':sorted({x['source_family'] for x in rows if x['split']==split})} for split in assignment.values()}
protocol={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'seed':20260913,'primary_task':'Source-family transfer: all primary-verified APP chemistry/assay records held out from fitting and development. Separately report internal held-out ENsiRNA source components. No unseen-target or independent-laboratory claim.','grouping':'Transitive union of complete source families and shared exact 13-mers on either strand, with T/U equivalence for grouping only; model retains T. All related APP patents share one source-family unit. Two largest ENsiRNA components in training; seeded shuffle of other components: two validation, two internal test, rest training.','summary':summary,'leakage':{'source_family':0,'study_component':0,'sequence_component':0},'prior_test_exposure':'Original ENsiRNA supplied split inspected in the recovered brief and again during identity auditing, never used for current model choice. APP raw outcomes seen only for provenance reconciliation before split. No candidate-model test scores inspected before final frozen predictions. Historical P3/P4 reporter datasets are not trained on here.','primary_model_selection':'Equal-study-component validation MSE; no APP labels in fitting, tuning, calibration or selection. Same grouping weights for all training objectives.','metrics':{'B1':'Unclipped activity fraction = percent inhibition /100. Report pooled and equal-assay MSE, Spearman, pairwise ordering, top-five mean midrank percentile (rank-1)/(n-1) on primary-verified pools n>=20; raw activity and top-five regret. Rank ties average; predicted ties resolved by record_id.','B2':'All pre-enumerated exact-sequence, same-table/assay/dose/conjugate chemistry pairs. Difference error, Spearman, sign accuracy at fixed 0.02 activity tie band. No mean-difference SE manufactured from reported SD. Cluster bootstrap preserves shared sequences across all doses and pairs.','B3':'Only complete measured rectangles of two distinct single-position chemistry changes on the same base sequences/conjugates/assay. Zero eligible at freeze.','uncertainty':'2000 paired bootstrap draws with seed 30260913 at sequence-component level within fixed APP cohort; reuse same groups for all models. Separately report source-component results. Conditional uncertainty does not establish new-study generalization; initialization seeds are not biological replicates.','practical_gain_threshold':None},'candidate_methods':['guide_ridge','guide_pairwise_ridge','chemistry_tree','token_cnn','nongraph','gnn','gnn_no_chemistry'],'stochastic_initialization_seeds':[1103,2207,3301],'novelty':'Conventional supervised chemistry-aware message passing; no distinct learning-rule or architectural novelty asserted.','theorem_coverage':False,'wet_lab_work':False}
write_json(r/'split_protocol.json',protocol)
outputs=['split_manifest.jsonl','split_protocol.json']+[f'{s}_{t}.jsonl' for s in ['train','validation','test_internal','test_APP'] for t in ['observations','inputs']];write_json(r/'split.outputs.json',outputs);print(json.dumps(protocol,indent=2))
