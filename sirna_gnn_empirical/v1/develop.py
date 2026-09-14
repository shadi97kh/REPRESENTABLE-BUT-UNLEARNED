import argparse,json,datetime,math,hashlib
from pathlib import Path
from common import write_json
from fit_engine import setup,fit_neural,fit_classical
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
profile=json.loads((r/'profile.json').read_text());base={'hidden':32,'batch_size':128,'max_epochs':80,'min_epochs':10,'patience':15}
configs=[dict(base,name='lr001_wd0001',lr=.001,weight_decay=.0001),dict(base,name='lr003_wd01',lr=.003,weight_decay=.01)]
plan={'frozen_before_development_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'authority':'Current explicit user authorization; task-local conservative stop conditions do not reuse or enlarge any historical allowance.','profile_seconds_per_update':profile['seconds_per_update'],'neural_development_fits':8,'neural_final_fits':12,'maximum_optimizer_updates_excluding_profile':20*80*math.ceil(2626/128),'profile_extrapolated_neural_active_wall_s':20*80*math.ceil(2626/128)*profile['seconds_per_update'],'neural_configs':configs,'development_seed':701,'final_seeds':[1103,2207,3301],'cpu_cores':2,'gpu_devices':[0] if profile['device'].startswith('cuda') else [],'gpu_memory_fraction_cap':.30,'development_stage_wall_cap_s':900,'final_training_stage_wall_cap_s':1800,'maximum_campaign_child_cpu_s':7200,'maximum_campaign_aggregate_command_wall_s':3600,'unknown_historical_balance_reused':False,'cost_limits':'Profile extrapolation excludes parsing/imports/tree/ridge fitting. Actual process and stage costs recorded separately. No paid/cloud services.','input_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [r/'split_protocol.json',r/'graph_manifest.json']}}
if not (r/'resource_plan.json').exists():write_json(r/'resource_plan.json',plan)
else:plan=json.loads((r/'resource_plan.json').read_text());configs=plan['neural_configs']
records=[]
for kind in ['guide_ridge','guide_pairwise_ridge']:
 for alpha in [.1,1.,10.,100.]:records.append(fit_classical(r,kind,{'name':'alpha'+str(alpha),'alpha':alpha},701,'development'))
for leaf,mf in [(2,.3),(2,1.),(8,.3),(8,1.)]:records.append(fit_classical(r,'chemistry_tree',{'name':f'leaf{leaf}_features{mf}','n_estimators':300,'min_samples_leaf':leaf,'max_features':mf},701,'development'))
env=setup(r)
for kind in ['token_cnn','nongraph','gnn','gnn_no_chemistry']:
 for config in configs:records.append(fit_neural(r,kind,config,701,'development',env))
write_json(r/'development_results.json',records);outputs=['resource_plan.json','development_results.json']+[str(p.relative_to(r)) for p in sorted((r/'fits/development').rglob('*')) if p.is_file()];write_json(r/'develop.outputs.json',outputs)
print('Completed',len(records),'development fits; no held-out predictions generated.')
