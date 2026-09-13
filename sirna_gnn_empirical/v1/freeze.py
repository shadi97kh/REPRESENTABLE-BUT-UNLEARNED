import argparse,collections,datetime,hashlib,json,shutil
from pathlib import Path
from common import write_json
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
records=json.loads((r/'development_results.json').read_text());kinds=sorted({x['kind'] for x in records});selected={kind:min((x for x in records if x['kind']==kind),key=lambda x:(x['validation_equal_study_mse'],x['config']['name'])) for kind in kinds};controls=[k for k in kinds if k not in ['gnn','gnn_no_chemistry']];reference=min(controls,key=lambda k:selected[k]['validation_equal_study_mse'])
code=Path(__file__).resolve().parent;snapshot=r/'frozen_source';snapshot.mkdir(exist_ok=True)
for p in code.glob('*.py'):shutil.copy2(p,snapshot/p.name)
shutil.copy2(code/'run_all.sh',snapshot/'run_all.sh')
result={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'selected_configs':{k:v['config'] for k,v in selected.items()},'development_scores':{k:v['validation_equal_study_mse'] for k,v in selected.items()},'development_selected_strongest_baseline':reference,'final_seeds':[1103,2207,3301],'selection_rule':'Minimum equal-study validation MSE; configuration-name lexical tie break. Final stochastic fits use independent seeds and the same fixed stopping protocol.','test_metrics_opened':False,'sources':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in snapshot.iterdir()}}
write_json(r/'frozen_selection.json',result);write_json(r/'freeze.outputs.json',['frozen_selection.json']+[str(p.relative_to(r)) for p in snapshot.iterdir()]);print(json.dumps(result,indent=2))
