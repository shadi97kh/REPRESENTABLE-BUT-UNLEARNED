import argparse,json
from pathlib import Path
from common import write_json
from fit_engine import setup,fit_neural,fit_classical
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run;selected=json.loads((r/'frozen_selection.json').read_text());records=[]
for kind in ['guide_ridge','guide_pairwise_ridge','chemistry_tree']:
 for seed in (selected['final_seeds'] if kind=='chemistry_tree' else [1103]):records.append(fit_classical(r,kind,selected['selected_configs'][kind],seed,'final'))
env=setup(r)
for kind in ['token_cnn','nongraph','gnn','gnn_no_chemistry']:
 for seed in selected['final_seeds']:records.append(fit_neural(r,kind,selected['selected_configs'][kind],seed,'final',env))
write_json(r/'final_fit_results.json',records);write_json(r/'train.outputs.json',['final_fit_results.json']+[str(p.relative_to(r)) for p in sorted((r/'fits/final').rglob('*')) if p.is_file()]);print('Final fits:',len(records),'optimizer updates:',sum(x['optimizer_updates'] for x in records))
