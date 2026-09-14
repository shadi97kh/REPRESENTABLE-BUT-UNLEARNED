"""Freeze held-out predictions before any outcome scoring. No test target reads."""
import argparse,datetime,hashlib,json,pickle
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from scipy import sparse
from common import write_json
from models import build,tensors,predict as nn_predict
from fit_engine import tabular
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
assert not (r/'evaluation_results.json').exists(),'Cannot generate fresh predictions after evaluation'
assert not (r/'prediction_commit.json').exists(),'Predictions already committed'
torch.set_num_threads(2);torch.use_deterministic_algorithms(True);device=json.loads((r/'profile.json').read_text())['device']
if device.startswith('cuda'):torch.cuda.set_per_process_memory_fraction(.30,0)
manifest=json.loads((r/'graph_manifest.json').read_text());fits=json.loads((r/'final_fit_results.json').read_text());selection=json.loads((r/'frozen_selection.json').read_text())
protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'ensemble':'Arithmetic mean of the three prespecified initialization predictions for each stochastic method. Deterministic ridge uses its one fitted prediction. No seed chosen by test scores.','seed_variability':'Sample SD across optimization seeds, not biological uncertainty.','reference':selection['development_selected_strongest_baseline'],'test_target_access':False,'checkpoint_reload_validation_tolerance':2e-6}
write_json(r/'prediction_protocol.json',protocol)
records=[];reload_checks=[]
for split in ['validation','test_internal','test_APP']:
 graph=tensors(r/f'{split}_graphs.npz',device);ids=json.loads((r/f'{split}_graph_ids.json').read_text());tab=None;sparse_cache={}
 for fit in fits:
  kind=fit['kind'];path=r/fit['checkpoint'];assert hashlib.sha256(path.read_bytes()).hexdigest()==fit['checkpoint_sha256']
  if path.suffix=='.pt':
   state=torch.load(path,map_location=device,weights_only=False);model=build(kind,manifest,state['config']).to(device);model.load_state_dict(state['model']);p=nn_predict(model,graph)*state['target_std']+state['target_mean'];del model
  else:
   with path.open('rb') as f:state=pickle.load(f)
   if kind=='chemistry_tree':
    if tab is None:tab=tabular(r,split)
    xx=tab[:,state['columns']]
   else:
    suffix='guide' if kind=='guide_ridge' else 'pairwise'
    if suffix not in sparse_cache:sparse_cache[suffix]=sparse.load_npz(r/f'{split}_{suffix}.npz')
    xx=sparse_cache[suffix]
   p=state['model'].predict(xx)
  assert len(p)==len(ids) and np.isfinite(p).all()
  if split=='validation':
   old=np.load(path.parent/'train_validation_predictions.npz')['validation'];err=float(np.max(np.abs(p-old)));assert err<2e-6,(kind,err);reload_checks.append({'kind':kind,'seed':fit['seed'],'max_error':err})
  else:records.extend({'record_id':rid,'split':split,'model':kind,'seed':fit['seed'],'prediction':float(v)} for rid,v in zip(ids,p))
 del graph
frame=pd.DataFrame(records);frame.to_csv(r/'heldout_predictions_by_seed.csv',index=False,float_format='%.17g')
ensemble=frame.groupby(['record_id','split','model'],sort=True).prediction.agg(['mean','std','count']).reset_index().rename(columns={'mean':'prediction','std':'optimization_seed_sd','count':'number_of_fits'});ensemble['optimization_seed_sd']=ensemble.optimization_seed_sd.fillna(0);ensemble.to_csv(r/'heldout_predictions.csv',index=False,float_format='%.17g')
write_json(r/'checkpoint_reload_checks.json',reload_checks)
files=['heldout_predictions_by_seed.csv','heldout_predictions.csv','prediction_protocol.json','checkpoint_reload_checks.json'];code=Path(__file__).resolve().parent
commit={'committed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'heldout_outcomes_scored':False,'files':{name:hashlib.sha256((r/name).read_bytes()).hexdigest() for name in files},'fit_checkpoints':{f['checkpoint']:f['checkpoint_sha256'] for f in fits},'prediction_rows_by_seed':len(frame),'ensemble_prediction_rows':len(ensemble),'prediction_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'evaluation_code_sha256':hashlib.sha256((code/'evaluate.py').read_bytes()).hexdigest(),'fit_source_hashes_unchanged':all(hashlib.sha256((code/name).read_bytes()).hexdigest()==h for name,h in selection['sources'].items() if name in ['models.py','fit_engine.py','graphs.py','split.py','common.py','train.py'])}
assert commit['fit_source_hashes_unchanged'];write_json(r/'prediction_commit.json',commit);write_json(r/'predict.outputs.json',files+['prediction_commit.json']);print(json.dumps(commit,indent=2))
