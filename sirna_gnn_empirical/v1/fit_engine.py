"""Training-only fitting helpers. This module never opens held-out targets or predictions."""
import copy,datetime,hashlib,json,os,pickle,random,time
from pathlib import Path
import numpy as np
import torch
from scipy import sparse
from sklearn.linear_model import Ridge
from sklearn.ensemble import ExtraTreesRegressor
from common import write_json
from models import build,tensors,take,predict

def atomic_torch(path,obj):
 temp=path.with_suffix(path.suffix+'.tmp');torch.save(obj,temp);temp.replace(path)
def setup(r):
 torch.set_num_threads(2);device=json.loads((r/'profile.json').read_text())['device']
 if device.startswith('cuda'):torch.cuda.set_per_process_memory_fraction(.30,0)
 torch.use_deterministic_algorithms(True)
 manifest=json.loads((r/'graph_manifest.json').read_text());data={s:tensors(r/f'{s}_graphs.npz',device) for s in ['train','validation']};targets={s:dict(np.load(r/f'{s}_targets.npz')) for s in data}
 mean=float(np.average(targets['train']['y'],weights=targets['train']['weight']));std=float(np.sqrt(np.average((targets['train']['y']-mean)**2,weights=targets['train']['weight'])));std=max(std,.05)
 return device,manifest,data,targets,mean,std

def fit_neural(r,kind,config,seed,phase,env):
 device,manifest,data,targets,mean,std=env;out=r/'fits'/phase/f'{kind}-seed{seed}-{config["name"]}'
 out.mkdir(parents=True,exist_ok=True);record=out/'fit.json'
 if record.exists():return json.loads(record.read_text())
 random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
 if device.startswith('cuda'):torch.cuda.manual_seed_all(seed);torch.cuda.reset_peak_memory_stats()
 model=build(kind,manifest,config).to(device);optimizer=torch.optim.AdamW(model.parameters(),lr=config['lr'],weight_decay=config['weight_decay']);gen=torch.Generator(device='cpu').manual_seed(seed+10000)
 y=torch.tensor((targets['train']['y']-mean)/std,device=device);w=torch.tensor(targets['train']['weight'],device=device);history=[];updates=0;best=float('inf');best_epoch=-1;start_epoch=0;accum_wall=0;accum_cpu=0
 before_message={n:p.detach().clone() for n,p in model.named_parameters() if 'message_weights' in n};grad_max=0.
 if (out/'last.pt').exists():
  state=torch.load(out/'last.pt',map_location=device,weights_only=False);model.load_state_dict(state['model']);optimizer.load_state_dict(state['optimizer']);history=state['history'];updates=state['optimizer_updates'];best=state['best'];best_epoch=state['best_epoch'];start_epoch=state['epoch']+1;torch.set_rng_state(state['torch_rng'].cpu());gen.set_state(state['permutation_rng'].cpu())
  if device.startswith('cuda'):torch.cuda.set_rng_state_all([s.cpu() for s in state['cuda_rng']])
  accum_wall=state['wall_s'];accum_cpu=state['cpu_s'];grad_max=state['message_gradient_max']
 if device.startswith('cuda'):torch.cuda.synchronize()
 start=time.perf_counter();cpu=time.process_time()
 for epoch in range(start_epoch,config['max_epochs']):
  model.train();perm=torch.randperm(len(y),generator=gen);loss_sum=0.
  for j in range(0,len(y),config['batch_size']):
   idx=perm[j:j+config['batch_size']].to(device);optimizer.zero_grad(set_to_none=True);prediction=model(**take(data['train'],idx));loss=(w[idx]*(prediction-y[idx])**2).mean();assert torch.isfinite(loss),'nonfinite training loss';loss.backward()
   gradient=sum(float(p.grad.norm()) for n,p in model.named_parameters() if 'message_weights' in n and p.grad is not None);grad_max=max(grad_max,gradient);torch.nn.utils.clip_grad_norm_(model.parameters(),5.);optimizer.step();updates+=1;loss_sum+=float(loss.detach())*len(idx)
  valpred=predict(model,data['validation'],config['batch_size'])*std+mean;val=float(np.average((valpred-targets['validation']['y'])**2,weights=targets['validation']['weight']));assert np.isfinite(val)
  history.append({'epoch':epoch+1,'optimizer_updates':updates,'train_weighted_normalized_mse':loss_sum/len(y),'validation_equal_study_mse':val})
  if val<best-1e-8:
   best=val;best_epoch=epoch;atomic_torch(out/'best.pt',{'model':model.state_dict(),'optimizer':optimizer.state_dict(),'optimizer_updates':updates,'epoch':epoch+1,'kind':kind,'config':config,'seed':seed,'target_mean':mean,'target_std':std,'validation_mse':val})
  if device.startswith('cuda'):torch.cuda.synchronize()
  state={'model':model.state_dict(),'optimizer':optimizer.state_dict(),'optimizer_updates':updates,'epoch':epoch,'kind':kind,'config':config,'seed':seed,'history':history,'best':best,'best_epoch':best_epoch,'torch_rng':torch.get_rng_state(),'permutation_rng':gen.get_state(),'cuda_rng':torch.cuda.get_rng_state_all() if device.startswith('cuda') else [],'wall_s':accum_wall+time.perf_counter()-start,'cpu_s':accum_cpu+time.process_time()-cpu,'message_gradient_max':grad_max};atomic_torch(out/'last.pt',state)
  if (epoch+1)%10==0:print(phase,kind,seed,config['name'],'epoch',epoch+1,'best validation MSE',best,flush=True)
  if epoch+1>=config['min_epochs'] and epoch-best_epoch>=config['patience']:break
 totalwall=accum_wall+time.perf_counter()-start;totalcpu=accum_cpu+time.process_time()-cpu
 beststate=torch.load(out/'best.pt',map_location=device,weights_only=False);model.load_state_dict(beststate['model']);preds={s:(predict(model,data[s],config['batch_size'])*std+mean) for s in data};np.savez_compressed(out/'train_validation_predictions.npz',**preds)
 changed=sum(float((p-before_message[n]).abs().sum()) for n,p in model.named_parameters() if n in before_message)
 if kind!='token_cnn':assert grad_max>0 and changed>0
 result={'phase':phase,'kind':kind,'seed':seed,'config':config,'epochs_executed':len(history),'optimizer_updates':updates,'selected_epoch':best_epoch+1,'selected_checkpoint_optimizer_updates':beststate['optimizer_updates'],'validation_equal_study_mse':best,'train_equal_study_mse':float(np.average((preds['train']-targets['train']['y'])**2,weights=targets['train']['weight'])),'parameters':sum(p.numel() for p in model.parameters()),'message_gradient_max':grad_max,'message_parameter_absolute_change_from_initialization':changed,'wall_s':totalwall,'cpu_s':totalcpu,'gpu_active_wall_s':totalwall if device.startswith('cuda') else 0,'device':device,'cuda_peak_allocated_bytes':torch.cuda.max_memory_allocated() if device.startswith('cuda') else 0,'checkpoint':str((out/'best.pt').relative_to(r)),'last_checkpoint':str((out/'last.pt').relative_to(r)),'history':history,'checkpoint_sha256':hashlib.sha256((out/'best.pt').read_bytes()).hexdigest(),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 write_json(record,result);return result

def tabular(r,split):
 z=np.load(r/f'{split}_graphs.npz');return np.concatenate([z['X'].reshape(len(z['X']),-1),z['A'].sum(-1).reshape(len(z['X']),-1),z['C']],axis=1)

def fit_classical(r,kind,config,seed,phase):
 out=r/'fits'/phase/f'{kind}-seed{seed}-{config["name"]}';out.mkdir(parents=True,exist_ok=True)
 if (out/'fit.json').exists():return json.loads((out/'fit.json').read_text())
 targets={s:dict(np.load(r/f'{s}_targets.npz')) for s in ['train','validation']}
 start=time.perf_counter();cpu=time.process_time();columns=None
 if kind=='chemistry_tree':
  xx=tabular(r,'train');columns=np.flatnonzero(np.ptp(xx,axis=0)>0);X={'train':xx[:,columns],'validation':tabular(r,'validation')[:,columns]};model=ExtraTreesRegressor(n_estimators=config['n_estimators'],min_samples_leaf=config['min_samples_leaf'],max_features=config['max_features'],random_state=seed,n_jobs=2)
 else:
  suffix='guide' if kind=='guide_ridge' else 'pairwise';X={s:sparse.load_npz(r/f'{s}_{suffix}.npz') for s in targets};model=Ridge(alpha=config['alpha'],solver='lsqr',tol=1e-7,max_iter=3000)
 model.fit(X['train'],targets['train']['y'],sample_weight=targets['train']['weight']);preds={s:model.predict(X[s]) for s in X}
 assert all(np.isfinite(p).all() for p in preds.values())
 with (out/'model.pkl').open('wb') as f:pickle.dump({'model':model,'columns':columns,'kind':kind,'config':config,'seed':seed},f)
 np.savez_compressed(out/'train_validation_predictions.npz',**preds)
 result={'phase':phase,'kind':kind,'seed':seed,'config':config,'optimizer_updates':0,'optimization':'weighted least squares LSQR' if kind!='chemistry_tree' else 'weighted randomized regression trees','validation_equal_study_mse':float(np.average((preds['validation']-targets['validation']['y'])**2,weights=targets['validation']['weight'])),'train_equal_study_mse':float(np.average((preds['train']-targets['train']['y'])**2,weights=targets['train']['weight'])),'wall_s':time.perf_counter()-start,'cpu_s':time.process_time()-cpu,'gpu_active_wall_s':0,'device':'cpu','features':X['train'].shape[1],'checkpoint':str((out/'model.pkl').relative_to(r)),'checkpoint_sha256':hashlib.sha256((out/'model.pkl').read_bytes()).hexdigest(),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 write_json(out/'fit.json',result);print(phase,kind,seed,config['name'],'validation MSE',result['validation_equal_study_mse'],flush=True);return result
