"""Fresh retrospective follow-up. Archived v1 is read-only; no pretrained external models."""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import json,hashlib,time,random,pickle,math,datetime
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
ROOT=Path(__file__).resolve().parents[2]
OLD=ROOT/'runs/sirna_gnn_empirical/v1-20260913T185027Z'
RUN=ROOT/(Path(__file__).parent/'active_run.txt').read_text().strip()
SEEDS=[1103,2207,3301]
SPLITS=['train','validation','test_internal','test_APP']
CONFIGS=[dict(name='lr001_wd0001',lr=.001,weight_decay=.0001),dict(name='lr003_wd01',lr=.003,weight_decay=.01)]
TREE_CONFIGS=[dict(name='leaf2',min_samples_leaf=2),dict(name='leaf8',min_samples_leaf=8)]
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x),allow_nan=False)+'\n');q.replace(p)
def readlines(p):return [json.loads(s) for s in Path(p).read_text().splitlines()]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def setup():
 torch.set_num_threads(2);torch.use_deterministic_algorithms(True)
 available=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,set(available[:2]))
 if torch.cuda.is_available():torch.cuda.set_per_process_memory_fraction(.30,0)
 return 'cuda:0' if torch.cuda.is_available() else 'cpu'
def load():
 rows=[];zs=[]
 for s in SPLITS:
  rr=readlines(OLD/f'{s}_observations.jsonl');ids=json.loads((OLD/f'{s}_graph_ids.json').read_text());assert [x['record_id'] for x in rr]==ids
  rows+=rr;zs.append(dict(np.load(OLD/f'{s}_graphs.npz')))
 data={k:np.concatenate([z[k] for z in zs]) for k in ['X','A','mask','C']}
 return rows,data
def weights(groups):
 g=np.array(groups);u,c=np.unique(g,return_counts=True);d=dict(zip(u,c));return np.array([len(g)/(len(u)*d[x]) for x in g],np.float32)
def rawcontext(rows):
 return np.array([[math.log1p(x['dose_reported']),float(x['dose_unit']!='nM'),(x['time_h'] or 0)/24,float(x['time_h'] is None),float(x['cell'] is not None),x['pairing_fraction'],len(x['guide']['sequence'])/32,len(x['passenger']['sequence'])/32] for x in rows],np.float32)
def transform(data,rows,train,R,C,original_frozen=False):
 out=dict(data);X=data['X'];valid=X[train][data['mask'][train].astype(bool)];zero=(valid==0).all(0);relations=(data['A'][train].sum((0,2,3))>0)
 if R:
  out['X']=X.copy();out['X'][:,:,zero]=0
 raw=rawcontext(rows)
 if C:
  cc=raw.copy();verified=np.array([x['dose_unit']=='nM' for x in rows]);cc[:,0]=np.where(verified,cc[:,0]/math.log(101),0)
  # A physically verified dose can activate only if some verified nonconstant dose exists in training.
  if not verified[train].any():cc[:,0]=0
  # All other physical coordinates are fractions/flags, not divided by small fitted SDs.
  constants=np.ptp(cc[train],axis=0)==0;cc[:,constants]=0;out['C']=cc
 elif original_frozen:out['C']=data['C']
 else:
  mu=raw[train].mean(0);sd=raw[train].std(0);sd[sd<1e-6]=1;out['C']=(raw-mu)/sd
 support={'zero_node_columns':np.flatnonzero(zero).tolist(),'relation_support':relations.tolist(),'R':R,'C':C,'context_constant_columns':np.flatnonzero(np.ptp(out['C'][train],axis=0)==0).tolist(),'train_record_ids':[rows[i]['record_id'] for i in train]}
 return out,support
class GraphModel(nn.Module):
 def __init__(self,kind='gnn',R=False,support=None):
  super().__init__();self.kind=kind;self.R=R;self.embed=nn.Linear(93,32);self.self_layers=nn.ModuleList([nn.Linear(32,32) for _ in range(2)]);self.norms=nn.ModuleList([nn.LayerNorm(32) for _ in range(2)]);self.message_weights=nn.ParameterList([nn.Parameter(torch.empty(8,32,32)) for _ in range(2)])
  for w in self.message_weights:nn.init.xavier_uniform_(w)
  self.dropout=nn.Dropout(.1);self.readout=nn.Sequential(nn.Linear(64*32+8,64),nn.ReLU(),nn.Dropout(.1),nn.Linear(64,1));self.support=support or [True]*8
 def forward(self,X,A,mask,C):
  if self.kind=='gnn_no_chemistry':
   X=X.clone();X[:,:,39:87]=0;X[:,:,89:93]=0
   A=A.clone();A[:,2]+=A[:,0]+A[:,1];A[:,5]+=A[:,3]+A[:,4];A[:,[0,1,3,4]]=0
  h=torch.relu(self.embed(X))*mask[:,:,None];degree=A.sum((1,3)).clamp(min=1)[:,:,None]
  for own,w,norm in zip(self.self_layers,self.message_weights,self.norms):
   if self.R:
    w=torch.stack([w[2]+w[0]*self.support[0],w[2]+w[1]*self.support[1],w[2],w[5]+w[3]*self.support[3],w[5]+w[4]*self.support[4],w[5],w[6]*self.support[6],w[7]*self.support[7]])
   agg=A.sum(-1)[:,:,:,None]*h[:,None,:,:] if self.kind=='nongraph' else torch.matmul(A,h[:,None,:,:])
   message=torch.einsum('brnh,rhk->bnk',agg,w)/degree
   h=norm(h+self.dropout(torch.relu(own(h)+message)))*mask[:,:,None]
  return self.readout(torch.cat([h.flatten(1),C],1)).squeeze(-1)
class TokenCNN(nn.Module):
 def __init__(self):
  super().__init__();self.embed=nn.Linear(93,32);self.layers=nn.Sequential(nn.Conv1d(32,32,3,padding=1),nn.ReLU(),nn.Dropout(.1),nn.Conv1d(32,32,3,padding=1),nn.ReLU());self.readout=nn.Sequential(nn.Linear(2056,64),nn.ReLU(),nn.Dropout(.1),nn.Linear(64,1))
 def forward(self,X,A,mask,C):
  h=torch.relu(self.embed(X))*mask[:,:,None];h=(h+self.layers(h.transpose(1,2)).transpose(1,2))*mask[:,:,None];return self.readout(torch.cat([h.flatten(1),C],1)).squeeze(-1)
def model(kind,R,support):return TokenCNN() if kind=='token_cnn' else GraphModel(kind,R,support)
def take(data,idx):return {k:v[idx] for k,v in data.items()}
def pred(m,data,idx,device,mean=0.,std=1.):
 m.eval();out=[]
 with torch.no_grad():
  for j in range(0,len(idx),128):
   dd={k:torch.as_tensor(v[idx[j:j+128]],device=device) for k,v in data.items()};out.append(m(**dd).cpu().numpy()*std+mean)
 return np.concatenate(out)
def atomic_torch(p,x):
 tmp=p.with_suffix('.tmp');torch.save(x,tmp);tmp.replace(p)
def fit(tag,kind,R,C,config,seed,rows,raw,train,val,test,groupkey='study_group',pairs=None,lam=0.,fixed_epochs=None):
 out=RUN/'fits'/tag/f'{kind}-{config["name"]}-s{seed}-l{lam:g}';out.mkdir(parents=True,exist_ok=True)
 if (out/'fit.json').exists():return json.loads((out/'fit.json').read_text()),dict(np.load(out/'predictions.npz'))
 data,support=transform(raw,rows,train,R,C,original_frozen=tag.startswith('factorial'))
 write(out/'support.json',support);y=np.array([x['activity'] for x in rows],np.float32);wt=weights([rows[i][groupkey] for i in train]);wv=weights([rows[i][groupkey] for i in val]) if len(val) else np.array([])
 mean=float(np.average(y[train],weights=wt));std=max(float(np.sqrt(np.average((y[train]-mean)**2,weights=wt))),.05)
 start=time.perf_counter();cpu=time.process_time();device=setup();history=[];updates=0;best=float('inf');best_epoch=0;params=0;grad_max=0.;selectedupdates=0
 random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
 if device.startswith('cuda'):torch.cuda.manual_seed_all(seed);torch.cuda.reset_peak_memory_stats()
 if kind in ['chemistry_tree','context_ridge']:
  degree_features=data['A'].sum(-1).copy()
  if R:
   original_degree=degree_features.copy();degree_features[:,2]=original_degree[:,:3].sum(1);degree_features[:,5]=original_degree[:,3:6].sum(1)
   for relation in [0,1,3,4,6,7]:degree_features[:,relation]*=support['relation_support'][relation]
  xx=data['C'] if kind=='context_ridge' else np.concatenate([data['X'].reshape(len(rows),-1),degree_features.reshape(len(rows),-1),data['C']],1)
  columns=np.flatnonzero(np.ptp(xx[train],axis=0)>0)
  estimator=Ridge(alpha=1.) if kind=='context_ridge' else ExtraTreesRegressor(n_estimators=300,min_samples_leaf=config['min_samples_leaf'],max_features=.3,random_state=seed,n_jobs=2)
  estimator.fit(xx[train][:,columns],y[train],sample_weight=wt);pr={s:estimator.predict(xx[ii][:,columns]) for s,ii in [('validation',val),('test',test)] if len(ii)}
  with (out/'model.pkl').open('wb') as f:pickle.dump(dict(model=estimator,columns=columns,support=support,config=config),f)
  best=float(np.average((pr['validation']-y[val])**2,weights=wv)) if len(val) else 0.;checkpoint=out/'model.pkl';parameters=sum(t.tree_.node_count for t in estimator.estimators_) if kind=='chemistry_tree' else len(columns)+1
 else:
  m=model(kind,R,support['relation_support']).to(device);parameters=sum(p.numel() for p in m.parameters());opt=torch.optim.AdamW(m.parameters(),lr=config['lr'],weight_decay=config['weight_decay']);gen=torch.Generator().manual_seed(seed+10000)
  dd={k:torch.as_tensor(v[train],device=device) for k,v in data.items()};yt=torch.as_tensor((y[train]-mean)/std,device=device);ww=torch.as_tensor(wt,device=device)
  initial={n:p.detach().clone() for n,p in m.named_parameters()};start_epoch=0
  if (out/'last.pt').exists():
   st=torch.load(out/'last.pt',map_location=device,weights_only=False);m.load_state_dict(st['model']);opt.load_state_dict(st['optimizer']);history=st['history'];updates=st['optimizer_updates'];best=st['best'];best_epoch=st['best_epoch'];start_epoch=st['epoch'];gen.set_state(st['permutation_rng'].cpu());torch.set_rng_state(st['torch_rng'].cpu());torch.cuda.set_rng_state_all([s.cpu() for s in st['cuda_rng']]) if device.startswith('cuda') else None
  trainpairs=[];valpairs=[]
  if pairs:
   lookup={rows[i]['record_id']:i for i in range(len(rows))};local={i:j for j,i in enumerate(train)};ts=set(train);vs=set(val)
   for p in pairs:
    a,b=lookup[p['reference_id']],lookup[p['changed_id']]
    assert (a in ts)==(b in ts) and (a in vs)==(b in vs)
    if a in ts:trainpairs.append((local[a],local[b],p['observed_difference'],rows[a]['sequence_group']))
    if a in vs:valpairs.append((a,b,p['observed_difference'],rows[a]['sequence_group']))
  if trainpairs:
   pa=torch.tensor([z[0] for z in trainpairs],device=device);pb=torch.tensor([z[1] for z in trainpairs],device=device);dy=torch.tensor([z[2]/std for z in trainpairs],device=device);pw=torch.tensor(weights([z[3] for z in trainpairs]),device=device)
  max_epochs=fixed_epochs or 80
  for ep in range(start_epoch,max_epochs):
   m.train();perm=torch.randperm(len(train),generator=gen);loss_sum=0.
   # Pair task uses one full-endpoint batch per epoch so both terms have exact group-normalized weight.
   batch=len(train) if pairs else 128
   for j in range(0,len(train),batch):
    ix=perm[j:j+batch].to(device);opt.zero_grad(set_to_none=True)
    if pairs:
     pp=m(**dd);loss=(ww*(pp-yt)**2).mean()
     if lam:loss=loss+lam*(pw*(pp[pb]-pp[pa]-dy)**2).mean()
    else:pp=m(**take(dd,ix));loss=(ww[ix]*(pp-yt[ix])**2).mean()
    assert torch.isfinite(loss);loss.backward();grad_max=max(grad_max,sum(float(p.grad.detach().norm()) for p in m.parameters() if p.grad is not None));torch.nn.utils.clip_grad_norm_(m.parameters(),5.);opt.step();updates+=1;loss_sum+=float(loss.detach())*len(ix)
   vp=pred(m,data,val,device,mean,std) if len(val) else np.array([])
   if valpairs:
    d=dict(zip(val,vp));ee=np.array([(d[b]-d[a]-dy)**2 for a,b,dy,g in valpairs]);score=float(np.average(ee,weights=weights([z[3] for z in valpairs])))
   else:score=float(np.average((vp-y[val])**2,weights=wv)) if len(val) else 0.
   history.append(dict(epoch=ep+1,optimizer_updates=updates,training_normalized_loss=loss_sum/len(train),validation_score=score))
   if fixed_epochs or score<best-1e-8:
    best=score;best_epoch=ep+1;selectedupdates=updates;atomic_torch(out/'best.pt',dict(model=m.state_dict(),optimizer=opt.state_dict(),epoch=ep+1,optimizer_updates=updates,target_mean=mean,target_std=std,kind=kind,R=R,C=C,support=support,config=config,seed=seed,lambda_pair=lam))
   atomic_torch(out/'last.pt',dict(model=m.state_dict(),optimizer=opt.state_dict(),epoch=ep+1,optimizer_updates=updates,history=history,best=best,best_epoch=best_epoch,permutation_rng=gen.get_state(),torch_rng=torch.get_rng_state(),cuda_rng=torch.cuda.get_rng_state_all() if device.startswith('cuda') else []))
   if not fixed_epochs and ep+1>=10 and ep+1-best_epoch>=15:break
  st=torch.load(out/'best.pt',map_location=device,weights_only=False);m.load_state_dict(st['model']);selectedupdates=st['optimizer_updates'];pr={s:pred(m,data,ii,device,mean,std) for s,ii in [('validation',val),('test',test)] if len(ii)};checkpoint=out/'best.pt';delta=sum(float((p.detach()-initial[n]).abs().sum()) for n,p in m.named_parameters());assert grad_max>0 and delta>0
 del data
 result=dict(tag=tag,kind=kind,R=R,C=C,config=config,seed=seed,lambda_pair=lam,training_n=len(train),validation_n=len(val),test_n=len(test),parameters=parameters,parameter_unit='tree nodes' if kind=='chemistry_tree' else 'scalars',optimizer_updates=updates,selected_updates=selectedupdates,selected_epoch=best_epoch,epochs_executed=len(history),validation_score=best,history=history,checkpoint=str(checkpoint.relative_to(RUN)),checkpoint_sha256=sha(checkpoint),wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,device=device if kind not in ['chemistry_tree','context_ridge'] else 'cpu',max_gradient_sum=grad_max,peak_gpu_allocated_bytes=torch.cuda.max_memory_allocated() if device.startswith('cuda') else 0)
 np.savez_compressed(out/'predictions.npz',**pr,test_indices=test,validation_indices=val);write(out/'fit.json',result);print(tag,kind,config['name'],seed,'lambda',lam,'best',best_epoch,'updates',updates,'score',round(best,6),flush=True)
 return result,pr
