"""Training-only preprocessing, grouped early stopping, immutable per-fit hashes."""
from common import *
import torch,random,pickle,math
from architecture import GraphModel,TokenCNN
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge

def setup():
 torch.set_num_threads(2);torch.use_deterministic_algorithms(True)
 available=sorted(os.sched_getaffinity(0));os.sched_setaffinity(0,set(available[:2]))
 if torch.cuda.is_available():torch.cuda.set_per_process_memory_fraction(.30,0)
 return 'cuda:0' if torch.cuda.is_available() else 'cpu'
def load():
 source_files=[OLD/f'{sp}_{part}' for sp in ['train','validation','test_internal','test_APP'] for part in ['observations.jsonl','graphs.npz']]+[PREV/'external_observations.jsonl',PREV/'external_graphs.npz']; baseline=json.loads((RUN/'preservation_before.json').read_text())
 for p in source_files:assert sha(p)==baseline[str(p.relative_to(ROOT))]['sha256'],p
 rows=rows_all();raws=[dict(np.load(OLD/f'{s}_graphs.npz')) for s in ['train','validation','test_internal','test_APP']];ext=readlines(PREV/'external_observations.jsonl');zz=dict(np.load(PREV/'external_graphs.npz'));zz['C']=np.zeros((len(ext),8),np.float32);raws.append(zz)
 return rows+ext,{k:np.concatenate([z[k] for z in raws]) for k in ['X','A','mask','C']}
def transform(raw,rows,tr):
 X=raw['X'];mask=raw['mask'];supported=np.any(X[tr][mask[tr].astype(bool)]!=0,axis=0);out=dict(raw);out['X']=X.copy();out['X'][:,:,~supported]=0
 cc=np.array([[math.log1p(x['dose_reported'])/math.log(101) if x['dose_unit']=='nM' else 0,float(x['dose_unit']!='nM'),(x['time_h'] or 0)/24,float(x['time_h'] is None),float(x['cell'] is not None),x['pairing_fraction'],len(x['guide']['sequence'])/32,len(x['passenger']['sequence'])/32] for x in rows],np.float32)
 varying=np.ptp(cc[tr],axis=0)>0;cc[:,~varying]=0;out['C']=cc
 counts=raw['A'][tr].sum((0,2,3));support=dict(node_support=supported.tolist(),context_varying=varying.tolist(),relation_support=(counts>0).tolist(),directed_relation_counts=counts.astype(int).tolist(),training_record_ids=[rows[i]['record_id'] for i in tr],support_threshold='strictly positive occurrence among active training nodes/edges; not adequate statistical support')
 return out,support

def network(method,support):return TokenCNN() if method=='token_cnn' else GraphModel('nongraph' if 'no_message' in method else 'gnn','corrected' in method,support['relation_support'])
def tp(data,idx,device):return {k:torch.as_tensor(v[idx],device=device) for k,v in data.items()}
def predict(m,data,idx,device,mu,sd):
 m.eval();zs=[]
 with torch.no_grad():
  for j in range(0,len(idx),128):zs.append(m(**tp(data,idx[j:j+128],device)).cpu().numpy()*sd+mu)
 return np.concatenate(zs) if zs else np.array([])
def savept(p,x):
 tmp=p.with_suffix('.tmp');torch.save(x,tmp);tmp.replace(p)
def classical_features(data,method):
 if method=='chemistry_tree':return np.concatenate([data['X'].reshape(len(data['X']),-1),data['A'].sum(-1).reshape(len(data['A']),-1),data['C']],1)
 bases=data['X'][:,:32,:5].reshape(len(data['X']),-1)
 if method=='pairwise_ridge':
  pairs=[(i,j) for i in range(24) for j in range(i+1,24)];xx=data['X'][:,:24,:5];pp=np.stack([(xx[:,i,:,None]*xx[:,j,None,:]).reshape(len(xx),25) for i,j in pairs],1).reshape(len(xx),-1);bases=np.concatenate([bases,pp],1)
 return np.concatenate([bases,data['C']],1)

def fit(name,method,cfg,seed,rows,raw,tr,va,te,groupkey='study_group',pairs=None,lam=0.,selection='activity'):
 out=RUN/'fits'/name;out.mkdir(parents=True,exist_ok=True)
 sig=dict(code={p:sha(HERE/p) for p in ['engine.py','architecture.py','common.py']},protocol=sha(RUN/'protocol.json'),input_manifest=sha(RUN/'preservation_before.json'),method=method,config=cfg,seed=seed,train=tr.tolist(),validation=va.tolist(),test=te.tolist(),lambda_pair=lam,selection=selection)
 if (out/'fit.json').exists():
  rec=json.loads((out/'fit.json').read_text());assert rec['signature']==sig
  for p,h in rec['hashes'].items():assert sha(out/p)==h
  return rec,dict(np.load(out/'predictions.npz'))
 data,sup=transform(raw,rows,tr);write(out/'preprocessing.json',sup);y=np.array([x['activity'] for x in rows],np.float32);w=weights([rows[i][groupkey] for i in tr]).astype(np.float32);wv=weights([rows[i][groupkey] for i in va]);mu=float(np.average(y[tr],weights=w));sd=max(float(np.sqrt(np.average((y[tr]-mu)**2,weights=w))),.05)
 start=time.perf_counter();cpu=time.process_time();device=setup();random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)
 if device.startswith('cuda'):torch.cuda.manual_seed_all(seed)
 history=[];updates=0;best=float('inf');bestep=0;parameters=0;status='completed';fail=None;maxgrad=0.;initial=None
 if method in ['chemistry_tree','guide_ridge','pairwise_ridge']:
  xx=classical_features(data,method);cols=np.flatnonzero(np.ptp(xx[tr],axis=0)>0);x=xx[:,cols];xm=x[tr].mean(0);xs=x[tr].std(0);xs[xs<1e-6]=1
  if method!='chemistry_tree':x=(x-xm)/xs
  m=ExtraTreesRegressor(n_estimators=300,min_samples_leaf=cfg['leaf'],max_features=.3,n_jobs=2,random_state=seed) if method=='chemistry_tree' else Ridge(alpha=cfg['alpha'],solver='lsqr')
  m.fit(x[tr],y[tr],sample_weight=w);pr={k:m.predict(x[ii]) for k,ii in [('validation',va),('test',te)] if len(ii)};best=float(np.average((pr['validation']-y[va])**2,weights=wv));parameters=sum(t.tree_.node_count for t in m.estimators_) if method=='chemistry_tree' else len(cols)+1
  with (out/'model.pkl').open('wb') as f:pickle.dump(dict(model=m,columns=cols,mean=xm,sd=xs,preprocessing=sup),f)
  checkpoint='model.pkl';selectedupdates=0;delta=None
 else:
  m=network(method,sup).to(device);parameters=sum(p.numel() for p in m.parameters());opt=torch.optim.AdamW(m.parameters(),lr=cfg['lr'],weight_decay=cfg['weight_decay']);initial={n:p.detach().clone() for n,p in m.named_parameters()};gen=torch.Generator().manual_seed(seed+10000)
  dd=tp(data,tr,device);yt=torch.tensor((y[tr]-mu)/sd,device=device);ww=torch.tensor(w,device=device);starep=0
  pairtr=[];pairva=[];lookup={x['record_id']:i for i,x in enumerate(rows)};local={int(i):j for j,i in enumerate(tr)};vset=set(va);tset=set(tr)
  for pair in pairs or []:
   a,b=lookup[pair['reference_id']],lookup[pair['changed_id']]
   assert (a in tset)==(b in tset) and (a in vset)==(b in vset)
   if a in tset:pairtr.append((local[a],local[b],pair['observed_difference'],rows[a]['sequence_group']))
   if a in vset:pairva.append((a,b,pair['observed_difference'],rows[a]['sequence_group']))
  if pairtr:
   pa=torch.tensor([p[0] for p in pairtr],device=device);pb=torch.tensor([p[1] for p in pairtr],device=device);py=torch.tensor([p[2]/sd for p in pairtr],dtype=torch.float32,device=device);pw=torch.tensor(weights([p[3] for p in pairtr]),dtype=torch.float32,device=device)
  if (out/'last.pt').exists():
   st=torch.load(out/'last.pt',map_location=device,weights_only=False);assert st['signature']==sig;m.load_state_dict(st['model']);opt.load_state_dict(st['optimizer']);history=st['history'];updates=st['updates'];best=st['best'];bestep=st['best_epoch'];starep=st['epoch'];gen.set_state(st['generator'].cpu());torch.set_rng_state(st['torch_rng'].cpu())
   if device.startswith('cuda'):torch.cuda.set_rng_state_all([a.cpu() for a in st['cuda_rng']])
  if device.startswith('cuda'):torch.cuda.synchronize();torch.cuda.reset_peak_memory_stats()
  fitregion=time.perf_counter()
  try:
   for ep in range(starep,500):
    m.train();perm=torch.randperm(len(tr),generator=gen);batch=len(tr) if pairs else 128;running=0.;gn=[]
    for j in range(0,len(tr),batch):
     ix=perm[j:j+batch].to(device);opt.zero_grad(set_to_none=True)
     if pairs:
      pp=m(**dd);loss=(ww*(pp-yt)**2).mean()
      if lam:loss=loss+lam*(pw*(pp[pb]-pp[pa]-py)**2).mean()
     else:
      pp=m(**{k:v[ix] for k,v in dd.items()});loss=(ww[ix]*(pp-yt[ix])**2).mean()
     if not torch.isfinite(loss):raise FloatingPointError('nonfinite training loss')
     loss.backward();gnorm=torch.nn.utils.clip_grad_norm_(m.parameters(),5.,error_if_nonfinite=True);gn.append(float(gnorm));maxgrad=max(maxgrad,float(gnorm));opt.step();updates+=1;running+=float(loss.detach())*len(ix)
    vp=predict(m,data,va,device,mu,sd);activity_score=float(np.average((vp-y[va])**2,weights=wv));score=activity_score
    if pairva and selection=='pair':
     v=dict(zip(va,vp));ee=np.array([(v[b]-v[a]-dy)**2 for a,b,dy,g in pairva]);score=float(np.average(ee,weights=weights([p[3] for p in pairva])))
    if not np.isfinite(score):raise FloatingPointError('nonfinite validation score')
    history.append(dict(epoch=ep+1,updates=updates,training_loss=running/len(tr),validation_score=score,validation_activity_mse=activity_score,max_batch_gradient_norm=max(gn)))
    if score<best-1e-8:
     best=score;bestep=ep+1;savept(out/'best.pt',dict(model=m.state_dict(),optimizer=opt.state_dict(),epoch=ep+1,updates=updates,target_mean=mu,target_sd=sd,signature=sig,preprocessing=sup))
    savept(out/'last.pt',dict(model=m.state_dict(),optimizer=opt.state_dict(),epoch=ep+1,updates=updates,history=history,best=best,best_epoch=bestep,signature=sig,generator=gen.get_state(),torch_rng=torch.get_rng_state(),cuda_rng=torch.cuda.get_rng_state_all() if device.startswith('cuda') else []))
    if ep+1>=50 and ep+1-bestep>=50:break
  except (FloatingPointError,RuntimeError) as e:
   status='failed';fail=repr(e);write(out/'failure.json',dict(reason=fail,updates=updates,history=history,policy='No restart. Training-fold mean prediction fallback; finite unfavorable seeds retained.'));best=None
  if device.startswith('cuda'):torch.cuda.synchronize()
  region=time.perf_counter()-fitregion
  if status=='completed':
   st=torch.load(out/'best.pt',map_location=device,weights_only=False);m.load_state_dict(st['model']);selectedupdates=st['updates'];delta=sum(float((p-initial[n]).detach().abs().sum()) for n,p in m.named_parameters());assert delta>0 and maxgrad>0
   pr={k:predict(m,data,ii,device,mu,sd) for k,ii in [('validation',va),('test',te)] if len(ii)};checkpoint='best.pt'
   # Data-only gradient audit at the selected checkpoint, dropout disabled; weight decay is not part of this loss.
   m.eval();opt.zero_grad(set_to_none=True)
   for j in range(0,len(tr),128):
    ii=tr[j:j+128];ad=tp(data,ii,device);loss=(torch.tensor(w[j:j+128],device=device)*(m(**ad)-torch.tensor((y[ii]-mu)/sd,device=device))**2).sum()/len(tr);loss.backward()
   audit=[]
   for n,p in m.named_parameters():
    if 'message_weights' in n:
     for r in range(8):audit.append(dict(parameter=n,relation=r,scalars=p[r].numel(),data_gradient_norm=float(p.grad[r].norm()) if p.grad is not None else 0.,parameter_l1_change=float((p[r]-initial[n][r]).detach().abs().sum()),training_directed_count=sup['directed_relation_counts'][r]))
   write(out/'gradient_audit.json',audit)
  else:pr={k:np.full(len(ii),mu) for k,ii in [('validation',va),('test',te)] if len(ii)};checkpoint='last.pt' if (out/'last.pt').exists() else None;selectedupdates=0;delta=None
  del m,opt,dd
 np.savez_compressed(out/'predictions.npz',**pr,test_indices=te,validation_indices=va)
 rec=dict(name=name,method=method,config=cfg,seed=seed,lambda_pair=lam,selection=selection,signature=sig,training_n=len(tr),validation_n=len(va),test_n=len(te),parameters=parameters,status=status,failure=fail,selected_epoch=bestep,selected_updates=selectedupdates,executed_updates=updates,epochs_executed=len(history),validation_score=best,target_mean=mu,target_sd=sd,training_row_mean=float(y[tr].mean()),training_equal_group_mean=mu,history=history,max_gradient=maxgrad,parameter_l1_change=delta,checkpoint=checkpoint,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,gpu_region_wall_s=region if method not in ['chemistry_tree','guide_ridge','pairwise_ridge'] and device.startswith('cuda') else 0.,hashes={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='fit.json'})
 write(out/'fit.json',rec);print(name,status,'epochs',len(history),'selected',bestep,'updates',updates,'val',best,flush=True);return rec,pr
