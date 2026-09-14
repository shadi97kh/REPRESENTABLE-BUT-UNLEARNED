"""Selected-checkpoint data-active scalar counts; read-only, no optimizer update."""
from common import *
from engine import load,transform,setup,network,tp
import torch
rows,raw=load();device=setup();y=np.array([r['activity'] for r in rows],np.float32);allrows=[];summary=[];t0=time.perf_counter()
for f in sorted((RUN/'fits').rglob('fit.json')):
 r=json.loads(f.read_text())
 if '/final/' not in r['name'] or r['executed_updates']<=0 or r['status']!='completed':continue
 start=time.perf_counter();tr=np.array(r['signature']['train']);data,sup=transform(raw,rows,tr);assert sup==json.loads((f.parent/'preprocessing.json').read_text());group='sequence_group' if r['name'].startswith('pair/') else 'study_group';w=weights([rows[i][group] for i in tr]).astype(np.float32);st=torch.load(f.parent/'best.pt',map_location=device,weights_only=False);m=network(r['method'],sup).to(device).eval();m.load_state_dict(st['model']);mu,sd=st['target_mean'],st['target_sd']
 m.zero_grad(set_to_none=True)
 for j in range(0,len(tr),128):
  ii=tr[j:j+128];loss=(torch.tensor(w[j:j+128],device=device)*(m(**tp(data,ii,device))-torch.tensor((y[ii]-mu)/sd,device=device))**2).sum()/len(tr);loss.backward()
 # CPU constructor follows the same initialization sequence as the original pre-device model.
 torch.manual_seed(r['seed']);initial=network(r['method'],sup);init={k:v.detach() for k,v in initial.named_parameters()};total=active=changed=0;gn2=0
 for name,p in m.named_parameters():
  size=p.numel();nz=int(torch.count_nonzero(p.grad)) if p.grad is not None else 0;norm=float(p.grad.norm()) if p.grad is not None else 0.;chg=int(torch.count_nonzero(p.detach().cpu()!=init[name]));total+=size;active+=nz;changed+=chg;gn2+=norm**2;allrows.append(dict(fit=r['name'],method=r['method'],seed=r['seed'],parameter=name,scalars=size,selected_activity_gradient_nonzero_scalars=nz,selected_activity_gradient_zero_scalars=size-nz,selected_activity_gradient_norm=norm,changed_from_initialized_scalars=chg))
 objective_active=active;objective_gn2=gn2
 if r.get('lambda_pair',0)>0:
  ppairs=readlines(OLD/'eligible_pairs.jsonl');local={int(v):i for i,v in enumerate(tr)};lookup={x['record_id']:i for i,x in enumerate(rows)};selected=[p for p in ppairs if lookup[p['reference_id']] in local];aa=torch.tensor([local[lookup[p['reference_id']]] for p in selected],device=device);bb=torch.tensor([local[lookup[p['changed_id']]] for p in selected],device=device);pw=torch.tensor(weights([rows[lookup[p['reference_id']]]['sequence_group'] for p in selected]),dtype=torch.float32,device=device);dy=torch.tensor([p['observed_difference']/sd for p in selected],dtype=torch.float32,device=device);m.zero_grad(set_to_none=True);z=m(**tp(data,tr,device));loss=(torch.tensor(w,device=device)*(z-torch.tensor((y[tr]-mu)/sd,device=device))**2).mean()+r['lambda_pair']*(pw*(z[bb]-z[aa]-dy)**2).mean();loss.backward();objective_active=sum(int(torch.count_nonzero(p.grad)) for p in m.parameters() if p.grad is not None);objective_gn2=sum(float(p.grad.norm())**2 for p in m.parameters() if p.grad is not None)
 parameter_map=dict(m.named_parameters())
 for record in allrows[-len(parameter_map):]:
  grad=parameter_map[record['parameter']].grad;record['selected_training_objective_nonzero_scalars']=int(torch.count_nonzero(grad)) if grad is not None else 0;record['selected_training_objective_gradient_norm']=float(grad.norm()) if grad is not None else 0.
 assert total==r['parameters'];assert sha(f.parent/'best.pt')==r['hashes']['best.pt'];summary.append(dict(fit=r['name'],method=r['method'],seed=r['seed'],parameters=total,selected_training_objective_nonzero_scalars=objective_active,selected_training_objective_gradient_norm=objective_gn2**.5,lambda_pair=r.get('lambda_pair',0),selected_activity_gradient_nonzero_scalars=active,selected_activity_gradient_zero_scalars=total-active,changed_from_initialized_scalars=changed,selected_activity_gradient_norm=gn2**.5,diagnostic_wall_s=time.perf_counter()-start,optimizer_updates=0,gradient_scope='Full training-partition normalized activity loss; dropout disabled, including pair-trained checkpoints. Exact zero comparison in recorded float arithmetic; not a count of all scalars ever data-updated over training.'))
 del m,initial
export(pd.DataFrame(allrows),'tables/selected_parameter_activity.csv');export(pd.DataFrame(summary),'tables/selected_parameter_activity_summary.csv');write(RUN/'parameter_audit_summary.json',dict(final_neural_checkpoints=len(summary),parameter_tensor_rows=len(allrows),new_fits=0,optimizer_updates=0,wall_s=time.perf_counter()-t0,scope='Read-only full-training activity and actual combined-objective gradients at all completed final neural checkpoints; dropout disabled. Parameter changes include weight decay. Ever-data-active counts throughout the full trajectory are not inferred.'))
print('Audited',len(summary),'final neural checkpoints;',len(allrows),'parameter tensors; zero fits/updates.')
