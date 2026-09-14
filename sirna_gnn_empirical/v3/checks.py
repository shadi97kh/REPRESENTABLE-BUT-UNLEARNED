"""Test actual failure modes on preserved measured inputs before final fitting."""
from common import *
from engine import *
rows,raw=load();pr=json.loads((RUN/'protocol.json').read_text());checks=[]
for pop in pr['activity_populations']:
 for f in pop['inner']:
  sets=[{x[k] for x in rows if x.get('study_group') in gg} for k in ['study_group','sequence_group','source_family'] for gg in [f['train_groups'],f['validation_groups'],pop['test_groups']]]
  for j in range(0,len(sets),3):assert not (sets[j]&sets[j+1] or sets[j]&sets[j+2] or sets[j+1]&sets[j+2])
checks.append('Zero protected study/source/sequence crossing in every activity fold')
tr=np.array([i for i,x in enumerate(rows) if x.get('split')=='train']);dd,sup=transform(raw,rows,tr);sample=np.array([tr[0]]);device=setup();batch=tp(dd,sample,device)
for method in ['original_gnn','corrected_gnn','original_no_message','corrected_no_message']:
 torch.manual_seed(7);m=network(method,sup).to(device).eval();h=torch.randn((1,64,32),device=device,requires_grad=True);A=batch['A'];agg=A.sum(-1)[:,:,:,None]*h[:,None,:,:] if 'no_message' in method else torch.matmul(A,h[:,None,:,:]);receiver=int(torch.nonzero(A.sum((1,3))[0])[0]);g=torch.autograd.grad(agg[0,:,receiver].sum(),h)[0];g[0,receiver]=0
 assert bool(g.abs().sum()>0)==('no_message' not in method)
 a=m(**batch);changed={k:v.clone() for k,v in batch.items()};changed['X'][changed['mask']==0]=500
 assert torch.allclose(m(**changed),a,atol=1e-6),method
checks+=['No-message aggregation has zero off-receiver hidden-state derivative; graph aggregation depends on neighbors','Padded node changes do not affect graph/no-message predictions']
# Transform scope test: arbitrary held-out chemistry values cannot alter fitted support.
other=dict(raw);other['X']=raw['X'].copy();other['X'][-1]=37;_,sup2=transform(other,rows,tr);assert sup==sup2
# Corrected sharing gives each unseen PO forward message exactly the learned common-base transform.
m=network('corrected_gnn',sup);w=m.message_weights[0];assert not sup['relation_support'][0];effective=w[2]+w[0]*sup['relation_support'][0];assert torch.equal(effective,w[2])
from metrics import ranking,metrics
z=pd.DataFrame(dict(activity=[1.,0.,.4,.7,.2,.6],prediction=[0.]*6,sequence_group=list('abcdef')));assert abs(ranking(z)['top5_mean_percentile']-.5)<1e-12
assert abs(metrics([1.,2.,3.],[2.,2.,2.])['r_squared'])<1e-12
pairs=readlines(OLD/'eligible_pairs.jsonl');lu={x['record_id']:x for x in rows}
for p in pairs:assert abs(lu[p['changed_id']]['activity']-lu[p['reference_id']]['activity']-p['observed_difference'])<1e-9
for f in pr['pair_populations']:
 for inner in f['inner']:assert not(set(inner['train_groups'])&set(inner['validation_groups']) or set(inner['train_groups'])&set(f['test_groups']) or set(inner['validation_groups'])&set(f['test_groups']))
checks+=['Held-out features do not change fitted masks','Unseen PO routes to trained directional base in corrected arm','Tied constant ranking is exactly random-tie expectation .5','Oracle-mean R-squared is zero','Measured pair orientation and grouped pair folds verified']
write(RUN/'scientific_checks.json',dict(passed=checks,relations=sup['directed_relation_counts'],parameters={m:sum(p.numel() for p in network(m,sup).parameters()) for m in pr['methods']}));print(checks)
