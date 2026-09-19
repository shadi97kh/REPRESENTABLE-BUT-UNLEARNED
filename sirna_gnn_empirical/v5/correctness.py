from common import *
from finite_range import *
import csv, math
protocol=json.loads((RUN/'correctness_protocol.json').read_text())
assert sha(RUN/'correctness_protocol.json')==json.loads((RUN/'correctness_protocol_freeze.json').read_text())['sha256']
q=[1,-1,-1,1];rows=[]
def case(name,result,expected):
 ok=all(result.get(k)==v for k,v in expected.items());rows.append(dict(case=name,passed=ok,result=result,expected=expected));assert ok,(name,result,expected)
def obs(values):return [([int(j==i) for j in range(4)],v,v) for i,v in enumerate(values)]
case('complete_nonlinear',exact_range([-2]*4,[2]*4,obs([0,0,0,1]),q),dict(status='exact',lower='1',upper='1'))
case('missing_joint',exact_range([-2]*4,[2]*4,obs([0,0,0]),q),dict(status='exact',lower='-2',upper='2'))
noise=[([int(j==i) for j in range(4)],str(F(v)-F(1,10)),str(F(v)+F(1,10))) for i,v in enumerate([0,0,0,1])]
case('bounded_noise',exact_range([-2]*4,[2]*4,noise,q),dict(status='exact',lower='3/5',upper='7/5'))
zero=exact_range([-2]*4,[2]*4,[],q,classes=[0]*4)
case('constant_encoder_false_zero_if_membership_unjustified',zero,dict(status='exact',lower='0',upper='0'))
case('constant_encoder_truth_incompatible',exact_range([-2]*4,[2]*4,obs([0,0,0,1]),q,classes=[0]*4),dict(status='empty'))
residual=exact_range([-2]*4,[2]*4,obs([0,0,0]),q,classes=[0]*4,eta='1/2')
case('approximation_error_retained',residual,dict(status='exact',lower='-1',upper='1'))
case('shared_endpoint_cancels_before_bounds',exact_range([-2]*3,[2]*3,[],[0,-1,1]),dict(status='exact',lower='-4',upper='4'))
case('encoder_alternatives_retained',union_range([zero,residual]),dict(status='exact',lower='-1',upper='1'))
case('empty_measurements',exact_range([-2]*4,[2]*4,obs([0,0,0,3]),q),dict(status='empty'))
case('budget_failure_not_zero',exact_range([-2]*4,[2]*4,[],q,budget=1),dict(status='unresolved'))
case('union_unresolved_not_zero',union_range([zero,{'status':'unresolved'}]),dict(status='unresolved'))
for name,v in [('nonfinite_nan','nan'),('nonfinite_inf','inf')]:
 try:exact_range([v]*4,[2]*4,[],q);raise AssertionError('nonfinite accepted')
 except ValueError:case(name,{'status':'rejected'},dict(status='rejected'))
# Neither missing joint support nor low-dimensional covariance forces a flat contrast.
cs=[-3,-1,0,2,5]
def fc(c,a,b):return 1+2*a+3*b+c*2*a*b
observed=[[fc(c,a,b) for a,b in [(0,0),(1,0),(0,1)]] for c in cs]
contrasts=[fc(c,1,1)-fc(c,1,0)-fc(c,0,1)+fc(c,0,0) for c in cs]
case('missing_joint_nonlinear_family',{'same_observed':len(set(map(tuple,observed)))==1,'contrast_values':contrasts},{'same_observed':True,'contrast_values':[-6,-2,0,4,10]})
# A local feasible extremum is an inner witness, not an outer certificate.
def nonconvex(x):return x*x-3*x**4+x**6
case('local_extrema_direction',{'second_derivative_at_zero':2,'local_min':nonconvex(0),'feasible_smaller_value':nonconvex(1)},{'second_derivative_at_zero':2,'local_min':0,'feasible_smaller_value':-1})
import torch
torch.set_num_threads(2)
z=torch.tensor([.2,-.1],dtype=torch.float64,requires_grad=True);target=.7;loss=((z[1]-z[0])-target)**2;loss.backward()
expected=torch.tensor([2.,-2.],dtype=torch.float64)
case('pair_endpoint_gradients',{'max_error':float((z.grad-expected).abs().max())},{'max_error':0.0})
# Independently absent route: data gradient zero; weight decay still changes it.
w=torch.tensor(2.,requires_grad=True);loss=(w*0+3)**2;loss.backward();before=w.item()
with torch.no_grad():w-=.1*(w.grad+.2*w)
case('absent_route_regularization_counterexample',{'data_gradient':float(w.grad),'weight_changed':w.item()!=before},{'data_gradient':0.0,'weight_changed':True})
import importlib.util
spec=importlib.util.spec_from_file_location('verified_architecture',ROOT/'sirna_gnn_empirical/v4/architecture.py');arch=importlib.util.module_from_spec(spec);spec.loader.exec_module(arch)
torch.manual_seed(1103);g=arch.GraphModel().double().eval();n=arch.GraphModel(kind='nongraph').double().eval();n.load_state_dict(g.state_dict())
X=torch.randn(1,64,93,dtype=torch.float64);mask=torch.ones(1,64,dtype=torch.float64);C=torch.zeros(1,8,dtype=torch.float64);A=torch.zeros(1,8,64,64,dtype=torch.float64);A[0,0,0,1]=1;B=A.clone();B[0,0,0,1]=0;B[0,0,0,2]=1
gn=float((g(X,A,mask,C)-g(X,B,mask,C)).abs().detach());nn=float((n(X,A,mask,C)-n(X,B,mask,C)).abs().detach())
case('actual_no_message_neighbor_invariance',{'parameters':sum(p.numel() for p in g.parameters()),'no_message_delta':nn,'gnn_delta_nonzero':gn>1e-12},{'parameters':153345,'no_message_delta':0.0,'gnn_delta_nonzero':True})
g.zero_grad();g(X,A,mask,C).sum().backward()
case('actual_absent_independent_route_data_gradient',{'absent_gradient':float(g.message_weights[0].grad[7].abs().max()),'present_gradient_nonzero':bool(g.message_weights[0].grad[0].abs().max()>0)},{'absent_gradient':0.0,'present_gradient_nonzero':True})
write(RUN/'correctness/results.json',dict(cases=rows,passed=len(rows),performance_pilot=False,new_fits=0,gradient_checks_not_training_fits=True))
with (RUN/'correctness/results.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=['case','passed','result']);writer.writeheader()
 for r in rows:writer.writerow({k:(json.dumps(r[k]) if k=='result' else r[k]) for k in writer.fieldnames})
print('PASSED',len(rows),'exact-reference/counterexample checks; zero performance fits')
