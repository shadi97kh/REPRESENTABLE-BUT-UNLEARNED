"""Strong basis control: nonlinear direct-target search within empirical CDF constraints."""
import time
import math
import numpy as np
import torch
from .profiles import ProfileModel,inverse,target,tensor
from .nonlinear_bias import ks_statistic

def differentiable_ks(model,samples,delta):
    values=[]
    for a,y in enumerate(samples):
        y=tensor(np.sort(y));n=len(y)
        p=torch.special.ndtr(inverse(model,a,y,delta,implicit=True))
        hi=torch.arange(1,n+1,dtype=p.dtype)/n
        lo=torch.arange(n,dtype=p.dtype)/n
        values.append(torch.maximum((hi-p).max(),(p-lo).max()))
    return torch.stack(values)

def direct_target_set_search(pooled_model,samples,delta,settings):
    start=time.process_time();eps=tensor([math.sqrt(math.log(10/.05)/(2*len(y))) for y in samples])
    endpoints=[]
    for sign in [1.,-1.]:
        model=ProfileModel.from_payload(pooled_model.payload())
        optimizer=torch.optim.Adam(model.parameters(),lr=.01)
        best=None;feasible_steps=0;trace=[]
        for step in range(settings['endpoint_steps']+1):
            optimizer.zero_grad();ks=differentiable_ks(model,samples,delta);I=target(model,delta,settings['geometry_nodes'])
            violation=torch.relu(ks-eps)
            if float(violation.max())<=1e-10:
                feasible_steps+=1
                if best is None or sign*float(I)<sign*best['target']:
                    best={'target':float(I),'model':model.payload(),'ks':ks.detach().tolist(),'step':step}
            if step%10==0:trace.append({'step':step,'target':float(I),'max_violation':float(violation.max())})
            if step==settings['endpoint_steps']:break
            loss=sign*I+5000*violation.sum();loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(),100.)
            optimizer.step();model.project()
        endpoints.append({'sign':sign,'best':best,'feasible_steps':feasible_steps,'trace':trace})
    success=all(x['best'] is not None for x in endpoints)
    return {'prediction':float(np.mean([x['best']['target'] for x in endpoints])) if success else None,
            'success':success,'endpoints':endpoints,'cpu_s':time.process_time()-start,
            'status':'approximate inner endpoint search in constrained basis sieve; NOT a confidence interval',
            'original_class_coverage':'DKW set covers truth, but sieve/optimizer endpoints do not bound its full target range'}
