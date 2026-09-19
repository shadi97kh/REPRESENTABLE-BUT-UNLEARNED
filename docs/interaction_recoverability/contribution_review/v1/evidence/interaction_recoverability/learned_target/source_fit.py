"""Source-only CRPS fitting; no imports or interfaces for evaluator truth."""
import copy
from functools import lru_cache
import time
import numpy as np
import torch
from .profiles import ProfileModel,tensor

@lru_cache(maxsize=8)
def source_quadrature(nodes):
    z,w=np.polynomial.hermite.hermgauss(nodes)
    z=tensor(np.sqrt(2)*z);w=tensor(w/np.sqrt(np.pi))
    return z,w,w*(2*torch.cumsum(w,0)-w-1)

def crps(model,samples,delta,nodes=64):
    z,w,pair_weight=source_quadrature(nodes)
    total=sum(len(y) for y in samples);loss=0.
    for a,y in enumerate(samples):
        y=tensor(y);f=model.response(a,z,delta)
        loss=loss+len(y)/total*((y[:,None]-f[None,:]).abs()@w).mean()-len(y)/total*(pair_weight@f)
    return loss

def fit(samples,selection,delta,kind,seed,settings):
    begin=time.process_time();wall=time.monotonic()
    model=ProfileModel(kind,settings['width'],seed)
    optimizer=torch.optim.Adam(model.parameters(),lr=settings['learning_rate'])
    checkpoints=[];best=None
    for step in range(1,max(settings['checkpoints'])+1):
        optimizer.zero_grad();loss=crps(model,samples,delta,settings['source_nodes'])
        if not torch.isfinite(loss):raise ArithmeticError('nonfinite source objective')
        loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),50.)
        optimizer.step();model.project()
        if step in settings['checkpoints']:
            with torch.no_grad():validation=float(crps(model,selection if selection is not None else samples,delta,settings['source_nodes']))
            checkpoints.append({'step':step,'training_crps':float(loss),'selection_crps':validation})
            if selection is None or best is None or validation<best[0]:best=(validation,copy.deepcopy(model.state_dict()),step)
    model.load_state_dict(best[1])
    with torch.no_grad():
        q64=float(crps(model,selection if selection is not None else samples,delta,settings['source_nodes']))
        q128=float(crps(model,selection if selection is not None else samples,delta,2*settings['source_nodes']))
    return model,dict(cpu_s=time.process_time()-begin,wall_s=time.monotonic()-wall,checkpoints=checkpoints,
                      chosen_step=best[2],source_quadrature_difference=abs(q64-q128),
                      source_objective='empirical CRPS; deterministic normal quadrature approximation',
                      success=True)
