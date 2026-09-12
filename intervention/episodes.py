"""Lazy function-level partitions; no model zoo or teacher table is generated."""
import numpy as np
from .tasks import AnalyticPredictor


def source_function(task,split,index,variant):
    boundaries={'train':(1000,64),'development':(2000,16),'test':(3000,32)}
    offset,count=boundaries[split]
    if index not in range(count) or variant not in (0,1): raise ValueError('undeclared source index')
    task_index={'edge_switch':0,'node_state':1}[task]
    rng=np.random.default_rng(np.random.SeedSequence([offset,index,task_index]))
    first=rng.normal(size=3)
    higher=rng.normal(size=(2,4))[variant]
    return AnalyticPredictor(task,tuple(np.r_[first,higher]),0.)
