"""Independent analytic graph response tasks, not trained GNN stand-ins."""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations
import numpy as np
import torch
from .transforms import BinaryGraphFamily


def all_actions(n=3):
    return [tuple(c) for k in range(n+1) for c in combinations(range(n),k)]


def family(task, seed=0):
    # Visible node-type tags distinguish locations and move with permutations.
    rng = np.random.default_rng(seed)
    x = torch.zeros((4,8),dtype=torch.float64)
    x[:,1:5] = torch.eye(4,dtype=torch.float64)
    x[:,5:] = torch.as_tensor(rng.normal(0,.1,(4,3)))
    adj = torch.zeros((4,4),dtype=torch.float64)
    if task == 'edge_switch':
        adj[0,1]=adj[1,0]=1
        switches=((0,2),(1,2),(2,3))
    elif task == 'node_state':
        for i in range(3): adj[i,i+1]=adj[i+1,i]=1
        switches=(0,1,2)
    else: raise ValueError('unknown task')
    return BinaryGraphFamily(task,x,adj,switches)


@dataclass(frozen=True)
class AnalyticPredictor:
    task: str
    coefficients: tuple
    intercept: float = 0.

    def __call__(self, state):
        # Resolve tagged nodes so the mathematical predictor is invariant to IDs.
        indices = state.x[:,1:5].argmax(0).tolist()
        if self.task == 'edge_switch':
            z = [float(state.adjacency[indices[i],indices[j]]) for i,j in ((0,2),(1,2),(2,3))]
        elif self.task == 'node_state':
            z = [float(state.x[indices[i],0]) for i in range(3)]
        else: raise ValueError('unknown task')
        terms = [z[0],z[1],z[2],z[0]*z[1],z[0]*z[2],z[1]*z[2],z[0]*z[1]*z[2]]
        return self.intercept + sum(w*t for w,t in zip(self.coefficients,terms))


def control_predictors(task):
    return {
        'additive': AnalyticPredictor(task,(.4,-.2,.7,0,0,0,0),.3),
        'pairwise': AnalyticPredictor(task,(.4,-.2,.7,.9,-.4,.2,0),.3),
        'higher_order': AnalyticPredictor(task,(.4,-.2,.7,.9,-.4,.2,1.2),.3),
        'indistinguishable_plus': AnalyticPredictor(task,(0,0,0,0,0,0,1),0),
        'indistinguishable_minus': AnalyticPredictor(task,(0,0,0,0,0,0,-1),0),
    }
