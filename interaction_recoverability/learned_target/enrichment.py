"""Residual-driven Galerkin enrichment; finite search optimum, no whole-class bound."""
import time
import numpy as np
from .geometry import DIM
from .correction import fixed_space,regularized_solve,orthonormal_span,selection_score

def construct(geo,lam,mode,rounds,seed,tau=.001):
    start=time.process_time();M=geo['M'];G=geo['G'];b=geo['b']
    V=np.eye(DIM) if mode=='rich' else fixed_space()
    trajectory=[];rng=np.random.default_rng(seed);witness=None
    for step in range(rounds if mode in ['adaptive','random'] else 0):
        solution=regularized_solve(geo,V,lam)
        residual=b-G@solution['c']
        rhs=residual if mode=='adaptive' else rng.normal(size=DIM)
        # Same stabilized solve and number of appended directions for matched random control.
        v=np.linalg.solve(G+tau*M,rhs)
        v/=np.sqrt(v@M@v)
        before=float((residual@v)**2/(v@G@v+tau))
        Q=orthonormal_span(V,M);orth=v-Q@(Q.T@M@v)
        novelty=float(np.sqrt(max(0.,orth@M@orth)))
        if novelty<1e-7:
            trajectory.append(dict(round=step,accepted=False,nonredundancy=novelty,witness=before));break
        V=np.column_stack([V,v]);after_solution=regularized_solve(geo,V,lam)
        after=float((after_solution['residual']@v)**2/(v@G@v+tau))
        trajectory.append(dict(round=step,accepted=True,nonredundancy=novelty,witness=before,
                               same_direction_residual_after=after,physical_direction=v.tolist(),
                               guarantee='global Rayleigh optimum only in declared finite numerical-metric sieve' if mode=='adaptive' else 'matched preconditioned random direction'))
        witness=v
    solution=regularized_solve(geo,V,lam)
    solution.update(trajectory=trajectory,witness_direction=witness,cpu_s=time.process_time()-start,mode=mode,lambda_=lam)
    return solution

def select(geo,selection_scores,evaluation_counts,mode,settings,seed):
    selection_start=time.process_time()
    candidates=[]
    for lam in settings['lambdas']:
        candidate=construct(geo,lam,mode,settings['enrichment_rounds'],seed,settings['search_tau'])
        proxy,bias=selection_score(candidate,geo,selection_scores,evaluation_counts)
        candidate['selection_proxy']=proxy;candidate['selection_bias_proxy']=bias
        candidates.append(candidate)
    best=min(candidates,key=lambda x:x['selection_proxy'])
    best['all_lambda_proxies']=[{'lambda':c['lambda_'],'proxy':c['selection_proxy'],'dimension':c['dimension']} for c in candidates]
    best['construction_cpu_all_lambdas']=sum(c['cpu_s'] for c in candidates)
    best['selection_and_construction_cpu_s']=time.process_time()-selection_start
    return best
