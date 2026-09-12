"""Ordinary regularized representer correction and source-only heuristic selection."""
import numpy as np
from .geometry import DIM,INITIAL,whitening,score_bank
from .profiles import inverse

def orthonormal_span(V,M,tol=1e-9):
    d,U=np.linalg.eigh((V.T@M@V+V.T@M.T@V)/2)
    keep=d>tol*max(1.,float(d.max()))
    return V@U[:,keep]/np.sqrt(d[keep])[None,:]

def regularized_solve(geo,V,lam):
    Q=orthonormal_span(V,geo['M']);G=Q.T@geo['G']@Q;b=Q.T@geo['b']
    A=(G+G.T)/2+lam*np.eye(len(b))
    small=np.linalg.solve(A,b);c=Q@small
    residual=geo['b']-geo['G']@c
    eigen,U=np.linalg.eigh((G+G.T)/2)
    null=eigen<1e-9*max(1.,float(eigen.max()))
    return dict(c=c,dimension=Q.shape[1],variance_model=float(c@geo['G']@c),
                residual=residual,solve_relative_residual=float(np.linalg.norm(A@small-b)/(1+np.linalg.norm(b))),
                regularized_condition=float(np.linalg.cond(A)),score_eigenvalues=eigen,
                score_null_target_norm=float(np.linalg.norm((U.T@b)[null])))

def score_observations(model,samples,delta,center):
    return [score_bank(model,a,inverse(model,a,y,delta),delta)-center[a] for a,y in enumerate(samples)]

def mean_scores(scores,weights=None):
    counts=np.array([len(s) for s in scores]);pi=counts/counts.sum() if weights is None else np.asarray(weights)
    return sum(p*s.mean(0) for p,s in zip(pi,scores))

def evaluation_correction(scores,c,weights=None):
    counts=np.array([len(s) for s in scores]);pi=counts/counts.sum() if weights is None else np.asarray(weights)
    action=[s@c for s in scores]
    correction=sum(p*x.mean() for p,x in zip(pi,action))
    variance=sum(p*p*np.var(x,ddof=1)/n for p,x,n in zip(pi,action,counts))
    return float(correction),float(variance)

def selection_score(solution,geo,selection_scores,evaluation_counts):
    """Heuristic local bias/variance proxy, never advertised as a risk certificate."""
    W=whitening(geo['M']);G=W.T@geo['G']@W
    mean=W.T@mean_scores(selection_scores)
    displacement=np.linalg.solve(G+.01*np.eye(DIM),mean)
    residual=W.T@solution['residual']
    proxy_bias=float(residual@displacement)
    proxy=solution['variance_model']/sum(evaluation_counts)+proxy_bias**2+.02**2*float(residual@residual)
    return float(proxy),proxy_bias

def fixed_space():return np.eye(DIM)[:,INITIAL]
