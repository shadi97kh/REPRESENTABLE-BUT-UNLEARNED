"""Deterministic constants for proved bounds, with explicit conventions."""
import math
import numpy as np

LIMIT=math.exp(.5)*(math.e+math.exp(2)-2*math.exp(1.5))

def score_information_bound():
    # P(t)=1.2 t^3+.3 t^2+2.6 t+1.15 bounds |delta-score| uniformly.
    square=np.polynomial.polynomial.polymul([1.15,2.6,.3,1.2],[1.15,2.6,.3,1.2])
    moments=[2**(k/2)*math.gamma((k+1)/2)/math.sqrt(math.pi) for k in range(7)]
    return float(np.dot(square,moments))

def experiment_bound(delta,counts):
    from .model import validate_delta
    validate_delta(delta)
    if len(counts)!=5 or any(type(n) is not int or n<0 for n in counts):
        raise ValueError('five nonnegative integer counts ordered [0,e1,e2,e3,e4] required')
    c=score_information_bound();n2=counts[2]
    h2_single=min(2.,c*delta**2)
    h2_product=2*(-math.expm1(n2*math.log1p(-h2_single/2))) if h2_single<2 else (2. if n2 else 0.)
    tv=min(1.,math.sqrt(max(0,h2_product)))
    return dict(hellinger_convention='H^2=int(sqrt(p)-sqrt(q))^2',C=c,n_informative=n2,
                h2_single_upper=h2_single,h2_experiment_upper=h2_product,tv_upper=tv,
                target_separation_lower=LIMIT,mse_lower=LIMIT**2*(1-tv)/8,
                status='STANDARD RESULT APPLIED',transport_bound_status='PROVED HERE')

def estimable_linear_functional(operator,target,tolerance=1e-10):
    """Known-operator row-space diagnostic, not an estimator fitted to samples."""
    a,l=np.asarray(operator,float),np.asarray(target,float)
    if a.ndim!=2 or l.shape!=(a.shape[1],):raise ValueError('operator/target dimensions')
    w=np.linalg.pinv(a.T)@l
    residual=float(np.linalg.norm(a.T@w-l))
    return dict(estimable=residual<=tolerance,weights=w.tolist(),residual=residual,
                iid_unit_noise_variance_factor=float(w@w) if residual<=tolerance else None)


def normalized_class_lower_bound(delta,n_per_action):
    from .model import validate_delta
    validate_delta(delta)
    if type(n_per_action) is not int or n_per_action<0:raise ValueError('nonnegative integer count required')
    m=math.sqrt(math.e);width=math.e-m;kappa=m*width
    eta=width if n_per_action==0 else min(width,m/(delta*math.sqrt(score_information_bound()*n_per_action)))
    return dict(coefficient_separation=eta,target_separation_lower=kappa*eta,
                mse_lower=kappa*kappa*eta*eta/16,
                rate='min(1,1/(n*delta^2)); constants uniform for fixed known delta in (0,.1]',
                status='STANDARD RESULT APPLIED',transport_bound_status='PROVED HERE')
