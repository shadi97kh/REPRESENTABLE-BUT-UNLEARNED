"""Known-warp two-quantile inversion; a standard minimum-distance control.

No fitting occurs on import. Sample estimation is explicitly gated and is never
called by the preparation CLI or its deterministic correctness tests.
"""
import numpy as np
from scipy.special import ndtri
from .model import h,validate_delta

PROBABILITIES=np.array([.25,.75])
LOW,HIGH=np.exp(.5),np.exp(1.)

def quantile_design(delta):
    validate_delta(delta)
    z=ndtri(PROBABILITIES)
    return np.column_stack([np.exp(z),np.exp(h(delta,z))])

def recover_coefficients_from_quantiles(quantiles,delta):
    """Exact two-by-two algebra; columns refer to anchored, known warp functions."""
    q=np.asarray(quantiles,float)
    if q.shape!=(2,2) or not np.isfinite(q).all():raise ValueError('rows=e1/e2, columns=.25/.75 quantiles')
    return np.clip(np.linalg.solve(quantile_design(delta),q.T).T,LOW,HIGH)

def primary_interaction(coefficients,m_delta):
    c=np.asarray(coefficients,float)
    if c.shape!=(2,2) or not np.isfinite(c).all() or not np.isfinite(m_delta):raise ValueError('invalid summaries')
    return float(np.exp(.5)*(c[0,0]-1)*(c[1,0]-1)+m_delta*(c[0,1]-1)*(c[1,1]-1))

def quantile_error_bound(delta,n_per_action,alpha=.05):
    """Union of two DKW events; returns None if central-range condition fails."""
    validate_delta(delta)
    if type(n_per_action) is not int or n_per_action<=0 or not 0<alpha<1:raise ValueError('invalid n/alpha')
    eps=np.sqrt(np.log(4/alpha)/(2*n_per_action))
    if eps>.125:return dict(informative=False,epsilon=float(eps))
    r=ndtri(.875)
    derivative_bound=HIGH*(np.exp(r)+np.exp(h(delta,r))*(1+delta*(1+2*r*r)/np.hypot(1,r)))
    quantile_lipschitz=derivative_bound/(np.exp(-r*r/2)/np.sqrt(2*np.pi))
    inverse=np.linalg.norm(np.linalg.inv(quantile_design(delta)),ord=np.inf)
    return dict(informative=True,epsilon=float(eps),quantile_lipschitz=float(quantile_lipschitz),
                inverse_norm=float(inverse),coefficient_sup_error=float(min(HIGH-LOW,inverse*quantile_lipschitz*eps)),
                status='STANDARD RESULT APPLIED; conservative finite-sample bound')

def estimate_from_samples(samples_e1,samples_e2,delta,m_delta,*,fitting_authorized=False):
    if not fitting_authorized:raise PermissionError('Sample fitting is outside preparation authorization')
    arrays=[np.asarray(v,float) for v in (samples_e1,samples_e2)]
    if any(a.ndim!=1 or not len(a) or not np.isfinite(a).all() or (a<=0).any() for a in arrays):raise ValueError('positive independent response arrays required')
    q=np.stack([np.quantile(a,PROBABILITIES,method='inverted_cdf') for a in arrays])
    c=recover_coefficients_from_quantiles(q,delta)
    return dict(coefficients=c.tolist(),interaction=primary_interaction(c,m_delta))


def protected_contrast_from_means(mean_e1,mean_zero):
    if not np.isfinite([mean_e1,mean_zero]).all():raise ValueError('finite means required')
    return float((np.e-1)*(mean_e1-mean_zero))
