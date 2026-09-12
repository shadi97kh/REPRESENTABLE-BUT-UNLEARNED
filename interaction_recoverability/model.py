"""Validated two-interaction fixtures; mathematical domain x in [0,1]^4."""
import numpy as np

DELTAS=(.1,.01,.001,.0001)
COEFFICIENTS={
    'T':np.array([[.5,.5,1.,0.],[1.,1.,0.,1.]]),
    'Q':np.array([[.5,1.,1.,0.],[1.,.5,0.,1.]])}

def validate_delta(delta, *, boundary=False):
    if not np.isscalar(delta) or not np.isfinite(delta) or not (0<=delta<=.1) or (delta==0 and not boundary):
        raise ValueError('delta must be in (0,.1]; zero needs boundary=True and is not identifiable')
    return float(delta)

def validate_action(x):
    x=np.asarray(x,dtype=float)
    if x.shape!=(4,) or not np.isfinite(x).all() or (x<0).any() or (x>1).any():
        raise ValueError('action must be a finite four-vector in [0,1]^4')
    return x

def s(z):
    z=np.asarray(z,dtype=float)
    if not np.isfinite(z).all():raise ValueError('finite latent coordinates required')
    with np.errstate(over='raise',invalid='raise'):return z*np.hypot(1.,z)

def s_prime(z):
    z=np.asarray(z,dtype=float)
    if not np.isfinite(z).all():raise ValueError('finite latent coordinates required')
    return (1+2*z*z)/np.hypot(1.,z)

def h(delta,z,*,boundary=False):
    validate_delta(delta,boundary=boundary)
    return np.asarray(z,dtype=float)+delta*s(z)

def coefficients(model,x):
    if model not in COEFFICIENTS:raise ValueError('model must be T or Q')
    return np.exp(COEFFICIENTS[model]@validate_action(x))

def log_response(model,x,z,delta,*,boundary=False):
    a,b=coefficients(model,x)
    return np.logaddexp(np.log(a)+np.asarray(z),np.log(b)+h(delta,z,boundary=boundary))

def response(model,x,z,delta,*,boundary=False):
    with np.errstate(over='raise',invalid='raise'):
        return np.exp(log_response(model,x,z,delta,boundary=boundary))

def log_derivative(model,x,z,delta,*,boundary=False):
    a,b=coefficients(model,x)
    return np.logaddexp(np.log(a)+z,np.log(b)+h(delta,z,boundary=boundary)+np.log1p(delta*s_prime(z)))

def log_density_at_latent(model,x,z,delta,*,boundary=False):
    """log p(y) evaluated at y=G_model(z), including its Jacobian."""
    return -.5*z*z-.5*np.log(2*np.pi)-log_derivative(model,x,z,delta,boundary=boundary)
