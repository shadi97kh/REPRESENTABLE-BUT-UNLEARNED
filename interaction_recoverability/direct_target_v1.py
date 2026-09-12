"""Deterministic finite direct-target weights. No sampled-data fitting interface.

Coefficients here are inputs to an identity, not oracle inputs to a learner.
The statistical construction estimating them is specified in the proof notes.
"""
import math
import numpy as np
from scipy.special import ndtr

SAMPLED_DATA_FITTING_ENABLED = False


def phi(x):
    x = np.asarray(x, dtype=float)
    return np.exp(-x*x/2)/math.sqrt(2*math.pi)


def h(z, delta):
    if not 0 < delta <= .1:
        raise ValueError('Original model requires 0 < delta <= .1')
    z = np.asarray(z, dtype=float)
    return z+delta*z*np.sqrt(1+z*z)


def hprime(z, delta):
    z = np.asarray(z, dtype=float)
    return 1+delta*(1+2*z*z)/np.sqrt(1+z*z)


def hinverse(x, delta):
    """Input-dependent certified real-arithmetic bracket; no fixed tail cutoff."""
    h(0., delta)
    x = np.asarray(x, dtype=float)
    if not np.isfinite(x).all():
        raise ValueError('Finite arguments required')
    lo, hi = -np.abs(x), np.abs(x)
    for _ in range(64):
        mid = (lo+hi)/2
        below = h(mid, delta) < x
        lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
    return (lo+hi)/2


def density(x, component, delta):
    if component == 1:
        return phi(x)
    if component == 2:
        z = hinverse(x, delta)
        return phi(z)/hprime(z, delta)
    raise ValueError('Component must be 1 or 2')


def cdf(x, component, delta):
    if component == 1:
        return ndtr(x)
    if component == 2:
        return ndtr(hinverse(x, delta))
    raise ValueError('Component must be 1 or 2')


def kernel(x, c1, c2, component, delta, primitive=False):
    if not .5 <= c1 <= 1 or not .5 <= c2 <= 1:
        raise ValueError('Coefficients must be in the original box')
    f = cdf if primitive else density
    x = np.asarray(x, dtype=float)
    return f(x-c1-c2,component,delta)-f(x-c1,component,delta)-f(x-c2,component,delta)+f(x,component,delta)


def weight(x, c1, c2, component, delta, truncation, primitive=False):
    if isinstance(truncation, bool) or int(truncation) != truncation or truncation < 1:
        raise ValueError('Positive integer truncation required')
    return sum(kernel(np.asarray(x)+k,c1,c2,component,delta,primitive)
               for k in range(1,int(truncation)+1))


def latent_weights(z, coefficients, delta, truncation):
    """Order (alpha1,alpha2,beta1,beta2); returns combined reference, e3, e4."""
    a1,a2,b1,b2 = coefficients
    w3 = weight(z,a1,a2,1,delta,truncation)
    w4 = hprime(z,delta)*weight(h(z,delta),b1,b2,2,delta,truncation)
    return -w3-w4, w3, w4


def taper(z, cutoff):
    if cutoff <= 0:
        raise ValueError('Positive cutoff required')
    t = np.clip(np.abs(np.asarray(z))-cutoff,0,1)
    return 1-3*t*t+2*t*t*t


def probability_weights_at_z(z, coefficients, delta, truncation, cutoff):
    z = np.asarray(z, dtype=float)
    active = np.abs(z) < cutoff+1
    safe = np.where(active,z,0.)
    raw = latent_weights(safe,coefficients,delta,truncation)
    factor = np.where(active,taper(safe,cutoff)/phi(safe),0.)
    return tuple(factor*w for w in raw)
