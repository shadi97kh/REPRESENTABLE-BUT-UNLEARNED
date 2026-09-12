"""Finite-difference profile reconstruction and conservative non-oracle estimator.

The sample interface is disabled by default. No sample calls occur in preparation.
The mathematical bound is exact-arithmetic; floating point is not certified.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.special import ndtr, log_ndtr
from ..model import h, s_prime, validate_delta, validate_action


def _positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name + ' must be a positive integer')


def inverse_h(delta, x):
    """Known warp inverse, bracketed by x and zero because h' >= 1."""
    validate_delta(delta)
    if not math.isfinite(x):
        raise ValueError('finite profile coordinate required')
    if x == 0:
        return 0.0
    return brentq(lambda z: float(h(delta, z)) - x, min(x, 0), max(x, 0), xtol=1e-13)


def separation(delta):
    validate_delta(delta)
    a0, b0 = .5 * math.exp(.5), 2 * math.e
    c = math.log(2 * b0 / a0)
    ratio = c / delta
    # hypot avoids squaring overflow for tiny delta; fail explicitly if unrepresentable.
    if not math.isfinite(ratio):
        raise ValueError('delta below supported numerical range')
    radius = math.sqrt((math.hypot(1, 2 * ratio) - 1) / 2)
    q = math.exp(-c)
    return dict(R=radius, hR=float(h(delta, radius)), q=q,
                gamma=a0-b0*q, J=b0*(1+q))


def reconstruct_profile(quantile_curve, component, x, delta, K):
    """Quantile curve signature (action index 0..4, latent z); no true profile input.

    Caller identifies whether curves are empirical or exact diagnostic fixtures.
    The K-term finite difference sum leaves the mathematical remainder g(x-K).
    """
    _positive_int(K, 'K')
    if component not in (1, 2) or not math.isfinite(x):
        raise ValueError('component 1/2 and finite x required')
    validate_delta(delta)
    terms = []
    for k in range(1, K+1):
        z = x-k if component == 1 else inverse_h(delta, x-k)
        anchor = 3 if component == 1 else 4
        terms.append(float(quantile_curve(anchor, z))-float(quantile_curve(0, z)))
    answer = math.fsum(terms)
    if not math.isfinite(answer):
        raise ArithmeticError('nonfinite profile reconstruction')
    return answer


def tail_envelope(delta, B):
    validate_delta(delta)
    if not math.isfinite(B) or B < 0:
        raise ValueError('nonnegative finite B required')
    a = .5-delta
    first = math.exp(.5)*(ndtr(1-B)+ndtr(-B-1))
    second = 2*math.exp(delta/2+1/(4*a))/math.sqrt(2*a)*ndtr(-math.sqrt(2*a)*(B-1/(2*a)))
    return float(first+second)


def target_range(delta):
    # tail_envelope(0) includes m0 exactly and an upper bound on m_delta.
    return 2*(math.e+1)**2*tail_envelope(delta, 0)


def preflight(delta, counts, K, B, failure_probability=.05):
    """No fitting or simulated observations. Check the proved DKW query domain."""
    _positive_int(K, 'K')
    if len(counts) != 5 or any(type(n) is not int or n < 1 for n in counts):
        raise ValueError('five positive integer stratum counts required')
    if not math.isfinite(B) or B <= 0 or not 0 < failure_probability < 1:
        raise ValueError('positive finite cutoff and probability in (0,1) required')
    sep = separation(delta)
    T = max(B, sep['R'])+K+2
    lo, hi = float(log_ndtr(-T)), float(log_ndtr(-T-1))
    log_gap = lo+math.log1p(-math.exp(hi-lo))
    log_eps = .5*(math.log(math.log(10/failure_probability)/2)-math.log(min(counts)))
    result = dict(**sep, T=T, U=T+1, log_epsilon=log_eps, log_cdf_gap=log_gap,
                  sufficient_log10_n=(math.log(math.log(10/failure_probability)/2)-2*log_gap)/math.log(10),
                  informative=log_eps <= log_gap, target_abs_bound=target_range(delta),
                  status='PROVED HERE bound evaluated in uncertified floating point')
    if result['informative']:
        U = T+1
        log_L = math.log(2*math.e)+float(np.logaddexp(U, float(h(delta, U))+math.log1p(delta*float(s_prime(U)))))+.5*U*U+.5*math.log(2*math.pi)
        if log_L+log_eps > 700:
            raise ArithmeticError('bound outside supported floating-point range')
        result['quantile_error'] = math.exp(log_L+log_eps)
    return result


def error_bound(delta, counts, K, B, grid_steps, integration_bins, failure_probability=.05):
    """Explicit P7 bound, including stochastic profile error and deterministic bias."""
    _positive_int(grid_steps, 'grid_steps')
    _positive_int(integration_bins, 'integration_bins')
    out = preflight(delta, counts, K, B, failure_probability)
    if not out['informative']:
        out['mse_bound'] = 4*out['target_abs_bound']**2
        return out
    EQ, R = out['quantile_error'], out['R']
    P = 4*K*math.exp(R)*EQ+2*math.e*(1+out['q'])*math.exp(-K)
    tau = .5/grid_steps  # conservative covering distance; grid endpoints included
    rho = min(.5, (2*(P+math.exp(R)*EQ)+out['J']*tau)/out['gamma'])
    edges = np.linspace(-B, B, integration_bins+1)
    z = (edges[:-1]+edges[1:])/2
    w = ndtr(edges[1:])-ndtr(edges[:-1])
    mG = float(w @ (np.exp(z)+np.exp(h(delta, z))))
    CH = 2*(math.e+1)**2*(math.exp(B)+math.exp(float(h(delta, B)))*(1+delta*float(s_prime(B))))
    terms = dict(stochastic_profile=16*K*EQ,
                 profile_truncation=2*(math.e+1)**2*mG*math.exp(-K),
                 coefficient_and_grid=4*(math.e**2+math.e)*mG*rho,
                 integration=CH*B/integration_bins,
                 outcome_tail=2*(math.e+1)**2*tail_envelope(delta, B))
    total = math.fsum(terms.values())
    out.update(coefficient_error=rho, error_terms=terms, absolute_error_bound=total,
               mse_bound=min(4*out['target_abs_bound']**2, total**2+4*out['target_abs_bound']**2*failure_probability),
               arithmetic_error='not certified; must be added for numerical theorem guarantees')
    return out


def estimate_from_samples(samples, delta, *, K, B, grid_steps, integration_bins,
                          failure_probability=.05, fitting_authorized=False):
    """Actual empirical-quantile/grid algorithm; preparation never calls this enabled.

    Input: five independent positive response arrays ordered 0,e1,e2,e3,e4.
    All tuning/truncation inputs are explicit and must be fixed before a future fit.
    No true profiles, coefficients, latent labels or held-out responses are inputs.
    """
    if not fitting_authorized:
        raise PermissionError('Sample fitting is not authorized in this preparation task')
    arrays = [np.asarray(x, dtype=float) for x in samples]
    if len(arrays) != 5 or any(x.ndim != 1 or len(x) == 0 or not np.isfinite(x).all() or (x <= 0).any() for x in arrays):
        raise ValueError('five nonempty finite positive one-dimensional arrays required')
    bounds = error_bound(delta, [len(x) for x in arrays], K, B, grid_steps, integration_bins, failure_probability)
    if not bounds['informative']:
        return dict(interaction=0., status='full-range fallback; sufficient quantile condition fails', bound=bounds)
    arrays = [np.sort(x) for x in arrays]

    def curve(a, z):
        if abs(z) > bounds['T']+1e-10:
            raise ValueError('query outside proved range')
        p = float(ndtr(z))
        if not 0 < p < 1:
            raise ArithmeticError('quantile probability not representable; higher precision required')
        index = max(0, min(len(arrays[a])-1, math.ceil(len(arrays[a])*p)-1))
        return float(arrays[a][index])

    def profile(j, x):
        return reconstruct_profile(curve, j, x, delta, K)

    grid = np.linspace(.5, 1, grid_steps+1)
    coefficients = []
    for a in (1, 2):
        best = None
        for alpha in grid:
            for beta in grid:
                residuals = []
                for z, weight in [(-bounds['R'], math.exp(bounds['R'])), (bounds['R'], math.exp(-bounds['hR']))]:
                    residuals.append(weight*(profile(1,z+alpha)+profile(2,float(h(delta,z))+beta)-curve(a,z)))
                residual = max(abs(x) for x in residuals)
                if best is None or residual < best[0]:
                    best = (residual, float(alpha), float(beta))
        coefficients.append(best[1:])
    edges = np.linspace(-B, B, integration_bins+1)
    terms = []
    for left, right in zip(edges[:-1], edges[1:]):
        z = (left+right)/2
        contrast = 0.
        for j, x in [(1,z),(2,float(h(delta,z)))]:
            c1, c2 = coefficients[0][j-1], coefficients[1][j-1]
            contrast += profile(j,x+c1+c2)-profile(j,x+c1)-profile(j,x+c2)+profile(j,x)
        terms.append(float(ndtr(right)-ndtr(left))*contrast)
    target = float(np.clip(math.fsum(terms), -bounds['target_abs_bound'], bounds['target_abs_bound']))
    return dict(interaction=target, coefficients=coefficients, bound=bounds,
                status='ordinary inverse-problem baseline; floating-point error uncertified')


def transport_score(z, F_prime, F_second, dot_F, dot_F_prime):
    """Observed-law score at y=F(z), including inverse-map and Jacobian terms."""
    if not np.isfinite([z,F_prime,F_second,dot_F,dot_F_prime]).all() or F_prime <= 0:
        raise ValueError('finite quantities and positive response derivative required')
    return z*dot_F/F_prime-dot_F_prime/F_prime+dot_F*F_second/F_prime**2


def sinusoidal_profile(amplitude, frequency):
    """Exact analytic diagnostic fixture only; never passed to a sample estimator."""
    if abs(amplitude) > .05 or not 1 <= frequency <= 2:
        raise ValueError('fixture outside explicitly checked subclass')
    def profile(x, derivative=0):
        return float(math.exp(x)*(1+amplitude*((1+1j*frequency)**derivative*complex(math.cos(frequency*x),math.sin(frequency*x))).imag))
    return profile


def response_derivatives(a, z, delta, coefficients, profiles):
    """Known-parameter fixture for checking the mathematical score; not a learner."""
    a = validate_action(a)
    c = np.asarray(coefficients, float)
    if c.shape != (2,2) or not np.isfinite(c).all() or (c < .5).any() or (c > 1).any():
        raise ValueError('coefficient rows e1/e2, columns alpha/beta in [.5,1]')
    x = z+c[0,0]*a[0]+c[1,0]*a[1]+a[2]
    y = float(h(delta,z))+c[0,1]*a[0]+c[1,1]*a[1]+a[3]
    hp = 1+delta*float(s_prime(z))
    hpp = delta*z*(3+2*z*z)/(1+z*z)**1.5
    g1,g2 = profiles
    return (g1(x)+g2(y),g1(x,1)+g2(y,1)*hp,
            g1(x,2)+g2(y,2)*hp**2+g2(y,1)*hpp,x,y,hp)
