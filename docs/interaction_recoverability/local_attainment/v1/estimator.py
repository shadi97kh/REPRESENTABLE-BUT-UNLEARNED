"""Finite operators for the local-attainment construction.

The public sampled-data entry is intentionally disabled for this continuation.
mpmath is an arbitrary-precision diagnostic backend, NOT a certified interval
backend. The statistical and certified-arithmetic specifications are in the
accompanying documents. No module-level computations or sample generation.
"""
from dataclasses import dataclass
from typing import Callable
import mpmath as mp

SAMPLED_EXECUTION_ENABLED = False


@dataclass(frozen=True)
class Schedule:
    N: int
    depth: int
    lattice: int
    radius: object
    mesh: object
    primitive_error: object
    precision_bits: int


def schedule(N: int) -> Schedule:
    if N < 10 or N % 5:
        raise ValueError("N=5n, with n>=2")
    lam = mp.mpf(575) / 613
    radius = mp.sqrt(mp.mpf(3) * mp.log(N) / 2)
    return Schedule(N, int(mp.ceil(2 * mp.log(N) / abs(mp.log(lam)))) + 2,
                    int(mp.ceil(4 * radius + 30)) + 2, radius,
                    mp.mpf(N) ** -6, mp.mpf(N) ** -40,
                    int(mp.ceil(10000 + 200 * mp.log(N, 2))))


def h(x):
    return x + x * mp.sqrt(1 + x*x) / 10


def hp(x):
    return 1 + (1 + 2*x*x) / (10 * mp.sqrt(1 + x*x))


def hinv(y, tolerance):
    if y == 0:
        return mp.mpf(0)
    sign = 1 if y > 0 else -1
    lo, hi = mp.mpf(0), abs(y)
    while hi - lo > tolerance:
        mid = (lo + hi) / 2
        if h(mid) < abs(y):
            lo = mid
        else:
            hi = mid
    return sign * (lo + hi) / 2


def phi(z):
    return mp.exp(-z*z/2) / mp.sqrt(2*mp.pi)


def Phi(z):
    return mp.erfc(-z/mp.sqrt(2)) / 2


def tail_lower(x):
    """Known Mills lower bound for Phi(-x), x>0."""
    return phi(x) * x / (1+x*x)


def domain_condition(s: Schedule):
    lower = min(tail_lower(mp.mpf(7)), tail_lower(s.radius + 1))
    effective_rank = (s.N // 5) * lower
    required = 64 * mp.log(s.N)**2
    return effective_rank >= required, effective_rank, required


def linear_quantile_from_sorted(values):
    """Polygon through (i/n, Y_(i)); constant on [0,1/n] and at 1.

    Deterministic order-statistic operator; no sampling, fitting or smoothing
    parameter is selected. Used by the disabled sample entry and arithmetic
    fixtures only.
    """
    y = tuple(mp.mpf(v) for v in values)
    if not y or any(y[i] > y[i+1] for i in range(len(y)-1)):
        raise ValueError("Nonempty sorted values required")
    n = len(y)
    def quantile(p):
        p = mp.mpf(p)
        if p <= mp.mpf(1)/n:
            return y[0]
        if p >= 1:
            return y[-1]
        rank = n*p
        k = int(mp.floor(rank))
        fraction = rank-k
        return (1-fraction)*y[k-1] + fraction*y[k]
    return quantile


class AnchorOperators:
    """Operators on response curves z -> Q(Phi(z)), not tangent q's."""
    def __init__(self, curves, depth, tolerance):
        if len(curves) != 5:
            raise ValueError("Exactly five source curves are required")
        self.curves = tuple(curves)
        self.depth = int(depth)
        self.tolerance = tolerance

    def transform(self, x, shift):
        return hinv(h(x)+shift, self.tolerance)

    def contract(self, x):
        return self.transform(x, 2)-1

    def E(self, x):
        f = self.curves
        t1, t2 = self.transform(x, 1), self.transform(x, 2)
        c = t2-1
        # Keep source-specific cancellations paired before summation.
        return (f[0](t2)-f[4](t1)) + (f[0](t1)-f[4](x)) + (f[0](c)-f[3](c))

    def difference(self, x, y=4):
        x, y = mp.mpf(x), mp.mpf(y)
        total = mp.mpf(0)
        for _ in range(self.depth):
            total -= self.E(x)-self.E(y)
            x, y = self.contract(x), self.contract(y)
        return total

    def g1(self, x):
        x = mp.mpf(x)
        shift = int(mp.ceil(4-x))
        center = x+shift
        f = self.curves
        base = 1 + mp.fsum(f[3](j)-f[0](j) for j in range(4))
        value = base+self.difference(center, 4)
        if shift >= 0:
            value -= mp.fsum(f[3](x+j)-f[0](x+j) for j in range(shift))
        else:
            value += mp.fsum(f[3](x+shift+j)-f[0](x+shift+j) for j in range(-shift))
        return value

    def g2(self, y):
        x = hinv(y, self.tolerance)
        return self.curves[0](x)-self.g1(x)

    def residual(self, pair_index, alpha, beta):
        return tuple(self.curves[pair_index](z)-self.g1(z+alpha)-self.g2(h(z)+beta)
                     for z in (mp.mpf(-1), mp.mpf(1)))


def signed_kernel(x, a, b, component, tolerance):
    def density(t):
        if component == 1:
            return phi(t)
        z = hinv(t, tolerance)
        return phi(z)/hp(z)
    return density(x-a-b)-density(x-a)-density(x-b)+density(x)


def target_weights(z, coeff, lattice, tolerance):
    alpha, beta = coeff[:2], coeff[2:]
    def A(x):
        return hp(x)*signed_kernel(h(x), *beta, 2, tolerance)
    def L(x):
        return signed_kernel(x, *alpha, 1, tolerance)-A(x)
    if z >= 4:
        W = mp.fsum(L(z+k) for k in range(1, lattice+1))
    else:
        W = -mp.fsum(L(z-k) for k in range(lattice))
    return A(z)-W, W


def periodized(t, coeff, lattice, tolerance):
    return mp.fsum(signed_kernel(t+k, *coeff[:2], 1, tolerance)
                   - hp(t+k)*signed_kernel(h(t+k), *coeff[2:], 2, tolerance)
                   for k in range(-lattice, lattice+1))


def midpoint(function: Callable, left, right, mesh):
    if right <= left:
        return mp.mpf(0)
    count = max(1, int(mp.ceil((right-left)/mesh)))
    step = (right-left)/count
    return step*mp.fsum(function(left+(j+mp.mpf(".5"))*step) for j in range(count))


def target_from_curves(operators, coeff, s: Schedule):
    """Finite, off-model functional; primary schedule is deliberately costly."""
    f = operators.curves
    def outside(z):
        u0, u3 = target_weights(z, coeff, s.lattice, s.primitive_error)
        return u0*f[0](z)+u3*f[3](z)
    cuts = sorted(set([-s.radius, s.radius] + ([mp.mpf(4)] if s.radius > 4 else [])))
    value = mp.fsum(midpoint(outside, a, b, s.mesh) for a,b in zip(cuts,cuts[1:]))
    value += midpoint(lambda t: periodized(t, coeff, s.lattice, s.primitive_error)
                      * operators.difference(t, 4), mp.mpf(4), mp.mpf(5), s.mesh)
    return value


def _grid_pair(operators, pair_index, N):
    """Finite lexicographic minimum-residual rule; never called in this task."""
    centers = [(mp.mpf(2)/3, mp.mpf(3)/5), (mp.mpf(3)/4, mp.mpf(4)/5)]
    ca, cb = centers[pair_index-1]
    radius = mp.mpf(1)/50
    count = (N+24)//25
    best, answer = mp.inf, None
    for j in range(count+1):
        a = ca-radius+2*radius*j/count
        for k in range(count+1):
            b = cb-radius+2*radius*k/count
            residual = operators.residual(pair_index, a, b)
            cost = mp.sqrt(mp.fsum(x*x for x in residual))
            if cost < best:
                best, answer = cost, (a,b)
    return answer


def estimate_from_sources(groups):
    """Inputs: only five balanced lists of responses in 0,e1,e2,e3,e4 order.

    Certified production arithmetic remains a separate backend requirement.
    This guard is intentionally unconditional under the authorized scope.
    """
    if not SAMPLED_EXECUTION_ENABLED:
        raise PermissionError("Sampled-data execution is disabled in this continuation")
    if len(groups) != 5 or len({len(g) for g in groups}) != 1:
        raise ValueError("Five balanced source groups required")
    n = len(groups[0])
    if n < 2:
        return {"estimate": 0, "fallback": "small_count"}
    N = 5*n
    s = schedule(N)
    with mp.workprec(s.precision_bits):
        s = schedule(N)
        # domain_condition(s) is a proof-onset diagnostic, not a data gate.
        if any(not mp.isfinite(y) or y <= 0 for g in groups for y in g):
            return {"estimate": 0, "fallback": "invalid_response"}
        if any(y > N*N for g in groups for y in g):
            return {"estimate": 0, "fallback": "response_cap"}
        quantiles = [linear_quantile_from_sorted(sorted(g)) for g in groups]
        curves = [(lambda z, q=q: q(Phi(z))) for q in quantiles]
        op = AnchorOperators(curves, s.depth, s.primitive_error)
        p1, p2 = _grid_pair(op, 1, N), _grid_pair(op, 2, N)
        coefficients = (p1[0],p2[0],p1[1],p2[1])
        raw = target_from_curves(op, coefficients, s)
        return {"estimate": min(mp.mpf(1000), max(mp.mpf(-1000), raw)),
                "coefficients": coefficients, "fallback": None,
                "arithmetic": "mpmath diagnostic backend; not a certified interval result"}
