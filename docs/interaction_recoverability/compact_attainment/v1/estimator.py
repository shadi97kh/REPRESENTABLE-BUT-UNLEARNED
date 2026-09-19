"""Compact finite operators; analytic fixtures only, no sampled execution.

mpmath is a diagnostic backend, not an interval certificate. The Gauss bound
below certifies truncation in exact arithmetic given its derivative bound;
rounding and knot enclosures have separate obligations in algorithm.md.
"""
from dataclasses import dataclass
import mpmath as mp

SAMPLED_EXECUTION_ENABLED = False


def h(x):
    return x + x * mp.sqrt(1+x*x)/10


def hp(x):
    return 1+(1+2*x*x)/(10*mp.sqrt(1+x*x))


def hinv(y, tol):
    if tol <= 0:
        raise ValueError("positive tolerance required")
    sign = 1 if y >= 0 else -1
    lo, hi = mp.mpf(0), abs(y)
    while hi-lo > tol:
        mid = (lo+hi)/2
        if h(mid) < abs(y):
            lo = mid
        else:
            hi = mid
    return sign*(lo+hi)/2


def phi(x):
    return mp.exp(-x*x/2)/mp.sqrt(2*mp.pi)


def Phi(x):
    return mp.erfc(-x/mp.sqrt(2))/2


@dataclass(frozen=True)
class Schedule:
    N: int
    depth: int
    lattice: int
    radius: object
    iterations: int


def schedule(N):
    if N < 10 or N % 5:
        raise ValueError("N=5n, n>=2")
    lam = mp.mpf(271)/290
    R = mp.sqrt(3*mp.log(N)/2)
    return Schedule(N, int(mp.ceil(2*mp.log(N)/abs(mp.log(lam))))+2,
                    int(mp.ceil(4*R+30))+2, R,
                    2*int(mp.ceil(mp.log(N)/abs(mp.log(mp.mpf(3)/4)))))


def polygon(values):
    """Deterministic polygon through supplied sorted heights at i/n."""
    y = tuple(map(mp.mpf, values))
    if not y or any(a > b for a, b in zip(y, y[1:])):
        raise ValueError("nonempty sorted heights required")
    n = len(y)
    def evaluate(p):
        if p <= mp.mpf(1)/n:
            return y[0]
        if p >= 1:
            return y[-1]
        t = n*p
        k = int(mp.floor(t))
        return y[k-1]+(t-k)*(y[k]-y[k-1])
    return evaluate


class CompactOperators:
    def __init__(self, curves, depth, tolerance):
        if len(curves) != 5 or depth < 0 or tolerance <= 0:
            raise ValueError("five curves, nonnegative depth, positive tolerance")
        self.curves, self.depth, self.tolerance = tuple(curves), depth, tolerance

    def T(self, x, j):
        return hinv(h(x)+j, self.tolerance)

    def C(self, x):
        return self.T(x, 5)-4

    def E(self, x):
        f = self.curves
        ts = [mp.mpf(x)]+[self.T(x, j) for j in range(1, 6)]
        return (mp.fsum(f[0](ts[j])-f[4](ts[j-1]) for j in range(1, 6))
                +mp.fsum(f[0](ts[5]-k)-f[3](ts[5]-k) for k in range(1, 5)))

    def difference(self, x, y=-1):
        x, y, total = mp.mpf(x), mp.mpf(y), mp.mpf(0)
        for _ in range(self.depth):
            total -= self.E(x)-self.E(y)
            x, y = self.C(x), self.C(y)
        return total

    def g1(self, x):
        x = mp.mpf(x)
        shift = int(mp.ceil(-mp.mpf(7)/4-x))
        f = self.curves
        value = 1-(f[3](-1)-f[0](-1))+self.difference(x+shift)
        if shift > 0:
            value -= mp.fsum(f[3](x+j)-f[0](x+j) for j in range(shift))
        elif shift < 0:
            value += mp.fsum(f[3](x-k)-f[0](x-k) for k in range(1, 1-shift))
        return value

    def g2(self, y):
        x = hinv(y, self.tolerance)
        return self.curves[0](x)-self.g1(x)

    def B(self, c):
        a, b = c
        return mp.matrix([self.g1(z+a)+self.g2(h(z)+b) for z in (-1, 1)])


def public_matrix(pair):
    centers = ((mp.mpf(2)/3, mp.mpf(3)/5), (mp.mpf(3)/4, mp.mpf(4)/5))
    if pair not in (1, 2):
        raise ValueError("pair 1 or 2")
    a, b = centers[pair-1]
    return mp.matrix([[mp.exp(-1+a), mp.exp(-h(1)+b)],
                      [mp.exp(1+a), mp.exp(h(1)+b)]]), (a, b)


def projected_iteration(B, y, pair, count):
    """Finite deterministic map. B must be an analytic fixture in this task."""
    if count < 0:
        raise ValueError("nonnegative iteration count")
    M, center = public_matrix(pair)
    A, c = M**-1, mp.matrix(center)
    radius = mp.mpf(1)/50
    for _ in range(count):
        proposed = c+A*(mp.matrix(y)-B(tuple(c)))
        c = mp.matrix([min(center[j]+radius, max(center[j]-radius, proposed[j]))
                       for j in range(2)])
    return tuple(c)


def kernel(x, a, b, component, tol):
    def density(y):
        if component == 1:
            return phi(y)
        t = hinv(y, tol)
        return phi(t)/hp(t)
    return density(x-a-b)-density(x-a)-density(x-b)+density(x)


def weights(z, coefficients, lattice, tol):
    alpha, beta = coefficients[:2], coefficients[2:]
    def A(x):
        return hp(x)*kernel(h(x), *beta, 2, tol)
    def L(x):
        return kernel(x, *alpha, 1, tol)-A(x)
    W = (mp.fsum(L(z+k) for k in range(1, lattice+1))
         if z >= -mp.mpf(7)/4 else -mp.fsum(L(z-k) for k in range(lattice)))
    return A(z)-W, W


def periodized(t, coefficients, lattice, tol):
    return mp.fsum(kernel(t+k, *coefficients[:2], 1, tol)
                   -hp(t+k)*kernel(h(t+k), *coefficients[2:], 2, tol)
                   for k in range(-lattice, lattice+1))


def monotone_knot_bracket(probability, knot, left, right, tolerance):
    """Diagnostic bisection; not a certified enclosure of inexact primitives."""
    if not probability(left) <= knot <= probability(right):
        raise ValueError("knot outside image")
    while right-left > tolerance:
        middle = (left+right)/2
        if probability(middle) < knot:
            left = middle
        else:
            right = middle
    return left, right


def gauss_two_piece(function, left, right, mesh, fourth_derivative_bound):
    """C4 piece only; returns value and exact-arithmetic quadrature bound.

    The caller must split at EVERY knot/jump and independently establish the
    fourth-derivative bound. This function does not estimate either obligation.
    """
    if mesh <= 0 or fourth_derivative_bound < 0 or right < left:
        raise ValueError("invalid panel parameters")
    if right == left:
        return mp.mpf(0), mp.mpf(0), 0
    panels = max(1, int(mp.ceil((right-left)/mesh)))
    width = (right-left)/panels
    total = mp.mpf(0)
    offset = width/(2*mp.sqrt(3))
    for j in range(panels):
        mid = left+(j+mp.mpf('.5'))*width
        total += width*(function(mid-offset)+function(mid+offset))/2
    return total, fourth_derivative_bound*(right-left)*width**4/4320, panels


def estimate_from_sources(groups):
    raise PermissionError("Sampled-data execution is disabled; analytic fixtures only")
