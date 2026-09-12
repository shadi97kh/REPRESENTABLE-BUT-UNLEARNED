"""Validated rounding: rigorous enclosures for every certifying operation.

WHY NOT EPSILON PADDING. Adding a fudge factor to a final float and declaring the result
safe proves nothing: the size of the accumulated error is exactly what is unknown. Here
every elementary operation returns an INTERVAL that is guaranteed to contain the exact
real result, and intervals compose, so the final enclosure is a consequence of the
arithmetic rather than an assumption about it.

HOW. IEEE-754 round-to-nearest commits an error of at most half an ulp, so for any
elementary operation the exact result lies within one ulp either side of the computed
value. `math.nextafter` gives exactly that neighbour, so

    exact(a op b)  in  [nextafter(fl(a op b), -inf), nextafter(fl(a op b), +inf)]

is rigorous without needing access to the hardware rounding mode, which neither Python
nor numpy exposes. The enclosures are slightly wider than directed rounding would give;
they are never wrong.

`tests/numerics` checks the enclosures against exact `Fraction` arithmetic, which is the
only way to know an enclosure is valid rather than merely plausible.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

INF = float("inf")


def down(x: float) -> float:
    """Largest representable value <= the exact result of an operation returning x."""
    if x == -INF:
        return x
    return math.nextafter(x, -INF)


def up(x: float) -> float:
    if x == INF:
        return x
    return math.nextafter(x, INF)


@dataclass(frozen=True)
class Interval:
    lo: float
    hi: float

    def __post_init__(self):
        if not (self.lo <= self.hi):
            raise ValueError(f"degenerate interval [{self.lo}, {self.hi}]")

    # -- constructors -----------------------------------------------------
    @staticmethod
    def exact(x: float) -> "Interval":
        """A float is exactly itself; no widening until an operation is applied."""
        return Interval(float(x), float(x))

    @staticmethod
    def hull(*vals) -> "Interval":
        los = [v.lo if isinstance(v, Interval) else float(v) for v in vals]
        his = [v.hi if isinstance(v, Interval) else float(v) for v in vals]
        return Interval(min(los), max(his))

    # -- rigorous arithmetic ---------------------------------------------
    def __add__(self, other) -> "Interval":
        o = other if isinstance(other, Interval) else Interval.exact(other)
        return Interval(down(self.lo + o.lo), up(self.hi + o.hi))

    __radd__ = __add__

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other) -> "Interval":
        o = other if isinstance(other, Interval) else Interval.exact(other)
        return self + (-o)

    def __mul__(self, other) -> "Interval":
        o = other if isinstance(other, Interval) else Interval.exact(other)
        ps = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Interval(down(min(ps)), up(max(ps)))

    __rmul__ = __mul__

    def __truediv__(self, other) -> "Interval":
        o = other if isinstance(other, Interval) else Interval.exact(other)
        if o.lo <= 0.0 <= o.hi:
            raise ZeroDivisionError("interval divisor straddles zero")
        qs = (self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi)
        return Interval(down(min(qs)), up(max(qs)))

    # -- predicates -------------------------------------------------------
    def contains(self, x) -> bool:
        if isinstance(x, Fraction):
            return Fraction(self.lo) <= x <= Fraction(self.hi)
        return self.lo <= x <= self.hi

    @property
    def width(self) -> float:
        return up(self.hi - self.lo)

    def __repr__(self):
        return f"[{self.lo!r}, {self.hi!r}]"


def isum(values) -> Interval:
    """Rigorously enclosed sum. Left-to-right, each step widened."""
    acc = Interval.exact(0.0)
    for v in values:
        acc = acc + (v if isinstance(v, Interval) else Interval.exact(v))
    return acc


def idot(a, b) -> Interval:
    """Rigorously enclosed dot product."""
    return isum((x if isinstance(x, Interval) else Interval.exact(x)) *
                (y if isinstance(y, Interval) else Interval.exact(y))
                for x, y in zip(a, b))


def imax(values) -> Interval:
    vals = [v if isinstance(v, Interval) else Interval.exact(v) for v in values]
    if not vals:
        return Interval(-INF, -INF)
    return Interval(max(v.lo for v in vals), max(v.hi for v in vals))


# ------------------------------------------------------------------ exact reference
def exact_affine(coef, const, z) -> Fraction:
    """Ground truth for the tests: the affine form evaluated in exact rationals."""
    total = Fraction(const)
    for c, zi in zip(coef, z):
        total += Fraction(c) * Fraction(zi)
    return total


# ------------------------------------------------------------------- status handling
FLOAT_UNVERIFIED = "float_unverified"
FLOAT_ENCLOSED = "float_enclosed"
EXACT_RATIONAL = "exact_rational"


def combine_status(*statuses) -> str:
    """The weakest status wins. A single unenclosed step contaminates the whole chain,
    which is the point: status is derived, never asserted."""
    order = [EXACT_RATIONAL, FLOAT_ENCLOSED, FLOAT_UNVERIFIED]
    worst = 0
    for s in statuses:
        if s in order:
            worst = max(worst, order.index(s))
        else:
            worst = len(order) - 1
    return order[worst]
