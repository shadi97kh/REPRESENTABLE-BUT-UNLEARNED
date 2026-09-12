"""Validated rounding: enclosures checked against exact rational arithmetic.

An enclosure is only worth something if it PROVABLY contains the exact result. Every
test here computes ground truth with `fractions.Fraction`, which is exact, and asserts
containment. Nothing is padded by a tolerance: a tolerance is what these tests exist to
avoid needing.
"""
import math
import random
from fractions import Fraction

import pytest

from scmp.numerics import (EXACT_RATIONAL, FLOAT_ENCLOSED, FLOAT_UNVERIFIED, Interval,
                           combine_status, down, exact_affine, idot, imax, isum, up)


def test_down_and_up_straddle_the_value():
    for x in [0.0, 1.0, -1.0, 1e-300, 1e300, math.pi]:
        assert down(x) <= x <= up(x)
        if x != 0.0:
            assert down(x) < x < up(x)


def test_infinities_are_absorbed():
    assert up(float("inf")) == float("inf")
    assert down(float("-inf")) == float("-inf")


@pytest.mark.parametrize("seed", range(20))
def test_sum_enclosure_contains_the_exact_rational_sum(seed):
    rng = random.Random(seed)
    vals = [rng.uniform(-1e3, 1e3) for _ in range(50)]
    enc = isum(vals)
    exact = sum((Fraction(v) for v in vals), Fraction(0))
    assert enc.contains(exact), (enc, float(exact))


@pytest.mark.parametrize("seed", range(20))
def test_dot_enclosure_contains_the_exact_rational_dot(seed):
    rng = random.Random(seed)
    a = [rng.uniform(-50, 50) for _ in range(30)]
    b = [rng.uniform(-50, 50) for _ in range(30)]
    enc = idot(a, b)
    exact = sum((Fraction(x) * Fraction(y) for x, y in zip(a, b)), Fraction(0))
    assert enc.contains(exact)


@pytest.mark.parametrize("seed", range(20))
def test_affine_enclosure_contains_the_exact_value(seed):
    rng = random.Random(seed)
    m = 12
    coef = [rng.uniform(-10, 10) for _ in range(m)]
    const = rng.uniform(-10, 10)
    z = [rng.choice([0.0, 1.0]) for _ in range(m)]
    enc = idot(coef, z) + Interval.exact(const)
    assert enc.contains(exact_affine(coef, const, z))


def test_products_enclose_all_four_endpoint_cases():
    for a in (Interval(-2.0, 3.0), Interval(1.0, 4.0), Interval(-5.0, -1.0)):
        for b in (Interval(-1.0, 2.0), Interval(0.5, 6.0), Interval(-7.0, -2.0)):
            p = a * b
            for x in (a.lo, a.hi, 0.5 * (a.lo + a.hi)):
                for y in (b.lo, b.hi, 0.5 * (b.lo + b.hi)):
                    if a.lo <= x <= a.hi and b.lo <= y <= b.hi:
                        assert p.contains(Fraction(x) * Fraction(y))


def test_division_by_a_straddling_interval_is_refused():
    with pytest.raises(ZeroDivisionError):
        Interval(1.0, 2.0) / Interval(-1.0, 1.0)


def test_degenerate_interval_is_refused():
    with pytest.raises(ValueError):
        Interval(1.0, 0.0)


def test_width_is_non_negative_and_grows_under_accumulation():
    a = isum([0.1] * 1000)
    assert a.width >= 0.0
    assert a.contains(Fraction(1, 10) * 1000)


def test_imax_encloses_the_true_maximum():
    vals = [Interval(-1.0, 0.0), Interval(2.0, 3.0), Interval(1.0, 5.0)]
    m = imax(vals)
    assert m.lo == 2.0 and m.hi == 5.0


def test_status_combination_takes_the_weakest():
    assert combine_status(EXACT_RATIONAL, EXACT_RATIONAL) == EXACT_RATIONAL
    assert combine_status(EXACT_RATIONAL, FLOAT_ENCLOSED) == FLOAT_ENCLOSED
    assert combine_status(FLOAT_ENCLOSED, FLOAT_UNVERIFIED) == FLOAT_UNVERIFIED
    assert combine_status(EXACT_RATIONAL, "something_unknown") == FLOAT_UNVERIFIED


def test_a_single_unenclosed_step_contaminates_the_chain():
    """Status is derived from what actually happened, never asserted."""
    assert combine_status(*([EXACT_RATIONAL] * 99 + [FLOAT_UNVERIFIED])) \
        == FLOAT_UNVERIFIED


def test_enclosure_is_not_a_tolerance():
    """The enclosure must be derived per operation, so widening a value by a fixed
    epsilon is NOT equivalent and must not be mistaken for it."""
    vals = [1e16, 1.0, -1e16]
    enc = isum(vals)
    exact = sum((Fraction(v) for v in vals), Fraction(0))
    assert enc.contains(exact)
    naive = 1e-9
    assert enc.width > naive          # cancellation: a fixed epsilon would be too small
