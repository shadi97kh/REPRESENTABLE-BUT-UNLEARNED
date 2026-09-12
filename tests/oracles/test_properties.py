"""Cross-family property tests.

These hold for any family behind the interface, so both oracles are run through the same
assertions. Finite-difference tests are explicitly floating-point diagnostics and are
kept separate from the exact-integer assertions elsewhere.
"""
import math
import random

import pytest

from scmp.oracles import (EXACT_INTEGER, FLOAT_DIAGNOSTIC, FLOAT_UNVERIFIED, INFEASIBLE,
                          NEG_INF, PROVED, LayeredDAGPathOracle, RNANonCrossingOracle)
from scmp.oracles.bruteforce import enumerate_dag_paths, enumerate_rna

RNA_CASES = ["GGAAACC", "GCGCAAAAGCGC", "GGGAAAUCC", "GGGGAAAACCCC"]
DAG_SIZES = [1, 3, 2, 2]
DAG_EDGES = [(0, 0, 0), (0, 0, 1), (0, 0, 2), (1, 0, 0), (1, 1, 0),
             (1, 1, 1), (1, 2, 1), (2, 0, 0), (2, 1, 1), (2, 0, 1)]


def oracles():
    out = [RNANonCrossingOracle(s, 3, True) for s in RNA_CASES]
    out.append(LayeredDAGPathOracle(DAG_SIZES, DAG_EDGES))
    return out


# ------------------------------------------------------- finite-difference marginals
@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_marginals_equal_gradient_of_log_partition(o):
    """FLOAT DIAGNOSTIC. marginal(e) must equal d logZ / d c_e.

    Central differences, h = 1e-5. This is an independent route to the marginals: the
    implementation computes Z(forced+{e})/Z, which never differentiates anything.
    """
    ground = o.ground_set()
    if not ground or o.count().value == 0:
        pytest.skip("degenerate family")
    rng = random.Random(20260905)
    coef = {e: rng.uniform(-1.0, 1.0) for e in ground}
    marg = o.marginals(coef).value
    h = 1e-5
    for e in ground:
        up = dict(coef); up[e] += h
        dn = dict(coef); dn[e] -= h
        zu = o.log_partition(up).value
        zd = o.log_partition(dn).value
        if zu == NEG_INF or zd == NEG_INF:
            continue
        fd = (zu - zd) / (2 * h)
        assert fd == pytest.approx(marg[e], abs=2e-5), (e, fd, marg[e])


@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_marginals_are_probabilities(o):
    if not o.ground_set() or o.count().value == 0:
        pytest.skip("degenerate family")
    rng = random.Random(1)
    coef = {e: rng.uniform(-2, 2) for e in o.ground_set()}
    for e, p in o.marginals(coef).value.items():
        assert -1e-12 <= p <= 1 + 1e-12


@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_forced_and_forbidden_marginals_are_pinned(o):
    ground = o.ground_set()
    if not ground or o.count().value == 0:
        pytest.skip("degenerate family")
    e = ground[0]
    coef = {g: 0.0 for g in ground}
    if o.feasibility(forced=(e,)):
        assert o.marginals(coef, forced=(e,)).value[e] == pytest.approx(1.0)
    if o.feasibility(forbidden=(e,)):
        assert o.marginals(coef, forbidden=(e,)).value[e] == pytest.approx(0.0)


# ------------------------------------------------------------------ conditioning algebra
@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_condition_composes(o):
    ground = o.ground_set()
    if len(ground) < 2:
        pytest.skip("need two edges")
    a, b = ground[0], ground[1]
    once = o.condition(forced=(a,)).condition(forbidden=(b,))
    direct = o.count(forced=(a,), forbidden=(b,))
    assert once.count().value == direct.value


@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_condition_is_equivalent_to_passing_the_mask(o):
    ground = o.ground_set()
    if not ground:
        pytest.skip("empty ground set")
    e = ground[0]
    rng = random.Random(2)
    coef = {g: rng.randint(-4, 4) for g in ground}
    assert o.condition(forced=(e,)).support(coef).value == o.support(coef, forced=(e,)).value
    assert o.condition(forced=(e,)).count().value == o.count(forced=(e,)).value


@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_conditioning_can_only_shrink_the_family(o):
    ground = o.ground_set()
    if not ground:
        pytest.skip("empty ground set")
    total = o.count().value
    for e in ground:
        assert o.count(forced=(e,)).value <= total
        assert o.count(forbidden=(e,)).value <= total
        # forced + forbidden on the same edge partitions the family
        assert o.count(forced=(e,)).value + o.count(forbidden=(e,)).value == total


# --------------------------------------------------------------- exact vs float status
@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_numerical_status_is_honest(o):
    ground = o.ground_set()
    int_coef = {e: 1 for e in ground}
    flt_coef = {e: 1.5 for e in ground}
    assert o.count().numerical_status == EXACT_INTEGER
    assert o.support(int_coef).numerical_status == EXACT_INTEGER
    if ground:
        assert o.support(flt_coef).numerical_status == FLOAT_UNVERIFIED
    assert o.log_partition(int_coef).numerical_status == FLOAT_DIAGNOSTIC
    if o.count().value > 0 and ground:
        assert o.marginals(int_coef).numerical_status == FLOAT_DIAGNOSTIC


@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_integer_support_returns_an_integer(o):
    ground = o.ground_set()
    if not ground or o.count().value == 0:
        pytest.skip("degenerate family")
    coef = {e: 3 for e in ground}
    r = o.support(coef)
    assert isinstance(r.value, int)          # not a float that happens to be integral


# ------------------------------------------------------------------ certificate fields
@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_every_result_carries_the_certificate_fields(o):
    coef = {e: 1 for e in o.ground_set()}
    for r in (o.count(), o.support(coef), o.log_partition(coef)):
        assert r.assumptions, "assumptions must be recorded"
        assert r.numerical_status in (EXACT_INTEGER, FLOAT_DIAGNOSTIC, FLOAT_UNVERIFIED)
        assert r.termination_status in (PROVED, INFEASIBLE)
        assert isinstance(r.family_hash, str) and len(r.family_hash) == 16
        assert r.bound_direction in ("exact", "upper", "lower", "none")


@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_family_hash_tracks_the_mask(o):
    ground = o.ground_set()
    if not ground:
        pytest.skip("empty ground set")
    e = ground[0]
    assert o.count().family_hash != o.count(forced=(e,)).family_hash
    assert o.count(forced=(e,)).family_hash != o.count(forbidden=(e,)).family_hash
    assert o.count(forced=(e,)).family_hash == o.condition(forced=(e,)).count().family_hash


@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_stats_report_chart_and_productions(o):
    s = o.count().stats
    assert s["chart_values"] > 0
    assert s["materialized_productions"] >= 0


# --------------------------------------------------------------- support/partition bounds
@pytest.mark.parametrize("o", oracles(), ids=lambda o: o.family_id + str(len(o.ground_set())))
def test_support_is_bounded_by_log_partition(o):
    """max score <= logZ <= max score + log|F| for any family."""
    ground = o.ground_set()
    n = o.count().value
    if n == 0 or not ground:
        pytest.skip("degenerate family")
    rng = random.Random(8)
    coef = {e: rng.uniform(-1, 1) for e in ground}
    mx = o.support(coef).value
    lz = o.log_partition(coef).value
    assert mx - 1e-9 <= lz <= mx + math.log(n) + 1e-9
