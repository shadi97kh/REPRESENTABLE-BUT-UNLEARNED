"""RNA non-crossing oracle against the independent enumerator.

Every assertion here compares the dynamic programme with `scmp.oracles.bruteforce`,
which never calls the recurrence. Exact-integer quantities are compared with `==`;
floating-point quantities are compared with a tolerance and are labelled diagnostics.
"""
import itertools
import math
import random

import pytest

from scmp.oracles import (EXACT_INTEGER, FLOAT_DIAGNOSTIC, INFEASIBLE, NEG_INF, PROVED,
                          RNANonCrossingOracle, derivation_edges, score_from_derivation)
from scmp.oracles.bruteforce import (brute_count, brute_log_partition, brute_marginals,
                                     brute_support, enumerate_rna, rna_is_valid)

TINY = ["", "A", "GC", "GGAAACC", "GCGCAAAAGCGC", "AUAUAUAUAU", "GGGAAAUCC",
        "ACGUACGUACGU", "GGGGAAAACCCC"]
SETTINGS = [(3, True), (3, False), (1, True), (0, False), (2, True)]


def _oracle(seq, min_loop, canon):
    return RNANonCrossingOracle(seq, min_loop=min_loop, canonical_only=canon)


# ------------------------------------------------------------------ counts / uniqueness
@pytest.mark.parametrize("seq", TINY)
@pytest.mark.parametrize("min_loop,canon", SETTINGS)
def test_count_matches_brute_force(seq, min_loop, canon):
    """If the grammar were ambiguous the count would exceed the true family size."""
    o = _oracle(seq, min_loop, canon)
    members = enumerate_rna(seq, min_loop, canon)
    c = o.count()
    assert c.numerical_status == EXACT_INTEGER
    assert c.value == brute_count(members)


@pytest.mark.parametrize("seq", TINY)
@pytest.mark.parametrize("min_loop,canon", SETTINGS)
def test_enumerated_members_are_valid_and_distinct(seq, min_loop, canon):
    members = enumerate_rna(seq, min_loop, canon)
    assert len(set(members)) == len(members)
    for m in members:
        assert rna_is_valid(m, seq, min_loop, canon)


def test_empty_sequence_and_empty_ground_set():
    """A family whose only member is the empty structure is size 1, not size 0."""
    for seq in ["", "A", "GC", "AAAA"]:
        o = _oracle(seq, 3, True)
        assert o.count().value == len(enumerate_rna(seq, 3, True))
    o = _oracle("AAAAAAAA", 3, True)          # canonical filter empties the ground set
    assert o.ground_set() == ()
    assert o.count().value == 1
    assert o.support({}).witness == ()
    assert o.log_partition({}).value == pytest.approx(0.0)


# ------------------------------------------------------------------------ support / max
@pytest.mark.parametrize("seq", TINY)
@pytest.mark.parametrize("min_loop,canon", SETTINGS)
def test_support_matches_brute_force_signed_integers(seq, min_loop, canon):
    o = _oracle(seq, min_loop, canon)
    members = enumerate_rna(seq, min_loop, canon)
    rng = random.Random(hash((seq, min_loop, canon)) & 0xFFFF)
    for _ in range(8):
        coef = {e: rng.randint(-9, 9) for e in o.ground_set()}
        r = o.support(coef)
        exp, _ = brute_support(members, coef)
        assert r.numerical_status == EXACT_INTEGER
        assert r.value == exp                      # exact integer comparison


@pytest.mark.parametrize("seq", ["GGAAACC", "GCGCAAAAGCGC", "GGGAAAUCC"])
def test_every_tiny_structure_is_selectable(seq):
    """Drive the objective at each member in turn; the argmax must be that member."""
    o = _oracle(seq, 3, True)
    members = enumerate_rna(seq, 3, True)
    ground = o.ground_set()
    for m in members:
        coef = {e: (10 if e in m else -10) for e in ground}
        r = o.support(coef)
        assert r.witness == tuple(sorted(m))
        assert r.value == 10 * len(m)


@pytest.mark.parametrize("seq", TINY)
@pytest.mark.parametrize("min_loop,canon", SETTINGS)
def test_witness_is_valid_and_scores_as_reported(seq, min_loop, canon):
    o = _oracle(seq, min_loop, canon)
    rng = random.Random(7)
    for _ in range(5):
        coef = {e: rng.randint(-5, 5) for e in o.ground_set()}
        r = o.support(coef)
        if r.termination_status == INFEASIBLE:
            continue
        assert rna_is_valid(r.witness, seq, min_loop, canon)
        assert sum(coef.get(e, 0) for e in r.witness) == r.value


# --------------------------------------------------- production-to-edge map P, uniqueness
@pytest.mark.parametrize("seq", ["GGAAACC", "GCGCAAAAGCGC", "GGGGAAAACCCC"])
def test_production_to_edge_map(seq):
    """P must recover exactly the witness, and the score must reconstruct through P."""
    o = _oracle(seq, 3, True)
    rng = random.Random(11)
    for _ in range(10):
        coef = {e: rng.randint(-6, 6) for e in o.ground_set()}
        r = o.support(coef)
        prods = r.stats["productions"]
        edges = derivation_edges(prods)
        assert len(edges) == len(set(edges))               # P injective on this derivation
        assert tuple(sorted(edges)) == r.witness           # edge multiset == structure
        assert score_from_derivation(prods, coef) == r.value
        for p in prods:
            assert p.kind in ("unpaired", "pair")
            assert (p.edge is None) == (p.kind == "unpaired")


@pytest.mark.parametrize("seq", ["GGAAACC", "GGGAAAUCC"])
def test_derivation_count_equals_member_count(seq):
    """Unambiguity, stated the other way round: one derivation per member."""
    o = _oracle(seq, 3, True)
    assert o.count().value == len(enumerate_rna(seq, 3, True))


# --------------------------------------------------------------------- log partition
@pytest.mark.parametrize("seq", TINY)
@pytest.mark.parametrize("min_loop,canon", SETTINGS)
def test_log_partition_matches_brute_force(seq, min_loop, canon):
    o = _oracle(seq, min_loop, canon)
    members = enumerate_rna(seq, min_loop, canon)
    rng = random.Random(3)
    for _ in range(5):
        coef = {e: rng.uniform(-2.0, 2.0) for e in o.ground_set()}
        r = o.log_partition(coef)
        assert r.numerical_status == FLOAT_DIAGNOSTIC        # never claimed exact
        assert r.value == pytest.approx(brute_log_partition(members, coef), rel=1e-9)


def test_log_partition_of_counts_recovers_count():
    """With all-zero coefficients, exp(logZ) is the family size."""
    for seq in ["GGAAACC", "GCGCAAAAGCGC"]:
        o = _oracle(seq, 3, True)
        z = o.log_partition({})
        assert math.exp(z.value) == pytest.approx(o.count().value, rel=1e-9)


# ---------------------------------------------------------------------------- marginals
@pytest.mark.parametrize("seq", ["GGAAACC", "GCGCAAAAGCGC", "GGGAAAUCC"])
def test_marginals_match_brute_force(seq):
    o = _oracle(seq, 3, True)
    members = enumerate_rna(seq, 3, True)
    rng = random.Random(5)
    coef = {e: rng.uniform(-1.5, 1.5) for e in o.ground_set()}
    got = o.marginals(coef)
    assert got.numerical_status == FLOAT_DIAGNOSTIC
    exp = brute_marginals(members, coef, o.ground_set())
    for e in o.ground_set():
        assert got.value[e] == pytest.approx(exp[e], abs=1e-9)


# ------------------------------------------------------------------- masks, exhaustively
@pytest.mark.parametrize("seq,min_loop,canon", [("GCGCAAAAGCGC", 3, True),
                                                ("GGGAAAUCC", 3, True),
                                                ("GGGGAAAACCCC", 3, True),
                                                ("GGAAACC", 1, True)])
def test_all_small_masks_consistent_and_inconsistent(seq, min_loop, canon):
    """Every forced/forbidden mask of size <= 2 each, consistent or not."""
    o = _oracle(seq, min_loop, canon)
    ground = o.ground_set()
    subsets = [()] + [(e,) for e in ground] + list(itertools.combinations(ground, 2))
    checked = 0
    for forced in subsets:
        for forbidden in subsets:
            dp = o.count(forced, forbidden)
            bf = len(enumerate_rna(seq, min_loop, canon, forced, forbidden))
            assert dp.value == bf, (forced, forbidden, dp.value, bf)
            feas = o.feasibility(forced, forbidden)
            assert bool(feas) == (bf > 0)
            assert feas.termination_status == (PROVED if bf > 0 else INFEASIBLE)
            checked += 1
    assert checked > 0


@pytest.mark.parametrize("seq", ["GCGCAAAAGCGC"])
def test_infeasible_masks_report_reason_and_status(seq):
    o = _oracle(seq, 3, True)
    g = o.ground_set()
    bad = [((g[0],), (g[0],), "both forced and forbidden"),
           (((0, 1),), (), "min_loop"),
           (((0, 99),), (), "out of range")]
    for forced, forbidden, needle in bad:
        f = o.feasibility(forced, forbidden)
        assert not f.feasible
        assert f.termination_status == INFEASIBLE
        assert needle in f.reason
        assert o.count(forced, forbidden).value == 0
        assert o.support({}, forced, forbidden).value == NEG_INF
        assert o.log_partition({}, forced, forbidden).value == NEG_INF


def test_forced_crossing_pairs_are_infeasible():
    seq = "GCGCAAAAGCGC"
    o = _oracle(seq, 3, True)
    g = o.ground_set()
    pairs = [(a, b) for a in g for b in g if a[0] < b[0] < a[1] < b[1]]
    assert pairs, "test needs a crossing pair to exist"
    a, b = pairs[0]
    f = o.feasibility((a, b))
    assert not f.feasible and "cross" in f.reason
    assert o.count((a, b)).value == 0
    assert len(enumerate_rna(seq, 3, True, (a, b))) == 0


def test_forced_shared_endpoint_is_zero_in_every_query_not_only_feasibility():
    """Regression: the recurrence built its partner map with a dict, so a second forced
    pair on the same index silently overwrote the first and the count came back
    non-zero. feasibility() caught it; count/support/log_partition did not."""
    seq = "GGGGAAAACCCC"
    o = _oracle(seq, 3, True)
    g = o.ground_set()
    shared = [(a, b) for a in g for b in g if a != b and (set(a) & set(b))]
    assert shared, "test needs two candidate pairs sharing an index"
    for a, b in shared[:12]:
        assert o.count((a, b)).value == 0
        assert o.support({}, (a, b)).value == NEG_INF
        assert o.log_partition({}, (a, b)).value == NEG_INF
        assert not o.feasibility((a, b)).feasible
        assert len(enumerate_rna(seq, 3, True, (a, b))) == 0


def test_forced_shared_endpoint_is_infeasible():
    seq = "GCGCAAAAGCGC"
    o = _oracle(seq, 3, True)
    g = o.ground_set()
    shared = [(a, b) for a in g for b in g if a != b and a[0] == b[0]]
    assert shared
    a, b = shared[0]
    f = o.feasibility((a, b))
    assert not f.feasible and "share index" in f.reason
    assert o.count((a, b)).value == 0


def test_forced_endpoint_partner_outside_interval():
    """A forced pair whose partner cannot be reached must kill the family, not silently
    produce a structure missing that pair."""
    seq = "GCGCAAAAGCGC"
    o = _oracle(seq, 3, True)
    e = o.ground_set()[0]
    r = o.support({}, forced=(e,))
    assert e in r.witness
    for m in enumerate_rna(seq, 3, True, forced=(e,)):
        assert e in m


# ------------------------------------------------------------------------- min-loop
@pytest.mark.parametrize("min_loop", [0, 1, 2, 3, 5])
def test_min_loop_is_respected(min_loop):
    seq = "GCGCAAAAGCGC"
    o = _oracle(seq, min_loop, True)
    for (i, j) in o.ground_set():
        assert j - i > min_loop
    assert o.count().value == len(enumerate_rna(seq, min_loop, True))
