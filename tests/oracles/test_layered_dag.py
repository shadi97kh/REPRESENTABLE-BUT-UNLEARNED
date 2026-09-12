"""Layered-DAG path oracle against the independent enumerator.

Second family behind the same interface. Its purpose is to catch assumptions that leaked
from RNA into `base.py`: here P is the identity, every member has the same cardinality,
and forced edges interact through layer chaining rather than nesting.
"""
import itertools
import math
import random

import pytest

from scmp.oracles import (EXACT_INTEGER, FLOAT_DIAGNOSTIC, INFEASIBLE, NEG_INF, PROVED,
                          LayeredDAGPathOracle, derivation_edges, score_from_derivation)
from scmp.oracles.bruteforce import (brute_log_partition, brute_marginals, brute_support,
                                     enumerate_dag_paths)


def make(seed, n_layers=4, width=3, density=0.7):
    rng = random.Random(seed)
    sizes = [1] + [rng.randint(1, width) for _ in range(n_layers - 1)]
    edges = []
    for l in range(len(sizes) - 1):
        for u in range(sizes[l]):
            for v in range(sizes[l + 1]):
                if rng.random() < density:
                    edges.append((l, u, v))
    return sizes, edges


CASES = [make(s) for s in range(12)]


@pytest.mark.parametrize("sizes,edges", CASES)
def test_count_matches_brute_force(sizes, edges):
    o = LayeredDAGPathOracle(sizes, edges)
    paths = enumerate_dag_paths(sizes, edges)
    c = o.count()
    assert c.numerical_status == EXACT_INTEGER
    assert c.value == len(paths)


@pytest.mark.parametrize("sizes,edges", CASES)
def test_support_matches_brute_force_signed_integers(sizes, edges):
    o = LayeredDAGPathOracle(sizes, edges)
    paths = enumerate_dag_paths(sizes, edges)
    rng = random.Random(hash((tuple(sizes), len(edges))) & 0xFFFF)
    for _ in range(6):
        coef = {e: rng.randint(-9, 9) for e in o.ground_set()}
        r = o.support(coef)
        exp, _ = brute_support(paths, coef)
        assert r.numerical_status == EXACT_INTEGER
        if not paths:
            assert r.value == NEG_INF and r.termination_status == INFEASIBLE
        else:
            assert r.value == exp


@pytest.mark.parametrize("sizes,edges", CASES)
def test_witness_is_a_real_path(sizes, edges):
    o = LayeredDAGPathOracle(sizes, edges)
    paths = set(enumerate_dag_paths(sizes, edges))
    if not paths:
        pytest.skip("no path in this instance")
    rng = random.Random(3)
    coef = {e: rng.randint(-4, 4) for e in o.ground_set()}
    r = o.support(coef)
    assert r.witness in paths
    assert sum(coef[e] for e in r.witness) == r.value
    assert len(r.witness) == len(sizes) - 1               # one edge per transition


@pytest.mark.parametrize("sizes,edges", CASES)
def test_production_to_edge_map_is_identity(sizes, edges):
    o = LayeredDAGPathOracle(sizes, edges)
    if not enumerate_dag_paths(sizes, edges):
        pytest.skip("no path in this instance")
    rng = random.Random(9)
    coef = {e: rng.randint(-5, 5) for e in o.ground_set()}
    r = o.support(coef)
    prods = r.stats["productions"]
    assert derivation_edges(prods) == list(r.witness)
    assert score_from_derivation(prods, coef) == r.value
    assert all(p.edge is not None for p in prods)


@pytest.mark.parametrize("sizes,edges", CASES)
def test_log_partition_matches_brute_force(sizes, edges):
    o = LayeredDAGPathOracle(sizes, edges)
    paths = enumerate_dag_paths(sizes, edges)
    rng = random.Random(4)
    coef = {e: rng.uniform(-2, 2) for e in o.ground_set()}
    r = o.log_partition(coef)
    assert r.numerical_status == FLOAT_DIAGNOSTIC
    if not paths:
        assert r.value == NEG_INF
    else:
        assert r.value == pytest.approx(brute_log_partition(paths, coef), rel=1e-9)


@pytest.mark.parametrize("sizes,edges", CASES[:6])
def test_marginals_match_brute_force(sizes, edges):
    o = LayeredDAGPathOracle(sizes, edges)
    paths = enumerate_dag_paths(sizes, edges)
    if not paths:
        pytest.skip("no path in this instance")
    rng = random.Random(6)
    coef = {e: rng.uniform(-1.5, 1.5) for e in o.ground_set()}
    got = o.marginals(coef).value
    exp = brute_marginals(paths, coef, o.ground_set())
    for e in o.ground_set():
        assert got[e] == pytest.approx(exp[e], abs=1e-9)


@pytest.mark.parametrize("sizes,edges", CASES[:6])
def test_all_small_masks(sizes, edges):
    o = LayeredDAGPathOracle(sizes, edges)
    g = o.ground_set()
    subsets = [()] + [(e,) for e in g] + list(itertools.combinations(g, 2))
    for forced in subsets:
        for forbidden in subsets[:len(g) + 1]:
            dp = o.count(forced, forbidden)
            bf = len(enumerate_dag_paths(sizes, edges, forced, forbidden))
            assert dp.value == bf, (forced, forbidden, dp.value, bf)
            assert bool(o.feasibility(forced, forbidden)) == (bf > 0)


def test_forced_edges_in_same_layer_are_infeasible():
    sizes = [1, 2, 2]
    edges = [(0, 0, 0), (0, 0, 1), (1, 0, 0), (1, 1, 1)]
    o = LayeredDAGPathOracle(sizes, edges)
    f = o.feasibility(((0, 0, 0), (0, 0, 1)))
    assert not f.feasible and "share layer" in f.reason
    assert o.count(((0, 0, 0), (0, 0, 1))).value == 0


def test_forced_edges_that_do_not_meet_are_infeasible():
    sizes = [1, 2, 2]
    edges = [(0, 0, 0), (0, 0, 1), (1, 0, 0), (1, 1, 1)]
    o = LayeredDAGPathOracle(sizes, edges)
    f = o.feasibility(((0, 0, 0), (1, 1, 1)))          # ends at 0, starts at 1
    assert not f.feasible and "do not meet" in f.reason
    assert o.count(((0, 0, 0), (1, 1, 1))).value == 0
    assert len(enumerate_dag_paths(sizes, edges, ((0, 0, 0), (1, 1, 1)))) == 0


def test_empty_family_when_no_edges_reach_the_sink():
    sizes = [1, 2, 2]
    edges = [(0, 0, 0)]                                 # nothing leaves layer 1
    o = LayeredDAGPathOracle(sizes, edges)
    assert o.count().value == 0
    assert len(enumerate_dag_paths(sizes, edges)) == 0
    assert o.support({}).value == NEG_INF
    assert o.support({}).termination_status == INFEASIBLE
    assert not o.feasibility()
