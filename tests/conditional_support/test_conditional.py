"""Deterministic tests for conditional message support.

No learned component is exercised here: everything below must hold before a policy is
allowed anywhere near the bound. Soundness is checked against exhaustive enumeration of
the feasible family, and the conditional bound is checked against the unconditional
baseline it must never be looser than.
"""
import numpy as np
import pytest
import torch

from scmp.bounds import INFEASIBLE, ResidualBounder, certified_upper_bound
from scmp.conditional import (BoundCache, CacheKey, ConditionalBounder, FamilyConc,
                              SCORERS, conditional_upper_bound)
from scmp.model import CertifiedModel, backbone_adjacency
from scmp.oracles import RNANonCrossingOracle
from scmp.oracles.bruteforce import enumerate_rna

SLACK = 1e-9
SEQS = ["GGGAAAUCC", "GCGCAAAAGCGC"]


def onehot(seq):
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return X


def build(seq, gamma=1.0, seed=0, d_hid=5, scale=1.0, pseed=0):
    o = RNANonCrossingOracle(seq, 3, True)
    m = CertifiedModel(4, o.ground_set(), len(seq), d_hid=d_hid, gamma=gamma, seed=seed)
    if scale:
        g = torch.Generator().manual_seed(pseed)
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype) * scale)
    return o, m, onehot(seq), backbone_adjacency(len(seq))


def cb(o, m, X, B, **kw):
    kw.setdefault("budget", 6)
    return ConditionalBounder(m, X, B, o, **kw)


# ------------------------------------------------------------------------- soundness
@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("pseed", [0, 1, 2])
@pytest.mark.parametrize("scorer", ["random", "largest_gap", "measured_gain"])
def test_conditional_bound_dominates_the_true_maximum(seq, pseed, scorer):
    o, m, X, B = build(seq, pseed=100 + pseed)
    members = enumerate_rna(seq, 3, True)
    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)
    r = conditional_upper_bound(m, X, o, B, budget=5, scorer=scorer, seed=pseed)
    assert r.upper >= true_max - SLACK


@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("pseed", [0, 1, 2])
def test_every_intermediate_bracket_still_holds(seq, pseed):
    o, m, X, B = build(seq, pseed=200 + pseed)
    bd = cb(o, m, X, B)
    lo, hi, tr, qb = bd.propagate()
    ref = ResidualBounder(m, X, B)
    for s in enumerate_rna(seq, 3, True):
        z = np.zeros(len(m.candidates))
        for e in s:
            z[m.index[e]] = 1.0
        truth = ref.reference_forward(z)
        for name, (L, H) in tr.layers.items():
            arr = truth[name]
            for i in range(len(L)):
                for k in range(len(L[i])):
                    assert L[i][k].eval(z) <= float(arr[i][k]) + SLACK
                    assert float(arr[i][k]) <= H[i][k].eval(z) + SLACK


@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("pseed", [0, 1, 2, 3])
def test_conditional_is_never_looser_than_the_baseline(seq, pseed):
    """The min of valid bounds is retained, so a conditional pass cannot hurt."""
    o, m, X, B = build(seq, pseed=300 + pseed)
    base = certified_upper_bound(m, X, o, B, mode="combined")
    cond = conditional_upper_bound(m, X, o, B, budget=6, seed=pseed)
    assert cond.upper <= base.upper + SLACK


def test_conditional_is_strictly_tighter_somewhere():
    seq = "GCGCAAAAGCGC"
    strict = 0
    for pseed in range(5):
        o, m, X, B = build(seq, pseed=400 + pseed)
        base = certified_upper_bound(m, X, o, B, mode="combined")
        cond = conditional_upper_bound(m, X, o, B, budget=8, seed=pseed)
        strict += cond.upper < base.upper - 1e-9
    assert strict > 0


# ---------------------------------------------------------- ignore only when proved
def test_edge_is_ignored_only_after_being_proved_impossible():
    seq = "GCGCAAAAGCGC"
    o, m, X, B = build(seq, pseed=7)
    g = o.ground_set()
    cross = [(a, b) for a in g for b in g if a[0] < b[0] < a[1] < b[1]]
    assert cross, "need a crossing pair"
    e0, e1 = cross[0]
    bd = cb(o, m, X, B, forced=(e0,), budget=len(g))
    uncond = bd._upto_z1(bd.forced, bd.forbidden)
    q = bd.conditional_message_support(e1[0], 0, e1, uncond)
    assert q.proved_impossible and q.source == "impossible"
    assert not o.feasibility(forced=(e0, e1)).feasible      # the proof itself


def test_unqueried_edges_keep_their_unconditional_bound():
    """An edge we did not spend budget on must not be dropped."""
    seq = "GCGCAAAAGCGC"
    o, m, X, B = build(seq, pseed=8)
    small = conditional_upper_bound(m, X, o, B, budget=1, seed=0)
    none = certified_upper_bound(m, X, o, B, mode="combined")
    assert small.upper <= none.upper + SLACK
    assert small.stats["conditional_passes"] <= 1 + len(o.ground_set())


def test_fallback_is_used_when_conditioning_does_not_help():
    seq = "GGGAAAUCC"
    seen = set()
    for pseed in range(6):
        o, m, X, B = build(seq, pseed=500 + pseed)
        bd = cb(o, m, X, B, budget=len(o.ground_set()))
        bd.propagate()
        if bd.report["fallback_used"]:
            seen.add("fallback")
        if bd.report["conditional_used"]:
            seen.add("conditional")
    assert "fallback" in seen and "conditional" in seen


# ------------------------------------------------------------- separation examples
def test_carrier_gating_separation():
    """Conditioning on e forbids everything that crosses it or shares an endpoint, so
    the carried message shrinks. Demonstrated as a strictly tighter q for some
    (node, direction, edge)."""
    seq = "GCGCAAAAGCGC"
    found = None
    for pseed in range(8):
        o, m, X, B = build(seq, pseed=600 + pseed)
        bd = cb(o, m, X, B, budget=len(o.ground_set()))
        _, _, _, qb = bd.propagate()
        tight = [q for q in qb.values() if q.source == "conditional"]
        if tight:
            found = max(tight, key=lambda q: q.unconditional - q.conditional)
            break
    assert found is not None, "no conditional tightening observed"
    assert found.conditional < found.unconditional
    assert found.value == found.conditional          # min of valid bounds


def test_relu_separation_conditioning_stabilises_units():
    """The tighter conditioned interval can turn an unstable ReLU stable, which is the
    other source of separation and is independent of gating."""
    seq = "GCGCAAAAGCGC"
    improved = 0
    for pseed in range(8):
        o, m, X, B = build(seq, pseed=700 + pseed)
        bd = cb(o, m, X, B, budget=0)
        uncond = bd._upto_z1(bd.forced, bd.forbidden)
        n_unstable = sum(r.count("unstable") for r in uncond["relu_kinds"])
        for e in o.ground_set():
            cond = bd._conditional_z1(e)
            if cond is None:
                continue
            if sum(r.count("unstable") for r in cond["relu_kinds"]) < n_unstable:
                improved += 1
                break
        if improved:
            break
    assert improved > 0, "conditioning never stabilised a ReLU"


# --------------------------------------------------------------------------- caching
def test_cache_key_covers_model_family_mask_layer_node_direction_condition():
    fields = CacheKey.__dataclass_fields__
    assert set(fields) == {"model_hash", "family_hash", "mask", "layer", "node",
                           "direction", "condition"}


def test_cache_hits_and_is_invalidated_by_a_model_change():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, pseed=9)
    cache = BoundCache()
    conditional_upper_bound(m, X, o, B, budget=4, seed=0, cache=cache)
    first = cache.stats()["entries"]
    conditional_upper_bound(m, X, o, B, budget=4, seed=0, cache=cache)
    assert cache.stats()["hits"] > 0
    assert cache.stats()["entries"] == first          # nothing recomputed
    with torch.no_grad():
        m.g.lin0.weight.add_(1.0)                      # different model
    conditional_upper_bound(m, X, o, B, budget=4, seed=0, cache=cache)
    assert cache.stats()["entries"] > first            # new keys, no stale reuse


def test_cache_separates_masks_and_conditions():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, pseed=10)
    bd = cb(o, m, X, B, budget=0)
    a = bd._upto_z1(bd.forced, bd.forbidden)
    e = o.ground_set()[0]
    b = bd._upto_z1((e,), bd.forbidden, cond=(e,))
    assert a is not b


# ------------------------------------------------------- recursion and bounded work
def test_conditional_recursion_reaches_earlier_layers():
    """The conditional Z1 must differ from the unconditional one: if it did not, no
    earlier layer was re-evaluated under the condition."""
    seq = "GCGCAAAAGCGC"
    o, m, X, B = build(seq, pseed=11)
    bd = cb(o, m, X, B, budget=0)
    uncond = bd._upto_z1(bd.forced, bd.forbidden)
    differs = 0
    for e in o.ground_set():
        cond = bd._conditional_z1(e)
        if cond is None:
            continue
        if any(cond["hi"][u][k] < uncond["hi"][u][k] - 1e-12
               for u in range(m.n_nodes) for k in range(len(uncond["hi"][0]))):
            differs += 1
    assert differs > 0


@pytest.mark.parametrize("budget", [0, 1, 3, 8])
def test_recursive_work_is_bounded_by_the_budget(budget):
    seq = "GCGCAAAAGCGC"
    o, m, X, B = build(seq, pseed=12)
    bd = cb(o, m, X, B, budget=budget, scorer="largest_gap")
    bd.propagate()
    assert bd.report["conditional_passes"] <= budget
    assert bd.report["budget"] == budget


def test_coordinate_directions_are_retained():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, pseed=13)
    bd = cb(o, m, X, B, budget=len(o.ground_set()))
    _, _, _, qb = bd.propagate()
    dirs = {q.direction for q in qb.values()}
    assert dirs == set(range(5))                      # every coordinate, not a mixture
    for (i, k, e), q in qb.items():
        assert q.node == i and q.direction == k and q.edge == e


def test_baseline_certificate_is_retained_alongside():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, pseed=14)
    bd = cb(o, m, X, B, budget=len(o.ground_set()))
    _, _, _, qb = bd.propagate()
    for q in qb.values():
        assert np.isfinite(q.unconditional)           # the inherited bound is kept
        assert q.value <= q.unconditional + SLACK
        assert q.assumptions and q.numerical_status == "float_unverified"


# ------------------------------------------------------------------ affine combining
def test_affine_term_is_combined_before_the_support_call():
    """One oracle call on the combined coefficient vector, not separate maxima."""
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, pseed=15)
    calls = {"n": 0}
    orig = o.support

    def counting(coef, forced=(), forbidden=()):
        calls["n"] += 1
        counting.last = dict(coef)
        return orig(coef, forced, forbidden)

    o.support = counting
    try:
        r = conditional_upper_bound(m, X, o, B, budget=0, seed=0)
    finally:
        o.support = orig
    a = m.additive_coefficients(X)
    combined_seen = any(abs(counting.last[e] - a[e]) > 1e-12 for e in m.candidates)
    assert combined_seen, "final support call used a, not a + gamma*c"
    assert r.stats["final_oracle_calls"] >= 1


# ----------------------------------------------------------------- scorer comparison
@pytest.mark.parametrize("scorer", list(SCORERS))
def test_all_scorers_are_sound_and_respect_budget(scorer):
    seq = "GCGCAAAAGCGC"
    o, m, X, B = build(seq, pseed=16)
    members = enumerate_rna(seq, 3, True)
    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)
    r = conditional_upper_bound(m, X, o, B, budget=4, scorer=scorer, seed=1)
    assert r.upper >= true_max - SLACK


def test_measured_gain_reports_its_teacher_cost():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, pseed=17)
    bd = cb(o, m, X, B, budget=3, scorer="measured_gain")
    bd.propagate()
    assert bd.report["teacher_passes"] > 0            # cost is counted, not hidden


def test_gamma_zero_conditional_bound_is_still_exact():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=0.0, pseed=18)
    members = enumerate_rna(seq, 3, True)
    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)
    r = conditional_upper_bound(m, X, o, B, budget=4)
    assert r.upper == pytest.approx(true_max, abs=1e-12)


def test_infeasible_mask_returns_no_bound():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, pseed=19)
    e = o.ground_set()[0]
    r = conditional_upper_bound(m, X, o, B, forced=(e,), forbidden=(e,))
    assert r.termination_status == INFEASIBLE and r.upper == float("-inf")
