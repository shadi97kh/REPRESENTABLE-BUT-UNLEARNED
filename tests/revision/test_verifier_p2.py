"""P2 mechanism tests: the ladder, the instrumentation and the v2 metrics.

These pin the properties the comparison depends on. In particular they check that
the difficulty ladder's generation rule cannot drift without its hash changing,
and that instrumentation does not alter any returned bound.
"""
import numpy as np
import pytest
import torch

from scmp.revision import ladder as L
from scmp.revision.instrument import COMPONENTS, Meter


# --------------------------------------------------------------- frozen ladder
def test_rule_hash_is_stable_and_covers_the_whole_rule():
    assert L.rule_hash() == L.RULE_HASH
    import copy
    saved = copy.deepcopy(L.LADDER)
    try:
        L.LADDER["rungs"][0]["d_hid"] += 1
        assert L.rule_hash() != L.RULE_HASH, "editing the rule must change its hash"
    finally:
        L.LADDER.clear()
        L.LADDER.update(saved)
    assert L.rule_hash() == L.RULE_HASH


def test_generation_is_deterministic_and_total():
    a, b = L.generate(), L.generate()
    assert a == b, "the rule must be deterministic"
    assert len({(c["rung"], c["family_tag"], c["seed"], c["margin"]) for c in a}) \
        == len(a), "cases must be distinct"
    assert all(c["rule_hash"] == L.RULE_HASH for c in a)


def test_rung0_preserves_the_v1_easy_cases():
    r0 = [c for c in L.generate() if c["rung"] == 0]
    assert r0, "rung 0 must exist"
    assert L.LADDER["rungs"][0]["name"] == "v1_easy_preserved"
    assert all(c["rung_name"] == "v1_easy_preserved" for c in r0)


def test_ladder_is_graded_in_candidate_count():
    cs = L.generate()
    by = {}
    for c in cs:
        by.setdefault(c["rung"], []).append(c["n_candidates"])
    peaks = [max(by[r]) for r in sorted(by)]
    assert peaks == sorted(peaks), f"candidate count must not decrease: {peaks}"
    assert peaks[-1] > peaks[0], "the ladder must actually get harder"


def test_every_case_builds_with_node_pair_edges():
    """CertifiedModel indexes node pairs; layered-DAG triples must be relabeled."""
    for c in L.generate(rungs=[0, 3], seeds_per_rung=1, margins=[0.1]):
        _, _, oracle, _ = L.build_case(c)
        assert all(len(e) == 2 for e in oracle.ground_set()), c["family_tag"]


# ------------------------------------------------------------ instrumentation
def test_meter_separates_cache_hits_from_misses_and_their_work():
    m = Meter()
    with m.cache_lookup() as c:
        c["hit"] = True
    with m.cache_lookup() as c:
        c["hit"] = False
        c["work"] = 7
    assert m.cache["hits"] == 1 and m.cache["misses"] == 1
    assert m.cache["work_on_hit"] == 0 and m.cache["work_on_miss"] == 7
    assert m.summary()["cache_hit_rate"] == pytest.approx(0.5)


def test_first_resolution_is_recorded_only_when_the_gap_closes():
    m = Meter(tolerance=1e-9)
    m.note_bound(-5.0, 5.0)
    assert m.first_resolution_s is None
    m.note_bound(1.0, 1.0)
    assert m.first_resolution_s is not None


def test_auc_charges_one_before_any_valid_bound():
    """An arm that has produced nothing must not be scored as if gap were zero."""
    m = Meter()
    auc, neg = m.anytime_auc(horizon=1.0, scale=1.0)
    assert auc == pytest.approx(1.0), "silence must cost the full integrand"
    assert neg == 0


def test_auc_records_negative_gaps_as_defects():
    m = Meter()
    m.note_bound(1.0, 0.5)          # inverted bound
    auc, neg = m.anytime_auc(horizon=0.1, scale=1.0)
    assert neg == 1, "a negative gap must be counted, not silently clipped"
    assert 0.0 <= auc <= 1.0


def test_instrumentation_does_not_change_the_bound():
    """A timed call must return exactly what an untimed one returns."""
    from scmp.bounds import certified_upper_bound
    c = L.generate(rungs=[1], seeds_per_rung=1, margins=[0.1])[0]
    model, X, oracle, B = L.build_case(c)
    plain = certified_upper_bound(model, X, oracle, B, mode="combined")
    m = Meter()
    with m.timed("bound_propagation"):
        timed = certified_upper_bound(model, X, oracle, B, mode="combined")
    assert timed.upper == plain.upper
    assert timed.incumbent_value == plain.incumbent_value


# ------------------------------------------------------- selective conditioning
def test_selective_skips_when_the_gap_is_already_closed():
    from scmp.revision.verifier_pilot import selective_decision
    c = L.generate(rungs=[0], seeds_per_rung=1, margins=[0.1])[0]
    model, X, oracle, B = L.build_case(c)
    sel = {"skip_if_gap_below": 1e18,            # force the skip branch
           "skip_if_stable_fraction": 1.1, "min_expected_gain": 0.0}
    action, _, _ = selective_decision(model, X, oracle, B, sel, Meter())
    assert action == "skip_already_resolved"


def test_selective_conditions_when_neither_skip_rule_fires():
    from scmp.revision.verifier_pilot import selective_decision
    c = L.generate(rungs=[3], seeds_per_rung=1, margins=[0.1])[0]
    model, X, oracle, B = L.build_case(c)
    sel = {"skip_if_gap_below": 0.0,             # never skip on gap
           "skip_if_stable_fraction": 1.1,       # never skip on stability
           "min_expected_gain": 0.0}
    action, _, _ = selective_decision(model, X, oracle, B, sel, Meter())
    assert action == "condition"


# --------------------------------------------------------- external references
def test_lp_polytope_bound_is_valid_against_exhaustive_optimum():
    """The LP relaxation must never fall below the true maximum."""
    from scmp.revision.verifier_pilot import exhaustive_optimum, lp_polytope_bound
    c = L.generate(rungs=[0], seeds_per_rung=1, margins=[0.1])[0]
    model, X, oracle, B = L.build_case(c)
    exact, info = exhaustive_optimum(model, X, oracle, B, 4000)
    ub, _ = lp_polytope_bound(model, X, oracle, B)
    if exact is None or ub is None:
        pytest.skip("family too large or LP unavailable")
    assert ub >= exact - 1e-6, f"LP bound {ub} below true optimum {exact}"


def test_exhaustive_optimum_is_attained_by_a_feasible_member():
    from scmp.revision.verifier_pilot import exhaustive_optimum
    c = L.generate(rungs=[0], seeds_per_rung=1, margins=[0.1])[0]
    model, X, oracle, B = L.build_case(c)
    exact, info = exhaustive_optimum(model, X, oracle, B, 4000)
    if exact is None:
        pytest.skip("family too large")
    z = [tuple(e) for e in info["argmax"]]
    assert oracle.feasibility(forced=tuple(z)).feasible
    with torch.no_grad():
        assert float(model(X, model.indicator(z), B)) == pytest.approx(exact)


# ---------------------------------------- ground truth must enumerate MEMBERS
@pytest.mark.parametrize("tag", ["layered_dag:3x3", "path_matching:6",
                                 "rna:GCGCAAAAGCGC"])
def test_enumeration_visits_exactly_the_family(tag):
    """Regression: 'feasible prefix' means EXTENDABLE, not 'is a member'.

    Layered-DAG paths are not downward-closed, so accepting extendable prefixes
    evaluated partial paths, inflated the maximum, and made every sound bound
    look unsound (66 false violations).
    """
    from scmp.revision.verifier_pilot import exhaustive_optimum
    c = [x for x in L.generate([0, 1], 1, [0.1]) if x["family_tag"] == tag][0]
    model, X, oracle, B = L.build_case(c)
    exact, info = exhaustive_optimum(model, X, oracle, B, 4000)
    assert exact is not None
    assert info["members_evaluated"] == oracle.count().value, (
        "enumeration must visit exactly the family the oracle counts")


@pytest.mark.parametrize("tag", ["layered_dag:3x3", "layered_dag:4x4",
                                 "path_matching:8", "rna:GCGCAAAAGCGC"])
def test_certified_bound_is_sound_against_exact_optimum(tag):
    from scmp.bounds import certified_upper_bound
    from scmp.revision.verifier_pilot import exhaustive_optimum
    c = [x for x in L.generate([0, 1, 2], 1, [0.1]) if x["family_tag"] == tag][0]
    model, X, oracle, B = L.build_case(c)
    exact, _ = exhaustive_optimum(model, X, oracle, B, 4000)
    if exact is None:
        pytest.skip("ground truth unavailable")
    r = certified_upper_bound(model, X, oracle, B, mode="combined")
    assert r.upper >= exact - 1e-9, f"upper {r.upper} below true optimum {exact}"
    assert r.incumbent_value <= exact + 1e-9, "incumbent above the true optimum"
