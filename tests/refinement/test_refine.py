"""Branch refinement: coverage, monotonicity, pruning discipline, and exactness.

Every property is checked against exhaustive enumeration of the feasible family, and
every certificate is re-checked by `scmp.check_certificates`, which does not import the
search it is checking.
"""
import json

import pytest
import torch

from scmp.check_certificates import Checker
from scmp.model import CertifiedModel, backbone_adjacency
from scmp.numerics import FLOAT_UNVERIFIED
from scmp.oracles import RNANonCrossingOracle
from scmp.oracles.bruteforce import enumerate_rna
from scmp.refine import INFEASIBLE, PROVED, UNRESOLVED, refine

TOL = 1e-9
SEQS = ["GGAAACC", "GGGAAAUCC", "GCGCAAAAGCGC"]


def build(seq, gamma=1.0, seed=0, d_hid=4, scale=1.0, pseed=0):
    o = RNANonCrossingOracle(seq, 3, True)
    m = CertifiedModel(4, o.ground_set(), len(seq), d_hid=d_hid, gamma=gamma, seed=seed)
    if scale:
        g = torch.Generator().manual_seed(pseed)
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype) * scale)
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return o, m, X, backbone_adjacency(len(seq))


def find_instance(seq, want_elims=False, max_pseed=16):
    """Locate a perturbation that actually exercises splitting (and optionally
    elimination). Many instances are proved at the root with zero splits -- correct
    behaviour, but useless for testing the branching machinery."""
    for pseed in range(max_pseed):
        o, m, X, B = build(seq, pseed=pseed)
        r = refine(m, X, o, B, max_nodes=200)
        splits = sum(1 for x in r.records if x["op"] == "split")
        elims = sum(1 for x in r.records if x["op"] == "eliminate")
        if splits and (elims or not want_elims):
            return pseed, o, m, X, B, r
    pytest.skip(f"no splitting instance found for {seq}")


def true_max(m, X, B, members):
    with torch.no_grad():
        return max(float(m(X, m.indicator(s), B)) for s in members)


# ------------------------------------------------------------------ core guarantees
@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("pseed", [0, 1, 2])
def test_global_bound_is_monotone_non_increasing(seq, pseed):
    o, m, X, B = build(seq, pseed=pseed)
    r = refine(m, X, o, B, max_nodes=60)
    hist = r.stats["history"]
    for a, b in zip(hist, hist[1:]):
        assert b <= a + TOL, (a, b)


@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("pseed", [0, 1, 2])
def test_global_bound_dominates_and_incumbent_is_attained(seq, pseed):
    o, m, X, B = build(seq, pseed=pseed)
    members = enumerate_rna(seq, 3, True)
    tm = true_max(m, X, B, members)
    r = refine(m, X, o, B, max_nodes=60)
    assert r.global_upper >= tm - TOL
    assert r.incumbent <= tm + TOL
    with torch.no_grad():
        assert float(m(X, m.indicator(r.incumbent_witness), B)) == r.incumbent


@pytest.mark.parametrize("seq", SEQS)
def test_eventual_exactness_on_tiny_families(seq):
    """With enough nodes the search must close the gap exactly."""
    for pseed in range(3):
        o, m, X, B = build(seq, pseed=pseed)
        members = enumerate_rna(seq, 3, True)
        tm = true_max(m, X, B, members)
        r = refine(m, X, o, B, max_nodes=500, time_limit=30.0)
        assert r.status == PROVED, r.stats
        assert r.global_upper == pytest.approx(tm, abs=1e-9)
        assert r.gap == pytest.approx(0.0, abs=1e-9)


def test_timeout_returns_unresolved_not_a_claim():
    """A node budget too small to close the tree must report unresolved, never a
    proved status with an unverified bound."""
    _, o, m, X, B, _ = find_instance("GGGAAAUCC")
    r = refine(m, X, o, B, max_nodes=1)
    assert r.status == UNRESOLVED
    members = enumerate_rna("GGGAAAUCC", 3, True)
    assert r.global_upper >= true_max(m, X, B, members) - TOL
    assert r.gap >= -TOL


def test_proving_at_the_root_is_a_valid_outcome():
    """When the root bound already equals an attained incumbent the search is finished
    with zero splits; that is optimality, not a missing expansion."""
    o, m, X, B = build("GCGCAAAAGCGC", pseed=8)
    members = enumerate_rna("GCGCAAAAGCGC", 3, True)
    r = refine(m, X, o, B, max_nodes=200)
    if r.stats["expanded"] == 0:
        assert r.status == PROVED
        assert r.global_upper == pytest.approx(true_max(m, X, B, members), abs=1e-9)


def test_zero_time_limit_returns_unresolved():
    o, m, X, B = build("GGGAAAUCC", pseed=4)
    r = refine(m, X, o, B, time_limit=0.0)
    assert r.status == UNRESOLVED


def test_infeasible_root_is_reported_not_bounded():
    o, m, X, B = build("GGGAAAUCC", pseed=5)
    e = o.ground_set()[0]
    r = refine(m, X, o, B, forced=(e,), forbidden=(e,))
    assert r.status == INFEASIBLE
    assert r.global_upper == float("-inf")


# ------------------------------------------------------------------------ structure
@pytest.mark.parametrize("seq", ["GGGAAAUCC", "GCGCAAAAGCGC"])
def test_splits_cover_the_parent(seq):
    _, o, m, X, B, r = find_instance(seq)
    members = enumerate_rna(seq, 3, True)
    rep = Checker(r.certificate(), m, X, o, B, members, seq).run()
    assert rep["accepted"], rep["failures"]
    assert any(x["op"] == "split" for x in r.records)


def test_child_cap_is_min_of_parent_and_new():
    _, o, m, X, B, r = find_instance("GCGCAAAAGCGC")
    caps = {}
    for rec in r.records:
        if rec["op"] != "bound":
            continue
        if rec.get("parent") is None:
            caps[rec["node"]] = rec["cap"]
            continue
        assert rec["cap"] == min(rec["parent_cap"], rec["raw"])
        assert rec["cap"] <= caps[rec["parent"]] + TOL
        caps[rec["node"]] = rec["cap"]


def test_empty_children_are_eliminated_with_a_proof():
    """Deeper in the tree a forced edge can make one side of a split empty; the
    elimination must carry the oracle's reason and be re-provable."""
    seq = "GGGAAAUCC"
    _, o, m, X, B, r = find_instance(seq, want_elims=True)
    elims = [x for x in r.records if x["op"] == "eliminate"]
    assert elims, "no elimination exercised"
    for e in elims:
        assert e["reason"]
    members = enumerate_rna(seq, 3, True)
    assert Checker(r.certificate(), m, X, o, B, members).run()["accepted"]


def test_pruning_only_when_the_bound_cannot_beat_the_incumbent():
    _, o, m, X, B, r = find_instance("GCGCAAAAGCGC")
    for rec in [x for x in r.records if x["op"] == "prune"]:
        assert rec["cap"] <= rec["incumbent"] + TOL


def test_no_feasible_branch_is_dropped_by_ranking():
    """A split policy that ranks adversarially must change the order, not the result."""
    seq = "GGGAAAUCC"
    members = enumerate_rna(seq, 3, True)
    o, m, X, B = build(seq, pseed=10)
    tm = true_max(m, X, B, members)
    a = refine(m, X, o, B, max_nodes=500)
    b = refine(m, X, o, B, max_nodes=500,
               split_policy=lambda und, node: und[-1])
    assert a.status == b.status == PROVED
    assert a.global_upper == pytest.approx(tm, abs=1e-9)
    assert b.global_upper == pytest.approx(tm, abs=1e-9)


def test_exact_singleton_evaluation_in_reference_mode():
    seq = "GGGAAAUCC"
    _, o, m, X, B, _ = find_instance(seq)
    r = refine(m, X, o, B, max_nodes=500, reference_mode=True)
    singles = [x for x in r.records if x["op"] == "singleton"]
    assert singles, "no singleton node reached"
    for s in singles:
        assert s["count"] == 1
        w = tuple(tuple(e) for e in s["witness"])
        with torch.no_grad():
            assert float(m(X, m.indicator(w), B)) == s["value"]


def test_fairness_eventually_serves_the_queue_fifo():
    _, o, m, X, B, _ = find_instance("GCGCAAAAGCGC")
    r = refine(m, X, o, B, max_nodes=60, fairness_period=2)
    assert r.stats["fairness_pops"] > 0
    r2 = refine(m, X, o, B, max_nodes=60, fairness_period=0)
    assert r2.stats["fairness_pops"] == 0


# ------------------------------------------------------------------- certificates
@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("nodes", [2, 500])
def test_certificate_is_accepted_by_the_independent_checker(seq, nodes):
    o, m, X, B = build(seq, pseed=13)
    members = enumerate_rna(seq, 3, True)
    r = refine(m, X, o, B, max_nodes=nodes)
    rep = Checker(r.certificate(), m, X, o, B, members, seq).run()
    assert rep["accepted"], rep["failures"]
    assert rep["checks"] > 0


def test_certificate_round_trips_through_json():
    o, m, X, B = build("GGGAAAUCC", pseed=14)
    r = refine(m, X, o, B, max_nodes=40)
    cert = json.loads(json.dumps(r.certificate()))
    members = enumerate_rna("GGGAAAUCC", 3, True)
    assert Checker(cert, m, X, o, B, members).run()["accepted"]


def test_stale_model_hash_is_rejected():
    o, m, X, B = build("GGGAAAUCC", pseed=15)
    r = refine(m, X, o, B, max_nodes=40)
    cert = r.certificate()
    with torch.no_grad():
        m.g.lin0.weight.add_(1.0)                      # certificate now describes a
    rep = Checker(cert, m, X, o, B).run()              # different model
    assert not rep["accepted"]
    assert any(f["kind"] == "stale_model_hash" for f in rep["failures"])


def test_stale_family_hash_is_rejected():
    o, m, X, B = build("GGGAAAUCC", pseed=16)
    r = refine(m, X, o, B, max_nodes=40)
    cert = r.certificate()
    other = RNANonCrossingOracle("GGGAAAUCC", 1, True)   # different min_loop
    rep = Checker(cert, m, X, other, B).run()
    assert not rep["accepted"]
    assert any(f["kind"] == "stale_family_hash" for f in rep["failures"])


def test_checker_catches_a_tampered_cap():
    _, o, m, X, B, r = find_instance("GCGCAAAAGCGC")
    cert = json.loads(json.dumps(r.certificate()))
    for rec in cert["records"]:
        if rec["op"] == "bound" and rec.get("parent") is not None:
            rec["cap"] = rec["cap"] - 1000.0            # claim a tighter cap
            break
    rep = Checker(cert, m, X, o, B).run()
    assert not rep["accepted"]


def test_checker_catches_a_tampered_incumbent():
    _, o, m, X, B, r = find_instance("GGGAAAUCC")
    cert = json.loads(json.dumps(r.certificate()))
    for rec in cert["records"]:
        if rec["op"] == "incumbent":
            rec["value"] = rec["value"] + 5.0
            break
    rep = Checker(cert, m, X, o, B).run()
    assert not rep["accepted"]
    assert any("incumbent" in f["kind"] for f in rep["failures"])


def test_checker_catches_an_eliminated_feasible_branch():
    _, o, m, X, B, r = find_instance("GCGCAAAAGCGC")
    cert = json.loads(json.dumps(r.certificate()))
    split = next(x for x in cert["records"] if x["op"] == "split")
    cert["records"].append({"op": "eliminate", "node": 9999, "parent": split["node"],
                            "side": "forced", "edge": split["edge"],
                            "reason": "fabricated"})
    rep = Checker(cert, m, X, o, B).run()
    assert not rep["accepted"]


# ---------------------------------------------------------------- numerical honesty
def test_numerical_status_is_inherited_not_upgraded():
    """Refinement consumes float64 bounds, so it cannot report better than they do."""
    o, m, X, B = build("GGGAAAUCC", pseed=20)
    r = refine(m, X, o, B, max_nodes=40)
    assert r.numerical_status == FLOAT_UNVERIFIED
    assert r.certificate()["numerical_status"] == FLOAT_UNVERIFIED
