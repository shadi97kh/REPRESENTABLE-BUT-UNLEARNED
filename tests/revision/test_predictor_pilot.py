"""P1 mechanism tests: the wiring the repair depends on.

These do not train to convergence. They pin the properties the A-initialised
residual relies on, so a refactor cannot silently reintroduce the defects P1
found in the historical pilot.
"""
import numpy as np
import pytest
import torch

from scmp.model import CertifiedModel, backbone_adjacency
from scmp.revision.predictor_pilot import (PlainMPNN, PlainMPNNSkip, _size_holdout,
                                           _state_digest, _zero_init_residual,
                                           feasible_mixed_differences, make_family,
                                           make_target, projection_residual_1z,
                                           stage_preserve, verify_preserved)


# ----------------------------------------------- the defect the pilot actually had
def test_empty_candidate_set_makes_gamma0_graph_blind():
    """certifiable_gamma0 with no candidates is b(X): it cannot see the adjacency."""
    n, d = 6, 4
    X = torch.rand(n, d, dtype=torch.float64,
                   generator=torch.Generator().manual_seed(0))
    cm = CertifiedModel(d, tuple(), n, d_hid=8, gamma=0.0, seed=0)
    z = torch.zeros(0, dtype=torch.float64)
    with torch.no_grad():
        a = float(cm(X, z, torch.eye(n, dtype=torch.float64)))
        b = float(cm(X, z, torch.ones(n, n, dtype=torch.float64)))
    assert a == b, "gamma0 must be graph-blind when there are no candidates"


def test_empty_candidate_set_means_a_dot_z_is_identically_zero():
    n, d = 6, 4
    cm = CertifiedModel(d, tuple(), n, d_hid=8, gamma=1.0, seed=0)
    X = torch.rand(n, d, dtype=torch.float64)
    assert len(cm.candidates) == 0
    assert cm.additive_coefficients(X) == {}


def test_skip_arm_restores_input_parity():
    """The historical baseline has no non-adjacency path; the parity arm does."""
    n, d, h = 6, 4, 5
    X = torch.rand(2, n, d, dtype=torch.float64)
    A = torch.zeros(2, n, n, dtype=torch.float64)      # adjacency with NO edges
    torch.manual_seed(0)
    plain, skip = PlainMPNN(d, h), PlainMPNNSkip(d, h)
    X2 = X + 3.0                       # a materially different input
    with torch.no_grad():
        p1, p2 = plain(X, A), plain(X2, A)
        s1, s2 = skip(X, A), skip(X2, A)
    # with no edges the plain MPNN emits only its readout bias: X cannot reach it
    assert torch.allclose(p1, p2), \
        "plain MPNN has no route from X to the output that bypasses the adjacency"
    assert not torch.allclose(s1, s2), \
        "the parity arm must still route X to the output when the graph is empty"


# ------------------------------------------------- A-initialised residual mechanics
def test_zero_init_leaves_output_zero_but_gradient_alive():
    n, d = 6, 3
    edges = ((0, 3), (1, 4))
    cm = CertifiedModel(d, edges, n, d_hid=6, gamma=1.0, seed=0)
    info = _zero_init_residual(cm)
    assert info["readout_zeroed"] and info["hidden_nonzero"]
    assert not info["gate_zeroed"], "zeroing the gate too would kill the gradient"

    X = torch.rand(n, d, dtype=torch.float64)
    B = backbone_adjacency(n)
    z = torch.tensor([1.0, 0.0], dtype=torch.float64)
    with torch.no_grad():
        anchor = float(cm.linear_value(X, z))
        full = float(cm(X, z, B))
    assert full == pytest.approx(anchor, abs=1e-12), \
        "at init the residual contributes exactly zero"

    cm(X, z, B).backward()
    assert cm.g.readout.weight.grad.abs().sum() > 0, "residual gradient must flow"


def test_gamma0_bypass_is_exact_under_a_poisoned_residual():
    n, d = 6, 3
    cm = CertifiedModel(d, ((0, 3),), n, d_hid=4, gamma=0.0, seed=0)
    with torch.no_grad():
        cm.g.readout.weight.fill_(float("nan"))
    X = torch.rand(n, d, dtype=torch.float64)
    z = torch.tensor([1.0], dtype=torch.float64)
    with torch.no_grad():
        v = float(cm(X, z, backbone_adjacency(n)))
    assert np.isfinite(v), "gamma=0 must never evaluate g"


# ------------------------------------------------------- interaction diagnostics
def test_projection_detects_additive_and_interacting_targets():
    n, edges, members = make_family({"name": "layered_dag", "layers": 3, "width": 3})
    Z = np.zeros((len(members), len(edges)))
    idx = {e: k for k, e in enumerate(edges)}
    for r, m in enumerate(members):
        for e in m:
            Z[r, idx[e]] = 1.0
    y_add, _ = make_target("additive", edges, members, seed=1)
    y_int, _ = make_target("pairwise", edges, members, seed=1)
    assert projection_residual_1z(Z, y_add)["additive_in_z"] is True
    assert projection_residual_1z(Z, y_int)["additive_in_z"] is False


def test_mixed_differences_vanish_only_for_additive_targets():
    n, edges, members = make_family({"name": "path_matching", "n_nodes": 6})
    look = {frozenset(m): i for i, m in enumerate(members)}
    for kind, expect_zero in (("additive", True), ("pairwise", False)):
        y, _ = make_target(kind, edges, members, seed=2)
        md = feasible_mixed_differences(members, edges,
                                        lambda S: float(y[look[S]]), limit=80)
        assert md, "the family must contain feasible rectangles"
        mx = max(abs(v) for v in md)
        assert (mx < 1e-9) is expect_zero, (kind, mx)


def test_target_generator_is_not_our_architecture():
    _, edges, members = make_family({"name": "path_matching", "n_nodes": 6})
    _, meta = make_target("pairwise", edges, members, seed=0)
    assert "independent of our MPNN" in meta["generator"]
    assert meta["n_pairs"] > 0


# -------------------------------------------------------------------- holdout
def test_size_holdout_separates_structure_sizes():
    _, edges, members = make_family({"name": "path_matching", "n_nodes": 8})
    tr, te, info = _size_holdout(members)
    assert info["kind"] == "size"
    assert not (set(tr) & set(te))
    sizes = np.array([len(m) for m in members])
    assert sizes[tr].max() < sizes[te].min(), "test must hold out LARGER structures"


# -------------------------------------------------------------------- preservation
def test_preservation_fingerprints_are_stable():
    a = stage_preserve(None)
    chk = verify_preserved(a)
    assert chk["preserved"] and not chk["drift"]


def test_state_digest_changes_when_a_parameter_changes():
    cm = CertifiedModel(3, ((0, 2),), 4, d_hid=4, gamma=0.0, seed=0)
    d0 = _state_digest(cm)
    with torch.no_grad():
        cm.anchor.a_out.bias.add_(1.0)
    assert _state_digest(cm) != d0


# ------------------------------------------------- batched forward must be exact
@pytest.mark.parametrize("gamma", [0.0, 1.0])
def test_batched_forward_is_bit_identical_to_the_loop(gamma):
    """The speedup that makes the family experiments affordable must change nothing."""
    from scmp.revision.predictor_pilot import (_fwd_batched, _fwd_members, _Zmat,
                                               _member_adjacencies)
    n, edges, members = make_family({"name": "path_matching", "n_nodes": 6})
    members = members[:40]
    dev = torch.device("cpu")
    X = torch.rand(n, 4, dtype=torch.float64,
                   generator=torch.Generator().manual_seed(0))
    B = backbone_adjacency(n)
    cm = CertifiedModel(4, edges, n, d_hid=5, gamma=gamma, seed=0)
    As = _member_adjacencies(cm, X, B, members, dev)
    Zt = torch.tensor(_Zmat(members, edges), dtype=torch.float64)
    with torch.no_grad():
        a = _fwd_batched(cm, X, B, Zt, As).numpy()
        b = _fwd_members(cm, X, B, members, dev).numpy()
    assert np.array_equal(a, b), "batched forward must be bit-identical"


def test_interaction_pairs_can_actually_co_occur():
    """Regression: keying interactions on vertex-sharing pairs made the
    'interacting' target silently additive on a matching family."""
    _, edges, members = make_family({"name": "path_matching", "n_nodes": 6})
    y, meta = make_target("pairwise", edges, members, seed=3)
    Z = np.zeros((len(members), len(edges)))
    idx = {e: k for k, e in enumerate(edges)}
    for r, m in enumerate(members):
        for e in m:
            Z[r, idx[e]] = 1.0
    assert projection_residual_1z(Z, y)["additive_in_z"] is False
    assert meta["n_pairs"] > 0
