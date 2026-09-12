"""Soundness of the affine bound propagation.

Every intermediate bracket is checked against exhaustive enumeration of the feasible
family, not only the final bound. A single violation anywhere is a hard failure: an
unsound intermediate invalidates everything downstream even when the final number
happens to look right.

All numbers here are float64 with no directed rounding, so they are DIAGNOSTIC. The
tests assert the real-arithmetic inequalities up to a stated slack; see
docs/bounds_proof.md for what would have to change to make them certificates.
"""
import numpy as np
import pytest
import torch

from scmp.bounds import (FLOAT_UNVERIFIED, INFEASIBLE, PROVED, Aff, ResidualBounder,
                         certified_upper_bound, mccormick, relu_bounds)
from scmp.model import CertifiedModel, backbone_adjacency
from scmp.oracles import RNANonCrossingOracle
from scmp.oracles.bruteforce import enumerate_rna

SLACK = 1e-9
SEQS = ["GGAAACC", "GGGAAAUCC", "GCGCAAAAGCGC"]


def onehot(seq):
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return X


def build(seq, gamma=1.0, seed=0, d_hid=5, scale=0.0, pseed=0):
    o = RNANonCrossingOracle(seq, 3, True)
    m = CertifiedModel(4, o.ground_set(), len(seq), d_hid=d_hid, gamma=gamma, seed=seed)
    if scale:
        g = torch.Generator().manual_seed(pseed)
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype) * scale)
    return o, m, onehot(seq), backbone_adjacency(len(seq))


def zvec(m, s):
    z = np.zeros(len(m.candidates))
    for e in s:
        z[m.index[e]] = 1.0
    return z


# ------------------------------------------------------------- envelope, in isolation
def test_all_four_mccormick_inequalities_are_formed_and_valid():
    """Each of the four is checked directly on a grid of (e, h)."""
    m = 1
    h_lo = Aff(np.array([-0.5]), 0.2)          # h = 0.2 - 0.5 z
    h_hi = Aff(np.array([-0.5]), 0.9)          # h = 0.9 - 0.5 z
    lowers, uppers = mccormick(0, h_lo, h_hi, m)
    assert len(lowers) == 2 and len(uppers) == 2
    for e in np.linspace(0.0, 1.0, 21):
        z = np.array([e])
        for h in np.linspace(h_lo.eval(z), h_hi.eval(z), 11):
            y = e * h
            for L in lowers:
                assert L.eval(z) <= y + SLACK
            for U in uppers:
                assert U.eval(z) >= y - SLACK


def test_relu_stable_cases_are_exact():
    m = 2
    pos_lo, pos_hi = Aff(np.array([1.0, 0.0]), 2.0), Aff(np.array([1.0, 0.0]), 3.0)
    lo, hi, kind = relu_bounds(pos_lo, pos_hi)
    assert kind == "stable_pos" and lo is pos_lo and hi is pos_hi
    neg_lo, neg_hi = Aff(np.array([0.0, 0.0]), -3.0), Aff(np.array([0.0, 0.0]), -1.0)
    lo, hi, kind = relu_bounds(neg_lo, neg_hi)
    assert kind == "stable_neg"
    assert lo.const == 0.0 and hi.const == 0.0 and not lo.coef.any()


def test_relu_unstable_brackets_the_true_value():
    lo = Aff(np.array([2.0]), -1.0)
    hi = Aff(np.array([2.0]), -1.0)
    low, up, kind = relu_bounds(lo, hi)
    assert kind == "unstable"
    for e in np.linspace(0, 1, 21):
        z = np.array([e])
        true = max(0.0, lo.eval(z))
        assert low.eval(z) <= true + SLACK <= up.eval(z) + 2 * SLACK


# --------------------------------------------- every intermediate, every feasible graph
@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("scale,pseed", [(0.0, 0), (0.5, 1), (1.0, 2)])
def test_every_intermediate_bracket_holds_on_every_feasible_graph(seq, scale, pseed):
    o, m, X, B = build(seq, gamma=1.0, seed=0, scale=scale, pseed=pseed)
    rb = ResidualBounder(m, X, B)
    lo_aff, hi_aff, tr = rb.propagate()
    members = enumerate_rna(seq, 3, True)
    checked = 0
    for s in members:
        z = zvec(m, s)
        truth = rb.reference_forward(z)
        for name, (L, H) in tr.layers.items():
            arr = truth[name]
            for i in range(len(L)):
                for k in range(len(L[i])):
                    lv, hv, tv = L[i][k].eval(z), H[i][k].eval(z), float(arr[i][k])
                    assert lv <= tv + SLACK, (name, i, k, lv, tv)
                    assert tv <= hv + SLACK, (name, i, k, tv, hv)
                    checked += 1
    assert checked > 0


@pytest.mark.parametrize("seq", SEQS)
def test_reference_forward_agrees_with_the_torch_model(seq):
    """The numpy reference must reproduce g exactly, or the intermediate checks above
    would be validating the wrong function."""
    o, m, X, B = build(seq, gamma=1.0, seed=3, scale=0.5, pseed=7)
    rb = ResidualBounder(m, X, B)
    for s in enumerate_rna(seq, 3, True)[:15]:
        z = zvec(m, s)
        with torch.no_grad():
            want = float(m.g(X, m.adjacency(torch.tensor(z, dtype=torch.float64), B)))
        assert rb.reference_forward(z)["out"][0][0] == pytest.approx(want, abs=1e-10)


# ------------------------------------------------------------------- the final bound
@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("gamma", [1.0, -1.0, 0.5])
@pytest.mark.parametrize("scale,pseed", [(0.0, 0), (0.8, 5)])
def test_upper_bound_dominates_the_true_maximum(seq, gamma, scale, pseed):
    o, m, X, B = build(seq, gamma=gamma, seed=1, scale=scale, pseed=pseed)
    members = enumerate_rna(seq, 3, True)
    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)
    for mode in ("combined", "terminal_only"):
        r = certified_upper_bound(m, X, o, B, mode=mode)
        assert r.upper >= true_max - SLACK, (mode, r.upper, true_max)


@pytest.mark.parametrize("seq", SEQS)
def test_combined_is_never_looser_than_separate_maxima(seq):
    """U = b + gamma d + support(a + gamma c) must dominate the terminal-only baseline
    that maximises the two pieces separately."""
    for scale, pseed in [(0.0, 0), (0.5, 1), (1.0, 2)]:
        o, m, X, B = build(seq, gamma=1.0, seed=2, scale=scale, pseed=pseed)
        rc = certified_upper_bound(m, X, o, B, mode="combined")
        rt = certified_upper_bound(m, X, o, B, mode="terminal_only")
        assert rc.upper <= rt.upper + SLACK


def test_combined_is_strictly_tighter_on_a_non_degenerate_instance():
    """When every coefficient is negative the family optimum is the empty structure and
    both modes coincide; the separation must be exhibited where signs are mixed."""
    seq = "GGGAAAUCC"
    strict = 0
    for pseed in range(6):
        o, m, X, B = build(seq, gamma=1.0, seed=pseed, scale=1.0, pseed=100 + pseed)
        rc = certified_upper_bound(m, X, o, B, mode="combined")
        rt = certified_upper_bound(m, X, o, B, mode="terminal_only")
        strict += rc.upper < rt.upper - 1e-9
    assert strict > 0


def test_baselines_share_identical_predictor_weights():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=1.0, seed=4, scale=0.7, pseed=9)
    rc = certified_upper_bound(m, X, o, B, mode="combined")
    rt = certified_upper_bound(m, X, o, B, mode="terminal_only")
    assert rc.model_hash == rt.model_hash


# ------------------------------------------------------------------------- incumbents
@pytest.mark.parametrize("seq", SEQS)
def test_incumbent_is_a_real_model_evaluation_not_an_envelope_value(seq):
    o, m, X, B = build(seq, gamma=1.0, seed=5, scale=0.6, pseed=3)
    r = certified_upper_bound(m, X, o, B)
    with torch.no_grad():
        direct = float(m(X, m.indicator(r.incumbent_witness), B))
    assert r.incumbent_value == pytest.approx(direct, abs=0.0, rel=0.0)
    assert r.incumbent_value <= r.upper + SLACK
    assert r.gap >= -SLACK
    members = {tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)}
    assert tuple(sorted(r.incumbent_witness)) in members       # feasible, not relaxed


# ------------------------------------------------------------------- gamma = 0 is exact
@pytest.mark.parametrize("seq", SEQS)
def test_gamma_zero_bound_is_exact(seq):
    """With no residual the model is additive and the oracle is exact, so the upper
    bound must equal the true maximum and the gap must close."""
    o, m, X, B = build(seq, gamma=0.0, seed=6, scale=0.5, pseed=4)
    members = enumerate_rna(seq, 3, True)
    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)
    r = certified_upper_bound(m, X, o, B)
    assert r.upper == pytest.approx(true_max, abs=1e-12)
    assert abs(r.gap) < 1e-12


# ----------------------------------------------------------------------- edge cases
def test_negative_anchor_coefficients_select_the_empty_structure():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=0.0, seed=0)
    coefs = m.additive_coefficients(X)
    assert all(c < 0 for c in coefs.values())
    r = certified_upper_bound(m, X, o, B)
    assert r.incumbent_witness == ()


def test_forced_crossing_edges_yield_infeasible_not_a_bound():
    seq = "GCGCAAAAGCGC"
    o, m, X, B = build(seq, gamma=1.0, seed=0, scale=0.5, pseed=1)
    g = o.ground_set()
    cross = [(a, b) for a in g for b in g if a[0] < b[0] < a[1] < b[1]]
    assert cross
    r = certified_upper_bound(m, X, o, B, forced=cross[0])
    assert r.termination_status == INFEASIBLE
    assert r.upper == float("-inf")
    assert r.incumbent_witness == ()


def test_incompatible_forced_and_forbidden_yield_infeasible():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=1.0, seed=0)
    e = o.ground_set()[0]
    r = certified_upper_bound(m, X, o, B, forced=(e,), forbidden=(e,))
    assert r.termination_status == INFEASIBLE


def test_empty_masks_match_the_unconditioned_bound():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=1.0, seed=0, scale=0.5, pseed=2)
    a = certified_upper_bound(m, X, o, B)
    b = certified_upper_bound(m, X, o, B, forced=(), forbidden=())
    assert a.upper == b.upper and a.incumbent_witness == b.incumbent_witness


def test_zero_weights_give_a_constant_residual_and_a_sound_bound():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=1.0, seed=0)
    with torch.no_grad():
        for p in m.g.parameters():
            p.zero_()
    members = enumerate_rna(seq, 3, True)
    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)
    r = certified_upper_bound(m, X, o, B)
    assert r.upper >= true_max - SLACK
    lo, hi, tr = ResidualBounder(m, X, B).propagate()
    assert not np.any(hi.coef) and hi.const == pytest.approx(0.0, abs=1e-12)


def test_stable_relus_are_detected_in_both_directions():
    seq = "GCGCAAAAGCGC"
    kinds = set()
    for pseed in range(8):
        o, m, X, B = build(seq, gamma=1.0, seed=0, scale=1.2, pseed=pseed)
        _, _, tr = ResidualBounder(m, X, B).propagate()
        for rows in tr.relu_kinds.values():
            for row in rows:
                kinds.update(row)
    assert {"stable_pos", "stable_neg", "unstable"} <= kinds


def test_tied_maxima_return_one_valid_argmax_and_a_sound_bound():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=0.0, seed=0)
    ties = {e: 1 for e in m.candidates}          # every edge worth the same
    r = o.support(ties)
    best = max(len(s) for s in enumerate_rna(seq, 3, True))
    assert r.value == best                        # exact integer
    argmaxes = {tuple(sorted(s)) for s in enumerate_rna(seq, 3, True) if len(s) == best}
    assert r.witness in argmaxes


def test_envelope_terms_are_counted_and_edge_variables_are_shared():
    """Four inequalities per incident-edge product term, and one variable per candidate
    edge in the final affine form -- not one per (edge, layer)."""
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=1.0, seed=0, scale=0.5, pseed=1)
    lo, hi, tr = ResidualBounder(m, X, B).propagate()
    assert tr.envelope_terms > 0 and tr.envelope_terms % 4 == 0
    assert len(hi.coef) == len(m.candidates)
    assert len(lo.coef) == len(m.candidates)


# ------------------------------------------------------------------ numerical honesty
def test_bounds_are_labelled_diagnostic_until_roundoff_is_bounded():
    seq = "GGGAAAUCC"
    o, m, X, B = build(seq, gamma=1.0, seed=0, scale=0.5, pseed=1)
    r = certified_upper_bound(m, X, o, B)
    assert r.numerical_status == FLOAT_UNVERIFIED
    assert any("float64" in a for a in r.assumptions)
    assert any("DIAGNOSTIC" in a for a in r.assumptions)
    assert r.termination_status == PROVED
    assert len(r.model_hash) == 16 and len(r.family_hash) == 16
