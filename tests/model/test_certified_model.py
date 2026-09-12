"""Tests for f(X, z) = b(X) + a(X)^T z + gamma * g(X, B + P z).

Exact-integer style assertions (bitwise equality, structural identity) are separated
from floating-point diagnostics, which carry an explicit tolerance.
"""
import itertools
import math

import pytest
import torch

from scmp.model import (AnchorEncoder, CertifiedModel, SignedResidualMPNN,
                        backbone_adjacency)
from scmp.oracles import RNANonCrossingOracle
from scmp.oracles.bruteforce import enumerate_rna

SEQS = ["GGAAACC", "GCGCAAAAGCGC", "GGGAAAUCC", "GGGGAAAACCCC"]


def onehot(seq):
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return X


def setup(seq, gamma=0.0, seed=0, use_residual=True, d_hid=8):
    o = RNANonCrossingOracle(seq, 3, True)
    cand = o.ground_set()
    m = CertifiedModel(4, cand, len(seq), d_hid=d_hid, gamma=gamma,
                       use_residual=use_residual, seed=seed)
    return o, m, onehot(seq)


# ------------------------------------------------------- gamma = 0 is an exact bypass
@pytest.mark.parametrize("seq", SEQS)
def test_gamma_zero_is_exactly_additive(seq):
    o, m, X = setup(seq, gamma=0.0, seed=1)
    assert not m.residual_active
    b, a = m.anchor_terms(X)
    for s in enumerate_rna(seq, 3, True):
        z = m.indicator(s)
        got = float(m(X, z))
        want = float(b) + sum(float(a[m.index[e]]) for e in sorted(s))
        assert got == want                       # bitwise: no residual contribution


@pytest.mark.parametrize("seq", SEQS)
def test_gamma_zero_matches_additive_oracle_inference(seq):
    """The oracle's argmax of a^T z must be the argmax of f itself when gamma = 0."""
    o, m, X = setup(seq, gamma=0.0, seed=2)
    coef = m.additive_coefficients(X)
    b = m.anchor_bias(X)
    members = enumerate_rna(seq, 3, True)
    r = o.support(coef)
    best_val = max(b + sum(coef[e] for e in s) for s in members)
    argmaxes = {tuple(sorted(s)) for s in members
                if b + sum(coef[e] for e in s) == best_val}
    assert r.witness in argmaxes                      # structural identity, exact
    assert float(m(X, m.indicator(r.witness))) == pytest.approx(best_val, rel=0, abs=1e-12)


@pytest.mark.parametrize("seq", ["GCGCAAAAGCGC"])
def test_gamma_zero_exact_with_integer_coefficients(seq):
    """With integer coefficients the agreement is exact integer arithmetic, not float."""
    o, m, X = setup(seq, gamma=0.0, seed=3)
    coef = {e: (7 * k % 11) - 5 for k, e in enumerate(m.candidates)}
    r = o.support(coef)
    assert r.is_exact
    best = max((sum(coef[e] for e in s), tuple(sorted(s)))
               for s in enumerate_rna(seq, 3, True))
    assert r.value == best[0]                         # exact integer equality


def test_residual_off_ablation_matches_gamma_zero():
    seq = "GCGCAAAAGCGC"
    _, m_off, X = setup(seq, gamma=1.0, seed=4, use_residual=False)
    _, m_zero, _ = setup(seq, gamma=0.0, seed=4, use_residual=True)
    assert m_off.g is None and not m_off.residual_active
    for s in enumerate_rna(seq, 3, True)[:12]:
        z = m_off.indicator(s)
        assert float(m_off(X, z)) == float(m_zero(X, z))


def test_gamma_switch_is_exact_not_small():
    """gamma = 0.0 must skip g entirely, so a NaN inside g cannot leak into f."""
    seq = "GCGCAAAAGCGC"
    _, m, X = setup(seq, gamma=0.0, seed=5)
    with torch.no_grad():
        m.g.lin0.weight.fill_(float("nan"))
    z = m.indicator(enumerate_rna(seq, 3, True)[3])
    assert math.isfinite(float(m(X, z)))
    m.set_gamma(1.0)
    assert math.isnan(float(m(X, z)))


# -------------------------------------------------------------- anchor independence
@pytest.mark.parametrize("seq", SEQS)
def test_anchor_does_not_depend_on_z(seq):
    o, m, X = setup(seq, gamma=1.0, seed=6)
    ref = m.additive_coefficients(X)
    for s in enumerate_rna(seq, 3, True)[:10]:
        _ = m(X, m.indicator(s))
        assert m.additive_coefficients(X) == ref


@pytest.mark.parametrize("seq", SEQS)
def test_anchor_does_not_depend_on_conditioning_mask(seq):
    """Coefficients handed to the oracle must be identical under every mask."""
    o, m, X = setup(seq, gamma=1.0, seed=7)
    ref = m.additive_coefficients(X)
    g = m.candidates
    for forced in [(), (g[0],)] + ([(g[0], g[1])] if len(g) > 1 else []):
        for forbidden in [(), (g[-1],)]:
            _ = o.count(forced, forbidden)
            assert m.additive_coefficients(X) == ref


def test_anchor_signature_excludes_adjacency():
    import inspect
    params = list(inspect.signature(AnchorEncoder.forward).parameters)
    assert params == ["self", "X", "candidates"]


# ------------------------------------------------------------- adjacency sensitivity
def test_residual_reads_the_adjacency():
    """g must be a function of A, not of X alone."""
    torch.manual_seed(0)
    X = torch.randn(6, 4, dtype=torch.float64).abs()
    g = SignedResidualMPNN(4, 8)
    A1 = torch.eye(6, dtype=torch.float64)
    A2 = torch.ones(6, 6, dtype=torch.float64)
    assert float(g(X, A1)) != float(g(X, A2))


@pytest.mark.parametrize("seq", ["GCGCAAAAGCGC", "GGGGAAAACCCC"])
def test_selected_edges_change_the_output_beyond_the_additive_part(seq):
    o, m, X = setup(seq, gamma=1.0, seed=8)
    members = enumerate_rna(seq, 3, True)
    seen = set()
    for s in members:
        z = m.indicator(s)
        with torch.no_grad():
            resid = float(m(X, z)) - float(m.linear_value(X, z))
        seen.add(round(resid, 12))
    assert len(seen) > 1, "residual took the same value on every structure"


def test_base_adjacency_changes_the_residual():
    seq = "GCGCAAAAGCGC"
    _, m, X = setup(seq, gamma=1.0, seed=9)
    z = m.indicator(enumerate_rna(seq, 3, True)[4])
    b1 = backbone_adjacency(len(seq), self_loops=True, chain=True)
    b2 = backbone_adjacency(len(seq), self_loops=True, chain=False)
    assert float(m(X, z, b1)) != float(m(X, z, b2))


def test_adjacency_is_symmetric_and_unnormalised():
    seq = "GCGCAAAAGCGC"
    _, m, X = setup(seq, gamma=1.0, seed=10)
    z = m.indicator(enumerate_rna(seq, 3, True)[5])
    base = backbone_adjacency(len(seq))
    A = m.adjacency(z, base)
    assert torch.equal(A, A.T)
    assert A.max() >= 1.0            # no row normalisation has been applied
    for e, k in m.index.items():
        if z[k] == 1.0:
            assert A[e[0], e[1]] >= 1.0


# ------------------------------------------------- relabelling with carried attributes
@pytest.mark.parametrize("seq", ["GGAAACC", "GCGCAAAAGCGC"])
def test_relabelling_with_carried_attributes(seq):
    """Permute the nodes, carry X and the candidate list, and f must not move."""
    n = len(seq)
    o = RNANonCrossingOracle(seq, 3, True)
    cand = o.ground_set()
    X = onehot(seq)
    m = CertifiedModel(4, cand, n, d_hid=8, gamma=1.0, seed=11)

    g = torch.Generator().manual_seed(4)
    perm = torch.randperm(n, generator=g).tolist()
    inv = [0] * n
    for new, old in enumerate(perm):
        inv[old] = new
    Xp = X[perm]                                  # attributes carried with the nodes
    cand_p = tuple(tuple(sorted((inv[i], inv[j]))) for (i, j) in cand)

    mp = CertifiedModel(4, cand_p, n, d_hid=8, gamma=1.0, seed=11)
    mp.load_state_dict(m.state_dict())            # identical parameters

    base = backbone_adjacency(n, chain=False)     # chain is not permutation invariant
    Pm = torch.zeros(n, n, dtype=torch.float64)
    for old, new in enumerate(inv):
        Pm[new, old] = 1.0
    base_p = Pm @ base @ Pm.T

    for s in enumerate_rna(seq, 3, True)[:12]:
        z = m.indicator(s)
        sp = [tuple(sorted((inv[i], inv[j]))) for (i, j) in s]
        zp = mp.indicator(sp)
        assert float(m(X, z, base)) == pytest.approx(float(mp(Xp, zp, base_p)), abs=1e-12)


# ------------------------------------------------------------------ tiny graph sweeps
@pytest.mark.parametrize("seq", ["GGAAACC", "GGGAAAUCC", "GCGCAAAAGCGC"])
@pytest.mark.parametrize("gamma", [0.0, 1.0, -0.5])
def test_all_tiny_graph_forwards(seq, gamma):
    """Every member of the family, both switch settings, must produce a finite value
    equal to the additive part plus gamma times the residual."""
    o, m, X = setup(seq, gamma=gamma, seed=12)
    base = backbone_adjacency(len(seq))
    for s in enumerate_rna(seq, 3, True):
        z = m.indicator(s)
        with torch.no_grad():
            f = float(m(X, z, base))
            lin = float(m.linear_value(X, z))
            assert math.isfinite(f)
            if gamma == 0.0:
                assert f == lin
            else:
                resid = float(m.g(X, m.adjacency(z, base)))
                assert f == pytest.approx(lin + gamma * resid, abs=1e-12)


def test_forward_on_the_empty_and_full_structures():
    seq = "GCGCAAAAGCGC"
    o, m, X = setup(seq, gamma=1.0, seed=13)
    empty = m.indicator(())
    assert math.isfinite(float(m(X, empty)))
    full_feasible = max(enumerate_rna(seq, 3, True), key=len)
    assert math.isfinite(float(m(X, m.indicator(full_feasible))))


def test_unknown_edge_is_rejected():
    seq = "GCGCAAAAGCGC"
    _, m, _ = setup(seq)
    with pytest.raises(KeyError):
        m.indicator([(0, 1)])


# --------------------------------------------------------------- parameter gradients
def test_parameter_gradients_match_finite_differences():
    """FLOAT DIAGNOSTIC. Central differences at h = 1e-6 in float64."""
    seq = "GCGCAAAAGCGC"
    o, m, X = setup(seq, gamma=0.7, seed=14)
    members = enumerate_rna(seq, 3, True)
    zs = [m.indicator(s) for s in members[:6]]
    base = backbone_adjacency(len(seq))

    def loss():
        return sum((m(X, z, base) - 0.3) ** 2 for z in zs)

    m.zero_grad()
    loss().backward()
    grads = {n: p.grad.clone() for n, p in m.named_parameters() if p.grad is not None}
    assert grads, "no parameter received a gradient"

    h = 1e-6
    checked = 0
    for name, p in m.named_parameters():
        flat = p.data.view(-1)
        for idx in range(0, flat.numel(), max(1, flat.numel() // 3)):
            orig = flat[idx].item()
            flat[idx] = orig + h
            with torch.no_grad():
                up = float(loss())
            flat[idx] = orig - h
            with torch.no_grad():
                dn = float(loss())
            flat[idx] = orig
            fd = (up - dn) / (2 * h)
            ad = float(grads[name].view(-1)[idx])
            assert fd == pytest.approx(ad, rel=1e-5, abs=1e-6), (name, idx, fd, ad)
            checked += 1
    assert checked >= 6


def test_anchor_and_residual_both_receive_gradients():
    seq = "GCGCAAAAGCGC"
    o, m, X = setup(seq, gamma=1.0, seed=15)
    z = m.indicator(enumerate_rna(seq, 3, True)[2])
    m.zero_grad()
    m(X, z).backward()
    anchor_grad = any(p.grad is not None and p.grad.abs().sum() > 0
                      for p in m.anchor.parameters())
    resid_grad = any(p.grad is not None and p.grad.abs().sum() > 0
                     for p in m.g.parameters())
    assert anchor_grad and resid_grad


def test_gamma_zero_gives_the_residual_no_gradient():
    seq = "GCGCAAAAGCGC"
    o, m, X = setup(seq, gamma=0.0, seed=16)
    z = m.indicator(enumerate_rna(seq, 3, True)[2])
    m.zero_grad()
    m(X, z).backward()
    assert all(p.grad is None or p.grad.abs().sum() == 0 for p in m.g.parameters())


# ------------------------------------------------- unsupported operations stay absent
def test_no_normalisation_attention_or_gating_in_the_certified_model():
    banned = (torch.nn.BatchNorm1d, torch.nn.LayerNorm, torch.nn.GroupNorm,
              torch.nn.MultiheadAttention, torch.nn.Softmax, torch.nn.Sigmoid,
              torch.nn.Dropout)
    _, m, _ = setup("GCGCAAAAGCGC", gamma=1.0)
    for mod in m.modules():
        assert not isinstance(mod, banned), f"unsupported layer present: {type(mod)}"


def test_residual_uses_only_relu_and_linear():
    _, m, _ = setup("GCGCAAAAGCGC", gamma=1.0)
    kinds = {type(x) for x in m.g.modules()} - {SignedResidualMPNN}
    assert kinds <= {torch.nn.Linear, torch.nn.ReLU}


def test_weights_are_signed_not_constrained():
    """The residual is allowed negative weights; that is the point of the ablation."""
    _, m, _ = setup("GCGCAAAAGCGC", gamma=1.0, seed=17)
    w = torch.cat([p.view(-1) for p in m.g.parameters()])
    assert (w < 0).any() and (w > 0).any()
