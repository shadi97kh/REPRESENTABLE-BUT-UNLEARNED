"""RNA ensemble machinery: covariates, windows, priors, expectations, brackets.

Ground truth is exhaustive enumeration of the family wherever a mean is claimed.
"""
import math, random
import numpy as np
import pytest
import torch

from scmp.model import CertifiedModel, backbone_adjacency
from scmp.numerics import Interval, idot
from scmp.oracles import RNANonCrossingOracle
from scmp.oracles.bruteforce import enumerate_rna, rna_is_valid
from scmp.rna.ensemble import (UNVALIDATED, VALIDATED, DeclaredGibbsPrior, ViennaPrior,
                               ensemble_bracket, expected_residual, mean_bracket,
                               predict_ensemble)
from scmp.rna.features import (ChemistryFoldingRefused, ContextSubstitutionRefused,
                               ExonUnionRefused, assert_chemistry_not_folded,
                               build_covariates, make_declared_window,
                               per_node_covariates, refuse_exon_union,
                               windows_for_isoforms)

CHEM = ["unmodified", "2p-OMe", "2p-F"]
ASSAY = ["reporter_yfp", "native_quantigene"]
SEQS = ["GGGAAAUCC", "GCGCAAAAGCGC", "GGGGAAAACCCC"]


def setup(seq, gamma=1.0, seed=0, d_hid=4):
    o = RNANonCrossingOracle(seq, 3, True)
    n = len(seq)
    X = build_covariates("GCGCAAAAGCGCAAAAGCGCA"[:21], "unmodified", "reporter_yfp",
                         CHEM, ASSAY)
    Xn = torch.tensor(per_node_covariates(X, n), dtype=torch.float64)
    m = CertifiedModel(Xn.shape[1], o.ground_set(), n, d_hid=d_hid, gamma=gamma,
                       seed=seed)
    g = torch.Generator().manual_seed(seed + 3)
    with torch.no_grad():
        for p in m.parameters():
            p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype))
    return o, m, Xn, backbone_adjacency(n)


def prior_for(o, scale=0.7):
    theta = {e: scale * ((k % 5) - 2) / 2.0 for k, e in enumerate(o.ground_set())}
    return DeclaredGibbsPrior(o, theta)


def exact_p(prior, members):
    c = prior.coefficients()
    lw = np.array([sum(c.get(e, 0.0) for e in m) for m in members])
    w = np.exp(lw - lw.max())
    return w / w.sum()


# ------------------------------------------------------- covariates / conditioning
def test_covariates_condition_on_guide_chemistry_and_assay():
    a = build_covariates("ACGUACGUACGUACGUACGUA", "unmodified", "reporter_yfp",
                         CHEM, ASSAY)
    b = build_covariates("ACGUACGUACGUACGUACGUA", "2p-OMe", "reporter_yfp", CHEM, ASSAY)
    c = build_covariates("ACGUACGUACGUACGUACGUA", "unmodified", "native_quantigene",
                         CHEM, ASSAY)
    d = build_covariates("GCGUACGUACGUACGUACGUA", "unmodified", "reporter_yfp",
                         CHEM, ASSAY)
    assert not np.array_equal(a, b)          # chemistry changes X
    assert not np.array_equal(a, c)          # assay changes X
    assert not np.array_equal(a, d)          # guide changes X
    assert a.min() >= 0.0                    # H4: non-negative covariates


def test_unknown_chemistry_is_an_error_not_a_silent_zero():
    with pytest.raises(KeyError):
        build_covariates("ACGUACGUACGUACGUACGUA", "LNA", "reporter_yfp", CHEM, ASSAY)
    with pytest.raises(KeyError):
        build_covariates("ACGUACGUACGUACGUACGUA", "unmodified", "qPCR", CHEM, ASSAY)


def test_modified_chemistry_is_never_folded():
    assert_chemistry_not_folded("unmodified", "target window")
    for chem in ("2p-OMe", "2p-F", "PS-backbone"):
        with pytest.raises(ChemistryFoldingRefused):
            assert_chemistry_not_folded(chem, "guide")


# ------------------------------------------------------------- transcript mapping
def test_isoforms_are_modelled_explicitly_one_window_each():
    site = "GCGCAAAAGCGC"
    isoforms = [
        {"accession": "NM_1", "isoform": "a", "species": "human", "strand": "+",
         "sequence": "AAAA" + site + "UUUUUUUU"},
        {"accession": "NM_2", "isoform": "b", "species": "human", "strand": "+",
         "sequence": "GGGGGG" + site + "CCCCCC"},
    ]
    ws = windows_for_isoforms(site, isoforms, 20, context_source="native_transcript")
    assert len(ws) == 2
    assert {w.accession for w in ws} == {"NM_1", "NM_2"}
    assert all(w.isoform is not None and w.species and w.strand for w in ws)
    assert ws[0].sequence != ws[1].sequence      # not collapsed into one window


def test_exon_unions_are_refused():
    with pytest.raises(ExonUnionRefused):
        refuse_exon_union()


def test_reporter_context_refuses_native_substitution_by_default():
    site = "GCGCAAAAGCGC"
    iso = [{"accession": "NM_1", "isoform": "a", "species": "human", "strand": "+",
            "sequence": "AAAA" + site + "UUUUUUUU"}]
    with pytest.raises(ContextSubstitutionRefused):
        windows_for_isoforms(site, iso, 20, context_source="reporter_construct")
    ws = windows_for_isoforms(site, iso, 20, context_source="reporter_construct",
                              allow_context_substitution=True,
                              substitution_justification="explicit override for a test")
    assert ws and ws[0].provenance["substituted"] is True
    assert ws[0].provenance["justification"]


@pytest.mark.parametrize("length", [50, 100, 150])
def test_window_length_sensitivity_configs(length):
    seq = ("GCGCAAAAGCGC" * 30)[:length + 40]
    site = "GCGCAAAAGCGC"
    iso = [{"accession": "NM_1", "isoform": "a", "species": "human", "strand": "+",
            "sequence": seq}]
    ws = windows_for_isoforms(site, iso, length, context_source="native_transcript")
    assert ws and all(w.length == length and len(w.sequence) == length for w in ws)


def test_declared_window_carries_no_transcript_claim():
    w = make_declared_window("GGGAAAUCC", 9)
    assert w.context_source == "declared"
    assert w.accession is None and w.isoform is None


# --------------------------------------------------------------------- marginals
@pytest.mark.parametrize("seq", SEQS)
def test_declared_marginals_are_exact_and_normalised(seq):
    o = RNANonCrossingOracle(seq, 3, True)
    prior = prior_for(o)
    M = prior.marginals()
    assert M.status == VALIDATED
    assert M.normalisation_ok
    assert np.all(M.mu >= -1e-12) and np.all(M.mu <= 1 + 1e-12)
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    p = exact_p(prior, members)
    exact = np.array([sum(pi for pi, m in zip(p, members) if e in m)
                      for e in o.ground_set()])
    assert np.max(np.abs(M.mu - exact)) < 1e-9


def test_vienna_marginals_are_flagged_unvalidated():
    o = RNANonCrossingOracle("GCGCAAAAGCGC", 3, True)
    M = ViennaPrior("GCGCAAAAGCGC", o).marginals()
    assert M.status == UNVALIDATED
    assert "not validated" in M.notes.get("note", "")
    assert M.normalisation_ok


def test_no_probability_floor_deletes_any_canonical_pair():
    seq = "GGGGAAAACCCC"
    o = RNANonCrossingOracle(seq, 3, True)
    CAN = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
    want = [(i, j) for i in range(len(seq)) for j in range(i + 4, len(seq))
            if (seq[i], seq[j]) in CAN]
    assert list(o.ground_set()) == want


# ------------------------------------------------------------------ valid samples
@pytest.mark.parametrize("seq", SEQS)
def test_every_declared_sample_is_a_valid_family_member(seq):
    o = RNANonCrossingOracle(seq, 3, True)
    prior = prior_for(o)
    members = {tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)}
    for s in prior.sample(40, random.Random(0)):
        assert rna_is_valid(s, seq, 3, True)
        assert tuple(sorted(s)) in members
        assert o.feasibility(forced=s).feasible


def test_sampler_frequencies_track_the_exact_distribution():
    seq = "GGGAAAUCC"
    o = RNANonCrossingOracle(seq, 3, True)
    prior = prior_for(o)
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    p = dict(zip(members, exact_p(prior, members)))
    draws = prior.sample(600, random.Random(1))
    from collections import Counter
    c = Counter(draws)
    for m, want in p.items():
        if want > 0.05:
            assert abs(c[m] / len(draws) - want) < 0.08


# ----------------------------------------------------------- tiny exact means
@pytest.mark.parametrize("seq", SEQS)
@pytest.mark.parametrize("gamma", [0.0, 1.0])
def test_yhat_equals_the_exact_ensemble_mean_of_f(seq, gamma):
    """yhat = b + a^T mu + gamma E[g] must equal E_p[f], by linearity of the anchor."""
    o, m, X, B = setup(seq, gamma=gamma)
    prior = prior_for(o)
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    p = exact_p(prior, members)
    mu = prior.marginals().mu
    with torch.no_grad():
        fs = np.array([float(m(X, m.indicator(s), B)) for s in members])
        gs = (np.array([float(m.g(X, m.adjacency(m.indicator(s), B)))
                        for s in members]) if m.residual_active
              else np.zeros(len(members)))
    assert predict_ensemble(m, X, mu, float(p @ gs)) == pytest.approx(float(p @ fs),
                                                                     abs=1e-9)


@pytest.mark.parametrize("seq", SEQS)
def test_deterministic_mean_bracket_contains_the_exact_mean(seq):
    o, m, X, B = setup(seq, gamma=1.0)
    prior = prior_for(o)
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    p = exact_p(prior, members)
    mu = prior.marginals().mu
    with torch.no_grad():
        gs = np.array([float(m.g(X, m.adjacency(m.indicator(s), B))) for s in members])
    lo, hi = mean_bracket(m, X, B, mu)
    assert lo - 1e-9 <= float(p @ gs) <= hi + 1e-9
    ylo, yhi = ensemble_bracket(m, X, mu, (lo, hi))
    with torch.no_grad():
        fs = np.array([float(m(X, m.indicator(s), B)) for s in members])
    assert ylo - 1e-9 <= float(p @ fs) <= yhi + 1e-9


def test_bracket_needs_no_sampling_and_gamma_zero_is_degenerate():
    o, m, X, B = setup("GGGAAAUCC", gamma=0.0)
    assert mean_bracket(m, X, B, prior_for(o).marginals().mu) == (0.0, 0.0)


@pytest.mark.parametrize("seq", ["GGGAAAUCC", "GCGCAAAAGCGC"])
def test_sampled_mean_agrees_with_exact_within_recorded_error(seq):
    o, m, X, B = setup(seq, gamma=1.0)
    prior = prior_for(o)
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    p = exact_p(prior, members)
    with torch.no_grad():
        gs = np.array([float(m.g(X, m.adjacency(m.indicator(s), B))) for s in members])
    exact = float(p @ gs)
    E = expected_residual(m, X, B, prior, 256, random.Random(2))
    assert E.n_valid == E.n_samples
    assert E.method == "monte_carlo"
    assert E.stderr == E.stderr and E.stderr >= 0.0     # recorded, not dropped
    assert abs(E.mean - exact) <= 6 * E.stderr + 1e-9


def test_no_latent_structure_is_labelled_with_the_ensemble_outcome():
    """The fitting target is a function of mu and E[g] only; nothing in the ensemble
    API accepts a per-structure label."""
    import inspect
    for fn in (predict_ensemble, ensemble_bracket, mean_bracket):
        params = set(inspect.signature(fn).parameters)
        assert not ({"y", "label", "target"} & params)


# ------------------------------------------------------------ numerical enclosure
def test_ensemble_prediction_has_a_validated_enclosure():
    o, m, X, B = setup("GGGAAAUCC", gamma=1.0)
    mu = prior_for(o).marginals().mu
    a = m.additive_coefficients(X)
    avec = [a[e] for e in m.candidates]
    enc = idot(avec, list(mu)) + Interval.exact(m.anchor_bias(X))
    naive = float(np.array(avec) @ mu + m.anchor_bias(X))
    assert enc.lo <= naive <= enc.hi
    assert enc.width >= 0.0


# ------------------------------------------------------------ stochastic backtrace
@pytest.mark.parametrize("seq", SEQS)
def test_backtrace_and_sequential_samplers_agree_with_the_exact_distribution(seq):
    """Two independent exact samplers, checked against enumeration, not each other."""
    o = RNANonCrossingOracle(seq, 3, True)
    prior = prior_for(o)
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    exact = dict(zip(members, exact_p(prior, members)))
    from collections import Counter
    for method in ("backtrace", "sequential"):
        n = 3000 if method == "backtrace" else 400
        draws = prior.sample(n, random.Random(4), method=method)
        assert all(rna_is_valid(s, seq, 3, True) for s in draws)
        c = Counter(draws)
        tv = 0.5 * sum(abs(c[m] / n - exact[m]) for m in members)
        assert tv < 0.15, (method, tv)


@pytest.mark.parametrize("length", [50, 100, 150])
def test_backtrace_sampling_is_tractable_at_the_declared_window_lengths(length):
    """The sequential sampler needs |E| partition calls per draw and is hopeless here;
    backtrace needs one inside pass."""
    seq = ("GCGCAAAAGCGCAUAUAUAU" * 20)[:length]
    o = RNANonCrossingOracle(seq, 3, True)
    prior = DeclaredGibbsPrior(o, {e: 0.0 for e in o.ground_set()})
    draws = prior.sample(64, random.Random(0))
    assert len(draws) == 64
    assert all(rna_is_valid(s, seq, 3, True) for s in draws)
    assert all(o.feasibility(forced=s).feasible for s in draws[:8])


def test_backtrace_respects_forced_and_forbidden_masks():
    seq = "GCGCAAAAGCGC"
    o = RNANonCrossingOracle(seq, 3, True)
    prior = prior_for(o)
    e_force, e_forbid = o.ground_set()[0], o.ground_set()[-1]
    for s in prior.sample(50, random.Random(5), forced=(e_force,),
                          forbidden=(e_forbid,)):
        assert e_force in s and e_forbid not in s
        assert rna_is_valid(s, seq, 3, True)
