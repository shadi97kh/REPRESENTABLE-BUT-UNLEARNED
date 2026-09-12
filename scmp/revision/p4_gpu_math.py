"""P4 execution-only acceleration, without changing CertifiedModel or its prior.

Inside/outside is the derivative of the existing unconditioned log partition.
Sampling uses its exact interval-production probabilities. No learned prior,
normalization layer, readout rescaling, or replacement graph features are used.
"""
from __future__ import annotations
import numpy as np
import torch
from ..oracles.rna_noncrossing import RNANonCrossingOracle
from ..model import backbone_adjacency


class FastPrior:
    def __init__(self, sequence, temperature=1.):
        self.oracle = RNANonCrossingOracle(sequence, 3, True)
        self.edges = tuple(self.oracle.ground_set())
        self.index = {e: k for k, e in enumerate(self.edges)}
        self.n = len(sequence)
        theta = .5 / temperature
        self.chart = self.oracle.log_inside({e: theta for e in self.edges})
        W = np.asarray(self.chart['W'])
        self.mu = np.zeros(len(self.edges))
        reach = np.zeros_like(W)
        reach[0, self.n] = 1.
        self.productions = {}
        pairable = {i: [k for a, k in self.edges if a == i] for i in range(self.n)}
        for length in range(self.n, 0, -1):
            for i in range(self.n-length+1):
                j = i+length  # exclusive endpoint
                ks = np.array([k for k in pairable[i] if k < j], dtype=int)
                pair_p = np.exp(theta + W[i+1, ks] + W[ks+1, j] - W[i, j])
                unpaired = np.exp(W[i+1, j]-W[i, j])
                probs = np.r_[unpaired, pair_p]
                if abs(probs.sum()-1.) > 1e-11:
                    raise ValueError('interval probabilities do not normalize')
                cdf = np.cumsum(probs)
                cdf[-1] = 1.
                self.productions[i, j] = (np.r_[-1, ks], cdf)
                r = reach[i, j]
                reach[i+1, j] += r*unpaired
                for k, p in zip(ks, pair_p):
                    v = r*p
                    self.mu[self.index[i, int(k)]] += v
                    reach[i+1, k] += v
                    reach[k+1, j] += v

    def sample(self, count, rng):
        samples = []
        for _ in range(count):
            stack, edges = [(0, self.n)], []
            while stack:
                i, j = stack.pop()
                if i >= j:
                    continue
                ks, cdf = self.productions[i, j]
                k = int(ks[np.searchsorted(cdf, rng.random_sample(), side='right')])
                if k < 0:
                    stack.append((i+1, j))
                else:
                    edges.append((i, k))
                    # Same inside-before-outside traversal as the original sampler.
                    stack.extend([(k+1, j), (i+1, k)])
            samples.append(tuple(sorted(edges)))
        return samples


def onehot(s):
    return np.array([[float(c == b) for b in 'ACGU'] for c in s])


def row_features(row, sequence):
    guide = onehot(row['guide']).reshape(-1)
    pos = np.arange(len(sequence))/max(1, len(sequence)-1)
    return np.column_stack([onehot(sequence), np.tile(guide, (len(sequence), 1)), pos, pos**2])


def pack(features, priors, device):
    X = torch.as_tensor(np.stack(features), dtype=torch.float64, device=device)
    max_edges = max(len(p.edges) for p in priors)
    ii = np.zeros((len(priors), max_edges), dtype=int)
    jj, mu = ii.copy(), np.zeros(ii.shape)
    for row, p in enumerate(priors):
        e = np.asarray(p.edges, dtype=int).reshape(-1, 2)
        ii[row, :len(e)], jj[row, :len(e)] = e[:, 0], e[:, 1]
        mu[row, :len(e)] = p.mu
    return (X, torch.as_tensor(ii, device=device), torch.as_tensor(jj, device=device),
            torch.as_tensor(mu, device=device))


def adjacencies(priors, count, rng, device):
    n = priors[0].n
    B = backbone_adjacency(n).numpy()
    As = np.broadcast_to(B, (len(priors), count, n, n)).copy()
    for row, prior in enumerate(priors):
        for k, s in enumerate(prior.sample(count, rng)):
            if s:
                a, b = np.asarray(s).T
                As[row, k, a, b] += 1.
                As[row, k, b, a] += 1.
    return torch.as_tensor(As, dtype=torch.float64, device=device)


def batched_forward(model, packed, As=None, return_draws=False):
    X, ii, jj, mu = packed
    row = torch.arange(len(X), device=X.device)[:, None]
    xi, xj = X[row, ii], X[row, jj]
    anchor = model.anchor
    b = anchor.psi(X).sum(dim=(1, 2))
    h = anchor.phi(torch.cat([xi, xj], -1)) + anchor.phi(torch.cat([xj, xi], -1))
    value = b + (anchor.a_out(h).squeeze(-1)*mu).sum(dim=1)
    if not model.residual_active:
        return (value, None) if return_draws else value
    g = model.g
    H = torch.relu(As @ g.lin0(X)[:, None])
    H = torch.relu(As @ g.lin1(H))
    draws = g.readout(H.sum(dim=2)).squeeze(-1)
    value = value + model.gamma*draws.mean(dim=1)
    return (value, draws) if return_draws else value


def validate(device):
    """Compare exact marginals, sample paths, outputs and gradients independently."""
    from ..rna.ensemble import DeclaredGibbsPrior
    from ..model import CertifiedModel
    from .p4_capability import ensemble_forward
    errors, same_samples = [], []
    for seq in ['GCGCAAAAGCGC', 'GGGGAAAACCCC', 'A'*12, 'ACGU'*6]:
        for t in [.5, 1., 2.]:
            p = FastPrior(seq, t)
            original = DeclaredGibbsPrior(p.oracle, {e: .5/t for e in p.edges})
            expected = np.asarray(original.marginals().mu)
            errors.append(float(np.max(abs(p.mu-expected), initial=0)))
            samples = p.sample(20, np.random.RandomState(17))
            ref = original.sample(20, np.random.RandomState(17))
            same_samples.append([set(s) for s in samples] == [set(s) for s in ref])
            for s in samples:
                assert p.oracle.feasibility(forced=s, forbidden=tuple(set(p.edges)-set(s))).feasible
    ps = [FastPrior('GCGCAAAAGCGC'), FastPrior('GGGGAAAACCCC')]
    fs = [row_features({'guide': 'ACGUACGUACGUACGUACGUA'}, s)
          for s in ['GCGCAAAAGCGC', 'GGGGAAAACCCC']]
    m = CertifiedModel(90, ps[0].edges, 12, d_hid=12, gamma=1., seed=417).to(device)
    pk = pack(fs, ps, device)
    As = adjacencies(ps, 8, np.random.RandomState(911), device)
    batched = batched_forward(m, pk, As)
    sequential = []
    for i, p in enumerate(ps):
        local = CertifiedModel(90, p.edges, 12, d_hid=12, gamma=1.).to(device)
        local.anchor, local.g = m.anchor, m.g
        sequential.append(ensemble_forward(local, pk[0][i], torch.as_tensor(p.mu, device=device), None, As[i]))
    sequential = torch.stack(sequential)
    value_error = float((batched-sequential).abs().max().detach().cpu())
    params = tuple(m.parameters())
    gb = torch.autograd.grad((batched-torch.tensor([.2, .7], device=device)).square().mean(), params)
    gs = torch.autograd.grad((sequential-torch.tensor([.2, .7], device=device)).square().mean(), params)
    gradient_error = max(float((a-b).abs().max().cpu()) for a, b in zip(gb, gs))
    m.set_gamma(0)
    with torch.no_grad():
        a = batched_forward(m, pk)
        for param in m.g.parameters():
            param.fill_(float('nan'))
        bypass = torch.equal(a, batched_forward(m, pk))
    assert max(errors) < 1e-10 and all(same_samples)
    assert value_error < 1e-10 and gradient_error < 1e-8 and bypass
    return {'max_marginal_error': max(errors), 'sample_sequences_match_existing_sampler': all(same_samples),
            'batched_output_max_abs_error': value_error, 'batched_gradient_max_abs_error': gradient_error,
            'gamma_zero_exact_nan_bypass': bypass, 'graph_membership_checks': 240,
            'device': str(device), 'optimizer_steps': 0}
