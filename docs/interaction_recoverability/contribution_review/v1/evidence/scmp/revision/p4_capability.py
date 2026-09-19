"""Actual-model ensemble objective; one prespecified independent control."""
from __future__ import annotations
import copy
import hashlib
import itertools
import time
import numpy as np
import torch
from ..model import CertifiedModel, backbone_adjacency
from ..oracles.bruteforce import enumerate_rna
from ..oracles.rna_noncrossing import RNANonCrossingOracle
from .predictor_pilot import _fwd_batched, _zero_init_residual


def array_hash(x):
    x = np.ascontiguousarray(x, dtype='<f8')
    return hashlib.sha256(str(x.shape).encode() + x.tobytes()).hexdigest()


def state_hash(m):
    h = hashlib.sha256()
    for name, value in sorted(m.state_dict().items()):
        h.update(name.encode())
        h.update(value.detach().cpu().numpy().tobytes())
    h.update(str(m.gamma).encode())
    return h.hexdigest()


def make_inputs(sequence, hidden, seed):
    oracle = RNANonCrossingOracle(sequence, 3, True)
    members = sorted({tuple(sorted(s)) for s in enumerate_rna(sequence, 3, True)})
    assert len(members) == oracle.count().value
    X = torch.zeros(len(sequence), 6, dtype=torch.float64)
    for i, c in enumerate(sequence):
        X[i, 'ACGU'.index(c)] = 1
        X[i, 4] = i/max(1, len(sequence)-1)
        X[i, 5] = X[i, 4]**2
    m = CertifiedModel(6, oracle.ground_set(), len(sequence), d_hid=hidden, gamma=1, seed=seed)
    B = backbone_adjacency(len(sequence))
    Z = torch.stack([m.indicator(s) for s in members])
    As = torch.stack([m.adjacency(z, B) for z in Z])
    return m, X, oracle, B, members, Z, As


def ensemble_forward(model, X, mu, Z, As, weights=None):
    """Differentiable b+a^T mu+gamma E[g], without modifying readout scale.

    As contains actual B+PZ graphs. mu may be exact while only g is sampled.
    The return value is one prediction for one ensemble observation.
    """
    b, a = model.anchor_terms(X)
    value = b + mu @ a
    if not model.residual_active:
        return value
    g = model.g
    H = torch.relu(As @ g.lin0(X))
    H = torch.relu(As @ g.lin1(H))
    vals = g.readout(H.sum(dim=1)).squeeze(-1)
    return value + model.gamma * (vals.mean() if weights is None else weights @ vals)


def additive_residual(Z, values):
    D = np.column_stack([np.ones(len(Z)), np.asarray(Z)])
    residual = values-D@np.linalg.lstsq(D, values, rcond=None)[0]
    return {'rmse': float(np.sqrt(np.mean(residual**2))),
            'max_abs': float(np.max(np.abs(residual))), 'rank': int(np.linalg.matrix_rank(D))}


def control_data(c):
    m, X, o, B, members, Z, As = make_inputs(c['family'], c['hidden'], c['seed'])
    z = Z.numpy()
    triple = next((p for p in itertools.combinations(range(z.shape[1]), 3)
                   if np.any(z[:, p].sum(axis=1) == 3)), None)
    if triple is None:
        raise ValueError('prespecified family has no feasible triple motif')
    h = .2*z.sum(axis=1)/z.sum(axis=1).max() + .8*(z[:, triple].sum(axis=1) >= 2)
    a, b = [m.candidates[i] for i in triple[:2]]
    inds = [members.index(s) for s in [(), (a,), (b,), tuple(sorted((a, b)))]]
    delta = np.zeros(len(members))
    delta[inds] = [1., -1., -1., 1.]
    assert abs(delta.sum()) < 1e-14 and np.max(abs(delta@z)) < 1e-14
    rng = np.random.RandomState(c['distribution_seed'])
    distributions, groups = [], []
    for j in range(c['n_instances']//2):
        p = rng.dirichlet(np.ones(len(members))*.2)*.5
        p[inds] += .125
        eps = min(p[inds])*.8
        distributions.extend([p+eps*delta, p-eps*delta])
        groups.extend([j, j])
    P = np.array(distributions)
    y = P@h + rng.normal(0, c['noise_sd'], len(P))
    from .p4_data import partition
    tr, dv, te = partition(groups, c['distribution_seed']+1, c['train_fraction'], c['dev_fraction'])
    pair_f = np.array([z[:, i]*z[:, j] for i, j in itertools.combinations(range(z.shape[1]), 2)]).T
    info = {'family_members': len(members), 'candidates': len(m.candidates),
            'motif_edges': [list(m.candidates[i]) for i in triple],
            'equal_marginal_pair_max_error': float(np.max(np.abs((P[::2]-P[1::2])@z))),
            'equal_marginal_label_difference_min': float(np.min(abs((P[::2]-P[1::2])@h))),
            'all_probabilities_nonnegative': bool(np.all(P >= 0)),
            'max_normalization_error': float(np.max(abs(P.sum(axis=1)-1))),
            'target_additive_projection': additive_residual(z, h),
            'target_pairwise_projection': additive_residual(np.column_stack([z, pair_f]), h),
            'label_noise_sd': c['noise_sd'], 'labels_sha256': array_hash(y),
            'distribution_sha256': array_hash(P),
            'input_contract': 'X fixed; observed distributions enter through member weights, not only first moments',
            'train_ids': tr.tolist(), 'dev_ids': dv.tolist(), 'test_ids': te.tolist(),
            'group_ids': groups, 'independent_groups': len(set(groups)),
            'split_is_fresh_synthetic': True, 'biological_evidence': False}
    return m, X, o, B, members, Z, As, P, y, tr, dv, te, pair_f, info


def smoke(c, budget):
    """No optimizer or parameter update: measured forward/backward diagnostic."""
    start = time.perf_counter()
    m, X, o, B, members, Z, As, P, y, tr, dv, te, pf, info = control_data(c)
    budget.check()
    before = state_hash(m)
    weights = torch.as_tensor(P[0], dtype=torch.float64)
    mu = weights@Z
    with torch.no_grad():
        loop = torch.stack([m(X, z, B) for z in Z])
        batched = _fwd_batched(m, X, B, Z, As)
        same = float(torch.max(torch.abs(loop-batched)))
        explicit_mean = weights@loop
        direct = ensemble_forward(m, X, mu, Z, As, weights)
    m.set_gamma(0)
    gamma_zero_exact = all(torch.equal(m(X, z, B), m.linear_value(X, z)) for z in Z)
    poisoned = copy.deepcopy(m)
    with torch.no_grad():
        for v in poisoned.g.parameters():
            v.fill_(float('nan'))
    poison_ok = bool(torch.isfinite(ensemble_forward(poisoned, X, mu, Z, As, weights)))
    m.set_gamma(1)
    m.zero_grad()
    step_start = time.perf_counter()
    pred = ensemble_forward(m, X, mu, Z, As, weights)
    ((pred-float(y[0]))**2).backward()
    measured = time.perf_counter()-step_start
    grads = {n: float(p.grad.norm()) if p.grad is not None else None for n, p in m.g.named_parameters()}
    with torch.no_grad():
        per_structure_loss = weights@((loop-y[0])**2)
        mean_loss = (explicit_mean-y[0])**2
    active = additive_residual(Z.numpy(), loop.numpy())
    assert gamma_zero_exact and poison_ok
    assert same < 1e-10 and abs(float(direct-explicit_mean)) < 1e-10
    assert any(v and v > 0 for v in grads.values())
    assert before == state_hash(m), 'smoke must never update parameters'
    info.update(gamma0_exact=gamma_zero_exact, gamma0_nan_bypass=poison_ok,
                forward_loop_vs_batched_max_abs=same,
                ensemble_vs_explicit_mean_abs=abs(float(direct-explicit_mean)),
                distinct_valid_adjacencies=len({array_hash(a.numpy()) for a in As}),
                fixed_model_nonadditivity=active, residual_grad_norms=grads,
                loss_of_mean=float(mean_loss), mean_of_losses=float(per_structure_loss),
                difference_is_structure_variance=float(per_structure_loss-mean_loss),
                parameters=sum(p.numel() for p in m.parameters()),
                parameters_unchanged=before == state_hash(m), model_sha256=before,
                readout_scale_modified=False, optimizer_steps=0,
                measured_one_observation_forward_backward_wall_s=measured,
                estimated_control_training_wall_s=measured*(c['anchor_epochs']+c['residual_epochs'])*4,
                estimate_note='4x measured full-family forward/backward; excludes imports and verifier teacher; hard cap required',
                total_wall_s=time.perf_counter()-start)
    return info


def fit_ridge(D, y, tr, lam=1.):
    D = np.column_stack([np.ones(len(D)), D])
    if D.shape[1] > len(tr):
        coef = D[tr].T@np.linalg.solve(D[tr]@D[tr].T + lam*np.eye(len(tr)), y[tr])
    else:
        coef = np.linalg.solve(D[tr].T@D[tr] + lam*np.eye(D.shape[1]), D[tr].T@y[tr])
    return D@coef, coef


def run_control(c, budget):
    m, X, o, B, members, Z, As, P, y, tr, dv, te, pf, info = control_data(c)
    Pt, yt = torch.as_tensor(P), torch.as_tensor(y)
    predictions = {}
    for name, D in [('additive_ridge', P@Z.numpy()),
                    ('pairwise_ridge', P@np.column_stack([Z.numpy(), pf]))]:
        budget.check()
        predictions[name] = fit_ridge(D, y, tr, 1e-6)[0]
    m.set_gamma(0)
    opt = torch.optim.Adam(m.anchor.parameters(), lr=c['lr_anchor'])
    for _ in range(c['anchor_epochs']):
        budget.check()
        opt.zero_grad()
        pred = Pt@_fwd_batched(m, X, B, Z, As)
        ((pred[tr]-yt[tr])**2).mean().backward()
        opt.step()
    with torch.no_grad():
        predictions['actual_A'] = (Pt@_fwd_batched(m, X, B, Z, As)).numpy()
    _zero_init_residual(m)  # P1 mechanism reused, not a second repair search.
    m.set_gamma(1)
    for p in m.anchor.parameters():
        p.requires_grad_(False)
    opt = torch.optim.Adam(m.g.parameters(), lr=c['lr_residual'])
    for _ in range(c['residual_epochs']):
        budget.check()
        opt.zero_grad()
        pred = Pt@_fwd_batched(m, X, B, Z, As)
        # Average predictions BEFORE the ensemble's loss.
        ((pred[tr]-yt[tr])**2).mean().backward()
        opt.step()
    with torch.no_grad():
        raw = (Pt@_fwd_batched(m, X, B, Z, As)).numpy()
    a = predictions['actual_A']
    grid = {str(g): float(np.mean((a[dv]+g*(raw[dv]-a[dv])-y[dv])**2)) for g in c['gamma_grid']}
    chosen = min(c['gamma_grid'], key=lambda g: (grid[str(g)], g))
    m.set_gamma(chosen)
    predictions['actual_A_plus_B'] = a+chosen*(raw-a)
    from .p4_metrics import metrics
    rows = [{'arm': name, 'train': metrics(p[tr], y[tr]), 'development': metrics(p[dv], y[dv]),
             'test': metrics(p[te], y[te]), 'prediction_sha256': array_hash(p),
             'predictions': p.tolist()} for name, p in predictions.items()]
    with torch.no_grad():
        projected = additive_residual(Z.numpy(), _fwd_batched(m, X, B, Z, As).numpy())
    by_name = {r['arm']: r for r in rows}
    fit_ok = by_name['actual_A_plus_B']['train']['r2'] >= c['fit_gate_train_r2']
    utility = by_name['actual_A_plus_B']['test']['mse'] <= by_name['pairwise_ridge']['test']['mse']+c['heldout_noninferiority_mse']
    info.update(rows=rows, chosen_gamma=chosen, development_gamma_grid=grid,
                trained_model_nonadditivity=projected, fit_gate=fit_ok, heldout_utility_gate=utility,
                capability_passed=fit_ok and utility and projected['max_abs'] > c['nonadditivity_tolerance'],
                model_sha256=state_hash(m), labels=y.tolist(), mc_error='zero: all valid members enumerated',
                interpretation='capacity/optimization control, not therapeutic efficacy or novelty')
    return info, m
