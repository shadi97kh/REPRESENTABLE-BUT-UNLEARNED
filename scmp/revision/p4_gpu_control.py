"""Unchanged preregistered capability optimization on one CUDA device."""
import numpy as np
import torch
from .p4_capability import control_data, fit_ridge, additive_residual, array_hash, state_hash
from .predictor_pilot import _fwd_batched, _zero_init_residual
from .p4_metrics import metrics


def run_control(c, budget, device):
    m, X, o, B, members, Z, As, P, y, tr, dv, te, pf, info = control_data(c)
    z = Z.numpy()
    predictions = {}
    for name, D in [('additive_ridge', P@z), ('pairwise_ridge', P@np.column_stack([z, pf]))]:
        budget.check()
        predictions[name] = fit_ridge(D, y, tr, 1e-6)[0]
    m = m.to(device)
    X, B, Z, As = [v.to(device) for v in (X, B, Z, As)]
    Pt, yt = torch.as_tensor(P, device=device), torch.as_tensor(y, device=device)
    train = torch.as_tensor(tr, device=device)
    optimizer_steps = 0
    for phase, epochs, lr in [('A', c['anchor_epochs'], c['lr_anchor']),
                              ('B', c['residual_epochs'], c['lr_residual'])]:
        if phase == 'A':
            m.set_gamma(0)
        else:
            with torch.no_grad():
                predictions['actual_A'] = (Pt@_fwd_batched(m, X, B, Z, As)).cpu().numpy()
            _zero_init_residual(m)
            m.set_gamma(1)
            for p in m.anchor.parameters():
                p.requires_grad_(False)
        opt = torch.optim.Adam(m.anchor.parameters() if phase == 'A' else m.g.parameters(), lr=lr)
        for epoch in range(epochs):
            budget.check()
            opt.zero_grad()
            pred = Pt@_fwd_batched(m, X, B, Z, As)
            loss = ((pred[train]-yt[train])**2).mean()
            loss.backward()
            opt.step()
            optimizer_steps += 1
            if (epoch+1) % 100 == 0:
                print(f'capability {phase} epoch {epoch+1}/{epochs} train_mse={float(loss.detach().cpu()):.8g}', flush=True)
    with torch.no_grad():
        raw = (Pt@_fwd_batched(m, X, B, Z, As)).cpu().numpy()
    a = predictions['actual_A']
    grid = {str(g): float(np.mean((a[dv]+g*(raw[dv]-a[dv])-y[dv])**2)) for g in c['gamma_grid']}
    chosen = min(c['gamma_grid'], key=lambda g: (grid[str(g)], g))
    m.set_gamma(chosen)
    predictions['actual_A_plus_B'] = a+chosen*(raw-a)
    rows = [{'arm': name, 'train': metrics(p[tr], y[tr]), 'development': metrics(p[dv], y[dv]),
             'test': metrics(p[te], y[te]), 'prediction_sha256': array_hash(p),
             'predictions': p.tolist()} for name, p in predictions.items()]
    with torch.no_grad():
        projected = additive_residual(z, _fwd_batched(m, X, B, Z, As).cpu().numpy())
    by_name = {r['arm']: r for r in rows}
    fit_ok = by_name['actual_A_plus_B']['train']['r2'] >= c['fit_gate_train_r2']
    utility = by_name['actual_A_plus_B']['test']['mse'] <= by_name['pairwise_ridge']['test']['mse']+c['heldout_noninferiority_mse']
    info.update(rows=rows, chosen_gamma=chosen, development_gamma_grid=grid,
                trained_model_nonadditivity=projected, fit_gate=fit_ok, heldout_utility_gate=utility,
                capability_passed=fit_ok and utility and projected['max_abs'] > c['nonadditivity_tolerance'],
                model_sha256=state_hash(m), labels=y.tolist(), mc_error='zero: all valid members enumerated',
                optimizer_steps=optimizer_steps, device=str(device), readout_scale_modified=False,
                interpretation='capacity/optimization control, not therapeutic efficacy or novelty')
    return info, m
