"""Conditional reporter comparison using the unchanged CertifiedModel modules.

This module is resource-gated by p4_decision. The historical 1,153 rows are
intersected with a consistent transcript/accession contract; no excluded row is
filled using a plausible sequence. Evaluation is explicitly exploratory.
"""
from __future__ import annotations
import copy
import numpy as np
import torch
from ..model import CertifiedModel, backbone_adjacency
from ..oracles.rna_noncrossing import RNANonCrossingOracle
from ..rna.ensemble import DeclaredGibbsPrior
from .p4_capability import ensemble_forward, fit_ridge, array_hash, state_hash
from .p4_data import audit_rows, components, partition
from .p4_metrics import metrics, paired_comparison
from .p4_resources import inventory
from .predictor_pilot import _zero_init_residual


def onehot(s):
    return np.array([[float(c == b) for b in 'ACGU'] for c in s])


def window_input(row, sequence, temperature, sample_count, seed, shared=None, hidden=12):
    oracle = RNANonCrossingOracle(sequence, 3, True)
    edges = oracle.ground_set()
    prior = DeclaredGibbsPrior(oracle, {e: .5/temperature for e in edges})
    mu = torch.as_tensor(prior.marginals().mu, dtype=torch.float64)
    guide = onehot(row['guide']).reshape(-1)
    pos = np.arange(len(sequence))/max(1, len(sequence)-1)
    X = torch.as_tensor(np.column_stack([onehot(sequence), np.tile(guide, (len(sequence),1)),
                                        pos, pos**2]), dtype=torch.float64)
    model = CertifiedModel(X.shape[1], edges, len(sequence), d_hid=hidden, gamma=0, seed=seed)
    if shared is not None:
        model.anchor, model.g = shared.anchor, shared.g
        model.set_gamma(shared.gamma)
    B = backbone_adjacency(len(sequence))
    samples = prior.sample(sample_count, np.random.RandomState(seed))
    for s in samples:
        if not oracle.feasibility(forced=s, forbidden=tuple(set(edges)-set(s))).feasible:
            raise ValueError('sample is not a complete valid member')
    Z = torch.stack([model.indicator(s) for s in samples])
    As = torch.stack([model.adjacency(z, B) for z in Z])
    return model, X, mu, Z, As, prior


def prediction(model, rows, temp, count, seed, budget, alternative=False):
    values, ses = [], []
    for i, r in enumerate(rows):
        budget.check()
        seq = r['alternative_window'] if alternative else r['window']
        m, X, mu, Z, As, prior = window_input(r, seq, temp, count, seed+i, model)
        with torch.no_grad():
            values.append(float(ensemble_forward(m, X, mu, Z, As)))
            if m.residual_active:
                vals = torch.stack([m.g(X, A) for A in As])
                ses.append(float(abs(m.gamma)*vals.std(unbiased=True)/np.sqrt(count)))
            else:
                ses.append(0.)
    return np.array(values), ses


def run_rna(cfg, budget):
    c = cfg['rna']
    audit = audit_rows(cfg)
    rows = [r for r in audit['rows'] if r['legacy_status'] == 'pilot_included' and
            r['context_status'] == 'inferred_contiguous_native_insert']
    allowance = inventory(cfg)
    # The old instance cap is not silently extended to a larger RNA dataset.
    required = len(rows)+cfg['capability']['n_instances']+2*(cfg['verifier']['development_instances']+cfg['verifier']['test_instances'])
    if not allowance['admitted'] or required > allowance.get('verified_spendable_instances', 0):
        return {'status': 'BLOCKED_RESOURCE_OR_INSTANCE_ALLOWANCE', 'required_observations': required,
                'conditional_rows': len(rows), 'training_executed': False}
    groups = components([r['guide'] for r in rows])
    for r, g in zip(rows, groups):
        r['global_group'] = g
    tr, dv, te = partition(groups, c['seed'], c['train_fraction'], c['dev_fraction'])
    if min(len(tr), len(dv), len(te)) < 20:
        return {'status': 'TOO_FEW_ROWS_AFTER_GLOBAL_GROUP_SPLIT', 'sizes': [len(tr),len(dv),len(te)]}
    y = np.array([r['efficacy'] for r in rows])
    G = np.stack([onehot(r['guide']).reshape(-1) for r in rows])
    W = np.stack([onehot(r['window']).reshape(-1) for r in rows])
    ij = np.triu_indices(G.shape[1], 1)
    designs = {'guide_ridge': G, 'guide_pairwise_ridge': np.column_stack([G,G[:,ij[0]]*G[:,ij[1]]]),
               'guide_window_ridge': np.column_stack([G,W])}
    data, summary = [], []
    for i, r in enumerate(rows):
        budget.check()
        m, X, mu, Z, As, prior = window_input(r, r['window'], 1., c['train_graph_samples'], c['seed']+i,
                                            hidden=c['hidden'])
        data.append((X, mu, prior, m.candidates))
        edges = m.candidates
        summary.append([float(mu.sum())/50, np.mean([b-a for a,b in edges])/50 if edges else 0.,
                        float(mu.mean()) if len(mu) else 0.,len(edges)/50])
    designs['summary_ridge'] = np.column_stack([G,np.array(summary)])
    preds, fits = {}, {}
    for arm, D in designs.items():
        choices = []
        for lam in c['ridge_grid']:
            budget.check()
            p, coef = fit_ridge(D,y,tr,lam)
            choices.append((float(np.mean((p[dv]-y[dv])**2)),lam,p,coef))
        _, lam, p, coef = min(choices,key=lambda x:(x[0],x[1]))
        preds[arm] = p
        fits[arm] = {'lambda':lam,'fitted_coefficients':len(coef),'prediction_sha256':array_hash(p),
                     'development_grid': [{'lambda':v[1],'mse':v[0]} for v in choices]}
    X, mu, prior, edges = data[0]
    model = CertifiedModel(X.shape[1], edges,50,d_hid=c['hidden'],gamma=0,seed=c['seed'])
    # Each row has its own candidates and actual sampled adjacency; only weights
    # are shared. The graph family is never replaced with summary features.
    train_rng = np.random.RandomState(c['seed']+20000)
    for phase, epochs in [('A',c['anchor_epochs']),('B',c['residual_epochs'])]:
        if phase == 'B':
            with torch.no_grad():
                preds['actual_A'], _ = prediction(model,rows,1.,16,c['seed']+30000,budget)
            anchor = copy.deepcopy(model)
            _zero_init_residual(model)
            model.set_gamma(1.)
            for p in model.anchor.parameters():
                p.requires_grad_(False)
        params = model.anchor.parameters() if phase == 'A' else model.g.parameters()
        opt = torch.optim.Adam(params,lr=c['learning_rate'])
        for epoch in range(epochs):
            order = train_rng.permutation(tr)
            for offset in range(0,len(order),c['batch_size']):
                budget.check()
                indices = order[offset:offset+c['batch_size']]
                opt.zero_grad()
                loss = 0.
                for i in indices:
                    X, mu, prior, edges = data[i]
                    m = CertifiedModel(X.shape[1],edges,50,d_hid=c['hidden'],gamma=model.gamma)
                    m.anchor, m.g = model.anchor, model.g
                    B = backbone_adjacency(50)
                    samples = prior.sample(c['train_graph_samples'],train_rng)
                    Z = torch.stack([m.indicator(s) for s in samples])
                    As = torch.stack([m.adjacency(z,B) for z in Z])
                    value = ensemble_forward(m,X,mu,Z,As)
                    loss = loss + (value-y[i])**2/len(indices)
                loss.backward()
                opt.step()
    raw, se = prediction(model,rows,1.,64,c['seed']+30000,budget)
    A = preds['actual_A']
    grid = {str(g):float(np.mean((A[dv]+g*(raw[dv]-A[dv])-y[dv])**2)) for g in cfg['capability']['gamma_grid']}
    gamma = min(cfg['capability']['gamma_grid'],key=lambda g:(grid[str(g)],g))
    model.set_gamma(gamma)
    preds['actual_A_plus_B'] = A+gamma*(raw-A)
    ref = min(['guide_ridge','guide_pairwise_ridge','guide_window_ridge'],
              key=lambda arm:np.mean((preds[arm][dv]-y[dv])**2))
    comparisons = {arm:paired_comparison(rows,te,p,preds[ref],cfg['metrics']) for arm,p in preds.items()}
    sensitivity = {}
    for temp in c['prior_temperatures']:
        p, err = prediction(model,rows,temp,64,c['seed']+30000,budget)
        sensitivity[str(temp)] = {'mean_prediction_se':float(np.mean(err)),
                                  'prediction_sha256':array_hash(p), 'predictions':p.tolist(),
                                  'max_prediction_change':float(np.max(abs(p-preds['actual_A_plus_B']))),
                                  'ranking':paired_comparison(rows,te,p,preds['actual_A_plus_B'],cfg['metrics'])}
    alt, alt_se = prediction(model,rows,1.,64,c['seed']+30000,budget,alternative=True)
    other_draw, other_se = prediction(model,rows,1.,64,c['seed']+40000,budget)
    identical = np.array([i for i in te if rows[i]['isoform_class']=='identical_local_windows'])
    return {'status':'COMPLETE_EXPLORATORY_CONDITIONAL_REPORTER','reference_selected_on_development':ref,
            'comparisons':comparisons,'ridge_fits':fits,'gamma_grid':grid,'gamma':gamma,
            'predictions':{k:v.tolist() for k,v in preds.items()},'row_ids':[r['row_id'] for r in rows],
            'split':{'train':tr.tolist(),'dev':dv.tolist(),'test':te.tolist(),'groups':groups},
            'prior_sensitivity':sensitivity,'isoform_alternative':paired_comparison(rows,te,alt,preds['actual_A_plus_B'],cfg['metrics']),
            'identical_window_only':paired_comparison(rows,identical,preds['actual_A_plus_B'],preds[ref],cfg['metrics']) if len(identical) else None,
            'eval_sampling_replicate_max_prediction_change':float(np.max(abs(other_draw-preds['actual_A_plus_B']))),
            'gradient_sampling':'independent valid backtraces resampled each minibatch; gradient variance not estimated',
            'sampling_se_units':'predicted inhibition fraction, including gamma','model_sha256':state_hash(model),
            'chemical_features':None,'assay_identity':'YFP / H1299 / 48 h','raw_labels_clipped':False,
            'context_status':'inferred_contiguous_native_insert','biological_order_certification':False}
