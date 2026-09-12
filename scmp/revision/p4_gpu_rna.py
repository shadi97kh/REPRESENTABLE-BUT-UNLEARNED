"""Execute all six preregistered RNA arms using one GPU and exact DP caches."""
from __future__ import annotations
import copy
import gc
import numpy as np
import torch
from ..model import CertifiedModel
from .p4_gpu_math import FastPrior, onehot, row_features, pack, adjacencies, batched_forward
from .p4_capability import fit_ridge, array_hash, state_hash
from .p4_data import components, partition
from .p4_metrics import paired_comparison, metrics
from .predictor_pilot import _zero_init_residual


def build_data(rows, temperature, budget, alternative=False):
    cache, features, priors = {}, [], []
    for i, row in enumerate(rows):
        budget.check()
        seq = row['alternative_window'] if alternative else row['window']
        if seq not in cache:
            cache[seq] = FastPrior(seq, temperature)
        priors.append(cache[seq])
        features.append(row_features(row, seq))
        if (i+1) % 250 == 0:
            print(f'RNA exact prior T={temperature} alternative={alternative}: {i+1}/{len(rows)}', flush=True)
    return features, priors


def predict(model, features, priors, count, seed, budget, device, batch_size=16):
    values, ses = [], []
    for offset in range(0, len(features), batch_size):
        budget.check()
        ps, fs = priors[offset:offset+batch_size], features[offset:offset+batch_size]
        pk = pack(fs, ps, device)
        # Independent row seeds make draws invariant to evaluation batch size.
        As = torch.cat([adjacencies([p], count, np.random.RandomState(seed+offset+j), device)
                        for j, p in enumerate(ps)], dim=0) if model.residual_active else None
        with torch.no_grad():
            value, draws = batched_forward(model, pk, As, return_draws=True)
            values.extend(value.cpu().tolist())
            ses.extend((abs(model.gamma)*draws.std(dim=1, unbiased=True)/np.sqrt(count)).cpu().tolist()
                       if draws is not None else [0.]*len(ps))
    return np.array(values), np.array(ses)


def run_rna(cfg, audited_rows, budget, device, save_stage, save_model):
    c = cfg['rna']
    rows = [copy.deepcopy(r) for r in audited_rows if r['legacy_status'] == 'pilot_included' and
            r['context_status'] == 'inferred_contiguous_native_insert']
    groups = components([r['guide'] for r in rows])
    for r, g in zip(rows, groups):
        r['global_group'] = g
    tr, dv, te = partition(groups, c['seed'], c['train_fraction'], c['dev_fraction'])
    assert min(len(tr), len(dv), len(te)) >= 20
    group_sets = [set(np.array(groups)[x]) for x in (tr, dv, te)]
    assert not any(group_sets[i] & group_sets[j] for i in range(3) for j in range(i))
    kmers = [{r['guide'][j:j+13] for i in part for r in [rows[i]] for j in range(len(r['guide'])-12)}
             for part in (tr, dv, te)]
    assert not any(kmers[i] & kmers[j] for i in range(3) for j in range(i))
    split = {'train': tr.tolist(), 'dev': dv.tolist(), 'test': te.tolist(), 'groups': groups,
             'group_leak': 0, 'shared_13mer_leak': 0, 'sizes': [len(tr),len(dv),len(te)]}
    print(f'RNA {len(rows)} rows; train/dev/test={split["sizes"]}; shared-13mer leak=0', flush=True)
    y = np.array([r['efficacy'] for r in rows])
    save_stage('rna_inputs', {'row_ids':[r['row_id'] for r in rows], 'split':split,
                            'labels':y.tolist(), 'labels_sha256':array_hash(y),
                            'rows':[{k:r[k] for k in ['row_id','gene','guide','window','alternative_window',
                                                     'isoform_class','efficacy','global_group']} for r in rows]})
    features, priors = build_data(rows, 1., budget)
    G = np.stack([onehot(r['guide']).reshape(-1) for r in rows])
    W = np.stack([onehot(r['window']).reshape(-1) for r in rows])
    ij = np.triu_indices(G.shape[1], 1)
    summary = [[float(p.mu.sum())/50, np.mean([b-a for a,b in p.edges])/50 if p.edges else 0.,
                float(p.mu.mean()) if len(p.mu) else 0., len(p.edges)/50] for p in priors]
    designs = {'guide_ridge':G, 'guide_pairwise_ridge':np.column_stack([G,G[:,ij[0]]*G[:,ij[1]]]),
               'guide_window_ridge':np.column_stack([G,W]), 'summary_ridge':np.column_stack([G,summary])}
    preds, fits = {}, {}
    for arm, D in designs.items():
        choices = []
        for lam in c['ridge_grid']:
            budget.check()
            p, coef = fit_ridge(D, y, tr, lam)
            choices.append((float(np.mean((p[dv]-y[dv])**2)),lam,p,coef))
        _, lam, p, coef = min(choices, key=lambda x:(x[0],x[1]))
        preds[arm] = p
        fits[arm] = {'lambda':lam, 'coefficients':coef.tolist(), 'prediction_sha256':array_hash(p),
                     'development_grid':[{'lambda':v[1], 'mse':v[0]} for v in choices]}
        print(f'RNA fitted {arm}; lambda={lam}', flush=True)
    save_stage('rna_ridges', {'fits':fits, 'predictions':{k:v.tolist() for k,v in preds.items()}, 'split':split})
    del designs, W
    model = CertifiedModel(90, priors[0].edges, 50, d_hid=c['hidden'], gamma=0, seed=c['seed']).to(device)
    train_rng = np.random.RandomState(c['seed']+20000)
    history, steps = [], 0
    for phase, epochs in [('A',c['anchor_epochs']),('B',c['residual_epochs'])]:
        if phase == 'B':
            preds['actual_A'], _ = predict(model,features,priors,16,c['seed']+30000,budget,device)
            save_model('rna_A', model)
            _zero_init_residual(model)
            model.set_gamma(1.)
            for p in model.anchor.parameters():
                p.requires_grad_(False)
        opt = torch.optim.Adam(model.anchor.parameters() if phase == 'A' else model.g.parameters(), lr=c['learning_rate'])
        for epoch in range(epochs):
            order, losses = train_rng.permutation(tr), []
            for offset in range(0,len(order),c['batch_size']):
                budget.check()
                indices = order[offset:offset+c['batch_size']]
                fs, ps = [features[i] for i in indices], [priors[i] for i in indices]
                pk = pack(fs, ps, device)
                # Consume the same sampler stream in A, preserving the fixed CPU plan's RNG trajectory.
                As = adjacencies(ps,c['train_graph_samples'],train_rng,device)
                opt.zero_grad()
                value = batched_forward(model,pk,As)
                target = torch.as_tensor(y[indices],device=device)
                loss = (value-target).square().mean()  # ensemble prediction BEFORE the observation loss
                loss.backward()
                opt.step()
                steps += 1
                losses.append(float(loss.detach().cpu()))
            history.append({'phase':phase,'epoch':epoch+1,'mean_minibatch_mse':float(np.mean(losses))})
            if (epoch+1) % 5 == 0:
                print(f'RNA {phase} epoch {epoch+1}/{epochs} train_mse={np.mean(losses):.8g}',flush=True)
        save_stage('rna_training_'+phase, {'history':history, 'optimizer_steps':steps, 'model_sha256':state_hash(model)})
    raw, raw_se = predict(model,features,priors,64,c['seed']+30000,budget,device)
    A = preds['actual_A']
    grid = {str(g):float(np.mean((A[dv]+g*(raw[dv]-A[dv])-y[dv])**2)) for g in cfg['capability']['gamma_grid']}
    gamma = min(cfg['capability']['gamma_grid'],key=lambda g:(grid[str(g)],g))
    model.set_gamma(gamma)
    preds['actual_A_plus_B'] = A+gamma*(raw-A)
    save_model('rna_A_plus_B',model)
    ref = min(['guide_ridge','guide_pairwise_ridge','guide_window_ridge'],key=lambda arm:np.mean((preds[arm][dv]-y[dv])**2))
    save_stage('rna_predictions', {'predictions':{k:v.tolist() for k,v in preds.items()}, 'gamma':gamma,
                'gamma_grid':grid, 'reference_selected_on_development':ref, 'row_ids':[r['row_id'] for r in rows]})
    print(f'RNA fitting complete; dev-selected gamma={gamma}, reference={ref}',flush=True)
    comparisons = {arm:paired_comparison(rows,te,p,preds[ref],cfg['metrics']) for arm,p in preds.items()}
    sensitivity = {}
    for temp in c['prior_temperatures']:
        budget.check()
        if temp == 1.:
            p, err = preds['actual_A_plus_B'], gamma*raw_se
        else:
            ff, pp = build_data(rows,temp,budget)
            p, err = predict(model,ff,pp,64,c['seed']+30000,budget,device)
            del ff, pp
            gc.collect()
        sensitivity[str(temp)] = {'mean_prediction_se':float(np.mean(err)), 'predictions':p.tolist(),
              'prediction_sha256':array_hash(p), 'max_prediction_change':float(np.max(abs(p-preds['actual_A_plus_B']))),
              'ranking':paired_comparison(rows,te,p,preds['actual_A_plus_B'],cfg['metrics'])}
        print(f'RNA prior sensitivity T={temp} complete',flush=True)
    ff, pp = build_data(rows,1.,budget,alternative=True)
    alt, alt_se = predict(model,ff,pp,64,c['seed']+30000,budget,device)
    del ff, pp
    gc.collect()
    other, other_se = predict(model,features,priors,64,c['seed']+40000,budget,device)
    p16, se16 = predict(model,features,priors,16,c['seed']+30000,budget,device)
    identical = np.array([i for i in te if rows[i]['isoform_class']=='identical_local_windows'])
    # Full residual, even when development chooses gamma zero, permits inspection of MC convergence.
    selected_gamma = model.gamma
    model.set_gamma(1.)
    raw16, raw16se = predict(model,features,priors,16,c['seed']+30000,budget,device)
    model.set_gamma(selected_gamma)
    result = {'status':'COMPLETE_EXPLORATORY_CONDITIONAL_REPORTER', 'reference_selected_on_development':ref,
        'comparisons':comparisons,'ridge_fits':fits,'gamma_grid':grid,'gamma':gamma,
        'predictions':{k:v.tolist() for k,v in preds.items()},'row_ids':[r['row_id'] for r in rows],
        'labels':y.tolist(), 'split':split, 'training_history':history, 'optimizer_steps':steps,
        'all_arm_metrics':{arm:{part:metrics(p[ix],y[ix]) for part,ix in [('train',tr),('dev',dv),('test',te)]} for arm,p in preds.items()},
        'prior_sensitivity':sensitivity,'isoform_alternative':paired_comparison(rows,te,alt,preds['actual_A_plus_B'],cfg['metrics']),
        'isoform_alternative_predictions':alt.tolist(),
        'identical_window_only':paired_comparison(rows,identical,preds['actual_A_plus_B'],preds[ref],cfg['metrics']) if len(identical) else None,
        'sampling':{'mean_selected_gamma_se_16':float(se16.mean()),'mean_selected_gamma_se_64':float((gamma*raw_se).mean()),
                    'mean_raw_gamma1_se_16':float(raw16se.mean()),'mean_raw_gamma1_se_64':float(raw_se.mean()),
                    'max_selected_prediction_change_16_to_64':float(np.max(abs(p16-preds['actual_A_plus_B']))),
                    'max_raw_prediction_change_16_to_64':float(np.max(abs(raw16-raw))),
                    'independent_64_draw_max_prediction_change':float(np.max(abs(other-preds['actual_A_plus_B']))),
                    'draw16_predictions':p16.tolist(),'independent_draw64_predictions':other.tolist()},
        'gradient_sampling':'independent valid backtraces resampled each minibatch; gradient variance not estimated',
        'sampling_se_units':'predicted inhibition fraction, including gamma','model_sha256':state_hash(model),
        'chemical_features':None,'assay_identity':'YFP / H1299 / 48 h','raw_labels_clipped':False,
        'context_status':'inferred_contiguous_native_insert','biological_order_certification':False,
        'device':str(device),'readout_scale_modified':False}
    return result, model
