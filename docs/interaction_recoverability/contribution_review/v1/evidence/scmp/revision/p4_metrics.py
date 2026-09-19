"""Paired row-aligned, group-resampled and target-specific measurements."""
from __future__ import annotations
import numpy as np
from scipy.stats import rankdata


def metrics(pred, y, k=5, threshold=.7):
    pred, y = np.asarray(pred), np.asarray(y)
    var = float(np.var(y)) if len(y) else 0
    rho = None
    if len(y) > 1 and np.std(pred) and np.std(y):
        rho = float(np.corrcoef(rankdata(pred), rankdata(y))[0, 1])
    mse = float(np.mean((pred-y)**2)) if len(y) else None
    out = {'n': len(y), 'spearman': rho, 'mse': mse, 'r2': 1-mse/var if var else None,
           'threshold_prevalence': float(np.mean(y >= threshold)) if len(y) else None,
           'k': k, 'ndcg_at_k': None, 'precision_at_k': None,
           'selected_activity': None, 'selection_regret': None}
    if len(y) < k:
        out['ranking_status'] = 'candidate_pool_smaller_than_k'
        return out
    chosen, ideal = np.argsort(-pred, kind='stable')[:k], np.argsort(-y, kind='stable')[:k]
    discount = 1/np.log2(np.arange(k)+2)
    dcg, idcg = [np.sum(np.expm1(np.log(2)*y[ix])*discount) for ix in (chosen, ideal)]
    out.update(ndcg_at_k=float(dcg/idcg) if idcg > 0 else None,
               precision_at_k=float(np.mean(y[chosen] >= threshold)),
               selected_activity=float(y[chosen].mean()),
               selection_regret=float(y[ideal].mean()-y[chosen].mean()),
               selected_local_indices=chosen.tolist(), ranking_status='measured')
    return out


def paired_comparison(rows, indices, p, q, cfg):
    """Resample global sequence components jointly across models, never seeds."""
    yy = np.array([r['efficacy'] for r in rows])[indices]
    pp, qq = np.asarray(p)[indices], np.asarray(q)[indices]
    rr = [rows[i] for i in indices]
    genes = sorted({r['gene'] for r in rr})
    a, b = [metrics(v, yy, cfg['k'], cfg['activity_threshold']) for v in (pp, qq)]
    by_gene = {}
    for g in genes:
        ix = [i for i, r in enumerate(rr) if r['gene'] == g]
        by_gene[g] = {'proposed': metrics(pp[ix], yy[ix], cfg['k'], cfg['activity_threshold']),
                      'reference': metrics(qq[ix], yy[ix], cfg['k'], cfg['activity_threshold'])}
    gids = sorted({r['global_group'] for r in rr})
    lookup = {g: [i for i, r in enumerate(rr) if r['global_group'] == g] for g in gids}
    rng, deltas = np.random.RandomState(cfg['bootstrap_seed']), []
    for _ in range(cfg['bootstrap_groups']):
        ix = [i for g in rng.choice(gids, len(gids), replace=True) for i in lookup[g]]
        ap, bp = [metrics(v[ix], yy[ix])['spearman'] for v in (pp, qq)]
        if ap is not None and bp is not None:
            deltas.append(ap-bp)
    leave_out = {}
    for g in genes:
        ix = [i for i, r in enumerate(rr) if r['gene'] != g]
        ap, bp = [metrics(v[ix], yy[ix])['spearman'] for v in (pp, qq)]
        leave_out[g] = ap-bp if ap is not None and bp is not None else None
    return {'proposed': a, 'reference': b, 'per_target': by_gene,
            'delta_spearman': a['spearman']-b['spearman'] if a['spearman'] is not None and b['spearman'] is not None else None,
            'group_paired_ci95': np.quantile(deltas, [.025,.975]).tolist() if deltas else None,
            'resampling_groups': len(gids), 'leave_one_target_out': leave_out,
            'interpretation': 'exploratory within observed targets; not new-target inference'}
