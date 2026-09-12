"""Frozen-predictor comparison with atomic anytime updates and a skip action.

The learned scalar chooses whether to condition. It never supplies a bound.
All numerical results remain float_unverified. No model parameters are trained
here. The caller must first pass the independent predictive capability gate.
"""
from __future__ import annotations
import copy
import itertools
import time
import numpy as np
import torch
from ..bounds import certified_upper_bound
from ..conditional import conditional_upper_bound
from .p4_capability import make_inputs, state_hash, fit_ridge
from .p4_resources import Budget, BudgetStop


def decision(lo, hi, tau):
    return bool(lo >= tau or hi < tau)


def policy_features(cheap, n_edges, tau):
    s = cheap.stats
    unstable = s.get('relu_unstable', 0)
    stable = s.get('relu_stable_pos', 0)+s.get('relu_stable_neg', 0)
    return np.array([cheap.gap, n_edges, unstable/max(1, stable+unstable),
                     tau-cheap.incumbent_value, cheap.upper-tau], dtype=float)


def run_one(frozen, cfg, case, arm, cap, policy=None):
    """Measure setup, inference, queries and completed updates inside one cap."""
    start = time.perf_counter()
    events, action, features = [], 'no_update', None
    counters = {'support': 0, 'count': 0, 'feasibility': 0, 'marginals': 0}
    times = {'setup_s': 0., 'policy_s': 0., 'bound_s': 0., 'enumeration_s': 0.}
    lo, hi = -float('inf'), float('inf')
    timed_out = False
    try:
        with Budget(cap, cap).enforce():
            setup_start = time.perf_counter()
            m, X, oracle, B, members, Z, As = make_inputs(case['sequence'], cfg['capability']['hidden'], 0)
            m.load_state_dict(frozen.state_dict())
            m.set_gamma(frozen.gamma)
            for p in m.parameters():
                p.requires_grad_(False)
            # Constrain only valid members; graph coordinates stay in this sequence.
            forbidden = tuple(tuple(e) for e in case['forbidden'])
            members = [s for s in members if not set(s) & set(forbidden)]
            for name in counters:
                orig = getattr(oracle, name)
                def wrapped(*a, _name=name, _fn=orig, **kw):
                    counters[_name] += 1
                    return _fn(*a, **kw)
                setattr(oracle, name, wrapped)
            times['setup_s'] = time.perf_counter()-setup_start
            t = time.perf_counter()
            cheap = certified_upper_bound(m, X, oracle, B, forbidden=forbidden)
            times['bound_s'] += time.perf_counter()-t
            lo, hi = float(cheap.incumbent_value), float(cheap.upper)
            events.append({'time_s': time.perf_counter()-start, 'lo': lo, 'hi': hi})
            features = policy_features(cheap, len(m.candidates)-len(forbidden), case['tau'])
            if decision(lo, hi, case['tau']):
                action = 'skip_already_decided'
            elif arm == 'combined':
                action = 'skip_fixed_no_conditioning'
            elif arm == 'exhaustive':
                action = 'enumerate'
                t = time.perf_counter()
                for s in members:
                    with torch.no_grad():
                        lo = max(lo, float(m(X, m.indicator(s), B)))
                    events.append({'time_s': time.perf_counter()-start, 'lo': lo, 'hi': hi})
                    if lo >= case['tau']:
                        break
                else:
                    hi = lo
                    events.append({'time_s': time.perf_counter()-start, 'lo': lo, 'hi': hi})
                times['enumeration_s'] = time.perf_counter()-t
            else:
                cond = arm == 'fixed_conditional'
                if arm == 'deterministic_selective':
                    cond = features[2] > 1-cfg['verifier']['skip_stable_fraction']
                elif arm == 'learned_selective':
                    t = time.perf_counter()
                    cond = float(np.r_[1., features]@policy) > 0
                    times['policy_s'] = time.perf_counter()-t
                if cond:
                    action = 'condition'
                    t = time.perf_counter()
                    r = conditional_upper_bound(m, X, oracle, B, forbidden=forbidden,
                                                budget=cfg['verifier']['condition_budget'], scorer='largest_gap')
                    times['bound_s'] += time.perf_counter()-t
                    # Incumbents and uppers improve independently; retain both.
                    lo, hi = max(lo, float(r.incumbent_value)), min(hi, float(r.upper))
                    events.append({'time_s': time.perf_counter()-start, 'lo': lo, 'hi': hi})
                else:
                    action = 'skip_policy' if arm == 'learned_selective' else 'skip_few_unstable'
    except BudgetStop:
        timed_out = True
    elapsed = time.perf_counter()-start
    accepted = [e for e in events if e['time_s'] <= cap]
    if accepted:
        lo, hi = accepted[-1]['lo'], accepted[-1]['hi']
    else:
        lo, hi = -float('inf'), float('inf')
    resolved = decision(lo, hi, case['tau'])
    return {'arm': arm, 'case': case, 'cap_s': cap, 'wall_s': elapsed,
            'lower': lo if np.isfinite(lo) else None, 'upper': hi if np.isfinite(hi) else None,
            'gap': hi-lo if np.isfinite(hi-lo) else None,
            'decision_coverage_indicator': resolved, 'timeout': timed_out,
            'outcome': 'decided' if resolved else 'unresolved',
            'action': action, 'oracle_calls': counters, 'component_times': times,
            'policy_features': features.tolist() if features is not None else None,
            'events': accepted, 'numerical_status': 'float_unverified',
            'negative_gap_defect': bool(hi < lo)}


def cases(cfg, n, seed):
    out = []
    for stratum, seq in [('original_family', cfg['capability']['family']),
                          ('wider_family', cfg['verifier']['wider_family'])]:
        _, _, o, _, _, _, _ = make_inputs(seq, cfg['capability']['hidden'], 0)
        edges = o.ground_set()
        masks = [()] + [(e,) for e in edges] + list(itertools.combinations(edges, 2))
        catalog = [(mask, tau) for mask in masks for tau in cfg['verifier']['thresholds']]
        rng = np.random.RandomState(cfg['verifier']['input_seeds'][0])
        order = rng.permutation(len(catalog))
        offset = 0 if seed == cfg['verifier']['input_seeds'][0] else cfg['verifier']['development_instances']
        if offset+n > len(catalog):
            raise ValueError('not enough distinct input/threshold groups')
        for idx in order[offset:offset+n]:
            mask, tau = catalog[idx]
            out.append({'id': f'{stratum}:{idx}', 'stratum': stratum, 'sequence': seq,
                        'forbidden': [list(e) for e in mask], 'tau': tau})
    return out


def compare_at_total_cost(rows, offline_s, cfg):
    results = []
    for cap in cfg['budgets_s']:
        for n in cfg['offline_amortization_instances']:
            coverage = {}
            for arm in cfg['arms']:
                rr = [r for r in rows if r['arm'] == arm and r['cap_s'] == cap and r['case']['stratum'] == 'original_family']
                effective = cap - (offline_s/n if arm == 'learned_selective' else 0.)
                resolved = 0
                for r in rr:
                    valid = [e for e in r['events'] if e['time_s'] <= effective]
                    resolved += bool(valid and decision(valid[-1]['lo'], valid[-1]['hi'], r['case']['tau']))
                coverage[arm] = resolved/len(rr) if rr else 0.
            strongest = max(v for a,v in coverage.items() if a != 'learned_selective')
            results.append({'budget_s':cap,'amortization_instances':n,'coverage':coverage,
                            'gain_over_strongest':coverage['learned_selective']-strongest,
                            'passes_minimum':coverage['learned_selective']-strongest >= cfg['min_coverage_gain']})
    return results


def run_verifier(model, cfg, budget):
    before = state_hash(model)
    c = cfg['verifier']
    development = cases(cfg, c['development_instances'], c['input_seeds'][0])
    heldout = cases(cfg, c['test_instances'], c['input_seeds'][1])
    teacher, feats, targets = [], [], []
    offline_start = time.perf_counter()
    for case in development:
        budget.check()
        pair = [run_one(model, cfg, case, a, max(c['budgets_s']))
                for a in ['combined', 'fixed_conditional']]
        teacher.append(pair)
        a, b = pair
        if a['policy_features'] is not None:
            # Downstream decision gain per TOTAL measured cost; never gap-only.
            gain = (int(b['decision_coverage_indicator'])-int(a['decision_coverage_indicator']))/max(b['wall_s'], 1e-9)
            cost = max(0., b['wall_s']-a['wall_s'])/max(b['wall_s'], 1e-9)
            feats.append(a['policy_features'])
            targets.append(gain-cost)
    if not feats:
        return {'status': 'NO_POLICY_TRAINING_FEATURES_WITHIN_CAP', 'useful_learned_gain': False}
    _, policy = fit_ridge(np.array(feats), np.array(targets), np.arange(len(feats)), 1.)
    offline = time.perf_counter()-offline_start
    rows, exact_audits = [], []
    for case in heldout:
        budget.check()
        audit_start = time.perf_counter()
        mm, XX, oo, BB, members, ZZ, AA = make_inputs(case['sequence'],cfg['capability']['hidden'],0)
        mm.load_state_dict(model.state_dict())
        mm.set_gamma(model.gamma)
        forbidden = {tuple(e) for e in case['forbidden']}
        with torch.no_grad():
            exact = max(float(mm(XX,mm.indicator(member),BB)) for member in members if not forbidden.intersection(member))
        exact_audits.append({'case_id':case['id'],'exact_float_optimum':exact,'audit_wall_s':time.perf_counter()-audit_start})
        for cap in c['budgets_s']:
            for arm in c['arms']:
                budget.check()
                r = run_one(model, cfg, case, arm, cap, policy)
                r['exact_float_optimum'] = exact
                r['upper_minus_exact'] = r['upper']-exact if r['upper'] is not None else None
                r['sound_vs_enumerated_float_tolerance_1e_10'] = r['upper'] >= exact-1e-10 if r['upper'] is not None else None
                rows.append(r)
    assert before == state_hash(model), 'verifier modified predictor'
    summary = []
    for cap in c['budgets_s']:
        for stratum in ['all', 'original_family', 'wider_family']:
            for arm in c['arms']:
                rr = [r for r in rows if r['arm'] == arm and r['cap_s'] == cap and
                      (stratum == 'all' or r['case']['stratum'] == stratum)]
                gaps = [r['gap'] for r in rr if r['gap'] is not None]
                summary.append({'arm': arm, 'stratum': stratum, 'budget_s': cap, 'n': len(rr),
                                'decision_coverage': sum(r['decision_coverage_indicator'] for r in rr)/len(rr),
                                'mean_gap_when_available': float(np.mean(gaps)) if gaps else None,
                                'gap_unavailable_n': len(rr)-len(gaps),
                                'timeouts': sum(r['timeout'] for r in rr),
                                'wall_quantiles': np.quantile([r['wall_s'] for r in rr], [.5,.9,1]).tolist()})
    return {'status': 'COMPLETE_FLOAT_DIAGNOSTIC', 'model_sha256': before, 'rows': rows,
            'summary': summary, 'offline_teacher_and_fit_wall_s': offline,
            'exact_optimum_audits': exact_audits,
            'audit_cost_contract': 'separate ground-truth audit, charged to cycle; each compared arm still pays its own setup and oracle work',
            'amortized_cost_per_instance_s': {str(n): offline/n for n in c['offline_amortization_instances']},
            'teacher': teacher, 'policy_coefficients': policy.tolist(),
            'useful_learned_gain': any(r['passes_minimum'] for r in compare_at_total_cost(rows,offline,c) if r['amortization_instances']==100),
            'equal_total_cost': compare_at_total_cost(rows,offline,c),
            'gain_status': 'primary amortization 100 instances; original family only; wider-family predictive utility not established',
            'numerical_status': 'float_unverified',
            'closest_external_baseline': 'exhaustive valid-member evaluation on these tiny families',
            'scope': 'model threshold decisions, not measured biological order'}
