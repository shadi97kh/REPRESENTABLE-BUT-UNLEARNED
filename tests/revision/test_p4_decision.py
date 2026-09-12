"""No-training tests for P4 scientific contracts and resource admission."""
import signal
import time
from pathlib import Path
import numpy as np
import pytest
import torch
import yaml

from scmp.revision.p4_capability import (control_data, ensemble_forward, make_inputs,
                                        state_hash)
from scmp.revision.p4_data import components, partition
from scmp.revision.p4_metrics import metrics, paired_comparison
from scmp.revision.p4_resources import Budget, BudgetStop, checked_record, inventory, write_record


@pytest.fixture
def cfg():
    return yaml.safe_load(Path('configs/revision_p4.yaml').read_text())


def test_budget_is_not_inferred_from_wall_proxy(cfg):
    ledger = inventory(cfg)
    assert ledger['authorization']['cpu_core_hours'] == 8
    assert ledger['counterexample_16_core_hours_pilots_plus_p3'] > 8
    assert not ledger['admitted']
    assert ledger['verified_spendable_cpu_seconds'] == 0
    assert ledger['balance_status'] == 'unverified_not_proven_exhausted'


def test_nested_timer_cannot_extend_enclosing_cap():
    began = time.monotonic()
    enclosing = signal.getitimer(signal.ITIMER_REAL)[0]
    handler = signal.getsignal(signal.SIGALRM)
    with pytest.raises(BudgetStop):
        with Budget(.02, .02).enforce():
            with Budget(1, 1).enforce():
                time.sleep(.1)
    assert time.monotonic()-began < .2
    remaining = signal.getitimer(signal.ITIMER_REAL)[0]
    assert signal.getsignal(signal.SIGALRM) == handler
    assert (remaining == 0 if not enclosing else 0 < remaining < enclosing)


def test_immutable_record_and_bad_sidecar(tmp_path):
    p = tmp_path/'record.json'
    write_record(p, {'preserved': True})
    assert checked_record(p)['preserved']
    with pytest.raises(FileExistsError):
        write_record(p, {})
    p.write_text('{}')
    with pytest.raises(ValueError):
        checked_record(p)


def test_distribution_observations_are_identifiable_and_paired(cfg):
    *_, info = control_data(cfg['capability'])
    assert info['equal_marginal_pair_max_error'] < 1e-14
    assert info['equal_marginal_label_difference_min'] > .1
    assert info['target_pairwise_projection']['max_abs'] > .1
    partitions = [set(info[k]) for k in ['train_ids','dev_ids','test_ids']]
    assert not any(a & b for i,a in enumerate(partitions) for b in partitions[i+1:])
    for i in range(0,96,2):
        assert any(i in s and i+1 in s for s in partitions)


def test_ensemble_mean_before_loss_and_exact_bypass(cfg):
    m, X, o, B, members, Z, As, P, y, *_ = control_data(cfg['capability'])
    before = state_hash(m)
    w = torch.as_tensor(P[0])
    vals = torch.stack([m(X,z,B) for z in Z])
    got = ensemble_forward(m,X,w@Z,Z,As,w)
    assert torch.allclose(got,w@vals,atol=1e-12,rtol=0)
    loss_of_mean, mean_loss = (got-y[0])**2, w@((vals-y[0])**2)
    assert float((mean_loss-loss_of_mean).detach()) > 0
    loss_of_mean.backward()
    assert any(p.grad is not None and torch.count_nonzero(p.grad) for p in m.g.parameters())
    assert before == state_hash(m)
    m.set_gamma(0)
    with torch.no_grad():
        m.g.readout.weight.fill_(float('nan'))
    b,a = m.anchor_terms(X)
    assert torch.equal(ensemble_forward(m,X,w@Z,Z,As,w), b+(w@Z)@a)


def test_global_shared_kmer_is_transitive_without_gene_key():
    seqs = ['AAAACCCC', 'CCCCGGGG', 'GGGGUUUU', 'ACGUACGU']
    g = components(seqs,k=4)
    assert g[0] == g[1] == g[2] and g[3] != g[0]
    folds = partition(g,3,.5,0)
    assert not ({g[i] for i in folds[0]} & {g[i] for i in folds[2]})


def test_ranking_uses_unclipped_activity_and_small_pool_status():
    y = np.array([1.341,1.1,.7,.2,.1])
    m = metrics(y,y,k=5)
    assert m['ndcg_at_k'] == 1
    assert m['selection_regret'] == 0
    assert m['threshold_prevalence'] == .6
    assert metrics(y[:2],y[:2])['precision_at_k'] is None
    assert metrics(y,y,k=1)['selected_activity'] == 1.341


def test_paired_bootstrap_keeps_same_rows_and_groups():
    rows = [{'efficacy':float(i)/10,'gene':str(i//4),'global_group':i//2} for i in range(12)]
    p = np.arange(12,dtype=float)
    c = dict(k=2,activity_threshold=.7,bootstrap_seed=2,bootstrap_groups=20)
    out = paired_comparison(rows,np.arange(12),p,p,c)
    assert out['group_paired_ci95'] == [0.,0.]
    assert out['delta_spearman'] == 0
    assert out['resampling_groups'] == 6


def test_rna_receives_valid_distinct_graphs_and_shared_weights():
    from scmp.revision.p4_rna import window_input
    row = {'guide':'A'*21}
    m,X,mu,Z,As,p = window_input(row,'GCGCAAAAGCGC',1.,16,41)
    assert len(torch.unique(Z,dim=0)) > 1
    m.set_gamma(1.)
    n,Y,nu,U,adj,q = window_input(row,'GGGGAAAACCCC',1.,16,42,m)
    assert n.anchor is m.anchor and n.g is m.g
    assert torch.isfinite(ensemble_forward(n,Y,nu,U,adj))


def test_verifier_does_not_mutate_predictor_and_records_timeout(cfg):
    from scmp.revision.p4_verifier import run_one
    m,*_ = make_inputs(cfg['capability']['family'],cfg['capability']['hidden'],1)
    before = state_hash(m)
    case = dict(sequence=cfg['capability']['family'],forbidden=[],tau=1e10,id='smoke',stratum='original_family')
    r = run_one(m,cfg,case,'deterministic_selective',1.)
    assert state_hash(m) == before
    assert r['action'] == 'skip_already_decided'
    assert r['numerical_status'] == 'float_unverified'
    t = run_one(m,cfg,case,'fixed_conditional',.000001)
    assert t['timeout'] and t['outcome'] == 'unresolved'
    assert t['events'] == []
