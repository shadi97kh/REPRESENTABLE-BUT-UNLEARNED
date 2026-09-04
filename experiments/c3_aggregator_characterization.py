"""C3: the aggregator characterisation, generalising the original four-row table.

Claim: endpoint exactness holds iff the aggregator is monotone under MULTISET INCLUSION,
in EITHER direction. Isotone gives min at the empty optional set and max at the full one;
antitone swaps the two. Both are certifiable in two passes. Neither direction means an
interior optimum is reachable and two passes are unsound.

This predicts membership for aggregators absent from the original table: logsumexp is
isotone and certifiable, min is antitone and certifiable, std is neither and is not.

POST-HOC. Not pre-registered; see the corrections amendment in PREREGISTRATION.md.
"""
import itertools
import numpy as np
from certmp.models import MonoMPNN
from certmp.aggregators import AGGREGATORS
from certmp.reach import build_A
from certmp.provenance import new_run, record

N, K, TRIALS, D = 7, 8, 25, 3


def trial(agg, seed):
    r = np.random.RandomState(seed)
    mand = [(0, 1), (1, 2), (2, 3)]
    cand = [(i, j) for i in range(N) for j in range(i + 1, N) if (i, j) not in mand]
    opt = [cand[t] for t in r.choice(len(cand), size=K, replace=False)]
    X = r.rand(N, D) * 2.0
    m = MonoMPNN(D, 8, 2, agg=agg, nonneg=True, act="relu", seed=seed)
    ys = {bits: m.forward(X, build_A(N, mand, opt, bits))
          for bits in itertools.product([0, 1], repeat=K)}
    empty, full = tuple([0] * K), tuple([1] * K)
    vals = list(ys.values()); ymax, ymin = max(vals), min(vals)
    sc = max(abs(ymax), abs(ymin), 1e-12)
    live = (ymax - ymin) > 1e-12 * sc
    endpoints = {ys[empty], ys[full]}
    exact = (min(abs(ymax - e) for e in endpoints) <= 1e-12 * sc and
             min(abs(ymin - e) for e in endpoints) <= 1e-12 * sc)
    # which endpoint carries the max, to confirm the predicted direction empirically
    at_full = abs(ymax - ys[full]) <= 1e-12 * sc
    return exact, live, at_full


def main():
    run = new_run("c3_aggregator_characterization", dict(n=N, k=K, trials=TRIALS))
    rows = []
    print(f"{'aggregator':>12s} {'predicted':>10s} {'exact':>10s} {'live':>6s} {'max at':>8s}  verdict")
    for agg, (_, pred) in AGGREGATORS.items():
        rs = [trial(agg, s) for s in range(TRIALS)]
        live = [(e, af) for e, l, af in rs if l]
        n_exact = sum(e for e, _ in live)
        should = pred in ("isotone", "antitone")
        got = (n_exact == len(live)) and len(live) > 0
        ok = got == should
        at = "full" if live and all(af for _, af in live) else (
             "empty" if live and not any(af for _, af in live) else "mixed")
        dir_ok = (not got) or (at == "full" and pred == "isotone") or \
                 (at == "empty" and pred == "antitone")
        print(f"{agg:>12s} {pred:>10s} {n_exact:4d}/{len(live):<5d} {len(live):6d} {at:>8s}  "
              f"{'PASS' if ok and dir_ok else '*** MISMATCH ***'}")
        rows.append(dict(agg=agg, predicted=pred, certifiable_predicted=should,
                         n_exact=n_exact, n_live=len(live), observed_exact=got,
                         max_endpoint=at, direction_agrees=dir_ok, agrees=ok and dir_ok))
    allok = all(r["agrees"] for r in rows)
    print(f"\n  characterisation {'holds' if allok else 'FAILS'} on all {len(rows)} aggregators")
    print("  logsumexp and min are certifiable and appear in no prior table here;")
    print("  min is ANTITONE, so its endpoint roles are swapped, not absent.")
    record(run, "c3_aggregator_characterization",
           {"rows": rows, "characterization_holds": allok})


if __name__ == "__main__": main()
