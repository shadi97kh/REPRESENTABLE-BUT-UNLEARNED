"""C7: the worst-case relaxation gap is a walk-count ratio.

Claim: for the certified class, gap = kappa_L exactly when biases vanish, and biases only
dilute, so kappa_L is a tight upper bound. Five families spanning the full range from
maximum exclusion (matchings) to none (join-closed), by EXHAUSTIVE enumeration so that
max over the feasible family is exact rather than sampled.

Registered expectations:
  K1 no-bias gap == kappa_L to machine precision on every family and depth
  K2 biased gap <= kappa_L always (tightness, never violated)
  K3 join-closed family gives kappa_L == 1.0 exactly at every depth (exactness criterion)
  K4 kappa_L is non-increasing along the NESTED chain
     matching subset deg<=2 subset deg<=3 subset all-subsets.
     Forests are deliberately excluded from the chain: a forest is neither a subset nor a
     superset of a degree-bounded family, so its kappa need not sit between theirs. An
     earlier version of this check assumed a total order and failed for exactly that
     reason; the corrected claim is monotonicity along nesting, not along a list.
"""
import numpy as np
from certmp.models import MonoMPNN
from certmp.kappa import (adjacency, kappa_L, measured_gap, enumerate_family, FAMILIES)
from certmp.provenance import new_run, record

N, N_CAND, SEEDS = 7, 12, 8


def zero_bias(model):
    """Return a copy with biases zeroed: the regime where kappa_L is attained exactly."""
    import copy
    m = copy.deepcopy(model)
    m.B = [np.zeros_like(b) for b in m.B]
    return m


def main():
    run = new_run("c7_kappa_characterization",
                  dict(n=N, n_candidate=N_CAND, seeds=SEEDS, depths=[1, 2, 3]))
    cand = [(i, j) for i in range(N) for j in range(i + 1, N)][:N_CAND]
    X = np.ones((N, 2))
    A_top = adjacency(N, cand)
    rows = []
    print(f"n={N}, |candidate|={len(cand)}, enumerating 2^{len(cand)} subsets exactly\n")
    print(f"{'family':>28s} {'L':>2s} {'kappa_L':>12s} {'gap(no bias)':>13s} "
          f"{'gap(bias)':>11s} {'rel err':>10s} {'tight':>6s}")
    for fname, pred in FAMILIES:
        feas = enumerate_family(N, cand, pred)
        A_feas = [adjacency(N, g) for g in feas]
        for L in (1, 2, 3):
            k, wtop, wbest = kappa_L(N, cand, feas, L)
            gnb, gb = [], []
            for s in range(SEEDS):
                m = MonoMPNN(2, 6, L, agg="sum", nonneg=True, act="relu", seed=s)
                gb.append(measured_gap(m, X, A_top, A_feas)[0])
                gnb.append(measured_gap(zero_bias(m), X, A_top, A_feas)[0])
            gnb_m, gb_m = float(np.mean(gnb)), float(np.mean(gb))
            err = abs(gnb_m - k) / k
            tight = gb_m <= k + 1e-9
            print(f"{fname:>28s} {L:2d} {k:12.6f} {gnb_m:13.6f} {gb_m:11.6f} "
                  f"{err:10.2e} {str(tight):>6s}")
            rows.append(dict(family=fname, L=L, kappa=k, w_top=wtop, w_best=wbest,
                             gap_nobias=gnb_m, gap_bias=gb_m, rel_err=err,
                             tight=bool(tight), n_feasible=len(feas)))
    k1 = max(r["rel_err"] for r in rows) < 1e-12
    k2 = all(r["tight"] for r in rows)
    k3 = all(abs(r["kappa"] - 1.0) < 1e-12 for r in rows if "join-closed" in r["family"])
    CHAIN = ["matching (deg<=1)", "deg<=2", "deg<=3", "all subsets (join-closed)"]
    by_L = {}
    for r in rows:
        by_L.setdefault(r["L"], {})[r["family"]] = r["kappa"]
    k4 = all(all(d[CHAIN[i]] >= d[CHAIN[i + 1]] - 1e-12 for i in range(len(CHAIN) - 1))
             for d in by_L.values())
    chain_vals = {L: [round(by_L[L][c], 4) for c in CHAIN] for L in sorted(by_L)}
    print(f"\n  K1 no-bias gap == kappa_L      {'PASS' if k1 else 'FAIL'}"
          f"   (max rel err {max(r['rel_err'] for r in rows):.2e})")
    print(f"  K2 kappa_L is a tight bound    {'PASS' if k2 else 'FAIL'}")
    print(f"  K3 join-closed -> kappa_L = 1  {'PASS' if k3 else 'FAIL'}")
    print(f"  K4 monotone along nested chain {'PASS' if k4 else 'FAIL'}")
    for L, v in chain_vals.items():
        print(f"     L={L}: " + " >= ".join(f"{x:.3f}" for x in v))
    print("     (forests are incomparable to degree bounds and are not in the chain)")
    record(run, "c7_kappa_characterization",
           {"rows": rows, "K1": k1, "K2": k2, "K3": k3, "K4": k4,
            "nested_chain": CHAIN, "chain_values": chain_vals})


if __name__ == "__main__":
    main()
