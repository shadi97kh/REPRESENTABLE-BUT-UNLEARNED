"""C6: closure of the certified class under common architectural additions.

A theorem restricted to one vanilla MPNN equation invites the objection that real
architectures do not look like that. Test whether endpoint exactness survives:
  - non-negative residual connections   h <- h + act(agg(...))
  - deeper stacks (1 to 5 layers)
  - each monotone activation in H3
  - both certifiable aggregation directions, isotone and antitone
and confirm it still FAILS for a non-monotone activation, so the test can discriminate.

Voids, where the output is constant across the whole lattice, are reported separately and
excluded rather than counted as confirmations.

POST-HOC. Not pre-registered.
"""
import itertools
import numpy as np
from certmp.models import MonoMPNN
from certmp.reach import build_A
from certmp.provenance import new_run, record

N, K, TRIALS, D = 6, 7, 15, 8      # d_in == d_hid so residual dims match


def exact_frac(agg, act, residual, L, seeds=TRIALS):
    n_exact = n_live = 0
    for s in range(seeds):
        r = np.random.RandomState(s)
        mand = [(0, 1), (1, 2)]
        cand = [(i, j) for i in range(N) for j in range(i + 1, N) if (i, j) not in mand]
        opt = [cand[t] for t in r.choice(len(cand), size=min(K, len(cand)), replace=False)]
        X = r.rand(N, D) * 2.0
        m = MonoMPNN(D, D, L, agg=agg, nonneg=True, act=act, residual=residual, seed=s)
        ys = {b: m.forward(X, build_A(N, mand, opt, b))
              for b in itertools.product([0, 1], repeat=len(opt))}
        vals = list(ys.values()); ymax, ymin = max(vals), min(vals)
        sc = max(abs(ymax), abs(ymin), 1e-12)
        if (ymax - ymin) <= 1e-12 * sc: continue
        n_live += 1
        eps = {ys[tuple([0] * len(opt))], ys[tuple([1] * len(opt))]}
        if (min(abs(ymax - e) for e in eps) <= 1e-12 * sc and
                min(abs(ymin - e) for e in eps) <= 1e-12 * sc): n_exact += 1
    return n_exact, n_live


def main():
    run = new_run("c6_closure", dict(n=N, k=K, trials=TRIALS, d=D))
    rows = []
    print("closure under residual connections, depth and monotone activations")
    print(f"{'agg':>10s} {'act':>9s} {'residual':>9s} {'L':>3s} {'exact/live':>13s}  verdict")
    for agg in ("sum", "max", "min"):
        for act in ("relu", "tanh", "sigmoid"):
            for residual in (False, True):
                for L in (1, 3, 5):
                    ne, nl = exact_frac(agg, act, residual, L)
                    void = (nl == 0)
                    ok = (ne == nl) and not void
                    if void or not ok:
                        tag = "VOID (all dead)" if void else "*** VIOLATION ***"
                        print(f"{agg:>10s} {act:>9s} {str(residual):>9s} {L:3d} "
                              f"{ne:5d}/{nl:<7d}  {tag}")
                    rows.append(dict(agg=agg, act=act, residual=residual, layers=L,
                                     n_exact=ne, n_live=nl, void=void, passed=ok))
    npass = sum(r["passed"] for r in rows)
    nvoid = sum(r["void"] for r in rows)
    viol = [r for r in rows if not r["passed"] and not r["void"]]
    print(f"\n  {npass}/{len(rows) - nvoid} live configurations retain endpoint exactness")
    print(f"  {nvoid} configurations VOID (output constant across the lattice), excluded")
    print(f"  {len(viol)} genuine violations")

    print("\n  discriminative control: a NON-monotone activation must BREAK exactness")
    ctrl = []
    for agg in ("sum", "max"):
        for residual in (False, True):
            ne, nl = exact_frac(agg, "sin", residual, 3)
            broke = ne < nl
            ctrl.append(dict(agg=agg, act="sin", residual=residual, n_exact=ne, n_live=nl,
                             broke_as_expected=broke))
            print(f"{agg:>10s} {'sin':>9s} {str(residual):>9s} {3:3d} {ne:5d}/{nl:<7d}  "
                  f"{'PASS (breaks)' if broke else '*** did not break ***'}")
    record(run, "c6_closure", {"rows": rows, "controls": ctrl, "n_void": nvoid,
                               "violations": viol, "closure_holds": len(viol) == 0,
                               "control_discriminates": all(c["broke_as_expected"] for c in ctrl)})


if __name__ == "__main__": main()
