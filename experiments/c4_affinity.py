"""C4: is the certified network affine?

With X >= 0, W > 0, b >= 0, every preactivation is non-negative, so ReLU never clips and
a linear aggregator collapses the network to an exactly AFFINE function of A. Test by
additivity: f(X1+X2) - f(0) == [f(X1)-f(0)] + [f(X2)-f(0)].

This matters three ways and all three belong in the paper:
  - it is an expressivity limitation of the ReLU instantiation a reviewer will find;
  - it explains the encoding result: compressed features starve a linear model;
  - it is what makes the slack scaling law of C5 explicable.
It is NOT a limitation of the certified class: tanh and sigmoid satisfy H3 and are
nonlinear on the positive orthant. max is nonlinear regardless.

POST-HOC. Not pre-registered.
"""
import numpy as np
from certmp.models import MonoMPNN
from certmp.provenance import new_run, record

TRIALS, N, D = 30, 6, 3


def violation(agg, act, seed):
    r = np.random.RandomState(seed)
    A = (r.rand(N, N) < .5).astype(float); A = np.maximum(A, A.T); np.fill_diagonal(A, 1.)
    m = MonoMPNN(D, 8, 2, agg=agg, nonneg=True, act=act, seed=seed)
    X1, X2 = r.rand(N, D), r.rand(N, D)
    f0 = m.forward(np.zeros((N, D)), A)
    lhs = m.forward(X1 + X2, A) - f0
    rhs = (m.forward(X1, A) - f0) + (m.forward(X2, A) - f0)
    return abs(lhs - rhs) / max(abs(lhs), 1e-12), m.is_affine()


def main():
    run = new_run("c4_affinity", dict(trials=TRIALS, n=N, d_in=D))
    rows = []
    print(f"{'aggregation':>12s} {'activation':>11s} {'max rel violation':>19s}  {'measured':>10s}  predicted")
    for agg in ("sum", "max", "logsumexp", "min"):
        for act in ("relu", "tanh", "sigmoid"):
            vs = [violation(agg, act, s) for s in range(TRIALS)]
            worst = max(v for v, _ in vs); pred = vs[0][1]
            measured_affine = worst < 1e-9
            ok = measured_affine == pred
            print(f"{agg:>12s} {act:>11s} {worst:19.3e}  "
                  f"{'AFFINE' if measured_affine else 'nonlinear':>10s}  "
                  f"{'affine' if pred else 'nonlinear':>9s} {'ok' if ok else '*** MISMATCH ***'}")
            rows.append(dict(agg=agg, act=act, max_rel_violation=worst,
                             measured_affine=measured_affine, predicted_affine=pred, agrees=ok))
    print("\n  sum+relu is exactly affine: ReLU never clips on the non-negative orthant.")
    print("  tanh and sigmoid satisfy H3 and restore nonlinearity. max is nonlinear throughout.")
    record(run, "c4_affinity", {"rows": rows, "all_agree": all(r["agrees"] for r in rows)})


if __name__ == "__main__": main()
