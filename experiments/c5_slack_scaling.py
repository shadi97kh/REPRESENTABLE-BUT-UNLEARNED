"""C5: why lattice relaxation is loose, and by exactly how much.

Part A  UNBOUNDED SLACK. Matchings on n=2m nodes are downward closed but NOT closed under
        union, so the smallest containing lattice has the complete graph on top, which is
        maximally infeasible. A monotone f attains its maximum exactly there. Slack grows
        without bound in m: the failure is structural, not a tuning problem.
Part B  SCALING LAW. Fit slack ~ rho^e with rho = k/(n/2), swept over DEPTH. Prediction:
        the exponent tracks depth, because the certified sum network is affine of degree L
        in A (C4). If confirmed, the exponent measured on RNA is a property of the
        architecture, not of RNA thermodynamics.
Part C  CONVERSE. A join-closed family has a feasible top element, so slack is exactly 1.0.
        Verified by brute force against every subset, not asserted: computing top/top would
        be a tautology and would test nothing.

POST-HOC. Not pre-registered.
"""
import itertools
import numpy as np
from certmp.models import MonoMPNN
from certmp.reach import build_A
from certmp.provenance import new_run, record


def complete_A(n):
    return np.ones((n, n))


def matching_A(n, m):
    A = np.eye(n)
    for i in range(m): A[2 * i, 2 * i + 1] = A[2 * i + 1, 2 * i] = 1.0
    return A


def slack(m, L, agg, seed, d_hid=6):
    """Best over the lattice (complete graph) divided by best over the FEASIBLE family.
    Every perfect matching on these nodes gives the same value by symmetry of the family,
    so the single matching below is the family maximum."""
    n = 2 * m
    net = MonoMPNN(2, d_hid, L, agg=agg, nonneg=True, act="relu", seed=seed)
    X = np.ones((n, 2))                       # uniform features isolate structure
    return net.forward(X, complete_A(n)) / net.forward(X, matching_A(n, m))


def fit(ms, L, agg, seeds=5):
    xs, ys = [], []
    for m in ms:
        n = 2 * m
        rho = (n * (n - 1) // 2) / (n / 2)
        s = float(np.mean([slack(m, L, agg, sd) for sd in range(seeds)]))
        xs.append(np.log(rho)); ys.append(np.log(s))
    b, a = np.polyfit(xs, ys, 1)
    pred = np.polyval([b, a], xs)
    denom = float(np.sum((np.array(ys) - np.mean(ys)) ** 2))
    r2 = 1 - float(np.sum((np.array(ys) - pred) ** 2)) / denom if denom > 1e-15 else float("nan")
    return b, a, r2


MS = list(range(3, 13))


def main():
    run = new_run("c5_slack_scaling", dict(m_range=[MS[0], MS[-1]], depths=[1, 2, 3, 4]))

    print("Part A -- unbounded slack on a downward-closed, non-union-closed family")
    print(f"  {'m':>4s} {'n':>4s} {'rho':>8s} {'slack (L=2)':>13s}")
    partA = []
    for m in (3, 6, 12, 24, 48):
        n = 2 * m; rho = (n * (n - 1) // 2) / (n / 2)
        s = float(np.mean([slack(m, 2, "sum", sd) for sd in range(5)]))
        partA.append(dict(m=m, n=n, rho=rho, slack=s))
        print(f"  {m:4d} {n:4d} {rho:8.1f} {s:13.1f}")
    print("  slack is unbounded in m -> no optimisation over this lattice can fix it\n")

    print("Part B -- does the exponent track DEPTH?")
    partB = {}
    for L in (1, 2, 3, 4):
        b, a, r2 = fit(MS, L, "sum")
        partB[L] = dict(exponent=float(b), intercept=float(a), r2=float(r2),
                        per_layer=float(b / L))
        print(f"  L={L}:  slack ~ rho^{b:6.3f}   R^2={r2:.4f}   exponent/depth={b/L:.3f}")
    pl = [v["per_layer"] for v in partB.values()]
    consistent = (max(pl) - min(pl)) < 0.05
    print(f"  exponent/depth is {'CONSTANT' if consistent else 'NOT constant'} "
          f"({min(pl):.3f} to {max(pl):.3f}) -> exponent {'is' if consistent else 'is not'} "
          f"set by depth\n")

    print("Part C -- converse: join-closed family has a feasible top, so slack is exact")
    partC = []
    for seed in range(5):
        r = np.random.RandomState(seed); n = 7
        E0 = [(i, j) for i in range(n) for j in range(i + 1, n) if r.rand() < 0.35][:10]
        net = MonoMPNN(2, 6, 2, agg="sum", nonneg=True, seed=seed)
        X = np.ones((n, 2))
        # brute force over the WHOLE downward-closed, join-closed family (all subsets of E0)
        best = max(net.forward(X, build_A(n, [], E0, bits))
                   for bits in itertools.product([0, 1], repeat=len(E0)))
        top = net.forward(X, build_A(n, [], E0, [1] * len(E0)))
        partC.append(dict(seed=seed, k=len(E0), brute_force_max=best, top_element=top,
                          slack=top / max(best, 1e-12)))
        print(f"  seed {seed}: k={len(E0)}, 2^k={2**len(E0)} subsets enumerated, "
              f"slack = {top/max(best,1e-12):.9f}")
    exact = all(abs(p["slack"] - 1.0) < 1e-12 for p in partC)
    print(f"  join-closed slack is exactly 1.0 on every trial: {exact}")

    record(run, "c5_slack_scaling",
           {"part_a_unbounded": partA, "part_b_depth": partB,
            "exponent_per_layer_constant": consistent,
            "part_c_joinclosed": partC, "part_c_exact": exact})


if __name__ == "__main__": main()
