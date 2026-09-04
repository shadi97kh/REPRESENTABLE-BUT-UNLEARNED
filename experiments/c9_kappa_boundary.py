"""C9: where does the kappa_L characterisation stop holding?

C7 establishes gap == kappa_L for sum aggregation, uniform features, ReLU (affine regime).
Each hypothesis is load-bearing, so each is removed in turn. Reporting the boundary is
worth more than pretending the theorem is universal.

  A  NON-UNIFORM features: f becomes a feature-WEIGHTED walk count, so equality should
     break. Question is whether kappa_L still upper-bounds the gap.
  B  MAX aggregation: nonlinear, so the walk-count derivation does not apply at all.
  C  TANH: monotone (H3 holds, endpoint exactness survives) but nonlinear, so affinity
     and hence the walk-count argument fail.

A bound that survives A/B/C is far more useful than an equality that only holds in the
affine corner. Whichever way it lands, it is stated.
"""
import copy
import numpy as np
from certmp.models import MonoMPNN
from certmp.kappa import adjacency, kappa_L, measured_gap, enumerate_family, FAMILIES
from certmp.provenance import new_run, record

N, N_CAND, SEEDS = 7, 12, 8


def zero_bias(m):
    m = copy.deepcopy(m); m.B = [np.zeros_like(b) for b in m.B]; return m


def run_case(label, agg, act, X, cand, feas, A_top, A_feas, depths=(1, 2, 3)):
    out = []
    for L in depths:
        k, _, _ = kappa_L(N, cand, feas, L)
        gaps = []
        for s in range(SEEDS):
            m = MonoMPNN(X.shape[1], 6, L, agg=agg, nonneg=True, act=act, seed=s)
            gaps.append(measured_gap(zero_bias(m), X, A_top, A_feas)[0])
        g = float(np.mean(gaps)); gmax = float(np.max(gaps))
        eq = abs(g - k) / k < 1e-12
        bounded = gmax <= k + 1e-9
        out.append(dict(case=label, agg=agg, act=act, L=L, kappa=k, gap_mean=g,
                        gap_max=gmax, equals_kappa=bool(eq), bounded_by_kappa=bool(bounded)))
        print(f"{label:>24s} {agg:>4s} {act:>5s} L={L}  kappa={k:9.4f}  "
              f"gap(mean)={g:9.4f}  gap(max)={gmax:9.4f}  "
              f"{'EQUAL' if eq else 'differs':>8s}  {'bounded' if bounded else '*** EXCEEDS ***'}")
    return out


def main():
    run = new_run("c9_kappa_boundary", dict(n=N, n_candidate=N_CAND, seeds=SEEDS))
    cand = [(i, j) for i in range(N) for j in range(i + 1, N)][:N_CAND]
    A_top = adjacency(N, cand)
    feas = enumerate_family(N, cand, dict(FAMILIES)["matching (deg<=1)"])
    A_feas = [adjacency(N, g) for g in feas]
    rng = np.random.RandomState(0)
    X_unif = np.ones((N, 2))
    X_var = rng.rand(N, 2) * 2.0

    rows = []
    print("control: the regime C7 established\n")
    rows += run_case("uniform X, sum, relu", "sum", "relu", X_unif, cand, feas, A_top, A_feas)
    print("\nA -- non-uniform features (weighted walk count)\n")
    rows += run_case("VARIED X, sum, relu", "sum", "relu", X_var, cand, feas, A_top, A_feas)
    print("\nB -- max aggregation (nonlinear, derivation does not apply)")
    print("   uniform X saturates max, so a gap of exactly 1.0 there is DEGENERATE,")
    print("   not evidence of tightness. Varied features are the informative case.\n")
    rows += run_case("uniform X, max, relu", "max", "relu", X_unif, cand, feas, A_top, A_feas)
    rows += run_case("VARIED X, max, relu", "max", "relu", X_var, cand, feas, A_top, A_feas)
    print("\nC -- tanh (monotone but nonlinear, so affinity fails)")
    print("   tanh saturates toward 1, which also compresses the ratio.\n")
    rows += run_case("uniform X, sum, tanh", "sum", "tanh", X_unif, cand, feas, A_top, A_feas)
    rows += run_case("VARIED X, sum, tanh", "sum", "tanh", X_var, cand, feas, A_top, A_feas)

    ctrl = [r for r in rows if r["case"].startswith("uniform X, sum, relu")]
    print(f"\n  control equality holds:            {all(r['equals_kappa'] for r in ctrl)}")
    for lab in ("VARIED X, sum, relu", "uniform X, max, relu", "VARIED X, max, relu",
                "uniform X, sum, tanh", "VARIED X, sum, tanh"):
        sub = [r for r in rows if r["case"] == lab]
        print(f"  {lab:>22s}: equality {all(r['equals_kappa'] for r in sub)}, "
              f"bounded by kappa {all(r['bounded_by_kappa'] for r in sub)}")
    allb = all(r["bounded_by_kappa"] for r in rows)
    print(f"\n  kappa_L upper-bounds the gap in EVERY case tested: {allb}")
    print("  equality holds only in the affine + uniform-feature regime (the control).")
    record(run, "c9_kappa_boundary", {"rows": rows, "kappa_bounds_all_cases": allb})


if __name__ == "__main__":
    main()
