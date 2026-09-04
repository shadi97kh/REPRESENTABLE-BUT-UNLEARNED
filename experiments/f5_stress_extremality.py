"""F5: robustness of the two surviving cases (sum, max) beyond the single configuration
used by F1, plus the max-saturation failure mode found while auditing F5's void rate.

NOT PRE-REGISTERED. This is a post-hoc robustness audit, added 2026-09-04 after F1
passed. It can only weaken confidence in the theorem, never establish it; the registered
test is F1. Reported separately for that reason.

Part A  randomise n, k, layers, hidden width, feature scale and the mandatory subgraph.
        A theorem must not depend on any of them.
Part B  a trial whose output is CONSTANT across the lattice is a void, not a confirmation.
        Part A shows max aggregation goes void often. Part B finds out why: a dense
        mandatory subgraph saturates the row-wise max, so optional edges change nothing.
Part C  does that failure mode reach the RNA application? Mandatory pairs come from
        p > 0.90 base pairs, and a base pairs with at most one partner.
"""
import itertools, statistics, numpy as np
from certmp.models import MonoMPNN
from certmp.reach import build_A
from certmp.provenance import new_run, record

TOL = 1e-12

def _sweep_values(model, n, mandatory, optional, X):
    return [model.forward(X, build_A(n, mandatory, optional, bits))
            for bits in itertools.product([0, 1], repeat=len(optional))]

def part_a_trial(agg, seed):
    r = np.random.RandomState(seed)
    n, layers = int(r.randint(4, 11)), int(r.randint(1, 5))
    d_in, d_hid = int(r.randint(1, 6)), int(r.randint(2, 24))
    allp = [(i, j) for i in range(n) for j in range(i + 1, n)]
    k = int(r.randint(1, min(len(allp), 11) + 1))
    idx = r.choice(len(allp), size=k, replace=False)
    optional = [allp[t] for t in idx]
    rest = [p for t, p in enumerate(allp) if t not in set(idx)]
    nf = int(r.randint(0, len(rest) + 1))
    mandatory = [rest[t] for t in r.choice(len(rest), size=nf, replace=False)] if nf else []
    X = r.rand(n, d_in) * float(r.choice([0.01, 1.0, 100.0]))
    m = MonoMPNN(d_in, d_hid, layers, agg=agg, nonneg=True, seed=seed)
    ys = _sweep_values(m, n, mandatory, optional, X)
    yhi = m.forward(X, build_A(n, mandatory, optional, [1] * k))
    ylo = m.forward(X, build_A(n, mandatory, optional, [0] * k))
    ymax, ymin = max(ys), min(ys)
    sc = max(abs(ymax), abs(ymin), 1e-12)
    return dict(gap_hi=(ymax - yhi) / sc, gap_lo=(ylo - ymin) / sc,
                live=(ymax - ymin) > TOL * sc,
                finite=bool(np.isfinite([ymax, ymin]).all()),
                n=n, k=k, layers=layers, d_hid=d_hid, n_mandatory=nf)

def part_b_live_frac(agg, n, k, n_mandatory, seeds=150):
    live = tried = 0
    for s in range(seeds):
        r = np.random.RandomState(20000 + s)
        allp = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if k + n_mandatory > len(allp): return None
        idx = r.choice(len(allp), size=k + n_mandatory, replace=False)
        optional, mandatory = [allp[t] for t in idx[:k]], [allp[t] for t in idx[k:]]
        X = r.rand(n, 3)
        m = MonoMPNN(3, 8, 2, agg=agg, nonneg=True, seed=s)
        ys = _sweep_values(m, n, mandatory, optional, X)
        sc = max(abs(max(ys)), abs(min(ys)), 1e-12)
        tried += 1; live += (max(ys) - min(ys)) > TOL * sc
    return live / tried

def part_c_rna(n_seq=20, length=60):
    import random
    from certmp.ensemble import bpp_lattice
    random.seed(1); degs, dens, ks = [], [], []
    for _ in range(n_seq):
        seq = "".join(random.choice("ACGU") for _ in range(length))
        lat = bpp_lattice(seq); n = lat["n"]; M = lat["mandatory"]
        d = {}
        for (i, j) in M:
            d[i] = d.get(i, 0) + 1; d[j] = d.get(j, 0) + 1
        degs.append(max(d.values()) if d else 0)
        dens.append(len(M) / (n * (n - 1) / 2)); ks.append(lat["k"])
    return dict(max_mandatory_degree=max(degs),
                median_density=statistics.median(dens),
                max_density=max(dens), median_k=statistics.median(ks),
                n_seq=n_seq, length=length)

TRIALS = 400

def main():
    run = new_run("f5_stress_extremality",
                  dict(trials=TRIALS, pre_registered=False, note="post-hoc robustness audit"))
    out = {"part_a": {}, "part_b": {}, "part_c": None}
    print("Part A -- randomised architecture and lattice, sum and max only")
    for agg in ("sum", "max"):
        rs = [part_a_trial(agg, s) for s in range(TRIALS)]
        live = [r for r in rs if r["live"]]
        bad = [r for r in live if r["gap_hi"] > TOL or r["gap_lo"] > TOL]
        out["part_a"][agg] = dict(
            trials=TRIALS, live=len(live), void_dead=TRIALS - len(live),
            non_finite=sum(not r["finite"] for r in rs), violations=len(bad),
            max_gap_hi=max((r["gap_hi"] for r in live), default=0.0),
            max_gap_lo=max((r["gap_lo"] for r in live), default=0.0),
            first_counterexample=bad[0] if bad else None)
        o = out["part_a"][agg]
        print(f"  {agg:4s} live {o['live']:3d}/{TRIALS}  void(dead) {o['void_dead']:3d}  "
              f"violations {o['violations']:3d}  max gap {max(o['max_gap_hi'], o['max_gap_lo']):.3e}")

    print("\nPart B -- void rate against mandatory-subgraph density (n=8, k=5)")
    print(f"  {'|mandatory|':>12s} {'sum live':>10s} {'max live':>10s}")
    for nf in (0, 1, 2, 4, 8, 12, 16, 20, 23):
        fs, fm = part_b_live_frac("sum", 8, 5, nf), part_b_live_frac("max", 8, 5, nf)
        if fs is None: break
        out["part_b"][nf] = dict(sum_live=fs, max_live=fm)
        print(f"  {nf:12d} {fs:10.2f} {fm:10.2f}")

    print("\nPart C -- does the saturation failure mode reach the RNA application?")
    c = part_c_rna(); out["part_c"] = c
    print(f"  mandatory subgraph: max node degree {c['max_mandatory_degree']}, "
          f"median density {c['median_density']*100:.2f}%, max {c['max_density']*100:.2f}%")
    print(f"  -> {'NO, mandatory core is a sparse matching' if c['max_mandatory_degree'] <= 1 else 'YES, investigate'}")

    out["violations_total"] = sum(v["violations"] for v in out["part_a"].values())
    record(run, "f5_stress_extremality", out)

if __name__ == "__main__": main()
