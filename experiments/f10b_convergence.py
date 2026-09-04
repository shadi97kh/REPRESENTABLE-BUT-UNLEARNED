"""F10b: does the maximal-structure set CONVERGE, or does it track the sample size?

f10 showed the maximal-element bound is nearly tight, ratio 1.008 to 1.045 against the
sampled ensemble maximum, where the edge lattice runs 21x to 217x. That is only useful if
the set of maximal structures needed to dominate 95% of mass is a property of the ENSEMBLE
rather than of the sample drawn.

If M_95 keeps growing with the number of samples, then "95% of mass" means 95% of this
particular sample, the bound says nothing about unseen structures, and the route has
degenerated into sampling with extra steps. This sweeps the sample size to find out.
"""
import numpy as np, RNA, statistics
from certmp.ensemble import generic_features
from certmp.maximal import saturate
from certmp.train import load_numpy_model
from certmp.provenance import new_run, record
from experiments._huesken import pairs_of, load_verified, gene_map
from experiments.f9_target_context import (GENCODE, MODEL_PATH, load_transcripts, rc,
                                           ALIASES, RNG_SEED)
from experiments.f10_maximal import coverage_curve, m_for

LENGTHS = (50, 150)
SAMPLE_SIZES = (1000, 5000, 20000)
N_WINDOWS = 4


def main():
    run = new_run("f10b_convergence", dict(lengths=LENGTHS, sample_sizes=SAMPLE_SIZES,
                                           n_windows=N_WINDOWS, model=MODEL_PATH))
    m = load_numpy_model(MODEL_PATH)
    records, rep, ver = load_verified(); gm = gene_map()
    tx = load_transcripts(GENCODE, {ALIASES.get(s.upper(), s.upper()) for s in set(gm.values())})
    sites = []
    for r in records:
        g = gm.get(r["sequence"])
        if not g: continue
        key = ALIASES.get(g.upper(), g.upper())
        if key not in tx: continue
        t = rc(r["sequence"])
        for ti, s in enumerate(tx[key]):
            p = s.find(t)
            if p >= 0: sites.append({"gene": key, "tx": ti, "pos": p}); break
    rng = np.random.RandomState(RNG_SEED)
    chosen = [sites[i] for i in rng.choice(len(sites), size=24, replace=False)][:N_WINDOWS]
    RNA.init_rand(RNG_SEED)

    out = {}
    print(f"{'L':>4s} {'window':>7s} {'N':>7s} {'M':>7s} {'M95':>7s} {'M95/N':>7s} {'ratio':>7s}")
    for L in LENGTHS:
        rows = []
        for wi, st in enumerate(chosen):
            seq = tx[st["gene"]][st["tx"]]
            lo = max(0, st["pos"] + 10 - L // 2); win = seq[lo:lo + L]
            if len(win) < L or set(win) - set("ACGU"): continue
            n = len(win); X = generic_features(win)
            backbone = [(i, i + 1) for i in range(n - 1)]
            md = RNA.md(); md.uniq_ML = 1
            fc = RNA.fold_compound(win, md); _, e = fc.mfe(); fc.exp_params_rescale(e); fc.pf()
            allS = [frozenset(pairs_of(x)) for x in fc.pbacktrack(max(SAMPLE_SIZES))]
            allT = [saturate(P, win, n) for P in allS]
            for N in SAMPLE_SIZES:
                S, T = allS[:N], allT[:N]
                freq = {}
                for t_ in T: freq[t_] = freq.get(t_, 0) + 1
                order = [t_ for t_, _ in sorted(freq.items(), key=lambda kv: -kv[1])]
                c = coverage_curve(S, order)
                m95 = m_for(c, 0.95)
                rows.append(dict(window=wi, gene=st["gene"], N=N, M=len(order), M95=m95,
                                 M95_over_N=(m95 / N) if m95 else None))
                print(f"{L:4d} {wi:7d} {N:7d} {len(order):7d} {m95 if m95 else -1:7d} "
                      f"{(m95/N if m95 else float('nan')):7.3f} "
                      f"{len(order)/N:7.3f}")
        out[L] = rows
        # growth exponent of M95 in N, per window
        exps = []
        for wi in sorted({r["window"] for r in rows}):
            rr = [r for r in rows if r["window"] == wi and r["M95"]]
            if len(rr) < 2: continue
            a = np.polyfit(np.log([r["N"] for r in rr]), np.log([r["M95"] for r in rr]), 1)[0]
            exps.append(float(a))
        out[f"{L}_M95_growth_exponents"] = exps
        print(f"  L={L}: M95 grows as N^{statistics.median(exps):.2f} "
              f"(per-window exponents {[round(e,2) for e in exps]})")
        print(f"  {'CONVERGED' if statistics.median(exps) < 0.15 else 'NOT CONVERGED'}: "
              f"exponent 0 would mean M95 is a property of the ensemble, 1 means it is "
              f"just the sample\n")
    record(run, "f10b_convergence", out)

if __name__ == "__main__": main()
