"""F10: test the MAXIMAL-ELEMENT route before committing to a reframing.

The edge lattice bounds everything, including edge subsets that are not valid secondary
structures, and f9 measured the resulting slack at 223x by 150 nt. The alternative is to
bound the union of downward closures of a set of MAXIMAL valid structures. Monotonicity
gives model(S) <= model(T) whenever S is contained in T, so a set {T_i} of maximal
structures bounds every structure inside any of them, and every configuration in that union
is realisable.

The trade is coverage. The lattice covers the ensemble provably and loosely. This covers a
MEASURED fraction of ensemble mass, tightly. This experiment measures both sides: how many
maximal structures are needed for 95% and 99% of sampled Boltzmann mass, how tight the
resulting bound is, and how large the covered set is.

Also fits the slack-scaling prediction from f9: regress log(slack) on log(k / (n/2)),
where k is the number of optional canonical pairs and n/2 bounds the pairs any single
structure can hold, so the ratio measures how much of the lattice is unrealisable.
"""
import glob, json, math, os, statistics
import numpy as np
import RNA
from certmp.ensemble import sound_lattice, generic_features
from certmp.maximal import saturate, is_valid, is_maximal
from certmp.reach import exact_interval, build_A
from certmp.train import load_numpy_model
from certmp.provenance import new_run, record
from experiments._huesken import pairs_of
from experiments.f9_target_context import (GENCODE, MODEL_PATH, load_transcripts, rc,
                                           ALIASES, RNG_SEED)
from experiments._huesken import load_verified, gene_map

LENGTHS = (50, 100, 150)
N_SAMPLES = 5000
N_WINDOWS = 12
TARGETS = (0.95, 0.99)


def coverage_curve(samples, order):
    """Incremental exact domination. Returns, for each prefix of `order`, the fraction of
    sampled structures contained in at least one maximal structure of that prefix.
    Exact, not a lower bound: every still-undominated sample is tested against each new
    maximal structure."""
    undom = list(range(len(samples)))
    covered, curve = 0, []
    for T in order:
        still = []
        for si in undom:
            if samples[si] <= T: covered += 1
            else: still.append(si)
        undom = still
        curve.append(covered / len(samples))
        if not undom: break
    return curve


def m_for(curve, target):
    for i, c in enumerate(curve):
        if c >= target: return i + 1
    return None


def log10_closure(sizes):
    """|union of downward closures| in log10. Lower bound: the largest single closure.
    Upper bound: the sum. Exact union needs inclusion-exclusion over 2^M terms."""
    if not sizes: return None, None
    lo = max(sizes) * math.log10(2)
    hi = math.log10(sum(2.0 ** min(s, 1000) for s in sizes)) if max(sizes) <= 1000 else \
         max(sizes) * math.log10(2) + math.log10(len(sizes))
    return lo, hi


def main():
    for p in (GENCODE, MODEL_PATH):
        if not os.path.exists(p):
            print(f"MISSING: {p}"); raise SystemExit(2)
    run = new_run("f10_maximal", dict(lengths=LENGTHS, n_samples=N_SAMPLES,
                                      n_windows=N_WINDOWS, targets=TARGETS,
                                      model=MODEL_PATH, rng_seed=RNG_SEED))
    m = load_numpy_model(MODEL_PATH)
    records, rep, ver = load_verified()
    gm = gene_map()
    wanted = {ALIASES.get(s.upper(), s.upper()) for s in set(gm.values())}
    tx = load_transcripts(GENCODE, wanted)

    sites = []
    for r in records:
        g = gm.get(r["sequence"])
        if not g: continue
        key = ALIASES.get(g.upper(), g.upper())
        if key not in tx: continue
        t = rc(r["sequence"])
        for ti, s in enumerate(tx[key]):
            pos = s.find(t)
            if pos >= 0:
                sites.append({"gene": key, "tx": ti, "pos": pos}); break
    rng = np.random.RandomState(RNG_SEED)
    chosen = [sites[i] for i in rng.choice(len(sites), size=24, replace=False)][:N_WINDOWS]
    print(f"{len(sites)} target sites located; using {len(chosen)} windows per length, "
          f"{N_SAMPLES} Boltzmann samples each")

    RNA.init_rand(RNG_SEED)
    out, slack_rows = {}, []
    for L in LENGTHS:
        per = []
        for st in chosen:
            seq = tx[st["gene"]][st["tx"]]
            lo = max(0, st["pos"] + 10 - L // 2)
            win = seq[lo:lo + L]
            if len(win) < L or set(win) - set("ACGU"): continue
            n = len(win)
            X = generic_features(win)
            backbone = [(i, i + 1) for i in range(n - 1)]

            md = RNA.md(); md.uniq_ML = 1
            fc = RNA.fold_compound(win, md)
            _, e = fc.mfe(); fc.exp_params_rescale(e); fc.pf()
            samples = [frozenset(pairs_of(x)) for x in fc.pbacktrack(N_SAMPLES)]
            ens_max = max(m.forward(X, build_A(n, backbone + sorted(P), [], [])) for P in samples)

            sats = [saturate(P, win, n) for P in samples]
            freq = {}
            for T in sats: freq[T] = freq.get(T, 0) + 1
            order = [T for T, _ in sorted(freq.items(), key=lambda kv: -kv[1])]
            M = len(order)
            curve = coverage_curve(samples, order)
            ys = {T: m.forward(X, build_A(n, backbone + sorted(T), [], [])) for T in order}

            rowd = dict(gene=st["gene"], M=M, n_samples=len(samples),
                        median_pairs=statistics.median(len(T) for T in order),
                        ens_max=ens_max, full_bound=max(ys.values()),
                        ratio_full=max(ys.values()) / max(ens_max, 1e-12))
            for tgt in TARGETS:
                mm = m_for(curve, tgt)
                if mm is None: continue
                sel = order[:mm]
                b = max(ys[T] for T in sel)
                clo, chi = log10_closure([len(T) for T in sel])
                rowd[f"M_{int(tgt*100)}"] = mm
                rowd[f"bound_{int(tgt*100)}"] = b
                rowd[f"ratio_{int(tgt*100)}"] = b / max(ens_max, 1e-12)
                rowd[f"closure_log10_lo_{int(tgt*100)}"] = clo
                rowd[f"closure_log10_hi_{int(tgt*100)}"] = chi
            # lattice comparison, same window
            lat = sound_lattice(win, floor=0.0, canonical_only=True)
            _, hi_lat, _ = exact_interval(m, n, backbone, lat["optional"], X)
            rowd.update(k=lat["k"], lattice_bound=hi_lat,
                        lattice_slack=hi_lat / max(ens_max, 1e-12),
                        lattice_log10=lat["k"] * math.log10(2))
            per.append(rowd)
            slack_rows.append(dict(n=n, k=lat["k"], slack=rowd["lattice_slack"]))

        med = lambda key: float(statistics.median([r[key] for r in per if key in r]))
        o = dict(n_windows=len(per), median_M=med("M"), max_M=max(r["M"] for r in per),
                 median_ratio_full=med("ratio_full"),
                 median_lattice_slack=med("lattice_slack"),
                 median_lattice_log10=med("lattice_log10"), rows=per)
        for tgt in TARGETS:
            t = int(tgt * 100)
            o[f"median_M_{t}"] = med(f"M_{t}")
            o[f"median_ratio_{t}"] = med(f"ratio_{t}")
            o[f"median_closure_log10_lo_{t}"] = med(f"closure_log10_lo_{t}")
            o[f"median_closure_log10_hi_{t}"] = med(f"closure_log10_hi_{t}")
        out[L] = o
        print(f"\nwindow {L} nt  ({o['n_windows']} windows, {N_SAMPLES} samples each)")
        print(f"  distinct maximal structures M      : median {o['median_M']:.0f}, max {o['max_M']}")
        for tgt in TARGETS:
            t = int(tgt * 100)
            print(f"  M for {t}% of sampled mass          : {o[f'median_M_{t}']:.0f}"
                  f"   bound/ensemble_max {o[f'median_ratio_{t}']:.3f}"
                  f"   covered closure 10^{o[f'median_closure_log10_lo_{t}']:.0f}"
                  f"-10^{o[f'median_closure_log10_hi_{t}']:.0f}")
        print(f"  all M, bound/ensemble_max          : {o['median_ratio_full']:.3f}")
        print(f"  edge lattice on same windows       : 10^{o['median_lattice_log10']:.0f} "
              f"subsets, slack {o['median_lattice_slack']:.1f}x")

    # ---- slack scaling: log(slack) ~ a * log(k / (n/2)) ----------------------
    x = np.log(np.array([r["k"] / (r["n"] / 2) for r in slack_rows]))
    yv = np.log(np.array([r["slack"] for r in slack_rows]))
    a, b = np.polyfit(x, yv, 1)
    pred = a * x + b
    ss_res = float(((yv - pred) ** 2).sum()); ss_tot = float(((yv - yv.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    print(f"\nslack scaling over {len(slack_rows)} windows across all three lengths")
    print(f"  log(slack) = {a:.3f} * log(k/(n/2)) + {b:.3f}   R^2 = {r2:.3f}")
    print(f"  exponent {a:.2f}: slack grows about as the {a:.2f} power of lattice redundancy")

    record(run, "f10_maximal", {"by_length": out,
                                "slack_scaling": {"exponent": float(a), "intercept": float(b),
                                                  "r2": float(r2), "n": len(slack_rows),
                                                  "rows": slack_rows}})

if __name__ == "__main__": main()
