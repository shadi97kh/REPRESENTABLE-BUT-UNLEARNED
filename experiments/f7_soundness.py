"""F7 (step 1): make the certificate SOUND, then price the soundness.

f3b measured two defects. Forcing mandatory pairs present made the lower bound invalid,
and discarding pairs below a probability floor left a median 26% of the sampled ensemble
outside the lattice, so the upper bound did not cover those structures at all.

The fix is structural, not statistical:
    mandatory = []      nothing is forced, so the empty endpoint is a subset of every
                        real structure's edge set
    floor     = 0.0     nothing is discarded, so the lattice is the full power set of
                        base pairs and every secondary structure is inside it
With both, monotonicity gives model(structure) <= model(all edges present) for EVERY
structure in the Boltzmann ensemble, with no coverage caveat and no appeal to sampling.

Raising the floor buys a smaller k and a tighter bound at the cost of coverage. This
experiment sweeps the floor and reports the coverage-versus-tightness curve so the trade
is explicit rather than assumed.

Model: the length-generalizing monotone model trained on real Huesken efficacy
(data/models/monotone_generic.json, held-out Spearman 0.385). Not an untrained network.
"""
import os, random, statistics
import numpy as np
import RNA
from certmp.ensemble import sound_lattice, generic_features
from certmp.reach import exact_interval, build_A
from certmp.train import load_numpy_model
from certmp.provenance import new_run, record
from experiments._huesken import pairs_of

N_SAMPLES, N_SEQ, LENGTH = 1000, 20, 60
RNG_SEED = 20260904
FLOORS = [0.0, 1e-6, 1e-4, 1e-3, 1e-2, 0.05, 0.10, 0.20]
MODEL_PATH = "data/models/monotone_generic.json"


def sample_structures(seq, n):
    md = RNA.md(); md.uniq_ML = 1
    fc = RNA.fold_compound(seq, md)
    _, mfe = fc.mfe(); fc.exp_params_rescale(mfe); fc.pf()
    return list(fc.pbacktrack(n))


def main():
    if not os.path.exists(MODEL_PATH):
        print(f"MODEL ABSENT: {MODEL_PATH}. Run `make f6` first."); raise SystemExit(2)
    run = new_run("f7_soundness", dict(n_samples=N_SAMPLES, n_seq=N_SEQ, length=LENGTH,
                                       floors=FLOORS, rng_seed=RNG_SEED,
                                       model=MODEL_PATH, mandatory="empty"))
    m = load_numpy_model(MODEL_PATH)
    random.seed(1)                       # same sequences as f3b
    RNA.init_rand(RNG_SEED)

    per_floor = {f: {"k": [], "cov": [], "tight": [], "viol": 0, "n": 0} for f in FLOORS}
    for t in range(N_SEQ):
        seq = "".join(random.choice("ACGU") for _ in range(LENGTH))
        X = generic_features(seq)
        n = len(seq)
        structs = sample_structures(seq, N_SAMPLES)
        pairsets = [pairs_of(s) for s in structs]
        ys = np.array([m.forward(X, build_A(n, list(P), [], [])) for P in pairsets])
        ens_max = float(ys.max())
        for f in FLOORS:
            lat = sound_lattice(seq, floor=f)
            _, hi, _ = exact_interval(m, n, lat["mandatory"], lat["optional"], X)
            inside = np.array([P <= set(lat["optional"]) for P in pairsets])
            cov = float(inside.mean())
            viol = int((ys > hi + 1e-9).sum())
            d = per_floor[f]
            d["k"].append(lat["k"]); d["cov"].append(cov)
            d["tight"].append(hi / max(ens_max, 1e-12))
            d["viol"] += viol; d["n"] += len(ys)

    print(f"{'floor':>8s} {'median k':>9s} {'coverage':>10s} {'worst cov':>10s} "
          f"{'lattice_hi/ens_max':>20s} {'bound violations':>18s}")
    rows = []
    for f in FLOORS:
        d = per_floor[f]
        row = dict(floor=f, median_k=statistics.median(d["k"]), max_k=max(d["k"]),
                   median_coverage=statistics.median(d["cov"]), min_coverage=min(d["cov"]),
                   median_tightness=statistics.median(d["tight"]),
                   max_tightness=max(d["tight"]),
                   violations=d["viol"], n_samples=d["n"])
        rows.append(row)
        print(f"{f:8.0e} {row['median_k']:9.0f} {row['median_coverage']*100:9.1f}% "
              f"{row['min_coverage']*100:9.1f}% {row['median_tightness']:20.2f} "
              f"{row['violations']:12d}/{row['n_samples']}")

    sound = rows[0]
    ok = (sound["min_coverage"] >= 1.0 - 1e-12 and sound["violations"] == 0)
    print(f"\nfloor=0, mandatory empty: coverage {sound['median_coverage']*100:.1f}% "
          f"(worst {sound['min_coverage']*100:.1f}%), {sound['violations']} violations in "
          f"{sound['n_samples']} sampled structures -> {'SOUND' if ok else 'NOT SOUND'}")
    print(f"price of soundness: k rises to a median {sound['median_k']:.0f} and the bound "
          f"loosens to {sound['median_tightness']:.2f}x the true sampled maximum")

    # coverage-versus-tightness curve
    png = os.path.join(run, "coverage_tightness.png")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        cov = [r["median_coverage"] * 100 for r in rows]
        tig = [r["median_tightness"] for r in rows]
        ax[0].plot(cov, tig, "o-", color="#2b6cb0")
        for r, c, tg in zip(rows, cov, tig):
            ax[0].annotate(f"{r['floor']:.0e}", (c, tg), textcoords="offset points",
                           xytext=(5, 4), fontsize=8)
        ax[0].set_xlabel("ensemble coverage (% of sampled structures inside lattice)")
        ax[0].set_ylabel("lattice max / sampled ensemble max")
        ax[0].set_title("coverage vs tightness (label = probability floor)")
        ax[0].grid(alpha=.3)
        ax[1].semilogx([max(r["floor"], 1e-7) for r in rows],
                       [r["median_k"] for r in rows], "s-", color="#276749")
        ax[1].set_xlabel("probability floor (1e-7 marks floor=0)")
        ax[1].set_ylabel("median k")
        ax[1].set_title("lattice size vs floor")
        ax[1].grid(alpha=.3)
        fig.tight_layout(); fig.savefig(png, dpi=140)
        print(f"[figure] {png}")
    except Exception as e:
        png = None; print(f"(figure skipped: {e})")

    record(run, "f7_soundness", {"rows": rows, "sound_at_floor_zero": ok, "figure": png,
                                 "model": MODEL_PATH})

if __name__ == "__main__": main()
