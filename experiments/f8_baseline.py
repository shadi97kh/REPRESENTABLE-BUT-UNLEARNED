"""F8 (step 2): the cost of certifiability, measured against a matched baseline.

The monotone model gives up signed weights. That is what buys the extremality theorem, so
the accuracy it costs belongs in the paper whichever way it comes out.

Matched by construction: the baseline differs from the monotone model in exactly one
respect, `nonneg=False`. Identical features, identical published split, identical
validation set, identical optimiser, learning rate, batch size, epoch budget, early
stopping and seeds, all through the same certmp.train.fit. Five seeds per configuration,
because a one-seed comparison cannot distinguish a real gap from initialisation noise.

Both feature sets are reported, since both appear elsewhere in the repo:
  positional  one-hot over (position, nucleotide), length-specific, used by f6
  generic     nucleotide plus end-proximity, length-free, used by f7 and f9
"""
import statistics
import numpy as np
from certmp import data as cdata, train as ctrain
from certmp.ensemble import generic_features
from certmp.provenance import new_run, record
from experiments._huesken import load_verified, build_graphs, published_split, gene_map
from experiments.f6_trained import positional_features

SEEDS = [0, 1, 2, 3, 4]
FEATURE_SETS = {"positional": positional_features, "generic": generic_features}


def evaluate(X, A, y, idx_tr, idx_va, idx_te, nonneg, seed):
    model, best_val, best_ep = ctrain.fit(X, A, y, idx_tr, idx_va, seed=seed, nonneg=nonneg)
    pred = ctrain.predict(model, X, A, idx_te)
    return (cdata.spearman(pred, y[idx_te]), cdata.pearson(pred, y[idx_te]), best_val)


def summarise(vals):
    return dict(mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)) if len(vals) > 1 else 0.0,
                min=float(min(vals)), max=float(max(vals)), values=[float(v) for v in vals])


def main():
    run = new_run("f8_baseline", dict(seeds=SEEDS, feature_sets=list(FEATURE_SETS),
                                      matched="identical code path, nonneg toggled"))
    records, rep, ver = load_verified()
    print(f"integrity verified={ver['verified']}, published split {rep['n_train']}/{rep['n_test']}")
    seqs = [r["sequence"] for r in records]
    y = np.array([r["efficacy"] for r in records], dtype=float)

    gm = gene_map()
    gseq = np.array([gm.get(s, "?") for s in seqs])

    out = {}
    for fname, ffn in FEATURE_SETS.items():
        print(f"\nfeatures = {fname}")
        graphs = build_graphs(seqs, ffn)
        X = np.stack([g["X"] for g in graphs]); A = np.stack([g["A_mfe"] for g in graphs])
        idx_tr, idx_va, idx_te = published_split(len(graphs), rep["n_train"])
        res = {}
        for label, nonneg in (("monotone (certifiable)", True), ("unconstrained (baseline)", False)):
            rho, r_, vl = zip(*[evaluate(X, A, y, idx_tr, idx_va, idx_te, nonneg, s) for s in SEEDS])
            res[label] = {"spearman": summarise(rho), "pearson": summarise(r_),
                          "val_mse": summarise(vl), "nonneg": nonneg}
            print(f"  {label:26s} Spearman {np.mean(rho):.3f} +/- {np.std(rho, ddof=1):.3f}   "
                  f"Pearson {np.mean(r_):.3f} +/- {np.std(r_, ddof=1):.3f}")
        m_, b_ = res["monotone (certifiable)"], res["unconstrained (baseline)"]
        res["cost_of_certifiability"] = {
            "spearman_delta": b_["spearman"]["mean"] - m_["spearman"]["mean"],
            "pearson_delta": b_["pearson"]["mean"] - m_["pearson"]["mean"]}
        print(f"  cost of certifiability: Spearman {res['cost_of_certifiability']['spearman_delta']:+.3f}  "
              f"Pearson {res['cost_of_certifiability']['pearson_delta']:+.3f}")
        out[fname] = res

    # gene-disjoint check on the headline feature set
    if gm:
        print("\ngene-disjoint split, positional features")
        graphs = build_graphs(seqs, positional_features)
        X = np.stack([g["X"] for g in graphs]); A = np.stack([g["A_mfe"] for g in graphs])
        genes = sorted(str(g) for g in set(gseq) - {"?"})
        grng = np.random.RandomState(20260904)
        held, acc = [], 0
        for gi in grng.permutation(len(genes)):
            g_ = genes[gi]; c = int((gseq == g_).sum())
            if acc + c > 400 and held: break
            held.append(g_); acc += c
            if acc >= 249: break
        te = np.where(np.isin(gseq, held))[0]
        rest = grng.permutation(np.where(~np.isin(gseq, held))[0])
        gres = {}
        for label, nonneg in (("monotone (certifiable)", True), ("unconstrained (baseline)", False)):
            rho, r_, _ = zip(*[evaluate(X, A, y, rest[200:], rest[:200], te, nonneg, s) for s in SEEDS])
            gres[label] = {"spearman": summarise(rho), "pearson": summarise(r_)}
            print(f"  {label:26s} Spearman {np.mean(rho):.3f} +/- {np.std(rho, ddof=1):.3f}   "
                  f"Pearson {np.mean(r_):.3f} +/- {np.std(r_, ddof=1):.3f}")
        gres["held_out_genes"] = held; gres["n_test"] = int(len(te))
        gres["cost_of_certifiability"] = {
            "spearman_delta": gres["unconstrained (baseline)"]["spearman"]["mean"]
                              - gres["monotone (certifiable)"]["spearman"]["mean"]}
        print(f"  cost of certifiability: Spearman "
              f"{gres['cost_of_certifiability']['spearman_delta']:+.3f}")
        out["gene_disjoint_positional"] = gres

    out["published_baseline_pearson"] = 0.66
    record(run, "f8_baseline", out)

if __name__ == "__main__": main()
