"""F6: replace F3's untrained-model numbers with a trained model.

Every magnitude in F3 and F3b came from a randomly initialised network, so they describe
the initialisation and say nothing about siRNA efficacy. This experiment trains the
monotone MPNN on real measured efficacy and repeats the comparison.

Pipeline
  1. Load a redistributed Huesken et al. 2005 table and REJECT it unless it matches the
     statistics published in the paper (certmp.data.verify). A dataset that does not
     verify is not silently used.
  2. Build a graph per siRNA: nodes are the 19 nucleotides, features are one-hot (H4),
     backbone edges (i, i+1) are mandatory, and ViennaRNA base-pair probabilities supply
     the mandatory/optional split exactly as in the RNA application.
  3. Train the torch monotone MPNN (sum aggregation, non-negative weights, non-negative
     features) on the MFE graph, which is what current pipelines condition on.
  4. Report held-out Spearman and Pearson against the published Pearson r = 0.66.
  5. Export to the numpy model and repeat the F3/F3b comparison with trained weights.

Honesty constraints
  - The PUBLISHED 2182/249 split is used, as redistributed with DSIR, not a random
     resplit. Gene annotations, when present, are used to report whether that split leaks
     targets across the train/test boundary rather than assuming it does not.
  - The published r = 0.66 came from a large feedforward ensemble on hand-built features.
     A sign-constrained monotone GNN is a strictly smaller hypothesis class. Any shortfall
     is the price of certifiability and is reported as such, not explained away.
"""
import os, sys, json, random, statistics
import numpy as np
import RNA
from certmp import data as cdata
from certmp.models import MonoMPNN, TorchMonoMPNN
from certmp.ensemble import bpp_lattice, onehot_features
from certmp.reach import exact_interval, build_A
from certmp.provenance import new_run, record
from experiments.f3b_relaxation_slack import pairs_of, sample_structures

TRAIN_PATH  = os.environ.get("HUESKEN_TRAIN", "data/TrainAll2182.txt")
TEST_PATH   = os.environ.get("HUESKEN_TEST",  "data/TestAll249.txt")
ANNOT_PATH  = os.environ.get("HUESKEN_ANNOT", "data/Huesken_2431_annotated.tsv")
VAL_SEED    = 20260904     # only for carving a validation set out of TRAIN
EPOCHS      = 250
LR          = 0.01
BATCH       = 128
D_HID       = 32
N_LAYERS    = 2
BAND_LO, BAND_HI = 0.05, 0.90
N_SAMPLES   = 200          # Boltzmann samples per sequence in the F3b repeat


def positional_features(seq):
    """Non-negative features (H4) that retain POSITION. Sum aggregation is permutation
    invariant, so plain one-hot nucleotides make the model a bag of nucleotides and it
    cannot express position-specific preferences, which is most of the siRNA signal.
    One-hot over (position, nucleotide) fixes that and stays non-negative."""
    idx = {c: t for t, c in enumerate("ACGU")}
    L = len(seq)
    X = np.zeros((L, L * 4))
    for t, c in enumerate(seq.replace("T", "U")):
        if c in idx: X[t, t * 4 + idx[c]] = 1.0
    return X


def build_graphs(seqs):
    """Per sequence: one-hot features, MFE adjacency for training, and the lattice."""
    out = []
    for s in seqs:
        lat = bpp_lattice(s, lo=BAND_LO, hi=BAND_HI)
        n = lat["n"]
        backbone = [(i, i + 1) for i in range(n - 1)]
        mfe_pairs = sorted(pairs_of(lat["mfe_structure"]))
        A_mfe = build_A(n, backbone + mfe_pairs, [], [])
        mandatory = sorted(set(backbone) | set(lat["mandatory"]))
        optional = [e for e in lat["optional"] if e not in set(mandatory)]
        out.append(dict(seq=s, n=n, X=positional_features(s), A_mfe=A_mfe,
                        mandatory=mandatory, optional=optional, k=len(optional),
                        mfe_structure=lat["mfe_structure"]))
    return out


def train(graphs, y, idx_tr, idx_va, seed=0, verbose=True):
    import torch
    torch.manual_seed(seed)
    X = torch.tensor(np.stack([g["X"] for g in graphs]), dtype=torch.float32)
    A = torch.tensor(np.stack([g["A_mfe"] for g in graphs]), dtype=torch.float32)
    Y = torch.tensor(np.asarray(y), dtype=torch.float32)
    tr = torch.tensor(idx_tr); va = torch.tensor(idx_va)
    model = TorchMonoMPNN(X.shape[2], D_HID, N_LAYERS, agg="sum", seed=seed)

    # Put the readout on the target's scale FIRST. Without this the raw output is ~1e3
    # against targets ~1, and the run spends every epoch fixing the scale instead of
    # learning; that is what produced the negative correlation on the first attempt.
    with torch.no_grad():
        raw = model(X[tr], A[tr])
        s = float(Y[tr].std() / raw.std().clamp(min=1e-12))
        model.set_output_affine(s, float(Y[tr].mean() - s * raw.mean()))
        if verbose:
            print(f"    readout initialised: scale {s:.3e}, "
                  f"init train MSE {float(torch.nn.functional.mse_loss(model(X[tr], A[tr]), Y[tr])):.4f}")

    opt = torch.optim.Adam(model.parameters(), lr=LR)
    sched = torch.optim.lr_scheduler.ReduceLROnPlateau(opt, factor=0.5, patience=20)
    best, best_state, best_ep = float("inf"), None, -1
    g = torch.Generator().manual_seed(seed)
    for ep in range(EPOCHS):
        model.train()
        perm = tr[torch.randperm(len(tr), generator=g)]
        for i in range(0, len(perm), BATCH):
            b = perm[i:i + BATCH]
            opt.zero_grad()
            torch.nn.functional.mse_loss(model(X[b], A[b]), Y[b]).backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            vl = float(torch.nn.functional.mse_loss(model(X[va], A[va]), Y[va]))
            tl = float(torch.nn.functional.mse_loss(model(X[tr], A[tr]), Y[tr]))
        sched.step(vl)
        if vl < best - 1e-7:
            best, best_ep = vl, ep
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        if verbose and ep % 25 == 0:
            print(f"    epoch {ep:4d}  train {tl:.4f}  val {vl:.4f}")
    model.load_state_dict(best_state)
    if verbose: print(f"    best val MSE {best:.4f} at epoch {best_ep}")
    return model, 0.0, 1.0, best_ep


def main():
    run = new_run("f6_trained", dict(train_path=TRAIN_PATH, test_path=TEST_PATH,
                                     val_seed=VAL_SEED, epochs=EPOCHS, lr=LR, d_hid=D_HID,
                                     n_layers=N_LAYERS, n_samples=N_SAMPLES, batch=BATCH,
                                     features="one-hot (position, nucleotide)"))
    # ---- 1. acquire + verify -------------------------------------------------
    for p_ in (TRAIN_PATH, TEST_PATH):
        if not os.path.exists(p_):
            print(f"DATASET ABSENT: {p_}")
            print("F6 cannot run without the real Huesken table. Nothing is synthesised.")
            record(run, "f6_trained", {"status": "dataset_absent", "path": p_})
            sys.exit(2)
    tr_rec, te_rec, split_report = cdata.load_dsir_split(TRAIN_PATH, TEST_PATH)
    records = tr_rec + te_rec
    print(f"loaded {split_report['n_train']} train + {split_report['n_test']} test "
          f"= {len(records)}  (published split, overlap {split_report['train_test_overlap']})")
    ver = cdata.verify("huesken", records)
    print(f"integrity: verified={ver['verified']}  problems={ver['problems']}")
    if not ver["verified"]:
        print("REJECTED. The file does not match the published statistics; F6 stops here.")
        record(run, "f6_trained", {"status": "integrity_rejected", "verification": ver})
        sys.exit(3)
    exp = cdata.EXPECTED["huesken"]
    assert split_report["n_train"] == exp["n_train"] and split_report["n_test"] == exp["n_test"]
    print(f"split sizes match the published 2182/249 exactly")

    # gene-level leakage across the published split, if annotations are available
    leak = {"checked": False}
    if os.path.exists(ANNOT_PATH):
        import csv
        gene = {}
        with open(ANNOT_PATH) as fh:
            rd = csv.reader(fh, delimiter="\t"); hdr = next(rd)
            gi = next((i for i, h in enumerate(hdr) if h.lower() in ("gene", "gene_symbol")), None)
            si = next((i for i, h in enumerate(hdr) if "21mer" in h.lower()), None)
            if gi is not None and si is not None:
                for row in rd:
                    if len(row) > max(gi, si):
                        gene[row[si].strip().upper().replace("T", "U")] = row[gi].strip()
        gtr = {gene[r["sequence"]] for r in tr_rec if r["sequence"] in gene}
        gte = {gene[r["sequence"]] for r in te_rec if r["sequence"] in gene}
        leak = {"checked": True, "annotated_rows": len(gene), "genes_train": len(gtr),
                "genes_test": len(gte), "genes_shared": len(gtr & gte),
                "shared": sorted(gtr & gte)[:10]}
        print(f"gene leakage across published split: {leak['genes_shared']} genes appear in "
              f"BOTH train and test (train {leak['genes_train']}, test {leak['genes_test']}, "
              f"annotated {leak['annotated_rows']}/{len(records)})")

    seqs = [r["sequence"] for r in records]
    y = np.array([r["efficacy"] for r in records], dtype=float)
    n_tr_pub = split_report["n_train"]

    # ---- 2. graphs -----------------------------------------------------------
    print("folding sequences and building graphs ...")
    graphs = build_graphs(seqs)
    print(f"  median optional pairs per siRNA k = {statistics.median([g['k'] for g in graphs]):.0f}")

    # ---- 3. split + train ----------------------------------------------------
    idx_te = np.arange(n_tr_pub, len(graphs))              # the published 249
    pub_tr = np.arange(n_tr_pub)
    rng = np.random.RandomState(VAL_SEED)                  # validation carved from TRAIN only
    shuf = rng.permutation(pub_tr)
    idx_va, idx_tr = shuf[:200], shuf[200:]
    print(f"published test {len(idx_te)}; train {len(idx_tr)} + val {len(idx_va)} "
          f"(val carved from published train, seed {VAL_SEED})")
    model, mu, sd, best_ep = train(graphs, y, idx_tr, idx_va, seed=0)

    # ---- 4. held-out correlation --------------------------------------------
    import torch  # noqa: F401  (also used by the gene-disjoint block below)
    Xt = torch.tensor(np.stack([graphs[i]["X"] for i in idx_te]), dtype=torch.float32)
    At = torch.tensor(np.stack([graphs[i]["A_mfe"] for i in idx_te]), dtype=torch.float32)
    with torch.no_grad():
        pred = model(Xt, At).numpy() * sd + mu
    truth = y[idx_te]
    rho = cdata.spearman(pred, truth); r = cdata.pearson(pred, truth)
    print(f"\nheld-out (n={len(idx_te)})  Spearman rho = {rho:.3f}   Pearson r = {r:.3f}")
    print(f"published Huesken 2005 baseline: Pearson r = 0.66 "
          f"(feedforward ensemble on hand-built features, target-wise split)")

    # ---- 4b. gene-disjoint split --------------------------------------------
    # The published split shares all 30 annotated genes between train and test, so the
    # held-out number above measures generalisation to new siRNAs against KNOWN targets.
    # A gene-disjoint split measures generalisation to a NEW target, which is the harder
    # and more honest question for a design tool. Report both, prefer neither silently.
    gene_eval = {"checked": False}
    if leak.get("checked") and leak["genes_shared"] > 0:
        import csv as _csv
        gmap = {}
        with open(ANNOT_PATH) as fh:
            rd = _csv.reader(fh, delimiter="\t"); hdr = next(rd)
            gi = next(i for i, h in enumerate(hdr) if h.lower() == "gene")
            si = next(i for i, h in enumerate(hdr) if "21mer" in h.lower())
            for row in rd:
                if len(row) > max(gi, si):
                    gmap[row[si].strip().upper().replace("T", "U")] = row[gi].strip()
        gseq = np.array([gmap.get(g_["seq"], "?") for g_ in graphs])
        genes = sorted(str(g_) for g_ in set(gseq) - {"?"})
        grng = np.random.RandomState(VAL_SEED)
        order = grng.permutation(len(genes))
        held, acc = [], 0
        for gi_ in order:                       # hold out whole genes until ~249 rows
            g_ = genes[gi_]; c = int((gseq == g_).sum())
            if acc + c > 400 and held: break
            held.append(str(g_)); acc += c
            if acc >= 249: break
        te_mask = np.isin(gseq, held)
        gidx_te = np.where(te_mask)[0]
        gidx_rest = np.where(~te_mask)[0]
        gshuf = grng.permutation(gidx_rest)
        gidx_va, gidx_tr = gshuf[:200], gshuf[200:]
        print(f"\ngene-disjoint split: held-out genes {held} -> "
              f"test {len(gidx_te)} rows, train {len(gidx_tr)}, val {len(gidx_va)}")
        gmodel, _, _, gbest = train(graphs, y, gidx_tr, gidx_va, seed=0, verbose=False)
        Xg = torch.tensor(np.stack([graphs[i]["X"] for i in gidx_te]), dtype=torch.float32)
        Ag = torch.tensor(np.stack([graphs[i]["A_mfe"] for i in gidx_te]), dtype=torch.float32)
        with torch.no_grad():
            gp = gmodel(Xg, Ag).numpy()
        gt = y[gidx_te]
        grho, gr = cdata.spearman(gp, gt), cdata.pearson(gp, gt)
        print(f"gene-disjoint held-out (n={len(gidx_te)})  Spearman rho = {grho:.3f}   "
              f"Pearson r = {gr:.3f}")
        gene_eval = {"checked": True, "held_out_genes": held, "n_test": len(gidx_te),
                     "n_train": len(gidx_tr), "spearman": grho, "pearson": gr,
                     "best_epoch": gbest}

    # ---- 5. certificates with the TRAINED model ------------------------------
    nm = model.export_numpy()
    with torch.no_grad():
        chk_t = float(model(Xt[:1], At[:1])[0])
    chk_n = nm.forward(graphs[idx_te[0]]["X"], graphs[idx_te[0]]["A_mfe"])
    parity = abs(chk_t - chk_n) / max(abs(chk_n), 1e-12)
    print(f"torch/numpy parity on exported model: rel diff {parity:.2e}")
    assert parity < 1e-4, "exported numpy model disagrees with the trained torch model"

    RNA.init_rand(20260904)
    rows = []
    n_zero_k = 0
    for i in idx_te:
        g = graphs[i]
        if g["k"] == 0:
            n_zero_k += 1; continue
        if len(rows) >= 40: break
        lo, hi, _ = exact_interval(nm, g["n"], g["mandatory"], g["optional"], g["X"])
        y_mfe = nm.forward(g["X"], g["A_mfe"])
        backbone = [(t, t + 1) for t in range(g["n"] - 1)]
        ys = []
        for s in sample_structures(g["seq"], N_SAMPLES):
            ys.append(nm.forward(g["X"], build_A(g["n"], backbone + sorted(pairs_of(s)), [], [])))
        rows.append(dict(k=g["k"], y_mfe=y_mfe, lattice_hi=hi, lattice_lo=lo,
                         ens_max=max(ys), ens_min=min(ys),
                         slack_ratio=hi / max(max(ys), 1e-12),
                         underestimate_lattice=(hi - y_mfe) / max(abs(y_mfe), 1e-12),
                         underestimate_ensemble=(max(ys) - y_mfe) / max(abs(y_mfe), 1e-12)))
    med = lambda k: float(statistics.median([r_[k] for r_ in rows]))
    print(f"\nTRAINED model, {len(rows)} held-out siRNAs, {N_SAMPLES} Boltzmann samples each "
          f"({n_zero_k} of the scanned held-out rows had k=0 and were skipped):")
    print(f"  MFE vs lattice worst case   : {med('underestimate_lattice')*100:6.1f}%")
    print(f"  MFE vs sampled ensemble max : {med('underestimate_ensemble')*100:6.1f}%")
    print(f"  relaxation slack ratio      : x{med('slack_ratio'):.2f}")

    record(run, "f6_trained", {
        "status": "ok", "verification": ver, "amendments": cdata.AMENDMENTS,
        "split": {"train": len(idx_tr), "val": len(idx_va), "test": len(idx_te),
                  "val_seed": VAL_SEED, "published_split_reproduced": True,
                  "split_report": split_report, "gene_leakage": leak},
        "best_epoch": best_ep, "held_out": {"n": len(idx_te), "spearman": rho, "pearson": r,
                                            "published_pearson_baseline": 0.66},
        "parity_rel_diff": parity, "gene_disjoint_eval": gene_eval,
        "trained_comparison": {"rows": rows,
                               "median_underestimate_lattice": med("underestimate_lattice"),
                               "median_underestimate_ensemble": med("underestimate_ensemble"),
                               "median_slack_ratio": med("slack_ratio"),
                               "n_zero_k_skipped": n_zero_k,
                               "median_k": float(statistics.median([r_["k"] for r_ in rows]))}})

if __name__ == "__main__": main()
