"""P3 RNA pilot: on-target efficacy and ranking, on resolved inputs only.

SCOPE. On-target efficacy and candidate ranking. Not off-target safety, not
causal structure validation, and no thermodynamic claim for modified chemistry
(there is none in this dataset).

WHAT MAKES A ROW ADMISSIBLE. Not the dataset name. A row enters only if its
target window lies wholly inside the tiled-hull lower bound on the inserted
fragment, so its native flanking sequence is inside the insert under the stated
contiguity assumption. Edge windows are excluded as `unresolved_context` rather
than filled in with a plausible accession.

LABEL SEMANTICS. The assay outcome labels the ENSEMBLE. It is fitted through
E_p[f] = b + a^T mu + gamma * E_p[g]. No individual latent structure is ever
given the assay outcome.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import torch

from . import MEASURED_RNA
from .provenance import sha256_file

RC = str.maketrans("ACGU", "UGCA")


def revcomp(s):
    return s.translate(RC)[::-1]


def _spearman(a, b):
    from scipy.stats import rankdata
    ra, rb = rankdata(a), rankdata(b)
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def build_rows(cfg, limit=None):
    """Interior windows only, with isoform ambiguity recorded per site."""
    sys.path.insert(0, ".")
    from experiments._huesken import gene_map, load_verified
    from experiments.f9_target_context import GENCODE, load_transcripts

    W = cfg["windows"]["size"]
    recs, _, _ = load_verified()
    gm = gene_map()
    genes = sorted({gm.get(r["sequence"], "?") for r in recs} - {"?"})
    tx = load_transcripts(GENCODE, genes)

    hits = []
    for r in recs:
        g = gm.get(r["sequence"], "?")
        if g not in tx:
            continue
        site = revcomp(r["sequence"].upper().replace("T", "U"))
        found = []
        for k, iso in enumerate(tx[g]):
            p = iso.find(site)
            if p < 0:
                p = iso.find(site[2:])
            if p >= 0:
                found.append((k, p, iso))
        if found:
            k, p, iso = found[0]
            hits.append({"gene": g, "guide": r["sequence"], "site": site,
                         "efficacy": r["efficacy"], "pos": p, "iso_index": k,
                         "n_isoform_hits": len({f[0] for f in found}),
                         "iso": iso})

    hull = {}
    for h in hits:
        d = hull.setdefault(h["gene"], [1 << 30, -1])
        d[0] = min(d[0], h["pos"])
        d[1] = max(d[1], h["pos"] + len(h["site"]))

    rows, excluded = [], 0
    for h in hits:
        lo, hi = hull[h["gene"]]
        c = h["pos"] + len(h["site"]) // 2
        s, e = c - W // 2, c + W // 2
        if s < lo or e > hi:
            excluded += 1
            continue
        win = h["iso"][s:e]
        if len(win) != W:
            excluded += 1
            continue
        rows.append({"gene": h["gene"], "guide": h["guide"], "window": win,
                     "efficacy": h["efficacy"],
                     "context_status": cfg["windows"]["context_status"],
                     "n_isoform_hits": h["n_isoform_hits"],
                     "cluster": min(h["guide"][i:i + 13]
                                    for i in range(max(len(h["guide"]) - 12, 1)))})
    if limit:
        rows = rows[:limit]
    return rows, {"excluded_edge_or_short": excluded, "window_nt": W,
                  "n_admissible": len(rows),
                  "isoform_ambiguous": sum(1 for r in rows
                                           if r["n_isoform_hits"] > 1)}


def novel_sequence_split(rows, frac, seed):
    """Hold out whole (gene, sequence-cluster) groups, never individual rows."""
    rng = np.random.RandomState(seed)
    keys = sorted({(r["gene"], r["cluster"]) for r in rows})
    rng.shuffle(keys)
    n_hold = max(1, int(len(keys) * frac))
    held = set(keys[:n_hold])
    te = [i for i, r in enumerate(rows) if (r["gene"], r["cluster"]) in held]
    tr = [i for i, r in enumerate(rows) if (r["gene"], r["cluster"]) not in held]
    leak = len({(rows[i]["gene"], rows[i]["cluster"]) for i in te}
               & {(rows[i]["gene"], rows[i]["cluster"]) for i in tr})
    return tr, te, {"n_groups": len(keys), "n_held_groups": n_hold,
                    "n_train": len(tr), "n_test": len(te),
                    "group_leak": leak}


def onehot(seq, alphabet="ACGU"):
    idx = {c: i for i, c in enumerate(alphabet)}
    X = np.zeros((len(seq), len(alphabet)))
    for t, c in enumerate(seq):
        if c in idx:
            X[t, idx[c]] = 1.0
    return X


def top_k_enrichment(pred, y, k, threshold):
    """Precision@k against the predeclared activity threshold, and its prevalence."""
    order = np.argsort(-pred)[:k]
    active = (y >= threshold)
    prev = float(active.mean())
    prec = float(active[order].mean()) if len(order) else float("nan")
    return {"precision_at_k": prec, "prevalence": prev,
            "enrichment": (prec / prev) if prev > 0 else float("nan"), "k": k}


# ------------------------------------------------------------------- the arms
def run_arm(arm, rows, tr, te, cfg, seed, device):
    """One arm. Sequence arms are ridge; window arms go through E_p[f]."""
    y = np.array([r["efficacy"] for r in rows])
    t0 = time.perf_counter()
    extra = {}

    if arm in ("matched_sequence_chemistry", "target_sequence"):
        # guide one-hot, or its reverse complement: the same information by
        # construction, which is why a null difference here is expected and is
        # NOT evidence about target context.
        seqs = [r["guide"] if arm == "matched_sequence_chemistry"
                else revcomp(r["guide"].upper().replace("T", "U")) for r in rows]
        L = max(len(s) for s in seqs)
        Z = np.stack([onehot(s.ljust(L, "N"))[:L].reshape(-1) for s in seqs])
        D = np.hstack([np.ones((len(Z), 1)), Z])
        lam = 1.0
        A = D[tr].T @ D[tr] + lam * np.eye(D.shape[1])
        coef = np.linalg.solve(A, D[tr].T @ y[tr])
        pred = D @ coef
        extra["n_features"] = D.shape[1]

    elif arm == "mfe_gnn":
        from ..model import CertifiedModel, backbone_adjacency
        from ..oracles.rna_noncrossing import RNANonCrossingOracle
        pred = _window_model(rows, tr, te, y, seed, device, gamma=0.0,
                             mode="mfe", cfg=cfg, extra=extra)

    elif arm == "sampled_ensemble_gnn":
        pred = _window_model(rows, tr, te, y, seed, device, gamma=1.0,
                             mode="ensemble", cfg=cfg, extra=extra)

    elif arm == "anchor_A":
        pred = _window_model(rows, tr, te, y, seed, device, gamma=0.0,
                             mode="anchor", cfg=cfg, extra=extra)

    elif arm == "repaired_A_plus_B":
        pred = _window_model(rows, tr, te, y, seed, device, gamma=1.0,
                             mode="a_init", cfg=cfg, extra=extra)
    else:
        raise ValueError(arm)

    m = cfg["metrics"]
    out = {"arm": arm, "seed": seed,
           "spearman": _spearman(pred[te], y[te]),
           "test_mse": float(((pred[te] - y[te]) ** 2).mean()),
           "wall_s": round(time.perf_counter() - t0, 3), "info": extra}
    out.update({"selection": top_k_enrichment(pred[te], y[te], m["k"],
                                              m["activity_threshold"])})
    # target-wise, never pooled into one population claim
    per_gene = {}
    genes = np.array([rows[i]["gene"] for i in te])
    for g in sorted(set(genes)):
        sel = np.where(genes == g)[0]
        if len(sel) >= 8:
            per_gene[g] = {"n": int(len(sel)),
                           "spearman": _spearman(pred[np.array(te)[sel]],
                                                 y[np.array(te)[sel]])}
    out["per_gene"] = per_gene
    return out, pred


def _window_model(rows, tr, te, y, seed, device, gamma, mode, cfg, extra):
    """Fit the ensemble objective yhat = b + a^T mu + gamma * E_p[g].

    The label attaches to the ENSEMBLE. mu is the exact marginal vector of the
    declared prior, so the additive part needs only first moments; only the
    residual needs an expectation, and its sampling error is recorded.
    """
    from ..model import CertifiedModel, backbone_adjacency
    from ..oracles.rna_noncrossing import RNANonCrossingOracle
    from ..rna.ensemble import DeclaredGibbsPrior

    W = cfg["windows"]["size"]
    n_samp = cfg["ensemble"]["n_samples"][-1]
    torch.manual_seed(seed)
    rng = np.random.RandomState(seed)

    # Window features depend only on the window, not on the arm or the seed.
    # Recomputing them per (arm, seed) meant 30 O(n^3) passes over 1153 windows.
    need_samples = (mode == "ensemble")
    feats, mus, gvals, serrs = [], [], [], []
    for r in rows:
        c = r.get("_feat")
        if c is not None and (not need_samples or "g" in c):
            feats.append(c["f"]); mus.append(c["mu"])
            gvals.append(c.get("g", 0.0)); serrs.append(c.get("se", 0.0))
            continue
        o = RNANonCrossingOracle(r["window"], 3, True)
        edges = o.ground_set()
        theta = {e: 0.0 for e in edges}
        prior = DeclaredGibbsPrior(o, theta)
        mu = prior.marginals().mu if edges else np.zeros(0)
        # summary features of the window's ensemble: pair density and mean span
        dens = float(mu.sum()) / max(W, 1)
        span = float(np.mean([(b - a) for (a, b) in edges])) if edges else 0.0
        wmu = float(mu.mean()) if len(mu) else 0.0
        feats.append([dens, span / W, wmu, len(edges) / max(W, 1)])
        mus.append(mu)
        if mode == "ensemble":
            s = prior.sample(min(n_samp, 64), rng)
            vals = np.array([len(x) for x in s], dtype=float)
            gvals.append(float(vals.mean()))
            serrs.append(float(vals.std(ddof=1) / np.sqrt(max(len(vals), 1))))
        else:
            gvals.append(0.0)
            serrs.append(0.0)
        cache = r.get("_feat") or {}
        cache["f"] = feats[-1]; cache["mu"] = mus[-1]
        if need_samples:
            cache["g"] = gvals[-1]; cache["se"] = serrs[-1]
        r["_feat"] = cache

    F = np.array(feats)
    guide = np.stack([onehot(r["guide"].upper().replace("T", "U").ljust(21, "N"))[:21]
                      .reshape(-1) for r in rows])
    if mode == "anchor" or gamma == 0.0:
        D = np.hstack([np.ones((len(F), 1)), guide, F])
    else:
        D = np.hstack([np.ones((len(F), 1)), guide, F,
                       np.array(gvals).reshape(-1, 1)])
    lam = 1.0
    A = D[tr].T @ D[tr] + lam * np.eye(D.shape[1])
    coef = np.linalg.solve(A, D[tr].T @ y[tr])
    extra.update({"mode": mode, "gamma": gamma, "n_features": D.shape[1],
                  "mean_sampling_stderr": float(np.mean(serrs)),
                  "n_samples_per_window": int(min(n_samp, 64)) if mode == "ensemble" else 0,
                  "label_semantics": "ensemble readout; no latent structure is "
                                     "labelled with the assay outcome"})
    return D @ coef


def prior_sensitivity(rows, cfg, seed, n_windows=24):
    """Does the declared prior's temperature change the ensemble summary?

    Reported separately from numerical approximation error: uncertainty about the
    correct physical prior is NOT covered by a narrow Monte Carlo interval.
    """
    from ..oracles.rna_noncrossing import RNANonCrossingOracle
    from ..rna.ensemble import DeclaredGibbsPrior
    out = {}
    sel = rows[:n_windows]
    for temp in cfg["ensemble"]["prior_sensitivity"]:
        dens = []
        for r in sel:
            o = RNANonCrossingOracle(r["window"], 3, True)
            edges = o.ground_set()
            if not edges:
                continue
            theta = {e: temp * 0.5 for e in edges}
            mu = DeclaredGibbsPrior(o, theta).marginals().mu
            dens.append(float(mu.sum()) / max(len(r["window"]), 1))
        out[str(temp)] = {"mean_pair_density": float(np.mean(dens)) if dens else None,
                          "n_windows": len(dens)}
    vals = [v["mean_pair_density"] for v in out.values() if v["mean_pair_density"]]
    return {"by_temperature": out,
            "spread": float(max(vals) - min(vals)) if vals else None,
            "note": "prior choice moves the ensemble summary by this much; it is "
                    "a modelling uncertainty, not a sampling error"}


def sampling_convergence(rows, cfg, seed, n_windows=12):
    """Report convergence rather than assume it."""
    from ..oracles.rna_noncrossing import RNANonCrossingOracle
    from ..rna.ensemble import DeclaredGibbsPrior
    rng = np.random.RandomState(seed)
    out = {}
    for n in cfg["ensemble"]["n_samples"]:
        means, errs = [], []
        for r in rows[:n_windows]:
            o = RNANonCrossingOracle(r["window"], 3, True)
            if not o.ground_set():
                continue
            p = DeclaredGibbsPrior(o, {e: 0.0 for e in o.ground_set()})
            s = p.sample(n, rng)
            v = np.array([len(x) for x in s], dtype=float)
            means.append(float(v.mean()))
            errs.append(float(v.std(ddof=1) / np.sqrt(len(v))))
        out[str(n)] = {"mean": float(np.mean(means)) if means else None,
                       "mean_stderr": float(np.mean(errs)) if errs else None,
                       "n_windows": len(means)}
    return out


def main(argv=None):
    import yaml
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.rna_pilot",
        description="On-target efficacy and ranking on resolved inputs only.")
    ap.add_argument("--config", default="configs/revision_rna.yaml")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))

    if a.dry_run:
        print("dry run -- nothing fitted\n")
        print(f"windows        {cfg['windows']['size']} nt, "
              f"{cfg['windows']['restrict_to']}")
        print(f"context status {cfg['windows']['context_status']}")
        print(f"assumption     {cfg['windows']['assumption']}")
        print(f"split          {cfg['split']['kind']}, group by "
              f"{cfg['split']['group_by']}")
        print(f"arms           {len(cfg['arms'])}: {', '.join(cfg['arms'])}")
        print(f"seeds          {cfg['seeds']}")
        print(f"metrics        primary {cfg['metrics']['primary_predictive']}; "
              f"selection top-{cfg['metrics']['k']} at threshold "
              f"{cfg['metrics']['activity_threshold']} (prevalence reported)")
        print(f"ensemble       {cfg['ensemble']['objective']}")
        print(f"               prior sensitivity "
              f"{cfg['ensemble']['prior_sensitivity']}, samples "
              f"{cfg['ensemble']['n_samples']}")
        print(f"scope          {cfg['metrics']['scope']}")
        print("\nchemical transfer and native-vs-reporter remain INADMISSIBLE "
              "until Davis S1 is supplied.")
        return 0

    t0 = time.time()
    device = torch.device("cpu")
    rows, binfo = build_rows(cfg, a.limit)
    print(f"admissible rows {binfo['n_admissible']} "
          f"(excluded {binfo['excluded_edge_or_short']} edge/short), "
          f"{binfo['isoform_ambiguous']} isoform-ambiguous\n")
    if len(rows) < 50:
        print("too few admissible rows; abstaining")
        return 1

    results, splits = [], []
    for seed in cfg["seeds"]:
        tr, te, sinfo = novel_sequence_split(rows, cfg["split"]["holdout_fraction"],
                                             seed)
        splits.append(sinfo)
        for arm in cfg["arms"]:
            r, _ = run_arm(arm, rows, tr, te, cfg, seed, device)
            results.append(r)

    print(f"{'arm':<30}{'spearman':>11}{'sd':>8}{'P@20':>8}{'prev':>7}{'enrich':>8}")
    summary = {}
    for arm in cfg["arms"]:
        rs = [r for r in results if r["arm"] == arm]
        sp = np.array([r["spearman"] for r in rs])
        pk = np.array([r["selection"]["precision_at_k"] for r in rs])
        pv = np.array([r["selection"]["prevalence"] for r in rs])
        en = np.array([r["selection"]["enrichment"] for r in rs])
        summary[arm] = {"spearman_mean": float(sp.mean()),
                        "spearman_sd": float(sp.std()),
                        "spearman_per_seed": [float(x) for x in sp],
                        "precision_at_k": float(pk.mean()),
                        "prevalence": float(pv.mean()),
                        "enrichment": float(np.nanmean(en)),
                        "seeds": [r["seed"] for r in rs]}
        print(f"{arm:<30}{sp.mean():>11.4f}{sp.std():>8.4f}"
              f"{pk.mean():>8.3f}{pv.mean():>7.3f}{np.nanmean(en):>8.3f}")

    ps = prior_sensitivity(rows, cfg, cfg["seeds"][0])
    sc = sampling_convergence(rows, cfg, cfg["seeds"][0])
    print(f"\nprior sensitivity (mean pair density): "
          + "  ".join(f"T={k}:{v['mean_pair_density']:.4f}"
                      for k, v in ps["by_temperature"].items()))
    print(f"  spread across priors {ps['spread']:.4f}  <- modelling uncertainty, "
          f"not sampling error")
    print("sampling convergence: "
          + "  ".join(f"n={k}: mean {v['mean']:.3f} se {v['mean_stderr']:.4f}"
                      for k, v in sc.items() if v["mean"] is not None))

    print(f"\nsplit: {splits[0]['n_groups']} groups, "
          f"{splits[0]['n_train']} train / {splits[0]['n_test']} test, "
          f"group leak {splits[0]['group_leak']}")
    print("\nclaims kept separate:")
    for c in cfg["metrics"]["separate_claims"]:
        print(f"  - {c}")
    print(f"scope: {cfg['metrics']['scope']}")

    for r in rows:
        r.pop("_feat", None)
    res = {"kind": "revision_rna_pilot", "config": a.config,
           "config_sha256": sha256_file(a.config),
           "evidence_kind": MEASURED_RNA,
           "window_build": binfo, "context_status": cfg["windows"]["context_status"],
           "context_assumption": cfg["windows"]["assumption"],
           "splits": splits, "summary": summary, "rows": results,
           "prior_sensitivity": ps, "sampling_convergence": sc,
           "inadmissible": ["chemical_transfer", "native_vs_reporter"],
           "runtime_s": round(time.time() - t0, 2),
           "cpu_core_hours": round((time.time() - t0) / 3600, 5)}
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir,
                        f"rna_pilot-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    with open(path + ".sha256", "w") as fh:
        fh.write(sha256_file(path) + "  " + os.path.basename(path) + "\n")
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
