"""RNA application evaluation: arm admissibility, then the comparisons that survive it.

ADMISSIBILITY IS CHECKED FIRST. Each arm declares the inputs it needs. An arm whose
inputs are not available for the dataset is reported as BLOCKED with the reason, and is
never quietly dropped or silently fed a substitute. For huesken2005 the assayed context
is a YFP reporter construct that is not in hand, so every arm needing a target window is
blocked; running them on native transcript windows would answer a different experiment.

WHAT IS AND IS NOT CERTIFIED. For lower-is-better mean intervals, `u_i < l_j` certifies
that the MODEL ranks candidate i above candidate j -- it says nothing about the measured
ordering. Agreement with held-out measurement is computed separately, reported
separately, and is an empirical quantity with a confidence interval, not a certificate.

Assay labels and scales are kept apart: every metric is computed within an assay and
never pooled across assays.

Run:  python -m scmp.rna.evaluate --protocol configs/rna_protocol.yaml
"""
from __future__ import annotations

import argparse, csv, hashlib, json, os, random, sys, time
from collections import defaultdict

import numpy as np
import yaml

from ..oracles import RNANonCrossingOracle

NUC = "ACGU"


# ------------------------------------------------------------------ availability
def dataset_capabilities(manifest_dir):
    """What the audited manifests say is actually available."""
    import glob
    fs = sorted(glob.glob(os.path.join(manifest_dir, "manifest-*.json")))
    if not fs:
        return {"available": set(), "reason": "no data manifest found"}
    man = json.load(open(fs[-1]))["manifests"]["huesken2005"]
    caps = {"guide", "target_site", "family"}
    ctx = man["context_resolution"]
    reason = ""
    if ctx["status"] == "UNRESOLVED":
        reason = ctx["reason"]
    else:
        caps |= {"target_window", "folding", "marginals", "sampling"}
    return {"available": caps, "reason": reason, "assay": man["assay"],
            "n_records": man["n_records"], "context": ctx}


def load_rows():
    tr = [l.split() for l in open("data/TrainAll2182.txt") if l.strip()]
    te = [l.split() for l in open("data/TestAll249.txt") if l.strip()]
    rows = [{"seq": r[0].upper().replace("T", "U"), "y": float(r[1])}
            for r in tr + te]
    for r in rows:
        r["assay"] = "reporter_yfp"          # single assay; never pooled with another
        r["target_site"] = r["seq"].translate(str.maketrans("ACGU", "UGCA"))[::-1]
    return rows


def kmer_clusters(seqs, k=13):
    parent = list(range(len(seqs)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    idx = defaultdict(list)
    for i, s in enumerate(seqs):
        for p in range(len(s) - k + 1):
            idx[s[p:p + k]].append(i)
    for _, mem in idx.items():
        for m in mem[1:]:
            a, b = find(mem[0]), find(m)
            if a != b:
                parent[a] = b
    return [find(i) for i in range(len(seqs))]


def onehot(seq):
    v = np.zeros(4 * len(seq))
    for t, c in enumerate(seq):
        if c in NUC:
            v[4 * t + NUC.index(c)] = 1.0
    return v


def ridge_fit(Xtr, ytr, lam):
    A = Xtr.T @ Xtr + lam * np.eye(Xtr.shape[1])
    return np.linalg.solve(A, Xtr.T @ ytr)


def spearman(a, b):
    def rank(v):
        v = np.asarray(v, float); o = v.argsort()
        r = np.empty(len(v)); r[o] = np.arange(len(v), dtype=float)
        _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
        for g in np.where(cnt > 1)[0]:
            m = inv == g; r[m] = r[m].mean()
        return r
    ra, rb = rank(a) - np.mean(rank(a)), rank(b) - np.mean(rank(b))
    d = np.sqrt((ra * ra).sum()) * np.sqrt((rb * rb).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")


def precision_at_k(pred, meas, k):
    """Lower measured value = less residual expression = better guide, so the best
    candidates are the SMALLEST measured values... here the label is inhibition, where
    LARGER is better. Stated explicitly rather than assumed."""
    k = min(k, len(pred))
    top_pred = set(np.argsort(-np.asarray(pred))[:k])
    top_meas = set(np.argsort(-np.asarray(meas))[:k])
    return len(top_pred & top_meas) / k


def cluster_bootstrap_diff(d, clusters, n=10000, seed=0):
    d = np.asarray(d, float); clusters = np.asarray(clusters)
    uniq = np.unique(clusters)
    if len(uniq) < 2:
        return float(np.mean(d)), (float("nan"), float("nan"))
    idx = {c: np.where(clusters == c)[0] for c in uniq}
    rng = np.random.default_rng(seed)
    draws = np.empty(n)
    for b in range(n):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        sel = np.concatenate([idx[c] for c in pick])
        draws[b] = d[sel].mean()
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(np.mean(d)), (float(lo), float(hi))


# ------------------------------------------------------- resolved comparisons demo
def resolved_comparison_demo(seed=0):
    """Engineering check of the u_i < l_j rule on TINY DECLARED windows.

    This is NOT an RNA result: the windows are declared sequences, not assayed context.
    It exists to show the rule is implemented correctly and to separate the two claims
    it can support -- certified model ranking, and empirical agreement with measurement.
    """
    import torch
    from ..model import CertifiedModel, backbone_adjacency
    from .ensemble import DeclaredGibbsPrior, ensemble_bracket, mean_bracket
    from .features import build_covariates, per_node_covariates
    from ..oracles.bruteforce import enumerate_rna

    seq = "GCGCAAAAGCGC"
    o = RNANonCrossingOracle(seq, 3, True)
    n = len(seq)
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    theta = {e: 0.7 * ((k % 5) - 2) / 2.0 for k, e in enumerate(o.ground_set())}
    prior = DeclaredGibbsPrior(o, theta)
    mu = prior.marginals().mu
    c = prior.coefficients()
    lw = np.array([sum(c.get(e, 0.0) for e in m) for m in members])
    p = np.exp(lw - lw.max()); p /= p.sum()

    cands = []
    for cand in range(12):
        X = build_covariates(("ACGU" * 6)[:21], "unmodified", "reporter_yfp",
                             ["unmodified"], ["reporter_yfp"])
        Xn = torch.tensor(per_node_covariates(X, n), dtype=torch.float64)
        m = CertifiedModel(Xn.shape[1], o.ground_set(), n, d_hid=4, gamma=1.0,
                           seed=seed + cand)
        g = torch.Generator().manual_seed(seed + cand + 100)
        with torch.no_grad():
            for pp in m.parameters():
                pp.add_(torch.randn(pp.shape, generator=g, dtype=pp.dtype))
        B = backbone_adjacency(n)
        br = mean_bracket(m, Xn, B, mu)
        lo, hi = ensemble_bracket(m, Xn, mu, br)
        with torch.no_grad():
            fs = np.array([float(m(Xn, m.indicator(s), B)) for s in members])
        cands.append({"id": cand, "l": lo, "u": hi, "true_mean": float(p @ fs)})

    pairs, resolved, correct = 0, 0, 0
    for i in range(len(cands)):
        for j in range(len(cands)):
            if i == j:
                continue
            pairs += 1
            if cands[i]["u"] < cands[j]["l"]:
                resolved += 1
                if cands[i]["true_mean"] < cands[j]["true_mean"]:
                    correct += 1
    return {"n_candidates": len(cands), "pairs": pairs, "resolved": resolved,
            "coverage": resolved / pairs if pairs else float("nan"),
            "model_ranking_agreement": (correct / resolved) if resolved else float("nan"),
            "median_width": float(np.median([c["u"] - c["l"] for c in cands])),
            "note": "declared windows; NOT an RNA result; certifies model ranking only"}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.rna.evaluate")
    ap.add_argument("--protocol", default="configs/rna_protocol.yaml")
    ap.add_argument("--manifest-dir", default="runs/rna_data")
    a = ap.parse_args(argv)
    proto = yaml.safe_load(open(a.protocol))
    t0 = time.time()

    print(f"RNA application evaluation: {proto['protocol']['name']}")
    print(f"task     : {proto['protocol']['task']}")
    print(f"primary  : {proto['metrics']['predictive_primary']['name']} "
          f"({proto['metrics']['predictive_primary']['scope']})")
    print(f"selection: {proto['metrics']['candidate_selection']['name']}"
          f"@{proto['metrics']['candidate_selection']['k']}\n")

    caps = dataset_capabilities(a.manifest_dir)
    print(f"dataset capabilities: {sorted(caps['available'])}")
    if caps.get("reason"):
        print(f"  blocked because: {caps['reason']}\n")

    # ---------------- arm admissibility ------------------------------
    admissible, blocked = [], []
    for arm in proto["arms"]:
        missing = [x for x in arm["needs"] if x not in caps["available"]]
        (blocked if missing else admissible).append(
            {**arm, "missing": missing})
    print(f"{'arm':>24} {'status':>10}  needs / missing")
    for arm in proto["arms"]:
        m = [x for x in arm["needs"] if x not in caps["available"]]
        print(f"{arm['name']:>24} {'BLOCKED' if m else 'admissible':>10}  "
              f"{','.join(arm['needs'])}"
              f"{'  missing: ' + ','.join(m) if m else ''}")
    print(f"\nadmissible: {len(admissible)}/{len(proto['arms'])}")

    # ---------------- run admissible arms ----------------------------
    rows = load_rows()
    assays = sorted({r["assay"] for r in rows})
    print(f"assays present: {assays} (metrics are computed within an assay, never "
          f"pooled)")
    cl = kmer_clusters([r["seq"] for r in rows])
    for r, c in zip(rows, cl):
        r["cluster"] = c

    results, excl = {}, {"non_finite": 0}
    for assay in assays:
        sub = [r for r in rows if r["assay"] == assay]
        clusters = sorted({r["cluster"] for r in sub})
        for arm in admissible:
            per_seed = []
            for seed in proto["design"]["seeds"]:
                rng = random.Random(seed)
                ks = list(clusters); rng.shuffle(ks)
                held = set(ks[:max(1, len(ks) // 4)])
                tr = [r for r in sub if r["cluster"] not in held]
                te = [r for r in sub if r["cluster"] in held]
                if arm["name"] == "linear_sequence_only":
                    feat = lambda r: onehot(r["seq"])
                else:
                    feat = lambda r: np.concatenate([onehot(r["seq"]),
                                                     onehot(r["target_site"])])
                Xtr = np.stack([feat(r) for r in tr]); ytr = np.array([r["y"] for r in tr])
                Xte = np.stack([feat(r) for r in te]); yte = np.array([r["y"] for r in te])
                w = ridge_fit(np.hstack([Xtr, np.ones((len(Xtr), 1))]), ytr, 1.0)
                pred = np.hstack([Xte, np.ones((len(Xte), 1))]) @ w
                if not np.all(np.isfinite(pred)):
                    excl["non_finite"] += 1
                    continue
                per_seed.append({"seed": seed, "n_test": len(te),
                                 "spearman": spearman(pred, yte),
                                 "precision_at_k": precision_at_k(
                                     pred, yte, proto["metrics"]
                                     ["candidate_selection"]["k"]),
                                 "clusters_held": len(held),
                                 "pred": pred.tolist(), "y": yte.tolist(),
                                 "cluster": [r["cluster"] for r in te]})
            results[(assay, arm["name"])] = per_seed

    print(f"\npredictive results, within assay (all seeds reported)")
    print(f"{'assay':>14} {'arm':>24} {'seeds':>6} {'Spearman':>18} "
          f"{'P@20':>16}")
    summary = {}
    for (assay, name), rs in results.items():
        sp = [r["spearman"] for r in rs]; pk = [r["precision_at_k"] for r in rs]
        summary[f"{assay}|{name}"] = {"spearman": sp, "precision_at_k": pk,
                                      "seeds": [r["seed"] for r in rs]}
        print(f"{assay:>14} {name:>24} {len(rs):6d} "
              f"{np.mean(sp):8.4f} +/-{np.std(sp):6.4f} "
              f"{np.mean(pk):8.4f} +/-{np.std(pk):5.4f}")

    # paired comparison between the two admissible arms
    keys = [k for k in results if k[1] == "sequence_target"]
    paired = None
    if keys and (keys[0][0], "linear_sequence_only") in results:
        assay = keys[0][0]
        A = results[(assay, "linear_sequence_only")]
        Bm = results[(assay, "sequence_target")]
        d, cls = [], []
        for ra, rb in zip(A, Bm):
            d.append(rb["spearman"] - ra["spearman"]); cls.append(ra["seed"])
        pt, ci = cluster_bootstrap_diff(d, cls,
                                        n=proto["design"]["confidence_interval"]
                                        ["resamples"])
        paired = {"contrast": "sequence_target - linear_sequence_only",
                  "mean_difference": pt, "ci95": list(ci), "n_seeds": len(d)}
        print(f"\npaired contrast: {paired['contrast']}")
        print(f"  mean difference {pt:+.4f}  95% CI [{ci[0]:+.4f}, {ci[1]:+.4f}]")
        print("  note: on this dataset the target site is the reverse complement of the "
              "guide, so the two arms carry the same information by construction")

    # ---------------- learned-layer / residual checks -----------------
    print(f"\nlearned-layer benefit and residual nonlinearity")
    print("  BLOCKED: both require the A+B arms, which need a target window")

    # ---------------- resolved comparisons (engineering) --------------
    print(f"\nresolved-comparison rule (engineering check, declared windows, "
          f"NOT an RNA result)")
    rc = resolved_comparison_demo()
    print(f"  candidates {rc['n_candidates']}, pairs {rc['pairs']}, "
          f"resolved {rc['resolved']} ({rc['coverage']:.1%})")
    print(f"  median ensemble-bracket width {rc['median_width']:.4f}")
    print(f"  resolved pairs whose MODEL ordering is correct: "
          f"{rc['model_ranking_agreement'] if rc['resolved'] else 'n/a'}")
    print("  agreement with held-out MEASUREMENT: not computable -- no measured "
          "ensemble readout exists for these declared windows")

    print(f"\nseparation of claims")
    for k, v in proto["separation"].items():
        print(f"  {k:>28}: {v}")

    out = "runs/rna_eval"
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"eval-{time.strftime('%Y%m%d-%H%M%S')}.json")
    payload = json.dumps({"protocol": a.protocol, "capabilities": sorted(caps["available"]),
                          "blocked_reason": caps.get("reason"),
                          "admissible": [x["name"] for x in admissible],
                          "blocked": [{"arm": x["name"], "missing": x["missing"]}
                                      for x in blocked],
                          "summary": summary, "paired": paired,
                          "resolved_comparison_demo": rc,
                          "exclusions": excl, "timeouts": 0,
                          "runtime_s": time.time() - t0}, indent=2, default=str)
    open(path, "w").write(payload)
    open(path + ".sha256", "w").write(hashlib.sha256(payload.encode()).hexdigest() + "\n")
    print(f"\nwritten: {path}")
    return 1 if blocked else 0


if __name__ == "__main__":
    sys.exit(main())
