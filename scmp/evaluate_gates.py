"""Preregistered gate table with confidence intervals.

Gates, thresholds and the split unit come from `configs/protocol.yaml`; none is chosen
here. A gate that fails is reported as failed. There is no path in this module that
selects seeds, retunes on the evaluation data, invents a gain, or converts a failure into
a negative-results framing -- a failed gate emits a redesign requirement instead.

Run:  python -m scmp.evaluate_gates --config configs/protocol.yaml
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np
import yaml


def cluster_bootstrap(values, clusters, stat=np.median, n=10000, level=0.95, seed=0):
    """Percentile CI, resampling whole clusters so within-cluster correlation is not
    mistaken for independent evidence."""
    values = np.asarray(values, dtype=float)
    clusters = np.asarray(clusters)
    uniq = np.unique(clusters)
    if len(uniq) < 2 or len(values) == 0:
        return float(stat(values)) if len(values) else float("nan"), (float("nan"),) * 2
    idx = {c: np.where(clusters == c)[0] for c in uniq}
    rng = np.random.default_rng(seed)
    draws = np.empty(n)
    for b in range(n):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        sel = np.concatenate([idx[c] for c in pick])
        draws[b] = stat(values[sel])
    lo, hi = np.percentile(draws, [(1 - level) / 2 * 100, (1 + level) / 2 * 100])
    return float(stat(values)), (float(lo), float(hi))


def paired(records, arm_a, arm_b, key="gap_rel"):
    """Pair arm_a against arm_b on identical (family, instance, seed)."""
    by = {}
    for r in records:
        if r.get("kind") != "verifier":
            continue
        by.setdefault((r["family"], r["instance"], r["seed"]), {})[r["arm"]] = r
    a_vals, b_vals, clusters = [], [], []
    for k, arms in by.items():
        if arm_a in arms and arm_b in arms:
            a_vals.append(arms[arm_a][key])
            b_vals.append(arms[arm_b][key])
            clusters.append(f"{k[0]}#{k[1]}")
    return np.array(a_vals), np.array(b_vals), np.array(clusters)


def relative_reduction(new, old):
    """(median(old) - median(new)) / median(old), guarded for a zero baseline."""
    mo = np.median(old)
    if abs(mo) < 1e-12:
        return float("nan")
    return (mo - np.median(new)) / abs(mo)


def gate_verifier(records, arm_new, arm_old, mmg, label, seed=0):
    a, b, cl = paired(records, arm_new, arm_old)
    if len(a) == 0:
        return {"gate": label, "status": "NO DATA", "n": 0}
    if abs(np.median(b)) < 1e-12:
        return {"gate": label, "n_pairs": int(len(a)),
                "median_new": float(np.median(a)), "median_old": float(np.median(b)),
                "status": "UNDETERMINED",
                "reason": "the baseline median is already 0 at matched wall time, so a "
                          "relative reduction is undefined; the instances are too easy "
                          "to discriminate the arms",
                "rule": "relative reduction requires a non-zero baseline"}
    point = relative_reduction(a, b)
    idx = np.arange(len(a))

    def stat(sel):
        return relative_reduction(a[sel], b[sel])

    uniq = np.unique(cl)
    rng = np.random.default_rng(seed)
    draws = []
    cidx = {c: np.where(cl == c)[0] for c in uniq}
    for _ in range(10000):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        sel = np.concatenate([cidx[c] for c in pick])
        draws.append(stat(sel))
    lo, hi = np.percentile(draws, [2.5, 97.5])
    passed = bool(lo > mmg)
    return {"gate": label, "n_pairs": int(len(a)), "n_clusters": int(len(uniq)),
            "median_new": float(np.median(a)), "median_old": float(np.median(b)),
            "relative_reduction": float(point), "ci95": [float(lo), float(hi)],
            "threshold": mmg, "status": "PASS" if passed else "FAIL",
            "rule": "pass iff lower 95% CI of the relative reduction exceeds the "
                    "minimum meaningful gain"}


def gate_noninferiority(records, arm_constrained, arm_free, margin, label, seed=0):
    rows = [r for r in records if r.get("kind") == "predictor"]
    c = {r["seed"]: r["spearman"] for r in rows if r["arm"] == arm_constrained}
    f = {r["seed"]: r["spearman"] for r in rows if r["arm"] == arm_free}
    common = sorted(set(c) & set(f))
    if not common:
        return {"gate": label, "status": "NO DATA", "n": 0}
    diff = np.array([f[s] - c[s] for s in common])
    point, (lo, hi) = cluster_bootstrap(diff, np.array(common), stat=np.mean, seed=seed)
    passed = bool(hi < margin)
    return {"gate": label, "n_seeds": len(common), "mean_difference": float(point),
            "ci95": [float(lo), float(hi)], "margin": margin,
            "status": "PASS" if passed else "FAIL",
            "rule": "noninferiority: pass iff the UPPER 95% CI of "
                    "(unconstrained - constrained) lies below the margin"}


def gate_positive_effect(records, arm_on, arm_off, label, seed=0):
    rows = [r for r in records if r.get("kind") == "predictor"]
    on = {r["seed"]: r["spearman"] for r in rows if r["arm"] == arm_on}
    off = {r["seed"]: r["spearman"] for r in rows if r["arm"] == arm_off}
    common = sorted(set(on) & set(off))
    if not common:
        return {"gate": label, "status": "NO DATA", "n": 0}
    diff = np.array([on[s] - off[s] for s in common])
    point, (lo, hi) = cluster_bootstrap(diff, np.array(common), stat=np.mean, seed=seed)
    passed = bool(lo > 0.0)
    return {"gate": label, "n_seeds": len(common), "mean_difference": float(point),
            "ci95": [float(lo), float(hi)],
            "status": "PASS" if passed else "FAIL",
            "rule": "pass iff the LOWER 95% CI of (residual on - residual off) "
                    "exceeds zero"}


REDESIGN = {
    "G1 conditional vs unconditional (UNDETERMINED)":
        "The pilot families are too easy: the unconditional baseline already closes "
        "98-100% of instances to gap exactly 0 within the matched wall time, so there "
        "is nothing for conditioning to improve and the metric is 0/0. Redesign means "
        "harder instances (longer windows, larger k), not a different metric.",
    "G1 conditional vs unconditional":
        "The conditioning mechanism is the contribution. If it does not clear the "
        "minimum meaningful gain at matched wall time, the bound is not worth its "
        "oracle cost and the relaxation must be tightened structurally "
        "(docs/bounds_proof.md, the nesting constraint the edge lattice discards) "
        "rather than by spending more queries.",
    "G2 learned vs deterministic allocation":
        "Learned allocation was already at the noise floor in development "
        "(random 0.0830 vs largest-gap 0.0830). If it fails here too, the honest move "
        "is to drop the learned allocator from the contribution and reposition on the "
        "conditioning itself. Do not retune the scorer on this evaluation data.",
    "G3 predictor noninferiority":
        "If the sign-constrained predictor is not noninferior, the certifiable model "
        "class costs real accuracy and that cost must be reported as the headline, not "
        "buried. Redesign means a richer non-negative encoding, not a weaker baseline.",
    "G4 residual predictive value":
        "If the nonlinear residual adds no held-out value on this label, note that the "
        "label is marginal over latent structures (preregistration section 6) and the "
        "test may be unable to detect a structure-conditional effect at all. Redesign "
        "means an ensemble-averaged objective or a structure-resolved label set, not a "
        "larger residual.",
}


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.evaluate_gates")
    ap.add_argument("--config", default="configs/protocol.yaml")
    ap.add_argument("--run", default="runs/pilot")
    a = ap.parse_args(argv)
    proto = yaml.safe_load(open(a.config))
    mmg = proto["effect_sizes"]["minimum_meaningful_gain"]["relative_reduction"]
    margin = proto["effect_sizes"]["predictive_noninferiority_margin"]["margin"]

    files = sorted(glob.glob(os.path.join(a.run, "pilot-*.json")))
    if not files:
        print(f"no pilot log under {a.run}")
        return 1
    records = json.load(open(files[-1]))
    meta = next((r for r in records if r.get("kind") == "meta"), {})
    final = next((r for r in records if r.get("kind") == "final"), {})

    print(f"preregistered gate evaluation")
    print(f"  protocol : {a.config}")
    print(f"  pilot log: {files[-1]}")
    hw = meta.get("hardware", {})
    print(f"  hardware : device={hw.get('device')} gpu={hw.get('gpu_name')} "
          f"index={hw.get('gpu_index')} frac={hw.get('gpu_memory_fraction')}")
    print(f"  budget   : {final.get('core_hours', float('nan')):.4f} core-hours used, "
          f"within={final.get('within_budget')}")
    nv = sum(1 for r in records if r.get("kind") == "verifier")
    nto = sum(1 for r in records if r.get("timeout"))
    nvoid = sum(1 for r in records if r.get("kind") == "void")
    print(f"  records  : {nv} verifier, {nto} timeouts, {nvoid} void (excluded), "
          f"{sum(1 for r in records if r.get('kind')=='predictor')} predictor")
    print(f"  matched wall time per arm: "
          f"{meta.get('config', {}).get('matched_wall_time_s')} s\n")

    gates = [
        gate_verifier(records, "fixed_conditional", "intermediate_unconditional",
                      mmg, "G1 conditional vs unconditional"),
        gate_verifier(records, "learned_conditional", "fixed_conditional",
                      mmg, "G2 learned vs deterministic allocation"),
        gate_noninferiority(records, "certifiable_gamma1", "unconstrained",
                            margin, "G3 predictor noninferiority"),
        gate_positive_effect(records, "certifiable_gamma1", "certifiable_gamma0",
                             "G4 residual predictive value"),
    ]

    print(f"{'gate':>38} {'point':>10} {'95% CI':>22} {'thresh':>8} {'status':>12}")
    for g in gates:
        if g["status"] in ("NO DATA", "UNDETERMINED"):
            print(f"{g['gate']:>38} {'-':>10} {'-':>22} {'-':>8} "
                  f"{g['status']:>12}")
            continue
        pt = g.get("relative_reduction", g.get("mean_difference"))
        ci = g["ci95"]
        th = g.get("threshold", g.get("margin", 0.0))
        print(f"{g['gate']:>38} {pt:10.4f} [{ci[0]:9.4f},{ci[1]:9.4f}] "
              f"{th:8.3f} {g['status']:>12}")

    failed = [g for g in gates if g["status"] in ("FAIL", "NO DATA", "UNDETERMINED")]
    print()
    for g in gates:
        print(f"  {g['gate']}: {g.get('rule','')}")
    if failed:
        print(f"\nFAILED GATES: {len(failed)}")
        for g in failed:
            print(f"\n  {g['gate']} -> {g['status']}")
            print(f"    REDESIGN REQUIRED: {REDESIGN.get(g['gate'], '')}")
        print("\nNo seed selection, no retuning on evaluation data, no invented gain,")
        print("and no negative-results pivot. The failing claim is withdrawn until a")
        print("redesign is proposed and preregistered.")

    out = "runs/gates"
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"gates-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump({"protocol": a.config, "pilot_log": files[-1], "gates": gates,
                   "hardware": hw, "budget": final,
                   "numerical_status": "float_unverified"}, fh, indent=2)
    print(f"\nwritten: {path}")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
